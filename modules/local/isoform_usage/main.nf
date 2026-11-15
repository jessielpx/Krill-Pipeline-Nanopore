/*
 * ============================================================
 * Krill Renew - 07 ISOFORM_USAGE
 *
 * Per-sample isoform fractions and delta IF
 *
 * Inputs:
 *   - Bambu uniqueCounts_transcript.txt
 *   - Samplesheet CSV / TSV
 *   - run_isoform_usage.R
 *
 * Modes:
 *   pilot  - original all-sample coverage filtering
 *   cohort - group-aware coverage filtering
 *
 * Optional Nextflow parameters:
 *   --isoform_mode pilot|cohort
 *   --isoform_min_gene_counts 10
 *   --isoform_group_fraction 0.25
 *
 * Output:
 *   07_isoform_usage/isoform_usage_results
 *
 * Native execution on Rorqual.
 * No containers.
 * ============================================================
 */

process ISOFORM_USAGE {

    tag "cohort"

    cpus 2
    memory '16 GB'
    time '4h'

    publishDir "${params.outdir}/07_isoform_usage",
        mode: 'copy'

    input:
    path bambu_results
    path samplesheet
    path usage_script

    output:
    path 'isoform_usage_results', emit: results

    shell:
    '''
    set -euo pipefail

    # ========================================================
    # 1. Load environment
    # ========================================================

    module --force purge

    module load StdEnv/2023
    module load r/4.6.1
    module load python/3.11.5

    unset PYTHONHOME PYTHONPATH
    hash -r

    export R_LIBS_USER="!{params.r_libs_user ?: ''}"

    export BIOCONDUCTOR_CONFIG_FILE="!{params.bioconductor_config_file ?: ''}"

    if [[ -z "$R_LIBS_USER" ]]; then
        echo "ERROR: Set r_libs_user in params.yaml" >&2
        exit 1
    fi

    if [[ ! -d "$R_LIBS_USER" ]]; then
        echo "ERROR: R library directory not found: $R_LIBS_USER" >&2
        exit 1
    fi

    if [[ -z "${BIOCONDUCTOR_CONFIG_FILE:-}" ]]; then
        echo "ERROR: Set bioconductor_config_file in params.yaml" >&2
        exit 1
    fi

    BIOC_CONFIG_PATH="${BIOCONDUCTOR_CONFIG_FILE#file://}"

    if [[ ! -s "$BIOC_CONFIG_PATH" ]]; then
        echo "ERROR: Bioconductor config not found: $BIOC_CONFIG_PATH" >&2
        exit 1
    fi


    # ========================================================
    # 2. Configure isoform analysis
    # ========================================================

    export KRILL_ISOFORM_MODE='!{params.isoform_mode ?: "pilot"}'

    export KRILL_ISOFORM_GROUP_FRACTION='!{params.isoform_group_fraction != null ? params.isoform_group_fraction : 0.25}'

    MIN_GENE_COUNTS='!{params.isoform_min_gene_counts != null ? params.isoform_min_gene_counts : 10}'

    echo "========================================"
    echo "Krill Renew - ISOFORM_USAGE"
    echo "Mode: $KRILL_ISOFORM_MODE"
    echo "Minimum gene counts: $MIN_GENE_COUNTS"
    echo "Group fraction: $KRILL_ISOFORM_GROUP_FRACTION"
    echo "========================================"

    if [[ "$KRILL_ISOFORM_MODE" != "pilot" &&
          "$KRILL_ISOFORM_MODE" != "cohort" ]]; then

        echo "ERROR: Invalid ISOFORM_USAGE mode" >&2
        exit 1
    fi

    # ========================================================
    # 3. Validate inputs
    # ========================================================

    COUNTS="!{bambu_results}/uniqueCounts_transcript.txt"
    SAMPLESHEET="!{samplesheet}"
    SCRIPT="!{usage_script}"

    for required in \
        "$COUNTS" \
        "$SAMPLESHEET" \
        "$SCRIPT"
    do
        if [[ ! -s "$required" ]]; then
            echo "ERROR: Missing input: $required" >&2
            exit 1
        fi
    done

    # ========================================================
    # 4. Prepare metadata
    #
    # Preserve all samplesheet columns:
    #   sample, condition, run, subtype, patient, etc.
    #
    # The R script currently uses sample and condition.
    # Keeping all columns allows future extensions.
    # ========================================================

    python3 - "$SAMPLESHEET" <<'PY'
import csv
import sys

source = sys.argv[1]

with open(source, newline="") as handle:
    first_line = handle.readline()

delimiter = chr(9) if chr(9) in first_line else ","

with open(source, newline="") as handle:
    reader = csv.DictReader(
        handle,
        delimiter=delimiter
    )

    columns = reader.fieldnames or []

    if not {"sample", "condition"}.issubset(columns):
        raise SystemExit(
            "ERROR: samplesheet requires "
            "sample and condition"
        )

    rows = list(reader)

if not rows:
    raise SystemExit("ERROR: Empty samplesheet")

if any(None in row for row in rows):
    raise SystemExit(
        "ERROR: Invalid samplesheet column count"
    )

samples = []
conditions = []

for row in rows:

    sample = (row.get("sample") or "").strip()
    condition = (row.get("condition") or "").strip()

    if not sample:
        raise SystemExit("ERROR: Empty sample ID")

    if condition not in {"control", "HR"}:
        raise SystemExit(
            f"ERROR: Invalid condition: {condition}"
        )

    samples.append(sample)
    conditions.append(condition)

if len(samples) != len(set(samples)):
    raise SystemExit("ERROR: Duplicate sample IDs")

for condition in ("control", "HR"):

    if conditions.count(condition) < 2:
        raise SystemExit(
            f"ERROR: At least two {condition} "
            "samples required"
        )

with open("metadata.tsv", "w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=columns,
        delimiter=chr(9),
        extrasaction="ignore"
    )

    writer.writeheader()

    for row in rows:
        writer.writerow({
            col: (row.get(col) or "").strip()
            for col in columns
        })

print("===== ISOFORM METADATA =====")
print("Samples:", len(samples))
print("HR:", conditions.count("HR"))
print("Control:", conditions.count("control"))
print("Columns:", ", ".join(columns))
PY

    # ========================================================
    # 5. Run ISOFORM_USAGE
    # ========================================================

    echo "===== RUNNING ISOFORM_USAGE ====="

    Rscript --vanilla "$SCRIPT" \
        "$COUNTS" \
        metadata.tsv \
        isoform_usage_results \
        "$MIN_GENE_COUNTS"

    # ========================================================
    # 6. Validate outputs
    # ========================================================

    for output in \
        isoform_usage_results/isoform_usage_supported.tsv \
        isoform_usage_results/isoform_usage_gene_summary.tsv \
        isoform_usage_results/isoform_gene_coverage.tsv \
        isoform_usage_results/isoform_usage_qc.tsv \
        isoform_usage_results/R_session_info.txt
    do

        if [[ ! -s "$output" ]]; then
            echo "ERROR: Missing output: $output" >&2
            exit 1
        fi

    done

    echo "========================================"
    echo "Krill ISOFORM_USAGE completed"
    echo "Mode: $KRILL_ISOFORM_MODE"
    echo "========================================"

    '''
}
