/*
 * ============================================================
 * Krill Renew - 06 DTU
 *
 * Differential Transcript Usage using DRIMSeq
 *
 * Input:
 *   - Bambu uniqueCounts_transcript.txt
 *   - Samplesheet (CSV or TSV)
 *   - run_dtu.R
 *
 * Supported modes:
 *   pilot  - original Batch17-compatible filtering
 *   cohort - group-aware coverage filtering
 *
 * Optional Nextflow parameters:
 *   --dtu_mode pilot|cohort
 *   --dtu_design "~ group"
 *   --dtu_coverage_min 10
 *   --dtu_group_fraction 0.25
 *
 * Metadata:
 *   Required: sample, condition
 *   Optional: run, subtype, patient, etc.
 *
 * Native execution on Rorqual.
 * No containers.
 * ============================================================
 */

process DTU {

    tag "cohort"

    cpus 4
    memory '32 GB'
    time '12h'

    publishDir "${params.outdir}/06_dtu",
        mode: 'copy'

    input:
    path bambu_results
    path samplesheet
    path dtu_script

    output:
    path 'dtu_results', emit: results

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
    # 2. DTU configuration
    # ========================================================

    export KRILL_DTU_MODE='!{params.dtu_mode ?: "pilot"}'

    export KRILL_DTU_DESIGN='!{params.dtu_design ?: "~ group"}'

    export KRILL_DTU_COVERAGE_MIN='!{params.dtu_coverage_min != null ? params.dtu_coverage_min : 10}'

    export KRILL_DTU_GROUP_FRACTION='!{params.dtu_group_fraction != null ? params.dtu_group_fraction : 0.25}'

    echo "========================================"
    echo "Krill Renew - DTU"
    echo "Mode: $KRILL_DTU_MODE"
    echo "Design: $KRILL_DTU_DESIGN"
    echo "Minimum gene coverage: $KRILL_DTU_COVERAGE_MIN"
    echo "Group fraction: $KRILL_DTU_GROUP_FRACTION"
    echo "========================================"

    if [[ "$KRILL_DTU_MODE" != "pilot" &&
          "$KRILL_DTU_MODE" != "cohort" ]]; then

        echo "ERROR: Invalid DTU mode" >&2
        exit 1
    fi

    # ========================================================
    # 3. Validate inputs
    # ========================================================

    COUNTS="!{bambu_results}/uniqueCounts_transcript.txt"
    SAMPLESHEET="!{samplesheet}"
    SCRIPT="!{dtu_script}"

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

    Rscript --vanilla -e '
        stopifnot(
            requireNamespace("DRIMSeq", quietly=TRUE),
            requireNamespace("BiocParallel", quietly=TRUE)
        )

        cat(
            "DRIMSeq version:",
            as.character(packageVersion("DRIMSeq")),
            "\\n"
        )
    '

    # ========================================================
    # 4. Prepare metadata
    #
    # Preserve ALL samplesheet columns.
    #
    # The original module only exported:
    #   sample, condition
    #
    # This version also retains:
    #   run, subtype, patient, etc.
    #
    # run_dtu.R uses selected covariates according
    # to KRILL_DTU_DESIGN.
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
            "ERROR: samplesheet must contain "
            "sample and condition"
        )

    rows = list(reader)

if not rows:
    raise SystemExit("ERROR: Empty samplesheet")

if any(None in row for row in rows):
    raise SystemExit(
        "ERROR: Samplesheet has more fields "
        "than its header"
    )

samples = []
conditions = []

for row in rows:
    sample = (row.get("sample") or "").strip()
    condition = (row.get("condition") or "").strip()

    if not sample:
        raise SystemExit("ERROR: Empty sample ID")

    if condition not in {"HR", "control"}:
        raise SystemExit(
            f"ERROR: Invalid condition: {condition}"
        )

    samples.append(sample)
    conditions.append(condition)

if len(samples) != len(set(samples)):
    raise SystemExit("ERROR: Duplicate sample IDs")

for condition in ("HR", "control"):
    if conditions.count(condition) < 2:
        raise SystemExit(
            "ERROR: At least two samples required "
            f"for {condition}"
        )

# Write TSV while preserving additional metadata.

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
            column: (
                row.get(column) or ""
            ).strip()
            for column in columns
        })

print("===== DTU METADATA =====")
print("Samples:", len(samples))
print("HR:", conditions.count("HR"))
print("Control:", conditions.count("control"))
print("Columns:", ", ".join(columns))
print("Metadata written to metadata.tsv")
PY

    # ========================================================
    # 5. Run DRIMSeq
    # ========================================================

    mkdir -p dtu_results

    Rscript --vanilla "$SCRIPT" \
        "$COUNTS" \
        metadata.tsv \
        dtu_results

    # ========================================================
    # 6. Validate output files
    # ========================================================

    for output in \
        dtu_results/dtu_gene_results.tsv \
        dtu_results/dtu_transcript_results.tsv \
        dtu_results/dtu_filtered_counts.tsv \
        dtu_results/dtu_gene_coverage_qc.tsv \
        dtu_results/dtu_input_qc.tsv \
        dtu_results/dtu_design_matrix.tsv \
        dtu_results/dtu_samples.tsv \
        dtu_results/drimseq_fit.rds \
        dtu_results/R_session_info.txt
    do

        if [[ ! -s "$output" ]]; then
            echo "ERROR: Missing output: $output" >&2
            exit 1
        fi

    done

    echo "========================================"
    echo "Krill DTU completed successfully"
    echo "Mode: $KRILL_DTU_MODE"
    echo "Design: $KRILL_DTU_DESIGN"
    echo "========================================"

    '''
}
