/*
 * Krill Renew 11 - Junction extraction
 *
 * Inputs:
 *   List of sorted BAM files from CORE
 *   Cohort samplesheet
 *   Reference GTF
 *   Python extraction script
 *
 * Only uses BAM files from the current workflow.
 * No reliance on previously published output files.
 */

process JUNCTION_EXTRACT {

    tag "junction_quantification"

    cpus 4
    memory '32 GB'
    time '24h'

    publishDir "${params.outdir}/11_alternative_splicing/junctions",
        mode: 'copy'

    input:
    path bam_files
    path samplesheet
    path reference_gtf
    path junction_script

    output:
    path 'junction_results', emit: results

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load python/3.11.5
    module load samtools/1.20

    unset PYTHONHOME PYTHONPATH
    hash -r

    mkdir -p junction_results

    # Construct the manifest using staged BAM files.
    # Confirm an exact match to each samplesheet ID.
    python3 - "!{samplesheet}" <<'PY'
import csv
import glob
import os
import sys
from pathlib import Path

sheet = sys.argv[1]

with open(sheet, newline="") as handle:
    first = handle.readline()
    handle.seek(0)
    delim = "\\t" if "\\t" in first else ","
    entries = list(csv.DictReader(handle, delimiter=delim))

if not entries or "sample" not in entries[0]:
    raise SystemExit("Invalid samplesheet")

samples = [x["sample"].strip() for x in entries]

if len(samples) != len(set(samples)):
    raise SystemExit("Duplicate sample IDs")

if any(not s for s in samples):
    raise SystemExit("Blank sample ID")

bam_paths = glob.glob("*.sorted.bam")

bam_map = {}

for bam in bam_paths:
    name = Path(bam).name

    if not name.endswith(".sorted.bam"):
        continue

    sample = name[:-len(".sorted.bam")]

    if sample in bam_map:
        raise SystemExit(
            "Duplicate BAM for sample: " + sample
        )

    bam_map[sample] = bam

if set(bam_map) != set(samples):
    raise SystemExit(
        "BAM/sample mismatch. Missing: "
        + str(sorted(set(samples) - set(bam_map)))
        + "; unexpected: "
        + str(sorted(set(bam_map) - set(samples)))
    )

with open("bam_manifest.tsv", "w") as out:
    out.write("sample\\tbam\\n")

    for sample in samples:
        out.write(
            sample + "\\t"
            + os.path.abspath(bam_map[sample])
            + "\\n"
        )

print("BAM manifest:", len(samples), "samples")
PY

    python3 "!{junction_script}" \
        --manifest bam_manifest.tsv \
        --reference-gtf "!{reference_gtf}" \
        --outdir junction_results \
        --min-mapq !{params.junction_min_mapq != null ? params.junction_min_mapq : 10} \
        --min-intron !{params.junction_min_intron != null ? params.junction_min_intron : 20}

    test -s junction_results/junction_counts.tsv
    test -s junction_results/junction_annotation.tsv
    test -s junction_results/junction_sample_qc.tsv
    test -s junction_results/junction_qc_summary.tsv

    cat junction_results/junction_qc_summary.tsv

    echo "JUNCTION_EXTRACT completed."
    '''
}


/*
 * Krill Renew 11 - Differential junction usage
 *
 * Exploratory gene-level and junction-level
 * DRIMSeq analysis.
 */

process JUNCTION_DJU {

    tag "junction_differential"

    cpus 2
    memory '32 GB'
    time '12h'

    publishDir "${params.outdir}/11_alternative_splicing/dju",
        mode: 'copy'

    input:
    path junction_results
    path samplesheet
    path dju_script

    output:
    path 'dju_results', emit: results

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load r/4.6.1

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


    export KRILL_DJU_DESIGN='!{params.dju_design ?: "~ group"}'

    export KRILL_DJU_MIN_COUNT='!{params.dju_min_count != null ? params.dju_min_count : 5}'

    export KRILL_DJU_GROUP_FRACTION='!{params.dju_group_fraction != null ? params.dju_group_fraction : 0.25}'

    COUNTS="!{junction_results}/junction_counts.tsv"

    test -s "$COUNTS"

    Rscript "!{dju_script}" \
        "$COUNTS" \
        "!{samplesheet}" \
        dju_results

    test -s dju_results/dju_gene_results.tsv
    test -s dju_results/dju_junction_results.tsv
    test -s dju_results/dju_input_qc.tsv

    cat dju_results/dju_input_qc.tsv

    echo "JUNCTION_DJU completed."
    '''
}
