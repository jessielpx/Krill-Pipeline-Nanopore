/*
 * Krill Renew 12 - Cohort fusion aggregation
 *
 * One CSV per sample.
 * Does not mix samples from prior analyses.
 */

process FUSION_MERGE {

    tag "cohort"

    cpus 1
    memory '8 GB'
    time '2h'

    publishDir "${params.outdir}/12_fusion_detection/cohort",
        mode: 'copy'

    input:
    path fusion_csvs
    path samplesheet
    path merge_script

    output:
    path 'fusion_results', emit: results

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load python/3.11.5

    unset PYTHONHOME PYTHONPATH
    hash -r

    # JAFFAL outputs have unique names:
    # sample14.jaffa_results.csv, etc.
    python3 - <<'PY'
import csv
from pathlib import Path

files = sorted(
    Path(".").glob("*.jaffa_results.csv")
)

if not files:
    raise SystemExit(
        "ERROR: No JAFFAL CSV files staged"
    )

suffix = ".jaffa_results.csv"

with open("fusion_manifest.tsv", "w", newline="") as out:
    writer = csv.writer(
        out,
        delimiter=chr(9)
    )

    writer.writerow(["sample", "csv"])

    for path in files:
        sample = path.name[:-len(suffix)]

        if not sample:
            raise SystemExit(
                "ERROR: Invalid JAFFAL filename"
            )

        writer.writerow([
            sample,
            str(path.resolve())
        ])

print("JAFFAL inputs:", len(files))
PY

    python3 "!{merge_script}" \
        --manifest fusion_manifest.tsv \
        --samplesheet "!{samplesheet}" \
        --outdir fusion_results \
        --min-spanning-reads !{params.fusion_min_spanning_reads != null ? params.fusion_min_spanning_reads : 2} \
        --min-recurrent-samples !{params.fusion_min_recurrent_samples != null ? params.fusion_min_recurrent_samples : 2}

    test -s fusion_results/all_fusions.tsv
    test -s fusion_results/filtered_fusions.tsv
    test -s fusion_results/recurrent_fusions.tsv
    test -s fusion_results/fusion_qc.tsv

    cat fusion_results/fusion_qc.tsv

    echo "FUSION_MERGE completed."
    '''
}
