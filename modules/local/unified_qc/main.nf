/*
 * Krill Renew - Unified Core QC
 *
 * Combine:
 *   FASTQ read statistics
 *   samtools alignment QC
 *   bambu quantification QC
 *   SQANTI3 structural classification
 *
 * Native execution on Rorqual.
 * No containers.
 */

process UNIFIED_QC {

    tag "cohort"

    cpus 1
    memory '4 GB'
    time '1h'

    publishDir "${params.outdir}/04_unified_qc",
        mode: 'copy'

    input:
    // Each FASTCAT output directory has the same name.
    // Stage in separate subdirectories to avoid collisions.
    path fastcat_stats_dirs, stageAs: 'fastcat/sample??/*'

    // Flagstat filenames contain unique sample IDs.
    path alignment_flagstats, stageAs: 'alignment/*'

    path bambu_results
    path sqanti3_results
    path qc_script

    output:
    tuple path('sample_qc_summary.tsv'),
          path('sqanti3_category_summary.tsv'),
          path('qc_summary.txt'),
          emit: results

    shell:
    '''
    set -euo pipefail

    # Use the previously verified Python environment.
    PYTHON="!{params.sqanti_env ?: ''}/bin/python"

    if [[ "$PYTHON" == "/bin/python" ]]; then
        echo "ERROR: Configure sqanti_env in params.yaml" >&2
        exit 1
    fi

    unset PYTHONHOME PYTHONPATH PYTHONUSERBASE

    test -x "$PYTHON"
    test -s "!{qc_script}"
    test -s "!{bambu_results}/bambu_sample_qc.tsv"
    test -s "!{sqanti3_results}/bambu_sqanti3_classification.txt"

    echo "===== UNIFIED CORE QC ====="

    "$PYTHON" "!{qc_script}" \
        --fastcat fastcat \
        --alignment alignment \
        --bambu "!{bambu_results}" \
        --sqanti "!{sqanti3_results}" \
        --out .

    test -s sample_qc_summary.tsv
    test -s sqanti3_category_summary.tsv
    test -s qc_summary.txt

    echo "Unified Core QC completed."
    '''
}
