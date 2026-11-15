/*
 * Krill Renew - 09 FUNCTIONAL_CONSEQUENCES
 *
 * Inputs:
 *   ISOFORM_USAGE results
 *   ISOFORM_STRUCTURE_ALL results
 *   DTU results
 *   Samplesheet
 *   Reference GTF
 *
 * Output:
 *   all_isoform_pairs.tsv
 *   pair_qc_pass_pairs.tsv
 *   functional_qc_summary.tsv
 *
 * Native execution on Rorqual.
 */

process FUNCTIONAL_CONSEQUENCES {

    tag "cohort"

    cpus 2
    memory '8 GB'
    time '4h'

    publishDir "${params.outdir}/09_functional_consequences",
        mode: 'copy'

    input:
    path usage_results
    path structure_results
    path dtu_results
    path samplesheet
    path reference_gtf
    path analysis_script

    output:
    path 'functional_results', emit: results

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load python/3.11.5

    unset PYTHONHOME PYTHONPATH
    hash -r

    USAGE="!{usage_results}/isoform_usage_supported.tsv"
    STRUCTURE="!{structure_results}/isoform_structure_annotations.tsv"
    DTU="!{dtu_results}/dtu_gene_results.tsv"
    SAMPLESHEET="!{samplesheet}"
    REFERENCE="!{reference_gtf}"
    SCRIPT="!{analysis_script}"

    for required in \
        "$USAGE" \
        "$STRUCTURE" \
        "$DTU" \
        "$SAMPLESHEET" \
        "$REFERENCE" \
        "$SCRIPT"
    do
        test -s "$required" || {
            echo "ERROR: Missing input: $required" >&2
            exit 1
        }
    done

    echo "===== KRILL FUNCTIONAL CONSEQUENCES ====="

    python3 "$SCRIPT" \
        --samplesheet "$SAMPLESHEET" \
        --usage "$USAGE" \
        --structure "$STRUCTURE" \
        --dtu "$DTU" \
        --reference-gtf "$REFERENCE" \
        --outdir functional_results \
        --min-reads !{params.functional_min_reads ?: 3} \
        --min-supported-samples !{params.functional_min_supported_samples ?: 2} \
        --group-fraction !{params.functional_group_fraction ?: 0.25}

    test -s functional_results/all_isoform_pairs.tsv

    test -s functional_results/pair_qc_pass_pairs.tsv

    test -s functional_results/functional_qc_summary.tsv

    echo "===== FUNCTIONAL QC ====="

    cat functional_results/functional_qc_summary.tsv

    echo "Krill FUNCTIONAL_CONSEQUENCES completed."
    '''
}
