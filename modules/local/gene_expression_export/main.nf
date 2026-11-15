/*
 * Krill Renew - Gene Expression Export
 *
 * Bambu gene counts
 * Bambu transcript-derived gene TPM
 * Ensembl gene ID -> gene symbol
 *
 * Independent of SQANTI3, DGE and DTU.
 */

process GENE_EXPRESSION_EXPORT {

    tag "cohort"

    cpus 2
    memory '16 GB'
    time '4h'

    publishDir "${params.outdir}/02_gene_expression",
        mode: 'copy'

    input:
    path bambu_results
    path reference_gtf
    path export_script

    output:
    path "gene_expression_results", emit: results

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load python/3.11.5

    unset PYTHONHOME PYTHONPATH PYTHONUSERBASE
    export PYTHONNOUSERSITE=1

    echo "===== GENE EXPRESSION EXPORT ====="

    python3 --version

    test -s "!{bambu_results}/counts_gene.txt"
    test -s "!{bambu_results}/counts_transcript.txt"
    test -s "!{bambu_results}/extended_annotations.gtf"
    test -s "!{reference_gtf}"
    test -s "!{export_script}"

    python3 "!{export_script}" \
        "!{bambu_results}" \
        "!{reference_gtf}" \
        gene_expression_results

    echo "===== VALIDATE EXPORT ====="

    for filename in \
        gene_counts_ENSG.tsv \
        gene_counts_symbols.tsv \
        gene_TPM_ENSG.tsv \
        gene_TPM_symbols.tsv \
        gene_id_mapping.tsv \
        export_qc_summary.tsv
    do
        test -s "gene_expression_results/$filename"
    done

    echo "GENE_EXPRESSION_EXPORT completed."
    '''
}
