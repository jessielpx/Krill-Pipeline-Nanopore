/*
 * Krill Renew 10 - GSEA
 *
 * DESeq2 Wald-statistic preranked fgsea.
 * Cached human MSigDB Ensembl GMT files.
 *
 * Independent of DTU, isoform usage and
 * functional consequences.
 */

process GSEA {

    tag "preranked"

    cpus 2
    memory '16 GB'
    time '8h'

    publishDir "${params.outdir}/10_gsea",
        mode: 'copy'

    input:
    path dge_results
    path gmt_dir
    path gsea_script

    output:
    path 'gsea_results', emit: results

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


    DGE="!{dge_results}/dge_HR_vs_control.tsv"
    GMT="!{gmt_dir}"
    SCRIPT="!{gsea_script}"

    test -s "$DGE"
    test -s "$GMT/hallmark.gmt"
    test -s "$GMT/reactome.gmt"
    test -s "$GMT/go_bp.gmt"
    test -s "$SCRIPT"

    Rscript "$SCRIPT" \
        --mode run \
        --dge "$DGE" \
        --gmt-dir "$GMT" \
        --outdir gsea_results \
        --min-size !{params.gsea_min_size != null ? params.gsea_min_size : 10} \
        --max-size !{params.gsea_max_size != null ? params.gsea_max_size : 500}

    test -s gsea_results/gsea_ranked_genes.tsv
    test -s gsea_results/gsea_all_results.tsv
    test -s gsea_results/gsea_qc.tsv

    cat gsea_results/gsea_qc.tsv

    echo "Krill GSEA completed."
    '''
}
