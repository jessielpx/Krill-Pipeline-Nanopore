process MSIGDB_PREPARE {

    tag "MSigDB"

    cpus 2
    memory '8 GB'
    time '4h'

    publishDir "${params.outdir}/00_resources",
        mode: 'copy'

    input:
    path gsea_script

    output:
    path "msigdb_human_ensembl", emit: gmt_dir

    shell:
    '''
    set -euo pipefail

    module --force purge
    module load StdEnv/2023
    module load r/4.6.1

    unset PYTHONPATH PYTHONHOME

    export R_LIBS_USER="!{params.r_libs_user}"

    Rscript --vanilla -e '
        stopifnot(requireNamespace("msigdbr", quietly=TRUE))
        cat("msigdbr:", as.character(packageVersion("msigdbr")), "\\n")
    '

    Rscript --vanilla "!{gsea_script}" \
        --mode prepare \
        --gmt-dir msigdb_human_ensembl

    for f in hallmark.gmt reactome.gmt go_bp.gmt; do
        test -s "msigdb_human_ensembl/$f"
    done
    '''
}
