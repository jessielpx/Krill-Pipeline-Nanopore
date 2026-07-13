process FASTCAT {

    tag "${barcode}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-common-sha72f3517dd994984e0e2da0b97cb3f23f8540be4b.sif'

    cpus 2
    memory '8 GB'
    time '6h'

    input:
    tuple val(barcode), path(read_dir)

    output:
    tuple val(barcode), path("${barcode}.fastq.gz"), emit: fastq
    tuple val(barcode), path("${barcode}.file_summary.tsv"), emit: summary

    script:
    """
    fastcat \
        --file ${barcode}.file_summary.tsv \
        ${read_dir} \
        | gzip -c > ${barcode}.fastq.gz
    """
}
