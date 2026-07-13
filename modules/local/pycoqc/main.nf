process PYCOQC {

    tag "${meta.id}"

    container 'quay.io/biocontainers/pycoqc:2.5.2--py_0'

    cpus 2
    memory '64 GB'
    time '2h'

    publishDir "${params.outdir}/pycoqc", mode: 'copy'

    input:
    tuple val(meta), path(summary_file)

    output:
    tuple val(meta), path("${meta.id}.pycoqc.html"), emit: html
    tuple val(meta), path("${meta.id}.pycoqc.json"), emit: json

    script:
    """
    pycoQC \
        --summary_file ${summary_file} \
        --html_outfile ${meta.id}.pycoqc.html \
        --json_outfile ${meta.id}.pycoqc.json \
        --report_title "${meta.id} pycoQC report"
    """
}
