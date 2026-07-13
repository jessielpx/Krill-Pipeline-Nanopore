process MINIMAP2_ALIGN {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 8
    memory '32 GB'
    time '12h'

    publishDir "${params.outdir}/alignment/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path minimap_index

    output:
    tuple val(meta), path("${meta.sample}.sam"), emit: sam

    script:
    """
    minimap2 \
        -ax splice \
        -uf \
        -k14 \
        -t ${task.cpus} \
        ${minimap_index} \
        ${reads} \
        > ${meta.sample}.sam
    """
}
