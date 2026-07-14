process SAMTOOLS_SORT {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 4
    memory '16 GB'
    time '6h'

    publishDir "${params.outdir}/bam/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.sample}.sorted.bam"), emit: bam

    script:
    """
    samtools sort \
        -@ ${task.cpus} \
        -m 3G \
        -o ${meta.sample}.sorted.bam \
        ${sam}
    """
}
