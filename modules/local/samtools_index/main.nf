process SAMTOOLS_INDEX {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 2
    memory '16 GB'
    time '2h'

    publishDir "${params.outdir}/bam/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
          path(bam),
          path("${bam}.bai"),
          emit: indexed_bam

    script:
    """
    samtools index \
        -@ ${task.cpus} \
        ${bam}
    """
}
