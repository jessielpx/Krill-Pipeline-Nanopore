process STRINGTIE_ASSEMBLY {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 4
    memory '32 GB'
    time '8h'

    publishDir "${params.outdir}/transcripts/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path annotation_gtf

    output:
    tuple val(meta),
          path("${meta.sample}.assembled.gtf"),
          emit: assembled_gtf

    script:
    """
    stringtie \
        ${bam} \
        -p ${task.cpus} \
        -G ${annotation_gtf} \
        -o ${meta.sample}.assembled.gtf
    """
}
