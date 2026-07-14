process MAP_TRANSCRIPTOME {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 8
    memory '32 GB'
    time '8h'

    publishDir "${params.outdir}/transcriptome_alignment/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path transcriptome_index

    output:
    tuple val(meta),
          path("${meta.sample}.transcriptome.sorted.bam"),
          emit: bam

    script:
    """
    set -euo pipefail

    minimap2 \
        -t ${task.cpus} \
        -ax map-ont \
        -p 1.0 \
        ${transcriptome_index} \
        ${reads} \
    | samtools view -Sb - \
    | samtools sort \
        -@ ${task.cpus} \
        -o ${meta.sample}.transcriptome.sorted.bam
    """
}
