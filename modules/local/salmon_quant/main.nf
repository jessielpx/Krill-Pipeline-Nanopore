process SALMON_QUANT {

    tag "${meta.sample}"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 8
    memory '32 GB'
    time '8h'

    publishDir "${params.outdir}/quantification/${meta.sample}", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path transcriptome_fasta

    output:
    tuple val(meta),
          path("${meta.sample}.transcript_counts.tsv"),
          emit: counts

    script:
    """
    salmon quant \
        --noErrorModel \
        -p ${task.cpus} \
        -t ${transcriptome_fasta} \
        -l SF \
        -a ${bam} \
        -o counts

    mv counts/quant.sf \
       ${meta.sample}.transcript_counts.tsv
    """
}
