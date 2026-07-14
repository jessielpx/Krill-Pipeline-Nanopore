process BUILD_TRANSCRIPTOME_INDEX {

    tag "transcriptome"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 4
    memory '32 GB'
    time '4h'

    publishDir "${params.outdir}/transcripts/merged", mode: 'copy'

    input:
    path transcriptome_fasta

    output:
    path "transcriptome.mmi", emit: index

    script:
    """
    minimap2 \
        -t ${task.cpus} \
        -I 1000G \
        -d transcriptome.mmi \
        ${transcriptome_fasta}
    """
}
