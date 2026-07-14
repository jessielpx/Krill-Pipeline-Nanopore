process PREPROCESS_ANNOTATION {

    tag "annotation"

    container '/project/def-lefranco/peixiliu/NXF_SINGULARITY_CACHEDIR/ontresearch-wf-transcriptomes-shaaaf20a5a0e76f9e18bad21af639a6b69e4a31a2f.sif'

    cpus 2
    memory '16 GB'
    time '2h'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path annotation_gtf
    path reference_fasta

    output:
    path "annotation.cleaned.gtf", emit: cleaned_gtf

    script:
    """
    gffread \
        ${annotation_gtf} \
        -g ${reference_fasta} \
        -T \
        -o annotation.cleaned.gtf
    """
}
