include { BUILD_MINIMAP2_INDEX }  from '../modules/local/build_minimap2_index/main'
include { PREPROCESS_ANNOTATION } from '../modules/local/preprocess_annotation/main'
include { MINIMAP2_ALIGN }        from '../modules/local/minimap2_align/main'
include { SAMTOOLS_SORT }         from '../modules/local/samtools_sort/main'
include { SAMTOOLS_INDEX }        from '../modules/local/samtools_index/main'

workflow GENOME_ALIGNMENT {

    take:
    merged_fastq_ch
    reference_fasta_ch
    annotation_gtf_ch

    main:

    BUILD_MINIMAP2_INDEX(
        reference_fasta_ch
    )

    PREPROCESS_ANNOTATION(
        annotation_gtf_ch,
        reference_fasta_ch
    )

    MINIMAP2_ALIGN(
        merged_fastq_ch,
        BUILD_MINIMAP2_INDEX.out
    )

    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.sam
    )

    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )

    emit:
    indexed_bam = SAMTOOLS_INDEX.out.indexed_bam
    cleaned_gtf = PREPROCESS_ANNOTATION.out.cleaned_gtf
    genome_index = BUILD_MINIMAP2_INDEX.out
}
