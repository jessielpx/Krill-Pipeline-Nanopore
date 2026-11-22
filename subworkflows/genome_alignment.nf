
/*
 * Krill Renew - Genome Alignment Subworkflow
 *
 * 1. Build minimap2 genome index
 * 2. Preprocess reference annotation
 * 3. Align cDNA reads to genome
 * 4. Index BAM and generate QC statistics
 */

include { BUILD_MINIMAP2_INDEX }  from '../modules/local/build_minimap2_index/main'
include { PREPROCESS_ANNOTATION } from '../modules/local/preprocess_annotation/main'
include { MINIMAP2_ALIGN }        from '../modules/local/minimap2_align/main'
include { SAMTOOLS_INDEX }        from '../modules/local/samtools_index/main'


workflow GENOME_ALIGNMENT {

    take:
    merged_fastq_ch
    reference_fasta_ch
    annotation_gtf_ch

    main:

    /*
     * Step 1: Build genome index
     */
    BUILD_MINIMAP2_INDEX(
        reference_fasta_ch
    )

    /*
     * Make reference index reusable across samples
     */
    genome_index_ch = BUILD_MINIMAP2_INDEX.out

    /*
     * Step 2: Prepare annotation
     */
    PREPROCESS_ANNOTATION(
        annotation_gtf_ch,
        reference_fasta_ch
    )

    /*
     * Step 3: Genome alignment
     */
    MINIMAP2_ALIGN(
        merged_fastq_ch,
        genome_index_ch,
        PREPROCESS_ANNOTATION.out.cleaned_gtf
    )

    /*
     * Step 4: BAM indexing and QC
     */
    SAMTOOLS_INDEX(
        MINIMAP2_ALIGN.out.bam
    )

    emit:
    indexed_bam = SAMTOOLS_INDEX.out.indexed_bam
    cleaned_gtf = PREPROCESS_ANNOTATION.out.cleaned_gtf
    genome_index = genome_index_ch
    alignment_qc = SAMTOOLS_INDEX.out.qc
    alignment_logs = MINIMAP2_ALIGN.out.log
    annotation_qc = PREPROCESS_ANNOTATION.out.qc
}
