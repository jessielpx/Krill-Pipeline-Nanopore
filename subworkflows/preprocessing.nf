
/*
 * Krill Renew - Stage 1: Preprocessing
 *
 * FASTCAT -> Genome Alignment -> BAM QC
 */

include { FASTCAT } from '../modules/local/fastcat/main'
include { GENOME_ALIGNMENT } from './genome_alignment'


workflow PREPROCESSING {

    take:
    samples_ch
    reference_fasta_ch
    annotation_gtf_ch

    main:

    /*
     * 1. Merge FASTQ files for each sample
     */
    FASTCAT(samples_ch)

    /*
     * 2. Extract merged FASTQ
     */
    merged_fastq_ch = FASTCAT.out.results.map {
        meta, fastq_file, stats_dir ->

        tuple(meta, fastq_file)
    }

    /*
     * 3. Genome alignment and BAM processing
     */
    GENOME_ALIGNMENT(
        merged_fastq_ch,
        reference_fasta_ch,
        annotation_gtf_ch
    )

    /*
     * 4. Expose outputs
     */
    emit:
    merged_fastq  = merged_fastq_ch
    fastcat_stats = FASTCAT.out.results
    indexed_bam   = GENOME_ALIGNMENT.out.indexed_bam
    cleaned_gtf   = GENOME_ALIGNMENT.out.cleaned_gtf
    genome_index  = GENOME_ALIGNMENT.out.genome_index
    alignment_qc  = GENOME_ALIGNMENT.out.alignment_qc
    annotation_qc = GENOME_ALIGNMENT.out.annotation_qc
}
