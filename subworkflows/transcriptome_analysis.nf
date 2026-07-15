include { STRINGTIE_ASSEMBLY }        from '../modules/local/stringtie_assembly/main'
include { STRINGTIE_MERGE }           from '../modules/local/stringtie_merge/main'
include { BUILD_TRANSCRIPTOME_INDEX } from '../modules/local/build_transcriptome_index/main'
include { MAP_TRANSCRIPTOME }         from '../modules/local/map_transcriptome/main'
include { SALMON_QUANT }              from '../modules/local/salmon_quant/main'
include { MERGE_TRANSCRIPT_COUNTS }   from '../modules/local/merge_transcript_counts/main'
include { ANNOTATE_AND_SUM_COUNTS }   from '../modules/local/annotate_and_sum_counts/main'

workflow TRANSCRIPTOME_ANALYSIS {

    take:
    indexed_bam_ch
    merged_fastq_ch
    cleaned_gtf_ch
    reference_fasta_ch
    merge_script_ch
    annotation_script_ch

    main:

    /*
     * Assemble transcripts for each sample
     */
    STRINGTIE_ASSEMBLY(
        indexed_bam_ch,
        cleaned_gtf_ch
    )

    /*
     * Collect all sample-level assembled GTF files
     */
    assembled_gtfs_ch = STRINGTIE_ASSEMBLY.out.assembled_gtf
        .map { meta, assembled_gtf ->
            assembled_gtf
        }
        .collect()

    /*
     * Build merged transcriptome
     */
    STRINGTIE_MERGE(
        assembled_gtfs_ch,
        cleaned_gtf_ch,
        reference_fasta_ch
    )

    /*
     * Build transcriptome Minimap2 index
     */
    BUILD_TRANSCRIPTOME_INDEX(
        STRINGTIE_MERGE.out.transcriptome_fasta
    )

    /*
     * Align reads to merged transcriptome
     */
    MAP_TRANSCRIPTOME(
        merged_fastq_ch,
        BUILD_TRANSCRIPTOME_INDEX.out.index
    )

    /*
     * Quantify transcript abundance
     */
    SALMON_QUANT(
        MAP_TRANSCRIPTOME.out.bam,
        STRINGTIE_MERGE.out.transcriptome_fasta
    )

    /*
     * Collect Salmon outputs
     */
    transcript_count_files_ch = SALMON_QUANT.out.counts
        .map { meta, count_file ->
            count_file
        }
        .collect()

    /*
     * Merge transcript counts and TPM
     */
    MERGE_TRANSCRIPT_COUNTS(
        transcript_count_files_ch,
        merge_script_ch
    )

    /*
     * Add gene IDs and construct gene-count matrix
     */
    ANNOTATE_AND_SUM_COUNTS(
        MERGE_TRANSCRIPT_COUNTS.out.counts,
        STRINGTIE_MERGE.out.merged_gtf,
        annotation_script_ch
    )

    emit:
    merged_gtf = STRINGTIE_MERGE.out.merged_gtf
    transcriptome_fasta = STRINGTIE_MERGE.out.transcriptome_fasta

    transcript_counts = MERGE_TRANSCRIPT_COUNTS.out.counts
    transcript_tpm = MERGE_TRANSCRIPT_COUNTS.out.tpm

    transcript_counts_with_genes =
        ANNOTATE_AND_SUM_COUNTS.out.transcript_counts_with_genes

    gene_counts =
        ANNOTATE_AND_SUM_COUNTS.out.gene_counts

    transcriptome_bam = MAP_TRANSCRIPTOME.out.bam
    salmon_counts = SALMON_QUANT.out.counts
}
