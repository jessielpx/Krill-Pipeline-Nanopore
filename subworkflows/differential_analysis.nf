include { FILTER_UNSTRANDED_ANNOTATION } from '../modules/local/filter_unstranded_annotation/main'
include { DIFFERENTIAL_EXPRESSION }      from '../modules/local/differential_expression/main'

workflow DIFFERENTIAL_ANALYSIS {

    take:
    sample_sheet_ch
    transcript_counts_ch
    merged_gtf_ch
    de_script_ch

    main:

    /*
     * Remove transcripts without a valid strand
     */
    FILTER_UNSTRANDED_ANNOTATION(
        merged_gtf_ch
    )

    /*
     * Run filtering, DGE and DTU
     */
    DIFFERENTIAL_EXPRESSION(
        sample_sheet_ch,
        transcript_counts_ch,
        FILTER_UNSTRANDED_ANNOTATION.out.stranded_gtf,
        de_script_ch
    )

    emit:
    stranded_gtf = FILTER_UNSTRANDED_ANNOTATION.out.stranded_gtf
    excluded_gtf = FILTER_UNSTRANDED_ANNOTATION.out.excluded_gtf

    dge = DIFFERENTIAL_EXPRESSION.out.dge
    dge_pdf = DIFFERENTIAL_EXPRESSION.out.dge_pdf
    dexseq = DIFFERENTIAL_EXPRESSION.out.dexseq

    dtu_gene = DIFFERENTIAL_EXPRESSION.out.dtu_gene
    dtu_transcript = DIFFERENTIAL_EXPRESSION.out.dtu_transcript
    dtu_stageR = DIFFERENTIAL_EXPRESSION.out.dtu_stageR
    dtu_pdf = DIFFERENTIAL_EXPRESSION.out.dtu_pdf

    cpm = DIFFERENTIAL_EXPRESSION.out.cpm
    filtered_counts = DIFFERENTIAL_EXPRESSION.out.filtered_counts
    gene_counts = DIFFERENTIAL_EXPRESSION.out.gene_counts
}
