include { GFFCOMPARE }       from '../modules/local/gffcompare/main'
include { PARSE_GFFCOMPARE } from '../modules/local/parse_gffcompare/main'

workflow ANNOTATION_ANALYSIS {

    take:
    merged_gtf_ch
    reference_gtf_ch
    parse_script_ch

    main:

    /*
     * Compare reconstructed transcripts against reference
     */
    GFFCOMPARE(
        merged_gtf_ch,
        reference_gtf_ch
    )

    /*
     * Create a transcript-level annotation table
     */
    PARSE_GFFCOMPARE(
        GFFCOMPARE.out.annotated_gtf,
        parse_script_ch
    )

    emit:
    annotated_gtf = GFFCOMPARE.out.annotated_gtf
    stats = GFFCOMPARE.out.stats
    tracking = GFFCOMPARE.out.tracking
    loci = GFFCOMPARE.out.loci

    annotation_table =
        PARSE_GFFCOMPARE.out.annotation_table
}
