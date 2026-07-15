nextflow.enable.dsl = 2

include { QC }                     from './subworkflows/qc'
include { GENOME_ALIGNMENT }       from './subworkflows/genome_alignment'
include { TRANSCRIPTOME_ANALYSIS } from './subworkflows/transcriptome_analysis'
include { DIFFERENTIAL_ANALYSIS }  from './subworkflows/differential_analysis'
include { ANNOTATION_ANALYSIS }    from './subworkflows/annotation_analysis'


workflow {

    /*
     * Check required parameters
     */
    if (!params.summary) {
        error "Please provide --summary /path/to/sequencing_summary.txt"
    }

    if (!params.fastq_dir) {
        error "Please provide --fastq_dir /path/to/barcode_directories"
    }

    if (!params.ref_genome) {
        error "Please provide --ref_genome /path/to/reference.fa"
    }

    if (!params.ref_annotation) {
        error "Please provide --ref_annotation /path/to/annotation.gtf"
    }


    /*
     * Input files
     */
    sample_sheet_file = file(
        "${projectDir}/assets/samples.csv",
        checkIfExists: true
    )

    summary_ch = Channel.value(
        tuple(
            [id: params.run_id],
            file(
                params.summary,
                checkIfExists: true
            )
        )
    )

    reference_fasta_ch = Channel.value(
        file(
            params.ref_genome,
            checkIfExists: true
        )
    )

    annotation_gtf_ch = Channel.value(
        file(
            params.ref_annotation,
            checkIfExists: true
        )
    )

    sample_sheet_ch = Channel.value(
        sample_sheet_file
    )


    /*
     * Read sample metadata and locate barcode directories
     */
    samples_ch = Channel
        .fromPath(
            sample_sheet_file,
            checkIfExists: true
        )
        .splitCsv(header: true)
        .map { row ->

            def meta = [
                barcode  : row.barcode,
                sample   : row.sample,
                alias    : row.alias,
                condition: row.condition
            ]

            def barcode_dir = file(
                "${params.fastq_dir}/${row.barcode}",
                checkIfExists: true
            )

            tuple(meta, barcode_dir)
        }


    /*
     * Local helper scripts
     */
    merge_script_ch = Channel.value(
        file(
            "${projectDir}/bin/merge_salmon_counts.py",
            checkIfExists: true
        )
    )

    annotation_script_ch = Channel.value(
        file(
            "${projectDir}/bin/add_gene_ids_and_sum_counts.py",
            checkIfExists: true
        )
    )

    parse_gffcompare_script_ch = Channel.value(
        file(
            "${projectDir}/bin/parse_gffcompare_annotation.py",
            checkIfExists: true
        )
    )

    de_script_ch = Channel.value(
        file(
            "${projectDir}/bin/de_analysis.R",
            checkIfExists: true
        )
    )


    /*
     * Run-level QC and FASTQ preparation
     */
    QC(
        summary_ch,
        samples_ch
    )


    /*
     * Extract merged FASTQ files
     */
    merged_fastq_ch = QC.out.fastcat_results.map {
        meta,
        fastq_files,
        stats_dir ->

        tuple(meta, fastq_files)
    }


    /*
     * Genome alignment and BAM processing
     */
    GENOME_ALIGNMENT(
        merged_fastq_ch,
        reference_fasta_ch,
        annotation_gtf_ch
    )


    /*
     * Transcript reconstruction and quantification
     */
    TRANSCRIPTOME_ANALYSIS(
        GENOME_ALIGNMENT.out.indexed_bam,
        merged_fastq_ch,
        GENOME_ALIGNMENT.out.cleaned_gtf,
        reference_fasta_ch,
        merge_script_ch,
        annotation_script_ch
    )


    /*
     * gffcompare transcript annotation
     */
    ANNOTATION_ANALYSIS(
        TRANSCRIPTOME_ANALYSIS.out.merged_gtf,
        GENOME_ALIGNMENT.out.cleaned_gtf,
        parse_gffcompare_script_ch
    )


    /*
     * Differential expression and transcript usage analysis
     */
    if (params.de_analysis) {

        DIFFERENTIAL_ANALYSIS(
            sample_sheet_ch,
            TRANSCRIPTOME_ANALYSIS.out.transcript_counts,
            TRANSCRIPTOME_ANALYSIS.out.merged_gtf,
            de_script_ch
        )
    }
}
