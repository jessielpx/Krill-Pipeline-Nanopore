nextflow.enable.dsl = 2

include { PYCOQC }               from './modules/local/pycoqc/main'
include { FASTCAT }              from './modules/local/fastcat/main'
include { BUILD_MINIMAP2_INDEX } from './modules/local/build_minimap2_index/main'
include { PREPROCESS_ANNOTATION } from './modules/local/preprocess_annotation/main'
include { MINIMAP2_ALIGN }       from './modules/local/minimap2_align/main'
include { SAMTOOLS_SORT }        from './modules/local/samtools_sort/main'
include { SAMTOOLS_INDEX }       from './modules/local/samtools_index/main'

workflow {

    /*
     * Check required inputs
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
     * Run-level quality control
     */
    summary_ch = Channel.value(
        tuple(
            [id: params.run_id],
            file(params.summary, checkIfExists: true)
        )
    )

    PYCOQC(summary_ch)

    /*
     * Read sample metadata and locate barcode directories
     */
    samples_ch = Channel
        .fromPath(
            "${projectDir}/assets/samples.csv",
            checkIfExists: true
        )
        .splitCsv(header: true)
        .map { row ->

            def meta = [
                barcode  : row.barcode,
                sample   : row.sample,
                condition: row.condition
            ]

            def barcode_dir = file(
                "${params.fastq_dir}/${row.barcode}",
                checkIfExists: true
            )

            tuple(meta, barcode_dir)
        }

    /*
     * Merge FASTQ files and generate Fastcat statistics
     */
    FASTCAT(samples_ch)

    /*
     * Reusable reference inputs
     */
    reference_fasta_ch = Channel.value(
        file(params.ref_genome, checkIfExists: true)
    )

    annotation_gtf_ch = Channel.value(
        file(params.ref_annotation, checkIfExists: true)
    )

    /*
     * Build Minimap2 reference index
     */
    BUILD_MINIMAP2_INDEX(reference_fasta_ch)

    /*
     * Clean and validate the reference annotation
     */
    PREPROCESS_ANNOTATION(
        annotation_gtf_ch,
        reference_fasta_ch
    )

    /*
     * Extract merged FASTQ files from FASTCAT output
     */
    merged_fastq_ch = FASTCAT.out.results.map {
        meta,
        fastq_files,
        stats_dir ->

        tuple(meta, fastq_files)
    }

    /*
     * Reuse the same Minimap2 index for every sample
     */
    minimap_index_ch = BUILD_MINIMAP2_INDEX.out.first()

    /*
     * Align merged FASTQ files to the reference
     */
    MINIMAP2_ALIGN(
        merged_fastq_ch,
        minimap_index_ch
    )

    /*
     * Sort SAM files into coordinate-sorted BAM files
     */
    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.sam
    )

    /*
     * Index sorted BAM files
     */
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
}
