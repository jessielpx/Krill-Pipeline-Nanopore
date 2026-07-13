nextflow.enable.dsl = 2

include { PYCOQC }               from './modules/local/pycoqc/main'
include { FASTCAT }              from './modules/local/fastcat/main'
include { BUILD_MINIMAP2_INDEX } from './modules/local/build_minimap2_index/main'
include { MINIMAP2_ALIGN }       from './modules/local/minimap2_align/main'

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

    /*
     * Run-level quality control
     */
    summary_ch = Channel.of(
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
     * Merge FASTQ files and generate fastcat statistics
     */
    FASTCAT(samples_ch)

    /*
     * Build Minimap2 reference index
     */
    reference_ch = Channel.of(
        file(params.ref_genome, checkIfExists: true)
    )

    BUILD_MINIMAP2_INDEX(reference_ch)

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
     * Align merged FASTQ files to the reference
     */
    minimap_index_ch = BUILD_MINIMAP2_INDEX.out.first()

    MINIMAP2_ALIGN(
        merged_fastq_ch,
        minimap_index_ch
    )
}
