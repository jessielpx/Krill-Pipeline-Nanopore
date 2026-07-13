nextflow.enable.dsl = 2

include { PYCOQC } from './modules/local/pycoqc/main'
include { FASTCAT } from './modules/local/fastcat/main'

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
     * Read sample metadata and locate each barcode directory
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
     * Barcode-level FASTQ concatenation
     */
    FASTCAT(samples_ch)
}
