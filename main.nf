nextflow.enable.dsl = 2

include { PYCOQC } from './modules/local/pycoqc/main'
include { FASTCAT } from './modules/local/fastcat/main'

workflow {

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
     * Barcode-level FASTQ concatenation
     */
    barcode_ch = Channel
        .fromPath(
            "${params.fastq_dir}/barcode*",
            type: 'dir',
            checkIfExists: true
        )
        .map { barcode_dir ->
            tuple(barcode_dir.baseName, barcode_dir)
        }

    FASTCAT(barcode_ch)
}
