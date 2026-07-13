nextflow.enable.dsl = 2

include { PYCOQC } from './modules/local/pycoqc/main'

workflow {

    if (!params.summary) {
        error "Please provide --summary /path/to/sequencing_summary.txt"
    }

    summary_ch = Channel.of(
        tuple(
            [id: params.run_id],
            file(params.summary, checkIfExists: true)
        )
    )

    PYCOQC(summary_ch)
}
