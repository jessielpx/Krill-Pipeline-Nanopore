include { PYCOQC } from '../modules/local/pycoqc/main'
include { FASTCAT } from '../modules/local/fastcat/main'

workflow QC {

    take:
    summary_ch
    samples_ch

    main:
    PYCOQC(summary_ch)
    FASTCAT(samples_ch)

    emit:
    fastcat_results = FASTCAT.out.results
    pycoqc_html      = PYCOQC.out.html
    pycoqc_json      = PYCOQC.out.json
}
