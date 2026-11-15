nextflow.enable.dsl = 2

include { JUNCTION_EXTRACT; JUNCTION_DJU } from '../modules/local/junction_analysis/main'
include { JAFFAL } from '../modules/local/jaffal/main'
include { FUSION_MERGE } from '../modules/local/fusion_merge/main'

/*
 * Shared validation for published-sample analysis.
 */
def publishedSampleRows(boolean requireComparisonGroups) {
    return Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .toList()
        .map { rows ->

            if (params.containsKey('expected_samples') && params.expected_samples != null && rows.size() != (params.expected_samples as Integer)) {
                error "Expected ${params.expected_samples} samples, found ${rows.size()}"
            }

            def samples = rows.collect {
                it.sample?.toString()?.trim()
            }

            if (samples.any { !it || !(it ==~ /[A-Za-z0-9._-]+/) }) {
                error "Invalid sample ID in samplesheet"
            }

            if (samples.size() != samples.toSet().size()) {
                error "Duplicate sample IDs"
            }

            def groups = rows.collect {
                it.condition?.toString()?.trim()
            }

            if (groups.any { !it }) {
                error "Missing condition in samplesheet"
            }

            if (requireComparisonGroups) {
                if (groups.toSet() != ['HR', 'control'].toSet()) {
                    error "Junction/DJU requires HR and control"
                }

                if (
                    groups.count('HR') < 2 ||
                    groups.count('control') < 2
                ) {
                    error "Junction/DJU requires at least two samples per group"
                }
            }

            log.info "Stage: ${params.stage}"
            log.info "Samples: ${rows.size()}"
            log.info "HR: ${groups.count('HR')}"
            log.info "Control: ${groups.count('control')}"

            rows.sort {
                a, b -> a.sample.toString() <=> b.sample.toString()
            }
        }
}

workflow JUNCTION_MODE {

    for (name in ['samplesheet', 'outdir', 'bam_dir', 'ref_annotation']) {
        if (!params[name])
            error "Missing --${name}"
    }

    samplesheet_ch = Channel.value(
        file(params.samplesheet, checkIfExists: true)
    )

    rows_ch = publishedSampleRows(true)


        bam_ch = rows_ch.map { rows ->

            rows.collect { row ->
                def sample = row.sample.toString().trim()

                file(
                    "${params.bam_dir}/${sample}/${sample}.sorted.bam",
                    checkIfExists: true
                )
            }
        }

        JUNCTION_EXTRACT(
            bam_ch,
            samplesheet_ch,
            Channel.value(file(params.ref_annotation, checkIfExists: true)),
            Channel.value(file(
                "${projectDir}/bin/run_junction_analysis.py",
                checkIfExists: true
            ))
        )

        JUNCTION_DJU(
            JUNCTION_EXTRACT.out.results,
            samplesheet_ch,
            Channel.value(file(
                "${projectDir}/bin/run_differential_junction_usage.R",
                checkIfExists: true
            ))
        )
}

workflow FUSION_MODE {

    def fastq_input_dir = params.merged_fastq_dir ?: params.fastq_dir

    for (name in ['samplesheet', 'outdir']) {
        if (!params[name])
            error "Missing --${name}"
    }

    if (!fastq_input_dir)
        error "Fusion stage requires --merged_fastq_dir (or legacy --fastq_dir)"

    samplesheet_ch = Channel.value(
        file(params.samplesheet, checkIfExists: true)
    )

    rows_ch = publishedSampleRows(false)


        fastq_ch = rows_ch.flatMap { rows ->

            rows.collect { row ->
                def sample = row.sample.toString().trim()

                def meta = [
                    sample: sample,
                    condition: row.condition.toString().trim()
                ]

                def fastq = file(
                    "${fastq_input_dir}/${sample}/merged_fastq/${sample}.fastq.gz",
                    checkIfExists: true
                )

                tuple(meta, fastq)
            }
        }

        JAFFAL(fastq_ch)

        fusion_csv_ch = JAFFAL.out.calls
            .map { meta, csv -> csv }
            .collect()

        FUSION_MERGE(
            fusion_csv_ch,
            samplesheet_ch,
            Channel.value(file(
                "${projectDir}/bin/merge_fusions.py",
                checkIfExists: true
            ))
        )
}
