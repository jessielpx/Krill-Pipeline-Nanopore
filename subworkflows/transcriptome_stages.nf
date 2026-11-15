nextflow.enable.dsl = 2

include { BAMBU_JOINT } from '../modules/local/bambu_joint/main'
include { SQANTI3_QC } from '../modules/local/sqanti3_qc/main'
include { GENE_EXPRESSION_EXPORT } from '../modules/local/gene_expression_export/main'

workflow BAMBU_MODE {

    if (!params.bam_dir)
        error "Missing --bam_dir"

    if (!params.ref_annotation)
        error "Missing --ref_annotation"

    if (!params.ref_genome)
        error "Missing --ref_genome"

    if (!params.outdir)
        error "Missing --outdir"

    // Find all published alignment BAMs.

    indexed_bams_ch = Channel
        .fromPath(
            "${params.bam_dir}/*/*.sorted.bam",
            checkIfExists: true
        )
        .toSortedList { a, b ->
            a.getFileName().toString() <=>
            b.getFileName().toString()
        }
        .map { bams ->

            if (params.expected_samples != null && bams.size() != (params.expected_samples as Integer)) {
                error(
                    "Expected ${params.expected_samples} BAMs, " +
                    "found ${bams.size()}"
                )
            }

            bams.each { bam ->
                if (
                    !java.nio.file.Files.isRegularFile(bam) ||
                    java.nio.file.Files.size(bam) == 0
                ) {
                    error "Missing or empty BAM: ${bam}"
                }
            }

            def sample_names = bams.collect {
                it.getFileName()
                  .toString()
                  .replaceFirst(/\.sorted\.bam$/, '')
            }

            if (sample_names.any {
                !it ||
                !(it ==~ /[A-Za-z0-9._-]+/) ||
                it in ['.', '..']
            }) {
                error "Invalid BAM sample ID"
            }

            if (sample_names.size() != sample_names.toSet().size()) {
                error "Duplicate BAM sample IDs"
            }

            log.info "===================================="
            log.info "Krill Renew - BAMBU JOINT"
            log.info "Samples: ${bams.size()}"
            log.info "BAM input: ${params.bam_dir}"
            log.info "===================================="

            bams
        }

    // Match each BAM to its existing BAI index.

    indexed_bais_ch = indexed_bams_ch.map { bams ->

        bams.collect { bam ->
            def bai = file(
                "${bam}.bai",
                checkIfExists: true
            )

            if (
                !java.nio.file.Files.isRegularFile(bai) ||
                java.nio.file.Files.size(bai) == 0
            ) {
                error "Missing or empty BAI: ${bai}"
            }

            bai
        }
    }

    annotation_ch = Channel.value(
        file(
            params.ref_annotation,
            checkIfExists: true
        )
    )

    genome_ch = Channel.value(
        file(
            params.ref_genome,
            checkIfExists: true
        )
    )

    bambu_script_ch = Channel.value(
        file(
            "${projectDir}/bin/run_bambu.R",
            checkIfExists: true
        )
    )

    BAMBU_JOINT(
        indexed_bams_ch,
        indexed_bais_ch,
        annotation_ch,
        genome_ch,
        bambu_script_ch
    )
}


workflow POST_BAMBU_MODE {

    if (!params.bambu_dir)
        error "Missing --bambu_dir"

    if (!params.ref_annotation)
        error "Missing --ref_annotation"

    if (!params.ref_genome)
        error "Missing --ref_genome"

    if (!params.outdir)
        error "Missing --outdir"

    def bambu_dir = file(
        params.bambu_dir,
        checkIfExists: true
    )

    // Validate completed Bambu outputs.
    [
        'bambu_transcripts.rds',
        'extended_annotations.gtf',
        'counts_transcript.txt',
        'counts_gene.txt',
        'sample_manifest.tsv'
    ].each { filename ->

        if (!file("${params.bambu_dir}/${filename}").exists()) {
            error "Missing Bambu output: ${filename}"
        }
    }

    bambu_ch = Channel.value(bambu_dir)

    annotation_ch = Channel.value(
        file(
            params.ref_annotation,
            checkIfExists: true
        )
    )

    genome_ch = Channel.value(
        file(
            params.ref_genome,
            checkIfExists: true
        )
    )

    export_script_ch = Channel.value(
        file(
            "${projectDir}/bin/export_gene_expression.py",
            checkIfExists: true
        )
    )

    log.info "===================================="
    log.info "Krill Renew - POST BAMBU"
    log.info "===================================="
    log.info "Bambu: ${params.bambu_dir}"
    log.info "Output: ${params.outdir}"
    log.info "Running SQANTI3 and Gene Expression Export"

    // Independent modules: Nextflow can run
    // them concurrently when resources permit.

    SQANTI3_QC(
        bambu_ch,
        annotation_ch,
        genome_ch
    )

    GENE_EXPRESSION_EXPORT(
        bambu_ch,
        annotation_ch,
        export_script_ch
    )
}
