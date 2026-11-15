nextflow.enable.dsl = 2

include { CORE } from '../subworkflows/core'
include { PREPROCESSING } from '../subworkflows/preprocessing'
include { GENE_EXPRESSION_EXPORT } from '../modules/local/gene_expression_export/main'
include { DGE } from '../modules/local/dge/main'
include { DTU } from '../modules/local/dtu/main'
include { ISOFORM_USAGE } from '../modules/local/isoform_usage/main'
include { ISOFORM_STRUCTURE } from '../modules/local/isoform_structure/main'
include { ISOFORM_STRUCTURE_ALL } from '../modules/local/isoform_structure_all/main'
include { FUNCTIONAL_CONSEQUENCES } from '../modules/local/functional_consequences/main'
include { GSEA } from '../modules/local/gsea/main'
include { JUNCTION_EXTRACT; JUNCTION_DJU } from '../modules/local/junction_analysis/main'
include { JAFFAL } from '../modules/local/jaffal/main'
include { FUSION_MERGE } from '../modules/local/fusion_merge/main'

/*
 * Shared raw FASTQ samplesheet parser.
 *
 * Used by CORE/ALL and ALIGNMENT.
 * Returns tuples: (meta, barcode_fastq_directory)
 */
def buildRawSamplesChannel(
    String samplesheetPath,
    String fastqRoot,
    boolean requireHRControl
) {
    return Channel
        .fromPath(samplesheetPath, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->

            def barcode = row.barcode?.toString()?.trim()
            def sample = row.sample?.toString()?.trim()
            def alias = row.alias?.toString()?.trim() ?: sample
            def condition = row.condition?.toString()?.trim()

            for (item in [
                barcode: barcode,
                sample: sample,
                alias: alias
            ]) {
                if (
                    !item.value ||
                    !(item.value ==~ /[A-Za-z0-9._-]+/) ||
                    item.value in ['.', '..']
                ) {
                    throw new IllegalArgumentException(
                        "Invalid ${item.key}: ${item.value}"
                    )
                }
            }

            if (
                requireHRControl &&
                !(condition in ['HR', 'control'])
            ) {
                throw new IllegalArgumentException(
                    "Invalid condition: ${condition}"
                )
            }

            // Optional per-sample FASTQ directory.
            // Allows barcode reuse across sequencing runs.
            def explicitPath = row.fastq_path?.toString()?.trim()

            if (!explicitPath && !fastqRoot) {
                throw new IllegalArgumentException(
                    "Sample ${sample} requires fastq_path " +
                    "or --raw_fastq_dir"
                )
            }

            def readDir = file(
                explicitPath ?: "${fastqRoot}/${barcode}",
                checkIfExists: true
            )

            if (!readDir.isDirectory()) {
                throw new IllegalArgumentException(
                    "Invalid FASTQ directory: ${readDir}"
                )
            }

            def meta = [
                barcode: barcode,
                sample: sample,
                alias: alias,
                condition: condition,
                run: row.run?.toString()?.trim(),
                subtype: row.subtype?.toString()?.trim(),
                patient: row.patient?.toString()?.trim()
            ]

            tuple(meta, readDir)
        }
        .toList()
        .flatMap { rows ->

            if (rows.isEmpty()) {
                throw new IllegalArgumentException(
                    "Samplesheet is empty"
                )
            }

            def samples = rows.collect { it[0].sample }
            def inputDirs = rows.collect {
                it[1].toAbsolutePath().normalize().toString()
            }

            if (samples.size() != samples.toSet().size()) {
                throw new IllegalArgumentException(
                    "Duplicate sample IDs"
                )
            }

            if (inputDirs.size() != inputDirs.toSet().size()) {
                throw new IllegalArgumentException(
                    "Duplicate FASTQ input directories"
                )
            }

            log.info "Validated samples: ${rows.size()}"

            rows
        }
}

workflow KRILL_CORE {

    // Full pipeline mode activates all standard downstream modules.
    def runAll = params.stage == 'all'
    def fastq_input_dir = params.raw_fastq_dir ?: params.fastq_dir


    /*
     * 1. Validate configuration
     */

    if (!(params.stage in ['core', 'all'])) {
        error "KRILL_CORE supports --stage core or all"
    }

    if (params.read_type != 'cdna') {
        error "Only --read_type cdna is supported"
    }

    if (!params.samplesheet) {
        error "Missing --samplesheet"
    }

    if (!params.ref_genome) {
        error "Missing --ref_genome"
    }

    if (!params.ref_annotation) {
        error "Missing --ref_annotation"
    }

    if (!params.outdir) {
        error "Missing --outdir"
    }

    /*
     * Structure annotation depends on
     * ISOFORM_USAGE output.
     */

    if (
        params.run_isoform_structure &&
        !runAll &&
        !params.run_isoform_usage
    ) {
        error(
            "ISOFORM_STRUCTURE requires " +
            "--run_isoform_usage true"
        )
    }


    /*
     * 2. Load input files
     */

    def samplesheet = file(
        params.samplesheet,
        checkIfExists: true
    )

    def reference_fasta_ch = Channel.value(
        file(
            params.ref_genome,
            checkIfExists: true
        )
    )

    def annotation_gtf_ch = Channel.value(
        file(
            params.ref_annotation,
            checkIfExists: true
        )
    )


    /*
     * 3. Parse samplesheet
     *
     * Required:
     *   barcode, sample, alias, condition
     *
     * Optional:
     *   run, subtype, patient
     */

    samples_ch = buildRawSamplesChannel(
        params.samplesheet,
        fastq_input_dir,
        false
    )

    /*
     * 4. Run unified CORE
     *
     * FASTQ -> Alignment -> Bambu -> QC
     */

    log.info "===================================="
    log.info "Krill Renew - BCC Long-read RNA-seq"
    log.info "===================================="

    log.info "Read type: ${params.read_type}"
    log.info "Output: ${params.outdir}"

    log.info "Run DGE: ${runAll || params.run_dge}"
    log.info "Run DTU: ${runAll || params.run_dtu}"
    log.info "Run ISOFORM_USAGE: ${runAll || params.run_isoform_usage}"
    log.info "Run ISOFORM_STRUCTURE: ${params.run_isoform_structure}"

    CORE(
        samples_ch,
        reference_fasta_ch,
        annotation_gtf_ch
    )



    /*
     * Optional gene expression export
     *
     * Bambu counts and transcript lengths
     * -> gene-level TPM
     * -> HGNC gene symbols
     *
     * Independent of SQANTI3 and DGE.
     */

    if (runAll || params.run_gene_expression_export) {

        GENE_EXPRESSION_EXPORT(
            CORE.out.bambu_results,
            annotation_gtf_ch,
            Channel.value(
                file(
                    "${projectDir}/bin/export_gene_expression.py",
                    checkIfExists: true
                )
            )
        )

    }

    /*
     * 5. Optional DGE
     *
     * DESeq2 gene-level analysis
     */

    if (runAll || params.run_dge) {

        DGE(
            CORE.out.bambu_results,
            Channel.value(samplesheet),
            Channel.value(
                file(
                    "${projectDir}/bin/run_dge.R",
                    checkIfExists: true
                )
            )
        )

    }


    /*
     * 6. Optional DTU
     *
     * DRIMSeq differential transcript usage
     *
     * Current version:
     * HR/control pilot model.
     */

    if (runAll || params.run_dtu) {

        DTU(
            CORE.out.bambu_results,
            Channel.value(samplesheet),
            Channel.value(
                file(
                    "${projectDir}/bin/run_dtu.R",
                    checkIfExists: true
                )
            )
        )

    }


    /*
     * 7. Optional ISOFORM_USAGE
     *
     * Bambu unique transcript counts
     *   -> isoform fractions
     *   -> group mean IF
     *   -> delta IF
     */

    if (runAll || params.run_isoform_usage) {

        ISOFORM_USAGE(
            CORE.out.bambu_results,
            Channel.value(samplesheet),
            Channel.value(
                file(
                    "${projectDir}/bin/run_isoform_usage.R",
                    checkIfExists: true
                )
            )
        )

    }


    /*
     * 8. Optional ISOFORM_STRUCTURE
     *
     * Input:
     *   - Bambu extended annotations
     *   - ISOFORM_USAGE results
     *   - Ensembl reference GTF
     *
     * Analysis:
     *   - Exon chain validation
     *   - CDS annotation
     *   - 5' and 3' UTR annotation
     *   - Structural annotation QC
     *
     * Automatically waits for ISOFORM_USAGE.
     */

    if (params.run_isoform_structure) {

        ISOFORM_STRUCTURE(
            CORE.out.bambu_results,
            ISOFORM_USAGE.out.results,
            annotation_gtf_ch,
            Channel.value(
                file(
                    "${projectDir}/bin/annotate_isoform_structure.py",
                    checkIfExists: true
                )
            )
        )

    }



    /*
     * 9. Optional functional consequences
     *
     * Depends on DTU, ISOFORM_USAGE,
     * and ISOFORM_STRUCTURE.
     */


    // 08_ALL: Bambu-wide structural annotation
    if (runAll || params.run_isoform_structure_all) {

        ISOFORM_STRUCTURE_ALL(
            CORE.out.bambu_results,
            annotation_gtf_ch,
            Channel.value(
                file(
                    "${projectDir}/bin/annotate_all_isoforms.py",
                    checkIfExists: true
                )
            )
        )

    }

    if (runAll || params.run_functional_consequences) {

        if (
            !runAll &&
            (
                !params.run_dtu ||
                !params.run_isoform_usage ||
                !params.run_isoform_structure_all
            )
        ) {
            error(
                "FUNCTIONAL_CONSEQUENCES requires " +
                "--run_dtu true, " +
                "--run_isoform_usage true, and " +
                "--run_isoform_structure_all true"
            )
        }

        FUNCTIONAL_CONSEQUENCES(
            ISOFORM_USAGE.out.results,
            ISOFORM_STRUCTURE_ALL.out.results,
            DTU.out.results,
            Channel.value(samplesheet),
            annotation_gtf_ch,
            Channel.value(
                file(
                    "${projectDir}/bin/run_functional_consequences.py",
                    checkIfExists: true
                )
            )
        )

    }









    /*
     * 10 - Preranked gene set enrichment
     *
     * Depends only on 05_DGE.
     */
    if (runAll || params.run_gsea) {

        if (!runAll && !params.run_dge) {
            error(
                "GSEA requires --run_dge true"
            )
        }

        def gmt_location = params.gsea_gmt_dir ?:
            "${projectDir}/assets/msigdb_human_ensembl"

        GSEA(
            DGE.out.results,
            Channel.value(
                file(gmt_location, checkIfExists: true)
            ),
            Channel.value(
                file(
                    "${projectDir}/bin/run_gsea.R",
                    checkIfExists: true
                )
            )
        )
    }



    /*
     * 11 - Alternative splicing / junction usage
     *
     * CORE indexed BAM
     *   -> junction extraction
     *   -> exploratory differential junction usage
     */

    if (runAll || params.run_alternative_splicing) {

        def junction_bams_ch = CORE.out.indexed_bam
            .map { record ->

                def items = (
                    record instanceof List
                    ? record
                    : [record]
                )

                def bams = items.findAll { item ->
                    item instanceof java.nio.file.Path &&
                    item.toString().endsWith('.bam')
                }

                if (bams.size() != 1) {
                    error(
                        "Expected exactly one BAM per " +
                        "CORE indexed_bam record; got: " +
                        record
                    )
                }

                return bams[0]
            }
            .collect()

        JUNCTION_EXTRACT(
            junction_bams_ch,
            Channel.value(samplesheet),
            annotation_gtf_ch,
            Channel.value(
                file(
                    "${projectDir}/bin/run_junction_analysis.py",
                    checkIfExists: true
                )
            )
        )

        JUNCTION_DJU(
            JUNCTION_EXTRACT.out.results,
            Channel.value(samplesheet),
            Channel.value(
                file(
                    "${projectDir}/bin/run_differential_junction_usage.R",
                    checkIfExists: true
                )
            )
        )
    }



    /*
     * 12 - Long-read fusion detection
     *
     * CORE merged FASTQ
     *   -> per-sample JAFFAL
     *   -> cohort fusion merge
     */

    if (runAll || params.run_fusion) {

        def fusion_fastq_ch = CORE.out.merged_fastq
            .map { record ->

                if (
                    !(record instanceof List) ||
                    record.size() < 2
                ) {
                    error(
                        "Expected CORE FASTQ tuple " +
                        "(meta, fastq[, ...]); got: " +
                        record
                    )
                }

                def meta = record[0]
                def fastq = record[1]

                if (
                    !(meta instanceof Map) ||
                    !meta.sample
                ) {
                    error(
                        "FASTQ metadata lacks sample ID: " +
                        meta
                    )
                }

                tuple(meta, fastq)
            }

        JAFFAL(fusion_fastq_ch)

        def fusion_csvs_ch = JAFFAL.out.calls
            .map { meta, csv -> csv }
            .collect()

        FUSION_MERGE(
            fusion_csvs_ch,
            Channel.value(samplesheet),
            Channel.value(
                file(
                    "${projectDir}/bin/merge_fusions.py",
                    checkIfExists: true
                )
            )
        )
    }

}

workflow ALIGNMENT_MODE {
    def fastq_input_dir = params.raw_fastq_dir ?: params.fastq_dir


    // ========================================================
    // 1. Validate parameters
    // ========================================================

    if (!params.samplesheet) {
        error "Missing --samplesheet"
    }

    if (!params.ref_genome) {
        error "Missing --ref_genome"
    }

    if (!params.ref_annotation) {
        error "Missing --ref_annotation"
    }

    if (!params.outdir) {
        error "Missing --outdir"
    }

    // ========================================================
    // 2. Reference inputs
    // ========================================================

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

    // ========================================================
    // 3. Read samplesheet
    // ========================================================

    samples_ch = buildRawSamplesChannel(
        params.samplesheet,
        fastq_input_dir,
        true
    )

    log.info "Krill Renew - Genome Alignment"
    log.info "FASTQ root: ${fastq_input_dir}"
    log.info "Output: ${params.outdir}"

    // ========================================================
    // 4. Run preprocessing and genome alignment
    // ========================================================

    PREPROCESSING(
        samples_ch,
        reference_fasta_ch,
        annotation_gtf_ch
    )

    // ========================================================
    // 5. Pipeline scope
    // ========================================================

    /*
     * PREPROCESSING includes:
     *
     * FASTCAT
     * BUILD_MINIMAP2_INDEX
     * PREPROCESS_ANNOTATION
     * MINIMAP2_ALIGN
     * SAMTOOLS_INDEX
     *
     * The following are NOT run:
     *
     * BAMBU_JOINT
     * SQANTI3_QC
     * UNIFIED_QC
     * DGE
     * DTU
     * ISOFORM_USAGE
     * ISOFORM_STRUCTURE
     * FUNCTIONAL_CONSEQUENCES
     * GSEA
     * JUNCTION_DJU
     * JAFFAL
     */
}
