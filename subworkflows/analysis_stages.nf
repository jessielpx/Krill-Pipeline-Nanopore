include { MSIGDB_PREPARE } from '../modules/local/msigdb_prepare/main'
nextflow.enable.dsl = 2

include { DGE } from '../modules/local/dge/main'
include { DTU } from '../modules/local/dtu/main'
include { ISOFORM_USAGE } from '../modules/local/isoform_usage/main'
include { ISOFORM_STRUCTURE_ALL } from '../modules/local/isoform_structure_all/main'
include { FUNCTIONAL_CONSEQUENCES } from '../modules/local/functional_consequences/main'
include { GSEA } from '../modules/local/gsea/main'

/*
 * Validate existing Bambu outputs before
 * launching an independent downstream stage.
 *
 * This checks file presence only.
 * Biological and statistical QC remain in
 * the corresponding analysis modules.
 */
def validateBambuFiles(bambuDir, requiredFiles) {

    requiredFiles.each { name ->

        def target = file("${bambuDir}/${name}")

        if (!target.isFile() || target.size() == 0) {
            error "Missing or empty Bambu file: ${target}"
        }
    }
}

workflow STATISTICS_MODE {

    if (!params.bambu_dir)
        error "Missing --bambu_dir"

    if (!params.samplesheet)
        error "Missing --samplesheet"

    if (!params.outdir)
        error "Missing --outdir"

    def bambu = file(
        params.bambu_dir,
        checkIfExists: true
    )

    def sheet = file(
        params.samplesheet,
        checkIfExists: true
    )

    validateBambuFiles(bambu, [
        'counts_gene.txt',
        'uniqueCounts_transcript.txt'
    ])

    DGE(
        Channel.value(bambu),
        Channel.value(sheet),
        Channel.value(file(
            "${projectDir}/bin/run_dge.R",
            checkIfExists: true
        ))
    )

    DTU(
        Channel.value(bambu),
        Channel.value(sheet),
        Channel.value(file(
            "${projectDir}/bin/run_dtu.R",
            checkIfExists: true
        ))
    )
}


/*
 * ============================================================
 * Stage 07-08: Isoform Usage + Structure
 * ============================================================
 */

workflow ISOFORMS_MODE {

    if (!params.bambu_dir)
        error "Missing --bambu_dir"

    if (!params.samplesheet)
        error "Missing --samplesheet"

    if (!params.ref_annotation)
        error "Missing --ref_annotation"

    if (!params.outdir)
        error "Missing --outdir"

    def bambu = file(
        params.bambu_dir,
        checkIfExists: true
    )

    def sheet = file(
        params.samplesheet,
        checkIfExists: true
    )

    def gtf = file(
        params.ref_annotation,
        checkIfExists: true
    )

    validateBambuFiles(bambu, [
        'uniqueCounts_transcript.txt',
        'extended_annotations.gtf'
    ])

    ISOFORM_USAGE(
        Channel.value(bambu),
        Channel.value(sheet),
        Channel.value(file(
            "${projectDir}/bin/run_isoform_usage.R",
            checkIfExists: true
        ))
    )

    ISOFORM_STRUCTURE_ALL(
        Channel.value(bambu),
        Channel.value(gtf),
        Channel.value(file(
            "${projectDir}/bin/annotate_all_isoforms.py",
            checkIfExists: true
        ))
    )
}


/*
 * ============================================================
 * Stage 09: Functional Consequences
 * ============================================================
 */

workflow FUNCTIONAL_MODE {

    for (name in [
        'usage_dir',
        'structure_dir',
        'dtu_dir',
        'samplesheet',
        'ref_annotation',
        'outdir'
    ]) {
        if (!params[name])
            error "Missing --${name}"
    }

    def usage = file(
        params.usage_dir,
        checkIfExists: true
    )

    def structure = file(
        params.structure_dir,
        checkIfExists: true
    )

    def dtu = file(
        params.dtu_dir,
        checkIfExists: true
    )

    def sheet = file(
        params.samplesheet,
        checkIfExists: true
    )

    def gtf = file(
        params.ref_annotation,
        checkIfExists: true
    )

    def required_files = [
        "${usage}/isoform_usage_supported.tsv",
        "${structure}/isoform_structure_annotations.tsv",
        "${dtu}/dtu_gene_results.tsv"
    ]

    required_files.each { name ->
        if (!file(name).exists())
            error "Missing input: ${name}"
    }

    FUNCTIONAL_CONSEQUENCES(
        Channel.value(usage),
        Channel.value(structure),
        Channel.value(dtu),
        Channel.value(sheet),
        Channel.value(gtf),
        Channel.value(file(
            "${projectDir}/bin/run_functional_consequences.py",
            checkIfExists: true
        ))
    )
}


/*
 * ============================================================
 * Stage 10: GSEA
 * ============================================================
 */

workflow GSEA_MODE {

    if (!params.dge_dir)
        error "Missing --dge_dir"

    if (!params.outdir)
        error "Missing --outdir"

    def dge = file(
        params.dge_dir,
        checkIfExists: true
    )

    if (!file("${dge}/dge_HR_vs_control.tsv").exists())
        error "Missing DGE results"

    def gsea_script = Channel.value(
        file("${projectDir}/bin/run_gsea.R",
             checkIfExists: true)
    )

    MSIGDB_PREPARE(gsea_script)

    GSEA(
        Channel.value(dge),
        MSIGDB_PREPARE.out.gmt_dir,
        gsea_script
    )
}
