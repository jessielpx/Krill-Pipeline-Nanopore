nextflow.enable.dsl = 2

/*
 * ==========================================
 * Krill Renew
 *
 * Nanopore long-read cDNA RNA-seq
 * Basal Cell Carcinoma (BCC)
 *
 * CORE:
 *   1. FASTQ preparation
 *   2. Genome alignment (minimap2)
 *   3. BAM indexing and alignment QC
 *   4. Joint transcript discovery (bambu)
 *   5. Gene/transcript quantification
 *   6. SQANTI3 structural QC
 *   7. Unified QC
 *
 * Optional downstream:
 *   - DGE (DESeq2)
 *   - DTU (DRIMSeq)
 *   - ISOFORM_USAGE (isoform fraction / delta IF)
 *   - ISOFORM_STRUCTURE (CDS / UTR annotation)
 *
 * Native execution on Rorqual.
 * No containers.
 *
 * Current downstream statistical design:
 *   HR vs control pilot.
 *   Cohort-ready analysis pending.
 * ==========================================
 */



// ==========================================
// Parameters
// ==========================================

params.stage = params.containsKey('stage') ? params.stage : 'core'
params.read_type = params.containsKey('read_type') ? params.read_type : 'cdna'

params.run_gene_expression_export = params.containsKey('run_gene_expression_export') ? params.run_gene_expression_export : false
params.run_dge = params.containsKey('run_dge') ? params.run_dge : false
params.run_dtu = params.containsKey('run_dtu') ? params.run_dtu : false
params.run_isoform_usage = params.containsKey('run_isoform_usage') ? params.run_isoform_usage : false
params.run_isoform_structure = params.containsKey('run_isoform_structure') ? params.run_isoform_structure : false
params.run_functional_consequences = params.containsKey('run_functional_consequences') ? params.run_functional_consequences : false


// ==========================================
params.run_isoform_structure_all = params.containsKey('run_isoform_structure_all') ? params.run_isoform_structure_all : false

// Main workflow
// ==========================================




params.run_alternative_splicing = params.containsKey('run_alternative_splicing') ? params.run_alternative_splicing : false



include {
    KRILL_CORE;
    ALIGNMENT_MODE
} from './subworkflows/raw_read_stages'


include {
    BAMBU_MODE;
    POST_BAMBU_MODE
} from './subworkflows/transcriptome_stages'




/*
 * ============================================================
 * Stage 05-06: DGE + DTU
 * ============================================================
 */


/*
 * Independent downstream analysis stages (05-10).
 */
include {
    STATISTICS_MODE;
    ISOFORMS_MODE;
    FUNCTIONAL_MODE;
    GSEA_MODE
} from './subworkflows/analysis_stages'


include {
    JUNCTION_MODE;
    FUSION_MODE
} from './subworkflows/remaining_stages'


/*
 * Unified pipeline entrypoint.
 *
 * Existing analysis implementations are unchanged.
 */

workflow {

    if (params.stage in ['core', 'all']) {
        KRILL_CORE()
    }
    else if (params.stage == 'alignment') {
        ALIGNMENT_MODE()
    }
    else if (params.stage == 'bambu') {
        BAMBU_MODE()
    }
    else if (params.stage == 'post_bambu') {
        POST_BAMBU_MODE()
    }
    else if (params.stage == 'statistics') {
        STATISTICS_MODE()
    }
    else if (params.stage == 'isoforms') {
        ISOFORMS_MODE()
    }
    else if (params.stage == 'functional') {
        FUNCTIONAL_MODE()
    }
    else if (params.stage == 'gsea') {
        GSEA_MODE()
    }
    else if (params.stage == 'junction') {
        JUNCTION_MODE()
    }
    else if (params.stage == 'fusion') {
        FUSION_MODE()
    }
    else if (params.stage == 'help') {
        log.info "Krill Renew available stages: all, core, alignment, bambu, post_bambu, statistics, isoforms, functional, gsea, junction, fusion"
    }
    else {
        error "Unsupported --stage: ${params.stage}"
    }
}
