/*
 * ============================================================
 * Krill Renew - 05 DGE
 *
 * Gene-level differential expression using DESeq2
 *
 * Inputs:
 *   - Bambu counts_gene.txt
 *   - Samplesheet CSV
 *   - run_dge.R
 *
 * Comparison:
 *   HR vs control
 *
 * Optional Nextflow parameters:
 *   --dge_design "~ condition"
 *   --dge_min_count 10
 *   --dge_min_samples 2
 *
 * Default design and filtering preserve the original
 * Batch17 pilot analysis.
 *
 * For full cohort:
 *   --dge_design "~ run + subtype + condition"
 *
 * Native execution on Rorqual.
 * No containers.
 * ============================================================
 */

process DGE {

    tag "HR_vs_control"

    cpus 2
    memory '16 GB'
    time '4h'

    publishDir "${params.outdir}/05_dge",
        mode: 'copy'

    input:
    path bambu_results
    path samplesheet
    path dge_script

    output:
    path 'dge_results', emit: results

    shell:
    '''
    set -euo pipefail

    # ========================================================
    # 1. Load R environment
    # ========================================================

    module --force purge

    module load StdEnv/2023
    module load r/4.6.1

    export R_LIBS_USER="!{params.r_libs_user ?: ''}"

    export BIOCONDUCTOR_CONFIG_FILE="!{params.bioconductor_config_file ?: ''}"

    if [[ -z "$R_LIBS_USER" ]]; then
        echo "ERROR: Set r_libs_user in params.yaml" >&2
        exit 1
    fi

    if [[ ! -d "$R_LIBS_USER" ]]; then
        echo "ERROR: R library directory not found: $R_LIBS_USER" >&2
        exit 1
    fi

    if [[ -z "${BIOCONDUCTOR_CONFIG_FILE:-}" ]]; then
        echo "ERROR: Set bioconductor_config_file in params.yaml" >&2
        exit 1
    fi

    BIOC_CONFIG_PATH="${BIOCONDUCTOR_CONFIG_FILE#file://}"

    if [[ ! -s "$BIOC_CONFIG_PATH" ]]; then
        echo "ERROR: Bioconductor config not found: $BIOC_CONFIG_PATH" >&2
        exit 1
    fi


    # ========================================================
    # 2. DGE configuration
    # ========================================================

    export KRILL_DGE_DESIGN='!{params.dge_design ?: "~ condition"}'

    export KRILL_DGE_MIN_COUNT='!{params.dge_min_count != null ? params.dge_min_count : 10}'

    DGE_MIN_SAMPLES='!{params.dge_min_samples != null ? params.dge_min_samples : ""}'

    if [[ -n "$DGE_MIN_SAMPLES" ]]; then
        export KRILL_DGE_MIN_SAMPLES="$DGE_MIN_SAMPLES"
    else
        unset KRILL_DGE_MIN_SAMPLES
    fi

    echo "========================================"
    echo "Krill Renew - DGE / DESeq2"
    echo "Design: $KRILL_DGE_DESIGN"
    echo "Minimum count: $KRILL_DGE_MIN_COUNT"
    echo "Minimum samples: ${KRILL_DGE_MIN_SAMPLES:-group minimum}"
    echo "========================================"

    # ========================================================
    # 3. Validate inputs
    # ========================================================

    COUNTS="!{bambu_results}/counts_gene.txt"
    SAMPLESHEET="!{samplesheet}"
    SCRIPT="!{dge_script}"

    for required in \
        "$COUNTS" \
        "$SAMPLESHEET" \
        "$SCRIPT"
    do
        if [[ ! -s "$required" ]]; then
            echo "ERROR: Missing input: $required" >&2
            exit 1
        fi
    done

    # ========================================================
    # 4. Validate DESeq2
    # ========================================================

    Rscript --vanilla -e '
        stopifnot(
            requireNamespace("DESeq2", quietly = TRUE)
        )

        cat(
            "DESeq2 version:",
            as.character(packageVersion("DESeq2")),
            "\\n"
        )
    '

    # ========================================================
    # 5. Run DGE
    # ========================================================

    echo "===== RUNNING KRILL DGE ====="

    Rscript --vanilla "$SCRIPT" \
        "$COUNTS" \
        "$SAMPLESHEET" \
        dge_results

    # ========================================================
    # 6. Validate output files
    # ========================================================

    for output in \
        dge_results/dge_HR_vs_control.tsv \
        dge_results/normalized_gene_counts.tsv \
        dge_results/dge_sample_qc.tsv \
        dge_results/dge_input_qc.tsv \
        dge_results/dge_design_matrix.tsv \
        dge_results/deseq2_dds.rds \
        dge_results/R_session_info.txt
    do

        if [[ ! -s "$output" ]]; then
            echo "ERROR: Missing output: $output" >&2
            exit 1
        fi

    done

    echo "========================================"
    echo "Krill DGE completed successfully"
    echo "Design: $KRILL_DGE_DESIGN"
    echo "========================================"

    '''
}
