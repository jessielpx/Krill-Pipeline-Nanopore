#!/usr/bin/env Rscript

# ============================================================
# Krill Renew - 05 DGE (DESeq2)
#
# Inputs:
#   1. Bambu counts_gene.txt
#   2. Samplesheet CSV
#   3. Output directory
#
# Comparison: HR vs control
#
# Optional environment settings:
#   KRILL_DGE_DESIGN
#     Default: ~ condition
#     Example: ~ run + subtype + condition
#
#   KRILL_DGE_MIN_COUNT
#     Default: 10
#
#   KRILL_DGE_MIN_SAMPLES
#     Default: minimum condition group size
#
# The default settings preserve the original
# Batch17 analysis for regression testing.
#
# Note:
#   Bambu estimated counts can be fractional.
#   This workflow rounds them for DESeq2.
#   This approximation is explicitly recorded.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
    stop(
        "Usage: Rscript run_dge.R ",
        "counts_gene.txt samplesheet.csv outdir"
    )
}

counts_file <- args[1]
samplesheet_file <- args[2]
outdir <- args[3]

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

suppressPackageStartupMessages({
    library(DESeq2)
})

# ============================================================
# 1. Configuration
# ============================================================

design_string <- Sys.getenv(
    "KRILL_DGE_DESIGN",
    unset = "~ condition"
)

min_count <- as.numeric(
    Sys.getenv(
        "KRILL_DGE_MIN_COUNT",
        unset = "10"
    )
)

min_samples_setting <- Sys.getenv(
    "KRILL_DGE_MIN_SAMPLES",
    unset = ""
)

if (
    length(min_count) != 1L ||
    !is.finite(min_count) ||
    min_count < 0
) {
    stop("Invalid KRILL_DGE_MIN_COUNT")
}

design_formula <- tryCatch(
    as.formula(design_string),
    error = function(e) {
        stop(
            "Invalid DGE design: ",
            conditionMessage(e)
        )
    }
)

if (length(design_formula) != 2L) {
    stop("Design must be a one-sided formula")
}

design_variables <- all.vars(design_formula)

if (!"condition" %in% design_variables) {
    stop("Design must include condition")
}

# Only additive models are supported in this version.
term_labels <- attr(
    terms(design_formula),
    "term.labels"
)

if (
    any(grepl(":", term_labels, fixed = TRUE)) ||
    any(grepl("*", term_labels, fixed = TRUE))
) {
    stop(
        "Interaction terms are not supported. ",
        "Use an additive design."
    )
}

message("========================================")
message("Krill Renew - DGE / DESeq2")
message("Design: ", design_string)
message("Minimum count: ", min_count)
message("========================================")

# ============================================================
# 2. Read gene counts
# ============================================================

if (!file.exists(counts_file)) {
    stop("Counts file not found")
}

count_df <- read.delim(
    counts_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!"GENEID" %in% names(count_df)) {
    stop("Missing GENEID column")
}

gene_ids <- as.character(count_df$GENEID)

if (
    anyNA(gene_ids) ||
    any(!nzchar(gene_ids)) ||
    anyDuplicated(gene_ids)
) {
    stop("Invalid or duplicate gene IDs")
}

sample_cols <- setdiff(
    names(count_df),
    "GENEID"
)

if (!length(sample_cols)) {
    stop("No sample columns found")
}

count_matrix <- as.matrix(
    count_df[, sample_cols, drop = FALSE]
)

rownames(count_matrix) <- gene_ids

suppressWarnings(
    storage.mode(count_matrix) <- "numeric"
)

if (
    any(!is.finite(count_matrix)) ||
    any(count_matrix < 0)
) {
    stop("Invalid gene count values")
}

# DESeq2 requires nonnegative integer counts.
# Keep the original rounding procedure.

fractional <- sum(
    abs(count_matrix - round(count_matrix)) > 1e-6
)

if (fractional > 0) {
    message(
        "Rounding ",
        fractional,
        " fractional estimated count values"
    )
}

count_matrix <- round(count_matrix)

if (any(count_matrix > .Machine$integer.max)) {
    stop("Count values exceed integer range")
}

storage.mode(count_matrix) <- "integer"

# ============================================================
# 3. Read sample metadata
# ============================================================

if (!file.exists(samplesheet_file)) {
    stop("Samplesheet file not found")
}

metadata <- read.csv(
    samplesheet_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (
    !all(c("sample", "condition") %in% names(metadata))
) {
    stop("Samplesheet requires sample and condition")
}

metadata$sample <- trimws(
    as.character(metadata$sample)
)

metadata$condition <- trimws(
    as.character(metadata$condition)
)

if (
    anyNA(metadata$sample) ||
    anyNA(metadata$condition) ||
    any(!nzchar(metadata$sample)) ||
    any(!nzchar(metadata$condition)) ||
    anyDuplicated(metadata$sample)
) {
    stop("Invalid sample metadata")
}

if (
    !setequal(
        colnames(count_matrix),
        metadata$sample
    )
) {
    stop(
        "Sample IDs do not match between ",
        "counts and metadata"
    )
}

metadata <- metadata[
    match(
        colnames(count_matrix),
        metadata$sample
    ),
    ,
    drop = FALSE
]

rownames(metadata) <- metadata$sample

if (
    !setequal(
        unique(metadata$condition),
        c("control", "HR")
    )
) {
    stop("Expected conditions: control and HR")
}

metadata$condition <- factor(
    metadata$condition,
    levels = c("control", "HR")
)

group_sizes <- table(metadata$condition)

if (any(group_sizes < 2L)) {
    stop(
        "At least two samples per condition required"
    )
}

# ============================================================
# 4. Prepare covariates
# ============================================================

additional_variables <- setdiff(
    design_variables,
    "condition"
)

missing_variables <- setdiff(
    additional_variables,
    names(metadata)
)

if (length(missing_variables)) {
    stop(
        "Missing design variables: ",
        paste(missing_variables, collapse = ", ")
    )
}

for (variable in additional_variables) {

    values <- metadata[[variable]]

    if (is.character(values)) {

        values <- trimws(values)

        values[
            is.na(values) | !nzchar(values)
        ] <- NA_character_

        metadata[[variable]] <- factor(values)
    }

    if (is.factor(metadata[[variable]])) {
        metadata[[variable]] <- droplevels(
            metadata[[variable]]
        )
    }

    if (anyNA(metadata[[variable]])) {
        stop(
            "Missing values in covariate: ",
            variable
        )
    }
}

if (any(colSums(count_matrix) == 0)) {
    stop("At least one sample has zero total counts")
}

# ============================================================
# 5. Validate design matrix
# ============================================================

design_matrix <- model.matrix(
    design_formula,
    data = metadata
)

if (
    qr(design_matrix)$rank != ncol(design_matrix)
) {
    stop(
        "Design matrix is not full rank. ",
        "Check confounding among run, subtype ",
        "and condition."
    )
}

if (nrow(design_matrix) <= ncol(design_matrix)) {
    stop(
        "Insufficient residual degrees of freedom ",
        "for the selected statistical design."
    )
}

if (
    !"conditionHR" %in% colnames(design_matrix)
) {
    stop(
        "Missing conditionHR coefficient. ",
        "Use condition as a main effect ",
        "with control as the reference."
    )
}

write.table(
    data.frame(
        sample = metadata$sample,
        design_matrix,
        check.names = FALSE
    ),
    file.path(
        outdir,
        "dge_design_matrix.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. Prefilter low-expression genes
# ============================================================

min_replicates <- min(group_sizes)

if (nzchar(min_samples_setting)) {

    min_replicates <- as.numeric(
        min_samples_setting
    )

    if (
        length(min_replicates) != 1L ||
        !is.finite(min_replicates) ||
        min_replicates < 1 ||
        min_replicates != floor(min_replicates) ||
        min_replicates > ncol(count_matrix)
    ) {
        stop("Invalid KRILL_DGE_MIN_SAMPLES")
    }
}

keep <- (
    rowSums(count_matrix >= min_count) >=
    min_replicates
)

message(
    "Genes before filtering: ",
    nrow(count_matrix)
)

message(
    "Genes after filtering: ",
    sum(keep)
)

if (!any(keep)) {
    stop("No genes remain after filtering")
}

# ============================================================
# 7. Export input QC
# ============================================================

input_qc <- data.frame(
    metric = c(
        "genes_input",
        "genes_after_filter",
        "samples",
        "control_samples",
        "HR_samples",
        "fractional_count_cells",
        "count_processing",
        "design",
        "minimum_count",
        "minimum_samples"
    ),
    value = as.character(c(
        nrow(count_matrix),
        sum(keep),
        ncol(count_matrix),
        group_sizes["control"],
        group_sizes["HR"],
        fractional,
        "round_estimated_counts",
        design_string,
        min_count,
        min_replicates
    ))
)

write.table(
    input_qc,
    file.path(
        outdir,
        "dge_input_qc.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. DESeq2 analysis
# ============================================================

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix[
        keep,
        ,
        drop = FALSE
    ],
    colData = metadata,
    design = design_formula
)

dds <- DESeq(dds)

res <- results(
    dds,
    contrast = c(
        "condition",
        "HR",
        "control"
    ),
    alpha = 0.05
)

# ============================================================
# 9. Differential expression results
# ============================================================

results_df <- data.frame(
    GENEID = rownames(res),
    as.data.frame(res),
    check.names = FALSE
)

results_df <- results_df[
    order(
        results_df$padj,
        results_df$pvalue,
        na.last = TRUE
    ),
    ,
    drop = FALSE
]

write.table(
    results_df,
    file.path(
        outdir,
        "dge_HR_vs_control.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# 10. Normalized gene counts
# ============================================================

norm_counts <- counts(
    dds,
    normalized = TRUE
)

write.table(
    data.frame(
        GENEID = rownames(norm_counts),
        norm_counts,
        check.names = FALSE
    ),
    file.path(
        outdir,
        "normalized_gene_counts.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 11. Sample normalization QC
# ============================================================

sample_qc <- data.frame(
    sample = colnames(dds),
    condition = metadata$condition,
    input_counts = colSums(count_matrix),
    size_factor = sizeFactors(dds)
)

write.table(
    sample_qc,
    file.path(
        outdir,
        "dge_sample_qc.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 12. Reproducibility
# ============================================================

saveRDS(
    dds,
    file.path(
        outdir,
        "deseq2_dds.rds"
    )
)

writeLines(
    capture.output(sessionInfo()),
    file.path(
        outdir,
        "R_session_info.txt"
    )
)

significant <- (
    !is.na(res$padj) &
    res$padj < 0.05
)

message("========================================")
message("Krill DGE completed")
message("Design: ", design_string)
message("Samples: ", ncol(dds))
message("Genes tested: ", nrow(dds))
message(
    "Genes FDR < 0.05: ",
    sum(significant)
)
message("Output: ", outdir)
message("========================================")
