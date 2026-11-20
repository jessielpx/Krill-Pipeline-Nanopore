#!/usr/bin/env Rscript

# ============================================================
# Krill Renew - Differential Transcript Usage (DRIMSeq)
#
# Inputs:
#   1. Bambu uniqueCounts_transcript.txt
#   2. Metadata TSV with sample and condition columns
#   3. Output directory
#
# Optional fourth argument:
#   --validate-only
#
# Modes:
#   pilot  - original all-sample coverage rule
#   cohort - group-aware coverage filtering
#
# Environment variables:
#   KRILL_DTU_MODE           pilot | cohort
#   KRILL_DTU_DESIGN         e.g. ~ run + subtype + group
#   KRILL_DTU_COVERAGE_MIN   default 10
#   KRILL_DTU_GROUP_FRACTION default 0.25
#
# Default pilot mode preserves the existing Nextflow interface.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L || length(args) > 4L) {
    stop(
        "Usage: run_dtu.R counts.tsv metadata.tsv ",
        "outdir [--validate-only]"
    )
}

validate_only <- (
    length(args) == 4L &&
    identical(args[4L], "--validate-only")
)

if (length(args) == 4L && !validate_only) {
    stop("Unknown option: ", args[4L])
}

counts_file <- args[1]
metadata_file <- args[2]
outdir <- args[3]

dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
)

# ============================================================
# 1. Configuration
# ============================================================

mode <- Sys.getenv(
    "KRILL_DTU_MODE",
    unset = "pilot"
)

if (!mode %in% c("pilot", "cohort")) {
    stop("KRILL_DTU_MODE must be pilot or cohort")
}

min_gene_coverage <- as.numeric(
    Sys.getenv(
        "KRILL_DTU_COVERAGE_MIN",
        unset = "10"
    )
)

group_fraction <- as.numeric(
    Sys.getenv(
        "KRILL_DTU_GROUP_FRACTION",
        unset = "0.25"
    )
)

if (
    length(min_gene_coverage) != 1L ||
    !is.finite(min_gene_coverage) ||
    min_gene_coverage < 0
) {
    stop("Invalid KRILL_DTU_COVERAGE_MIN")
}

if (
    length(group_fraction) != 1L ||
    !is.finite(group_fraction) ||
    group_fraction <= 0 ||
    group_fraction > 1
) {
    stop("KRILL_DTU_GROUP_FRACTION must be in (0,1]")
}

design_string <- Sys.getenv(
    "KRILL_DTU_DESIGN",
    unset = "~ group"
)

design_formula <- tryCatch(
    as.formula(design_string),
    error = function(e) {
        stop(
            "Invalid KRILL_DTU_DESIGN: ",
            conditionMessage(e)
        )
    }
)

if (length(design_formula) != 2L) {
    stop("Design must be a one-sided formula")
}

design_variables <- all.vars(design_formula)

if (!any(c("condition", "group") %in% design_variables)) {
    stop("Design must contain condition or group")
}

message("========================================")
message("Krill Renew - DRIMSeq")
message("Mode: ", mode)
message("Design: ", design_string)
message("Minimum gene coverage: ", min_gene_coverage)
message("Group fraction: ", group_fraction)
message("========================================")

# ============================================================
# 2. Read transcript counts
# ============================================================

if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
}

x <- read.delim(
    counts_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!all(c("TXNAME", "GENEID") %in% names(x))) {
    stop("Missing TXNAME or GENEID")
}

if (
    anyNA(x$TXNAME) ||
    anyNA(x$GENEID) ||
    any(!nzchar(x$TXNAME)) ||
    any(!nzchar(x$GENEID)) ||
    anyDuplicated(x$TXNAME)
) {
    stop("Invalid transcript or gene IDs")
}

sample_cols <- setdiff(
    names(x),
    c("TXNAME", "GENEID")
)

if (!length(sample_cols)) {
    stop("No sample columns in counts")
}

# ============================================================
# 3. Read metadata
# ============================================================

if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
}

meta <- read.delim(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!all(c("sample", "condition") %in% names(meta))) {
    stop("Metadata requires sample and condition")
}

meta$sample <- trimws(as.character(meta$sample))
meta$condition <- trimws(as.character(meta$condition))

if (
    anyNA(meta$sample) ||
    anyNA(meta$condition) ||
    any(!nzchar(meta$sample)) ||
    any(!nzchar(meta$condition)) ||
    anyDuplicated(meta$sample)
) {
    stop("Invalid sample metadata")
}

if (!setequal(sample_cols, meta$sample)) {
    stop("Count columns do not match metadata samples")
}

meta <- meta[
    match(sample_cols, meta$sample),
    ,
    drop = FALSE
]

if (!setequal(unique(meta$condition), c("control", "HR"))) {
    stop("Expected conditions: control and HR")
}

meta$condition <- factor(
    meta$condition,
    levels = c("control", "HR")
)

group_sizes <- table(meta$condition)

if (any(group_sizes < 2L)) {
    stop("At least two samples per group required")
}

# ============================================================
# 4. Prepare DRIMSeq sample table
# ============================================================

samples <- data.frame(
    sample_id = meta$sample,
    condition = factor(
        meta$condition,
        levels = c("control", "HR")
    ),
    group = factor(
        meta$condition,
        levels = c("control", "HR")
    ),
    stringsAsFactors = FALSE
)

# Preserve additional metadata columns for cohort models.

covariates <- setdiff(
    names(meta),
    c("sample", "condition", "group", "sample_id")
)

for (name in covariates) {
    values <- trimws(as.character(meta[[name]]))

    values[is.na(values) | !nzchar(values)] <- NA_character_

    samples[[name]] <- factor(values)
}

missing_design <- setdiff(
    design_variables,
    names(samples)
)

if (length(missing_design)) {
    stop(
        "Missing design covariates: ",
        paste(missing_design, collapse = ", ")
    )
}

for (name in design_variables) {
    if (anyNA(samples[[name]])) {
        stop("Missing values in design covariate: ", name)
    }
}

write.table(
    samples,
    file.path(outdir, "dtu_samples.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# 5. Validate integer counts
# ============================================================

count_mat <- as.matrix(
    x[, sample_cols, drop = FALSE]
)

suppressWarnings(
    storage.mode(count_mat) <- "numeric"
)

if (
    any(!is.finite(count_mat)) ||
    any(count_mat < 0)
) {
    stop("Counts must be finite and nonnegative")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
    stop(
        "Non-integer counts detected. ",
        "Use Bambu uniqueCounts_transcript.txt"
    )
}

if (any(count_mat > .Machine$integer.max)) {
    stop("Counts exceed integer range")
}

count_df <- as.data.frame(
    matrix(
        as.integer(round(count_mat)),
        nrow = nrow(count_mat),
        dimnames = dimnames(count_mat)
    ),
    check.names = FALSE
)

input <- data.frame(
    gene_id = x$GENEID,
    feature_id = x$TXNAME,
    count_df,
    check.names = FALSE
)

input <- input[
    order(input$gene_id, input$feature_id),
    ,
    drop = FALSE
]

message(
    "Input: ",
    nrow(input),
    " transcripts / ",
    length(unique(input$gene_id)),
    " genes"
)

# ============================================================
# 6. Initial QC
# ============================================================

qc <- data.frame(
    metric = c(
        "transcripts_input",
        "genes_input",
        "samples",
        "genes_with_multiple_isoforms",
        "noninteger_cells",
        "mode",
        "design",
        "minimum_gene_coverage",
        "group_fraction"
    ),
    value = as.character(c(
        nrow(input),
        length(unique(input$gene_id)),
        nrow(samples),
        sum(table(input$gene_id) >= 2),
        0,
        mode,
        design_string,
        min_gene_coverage,
        group_fraction
    ))
)

write.table(
    qc,
    file.path(outdir, "dtu_input_qc.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. Check design matrix
# ============================================================

full <- model.matrix(
    design_formula,
    data = samples
)

if (qr(full)$rank != ncol(full)) {
    stop(
        "Design matrix is not full rank. ",
        "Check confounding between run, subtype and group."
    )
}

if (nrow(full) <= ncol(full)) {
    stop(
        "Insufficient residual degrees of freedom ",
        "for the selected design"
    )
}

comparison_coef <- if ("condition" %in% design_variables) {
    "conditionHR"
} else {
    "groupHR"
}

if (!comparison_coef %in% colnames(full)) {
    stop(
        "Missing comparison coefficient: ",
        comparison_coef
    )
}

write.table(
    data.frame(
        sample = samples$sample_id,
        full,
        check.names = FALSE
    ),
    file.path(outdir, "dtu_design_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

if (validate_only) {
    message("Input validation completed.")
    quit(save = "no", status = 0)
}

# ============================================================
# 8. Load DRIMSeq
# ============================================================

if (!requireNamespace("DRIMSeq", quietly = TRUE)) {
    stop("DRIMSeq not installed")
}

if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel not installed")
}

d <- DRIMSeq::dmDSdata(
    counts = input,
    samples = samples
)

# ============================================================
# 9. Initial expression filtering
# ============================================================

message("Applying DRIMSeq expression filtering")

d <- DRIMSeq::dmFilter(
    d,
    min_samps_gene_expr = 3,
    min_gene_expr = 10,
    min_samps_feature_expr = 2,
    min_feature_expr = 5,
    min_samps_feature_prop = 0,
    min_feature_prop = 0
)

filtered <- DRIMSeq::counts(d)

if (!nrow(filtered)) {
    stop("No transcripts survived initial filtering")
}

message(
    "Initial filter: ",
    nrow(filtered),
    " transcripts / ",
    length(unique(filtered$gene_id)),
    " genes"
)

# ============================================================
# 10. Gene-level coverage
# ============================================================

gene_counts <- aggregate(
    filtered[, sample_cols, drop = FALSE],
    by = list(gene_id = filtered$gene_id),
    FUN = sum
)

coverage_matrix <- as.matrix(
    gene_counts[, sample_cols, drop = FALSE]
)

supported_matrix <- (
    coverage_matrix >= min_gene_coverage
)

gene_counts$n_samples_supported <- rowSums(
    supported_matrix
)

if (mode == "pilot") {

    # Preserve original Batch17 coverage rule.

    gene_counts$coverage_pass <- (
        gene_counts$n_samples_supported ==
        length(sample_cols)
    )

} else {

    # Cohort-aware filtering:
    # require support independently in each group.

    group_requirements <- pmax(
        2L,
        ceiling(
            as.numeric(group_sizes) * group_fraction
        )
    )

    names(group_requirements) <- names(group_sizes)

    for (group_name in names(group_sizes)) {

        group_samples <- meta$sample[
            meta$condition == group_name
        ]

        n_supported <- rowSums(
            supported_matrix[
                ,
                group_samples,
                drop = FALSE
            ]
        )

        gene_counts[[
            paste0("n_", group_name, "_supported")
        ]] <- n_supported
    }

    gene_counts$coverage_pass <- (
        gene_counts$n_control_supported >=
            group_requirements["control"] &
        gene_counts$n_HR_supported >=
            group_requirements["HR"]
    )

    message(
        "Cohort requirements: control >= ",
        group_requirements["control"],
        ", HR >= ",
        group_requirements["HR"]
    )
}

write.table(
    gene_counts,
    file.path(outdir, "dtu_gene_coverage_qc.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

keep_ids <- gene_counts$gene_id[
    gene_counts$coverage_pass
]

filtered <- filtered[
    filtered$gene_id %in% keep_ids,
    ,
    drop = FALSE
]

# Keep genes with at least two isoforms.

isoform_number <- table(filtered$gene_id)

keep_multi <- names(isoform_number)[
    isoform_number >= 2L
]

filtered <- filtered[
    filtered$gene_id %in% keep_multi,
    ,
    drop = FALSE
]

n_genes <- length(unique(filtered$gene_id))

if (n_genes < 2L) {
    stop("Too few genes after coverage filtering")
}

message(
    "Final filter: ",
    nrow(filtered),
    " transcripts / ",
    n_genes,
    " genes"
)

write.table(
    filtered,
    file.path(outdir, "dtu_filtered_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 11. Fit DRIMSeq model
# ============================================================

d <- DRIMSeq::dmDSdata(
    counts = filtered,
    samples = samples
)

set.seed(42)

bp <- BiocParallel::SerialParam()

message("Estimating precision")

d <- DRIMSeq::dmPrecision(
    d,
    design = full,
    BPPARAM = bp
)

message("Fitting DRIMSeq model")

d <- DRIMSeq::dmFit(
    d,
    design = full,
    BPPARAM = bp
)

message("Testing HR vs control")

d <- DRIMSeq::dmTest(
    d,
    coef = comparison_coef,
    BPPARAM = bp
)

# ============================================================
# 12. Export results
# ============================================================

gene_res <- DRIMSeq::results(
    d,
    level = "gene"
)

feature_res <- DRIMSeq::results(
    d,
    level = "feature"
)

write.table(
    gene_res,
    file.path(outdir, "dtu_gene_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

write.table(
    feature_res,
    file.path(outdir, "dtu_transcript_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

saveRDS(
    d,
    file.path(outdir, "drimseq_fit.rds")
)

writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "R_session_info.txt")
)

message("========================================")
message("Krill DTU completed")
message("Mode: ", mode)
message("Design: ", design_string)
message("Genes tested: ", nrow(gene_res))
message("Transcripts tested: ", nrow(feature_res))
message("========================================")
