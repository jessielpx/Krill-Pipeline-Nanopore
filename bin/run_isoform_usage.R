#!/usr/bin/env Rscript

# ============================================================
# Krill Renew - 07 ISOFORM_USAGE
#
# Per-sample transcript isoform fractions
# and descriptive HR vs control delta IF.
#
# Usage:
# Rscript run_isoform_usage.R counts.tsv metadata.tsv outdir [min_gene_counts]
#
# Environment:
# KRILL_ISOFORM_MODE=pilot|cohort
# KRILL_ISOFORM_GROUP_FRACTION=0.25
#
# Default mode: pilot (backward compatible)
#
# No differential statistical test is performed here.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L || length(args) > 4L) {
    stop(
        "Usage: run_isoform_usage.R ",
        "counts.tsv metadata.tsv outdir [min_gene_counts]"
    )
}

counts_file <- args[1]
metadata_file <- args[2]
outdir <- args[3]

min_gene_counts <- if (length(args) == 4L) {
    as.numeric(args[4])
} else {
    10
}

mode <- Sys.getenv(
    "KRILL_ISOFORM_MODE",
    unset = "pilot"
)

group_fraction <- as.numeric(
    Sys.getenv(
        "KRILL_ISOFORM_GROUP_FRACTION",
        unset = "0.25"
    )
)

if (!mode %in% c("pilot", "cohort")) {
    stop("KRILL_ISOFORM_MODE must be pilot or cohort")
}

if (
    length(min_gene_counts) != 1L ||
    !is.finite(min_gene_counts) ||
    min_gene_counts < 0
) {
    stop("min_gene_counts must be >= 0")
}

if (
    length(group_fraction) != 1L ||
    !is.finite(group_fraction) ||
    group_fraction <= 0 ||
    group_fraction > 1
) {
    stop("Invalid KRILL_ISOFORM_GROUP_FRACTION")
}

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

message("========================================")
message("Krill Renew - ISOFORM_USAGE")
message("Mode: ", mode)
message("Minimum gene counts: ", min_gene_counts)
message("Group fraction: ", group_fraction)
message("========================================")

# ============================================================
# 1. Read Bambu transcript counts
# ============================================================

x <- read.delim(
    counts_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!all(c("TXNAME", "GENEID") %in% names(x))) {
    stop("Missing TXNAME / GENEID")
}

if (
    anyNA(x$TXNAME) ||
    anyNA(x$GENEID) ||
    any(!nzchar(x$TXNAME)) ||
    any(!nzchar(x$GENEID)) ||
    anyDuplicated(x$TXNAME)
) {
    stop("Invalid transcript/gene IDs")
}

samp <- setdiff(
    names(x),
    c("TXNAME", "GENEID")
)

if (!length(samp)) {
    stop("No sample count columns found")
}

# ============================================================
# 2. Read metadata
# ============================================================

meta <- read.delim(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (!all(c("sample", "condition") %in% names(meta))) {
    stop("Missing sample / condition columns")
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

if (!setequal(samp, meta$sample)) {
    stop("Sample names do not match metadata")
}

meta <- meta[
    match(samp, meta$sample),
    ,
    drop = FALSE
]

if (!setequal(unique(meta$condition), c("HR", "control"))) {
    stop("Expected HR and control")
}

hr <- meta$sample[
    meta$condition == "HR"
]

ctl <- meta$sample[
    meta$condition == "control"
]

if (length(hr) < 2L || length(ctl) < 2L) {
    stop("At least two samples per group required")
}

# ============================================================
# 3. Validate numeric counts
# ============================================================

mat <- as.matrix(
    x[, samp, drop = FALSE]
)

suppressWarnings(
    storage.mode(mat) <- "numeric"
)

if (
    any(!is.finite(mat)) ||
    any(mat < 0)
) {
    stop(
        "Counts must be finite, numeric, and nonnegative"
    )
}

# ============================================================
# 4. Calculate gene totals and isoform fractions
#
# Denominators always include ALL transcripts
# belonging to the gene.
# ============================================================

genes <- unique(x$GENEID)

gene_index <- match(
    x$GENEID,
    genes
)

gene_totals <- rowsum(
    mat,
    group = x$GENEID,
    reorder = FALSE
)

gene_totals <- gene_totals[
    match(genes, rownames(gene_totals)),
    ,
    drop = FALSE
]

stopifnot(
    identical(rownames(gene_totals), genes)
)

denoms <- gene_totals[
    gene_index,
    ,
    drop = FALSE
]

# Gene total == 0 means the isoform fraction
# is undefined, not zero.

fractions <- mat / ifelse(
    denoms > 0,
    denoms,
    NA_real_
)

colnames(fractions) <- paste0(
    samp,
    "_IF"
)

# ============================================================
# 5. Gene-level coverage filtering
# ============================================================

support_mat <- (
    gene_totals[, samp, drop = FALSE] >=
    min_gene_counts
)

n_supported <- rowSums(support_mat)

n_control_supported <- rowSums(
    support_mat[, ctl, drop = FALSE]
)

n_HR_supported <- rowSums(
    support_mat[, hr, drop = FALSE]
)

if (mode == "pilot") {

    # Original Batch17 criterion:
    # each sample must pass.

    supported <- (
        n_supported == length(samp)
    )

    required_control <- length(ctl)
    required_HR <- length(hr)

} else {

    # Cohort criterion:
    # adequate support independently
    # within each condition.

    required_control <- max(
        2L,
        ceiling(
            length(ctl) * group_fraction
        )
    )

    required_HR <- max(
        2L,
        ceiling(
            length(hr) * group_fraction
        )
    )

    supported <- (
        n_control_supported >= required_control &
        n_HR_supported >= required_HR
    )
}

message(
    "Required supported samples: control=",
    required_control,
    ", HR=",
    required_HR
)

# Keep genes with >=2 isoform features.

n_iso <- tabulate(
    gene_index,
    nbins = length(genes)
)

eligible <- (
    supported &
    n_iso >= 2L
)

# ============================================================
# 6. Export coverage QC
# ============================================================

qc <- data.frame(
    gene_id = genes,
    n_transcripts = n_iso,
    gene_totals,
    check.names = FALSE
)

# Preserve original column names.

qc$all_samples_ge_threshold <- (
    n_supported == length(samp)
)

qc$eligible_for_delta_IF <- eligible

# Additional cohort diagnostics.

qc$n_samples_supported <- n_supported
qc$n_control_supported <- n_control_supported
qc$n_HR_supported <- n_HR_supported
qc$coverage_pass <- supported

write.table(
    qc,
    file.path(
        outdir,
        "isoform_gene_coverage.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. Build per-transcript output
# ============================================================

out <- data.frame(
    gene_id = x$GENEID,
    feature_id = x$TXNAME,
    x[, samp, drop = FALSE],
    as.data.frame(
        fractions,
        check.names = FALSE
    ),
    check.names = FALSE
)

out$gene_eligible <- eligible[
    gene_index
]

# ============================================================
# 8. HR/control isoform fractions
#
# Pilot:
#   requires all samples to have sufficient counts.
#
# Cohort:
#   averages samples with defined fractions.
#   Missing fractions are not replaced with zero.
# ============================================================

hr_mat <- fractions[
    ,
    match(
        paste0(hr, "_IF"),
        colnames(fractions)
    ),
    drop = FALSE
]

ctl_mat <- fractions[
    ,
    match(
        paste0(ctl, "_IF"),
        colnames(fractions)
    ),
    drop = FALSE
]

if (mode == "pilot") {

    out$HR_mean_IF <- rowMeans(hr_mat)
    out$Control_mean_IF <- rowMeans(ctl_mat)

} else {

    safe_row_mean <- function(z) {

        observed <- rowSums(
            is.finite(z)
        )

        result <- rowMeans(
            z,
            na.rm = TRUE
        )

        result[observed == 0] <- NA_real_

        result
    }

    out$HR_mean_IF <- safe_row_mean(hr_mat)

    out$Control_mean_IF <- safe_row_mean(ctl_mat)
}

out$delta_IF <- (
    out$HR_mean_IF -
    out$Control_mean_IF
)

out$HR_mean_IF[
    !out$gene_eligible
] <- NA_real_

out$Control_mean_IF[
    !out$gene_eligible
] <- NA_real_

out$delta_IF[
    !out$gene_eligible
] <- NA_real_

# ============================================================
# 9. Export supported isoforms
#
# Preserve all transcripts of eligible genes.
# ============================================================

usable <- out[
    out$gene_eligible,
    ,
    drop = FALSE
]

usable <- usable[
    order(
        usable$gene_id,
        -abs(usable$delta_IF)
    ),
    ,
    drop = FALSE
]

write.table(
    usable,
    file.path(
        outdir,
        "isoform_usage_supported.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
)

# ============================================================
# 10. Gene summary
# ============================================================

if (!nrow(usable)) {
    stop("No eligible genes after coverage filtering")
}

chunks <- split(
    usable,
    usable$gene_id
)

summary <- do.call(
    rbind,
    lapply(chunks, function(z) {

        j <- which.max(
            abs(z$delta_IF)
        )

        data.frame(
            gene_id = z$gene_id[j],
            top_transcript = z$feature_id[j],
            max_abs_delta_IF = abs(
                z$delta_IF[j]
            ),
            top_delta_IF = z$delta_IF[j],
            n_transcripts = nrow(z)
        )
    })
)

summary <- summary[
    order(
        -summary$max_abs_delta_IF
    ),
    ,
    drop = FALSE
]

write.table(
    summary,
    file.path(
        outdir,
        "isoform_usage_gene_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 11. Summary QC
# ============================================================

summary_qc <- data.frame(
    metric = c(
        "total_transcripts",
        "total_genes",
        "eligible_genes",
        "eligible_transcripts",
        "minimum_gene_counts_per_sample",
        "source",
        "mode",
        "group_fraction",
        "required_control_samples",
        "required_HR_samples"
    ),
    value = as.character(c(
        nrow(x),
        length(genes),
        sum(eligible),
        nrow(usable),
        min_gene_counts,
        basename(counts_file),
        mode,
        group_fraction,
        required_control,
        required_HR
    ))
)

write.table(
    summary_qc,
    file.path(
        outdir,
        "isoform_usage_qc.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

writeLines(
    capture.output(sessionInfo()),
    file.path(
        outdir,
        "R_session_info.txt"
    )
)

message("========================================")
message("ISOFORM_USAGE completed")
message("Mode: ", mode)
message("Eligible genes: ", sum(eligible))
message("Eligible transcripts: ", nrow(usable))
message("No statistical tests performed")
message("========================================")
