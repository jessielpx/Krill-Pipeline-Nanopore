#!/usr/bin/env Rscript

# Krill Renew 11 - Differential Junction Usage
# Exploratory DRIMSeq-based analysis of annotated junctions.
#
# Arguments:
#   1 junction_counts.tsv
#   2 samplesheet CSV/TSV
#   3 output directory
#
# Environment:
#   KRILL_DJU_DESIGN
#   KRILL_DJU_MIN_COUNT
#   KRILL_DJU_GROUP_FRACTION

suppressPackageStartupMessages(
    library(DRIMSeq)
)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3L) {
    stop(
        "Usage: run_differential_junction_usage.R ",
        "<junction_counts.tsv> <samplesheet> <outdir>"
    )
}

counts_path <- args[1]
metadata_path <- args[2]
outdir <- args[3]

dir.create(
    outdir, recursive=TRUE,
    showWarnings=FALSE
)

min_count <- as.integer(
    Sys.getenv("KRILL_DJU_MIN_COUNT", "5")
)

fraction <- as.numeric(
    Sys.getenv("KRILL_DJU_GROUP_FRACTION", "0.25")
)

design_text <- Sys.getenv(
    "KRILL_DJU_DESIGN", "~ group"
)

if (
    is.na(min_count) || min_count < 1 ||
    is.na(fraction) || fraction <= 0 ||
    fraction > 1
) {
    stop("Invalid DJU filtering parameters")
}

read_metadata <- function(path) {
    header <- readLines(path, n=1L)

    sep <- if (grepl("\t", header)) "\t" else ","

    read.table(
        path, header=TRUE, sep=sep,
        check.names=FALSE,
        stringsAsFactors=FALSE,
        comment.char="",
        quote="\""
    )
}

write_tsv <- function(x, path) {
    write.table(
        x, file=path, sep="\t",
        quote=FALSE, row.names=FALSE,
        col.names=TRUE, na="NA"
    )
}

metadata <- read_metadata(metadata_path)

if (!all(
    c("sample", "condition") %in% names(metadata)
)) {
    stop("Samplesheet requires sample and condition")
}

metadata$sample <- trimws(metadata$sample)
metadata$condition <- trimws(metadata$condition)

if (
    anyNA(metadata$sample) ||
    any(!nzchar(metadata$sample)) ||
    anyDuplicated(metadata$sample)
) {
    stop("Invalid or duplicated sample IDs")
}

if (
    anyNA(metadata$condition) ||
    any(!nzchar(metadata$condition))
) {
    stop("Invalid condition values")
}

groups <- unique(metadata$condition)

if (length(groups) != 2L) {
    stop(
        "DJU currently requires exactly two groups"
    )
}

if (all(c("control", "HR") %in% groups)) {
    levels <- c("control", "HR")
} else {
    levels <- sort(groups)
}

metadata$group <- factor(
    metadata$condition,
    levels=levels
)

group_sizes <- table(metadata$group)

if (any(group_sizes < 2L)) {
    stop("At least two samples per group required")
}

# Preserve any additional samplesheet covariates.
samples <- metadata

names(samples)[
    names(samples) == "sample"
] <- "sample_id"

for (name in names(samples)) {
    if (
        name != "sample_id" &&
        is.character(samples[[name]])
    ) {
        samples[[name]] <- factor(
            samples[[name]]
        )
    }
}

samples$group <- metadata$group

design_formula <- as.formula(design_text)

design <- model.matrix(
    design_formula,
    data=samples
)

if (qr(design)$rank != ncol(design)) {
    stop(
        "DJU design matrix is not full rank: ",
        design_text
    )
}

group_coef <- paste0(
    "group", levels[2]
)

if (!group_coef %in% colnames(design)) {
    stop(
        "Design must contain the comparison ",
        "coefficient: ", group_coef
    )
}

junction <- read.delim(
    counts_path,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

required <- c(
    "junction_id",
    "gene_id",
    "reference_known"
)

if (!all(required %in% names(junction))) {
    stop("Missing junction matrix columns")
}

sample_ids <- as.character(samples$sample_id)

if (!all(sample_ids %in% names(junction))) {
    stop("Missing sample columns in junction matrix")
}

initial_n <- nrow(junction)

# Only uniquely assigned, reference-known junctions.
keep <- (
    !is.na(junction$gene_id) &
    nzchar(junction$gene_id) &
    junction$reference_known == "TRUE"
)

junction <- junction[
    keep,
    c(required, sample_ids),
    drop=FALSE
]

count_matrix <- as.matrix(
    junction[, sample_ids, drop=FALSE]
)

storage.mode(count_matrix) <- "numeric"

if (
    anyNA(count_matrix) ||
    any(!is.finite(count_matrix)) ||
    any(count_matrix < 0) ||
    any(abs(count_matrix - round(count_matrix)) > 1e-8)
) {
    stop(
        "Junction counts must be nonnegative integers"
    )
}

storage.mode(count_matrix) <- "integer"

# Require support in each comparison group.
# At least two replicates per group.
support <- rep(TRUE, nrow(junction))

for (group in levels) {
    members <- which(
        samples$group == group
    )

    minimum <- max(
        2L,
        ceiling(fraction * length(members))
    )

    n_supported <- rowSums(
        count_matrix[, members, drop=FALSE] >=
            min_count
    )

    support <- support &
        (n_supported >= minimum)
}

junction <- junction[
    support, , drop=FALSE
]

count_matrix <- count_matrix[
    support, , drop=FALSE
]

# DRIMSeq requires genes with >=2 retained features.
features_per_gene <- table(
    junction$gene_id
)

valid_genes <- names(
    features_per_gene[
        features_per_gene >= 2L
    ]
)

junction <- junction[
    junction$gene_id %in% valid_genes,
    ,
    drop=FALSE
]

junction <- junction[
    order(junction$gene_id, junction$junction_id),
    ,
    drop=FALSE
]

count_matrix <- as.matrix(
    junction[, sample_ids, drop=FALSE]
)
storage.mode(count_matrix) <- "integer"


qc <- data.frame(
    metric=c(
        "input_junctions",
        "eligible_reference_junctions",
        "tested_genes",
        "tested_junctions",
        "samples",
        "comparison",
        "design",
        "min_count",
        "group_fraction"
    ),
    value=as.character(c(
        initial_n,
        sum(keep),
        length(unique(junction$gene_id)),
        nrow(junction),
        nrow(samples),
        paste(levels, collapse="_vs_"),
        design_text,
        min_count,
        fraction
    ))
)

write_tsv(
    qc,
    file.path(
        outdir,
        "dju_input_qc.tsv"
    )
)

write_tsv(
    data.frame(
        gene_id=junction$gene_id,
        feature_id=junction$junction_id,
        count_matrix,
        check.names=FALSE
    ),
    file.path(
        outdir,
        "dju_filtered_counts.tsv"
    )
)

empty_gene <- data.frame(
    gene_id=character(),
    lr=numeric(),
    df=numeric(),
    pvalue=numeric(),
    adj_pvalue=numeric()
)

empty_feature <- data.frame(
    gene_id=character(),
    feature_id=character(),
    lr=numeric(),
    df=numeric(),
    pvalue=numeric(),
    adj_pvalue=numeric()
)

if (nrow(junction) == 0L) {

    write_tsv(
        empty_gene,
        file.path(outdir, "dju_gene_results.tsv")
    )

    write_tsv(
        empty_feature,
        file.path(outdir, "dju_junction_results.tsv")
    )

    message("No junctions passed filters")
    quit(save="no", status=0)
}

counts_long <- data.frame(
    gene_id=junction$gene_id,
    feature_id=junction$junction_id,
    count_matrix,
    check.names=FALSE
)

set.seed(20261110)

d <- dmDSdata(
    counts=counts_long,
    samples=samples
)

d <- dmPrecision(
    d,
    design=design
)

d <- dmFit(
    d,
    design=design
)

d <- dmTest(
    d,
    coef=group_coef
)

gene_results <- as.data.frame(
    results(d, level="gene")
)

junction_results <- as.data.frame(
    results(d, level="feature")
)

write_tsv(
    gene_results,
    file.path(outdir, "dju_gene_results.tsv")
)

write_tsv(
    junction_results,
    file.path(outdir, "dju_junction_results.tsv")
)

writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "R_session_info.txt")
)

message(
    "DJU complete: ",
    nrow(gene_results),
    " genes; ",
    nrow(junction_results),
    " junctions"
)

message(
    "Gene FDR <0.05: ",
    sum(
        gene_results$adj_pvalue < 0.05,
        na.rm=TRUE
    )
)
