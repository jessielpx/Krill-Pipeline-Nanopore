if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org"
    )
}

if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {

    cat("\nIsoformSwitchAnalyzeR not found.\n")
    cat("Installing from Bioconductor...\n\n")

    BiocManager::install(
        "IsoformSwitchAnalyzeR",
        ask = FALSE,
        update = FALSE
    )

    if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
        stop("Failed to install IsoformSwitchAnalyzeR.")
    }
}

library(IsoformSwitchAnalyzeR)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop(
        paste0(
            "Usage:\n",
            "Rscript run_isoform_switch.R ",
            "<transcript_counts.tsv> ",
            "<transcript_tpm.tsv> ",
            "<samples.csv> ",
            "<merged_transcriptome.gtf>"
        )
    )
}

counts_file <- args[1]
tpm_file <- args[2]
design_file <- args[3]
gtf_file <- args[4]

cat("Counts file:", counts_file, "\n")
cat("TPM file:", tpm_file, "\n")
cat("Design file:", design_file, "\n")
cat("GTF file:", gtf_file, "\n\n")

for (input_file in c(
    counts_file,
    tpm_file,
    design_file,
    gtf_file
)) {
    if (!file.exists(input_file)) {
        stop("Input file does not exist: ", input_file)
    }
}

counts <- read.delim(
    counts_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

tpm <- read.delim(
    tpm_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

design <- read.csv(
    design_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]

rownames(tpm) <- tpm[[1]]
tpm <- tpm[, -1, drop = FALSE]

counts[] <- lapply(counts, as.numeric)
tpm[] <- lapply(tpm, as.numeric)

cat("Count matrix dimensions:\n")
print(dim(counts))

cat("\nTPM matrix dimensions:\n")
print(dim(tpm))

cat("\nDesign matrix:\n")
print(design)

if (!all(c("sample", "condition") %in% colnames(design))) {
    stop(
        "Design file must contain columns named 'sample' and 'condition'."
    )
}

if (!identical(colnames(counts), colnames(tpm))) {
    stop("Count and TPM sample columns are not identical.")
}

if (!identical(rownames(counts), rownames(tpm))) {
    stop("Count and TPM transcript IDs are not identical.")
}

if (!setequal(colnames(counts), design$sample)) {
    stop(
        "Sample names in the expression matrices and design file do not match."
    )
}

design <- design[
    match(colnames(counts), design$sample),
    ,
    drop = FALSE
]

cat("\nSample order after matching:\n")
print(design[, c("sample", "condition")])

design_for_isa <- design[, c("sample", "condition")]
colnames(design_for_isa)[1] <- "sampleID"

switchList <- importRdata(
    isoformCountMatrix = counts,
    isoformRepExpression = tpm,
    designMatrix = design_for_isa,
    isoformExonAnnoation = gtf_file,
    ignoreAfterPeriod = FALSE,
    showProgress = TRUE
)

cat("\nSuccessfully created IsoformSwitchAnalyzeR object.\n")
print(switchList)

cat("\nRunning differential isoform usage analysis with DEXSeq...\n")

switchList <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = switchList,
    reduceToSwitchingGenes = FALSE
)

cat("\nDEXSeq analysis completed.\n")

saveRDS(
    switchList,
    file = "isoform_switch_results.rds"
)

results_all <- switchList$isoformFeatures

write.csv(
    results_all,
    file = "isoform_switch_all.csv",
    row.names = FALSE
)

results_significant <- results_all[
    !is.na(results_all$isoform_switch_q_value) &
    results_all$isoform_switch_q_value < 0.05 &
    abs(results_all$dIF) >= 0.10,
    ,
    drop = FALSE
]

write.csv(
    results_significant,
    file = "isoform_switch_significant.csv",
    row.names = FALSE
)

summary_table <- data.frame(
    comparison = paste(
        unique(results_all$condition_1),
        "vs",
        unique(results_all$condition_2)
    ),
    total_isoforms_tested = nrow(results_all),
    total_genes_tested = length(unique(results_all$gene_ref)),
    significant_isoforms_q0.05 = sum(
        !is.na(results_all$isoform_switch_q_value) &
        results_all$isoform_switch_q_value < 0.05
    ),
    significant_isoforms_q0.05_dIF0.10 = nrow(
        results_significant
    ),
    significant_genes_q0.05_dIF0.10 = length(
        unique(results_significant$gene_ref)
    )
)

write.csv(
    summary_table,
    file = "isoform_switch_summary.csv",
    row.names = FALSE
)

cat("\nExported result files:\n")
cat("  isoform_switch_all.csv\n")
cat("  isoform_switch_significant.csv\n")
cat("  isoform_switch_summary.csv\n")

cat(
    "\nSignificant isoforms with isoform switch q-value < 0.05 ",
    "and |dIF| >= 0.10:\n"
)
cat(nrow(results_significant), "\n")

cat(
    "\nSaved IsoformSwitchAnalyzeR results to:",
    normalizePath("isoform_switch_results.rds"),
    "\n"
)
