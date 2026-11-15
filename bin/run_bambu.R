# Krill Renew - Joint bambu transcriptome analysis

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5L) {
    stop(
        "Usage: run_bambu.R ",
        "<manifest.tsv> <annotation.gtf> ",
        "<genome.fa> <output_dir> <cores>"
    )
}

manifest_file <- args[[1]]
annotation_file <- args[[2]]
genome_file <- args[[3]]
output_dir <- args[[4]]
cores <- as.integer(args[[5]])

if (is.na(cores) || cores < 1L) {
    stop("Invalid CPU count")
}

for (f in c(manifest_file, annotation_file, genome_file)) {
    if (!file.exists(f)) {
        stop("Missing input file: ", f)
    }
}

# Load required genomic functions before starting workers.
suppressPackageStartupMessages({
    library(Seqinfo)
    library(GenomeInfoDb)
    library(BiocParallel)
})

for (pkg in c("bambu", "Rsamtools", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required R package missing: ", pkg)
    }
}

# Read and validate sample manifest.
manifest <- read.delim(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!all(c("sample", "bam") %in% names(manifest))) {
    stop("Manifest requires sample and bam columns")
}

if (nrow(manifest) < 1L) {
    stop("No samples found")
}

if (
    anyNA(manifest$sample) ||
    anyNA(manifest$bam) ||
    any(!nzchar(manifest$sample)) ||
    any(!nzchar(manifest$bam)) ||
    anyDuplicated(manifest$sample)
) {
    stop("Invalid or duplicated sample IDs/BAM paths")
}

manifest <- manifest[
    order(manifest$sample),
    ,
    drop = FALSE
]

bam_paths <- manifest$bam

for (bam in bam_paths) {
    if (!file.exists(paste0(bam, ".bai"))) {
        stop("Missing BAM index: ", bam, ".bai")
    }
}

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# Explicitly configure BiocParallel on Rorqual.
# Use no more workers than available samples.
workers <- min(cores, nrow(manifest))

BiocParallel::register(
    BiocParallel::MulticoreParam(
        workers = workers,
        timeout = 86400,
        stop.on.error = TRUE
    ),
    default = TRUE
)

cat("===== KRILL BAMBU =====\n")
cat("Samples:", nrow(manifest), "\n")
cat("Workers:", workers, "\n")
cat("bambu version:", as.character(
    packageVersion("bambu")
), "\n")
cat("BiocParallel timeout: 86400 seconds\n")

print(manifest)

# Build local FASTA index if it is not staged.
# Passing FaFile avoids an online BSgenome repository lookup.
fai_file <- paste0(genome_file, ".fai")

if (!file.exists(fai_file)) {
    message("Creating local FASTA index")
    Rsamtools::indexFa(genome_file)
}

genome <- Rsamtools::FaFile(genome_file)

if (length(Seqinfo::seqlengths(genome)) == 0L) {
    stop("Could not read reference genome sequence lengths")
}

# Prepare reference annotation.
message("--- Preparing bambu annotations ---")

annotation <- bambu::prepareAnnotations(
    annotation_file
)

# Create indexed BAM collection.
bam_files <- Rsamtools::BamFileList(bam_paths)
names(bam_files) <- manifest$sample

# Joint discovery and quantification.
message("--- Starting bambu joint analysis ---")

se <- bambu::bambu(
    reads = bam_files,
    annotations = annotation,
    genome = genome,
    discovery = TRUE,
    quant = TRUE,
    ncore = workers,
    stranded = FALSE
)

if (ncol(se) != nrow(manifest)) {
    stop("Bambu sample count does not match manifest")
}

colnames(se) <- manifest$sample

# Save complete SummarizedExperiment.
saveRDS(
    se,
    file.path(output_dir, "bambu_transcripts.rds")
)

# Export extended annotations and count matrices.
bambu::writeBambuOutput(
    se,
    path = output_dir
)

# Per-sample basic QC.
counts <- SummarizedExperiment::assay(
    se,
    "counts"
)

sample_qc <- data.frame(
    sample = colnames(se),
    reference_transcripts = nrow(se),
    detected_transcripts = colSums(counts > 0),
    transcript_estimated_counts = colSums(counts),
    check.names = FALSE
)

write.table(
    sample_qc,
    file.path(output_dir, "bambu_sample_qc.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    manifest,
    file.path(output_dir, "sample_manifest.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

capture.output(
    sessionInfo(),
    file = file.path(output_dir, "R_session_info.txt")
)

message("--- Finished KRILL BAMBU ---")
