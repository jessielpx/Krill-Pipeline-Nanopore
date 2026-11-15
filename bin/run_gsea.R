#!/usr/bin/env Rscript

# Krill Renew 10 - Preranked GSEA
#
# Modes:
#   prepare: cache MSigDB Ensembl-ID GMT files
#   run:     DESeq2 Wald stat -> fgseaMultilevel
#
# No gene symbol conversion.
# No significant-DEG prefiltering.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
    position <- match(name, args)
    if (is.na(position)) return(default)
    if (position == length(args)) {
        stop("Missing value for ", name)
    }
    args[position + 1L]
}

write_tsv <- function(x, path) {
    write.table(
        x, file = path, sep = "\t",
        quote = FALSE, row.names = FALSE,
        col.names = TRUE, na = "NA"
    )
}

safe_id <- function(x) {
    sub("\\.[0-9]+$", "", trimws(as.character(x)))
}

mode <- get_arg("--mode", "run")
gmt_dir <- get_arg("--gmt-dir")

if (is.null(gmt_dir)) {
    stop("Required: --gmt-dir")
}

collections <- c("Hallmark", "Reactome", "GO_BP")

gmt_names <- c(
    Hallmark = "hallmark.gmt",
    Reactome = "reactome.gmt",
    GO_BP = "go_bp.gmt"
)

if (mode == "prepare") {

    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("Install msigdbr before preparing GMT files")
    }

    dir.create(
        gmt_dir, recursive = TRUE,
        showWarnings = FALSE
    )

    get_msigdb <- function(collection, subcollection = NULL) {
        fn <- msigdbr::msigdbr

        arguments <- list(
            species = "Homo sapiens"
        )

        formal_names <- names(formals(fn))

        if ("db_species" %in% formal_names) {
            arguments$db_species <- "HS"
        }

        if ("collection" %in% formal_names) {
            arguments$collection <- collection
            if (!is.null(subcollection)) {
                arguments$subcollection <- subcollection
            }
        } else {
            arguments$category <- collection
            if (!is.null(subcollection)) {
                arguments$subcategory <- subcollection
            }
        }

        do.call(fn, arguments)
    }

    configs <- list(
        Hallmark = c("H", NA_character_),
        Reactome = c("C2", "CP:REACTOME"),
        GO_BP = c("C5", "GO:BP")
    )

    manifest <- list()

    for (label in collections) {

        cfg <- configs[[label]]

        subcollection <- if (is.na(cfg[2])) {
            NULL
        } else {
            cfg[2]
        }

        message("Fetching ", label)

        msig <- get_msigdb(cfg[1], subcollection)

        if (!all(c(
            "gs_name", "ensembl_gene"
        ) %in% names(msig))) {
            stop("MSigDB missing required fields")
        }

        msig$ensembl_gene <- safe_id(
            msig$ensembl_gene
        )

        msig <- msig[
            !is.na(msig$ensembl_gene) &
            grepl("^ENSG[0-9]+$", msig$ensembl_gene),
            ,
            drop = FALSE
        ]

        pathways <- split(
            msig$ensembl_gene,
            msig$gs_name
        )

        pathways <- lapply(pathways, unique)
        pathways <- pathways[lengths(pathways) > 0]

        if (!length(pathways)) {
            stop("No Ensembl pathways for ", label)
        }

        outfile <- file.path(
            gmt_dir, gmt_names[[label]]
        )

        con <- file(outfile, open = "wt")

        for (name in sort(names(pathways))) {
            line <- paste(
                c(name, "na", pathways[[name]]),
                collapse = "\t"
            )
            writeLines(line, con)
        }

        close(con)

        version <- if ("db_version" %in% names(msig)) {
            paste(unique(msig$db_version), collapse = ";")
        } else {
            "not_reported"
        }

        manifest[[label]] <- data.frame(
            collection = label,
            pathways = length(pathways),
            db_version = version,
            msigdbr_version = as.character(
                packageVersion("msigdbr")
            ),
            file = gmt_names[[label]]
        )

        message(
            label, ": ", length(pathways),
            " pathways"
        )
    }

    write_tsv(
        do.call(rbind, manifest),
        file.path(gmt_dir, "manifest.tsv")
    )

    message("MSigDB GMT preparation completed")
    quit(save = "no", status = 0)
}

if (mode != "run") {
    stop("Unknown mode: ", mode)
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Missing R package: fgsea")
}

dge_file <- get_arg("--dge")
outdir <- get_arg("--outdir")

if (is.null(dge_file) || is.null(outdir)) {
    stop("Required: --dge and --outdir")
}

min_size <- as.integer(get_arg("--min-size", "10"))
max_size <- as.integer(get_arg("--max-size", "500"))

if (
    is.na(min_size) || is.na(max_size) ||
    min_size < 1 || max_size < min_size
) {
    stop("Invalid gene set size limits")
}

dir.create(
    outdir, recursive = TRUE,
    showWarnings = FALSE
)

dge <- read.delim(
    dge_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required <- c("GENEID", "stat")
missing <- setdiff(required, names(dge))

if (length(missing)) {
    stop(
        "DGE missing: ",
        paste(missing, collapse = ", ")
    )
}

dge$GENEID <- safe_id(dge$GENEID)
dge$stat <- suppressWarnings(
    as.numeric(dge$stat)
)

valid <- (
    grepl("^ENSG[0-9]+$", dge$GENEID) &
    is.finite(dge$stat)
)

rank_table <- dge[valid, c("GENEID", "stat")]

if (!nrow(rank_table)) {
    stop("No valid Ensembl-ranked genes")
}

# Deterministic handling of duplicate Ensembl IDs:
# keep the transcript/gene row with largest |stat|.
rank_table <- rank_table[
    order(
        -abs(rank_table$stat),
        rank_table$GENEID
    ),
    ,
    drop = FALSE
]

rank_table <- rank_table[
    !duplicated(rank_table$GENEID),
    ,
    drop = FALSE
]

ranks <- rank_table$stat
names(ranks) <- rank_table$GENEID

ranks <- sort(ranks, decreasing = TRUE)

if (length(ranks) < 10) {
    stop("Fewer than 10 ranked genes")
}

write_tsv(
    data.frame(
        GENEID = names(ranks),
        stat = unname(ranks)
    ),
    file.path(outdir, "gsea_ranked_genes.tsv")
)

result_columns <- c(
    "collection", "pathway", "pval", "padj",
    "ES", "NES", "size", "log2err", "leadingEdge"
)

empty_result <- function() {
    data.frame(
        collection = character(),
        pathway = character(),
        pval = numeric(),
        padj = numeric(),
        ES = numeric(),
        NES = numeric(),
        size = integer(),
        log2err = numeric(),
        leadingEdge = character()
    )
}

all_results <- list()
qc <- list()
universe <- names(ranks)

for (label in collections) {

    gmt_file <- file.path(
        gmt_dir, gmt_names[[label]]
    )

    if (!file.exists(gmt_file)) {
        stop("Missing GMT file: ", gmt_file)
    }

    pathways <- fgsea::gmtPathways(gmt_file)
    pathways <- lapply(pathways, unique)

    all_genes <- unique(
        unlist(pathways, use.names = FALSE)
    )

    overlap <- intersect(universe, all_genes)

    # Explicit intersection guarantees that pathway
    # sizes are evaluated in the ranked-gene universe.
    restricted <- lapply(
        pathways, intersect, y = universe
    )

    sizes <- lengths(restricted)

    eligible <- restricted[
        sizes >= min_size & sizes <= max_size
    ]

    result <- empty_result()

    if (length(eligible) > 0) {

        set.seed(20261110)

        fg <- fgsea::fgseaMultilevel(
            pathways = eligible,
            stats = ranks,
            minSize = min_size,
            maxSize = max_size,
            eps = 0,
            nproc = 1
        )

        fg <- as.data.frame(fg)

        if (nrow(fg)) {

            fg$leadingEdge <- vapply(
                fg$leadingEdge,
                function(x) paste(x, collapse = ";"),
                character(1)
            )

            fg$collection <- label

            result <- fg[, result_columns]

            result <- result[
                order(
                    result$padj,
                    -abs(result$NES),
                    na.last = TRUE
                ),
                ,
                drop = FALSE
            ]
        }
    }

    all_results[[label]] <- result

    qc[[label]] <- data.frame(
        collection = label,
        input_dge_rows = nrow(dge),
        ranked_genes = length(ranks),
        pathways_total = length(pathways),
        pathways_eligible = length(eligible),
        genes_overlapping_msigdb = length(overlap),
        pathways_tested = nrow(result),
        pathways_fdr_0_05 = sum(
            result$padj < 0.05,
            na.rm = TRUE
        ),
        min_size = min_size,
        max_size = max_size
    )

    write_tsv(
        result,
        file.path(
            outdir,
            paste0(
                "gsea_",
                tolower(label),
                ".tsv"
            )
        )
    )

    message(
        label, ": ", nrow(result),
        " tested, ",
        sum(result$padj < 0.05, na.rm = TRUE),
        " FDR < 0.05"
    )
}

combined <- do.call(
    rbind, all_results
)

# Preserve an empty but valid table if
# no pathways pass size filters.
if (is.null(combined)) {
    combined <- empty_result()
}

write_tsv(
    combined,
    file.path(outdir, "gsea_all_results.tsv")
)

write_tsv(
    do.call(rbind, qc),
    file.path(outdir, "gsea_qc.tsv")
)

writeLines(
    capture.output(sessionInfo()),
    file.path(outdir, "R_session_info.txt")
)

message("Krill GSEA completed")
