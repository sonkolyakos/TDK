#!/usr/bin/env Rscript

###############################################################################
# Multi-database ORA of shared (common) splicing genes -> Figure 11
#
# Databases (MSigDB):
#   CANONICAL_ALL = C2:CP  Reactome + WikiPathways + KEGG Medicus + BioCarta + PID
#   GO_BP         = C5:GO:BP
#   HALLMARK      = H
#
# Usage:
#   Rscript ora_multidb.R common_genes_MCF7_T47D_core_pan.csv ORA_out tested_universe.tsv
#
# Universe: all genes testable for splicing in BOTH cell lines.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(clusterProfiler)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)

# ---- Arguments --------------------------------------------------------------
ARGS <- commandArgs(trailingOnly = TRUE)
INPUT <- if (length(ARGS) >= 1) ARGS[1] else
  "/home/akos/Desktop/lab/TDK/common_genes_MCF7_T47D_core_pan.csv"
OUT <- if (length(ARGS) >= 2) ARGS[2] else "ORA_common_genes_multiDB"
UNIVERSE_FILE <- if (length(ARGS) >= 3 && nzchar(ARGS[3])) ARGS[3] else NA_character_

FIGDIR <- file.path(OUT, "figures")
TABDIR <- file.path(OUT, "tables")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Parameters -------------------------------------------------------------
ORA_FDR         <- 0.05
MIN_GS_SIZE     <- 15
MAX_GS_SIZE     <- 500
TOP_N_FACET     <- 10
TEXT_WRAP_WIDTH <- 40

# ---- Helpers ----------------------------------------------------------------
pick_column <- function(x, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) {
    lower_names <- tolower(names(x))
    for (candidate in candidates) {
      idx <- which(lower_names == tolower(candidate))
      if (length(idx)) return(names(x)[idx[1]])
    }
    if (required) {
      stop("Could not identify required column. Candidates: ",
           paste(candidates, collapse = ", "),
           "\nAvailable columns: ", paste(names(x), collapse = ", "))
    }
    return(NA_character_)
  }
  hit[1]
}

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", ".")] <- NA_character_
  unique(x[!is.na(x)])
}

ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2) return(NA_real_)
    a <- suppressWarnings(as.numeric(z[1]))
    b <- suppressWarnings(as.numeric(z[2]))
    if (!is.finite(a) || !is.finite(b) || b == 0) return(NA_real_)
    a / b
  }, numeric(1))
}

get_msig <- function(collection, subcollection = NULL) {
  fml <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% fml) {
    args <- list(species = "Homo sapiens", collection = collection)
    if (!is.null(subcollection) && nzchar(subcollection)) args$subcollection <- subcollection
  } else {  # older msigdbr releases
    args <- list(species = "Homo sapiens", category = collection)
    if (!is.null(subcollection) && nzchar(subcollection)) args$subcategory <- subcollection
  }
  X <- as.data.table(do.call(msigdbr::msigdbr, args))
  gene_col <- pick_column(X, c("gene_symbol", "human_gene_symbol"))
  set_col  <- pick_column(X, c("gs_name"))
  X[, gene := as.character(get(gene_col))]
  X[, term := as.character(get(set_col))]
  unique(X[!is.na(gene) & nzchar(gene) & !is.na(term) & nzchar(term), .(term, gene)])
}

safe_get_msig <- function(collection, subcollection = NULL, label) {
  cat("\nLoading MSigDB: ", label, "\n", sep = "")
  out <- tryCatch(get_msig(collection, subcollection), error = function(e) {
    warning("Could not load ", label, ": ", conditionMessage(e))
    data.table(term = character(), gene = character())
  })
  cat("  gene sets = ", uniqueN(out$term), "; mappings = ", nrow(out), "\n", sep = "")
  out
}

run_ora <- function(foreground, universe, term2gene, database) {
  term2gene <- unique(term2gene[!is.na(term) & !is.na(gene)])
  if (!nrow(term2gene)) return(data.table())

  db_genes <- unique(term2gene$gene)
  fg <- intersect(foreground, db_genes)
  bg <- intersect(universe, db_genes)
  fg <- intersect(fg, bg)   # foreground must be inside the background

  qc <- data.table(Database = database,
                   Foreground_total = length(foreground),
                   Foreground_in_database = length(fg),
                   Universe_total = length(universe),
                   Universe_in_database = length(bg),
                   Gene_sets = uniqueN(term2gene$term))
  fwrite(qc, file.path(TABDIR, paste0(database, "__input_QC.tsv")), sep = "\t")

  if (length(fg) < 2 || length(bg) < 10) {
    warning(database, ": insufficient genes after universe/database intersection.")
    return(data.table())
  }

  E <- clusterProfiler::enricher(
    gene = fg, universe = bg,
    TERM2GENE = term2gene[, .(term, gene)],
    pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
    minGSSize = MIN_GS_SIZE, maxGSSize = MAX_GS_SIZE
  )
  if (is.null(E) || nrow(as.data.frame(E)) == 0) return(data.table())

  R <- as.data.table(as.data.frame(E))
  R[, Database := database]
  R[, GeneRatio_numeric := ratio_to_numeric(GeneRatio)]
  R[, minus_log10_FDR := -log10(pmax(p.adjust, 1e-300))]
  setorder(R, p.adjust, pvalue, -Count)
  R
}

# ---- Foreground -------------------------------------------------------------
FG_TABLE <- fread(INPUT)
FG_COL <- pick_column(FG_TABLE, c("Gene_Symbol", "gene", "gene_symbol", "SYMBOL", "symbol", "Gene"))
FOREGROUND <- clean_symbols(FG_TABLE[[FG_COL]])
cat("Unique foreground symbols: ", length(FOREGROUND), "\n", sep = "")

# ---- Universe ---------------------------------------------------------------
HAS_REAL_UNIVERSE <- FALSE
if (!is.na(UNIVERSE_FILE) && file.exists(UNIVERSE_FILE)) {
  U <- fread(UNIVERSE_FILE)
  U_COL <- pick_column(U, c("Gene_Symbol", "gene", "gene_symbol", "SYMBOL", "symbol", "Gene"))
  cell_col <- pick_column(U, c("cellline", "cell_line", "CellLine"), required = FALSE)

  if (!is.na(cell_col)) {
    # long table (one cell line per row): keep genes testable in >= 2 cell lines
    TMP <- unique(U[, .(gene = as.character(get(U_COL)), cellline = as.character(get(cell_col)))])
    TMP <- TMP[!is.na(gene) & nzchar(gene)]
    UNIVERSE <- clean_symbols(TMP[, .(n = uniqueN(cellline)), by = gene][n >= 2, gene])
  } else {
    UNIVERSE <- clean_symbols(U[[U_COL]])
  }
  HAS_REAL_UNIVERSE <- TRUE
  cat("Using explicit tested universe: ", length(UNIVERSE), " genes\n", sep = "")
} else {
  UNIVERSE <- NULL
  cat("WARNING: no tested-gene universe supplied; falling back to MSigDB-annotated genes.\n")
}

# ---- MSigDB databases -------------------------------------------------------
DB <- list()
DB[["HALLMARK"]] <- safe_get_msig("H", NULL, "Hallmark")
DB[["GO_BP"]]    <- safe_get_msig("C5", "GO:BP", "GO Biological Process")

canonical_src <- list(
  REACTOME     = safe_get_msig("C2", "CP:REACTOME",     "Reactome"),
  WIKIPATHWAYS = safe_get_msig("C2", "CP:WIKIPATHWAYS", "WikiPathways"),
  KEGG_MEDICUS = safe_get_msig("C2", "CP:KEGG_MEDICUS", "KEGG Medicus"),
  BIOCARTA     = safe_get_msig("C2", "CP:BIOCARTA",     "BioCarta"),
  PID          = safe_get_msig("C2", "CP:PID",          "PID")
)

# Source is kept in the term name so identical names from different sources cannot collide
DB[["CANONICAL_ALL"]] <- unique(rbindlist(lapply(names(canonical_src), function(nm) {
  X <- copy(canonical_src[[nm]])
  if (nrow(X)) X[, term := paste0(nm, "::", term)]
  X
}), fill = TRUE))

if (!HAS_REAL_UNIVERSE) {
  UNIVERSE <- unique(unlist(lapply(DB[c("HALLMARK", "CANONICAL_ALL", "GO_BP")],
                                   function(x) x$gene), use.names = FALSE))
}

FOREGROUND_IN_UNIVERSE <- intersect(FOREGROUND, UNIVERSE)
cat("Foreground genes in universe: ", length(FOREGROUND_IN_UNIVERSE), " / ",
    length(FOREGROUND), "\n", sep = "")

# ---- ORA per database -------------------------------------------------------
RUN_ORDER <- c("CANONICAL_ALL", "GO_BP", "HALLMARK")
RESULTS <- list()

for (name in RUN_ORDER) {
  cat("\nORA: ", name, "\n", sep = "")
  R <- run_ora(FOREGROUND_IN_UNIVERSE, UNIVERSE, DB[[name]], name)
  RESULTS[[name]] <- R
  if (nrow(R)) {
    fwrite(R, file.path(TABDIR, paste0(name, "__ORA_all.tsv")), sep = "\t")
    fwrite(R[p.adjust < ORA_FDR], file.path(TABDIR, paste0(name, "__ORA_FDR05.tsv")), sep = "\t")
    cat("Tested terms: ", nrow(R), "; FDR < 0.05: ", sum(R$p.adjust < ORA_FDR), "\n", sep = "")
  }
}

# ---- Figure 11: faceted multi-database summary ------------------------------
FACET <- rbindlist(lapply(RUN_ORDER, function(name) {
  R <- RESULTS[[name]]
  if (is.null(R) || !nrow(R)) return(NULL)
  X <- R[p.adjust < ORA_FDR]
  if (!nrow(X)) return(NULL)
  head(X, TOP_N_FACET)
}), fill = TRUE)

if (nrow(FACET)) {
  # clean term names and wrap long labels
  FACET[, clean_term := sub("^.*::", "", Description)]
  FACET[, clean_term := sub("^[A-Z]+_", "", clean_term)]
  FACET[, clean_term := gsub("_", " ", clean_term)]
  FACET[, clean_term := vapply(clean_term, function(x)
    paste(strwrap(x, width = TEXT_WRAP_WIDTH), collapse = "\n"), character(1))]

  setorder(FACET, Database, GeneRatio_numeric)
  FACET[, plot_label := factor(clean_term, levels = unique(clean_term))]

  p <- ggplot(FACET, aes(GeneRatio_numeric, plot_label)) +
    geom_point(aes(size = Count, color = minus_log10_FDR), alpha = 0.9) +
    facet_grid(Database ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "grey70", high = "black") +
    labs(title = "Shared splicing genes: multi-database ORA",
         x = "Gene ratio", y = NULL, size = "Genes", color = "-log10(FDR)") +
    theme_classic(base_size = 14) +
    theme(strip.text.y = element_text(face = "bold", angle = 0, size = 12),
          axis.text.y  = element_text(size = 11),
          plot.title   = element_text(face = "bold", size = 16))

  h <- max(6, 0.45 * nrow(FACET) + 2)
  ggsave(file.path(FIGDIR, "MULTIDB_ORA_top_significant.png"), p,
         width = 10, height = h, dpi = 300, bg = "white")
  ggsave(file.path(FIGDIR, "MULTIDB_ORA_top_significant.pdf"), p,
         width = 10, height = h, bg = "white")
}
