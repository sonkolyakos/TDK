#!/usr/bin/env Rscript

###############################################################################
# MULTI-DATABASE ORA FOR SHARED / COMMON SPLICING GENES
#
# Primary pathway analysis:
#   MSigDB C2:CP canonical pathways
#     - Reactome
#     - WikiPathways
#     - KEGG Medicus
#     - BioCarta
#     - PID
#     - KEGG Legacy (optional historical)
#
# Complementary:
#   H        = Hallmark
#   C5:GO:BP = Gene Ontology Biological Process
#   C6       = Oncogenic signatures (optional, cancer-focused)
#
# IMPORTANT:
# For a statistically correct splicing ORA, the preferred universe is:
#   ALL GENES THAT WERE TESTABLE FOR SPLICING IN BOTH MCF7 AND T47D
#
# Usage:
#
# Rscript 24_common_genes_multiDB_ORA.R \
#   common_genes_MCF7_T47D_core_pan.csv \
#   ORA_common_genes_multiDB \
#   tested_splicing_universe.tsv
#
# If the 3rd argument is omitted, the script will run using each MSigDB
# collection's annotated genes as the fallback universe, but it will print a
# strong warning because that is less appropriate for the final TDK.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(clusterProfiler)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)

###############################################################################
# ARGUMENTS
###############################################################################

ARGS <- commandArgs(trailingOnly = TRUE)
INPUT <- if (length(ARGS) >= 1) {
  ARGS[1]
} else {
  "/home/akos/Desktop/lab/TDK/common_genes_MCF7_T47D_core_pan.csv"
}

OUT <- if (length(ARGS) >= 2) {
  ARGS[2]
} else {
  "ORA_common_genes_multiDB"
}

UNIVERSE_FILE <- if (length(ARGS) >= 3 && nzchar(ARGS[3])) {
  ARGS[3]
} else {
  NA_character_
}

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

FIGDIR <- file.path(
  OUT,
  "figures"
)

TABDIR <- file.path(
  OUT,
  "tables"
)

dir.create(
  FIGDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# PARAMETERS
###############################################################################

ORA_FDR <- 0.05

MIN_GS_SIZE <- 15
MAX_GS_SIZE <- 500

TOP_N_PLOT <- 25
TOP_N_FACET <- 10

# Presentation-only redundancy filter based on overlap among foreground hits.
# This does NOT change statistical testing.
JACCARD_THRESHOLD <- 0.55

RUN_ONCOGENIC <- TRUE
RUN_KEGG_LEGACY <- FALSE

# Hány karakternél törjenek a sorok a plotokon?
TEXT_WRAP_WIDTH <- 40 

###############################################################################
# HELPERS
###############################################################################

pick_column <- function(
    x,
    candidates,
    required = TRUE
) {
  
  hit <- candidates[
    candidates %in% names(x)
  ]
  
  if (!length(hit)) {
    
    # case-insensitive fallback
    lower_names <- tolower(names(x))
    
    for (candidate in candidates) {
      
      idx <- which(
        lower_names ==
          tolower(candidate)
      )
      
      if (length(idx)) {
        return(
          names(x)[idx[1]]
        )
      }
    }
    
    if (required) {
      stop(
        "Could not identify required column. Candidates: ",
        paste(
          candidates,
          collapse = ", "
        ),
        "\nAvailable columns: ",
        paste(
          names(x),
          collapse = ", "
        )
      )
    }
    
    return(NA_character_)
  }
  
  hit[1]
}


clean_symbols <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "NULL",
      "."
    )
  ] <- NA_character_
  
  unique(
    x[
      !is.na(x)
    ]
  )
}


ratio_to_numeric <- function(x) {
  
  vapply(
    strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    ),
    function(z) {
      
      if (length(z) != 2) {
        return(NA_real_)
      }
      
      a <- suppressWarnings(
        as.numeric(z[1])
      )
      
      b <- suppressWarnings(
        as.numeric(z[2])
      )
      
      if (
        !is.finite(a) ||
        !is.finite(b) ||
        b == 0
      ) {
        return(NA_real_)
      }
      
      a / b
    },
    numeric(1)
  )
}


get_msig <- function(
    collection,
    subcollection = NULL
) {
  
  fml <- names(
    formals(
      msigdbr::msigdbr
    )
  )
  
  if ("collection" %in% fml) {
    
    args <- list(
      species = "Homo sapiens",
      collection = collection
    )
    
    if (
      !is.null(subcollection) &&
      nzchar(subcollection)
    ) {
      args$subcollection <- subcollection
    }
    
  } else {
    
    # compatibility with older msigdbr releases
    args <- list(
      species = "Homo sapiens",
      category = collection
    )
    
    if (
      !is.null(subcollection) &&
      nzchar(subcollection)
    ) {
      args$subcategory <- subcollection
    }
  }
  
  X <- do.call(
    msigdbr::msigdbr,
    args
  )
  
  X <- as.data.table(
    X
  )
  
  gene_col <- pick_column(
    X,
    c(
      "gene_symbol",
      "human_gene_symbol"
    )
  )
  
  set_col <- pick_column(
    X,
    c(
      "gs_name"
    )
  )
  
  X[
    ,
    gene :=
      as.character(
        get(
          gene_col
        )
      )
  ]
  
  X[
    ,
    term :=
      as.character(
        get(
          set_col
        )
      )
  ]
  
  unique(
    X[
      !is.na(gene) &
        nzchar(gene) &
        !is.na(term) &
        nzchar(term),
      .(
        term,
        gene
      )
    ]
  )
}


safe_get_msig <- function(
    collection,
    subcollection = NULL,
    label
) {
  
  cat(
    "\nLoading MSigDB: ",
    label,
    " [",
    collection,
    ifelse(
      is.null(subcollection),
      "",
      paste0(
        " / ",
        subcollection
      )
    ),
    "]\n",
    sep = ""
  )
  
  out <- tryCatch(
    get_msig(
      collection,
      subcollection
    ),
    error = function(e) {
      
      warning(
        "Could not load ",
        label,
        ": ",
        conditionMessage(e)
      )
      
      data.table(
        term = character(),
        gene = character()
      )
    }
  )
  
  cat(
    "  gene sets = ",
    uniqueN(out$term),
    "; mappings = ",
    nrow(out),
    "\n",
    sep = ""
  )
  
  out
}


run_ora <- function(
    foreground,
    universe,
    term2gene,
    database
) {
  
  term2gene <- unique(
    term2gene[
      !is.na(term) &
        !is.na(gene)
    ]
  )
  
  if (!nrow(term2gene)) {
    return(
      data.table()
    )
  }
  
  db_genes <- unique(
    term2gene$gene
  )
  
  fg <- intersect(
    foreground,
    db_genes
  )
  
  bg <- intersect(
    universe,
    db_genes
  )
  
  # Foreground must be inside background.
  fg <- intersect(
    fg,
    bg
  )
  
  qc <- data.table(
    Database = database,
    Foreground_total = length(foreground),
    Foreground_in_database = length(fg),
    Universe_total = length(universe),
    Universe_in_database = length(bg),
    Gene_sets = uniqueN(term2gene$term)
  )
  
  fwrite(
    qc,
    file.path(
      TABDIR,
      paste0(
        database,
        "__input_QC.tsv"
      )
    ),
    sep = "\t"
  )
  
  if (
    length(fg) < 2 ||
    length(bg) < 10
  ) {
    
    warning(
      database,
      ": insufficient genes after universe/database intersection."
    )
    
    return(
      data.table()
    )
  }
  
  E <- clusterProfiler::enricher(
    gene = fg,
    universe = bg,
    TERM2GENE = term2gene[
      ,
      .(
        term,
        gene
      )
    ],
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = MIN_GS_SIZE,
    maxGSSize = MAX_GS_SIZE
  )
  
  if (
    is.null(E) ||
    nrow(
      as.data.frame(E)
    ) == 0
  ) {
    return(
      data.table()
    )
  }
  
  R <- as.data.table(
    as.data.frame(E)
  )
  
  R[
    ,
    Database := database
  ]
  
  R[
    ,
    GeneRatio_numeric :=
      ratio_to_numeric(
        GeneRatio
      )
  ]
  
  R[
    ,
    BgRatio_numeric :=
      ratio_to_numeric(
        BgRatio
      )
  ]
  
  R[
    ,
    minus_log10_FDR :=
      -log10(
        pmax(
          p.adjust,
          1e-300
        )
      )
  ]
  
  setorder(
    R,
    p.adjust,
    pvalue,
    -Count
  )
  
  R
}


reduce_redundancy <- function(
    R,
    threshold = JACCARD_THRESHOLD,
    max_terms = 100
) {
  
  if (!nrow(R)) {
    return(
      R
    )
  }
  
  X <- copy(R)
  
  X <- X[
    is.finite(p.adjust)
  ]
  
  setorder(
    X,
    p.adjust,
    pvalue,
    -Count
  )
  
  selected <- integer()
  selected_genes <- list()
  
  for (i in seq_len(nrow(X))) {
    
    gi <- unique(
      strsplit(
        as.character(
          X$geneID[i]
        ),
        "/",
        fixed = TRUE
      )[[1]]
    )
    
    keep <- TRUE
    
    if (length(selected_genes)) {
      
      sims <- vapply(
        selected_genes,
        function(gj) {
          
          denom <- length(
            union(
              gi,
              gj
            )
          )
          
          if (!denom) {
            return(0)
          }
          
          length(
            intersect(
              gi,
              gj
            )
          ) / denom
        },
        numeric(1)
      )
      
      if (
        any(
          sims >= threshold
        )
      ) {
        keep <- FALSE
      }
    }
    
    if (keep) {
      
      selected <- c(
        selected,
        i
      )
      
      selected_genes[[length(selected_genes) + 1]] <- list(
        gi
      )
    }
    
    if (
      length(selected) >=
      max_terms
    ) {
      break
    }
  }
  
  X[
    selected
  ]
}


plot_ora <- function(
    R,
    database,
    filename,
    top_n = TOP_N_PLOT
) {
  
  if (!nrow(R)) {
    return(
      invisible(NULL)
    )
  }
  
  PLOT <- copy(
    R[
      p.adjust < ORA_FDR
    ]
  )
  
  if (!nrow(PLOT)) {
    
    PLOT <- copy(
      head(
        R,
        top_n
      )
    )
    
    subtitle <- paste0(
      "No terms at FDR < ",
      ORA_FDR,
      "; showing best nominal results"
    )
    
  } else {
    
    PLOT <- head(
      PLOT,
      top_n
    )
    
    subtitle <- paste0(
      "Significant terms: BH FDR < ",
      ORA_FDR
    )
  }
  
  # Sortörés alkalmazása a címkéken
  PLOT[
    ,
    clean_desc := vapply(
      as.character(Description), 
      function(x) { paste(strwrap(x, width = TEXT_WRAP_WIDTH), collapse = "\n") }, 
      character(1)
    )
  ]
  
  PLOT[
    ,
    clean_desc :=
      factor(
        clean_desc,
        levels = rev(
          unique(clean_desc)
        )
      )
  ]
  
  p <- ggplot(
    PLOT,
    aes(
      x = GeneRatio_numeric,
      y = clean_desc
    )
  ) +
    
    geom_point(
      aes(
        size = Count,
        color = minus_log10_FDR
      ),
      alpha = 0.9
    ) +
    
    scale_color_gradient(
      low = "grey70",
      high = "black"
    ) +
    
    labs(
      title = paste0(
        database,
        " ORA"
      ),
      subtitle = subtitle,
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = "-log10(FDR)"
    ) +
    
    theme_classic(
      base_size = 14 # Megnövelt betűméret
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 16 # Megnövelt cím méret
        ),
      axis.text.y =
        element_text(
          size = 12 # Megnövelt Y tengely betűméret
        )
    )
  
  ggsave(
    file.path(
      FIGDIR,
      paste0(
        filename,
        ".png"
      )
    ),
    p,
    width = 10,
    height = max(
      6,
      0.45 *
        nrow(PLOT) +
        2
    ),
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(
      FIGDIR,
      paste0(
        filename,
        ".pdf"
      )
    ),
    p,
    width = 10,
    height = max(
      6,
      0.45 *
        nrow(PLOT) +
        2
    ),
    bg = "white"
  )
}


###############################################################################
# READ FOREGROUND
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "READING FOREGROUND GENES\n"
)

cat(
  "============================================================\n"
)

FG_TABLE <- fread(
  INPUT
)

FG_COL <- pick_column(
  FG_TABLE,
  c(
    "Gene_Symbol",
    "gene",
    "gene_symbol",
    "SYMBOL",
    "symbol",
    "Gene"
  )
)

FOREGROUND <- clean_symbols(
  FG_TABLE[[FG_COL]]
)

cat(
  "Input file: ",
  INPUT,
  "\n",
  sep = ""
)

cat(
  "Gene column: ",
  FG_COL,
  "\n",
  sep = ""
)

cat(
  "Unique foreground symbols: ",
  length(FOREGROUND),
  "\n",
  sep = ""
)

if (
  "In_Tested_Universe" %in%
  names(FG_TABLE)
) {
  
  cat(
    "In_Tested_Universe distribution:\n"
  )
  
  print(
    table(
      FG_TABLE$In_Tested_Universe,
      useNA = "ifany"
    )
  )
}

fwrite(
  data.table(
    gene = FOREGROUND
  ),
  file.path(
    TABDIR,
    "foreground_genes_used.tsv"
  ),
  sep = "\t"
)

###############################################################################
# READ / DEFINE UNIVERSE
###############################################################################

HAS_REAL_UNIVERSE <- FALSE

if (
  !is.na(UNIVERSE_FILE) &&
  file.exists(UNIVERSE_FILE)
) {
  
  U <- fread(
    UNIVERSE_FILE
  )
  
  U_COL <- pick_column(
    U,
    c(
      "Gene_Symbol",
      "gene",
      "gene_symbol",
      "SYMBOL",
      "symbol",
      "Gene"
    )
  )
  
  # If a long table contains one cell line per row, derive genes present
  # in both MCF7 and T47D.
  cell_col <- pick_column(
    U,
    c(
      "cellline",
      "cell_line",
      "CellLine"
    ),
    required = FALSE
  )
  
  if (!is.na(cell_col)) {
    
    TMP <- unique(
      U[
        ,
        .(
          gene =
            as.character(
              get(
                U_COL
              )
            ),
          cellline =
            as.character(
              get(
                cell_col
              )
            )
        )
      ]
    )
    
    TMP <- TMP[
      !is.na(gene) &
        nzchar(gene)
    ]
    
    COMMON_UNIVERSE <- TMP[
      ,
      .(
        n_celllines =
          uniqueN(
            cellline
          )
      ),
      by = gene
    ][
      n_celllines >= 2,
      gene
    ]
    
    UNIVERSE <- clean_symbols(
      COMMON_UNIVERSE
    )
    
    cat(
      "\nUniverse file contains cell-line labels.\n"
    )
    
    cat(
      "Using genes present/testable in >=2 cell lines.\n"
    )
    
  } else {
    
    UNIVERSE <- clean_symbols(
      U[[U_COL]]
    )
  }
  
  HAS_REAL_UNIVERSE <- TRUE
  
  cat(
    "\nUsing explicit tested universe: ",
    UNIVERSE_FILE,
    "\n",
    sep = ""
  )
  
  cat(
    "Universe genes: ",
    length(UNIVERSE),
    "\n",
    sep = ""
  )
  
} else {
  
  UNIVERSE <- NULL
  
  cat(
    "\n************************************************************\n"
  )
  
  cat(
    "WARNING: NO TESTED-GENE UNIVERSE WAS PROVIDED.\n"
  )
  
  cat(
    "The script will fall back to genes annotated in each MSigDB\n"
  )
  
  cat(
    "collection. This is acceptable for exploration, but for the\n"
  )
  
  cat(
    "final TDK ORA you should use the genes testable for splicing\n"
  )
  
  cat(
    "in BOTH MCF7 and T47D as the background universe.\n"
  )
  
  cat(
    "************************************************************\n\n"
  )
}

###############################################################################
# MSIGDB DATABASES
###############################################################################

DB <- list()

DB[["HALLMARK"]] <- safe_get_msig(
  "H",
  NULL,
  "Hallmark"
)

DB[["REACTOME"]] <- safe_get_msig(
  "C2",
  "CP:REACTOME",
  "Reactome"
)

DB[["WIKIPATHWAYS"]] <- safe_get_msig(
  "C2",
  "CP:WIKIPATHWAYS",
  "WikiPathways"
)

DB[["KEGG_MEDICUS"]] <- safe_get_msig(
  "C2",
  "CP:KEGG_MEDICUS",
  "KEGG Medicus"
)

DB[["BIOCARTA"]] <- safe_get_msig(
  "C2",
  "CP:BIOCARTA",
  "BioCarta"
)

DB[["PID"]] <- safe_get_msig(
  "C2",
  "CP:PID",
  "PID"
)

if (RUN_KEGG_LEGACY) {
  
  DB[["KEGG_LEGACY"]] <- safe_get_msig(
    "C2",
    "CP:KEGG_LEGACY",
    "KEGG Legacy"
  )
}

DB[["GO_BP"]] <- safe_get_msig(
  "C5",
  "GO:BP",
  "GO Biological Process"
)

if (RUN_ONCOGENIC) {
  
  DB[["ONCOGENIC_C6"]] <- safe_get_msig(
    "C6",
    "NULL",
    "Oncogenic signatures"
  )
}

###############################################################################
# COMBINED CANONICAL PATHWAY DATABASE
###############################################################################

CANONICAL_NAMES <- intersect(
  c(
    "REACTOME",
    "WIKIPATHWAYS",
    "KEGG_MEDICUS",
    "BIOCARTA",
    "PID",
    "KEGG_LEGACY"
  ),
  names(DB)
)

CANONICAL_LIST <- lapply(
  CANONICAL_NAMES,
  function(name) {
    
    X <- copy(
      DB[[name]]
    )
    
    if (!nrow(X)) {
      return(X)
    }
    
    # Keep source visible in term name so identical names from different
    # sources cannot collide.
    X[
      ,
      term :=
        paste0(
          name,
          "::",
          term
        )
    ]
    
    X
  }
)

DB[["CANONICAL_ALL"]] <- unique(
  rbindlist(
    CANONICAL_LIST,
    fill = TRUE
  )
)

cat(
  "\nCombined canonical pathway sets: ",
  uniqueN(
    DB[["CANONICAL_ALL"]]$term
  ),
  "\n",
  sep = ""
)

###############################################################################
# EXPERIMENT-WIDE EXPLORATORY COLLECTION
#
# This combines Hallmark + canonical pathways + GO BP.
# Use for discovery, but report the focused database-specific results as the
# primary analyses.
###############################################################################

BROAD_LIST <- lapply(
  c(
    "HALLMARK",
    "CANONICAL_ALL",
    "GO_BP"
  ),
  function(name) {
    
    X <- copy(
      DB[[name]]
    )
    
    if (!nrow(X)) {
      return(X)
    }
    
    X[
      ,
      term :=
        paste0(
          name,
          "::",
          term
        )
    ]
    
    X
  }
)

DB[["BROAD_EXPLORATORY"]] <- unique(
  rbindlist(
    BROAD_LIST,
    fill = TRUE
  )
)

###############################################################################
# DEFINE FALLBACK UNIVERSE IF NEEDED
###############################################################################

if (!HAS_REAL_UNIVERSE) {
  
  UNIVERSE <- unique(
    unlist(
      lapply(
        DB[
          c(
            "HALLMARK",
            "CANONICAL_ALL",
            "GO_BP"
          )
        ],
        function(x) {
          x$gene
        }
      ),
      use.names = FALSE
    )
  )
  
  cat(
    "Fallback MSigDB universe genes: ",
    length(UNIVERSE),
    "\n",
    sep = ""
  )
}

# Make sure all foreground genes that are testable are inside the universe.
FOREGROUND_IN_UNIVERSE <- intersect(
  FOREGROUND,
  UNIVERSE
)

cat(
  "Foreground genes in supplied/fallback universe: ",
  length(FOREGROUND_IN_UNIVERSE),
  " / ",
  length(FOREGROUND),
  "\n",
  sep = ""
)

fwrite(
  data.table(
    gene = FOREGROUND,
    in_universe =
      FOREGROUND %in%
      UNIVERSE
  ),
  file.path(
    TABDIR,
    "foreground_universe_QC.tsv"
  ),
  sep = "\t"
)

###############################################################################
# RUN ALL DATABASES
###############################################################################

RUN_ORDER <- c(
  "HALLMARK",
  "CANONICAL_ALL",
  "REACTOME",
  "WIKIPATHWAYS",
  "KEGG_MEDICUS",
  "BIOCARTA",
  "PID",
  "GO_BP"
)

if (RUN_ONCOGENIC) {
  RUN_ORDER <- c(
    RUN_ORDER,
    "ONCOGENIC_C6"
  )
}

RUN_ORDER <- c(
  RUN_ORDER,
  "BROAD_EXPLORATORY"
)

RUN_ORDER <- RUN_ORDER[
  RUN_ORDER %in%
    names(DB)
]

RESULTS <- list()

SUMMARY <- list()

for (name in RUN_ORDER) {
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "ORA: ",
    name,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================================\n"
  )
  
  R <- run_ora(
    foreground = FOREGROUND_IN_UNIVERSE,
    universe = UNIVERSE,
    term2gene = DB[[name]],
    database = name
  )
  
  RESULTS[[name]] <- R
  
  if (nrow(R)) {
    
    fwrite(
      R,
      file.path(
        TABDIR,
        paste0(
          name,
          "__ORA_all.tsv"
        )
      ),
      sep = "\t"
    )
    
    SIG <- R[
      p.adjust < ORA_FDR
    ]
    
    fwrite(
      SIG,
      file.path(
        TABDIR,
        paste0(
          name,
          "__ORA_FDR05.tsv"
        )
      ),
      sep = "\t"
    )
    
    NONRED <- reduce_redundancy(
      SIG,
      threshold = JACCARD_THRESHOLD,
      max_terms = 100
    )
    
    fwrite(
      NONRED,
      file.path(
        TABDIR,
        paste0(
          name,
          "__ORA_FDR05_nonredundant.tsv"
        )
      ),
      sep = "\t"
    )
    
    plot_ora(
      R,
      name,
      paste0(
        name,
        "__ORA_dotplot"
      )
    )
    
    best_term <- if (nrow(R)) {
      R$Description[1]
    } else {
      NA_character_
    }
    
    best_fdr <- if (nrow(R)) {
      R$p.adjust[1]
    } else {
      NA_real_
    }
    
    SUMMARY[[name]] <- data.table(
      Database = name,
      Tested_terms = nrow(R),
      Significant_FDR05 = nrow(SIG),
      Significant_nonredundant = nrow(NONRED),
      Best_term = best_term,
      Best_FDR = best_fdr
    )
    
    cat(
      "Tested terms: ",
      nrow(R),
      "\n",
      sep = ""
    )
    
    cat(
      "FDR < 0.05: ",
      nrow(SIG),
      "\n",
      sep = ""
    )
    
    if (nrow(SIG)) {
      
      cat(
        "\nTop significant terms:\n"
      )
      
      print(
        SIG[
          1:min(
            10,
            .N
          ),
          .(
            Description,
            GeneRatio,
            Count,
            p.adjust
          )
        ]
      )
    }
    
  } else {
    
    SUMMARY[[name]] <- data.table(
      Database = name,
      Tested_terms = 0L,
      Significant_FDR05 = 0L,
      Significant_nonredundant = 0L,
      Best_term = NA_character_,
      Best_FDR = NA_real_
    )
  }
}

###############################################################################
# COMBINED RESULT TABLE
###############################################################################

ALL_RESULT <- rbindlist(
  RESULTS,
  fill = TRUE
)

fwrite(
  ALL_RESULT,
  file.path(
    TABDIR,
    "ALL_DATABASES_ORA_all.tsv"
  ),
  sep = "\t"
)

ALL_SIG <- ALL_RESULT[
  p.adjust < ORA_FDR
]

fwrite(
  ALL_SIG,
  file.path(
    TABDIR,
    "ALL_DATABASES_ORA_FDR05.tsv"
  ),
  sep = "\t"
)

SUMMARY_TABLE <- rbindlist(
  SUMMARY,
  fill = TRUE
)

setorder(
  SUMMARY_TABLE,
  -Significant_FDR05,
  Best_FDR
)

fwrite(
  SUMMARY_TABLE,
  file.path(
    TABDIR,
    "ORA_database_summary.tsv"
  ),
  sep = "\t"
)

###############################################################################
# FACETED MULTI-DATABASE SUMMARY
#
# Exclude BROAD_EXPLORATORY and redundant sub-databases to avoid duplication.
###############################################################################

# 1. MÓDOSÍTÁS: Csak a 3 fő, nem redundáns kategóriát tartjuk meg a plothoz
FACET_DB <- intersect(
  c("CANONICAL_ALL", "GO_BP", "HALLMARK"), 
  RUN_ORDER
)

FACET_LIST <- list()

for (name in FACET_DB) {
  
  R <- RESULTS[[name]]
  
  if (!nrow(R)) {
    next
  }
  
  X <- R[
    p.adjust < ORA_FDR
  ]
  
  if (!nrow(X)) {
    next
  }
  
  X <- head(
    X,
    TOP_N_FACET
  )
  
  FACET_LIST[[name]] <- X
}

FACET <- rbindlist(
  FACET_LIST,
  fill = TRUE
)

if (nrow(FACET)) {
  
  # 2. MÓDOSÍTÁS: Névtisztítás (A tolower és az első betűs trükkök eltávolítva)
  FACET[, clean_term := sub("^.*::", "", Description)]
  FACET[, clean_term := sub("^[A-Z]+_", "", clean_term)]
  FACET[, clean_term := gsub("_", " ", clean_term)]
  
  # Hosszú nevek tördelése adott karakterszámnál (strwrap)
  FACET[
    ,
    clean_term := vapply(
      clean_term, 
      function(x) { paste(strwrap(x, width = TEXT_WRAP_WIDTH), collapse = "\n") }, 
      character(1)
    )
  ]
  
  # Rendezés és faktorálás a tiszta nevek alapján
  setorder(
    FACET,
    Database,
    GeneRatio_numeric
  )
  
  FACET[
    ,
    plot_label :=
      factor(
        clean_term,
        levels = unique(clean_term)
      )
  ]
  
  p <- ggplot(
    FACET,
    aes(
      GeneRatio_numeric,
      plot_label  
    )
  ) +
    
    geom_point(
      aes(
        size = Count,
        color = minus_log10_FDR
      ),
      alpha = 0.9
    ) +
    
    facet_grid(
      Database ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    
    scale_color_gradient(
      low = "grey70",
      high = "black"
    ) +
    
    labs(
      title =
        "Shared splicing genes: multi-database ORA",
      x =
        "Gene ratio",
      y =
        NULL,
      size =
        "Genes",
      color =
        "-log10(FDR)"
    ) +
    
    theme_classic(
      base_size = 14 # Megnövelt alap betűméret (korábban 10)
    ) +
    
    theme(
      strip.text.y =
        element_text(
          face = "bold",
          angle = 0,
          size = 12 # Megnövelt kategória felirat
        ),
      axis.text.y =
        element_text(
          size = 11 # Megnövelt tengely felirat
        ),
      plot.title =
        element_text(
          face = "bold",
          size = 16 # Megnövelt cím méret
        )
    )
  
  ggsave(
    file.path(
      FIGDIR,
      "MULTIDB_ORA_top_significant.png"
    ),
    p,
    width = 10,
    height = max(
      6,
      0.45 *
        nrow(FACET) +
        2
    ),
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(
      FIGDIR,
      "MULTIDB_ORA_top_significant.pdf"
    ),
    p,
    width = 10,
    height = max(
      6,
      0.45 *
        nrow(FACET) +
        2
    ),
    bg = "white"
  )
}

###############################################################################
# DATABASE COVERAGE PLOT
###############################################################################

QC_FILES <- list.files(
  TABDIR,
  pattern = "__input_QC\\.tsv$",
  full.names = TRUE
)

QC_ALL <- rbindlist(
  lapply(
    QC_FILES,
    fread
  ),
  fill = TRUE
)

fwrite(
  QC_ALL,
  file.path(
    TABDIR,
    "ORA_database_gene_coverage.tsv"
  ),
  sep = "\t"
)

if (nrow(QC_ALL)) {
  
  QC_ALL[
    ,
    foreground_coverage :=
      Foreground_in_database /
      Foreground_total
  ]
  
  p_cov <- ggplot(
    QC_ALL[
      Database !=
        "BROAD_EXPLORATORY"
    ],
    aes(
      reorder(
        Database,
        foreground_coverage
      ),
      foreground_coverage
    )
  ) +
    
    geom_col() +
    
    coord_flip() +
    
    scale_y_continuous(
      limits = c(
        0,
        1
      ),
      labels = function(x) {
        paste0(
          round(
            100 * x
          ),
          "%"
        )
      }
    ) +
    
    labs(
      title =
        "Foreground gene coverage by enrichment database",
      x =
        NULL,
      y =
        "Common genes represented"
    ) +
    
    theme_classic(
      base_size = 14 # Megnövelt alap betűméret
    )
  
  ggsave(
    file.path(
      FIGDIR,
      "database_gene_coverage.png"
    ),
    p_cov,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

###############################################################################
# README
###############################################################################

README <- c(
  "# Multi-database ORA",
  "",
  paste0(
    "Input foreground: ",
    normalizePath(
      INPUT,
      mustWork = FALSE
    )
  ),
  "",
  paste0(
    "Foreground genes: ",
    length(
      FOREGROUND
    )
  ),
  "",
  paste0(
    "Foreground genes in universe: ",
    length(
      FOREGROUND_IN_UNIVERSE
    )
  ),
  "",
  paste0(
    "Explicit tested universe supplied: ",
    HAS_REAL_UNIVERSE
  ),
  "",
  "## Recommended interpretation order",
  "",
  "1. CANONICAL_ALL — primary pathway ORA.",
  "2. HALLMARK — compact high-level biology.",
  "3. GO_BP — broad biological-process discovery.",
  "4. REACTOME / WIKIPATHWAYS / KEGG_MEDICUS — database-specific confirmation.",
  "5. ONCOGENIC_C6 — optional cancer-specific signatures.",
  "",
  "## Important",
  "",
  "For the final TDK, use the actual splicing-testable background universe.",
  "Do not use the 265 foreground/common genes themselves as the ORA universe.",
  "",
  "The BROAD_EXPLORATORY analysis pools Hallmark, canonical pathways and GO BP",
  "for discovery. The focused database-specific analyses should be used for",
  "formal reporting because their multiple-testing families are clearer.",
  ""
)

writeLines(
  README,
  file.path(
    OUT,
    "README_ORA.md"
  )
)

###############################################################################
# SESSION INFO
###############################################################################

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "sessionInfo.txt"
  )
)

###############################################################################
# FINISHED
###############################################################################

cat(
  "\n============================================================\n"
)

cat(
  "MULTI-DATABASE ORA FINISHED\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nForeground genes: ",
  length(FOREGROUND),
  "\n",
  sep = ""
)

cat(
  "Foreground in universe: ",
  length(FOREGROUND_IN_UNIVERSE),
  "\n",
  sep = ""
)

cat(
  "Explicit tested universe: ",
  HAS_REAL_UNIVERSE,
  "\n",
  sep = ""
)

cat(
  "\nDatabase summary:\n"
)

print(
  SUMMARY_TABLE
)

cat(
  "\nOutput:\n",
  normalizePath(
    OUT,
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)