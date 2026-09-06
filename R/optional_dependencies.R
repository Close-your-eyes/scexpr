.scexpr_bioconductor_packages <- c(
  "AnnotationDbi",
  "BiocParallel",
  "DropletUtils",
  "GO.db",
  "MAST",
  "SingleCellExperiment",
  "SingleR",
  "SummarizedExperiment",
  "UCell",
  "biomaRt",
  "celda",
  "fgsea",
  "limma",
  "metapod",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "scrapper",
  "scuttle"
)

.scexpr_pak_packages <- c(
  CellChat = "jinworks/CellChat",
  ProjecTILs = "carmonalab/ProjecTILs",
  STACAS = "carmonalab/STACAS",
  SeuratWrappers = "satijalab/seurat-wrappers",
  brathering = "close-your-eyes/brathering",
  colrr = "close-your-eyes/colrr",
  fcexpr = "close-your-eyes/fcexpr",
  monocle3 = "cole-trapnell-lab/monocle3",
  presto = "immunogenomics/presto",
  scDblFinder = "plger/scDblFinder",
  scGate = "carmonalab/scGate"
)

.scexpr_cran_packages <- c(
  "FactoMineR", "MASS", "RANN", "SoupX", "caret", "clue", "cluster",
  "cowplot", "diptest", "factoextra", "forcats", "ggdendro", "ggforce",
  "ggnewscale", "ggpubr", "ggrepel", "ggtext", "ggraph", "glue", "gt",
  "harmony", "igraph", "knitr", "lobstr", "matrixTests", "mclust",
  "msigdbr", "pROC", "patchwork", "psych", "readr",
  "rmarkdown", "sp", "stringdist", "uwot", "vroom", "xgboost"
)

.scexpr_optional_packages <- c(
  .scexpr_bioconductor_packages,
  names(.scexpr_pak_packages),
  .scexpr_cran_packages
)

.ensure_package <- function(package) {
  if (!package %in% .scexpr_optional_packages) {
    stop("Optional package '", package, "' is not registered by scexpr.", call. = FALSE)
  }

  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  message("Installing optional package '", package, "'.")

  tryCatch(
    {
      if (package %in% .scexpr_bioconductor_packages) {
        BiocManager::install(package, ask = FALSE, update = FALSE)
      } else if (package %in% names(.scexpr_pak_packages)) {
        pak::pak(unname(.scexpr_pak_packages[[package]]))
      } else {
        utils::install.packages(package)
      }
    },
    error = function(error) {
      stop(
        "Installation of optional package '", package, "' failed: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Optional package '", package,
      "' is required for this operation but could not be installed.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.ensure_packages <- function(packages) {
  invisible(lapply(unique(packages), .ensure_package))
}
