#' Copy dimension reductions from a cell data set (cds) object to respective slots of a Seurat object (SO)
#'
#' cds may have been created by monocle3.
#'
#' @param cds a cell data set
#' @param SO an Seurat object
#' @param reductions which reductions to copy; must be one or multiple of names(SingleCellExperiment::reducedDims(cds))
#' @param suffix which suffix to add to the reduction names in the Seurat object
#' @param assay which assay the dimension reduction were originally based on; must be one of names(SO@assays)
#'
#' @return
#' @export
#'
#' @examples
cds_dimred_to_seurat <- function(cds,
                                 SO,
                                 reductions = NULL,
                                 suffix = "monocle",
                                 assay = "RNA") {


  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }

  if (!requireNamespace("SingleCellExperiment", quietly = T)) {
    BiocManager::install("SingleCellExperiment")
  }

  if (grepl("_", suffix)) {
    stop("suffix should not contain an underscore.")
  }

  if (missing(cds)) {
    stop("cds is required.")
  }

  if (missing(SO)) {
    stop("SO is required.")
  }

  assay <- match.arg(assay, names(SO@assays))

  if (is.null(reductions)) {
    reductions <- names(SingleCellExperiment::reducedDims(cds))
  } else {
    reductions <- intersect(reductions, names(SingleCellExperiment::reducedDims(cds)))
    if (length(reductions) == 0) {
      stop("Non of reductions found in cds.")
    }
  }

  if (length(colnames(SO)) != length(colnames(cds))) {
    stop("Different number of cells in cds and SO.")
  }

  if (length(intersect(colnames(cds), colnames(SO))) != length(colnames(SO))) {
    stop("Different cells or cellnames in cds and SO. Cannot match dimension reduction information.")
  }

  for (i in names(SingleCellExperiment::reducedDims(cds))) {
    emb <- SingleCellExperiment::reducedDims(cds)[[i]][colnames(SO),]
    colnames(emb) <- paste0(i, "_", seq(1,ncol(emb)))
    SO[[paste0(i, suffix)]] <- Seurat::CreateDimReducObject(embeddings = emb,
                                                            key = paste0(i, suffix, "_"),
                                                            assay = assay)
  }

  return(SO)
}
