#' Subset an Object and Remove Empty Assay Entries
#'
#' Subsets an object using its \code{subset()} method and restricts each assay's
#' cell and feature membership maps to entries present in the \code{counts}
#' layer. Optionally removes groups represented by fewer than a specified
#' number of cells.
#'
#' @param obj An object containing \code{assays} and \code{meta.data} slots.
#'   Assays are expected to contain \code{cells} and \code{features} membership
#'   maps with a \code{counts} column.
#' @param ... Additional arguments passed to \code{\link[base]{subset}}.
#' @param min_obs A numeric value specifying the minimum number of cells
#'   required for each group defined by \code{by}. Group filtering is performed
#'   only when \code{min_obs > 1}. Defaults to \code{1}.
#' @param by A character string naming the column in \code{obj@meta.data} used
#'   to group cells when applying \code{min_obs}. Defaults to
#'   \code{"orig.ident"}.
#'
#' @return The subsetted object with assay cell and feature membership maps
#'   restricted to entries present in the \code{counts} layer. When
#'   \code{min_obs > 1}, only groups containing at least \code{min_obs} cells
#'   are retained. If \code{by} is not found, or if no groups meet the
#'   threshold, the object is returned without group filtering.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj_subset <- subset2(obj, subset = cell_type == "T cell")
#'
#' obj_subset <- subset2(
#'   obj,
#'   subset = cell_type == "T cell",
#'   min_obs = 10,
#'   by = "orig.ident"
#' )
#'
#' # manual subset:
#' so@assays$RNA@layers$data <- so@assays$RNA@layers$data[, inds]
#' so@assays$RNA@cells <- so@assays$RNA@cells[inds,]
#' so@assays$ADT <- NULL
#' so@assays$HTO <- NULL
#' so@graphs <- list()
#' so@meta.data <- so@meta.data[inds,]
#' so@reductions$umap_harmony@cell.embeddings <- so@reductions$umap_harmony@cell.embeddings[inds,]
#' so@reductions$umap_integrated@cell.embeddings <- so@reductions$umap_integrated@cell.embeddings[inds,]
#' so@reductions$umap_simple@cell.embeddings <- so@reductions$umap_simple@cell.embeddings[inds,]
#' Seurat::Idents(so) <- so@meta.data$orig.ident
#' }
subset2 <- function(obj, ..., min_obs = 1, by = "orig.ident") {
  obj <- subset(obj, ...)

  for (i in names(obj@assays)) {

    try(obj@assays[[i]]@cells <- obj@assays[[i]]@cells[which(obj@assays[[i]]@cells@.Data[,"counts"]),,drop = FALSE], silent = T)
    try(obj@assays[[i]]@features <- obj@assays[[i]]@features[which(obj@assays[[i]]@features@.Data[,"counts"]),,drop = FALSE], silent = T)

  }

  if (min_obs > 1) {
    if (!by %in% names(obj@meta.data)) {
      message("'by' not found in meta.data.")
      return(obj)
    }

    nobs <- dplyr::summarise(obj@meta.data, n = dplyr::n(), .by = !!rlang::sym(by))
    filt <- sum(nobs$n<min_obs)

    nobs2 <- dplyr::filter(nobs, n >= min_obs)
    cells <- rownames(dplyr::filter(obj@meta.data, !!rlang::sym(by) %in% nobs2[[by]]))

    if (!length(cells)) {
      message("no cells left. will not filter.")
      return(obj)
    } else if (filt>0) {
      message("filtering ", nrow(obj@meta.data)-length(cells), " cells from ", filt, " individuals from ", by, ".")
    }

    obj <- subset(obj, cells = cells)
    for (i in names(obj@assays)) {

      try(obj@assays[[i]]@cells <- obj@assays[[i]]@cells[which(obj@assays[[i]]@cells@.Data[,"counts"]),,drop = FALSE], silent = T)
      try(obj@assays[[i]]@features <- obj@assays[[i]]@features[which(obj@assays[[i]]@features@.Data[,"counts"]),,drop = FALSE], silent = T)

    }
  }

  return(obj)
}
