#' Return Cell Names from a Seurat Object
#'
#' Extracts the cell identifiers (rownames of the metadata table) from a
#' Seurat object.
#'
#' Seurat::Cells sometimes fails after subset.
#'
#' @param obj A Seurat object.
#'
#' @return A character vector containing the cell names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Create or load a Seurat object
#' data("pbmc_small")
#'
#' # Get cell names
#' cells2(pbmc_small)
#' }
cells2 <- function(obj) {
  stopifnot(inherits(obj, "Seurat"))
  return(rownames(obj@meta.data))
}
