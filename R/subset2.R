#' Subset an Object and Remove Empty Assay Entries
#'
#' Subsets an object using its \code{subset()} method, then removes cells and
#' features from each assay when they are absent from the assay's
#' \code{counts} layer.
#'
#' @param obj An object containing an \code{assays} slot. Each assay must
#'   contain \code{cells} and \code{features} membership maps with a
#'   \code{counts} column.
#' @param ... Additional arguments passed to \code{\link[base]{subset}}.
#'
#' @return The subsetted object with assay cell and feature membership maps
#'   restricted to entries present in the \code{counts} layer.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj_subset <- subset2(obj, subset = cell_type == "T cell")
#' }
subset2 <- function(obj, ...) {
  obj <- subset(obj, ...)

  for (i in names(obj@assays)) {

    try(obj@assays[[i]]@cells <- obj@assays[[i]]@cells[which(obj@assays[[i]]@cells@.Data[,"counts"]),,drop = FALSE], silent = T)
    try(obj@assays[[i]]@features <- obj@assays[[i]]@features[which(obj@assays[[i]]@features@.Data[,"counts"]),,drop = FALSE], silent = T)

  }

  return(obj)
}
