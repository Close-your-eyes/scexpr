#' Add Feature Values as a Dimensional Reduction
#'
#' Creates a two-dimensional Seurat dimensional reduction from the expression
#' values of two features. Duplicate feature names are removed, and only the
#' first two unique features are used.
#'
#' @param obj A Seurat object.
#' @param features A character vector of feature names. At least two unique
#'   features must be supplied.
#'
#' @return The Seurat object with a new dimensional reduction added to
#'   `obj@reductions`. The reduction name is formed by concatenating the two
#'   feature names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- add_feature_dimreduc(
#'   obj,
#'   features = c("CD3D", "CD8A")
#' )
#'
#' Embeddings(obj, reduction = "CD3DCD8A")
#' }
add_feature_dimreduc <- function(obj, features) {
  features <- unique(features)[c(1,2)]
  key <- paste0(features, collapse = "")
  mat <- get_layer(obj, features = features, transpose = T, as = "dense")
  colnames(mat) <- paste0(key, "_", c(1,2))
  # colnames(mat) <- paste0(colnames(mat), "_", c(1,2))
  obj@reductions[[key]] <- SeuratObject::CreateDimReducObject(embeddings = mat, assay = "RNA", key = paste0(key, "_"))
  # colnames(obj@reductions[[key]]@cell.embeddings) <- colnames(mat)
  return(obj)
}
