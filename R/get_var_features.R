#' Get Variable Features from a Seurat Object
#'
#' Retrieves variable features from a Seurat object. The function first checks
#' the current default assay for variable features. If none are found, it
#' iterates through all other assays until variable features are identified.
#'
#' @param obj A \code{Seurat} object.
#'
#' @returns A character vector of variable feature names if found. Returns
#' \code{NULL} if no variable features are present in any assay.
#'
#' @details
#' The function searches for variable features in the default assay using
#' \code{SeuratObject::VariableFeatures()}. If no features are found, it
#' sequentially switches to other assays in the object and repeats the search.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#'
#' # Create a small example Seurat object
#' data <- matrix(rpois(2000, lambda = 5), nrow = 100)
#' rownames(data) <- paste0("Gene", 1:100)
#' colnames(data) <- paste0("Cell", 1:20)
#'
#' obj <- CreateSeuratObject(counts = data)
#' obj <- FindVariableFeatures(obj)
#'
#' # Retrieve variable features
#' vf <- get_var_features(obj)
#' head(vf)
#' }
get_var_features <- function(obj) {
  varfeat <- SeuratObject::VariableFeatures(obj)
  if (!is.null(varfeat)) {
    message("var features from assay ", Seurat::DefaultAssay(obj))
    return(varfeat)
  } else {
    assays <- setdiff(names(obj@assays), Seurat::DefaultAssay(obj))
    for (i in assays) {
      Seurat::DefaultAssay(obj) <- i
      varfeat <- get_var_features(obj)
      if (!is.null(varfeat)) return(varfeat)
    }
  }
  if (is.null(varfeat)) {
    message("no var features in any assay found.")
    return(NULL)
  }
}
