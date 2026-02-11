#' Rename a reduction in Seurat at all necessary positions
#'
#' @param obj seurat object
#' @param reduction reduction name (current)
#' @param new_name new name to assign
#'
#' @returns seurat object
#' @export
#'
#' @examples
rename_reduction <- function(obj,
                             reduction,
                             new_name) {

  reduction <- scexpr:::check.reduction(SO = obj,
                                        reduction = reduction)

  if (new_name %in% names(obj@reduction)) {
    stop("new name already exists in obj@reduction")
  }

  names(obj@reductions)[reduction] <- new_name
  obj@reductions[[new_name]]@key <- paste0(new_name, "_")
  colnames(obj@reductions[[new_name]]@cell.embeddings) <- paste0(new_name, "_", seq(1,ncol(obj@reductions[[new_name]]@cell.embeddings),1))

  return(obj)
}
