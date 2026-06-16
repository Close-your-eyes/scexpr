#' Downsample a Seurat object
#'
#' @param obj A Seurat object.
#' @param downsample Either:
#'   - a proportion in (0, 1], or
#'   - the number of cells to retain (> 1).
#' @param group Optional metadata column used for stratified sampling.
#'
#' @returns A downsampled Seurat object.
#' @export
downsample_seurat <- function(obj, downsample = 1, group = NULL) {

  if (downsample == 1) {
    return(obj)
  }
  stopifnot(inherits(obj, "Seurat"))

  if (!is.null(group)) {
    group <- rlang::sym(group)
  }
  if (downsample < 1) {
    cells <- rownames(dplyr::slice_sample(obj@meta.data, prop = downsample, by = group))
  } else if (downsample > 1) {
    cells <- rownames(dplyr::slice_sample(obj@meta.data, n = downsample, by = group))
  }
  obj <- subset(obj, cells = cells)
  return(obj)
}
