#' Compute Seurat-like Average Log2 Fold Change
#'
#' Computes the difference in average log-transformed expression between two
#' groups of cells using the same mean expression calculation employed by
#' Seurat's differential expression workflow. Expression values are assumed to
#' be log1p-transformed, are converted back to linear space, averaged with a
#' pseudocount, and transformed back to log2 space.
#'
#' @param mat1 A numeric matrix of log1p-transformed gene expression values for
#'   the first group of cells, with genes in rows and cells in columns.
#' @param mat2 A numeric matrix of log1p-transformed gene expression values for
#'   the second group of cells, with genes in rows and cells in columns.
#'
#' @returns A numeric vector containing the average log2 fold change for each
#'   gene (`mat1 - mat2`).
#' @export
#'
#' @examples
#' \dontrun{
#' avg_log2fc <- avg_log2fc_seuratlike(group1_matrix, group2_matrix)
#' }
avg_log2fc_seuratlike <- function(mat1, mat2) {

  ##  avg exprs calc
  pseudocount.use <- 1
  base <- 2
  ## no checking

  return(default.mean.fxn(mat1) - default.mean.fxn(mat2))

}

default.mean.fxn <- function(x) {
  return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = base))
}
