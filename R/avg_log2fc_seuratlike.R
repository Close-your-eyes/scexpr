#' Title
#'
#' @param mat1
#' @param mat2
#'
#' @returns
#' @export
#'
#' @examples
avg_log2fc_seuratlike <- function(mat1, mat2) {

  ##  avg exprs calc
  pseudocount.use <- 1
  base <- 2
  default.mean.fxn <- function(x) {
    return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = base))
  }

  ## no checking

  return(default.mean.fxn(mat1) - default.mean.fxn(mat2))

}
