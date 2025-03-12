#' Recycle vectors or lists
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
recycle <- function(x, y) {
  if (length(y) > length(x)) {
    recycle2(y, x)
  } else {
    recycle2(x, y)
  }
}

recycle2 <- function(x, y) {
  if (length(x) %% length(y) != 0) {
    message("recycling: longer element is not a multiple of the shorter.")
  }
  return(rep(y, length.out = length(x)))
}
