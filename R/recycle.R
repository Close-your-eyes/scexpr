#' Recycle vectors or lists
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
recycle <- function(short, long) {

  if (length(long) < length(short)) {
    message("recycling: short is longer than long and will be cut to length(long).")
  } else if (length(long) %% length(short) != 0) {
    message("recycling: longer element is not a multiple of the shorter.")
  }

  return(rep(short, length.out = length(long)))
}
