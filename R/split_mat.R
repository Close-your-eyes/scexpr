#' Split matrices just like splitting a data frame.
#'
#' split.data.frame method from base R actually does the job.
#' The name is very intuitive though and by default it
#' only splits along rows. For splitting by column two
#' transposations would be needed which is not very elegant
#' in case of large matrices. So here you can define if
#' you want to split the matrix along its rows or columns.
#'
#' @param x a matrix
#' @param f a character or factor vector for splitting,
#' must be same length as rows or cols of x
#' @param byrow split by rows (T) or by columns (F)
#' @param ... arguments passed to split
#'
#' @return
#' @export
#'
#' @examples
split_mat <- function(x,
                      f,
                      n_chunks,
                      chunk_size,
                      byrow = T, ...) {


  ## https://stackoverflow.com/questions/62161916/is-there-a-function-in-r-that-splits-a-matrix-along-a-margin-using-a-factor-or-c
  ## modified from base function split.data.frame (to avoid 2 x transposation)
  ## multi dirs: split count matrix by orig.idents

  if (!is.logical(byrow)) {
    stop("byrow should be logical, TRUE or FALSE.")
  }

  if (!is.matrix(x)) {
    stop("x should be a matrix.")
    # sparse?
  }

  if (missing(f) && missing(n_chunks) && missing(chunk_size)) {
    stop("f, n_chunks or chunk_size must be provided.")
  }

  if (sum(c(!missing(f), !missing(n_chunks), !missing(chunk_size))) > 1) {
    message("more than one argument of f, n_chunks and chunk_size provided. will ignore f.")
    f <- NULL
    if (sum(c(!missing(n_chunks), !missing(chunk_size))) > 1) {
      message("n_chunks and chunk_size provided. will ignore n_chunks.")
      n_chunks <- NULL
    }
  }


  if (missing(f) || is.null(f)) {

    if (!missing(n_chunks) && !is.null(n_chunks)) {
      if (byrow) {
        chunk_size <- ceiling(nrow(x)/n_chunks)
      } else {
        chunk_size <- ceiling(ncol(x)/n_chunks)
      }
    }

    if (byrow) {
      d <- 1:nrow(x)
    } else {
      d <- 1:ncol(x)
    }
    if (chunk_size > length(d)) {
      stop("chunk_size is larger than nrow or ncol. This cannot be.")
    }
    f <- split(d, ceiling(seq_along(d)/chunk_size))
    f <- rep(names(f), lengths(f))
  }

  if (byrow) {
    if (length(f) != nrow(x)) {
      stop("length(f) should be equal to nrow(x).")
    }
    split_x <- lapply(split(x = seq_len(nrow(x)), f = f, ...), function(ind) x[ind,,drop = FALSE])
  } else {
    if (length(f) != ncol(x)) {
      stop("length(f) should be equal to ncol(x).")
    }
    split_x <- lapply(split(x = seq_len(ncol(x)), f = f, ...), function(ind) x[,ind,drop = FALSE])
  }

  return(split_x)
}

