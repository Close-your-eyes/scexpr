#' Compare two UMI count matrices
#'
#' Compares two UMI count matrices after aligning them to their shared
#' rownames (features) and colnames (cells). Rows and columns that are not
#' present in both matrices are discarded before comparison.
#'
#' The function computes several measures of disagreement between matrices,
#' both globally and for each feature. Metrics include root mean squared
#' error (RMSE), absolute differences in UMI counts, and relative differences
#' normalized by total abundance.
#'
#' Sparse matrices from the Matrix package are supported and processed
#' without densification, making the function suitable for large
#' single-cell expression matrices.
#'
#' @param mat1 First UMI count matrix. Rows should represent features and
#'   columns should represent cells. Must contain rownames and colnames.
#'
#' @param mat2 Second UMI count matrix. Rows should represent features and
#'   columns should represent cells. Must contain rownames and colnames.
#'
#' @return A list with two data frames:
#'
#' \describe{
#'   \item{matrixwise_error}{
#'     One-row summary of differences between the aligned matrices.
#'     Contains:
#'     \describe{
#'       \item{RMSE}{Root mean squared error across all matrix entries.}
#'       \item{total_diff}{Sum of absolute differences,
#'         \eqn{\sum |mat1 - mat2|}.}
#'       \item{rel_diff}{Total absolute difference relative to the average
#'         total UMI count,
#'         \eqn{\sum |mat1 - mat2| / ((\sum mat1 + \sum mat2)/2)}.}
#'       \item{abs_unequal_entry}{Number of matrix entries with different
#'         values between matrices.}
#'       \item{rel_unequal_entry}{Fraction of entries that differ,
#'         \eqn{abs\_unequal\_entry / (nrow \times ncol)}.}
#'     }
#'   }
#'
#'   \item{rowwise_error}{
#'     Feature-wise summary of differences, with one row per feature
#'     (matrix row). Contains:
#'     \describe{
#'       \item{rowname}{Feature identifier.}
#'       \item{RMSE}{Root mean squared error across all cells for the feature.}
#'       \item{abs_umi}{Average total abundance of the feature across the two
#'         matrices,
#'         \eqn{(rowSums(mat1) + rowSums(mat2))/2}.}
#'       \item{abs_err}{Signed difference in abundance,
#'         \eqn{rowSums(mat1 - mat2)}. Positive values indicate larger counts
#'         in \code{mat1}.}
#'       \item{total_diff}{Sum of absolute per-cell differences for the
#'         feature,
#'         \eqn{rowSums(|mat1 - mat2|)}.}
#'       \item{rel_total_diff}{Feature-level total absolute difference relative
#'         to feature abundance,
#'         \eqn{total\_diff / abs\_umi}.}
#'       \item{abs_umi_group}{Abundance bin derived from \code{abs_umi} for
#'         stratified analyses and plotting.}
#'     }
#'   }
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' mat1 <- random_matrix_normal <- matrix(rnorm(5 * 6), nrow = 5, ncol = 6)
#' mat2 <- random_matrix_normal <- matrix(rnorm(5 * 6), nrow = 5, ncol = 6)
#' colnames(mat1) <- colnames(mat2) <- 1:ncol(mat1)
#' rownames(mat1) <- rownames(mat2) <- 1:nrow(mat1)
#' mat_diff <- compare_UMI_matrices(mat1 = mat1, mat2 = mat2)
#' }
compare_UMI_matrices <- function(mat1, mat2) {

  if (is.null(rownames(mat1)) || is.null(rownames(mat2)) ||
      is.null(colnames(mat1)) || is.null(colnames(mat2))) {
    stop("mat1 and mat2 must have rownames and colnames.")
  }

  if (nrow(mat1) != nrow(mat2)) {
    message("mat1 and mat2 have different numbers of rows.")
  }

  if (ncol(mat1) != ncol(mat2)) {
    message("mat1 and mat2 have different numbers of columns.")
  }

  matching_rows <- intersect(rownames(mat1), rownames(mat2))
  matching_cols <- intersect(colnames(mat1), colnames(mat2))

  if (length(matching_rows) == 0) {
    stop("No matching rownames found.")
  }

  if (length(matching_cols) == 0) {
    stop("No matching colnames found.")
  }

  if (length(matching_rows) < nrow(mat1)) {
    message(
      nrow(mat1) - length(matching_rows),
      " rows from mat1 removed."
    )
  }

  if (length(matching_rows) < nrow(mat2)) {
    message(
      nrow(mat2) - length(matching_rows),
      " rows from mat2 removed."
    )
  }

  if (length(matching_cols) < ncol(mat1)) {
    message(
      ncol(mat1) - length(matching_cols),
      " columns from mat1 removed."
    )
  }

  if (length(matching_cols) < ncol(mat2)) {
    message(
      ncol(mat2) - length(matching_cols),
      " columns from mat2 removed."
    )
  }

  message(length(matching_cols), " columns retained.")

  mat1 <- mat1[matching_rows, matching_cols, drop = FALSE]
  mat2 <- mat2[matching_rows, matching_cols, drop = FALSE]
  n_entries <- nrow(mat1) * ncol(mat1)
  diff_mat <- mat1 - mat2

  ## ------------------------------------------------------------
  ## Matrix-wide metrics
  ## ------------------------------------------------------------

  if (inherits(diff_mat, "sparseMatrix")) {

    #table(diff_mat@x)
    mse <- sum(diff_mat@x^2) / n_entries
    total_diff <- sum(abs(diff_mat@x))
    abs_entry_diff <- sum(diff_mat@x != 0)
    abs_diff_mat <- diff_mat
    abs_diff_mat@x <- abs(abs_diff_mat@x)
    sq_diff_mat <- diff_mat
    sq_diff_mat@x <- sq_diff_mat@x^2

  } else {

    mse <- mean(diff_mat^2)
    total_diff <- sum(abs(diff_mat))
    abs_entry_diff <- sum(diff_mat != 0)
    abs_diff_mat <- abs(diff_mat)
    sq_diff_mat <- diff_mat^2
  }

  rmse <- sqrt(mse)

  total <- (sum(mat1) + sum(mat2)) / 2

  rel_diff <- if (total == 0) {
    NA_real_
  } else {
    total_diff / total
  }

  rel_entry_diff <- abs_entry_diff / n_entries

  ## ------------------------------------------------------------
  ## Row-wise metrics
  ## ------------------------------------------------------------

  abs_umi_r <- (Matrix::rowSums(mat1) + Matrix::rowSums(mat2)) / 2
  abs_err_r <- Matrix::rowSums(diff_mat)
  total_diff_r <- Matrix::rowSums(abs_diff_mat)
  rel_total_diff_r <- total_diff_r / abs_umi_r
  rel_total_diff_r[abs_umi_r == 0] <- NA_real_
  rmse_r <- sqrt(Matrix::rowSums(sq_diff_mat) / ncol(mat1))

  rowwise_error_df <- data.frame(
    rowname = rownames(mat1),
    RMSE = as.numeric(rmse_r),
    abs_umi = as.numeric(abs_umi_r),
    abs_err = as.numeric(abs_err_r),
    total_diff = as.numeric(total_diff_r),
    rel_total_diff = as.numeric(rel_total_diff_r),
    stringsAsFactors = FALSE
  )

  ## ------------------------------------------------------------
  ## UMI abundance groups
  ## ------------------------------------------------------------

  rowwise_error_df$abs_umi_group <- cut(
    rowwise_error_df$abs_umi,
    include.lowest = TRUE,
    breaks = c(
      1e0, 1e1, 1e2, 1e3,
      1e4, 1e5, 1e6, 1e7, 1e8
    )
  )

  umi_levels <- c("0", levels(rowwise_error_df$abs_umi_group))
  rowwise_error_df$abs_umi_group <- as.character(rowwise_error_df$abs_umi_group)
  rowwise_error_df$abs_umi_group[is.na(rowwise_error_df$abs_umi_group)] <- "0"

  rowwise_error_df$abs_umi_group <- factor(
    rowwise_error_df$abs_umi_group,
    levels = umi_levels
  )

  ## ------------------------------------------------------------
  ## Matrix-wise summary
  ## ------------------------------------------------------------

  matrixwise_error_df <- data.frame(
    RMSE = rmse,
    total_diff = total_diff,
    rel_diff = rel_diff,
    abs_unequal_entry = abs_entry_diff,
    rel_unequal_entry = rel_entry_diff
  )

  return(list(matrixwise_error = matrixwise_error_df,
              rowwise_error = rowwise_error_df))
}



## obsolete funs: very slow.

rsme <- function(x, y) {
  squared_diff <- (x - y)^2
  mse <- mean(squared_diff)
  rmse <- sqrt(mse)
  return(rmse)
}

abs_umi <- function(x, y) {
  return(sum(x+y))
}
abs_err <- function(x, y) {
  return(sum(x-y))
}

abs_abs_err <- function(x, y) {
  return(sum(abs(x-y)))
}

rel_abs_err <- function(x, y) {
  return(sum(abs(x-y))/(sum(x+y)))
}
