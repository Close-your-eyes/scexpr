#' Compare entries of two compatible UMI matrices
#'
#' Matrices will be reduces to matching rownames and columnnames.
#' Then root mean square error (RSME), total and relative differences in UMIs counts
#' are computed. This is done for matrix as a whole as well as rowwise. Rows are expected
#' to contain features. It is expected that matrices contain positive values only.
#' An error will not be thrown if this is violated, though.
#'
#' @param mat1 matrix 1 with rownames and colnames
#' @param mat2 matrix 2 with rownames and colnames
#'
#' @return list with measures of differences or errors, as a whole and rowwise
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

  # how to work completely on sparse matrices?

  if (is.null(colnames(mat1)) || is.null(colnames(mat2)) || is.null(rownames(mat1)) || is.null(rownames(mat2))) {
    stop("mat1 and mat2 need to have rows.")
  }

  if (nrow(mat1) != nrow(mat2)) {
    message("mat1 and mat2 have different number of rows.")
  }

  if (ncol(mat1) != ncol(mat2)) {
    message("mat1 and mat2 have different number of columns.")
  }

  matching_cols <- intersect(colnames(mat1), colnames(mat2))
  matching_rows <- intersect(rownames(mat1), rownames(mat2))

  if (length(matching_rows) != nrow(mat1) || length(matching_rows) != nrow(mat2)) {

  }

  if (length(matching_cols) != ncol(mat1) || length(matching_cols) != ncol(mat2)) {
    if (length(setdiff(matching_cols, colnames(mat1))) > 0) {
      message(x, " of matching columns")
    }
    if (length(x <- setdiff(colnames(mat1), matching_cols)) > 0) {
      message(x, " of matching columns")
    }
    setdiff(colnames(mat1), matching_cols)
  }

  mat1 <- as.matrix(mat1[matching_rows,matching_cols])
  mat2 <- as.matrix(mat2[matching_rows,matching_cols])

  abs_diff <- mat1 - mat2
  squ_diff <- abs_diff^2
  rel_diff <- (abs(abs_diff)+abs(abs_diff))/(mat1 + mat2)
  # mean squared error

  mse <- Matrix::mean(squ_diff)

  # this avoids coercion to dense matrix, but is it always correct?
  # mean(c(mean(1:10), mean(11:20))) == mean(1:20)
  #mse <- mean(Matrix::colMeans(squ_diff))

  # root mean squared error (RMSE)
  rmse <- sqrt(mse)

  ## summarise the overall difference; difference in UMI to total number of UMI
  total_diff <- sum(abs(abs_diff))
  total <- sum((mat1+mat2)/2) # same sum((mat1+mat2))/2
  rel_diff <- total_diff/total

  ## how many entries are unequal
  abs_entry_diff <- sum(abs_diff != 0)
  #total_entries <- ncol(mat1)*nrow(mat1) # same length(mat1)
  rel_entry_diff <- abs_entry_diff/length(mat1)


  # _r for row
  abs_umi_r <- utils::stack(rowSums(mat1+mat2)/2)
  names(abs_umi_r)[1] <- "abs_umi"

  abs_err_r <- utils::stack(rowSums(mat1-mat2))
  names(abs_err_r)[1] <- "abs_err"

  abs_abs_err_r <- utils::stack(rowSums(abs(mat1-mat2)))
  names(abs_abs_err_r)[1] <- "total_diff"

  rel_abs_err_r <- utils::stack(stats::setNames(abs_abs_err_r$total_diff/abs_umi_r$abs_umi, abs_abs_err_r$ind))
  names(rel_abs_err_r)[1] <- "rel_total_diff"

  # purrr::map_dbl(setNames(1:nrow(mat1), rownames(mat1)), function(x) rsme(mat1[x,], mat2[x,]), .progress = T)
  rsme_r <- utils::stack(sqrt(Matrix::rowMeans((mat1-mat2)^2)))
  names(rsme_r)[1] <- "RSME"


  rowwise_error_df <- purrr::reduce(list(rsme_r,
                                         abs_umi_r,
                                         abs_err_r,
                                         abs_abs_err_r,
                                         rel_abs_err_r),
                                    dplyr::left_join, by = "ind")
  names(rowwise_error_df)[which(names(rowwise_error_df) == "ind")] <- "rowname"
  rowwise_error_df <- rowwise_error_df[,c(2,1,3,4,5,6)]
  rowwise_error_df$rowname <- as.character(rowwise_error_df$rowname)

  rowwise_error_df$abs_umi_group = cut(rowwise_error_df$abs_umi,
                                       include.lowest = T,
                                       breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8))
  levels <- c(0, levels(rowwise_error_df$abs_umi_group))
  rowwise_error_df$abs_umi_group <- as.character(rowwise_error_df$abs_umi_group)
  rowwise_error_df$abs_umi_group[which(is.na(rowwise_error_df$abs_umi_group))] <- "0"
  rowwise_error_df$abs_umi_group <- factor(rowwise_error_df$abs_umi_group, levels = levels)

  matrixwise_error_df <- data.frame(RSME = rmse,
                                    #MSE = mse,
                                    total_diff = total_diff,
                                    rel_diff = rel_diff,
                                    abs_unequal_entry = abs_entry_diff,
                                    rel_unequal_entry = rel_entry_diff)

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
