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

  matching_cols <- intersect(colnames(mat1), colnames(mat2))
  matching_rows <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- as.matrix(mat1[matching_rows,matching_cols])
  mat2 <- as.matrix(mat2[matching_rows,matching_cols])

  abs_diff <- mat1 - mat2
  squ_diff <- abs_diff^2
  rel_diff <- (abs(abs_diff)+abs(abs_diff))/(mat1+mat2)
  # mean squared error
  mse <- mean(squ_diff)
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
  rsme_r <- utils::stack(purrr::map_dbl(setNames(1:nrow(mat1), rownames(mat1)), function(x) rsme(mat1[x,], mat2[x,])))
  names(rsme_r)[1] <- "RSME"

  abs_umi_r <- utils::stack(purrr::map_dbl(setNames(1:nrow(mat1), rownames(mat1)), function(x) abs_umi(mat1[x,], mat2[x,])))
  names(abs_umi_r)[1] <- "abs_umi"

  abs_err_r <- utils::stack(purrr::map_dbl(setNames(1:nrow(mat1), rownames(mat1)), function(x) abs_err(mat1[x,], mat2[x,])))
  names(abs_err_r)[1] <- "abs_err"

  abs_abs_err_r <- utils::stack(purrr::map_dbl(setNames(1:nrow(mat1), rownames(mat1)), function(x) abs_abs_err(mat1[x,], mat2[x,])))
  names(abs_abs_err_r)[1] <- "total_diff"

  rel_abs_err_r <- utils::stack(purrr::map_dbl(setNames(1:nrow(mat1), rownames(mat1)), function(x) rel_abs_err(mat1[x,], mat2[x,])))
  names(rel_abs_err_r)[1] <- "rel_total_diff"

  rowwise_error_df <- purrr::reduce(list(rsme_r,
                                         abs_umi_r,
                                         abs_err_r,
                                         abs_abs_err_r,
                                         rel_abs_err_r),
                                    dplyr::left_join, by = "ind")
  names(rowwise_error_df)[which(names(rowwise_error_df) == "ind")] <- "rowname"
  rowwise_error_df <- rowwise_error_df[,c(2,1,3,4,5,6)]
  rowwise_error_df$rowname <- as.character(rowwise_error_df$rowname)

  return(list(RSME = rmse,
              MSE = mse,
              total_diff = total_diff,
              rel_diff = rel_diff,
              abs_unequal_entry = abs_entry_diff,
              rel_unequal_entry = rel_entry_diff,
              rowwise_error = rowwise_error_df))

}



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
