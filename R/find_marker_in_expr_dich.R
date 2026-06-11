#' Find marker features among expressing cells between two groups
#'
#' Performs differential expression testing between two groups using only
#' cells with non-zero expression for each feature. For every feature (gene),
#' the function compares expression values among expressing cells in the two
#' groups using a Wilcoxon rank-sum test.
#'
#' Unlike standard marker detection methods that include zero-expression cells,
#' this function evaluates differences conditional on expression, making it
#' useful for identifying changes in expression magnitude among expressing
#' cells rather than changes in detection rate.
#'
#' @param obj A Seurat object or a numeric matrix with features (e.g. genes)
#'   in rows and cells in columns.
#' @param group Either:
#'   \itemize{
#'     \item a character string giving the name of a metadata column in a
#'       Seurat object, or
#'     \item a vector of length \code{ncol(obj)} defining two cell groups.
#'   }
#'   Exactly two groups must be present.
#' @param get_layer_args Named list of additional arguments passed to
#'   \code{get_layer()} when \code{obj} is a Seurat object.
#' @param mc.cores Number of cores passed to
#'   \code{\link[parallel]{mclapply}}.
#'
#' @return A data frame with one row per tested feature containing:
#'   \describe{
#'     \item{feature}{Feature (gene) name.}
#'     \item{pval}{Wilcoxon rank-sum test p-value.}
#'     \item{mean_diff}{Difference in mean expression among expressing cells
#'       (\code{group.1 - group.2}).}
#'     \item{log2fc}{Log2 fold-change of mean expression among expressing
#'       cells (\code{group.1 / group.2}).}
#'     \item{pct.1}{Fraction of cells expressing the feature in group 1.}
#'     \item{pct.2}{Fraction of cells expressing the feature in group 2.}
#'     \item{padj}{Benjamini-Hochberg adjusted p-value.}
#'     \item{group.1}{Name of the first group.}
#'     \item{group.2}{Name of the second group.}
#'   }
#'
#' Features with no expressing cells in either group are excluded from the
#' results.
#'
#' @details
#' For each feature, cells with expression values greater than zero are
#' extracted separately within each group. A Wilcoxon rank-sum test is then
#' performed on these non-zero values only. Features with no expressing cells
#' in either group are skipped.
#'
#' This approach is useful when the biological question concerns expression
#' levels among expressing cells rather than overall abundance driven by
#' differences in detection frequency.
#'
#' @export
#'
#' @examples
#' # Using a Seurat object
#' res <- find_marker_in_expr_dich(
#'   obj = seu,
#'   group = "cell_type"
#' )
#'
#' # Using an expression matrix
#' grp <- rep(c("A", "B"), each = ncol(mat) / 2)
#' res <- find_marker_in_expr_dich(
#'   obj = mat,
#'   group = grp
#' )
find_marker_in_expr_dich <- function(obj,
                                     group,
                                     get_layer_args = list(),
                                     mc.cores = 4) {

  # dich  =  dichotomous

  if (length(group) == 1 && !methods::is(obj, "Seurat")) {
    stop("when obj is a matrix, group has to be a vector.")
  } else if (length(group) == 1) {
    group <- obj@meta.data[[group]]
  }
  group <- as.character(group)

  if (methods::is(obj, "Seurat")) {
    obj <- Gmisc::fastDoCall(get_layer, args = c(list(obj = obj), get_layer_args))
  }

  if (ncol(obj) != length(group)) {
    stop("column number of obj not equal to length of group.")
  }

  if (length(unique(group)) != 2) {
    stop("group must have two levels.")
  }

  obj <- brathering::split_mat(x = obj, f = group, byrow = F)

  idx <- stats::setNames(seq_len(nrow(obj[[1]])), rownames(obj[[1]]))
  out <- parallel::mclapply(idx, function(i) {

    v1 <- obj[[1]][i, ]
    v2 <- obj[[2]][i, ]

    z1 <- v1[v1 > 0]
    z2 <- v2[v2 > 0]

    if (length(z1) == 0L || length(z2) == 0L) return(NULL)
    m1 <- mean(z1)
    m2 <- mean(z2)

    df <- data.frame(
      pval = wilcox.test(z1, z2, exact = FALSE)$p.value,
      mean_diff = m1 - m2,
      log2fc = log2(m1) - log2(m2),
      pct.1 = length(z1)/length(v1),
      pct.2 = length(z2)/length(v2)
    )
    return(df)
  }, mc.cores = mc.cores)

  out <- purrr::discard(out, is.null)
  out <- dplyr::bind_rows(out, .id = "feature")
  out$padj <- stats::p.adjust(out$pval, method = "BH")
  out$group.1 <- names(obj)[1]
  out$group.2 <- names(obj)[2]

  return(out)
}


#' Find marker features among expressing cells using one-vs-rest comparisons
#'
#' Performs differential expression testing for each group against all
#' remaining cells, using only cells with non-zero expression for each
#' feature. For every unique group level, the function calls
#' \code{find_marker_in_expr_dich()} to compare that group against the union
#' of all other groups.
#'
#' This function is useful for identifying group-specific marker features
#' when more than two groups are present.
#'
#' @param obj A Seurat object or a numeric matrix with features (e.g. genes)
#'   in rows and cells in columns.
#' @param group Either:
#'   \itemize{
#'     \item a character string giving the name of a metadata column in a
#'       Seurat object, or
#'     \item a vector of length \code{ncol(obj)} defining cell group
#'       assignments.
#'   }
#' @param get_layer_args Named list of additional arguments passed to
#'   \code{get_layer()} when \code{obj} is a Seurat object.
#' @param mc.cores Number of cores passed to
#'   \code{\link[parallel]{mclapply}}.
#'
#' @return A data frame containing marker statistics for all one-vs-rest
#'   comparisons. The output includes the columns returned by
#'   \code{find_marker_in_expr_dich()}:
#'   \describe{
#'     \item{feature}{Feature (gene) name.}
#'     \item{pval}{Wilcoxon rank-sum test p-value.}
#'     \item{mean_diff}{Difference in mean expression among expressing cells
#'       (\code{group.1 - group.2}).}
#'     \item{log2fc}{Log2 fold-change of mean expression among expressing
#'       cells (\code{group.1 / group.2}).}
#'     \item{pct.1}{Fraction of cells expressing the feature in the target
#'       group.}
#'     \item{pct.2}{Fraction of cells expressing the feature in all remaining
#'       cells.}
#'     \item{padj}{Benjamini-Hochberg adjusted p-value.}
#'     \item{group.1}{Target group being tested.}
#'     \item{group.2}{Reference population (all other cells).}
#'   }
#'
#' @details
#' For each unique value in \code{group}, a binary grouping variable is
#' created consisting of:
#' \itemize{
#'   \item the focal group, and
#'   \item all remaining cells combined.
#' }
#'
#' Differential expression is then performed using
#' \code{find_marker_in_expr_dich()}, which evaluates expression differences
#' only among expressing cells (expression > 0).
#'
#' @seealso
#' \code{\link{find_marker_in_expr_dich}}
#'
#' @export
#'
#' @examples
#' # One-vs-rest marker analysis from a Seurat object
#' res <- find_marker_in_expr_all(
#'   obj = seu,
#'   group = "cell_type"
#' )
#'
#' # One-vs-rest marker analysis from a matrix
#' res <- find_marker_in_expr_all(
#'   obj = mat,
#'   group = cluster_ids
#' )
find_marker_in_expr_all <- function(obj,
                                    group,
                                    get_layer_args = list(),
                                    mc.cores = 4) {


  if (length(group) == 1 && !methods::is(obj, "Seurat")) {
    stop("when obj is a matrix, group has to be a vector.")
  } else if (length(group) == 1) {
    group <- obj@meta.data[[group]]
  }
  group <- as.character(group)

  if (methods::is(obj, "Seurat")) {
    obj <- Gmisc::fastDoCall(get_layer, args = c(list(obj = obj), get_layer_args))
  }

  if (ncol(obj) != length(group)) {
    stop("column number of obj not equal to length of group.")
  }

  out <- purrr::map_dfr(unique(group), function(x) {
    group_dich <- ifelse(group == x, x, paste0("not_", x))
    y <- find_marker_in_expr_dich(
      obj,
      group = group_dich,
      get_layer_args = get_layer_args,
      mc.cores = mc.cores
    )
    if (y$group.1[1] != x) {
      y <- dplyr::rename(y, "group.1" = group.2, "group.2" = group.1)
    }
    return(y)
  })

  return(out)
}

