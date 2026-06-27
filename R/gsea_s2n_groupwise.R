#' Computes the signal-to-noise (S2N) statistic for every gene in a one-vs-rest
#' comparison for each group (e.g. cluster or cell type). For a given group
#' \eqn{k}, the score is calculated as
#'
#' \deqn{
#' S2N(g, k) = \frac{\mu_k - \mu_{rest}}{\sigma_k + \sigma_{rest} + \epsilon}
#' }
#'
#' where \eqn{\mu} and \eqn{\sigma} denote the mean and standard deviation of
#' gene expression within the target group and all remaining cells,
#' respectively. Positive scores indicate genes enriched in the target group,
#' whereas negative scores indicate higher expression outside the group.
#'
#' @param obj A Seurat object or a gene-by-cell expression matrix (dense or
#' sparse `Matrix`).
#'
#' @param group Group assignments. Either:
#' \itemize{
#'   \item a character scalar specifying a column in `obj@meta.data` (when
#'   `obj` is a Seurat object), or
#'   \item a vector of length `ncol(obj)` giving the group assignment for each
#'   cell.
#' }
#'
#' @param eps Small positive constant added to the denominator to prevent
#' division by zero.
#'
#' @param get_layer_args Named list of additional arguments passed to
#' `scexpr::get_layer()` when extracting the expression matrix from a Seurat
#' object.
#'
#' @return A numeric matrix with one row per gene and one column per group.
#' Matrix entries contain the S2N statistic for each gene in each one-vs-rest
#' comparison.
#'
#' @details
#' Standard deviations are computed independently within each group and its
#' complement using a row-wise second-moment estimator that supports both dense
#' and sparse matrices efficiently.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## From a Seurat object
#' s2n <- gsea_s2n_groupwise(
#'   obj = pbmc,
#'   group = "celltype"
#' )
#'
#' ## From an expression matrix
#' expr <- scexpr::get_layer(pbmc)
#' s2n <- gsea_s2n_groupwise(
#'   obj = expr,
#'   group = pbmc$celltype
#' )
#'
#' ## Top marker genes for one group
#' head(sort(s2n[, "B cells"], decreasing = TRUE))
#' }
gsea_s2n_groupwise <- function(obj,
                               group,
                               eps = 1e-6,
                               get_layer_args = list()) {

  if (length(group) == 1 && methods::is(obj, "Seurat") && group %in% names(obj@meta.data)) {
    group <- obj@meta.data[[group]]
  }

  if (methods::is(obj, "Seurat")) {
    obj <- Gmisc::fastDoCall(get_layer, args = c(list(obj = obj), get_layer_args))
  }

  stopifnot(ncol(obj) == length(group))
  if (!inherits(group, "factor")) group <- factor(group)
  K <- levels(group)

  # Helper: row-wise SD using second-moment trick (works for dense & sparse)
  row_sd <- function(x) {
    m  <- Matrix::rowMeans(x)
    m2 <- Matrix::rowMeans(x^2)
    sqrt(pmax(m2 - m^2, 0))
  }

  # Precompute per-cluster indexes
  idx_list <- lapply(K, function(k) which(group == k))
  names(idx_list) <- K

  # Allocate result
  S2N <- matrix(NA_real_, nrow = nrow(obj), ncol = length(K),
                dimnames = list(rownames(obj), K))

  n_total <- ncol(obj)

  for (k in K) {
    i1 <- idx_list[[k]]
    i2 <- setdiff(seq_len(n_total), i1)
    X1 <- obj[, i1, drop = FALSE]
    X2 <- obj[, i2, drop = FALSE]

    m1 <- Matrix::rowMeans(X1)
    m2 <- Matrix::rowMeans(X2)
    s1 <- row_sd(X1)
    s2 <- row_sd(X2)

    S2N[, k] <- (m1 - m2) / (s1 + s2 + eps)
  }

  return(S2N)
}
