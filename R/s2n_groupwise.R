#' Compute S2N per gene for each group (cluster) (one-vs-rest)
#'
#' S2N(g, k) = (mean_k - mean_rest) / (sd_k + sd_rest + eps)
#'
#' @param obj Seurat object or gene x cell matrix
#' @param group grouping column in obj meta.data or vector of length ncol(obj)
#' @param eps epsilon to avoid division by 0 in S2N formula
#' @param get_layer_args arguments to scexpr::get_layer
#'
#' @returns matrix
#' @export
#'
#' @examples
s2n_groupwise <- function(obj,
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
