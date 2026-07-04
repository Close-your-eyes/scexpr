#' Find nearest-neighbor and shared-nearest-neighbor graphs
#'
#' Computes a k-nearest-neighbor (KNN) graph, and optionally a
#' shared-nearest-neighbor (SNN) graph, from an input matrix or distance matrix.
#' This is a modified neighbor-finding helper based on Seurat internals.
#' It always returns the ranked nearest neighbor graph.
#'
#' @param object A matrix-like object with cells in rows and features in columns,
#'   or a precomputed distance matrix if `distance.matrix = TRUE`. Row names are
#'   required and are used as cell names.
#' @param query Optional query matrix. If `NULL`, `object` is used as the query.
#' @param distance.matrix Logical. Whether `object` is a precomputed distance
#'   matrix. Default is `FALSE`.
#' @param k.param Integer. Number of nearest neighbors to compute. Default is `20`.
#' @param return.neighbor Logical. Whether to return the nearest-neighbor object
#'   directly instead of graph objects. Default is `FALSE`.
#' @param compute.SNN Logical. Whether to compute the SNN graph. Defaults to
#'   `!return.neighbor`.
#' @param prune.SNN Numeric. SNN pruning threshold. Edges with values less than
#'   or equal to this threshold are removed. Default is `1/15`.
#' @param nn.method Character. Nearest-neighbor search method. Default is
#'   `"annoy"`.
#' @param n.trees Integer. Number of trees to use when `nn.method = "annoy"`.
#'   Default is `50`.
#' @param annoy.metric Character. Distance metric used by Annoy. Default is
#'   `"euclidean"`.
#' @param nn.eps Numeric. Error bound for nearest-neighbor search. Default is `0`.
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#' @param l2.norm Logical. Whether to L2-normalize `object` and `query` before
#'   neighbor search. Default is `FALSE`.
#' @param cache.index Logical. Whether to cache the Annoy index. Default is
#'   `FALSE`.
#' @param index Optional precomputed nearest-neighbor index.
#' @param ... Additional arguments. Currently checked by `Seurat:::CheckDots`.
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{graphs}{A list containing the KNN graph as `nn` and, if requested,
#'   the SNN graph as `snn`.}
#'   \item{nn.ranked}{The nearest-neighbor object returned by Seurat's internal
#'   neighbor search helper.}
#' }
#'
#' If `return.neighbor = TRUE`, the nearest-neighbor object is returned directly
#' and graph construction is skipped.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- matrix(rnorm(1000), nrow = 100)
#' rownames(mat) <- paste0("cell_", seq_len(nrow(mat)))
#'
#' neighbors <- FindNeighbors2(mat, k.param = 20)
#' names(neighbors$graphs)
#' }
FindNeighbors2 <- function(
    object,
    query = NULL,
    distance.matrix = FALSE,
    k.param = 20,
    return.neighbor = FALSE,
    compute.SNN = !return.neighbor,
    prune.SNN = 1/15,
    nn.method = "annoy",
    n.trees = 50,
    annoy.metric = "euclidean",
    nn.eps = 0,
    verbose = TRUE,
    l2.norm = FALSE,
    cache.index = FALSE,
    index = NULL,
    ...) {

  Seurat:::CheckDots(...)
  if (is.null(x = dim(x = object))) {
    warning(
      "Object should have two dimensions, attempting to coerce to matrix",
      call. = FALSE
    )
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning(
      "k.param set larger than number of cells. Setting k.param to number of cells - 1.",
      call. = FALSE
    )
    k.param <- n.cells - 1
  }
  if (l2.norm) {
    object <- Seurat:::L2Norm(mat = object)
    query <- query %iff% Seurat:::L2Norm(mat = query)
  }
  query <- query %||% object
  # find the k-nearest neighbors for each single cell
  if (!distance.matrix) {
    if (verbose) {
      if (return.neighbor) {
        message("Computing nearest neighbors")
      } else {
        message("Computing nearest neighbor graph")
      }
    }
    nn.ranked_keep <- Seurat:::NNHelper(
      data = object,
      query = query,
      k = k.param,
      method = nn.method,
      n.trees = n.trees,
      searchtype = "standard",
      eps = nn.eps,
      metric = annoy.metric,
      cache.index = cache.index,
      index = index
    )
    if (return.neighbor) {
      if (compute.SNN) {
        warning("The SNN graph is not computed if return.neighbor is TRUE.", call. = FALSE)
      }
      return(nn.ranked_keep)
    }

    nn.ranked <- SeuratObject::Indices(object = nn.ranked_keep)
  } else {
    if (verbose) {
      message("Building SNN based on a provided distance matrix")
    }
    knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
    knd.mat <- knn.mat
    for (i in 1:n.cells) {
      knn.mat[i, ] <- order(object[i, ])[1:k.param]
      knd.mat[i, ] <- object[i, knn.mat[i, ]]
    }
    nn.ranked <- knn.mat[, 1:k.param]
  }
  # convert nn.ranked into a Graph
  j <- as.numeric(x = t(x = nn.ranked))
  i <- ((1:length(x = j)) - 1) %/% k.param + 1
  nn.matrix <- as(object = Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(x = object), nrow(x = object))), Class = "Graph")
  rownames(x = nn.matrix) <- rownames(x = object)
  colnames(x = nn.matrix) <- rownames(x = object)
  neighbor.graphs <- list(nn = nn.matrix)
  if (compute.SNN) {
    if (verbose) {
      message("Computing SNN")
    }
    snn.matrix <- Seurat:::ComputeSNN(
      nn_ranked = nn.ranked,
      prune = prune.SNN
    )
    rownames(x = snn.matrix) <- rownames(x = object)
    colnames(x = snn.matrix) <- rownames(x = object)
    snn.matrix <- SeuratObject::as.Graph(x = snn.matrix)
    neighbor.graphs[["snn"]] <- snn.matrix
  }
  return(list(graphs = neighbor.graphs, nn.ranked = nn.ranked_keep))
}
