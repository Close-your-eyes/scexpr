#' @param query Matrix of data to query against object. If missing, defaults to
#' object.
#' @param distance.matrix Boolean value of whether the provided matrix is a
#' distance matrix; note, for objects of class \code{dist}, this parameter will
#' be set automatically
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param return.neighbor Return result as \code{\link[SeuratObject]{Neighbor}} object. Not
#' used with distance matrix input.
#' @param compute.SNN also compute the shared nearest neighbor graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param nn.method Method for nearest neighbor finding. Options include: rann,
#' annoy
#' @param annoy.metric Distance metric for annoy. Options include: euclidean,
#' cosine, manhattan, and hamming
#' @param n.trees More trees gives higher precision when using annoy approximate
#' nearest neighbor search
#' @param nn.eps Error bound when performing nearest neighbor search using RANN;
#' default of 0.0 implies exact nearest neighbor search
#' @param verbose Whether or not to print output to the console
#' @param l2.norm Take L2Norm of the data
#' @param cache.index Include cached index in returned Neighbor object
#' (only relevant if return.neighbor = TRUE)
#' @param index Precomputed index. Useful if querying new data against existing
#' index to avoid recomputing.
#'
#' @importFrom RANN nn2
#' @importFrom methods as
#'
#' @export
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
