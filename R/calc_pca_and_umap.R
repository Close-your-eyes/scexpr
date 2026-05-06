#' Compute PCA and UMAP embeddings
#'
#' Performs principal component analysis (PCA) on the input data and computes
#' a UMAP embedding based on the PCA coordinates.
#'
#' @param x A numeric matrix or data frame containing expression values
#'   (rows = features, columns = samples or cells).
#' @param npc Integer. Number of principal components to compute. Default is 10.
#' @param rows_to Character. Optional column name to store row names when
#'   converting to a data frame.
#' @param umap_args List. Additional arguments passed to
#'   \code{fcexpr::ff_calc_umap_tsne} for UMAP computation. Defaults include
#'   \code{n_neighbors = 10}, \code{metric = "euclidean"},
#'   \code{verbose = TRUE}, and \code{scale = FALSE}.
#'
#' @returns A data frame containing PCA coordinates and UMAP embeddings.
#'   If \code{rows_to} is provided, row names are added as a column.
#'
#' @details
#' PCA is computed using \code{FactoMineR::PCA}, and the resulting coordinates
#' are passed to \code{fcexpr::ff_calc_umap_tsne} to generate UMAP embeddings.
#'
#' @export
#' @examples
calc_pca_and_umap <- function(x,
                              npc = 10,
                              rows_to = NULL,
                              umap_args = list(n_neighbors = 10,
                                               metric = "euclidean",
                                               verbose = T,
                                               scale = F)) {
  y <- factoextra::get_pca_ind(FactoMineR::PCA(x, ncp = npc, graph = F))$coord
  um <- fcexpr::ff_calc_umap_tsne(exprs = y, fun_args = umap_args)
  df <- as.data.frame(cbind(y, um))

  if (!is.null(rows_to)) {
    df <- df |> tibble::rownames_to_column(rows_to)
  }
  return(df)
}
