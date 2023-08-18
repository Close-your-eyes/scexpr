#' Title
#'
#' @param gsea_data_df
#' @param umap_args
#' @param FindNeighbors_args
#' @param FindClusters_args
#'
#' @return
#' @export
#'
#' @examples
cluster_gene_sets_by_leadingEdge <- function(gsea_data_df,
                                             umap_args = list(),
                                             FindNeighbors_args = list(),
                                             FindClusters_args = list(resolution = 0.6),
                                             seed = 42) {

  LE_genes <- unique(unlist(gsea_data_df$leadingEdge))
  LE_elements <- lapply(gsea_data_df$leadingEdge, function(x) LE_genes %in% x)

  # chunk-wise rbind is faster !!!
  # make that a separate function somewhen
  while (length(LE_elements) > 20) {
    LE_elements <- purrr::map(split(c(1:length(LE_elements)), ceiling(seq_along(c(1:length(LE_elements)))/10)), function(x) purrr::reduce(LE_elements[x], rbind))
  }
  LE_elements <- Reduce(rbind, LE_elements)
  LE_elements <- apply(LE_elements, 2, function(x) as.integer(as.logical(x)))
  colnames(LE_elements) <- LE_genes
  rownames(LE_elements) <- gsea_data_df$pathway
  # find groups of gene set based on leading edges
  set.seed(seed)
  dr <- Gmisc::fastDoCall(what = uwot::umap, args = c(umap_args, list(X = LE_elements)))
  colnames(dr) <- c("UMAP_1", "UMAP_2")
  set.seed(seed)
  snn <- Gmisc::fastDoCall(what = Seurat::FindNeighbors, args = c(FindNeighbors_args, list(object = LE_elements)))
  set.seed(seed)
  clust <- Gmisc::fastDoCall(what = Seurat::FindClusters, args = c(FindClusters_args, list(object = snn$snn)))

  dr <- as.data.frame(dr)
  dr$cluster <- clust[,1]
  dr$pathway <- rownames(dr)

  gsea_data_df <- dplyr::left_join(gsea_data_df, dr, by = "pathway")

  # which genes are shared by clusters
  LE_elements_df <- as.data.frame(LE_elements)
  LE_elements_df$pathway <- rownames(LE_elements_df)
  LE_elements_df <- tidyr::pivot_longer(LE_elements_df, cols = -pathway, values_to = "in_pathway", names_to = "gene")
  LE_elements_df$cluster <- clust[LE_elements_df$pathway, 1]
  LE_elements_df$in_pathway <- as.logical(LE_elements_df$in_pathway)
  LE_elements_df <-
    LE_elements_df %>%
    dplyr::group_by(cluster, gene) %>%
    dplyr::summarise(pathway_freq = sum(in_pathway)/dplyr::n(), .groups = "drop") %>%
    as.data.frame()

  return(list(data = gsea_data_df, LE_elements = LE_elements, gene_freq_by_cluster = LE_elements_df))

}
