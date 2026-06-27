#' Make a silhouette plot for clustering
#'
#' Computes silhouette widths for clusters in a Seurat object using cell embeddings
#' from a specified dimensionality reduction, and returns both the silhouette data
#' frame and a ggplot object.
#'
#' @param obj A Seurat object.
#' @param clustering Character string. Column name in `obj@meta.data` containing
#'   cluster assignments.
#' @param reduction Character string. Name of the dimensionality reduction to use.
#'   Default is `"pca"`.
#' @param theme A ggplot2 theme. Default is `colrr::theme_material(white = TRUE)`.
#' @param col_pal Named or unnamed vector of colors used for cluster outlines and fills.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{data}{A data frame containing silhouette widths and cluster assignments.}
#'   \item{plot}{A ggplot2 silhouette plot.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' silhouette_plot(so, clustering = "celltype")
#' }
silhouette_plot <- function(obj,
                           clustering,
                           reduction = "pca",
                           theme = colrr::theme_material(white = T),
                           col_pal = colrr::col_pal("custom")) {

  sil <- cluster::silhouette(x = as.numeric(as.character(obj@meta.data[[clustering]])),
                             dist = stats::dist(obj@reductions[[reduction]]@cell.embeddings))
  # factoextra::fviz_silhouette(sil)
  sil_df <- as.data.frame(sil) |>
    dplyr::arrange(cluster, dplyr::desc(sil_width)) |>
    dplyr::mutate(row = dplyr::row_number()) |>
    dplyr::mutate(cluster = as.character(cluster))
  rownames(sil_df) <- rownames(obj@meta.data)
  #obj@misc[[paste0(clustering, "_silhoutte")]] <- sil_df

  plot <- ggplot2::ggplot(sil_df, ggplot2::aes(x = row, y = sil_width)) +
    ggplot2::geom_col(ggplot2::aes(color = cluster, fill = cluster)) +
    ggplot2::geom_hline(yintercept = mean(sil_df$sil_width), linetype = "dashed") +
    theme +
    ggplot2::scale_color_manual(values = col_pal) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(y = "silhouette width", color = clustering, fill = clustering)

  return(list(data = sil_df, plot = plot))
}
