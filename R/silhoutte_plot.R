#' Make silhouette plot for clustering
#'
#' From Seurat object. Data frame added to misc slot.
#'
#' @param obj seurat object
#' @param clustering column in meta.data
#' @param reduction reduction name
#' @param theme ggplot theme
#' @param col_pal color and fill palette
#'
#' @returns ggplot and data frame to misc slot
#' @export
#'
#' @examples
silhoutte_plot <- function(obj,
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
  obj@misc[[paste0(clustering, "_silhoutte")]] <- sil_df

  ggplot2::ggplot(sil_df, ggplot2::aes(x = row, y = sil_width)) +
    ggplot2::geom_col(ggplot2::aes(color = cluster, fill = cluster)) +
    ggplot2::geom_hline(yintercept = mean(sil_df$sil_width), linetype = "dashed") +
    theme +
    ggplot2::scale_color_manual(values = col_pal) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::labs(y = "silhouette width", color = clustering, fill = clustering)

}
