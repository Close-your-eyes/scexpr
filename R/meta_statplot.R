#' Plot average feature expression across metadata groups
#'
#' Computes average expression for selected features grouped by one metadata
#' column, joins those averages to a second metadata grouping column, and creates
#' faceted statistical plots for each feature.
#'
#' @param obj A Seurat object containing expression data and metadata.
#' @param group1 Character string giving the metadata column used to calculate
#'   average expression with `scexpr::avg_expression()`.
#' @param group2 Character string giving the metadata column used on the x-axis
#'   for statistical comparison.
#' @param features Character vector of feature names, typically genes, to average
#'   and plot.
#'
#' @returns A `ggplot` object showing average expression values by `group2`,
#'   faceted by feature.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' meta_statplot(
#'   obj = SO,
#'   group1 = "seurat_clusters",
#'   group2 = "orig.ident",
#'   features = c("CD3D", "MS4A1", "LYZ")
#' )
#' }
meta_statplot <- function(obj, group1, group2, features) {

  avgexpr <- scexpr::avg_expression(obj, group = group1, features = features)[[1]] |>
    as.data.frame() |>
    brathering::mat_to_df_long(values_to = "avg_expr", colnames_to = group1, rownames_to = "feature") |>
    dplyr::left_join(dplyr::distinct(obj@meta.data, !!rlang::sym(group1), !!rlang::sym(group2)), by = group1)
  plot <- ggplot2::ggplot(avgexpr, ggplot2::aes(!!rlang::sym(group2), avg_expr)) +
    ggplot2::geom_point() +
    colrr::theme_material(white = T, style = "prismy") +
    ggpubr::geom_pwc() +
    ggplot2::facet_wrap(vars(feature), scales = "free_y")

  return(plot)
}
