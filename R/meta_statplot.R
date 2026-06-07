#' Title
#'
#' @param obj
#' @param group1
#' @param group2
#' @param features
#'
#' @returns
#' @export
#'
#' @examples
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
