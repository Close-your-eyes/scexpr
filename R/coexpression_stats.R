#' Inspect reference feature expression in feature coexpression groups
#'
#' Calculates coexpression groups for all feature combinations of length two or
#' greater, then summarizes and plots expression of a reference feature across
#' those coexpression groups.
#'
#' @param obj A Seurat object.
#' @param features Character vector of features used to define coexpression
#'   groups.
#' @param ref_feature Reference feature whose expression is inspected across
#'   coexpression groups.
#' @param metric Name of the metric returned by `scexpr::coexpression_metrics()`
#'   to use for grouping and plotting. Defaults to `"coexpr_sum"`.
#'
#' @return
#' A named list with three elements:
#' \describe{
#'   \item{`df`}{Cell-level data containing coexpression metrics and reference
#'   feature expression.}
#'   \item{`df2`}{Summary table with coexpression frequencies and reference
#'   feature expression frequencies per group.}
#'   \item{`plot`}{A faceted `ggplot2` object showing reference feature
#'   expression across coexpression metric groups.}
#' }
#'
#' @details
#' Feature combinations are generated with `brathering::combnn()` using all
#' unique input features and a minimum combination length of two. For each
#' combination, `scexpr::coexpression_metrics()` is called and the selected
#' `metric` is extracted.
#'
#' The returned plot shows individual cells as jittered points. Text labels
#' indicate the fraction of cells in each coexpression group and the fraction
#' expressing `ref_feature`.
#'
#' @examples
#' \dontrun{
#' coexpression_stats(
#'   obj = so,
#'   features = c("CD53", "CD70", "CD27"),
#'   ref_feature = "KLF2"
#' )
#' }
#'
#' @export
coexpression_stats <- function(obj,
                               features,
                               ref_feature,
                               metric = "coexpr_sum") {

  combs <- purrr::list_flatten(
    brathering::combnn(unique(features), min_len = 2, return_numeric = F),
    name_spec = "{inner}"
  )
  df <- purrr::map_dfr(combs,
                       ~scexpr::coexpression_metrics(obj, features = .x)[[1]] |>
                         tibble::rownames_to_column("id"), .id = "features") |>
    dplyr::select(id, features, !!rlang::sym(metric)) |>
    tidyr::pivot_wider(names_from = features, values_from = !!rlang::sym(metric)) |>
    dplyr::left_join(scexpr::get_data(obj, ref_feature, reduction = NULL, try_df = T)[,c(1,2)], by = "id") |>
    tidyr::pivot_longer(cols = names(combs),
                        names_to = "features",
                        values_to = metric)

  feature_coexpr_freq <- dplyr::summarise(df,
                                          coexpr_pct = round(dplyr::n()/ncol(obj), 2),
                                          .by = c(features, !!rlang::sym(metric)))
  ref_feature_expr_freq <- dplyr::summarise(df,
                                            !!rlang::sym(ref_feature) := round(mean(!!rlang::sym(ref_feature)>0), 2),
                                            .by = c(features, !!rlang::sym(metric)))
  df2 <- dplyr::left_join(feature_coexpr_freq, ref_feature_expr_freq, by = dplyr::join_by(features, coexpr_sum))

  coexpr_plot <- ggplot2::ggplot(df, ggplot2::aes(x = !!rlang::sym(metric), y = !!rlang::sym(ref_feature))) +
    ggplot2::geom_jitter(width = 0.1) +
    ggplot2::geom_text(data = df2, y = max(df[[ref_feature]])*1.3, ggplot2::aes(label = coexpr_pct)) +
    ggplot2::geom_text(data = df2, y = max(df[[ref_feature]])*1.15, ggplot2::aes(label = !!rlang::sym(ref_feature)), color = "tomato2") +
    # ylim(c(NA, 4.2)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.5))) +
    colrr::theme_material(white = T) +
    ggplot2::facet_wrap(ggplot2::vars(features), axes = "all", axis.labels = "margins")

  return(list(df = df, df2 = df2, plot = coexpr_plot))
}
