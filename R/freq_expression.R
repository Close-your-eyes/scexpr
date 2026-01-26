#' Title
#'
#' @param obj
#' @param features
#' @param meta_col
#' @param assay
#' @param cells
#' @param split_feature
#' @param format_pct
#'
#' @returns
#' @export
#'
#' @examples
freq_expression <- function(obj,
                            features,
                            meta_col,
                            assay = "RNA",
                            cells = NULL,
                            split_feature = NULL,
                            format_pct = F) {

  data <- get_data(obj,
                   feature = features,
                   reduction = NULL,
                   assay = assay,
                   layer = "data",
                   cells = cells,
                   meta_col = meta_col,
                   split_feature = split_feature)
  data <- dplyr::bind_rows(data, .id = "feature_split")

  stat <-
    data %>%
    dplyr::mutate(max.feat.expr = max(feature), .by = feature_split) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("feature_split", meta_col, "SO.split", "max.feat.expr")))) %>%
    dplyr::summarise(pct.expr = sum(feature > 0)/dplyr::n(), .groups = "drop")

  if (format_pct) {
    stat <- stat |>
      dplyr::mutate(pct.expr.adjust.pct = dplyr::case_when(pct.expr == 0 ~ "0 %",
                                                           pct.expr > 0 & pct.expr < 0.01 ~ "> 1 %",
                                                           pct.expr >= 0.01 ~ paste0(round(pct.expr*100, expr.freq.decimals), " %"))) %>%
      dplyr::mutate(pct.expr.adjust = dplyr::case_when(pct.expr == 0 ~ "0",
                                                       pct.expr > 0 & pct.expr < 0.01 ~ "> 0.01",
                                                       pct.expr >= 0.01 ~ as.character(round(pct.expr, expr.freq.decimals))))
  }
  return(stat)
}
