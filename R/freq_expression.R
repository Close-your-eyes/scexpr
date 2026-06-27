#' Calculate the fraction of cells expressing one or more features
#'
#' Computes the fraction of cells with non-zero expression for one or more
#' features across groups defined by a metadata column. Optionally, cells can be
#' further split by an additional metadata variable before summarization.
#'
#' @param obj A Seurat object or list of Seurat objects.
#'
#' @param features Character vector of feature (gene) names whose expression
#' frequency should be calculated.
#'
#' @param meta_col Metadata column used to define the groups for which
#' expression frequencies are computed.
#'
#' @param assay Assay from which expression values are extracted.
#'
#' @param cells Optional vector of cell names to subset prior to computing
#' expression frequencies. If `NULL`, all cells are used.
#'
#' @param split_feature Optional metadata column used to split the data before
#' summarization. Frequencies are calculated independently for each level of
#' this variable.
#'
#' @param format_pct Logical. If `TRUE`, add formatted character columns for the
#' expression frequencies (e.g. `"25 %"` and `"0.25"`).
#'
#' @return A data frame containing one row per feature (or feature split) and
#' group. The output includes:
#' \describe{
#'   \item{feature_split}{Feature identifier (or feature/split combination).}
#'   \item{meta_col}{Grouping variable supplied by `meta_col`.}
#'   \item{SO.split}{Identifier of the originating Seurat object when multiple
#'   objects are supplied.}
#'   \item{max.feat.expr}{Maximum observed expression of the feature within the
#'   subset.}
#'   \item{pct.expr}{Fraction of cells with expression greater than zero.}
#'   \item{pct.expr.adjust.pct}{Formatted percentage (only when
#'   `format_pct = TRUE`).}
#'   \item{pct.expr.adjust}{Formatted fraction (only when
#'   `format_pct = TRUE`).}
#' }
#'
#' @details
#' Expression is considered positive when the extracted expression value is
#' greater than zero. Data are retrieved using `get_data()` with
#' `layer = "data"`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Fraction of cells expressing selected genes by cell type
#' freq_expression(
#'   obj = pbmc,
#'   features = c("MS4A1", "CD3D", "LYZ"),
#'   meta_col = "celltype"
#' )
#'
#' ## Split by sample
#' freq_expression(
#'   obj = pbmc,
#'   features = "MS4A1",
#'   meta_col = "celltype",
#'   split_feature = "sample"
#' )
#' }
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
