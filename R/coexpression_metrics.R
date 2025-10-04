#' Measures of coexpression for a set of features
#'
#' Metrics:
#' sum = number of expressed genes per cell.
#' rel = fraction of genes expressed per cell.
#' abs = total (summed) expression
#' norm = sum of rescaled [0,3] expression, max is 3*length(features), rounded
#' lognorm = log2 of norm (dapening extreme values)
#' Check out: UCell::AddModuleScore_UCell
#' Check out: scexpr::fgsea_groupwise
#'
#' @param obj Seurat object
#' @param features gene features
#' @param assay assay in obj
#' @param layer layer in assay
#' @param nameprefix prefix to metrics
#'
#' @returns data frame
#' @export
#'
#' @examples
#'\dontrun{
#' coexpr <- coexpression_metrics(obj = so,
#'                                features = c("NKG7", "PRF1", "KLRG1"),
#'                                nameprefix = "killer")
#' so <- Seurat::AddMetaData(so, coexpr)
#' }
coexpression_metrics <- function(obj,
                                 features,
                                 assay = "RNA",
                                 layer = "data",
                                 nameprefix = "coexpr",
                                 group = NULL) {

  # weights?

  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  lay <- scexpr::get_layer(
    obj = obj,
    assay = assay,
    layer = layer,
    features = features
  )

  coexpr_sum <- Matrix::colSums(lay>0)
  coexpr_rel <- coexpr_sum/length(features)
  coexpr_abs <- Matrix::colSums(lay)
  coexpr_scale <- Matrix::rowSums(brathering::scale2(as.matrix(lay), 0, 3, margin = 1))
  coexpr_scale2 <- log2(coexpr_scale+1)

  df <- data.frame(
    coexpr_sum,
    coexpr_rel,
    coexpr_abs,
    round(coexpr_scale),
    coexpr_scale2,
    row.names = names(coexpr_sum)
  )

  names(df) <- paste0(nameprefix, "_", c("sum", "rel", "abs", "norm", "lognorm"))

  if (!is.null(group)) {
    group <- purrr::map_dfr(split(df, so@meta.data[[group]]), ~Matrix::colMeans(.x), .id = group)
  }

  # df <- df |>
  #   dplyr::mutate(dplyr::across(dplyr::everything(), ~dplyr::min_rank(.x), .names = "rank_{col}")) |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(avg_rank = mean(dplyr::c_across(dplyr::starts_with("rank_")))) |>
  #   dplyr::ungroup() |>
  #   dplyr::select(-dplyr::starts_with("rank")) |>
  #   dplyr::mutate(rank = dplyr::min_rank(-avg_rank)) |>
  #   dplyr::select(-avg_rank)

  return(list(bycell = df, grouped = group))
}
