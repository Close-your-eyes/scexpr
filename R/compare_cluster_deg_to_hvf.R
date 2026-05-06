#' Compare cluster marker genes to highly variable features (HVFs)
#'
#' Identifies cluster-specific marker genes from a Seurat object and compares them
#' to the set of highly variable features (HVFs). The function filters markers
#' based on log fold-change, expression percentage, and adjusted p-value, then
#' selects the top markers per group. It returns a combined data frame indicating
#' whether each feature is a marker and/or an HVF.
#'
#' @param obj A \code{Seurat} object containing expression data and metadata.
#' @param meta_col Character string specifying the metadata column used to define clusters.
#'   If \code{NULL}, a default clustering identity in \code{obj} is used.
#' @param avg_log2FC_cut Numeric. Minimum average log2 fold-change threshold for marker filtering.
#' @param pct_in_cut Numeric. Minimum percentage of cells expressing the gene in the cluster.
#' @param padj_cut Numeric. Maximum adjusted p-value threshold for marker filtering.
#' @param topn_per_group Integer. Number of top marker genes to retain per group based on AUC.
#' @param verbose Logical. Whether to print progress messages and summaries.
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{df}{A data frame containing all features with binary indicators:
#'     \code{hvf} (highly variable feature) and \code{marker} (selected marker gene).}
#'   \item{marker}{A data frame of filtered marker genes used in the comparison.}
#' }
#'
#' @details
#' The function first identifies marker genes using \code{find_all_marker}, then applies
#' filtering thresholds. It compares these markers with the set of highly variable features
#' extracted via \code{Seurat::VariableFeatures}. Features not overlapping between the two
#' sets are reported for diagnostic purposes.
#'
#' @export
#'
#' @examples
compare_cluster_deg_to_hvf <- function(obj,
                                       meta_col = NULL,
                                       avg_log2FC_cut = 1,
                                       pct_in_cut = 30,
                                       padj_cut = 1e-5,
                                       topn_per_group = 10,
                                       verbose = T) {

  rawmarker <- find_all_marker(obj, meta_col = meta_col, na_rm = T)
  marker <- rawmarker |>
    dplyr::filter(avg_log2FC>avg_log2FC_cut & pct_in>pct_in_cut & padj<padj_cut) |>
    dplyr::slice_max(auc, n = topn_per_group, by = group, with_ties = F)

  print(table(marker$group))
  if (verbose) {
    if (all(marker$group %in% rawmarker$group)) {
      message("all groups have markers.")
    } else {
      message("!!! not all groups have markers.")
    }
  }

  hvfdf <- data.frame(feature = Seurat::VariableFeatures(obj), hvf = 1)
  if (nrow(hvfdf) == 0) {
    stop("no hvf in obj.")
  } else {
    if (verbose) {
      message("nhvf: ", nrow(hvfdf))
    }
  }

  hvfmarkerdf <- hvfdf |>
    dplyr::full_join(marker |> dplyr::distinct(feature) |> dplyr::mutate(marker = 1), by = "feature") |>
    brathering::replace_na_all()

  marker_but_no_hvf <- hvfmarkerdf |>
    dplyr::filter(hvf == 0 & marker == 1)
  if (verbose) {
    message("n marker which are not hvf: ", nrow(marker_but_no_hvf))
  }

  return(list(df = hvfmarkerdf, marker = marker))
}
