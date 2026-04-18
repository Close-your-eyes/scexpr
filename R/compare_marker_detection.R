#' Compare Marker Detection Methods
#'
#' Compare markers identified by a global marker detection method and a
#' pairwise marker detection approach. The function returns markers that are
#' uniquely identified by the pairwise method but not by the global method.
#'
#' @param obj A Seurat object
#' @param meta_col A character string specifying the metadata column used to define groups/clusters.
#'
#' @return A \code{data.frame} containing markers detected only by the pairwise
#' method, after filtering by log-fold change, adjusted p-value, and percentage expression.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Identifies markers using a global method via \code{find_all_marker()}.
#'   \item Identifies markers using a pairwise comparison method via \code{find_markers_pairwise()}.
#'   \item Filters both marker sets using:
#'     \itemize{
#'       \item \code{avg_log2FC > 1}
#'       \item \code{padj < 1e-5}
#'       \item \code{pct_in > 20}
#'     }
#'   \item Returns markers present only in the pairwise results.
#' }
#'
#' @seealso \code{\link{find_all_marker}}, \code{\link{find_markers_pairwise}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with a Seurat object
#' result <- compare_marker_detection(
#'   obj = seurat_object,
#'   meta_col = "cell_type"
#' )
#'
#' head(result)
#' }
compare_marker_detection <- function(obj,
                                     meta_col = NULL) {

  ## which markers are additionally discovered by pairwise method

  marker1 <- find_all_marker(obj, meta_col = meta_col)
  marker2 <- find_markers_pairwise(obj, group = meta_col,
                                   mc.cores = 4, method = "scexpr", redundant = T)

  marker2 <- marker2 |>
    dplyr::filter(avg_log2FC>1 & padj<1e-5 & pct_in>20)
  marker1 <- marker1 |>
    dplyr::filter(avg_log2FC>1 & padj<1e-5 & pct_in>20)

  marker_add <- marker2 |> dplyr::anti_join(marker1, by = c("feature" = "feature", "ident.1" = "group"))

  return(list(all = marker1, pairwise = marker2, pairwise_add = marker_add))
}


