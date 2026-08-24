#' Add Quantile-Based Feature Bins to Seurat Metadata
#'
#' Divides a feature into quantile-based intervals and adds the resulting
#' categorical variable to a Seurat object's metadata.
#'
#' @param so A Seurat object.
#' @param feature A single character string specifying the feature to bin,
#'   such as `"nCount_RNA"`.
#' @param probs A numeric vector of probabilities used to calculate the
#'   quantile breakpoints. Defaults to quintiles:
#'   `seq(0, 1, 0.2)`.
#'
#' @return The Seurat object with a new metadata column named
#'   `<feature>_bins`.
#'
#' @examples
#' \dontrun{
#' so <- add_feature_bins(
#'   so,
#'   feature = "nCount_RNA",
#'   probs = seq(0, 1, 0.2)
#' )
#' }
#'
#' @export
add_feature_bins <- function(
    so,
    feature,
    probs = seq(0, 1, 0.2)
) {

  if (!is.character(feature) || length(feature) != 1L || is.na(feature)) {
    stop("`feature` must be a single, non-missing character string.")
  }

  vals <- scexpr::get_data(
    so,
    feature,
    reduction = NULL
  )[[1]]

  if (!is.numeric(vals$feature)) {
    stop(
      sprintf(
        "`%s` must be numeric; it is <%s>.",
        feature,
        paste(class(vals$feature), collapse = "/")
      )
    )
  }

  breaks <- unique(stats::quantile(
    vals$feature,
    probs = probs,
    na.rm = TRUE
  ))

  vals <- vals |>
    dplyr::mutate(
      bins = cut(
        feature,
        breaks = breaks,
        include.lowest = TRUE,
        dig.lab = 15
      )
    ) |>
    dplyr::select(id, bins) |>
    tibble::remove_rownames() |>
    tibble::column_to_rownames("id")

  names(vals) <- paste0(feature, "_bins")

  Seurat::AddMetaData(so, vals)
}
