#' Transfer cell type labels based on relative frequencies
#'
#' Creates a named character vector mapping cluster IDs to cell type labels
#' based on relative frequency thresholds. For each group, if any label exceeds
#' the upper cutoff, only those labels are used; otherwise, all labels exceeding
#' the lower cutoff are retained and collapsed.
#'
#' @param compbar_data_df A data frame or tibble containing cluster IDs,
#'   cell type labels, and relative frequencies.
#' @param cutoff_top Numeric scalar. Upper threshold for selecting dominant
#'   labels within each group.
#' @param cutoff_bottom Numeric scalar. Lower threshold used when no label
#'   exceeds `cutoff_top`.
#' @param collapse_str Character string used to collapse multiple labels
#'   into a single string.
#' @param col_freq Character string giving the name of the column containing
#'   relative frequencies.
#' @param col_group Character string giving the name of the grouping column
#'   (e.g., cluster or resolution ID).
#' @param col_label Character string giving the name of the label column
#'   (e.g., cell type annotations).
#'
#' @details
#' For each group defined by `col_group`:
#' \itemize{s
#'   \item If any value in `col_freq` exceeds `cutoff_top`, only labels
#'   associated with those values are retained.
#'   \item Otherwise, all labels with values exceeding `cutoff_bottom`
#'   are retained.
#' }
#' When multiple labels are retained, they are concatenated using
#' `collapse_str`.
#'
#' @return
#' A named character vector where names correspond to the unique values
#' of `col_group` and values are the assigned (possibly collapsed) labels.
#'
#' @examples
#' \dontrun{
#' compbar_main <- scexpr::composition_barplot(test,
#'                                             x_cat = "SCT_snn_res.2",
#'                                             fill_cat = "celltype__labels__single__main")
#' compbar_main_data <- compbar_main[["data"]]
#' transfer_labels2(compbar_main_data)
#' }
#'
#' @export
transfer_labels2 <- function(compbar_data_df,
                             cutoff_top = 0.75,
                             cutoff_bottom = 0.2,
                             collapse_str = "; ",
                             col_freq = "rel_x",
                             col_group = names(compbar_data_df)[1],
                             col_label = names(compbar_data_df)[2]) {

  out <- compbar_data_df |>
    dplyr::group_by(.data[[col_group]]) |>
    dplyr::summarise(
      label = {
        freq  <- .data[[col_freq]]
        label <- .data[[col_label]]

        if (any(freq > cutoff_top)) {
          paste(label[freq > cutoff_top], collapse = collapse_str)
        } else {
          paste(label[freq > cutoff_bottom], collapse = collapse_str)
        }
      },
      .groups = "drop"
    )

  return(stats::setNames(out$label, as.character(out[[col_group]])))

}
