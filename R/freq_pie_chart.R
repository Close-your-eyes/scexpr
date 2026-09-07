#' Frequency pie chart from metadata
#'
#' Creates a pie chart showing the frequency distribution of a metadata column
#' from a Seurat object or a data frame. This function is a wrapper around
#' \code{brathering::piechart()} with convenient defaults for single-cell workflows.
#'
#' @param SO A Seurat object or a data frame. If a Seurat object is provided,
#'   its \code{meta.data} slot will be used.
#' @param meta_col Character scalar. Column name in metadata used to compute frequencies.
#' @param order Optional character vector specifying the order of categories in the plot.
#' @param fill A vector of fill colors or a palette function. Defaults to
#'   \code{colrr::col_pal("custom")}.
#' @param fill_na Color used for missing values. Defaults to \code{"grey50"}.
#' @param color Border color of pie slices. Defaults to \code{"white"}.
#' @param radius_inside Numeric. Inner radius of the pie chart (for donut charts).
#' @param label_outside Character. Type of labels outside the pie. One of
#'   \code{"none"}, \code{"abs"}, or \code{"rel"}.
#' @param label_inside Character. Type of labels inside the pie. One of
#'   \code{"rel"}, \code{"abs"}, or \code{"none"}.
#' @param label_rel_cutoff Numeric. Minimum relative fraction required to display labels.
#' @param label_size Numeric. Text size of labels.
#' @param label_radius_inside Numeric. Radial position for inside labels.
#' @param label_radius_outside Numeric. Radial position for outside labels.
#' @param label_angle_inside Inside-label orientation. Supply
#'   `"radial_readable"`, `"radial"`, `"tangent_readable"`, or `"tangent"`;
#'   alternatively, supply one or more numeric angles in degrees as interpreted
#'   by [ggplot2::geom_text()]. See **Label angles**.
#' @param label_angle_outside Outside-label orientation. Accepts the same values
#'   as `label_angle_inside`.
#' @param label_overlap Character. Strategy for overlapping labels. One of
#'   \code{"ignore"}, \code{"alternate"}, or \code{"outside"}.
#' @param overlap_outside_radius Numeric. Radius used when resolving label overlap outside.
#' @param label_rel_pct Logical. Whether to display relative values as percentages.
#' @param label_rel_dec Integer. Number of decimal places for relative labels.
#' @param legend_title Optional character. Title of the legend.
#' @param theme A ggplot2 theme object. Defaults to \code{ggplot2::theme_classic()}.
#' @param theme_args A named list of additional theme modifications applied via
#'   \code{ggplot2::theme()}.
#' @param col_pal_args List of arguments passed to \code{colrr::col_pal()}.
#' @param split split seurat object and produce multiple pie charts?
#' @param label_color_inside label color
#' @param label_color_outside label color
#'
#' @section Label angles:
#' Slice angles used by [ggforce::geom_arc_bar()] are measured in radians,
#' beginning at 12 o'clock and increasing clockwise. Text angles passed to
#' [ggplot2::geom_text()] are measured in degrees, with positive values rotating
#' counterclockwise.
#'
#' The named label-angle modes are:
#'
#' * `"radial"`: align the text baseline with the slice's radial direction.
#' * `"radial_readable"`: use radial alignment and flip labels on the opposite
#'   half of the circle to keep them upright.
#' * `"tangent"`: align the text baseline with the circle's tangent, making it
#'   orthogonal to the slice's radial direction.
#' * `"tangent_readable"`: use tangential alignment and flip equivalent
#'   orientations to keep labels upright.
#'
#' Modes without the `_readable` suffix preserve the uncorrected geometric
#' orientation, so some labels may appear upside down. Numeric angles bypass
#' automatic calculation; for example, `0` makes labels horizontal.
#'
#' @return A ggplot2 object representing the pie chart.
#'
#' @details
#' This function automatically installs the required packages
#' \code{colrr} and \code{brathering} if they are not already available.
#'
#' Internally, it extracts the selected metadata column and forwards all
#' arguments to \code{brathering::piechart()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot <- freq_pie_chart(so, "celltype")
#' }
freq_pie_chart <- function(SO,
                           meta_col,
                           split = NULL,
                           order = NULL,
                           fill = "..auto..",
                           fill_na = "grey50",
                           color = "white",
                           label_color_inside = "..auto..",
                           label_color_outside = "..auto..",
                           radius_inside = 0.3,
                           label_outside = c("none", "abs", "rel"),
                           label_inside = c("rel", "abs", "none"),
                           label_rel_cutoff = 0,
                           label_size = 5,
                           label_radius_inside = 0.75,
                           label_radius_outside = 1.1,
                           label_angle_inside = "radial_readable",
                           label_angle_outside =  "radial_readable",
                           label_overlap = c("ignore", "alternate", "outside"),
                           overlap_outside_radius = 1.1,
                           label_rel_pct = F,
                           label_rel_dec = 2,
                           legend_title = NULL,
                           theme = ggplot2::theme_classic(),
                           theme_args = list(panel.grid = ggplot2::element_blank(),
                                             axis.title.x = ggplot2::element_blank(),
                                             axis.title.y= ggplot2::element_blank(),
                                             axis.text.x = ggplot2::element_blank(),
                                             axis.text.y = ggplot2::element_blank(),
                                             axis.ticks.x = ggplot2::element_blank(),
                                             axis.ticks.y = ggplot2::element_blank(),
                                             axis.line.x = ggplot2::element_blank(),
                                             axis.line.y = ggplot2::element_blank()),
                           col_pal_args = list(missing_fct_to_na = T)) {

  #  panel.background = ggplot2::element_rect(fill = "white")
  scexpr:::.ensure_packages(c("colrr", "brathering",))


  if (fill[1] == "..auto..") {
    if (methods::is(SO, "Seurat")) {
      if ("metacolors" %in% names(SO@misc) && is.list(SO@misc[["metacolors"]])) {
        if (meta_col %in% names(SO@misc[["metacolors"]])) {
          fill <- SO@misc[["metacolors"]][[meta_col]]
        } else {
          fill <- colrr::col_pal("custom")
        }
      } else {
        fill <- colrr::col_pal("custom")
      }
    } else {
      fill <- colrr::col_pal("custom")
    }
  }

  if (methods::is(SO, "Seurat")) {
    if (!is.null(split)) {
      SO <- purrr::map(Seurat::SplitObject(SO, split.by = split), ~.x@meta.data)
    } else {
      SO <- list(SO@meta.data)
    }
  } else {
    if (is.data.frame(SO)) {
      SO <- list(SO)
    }
  }

  # brathering::
  pies <- purrr::map(SO, ~brathering::piechart(x = .x[[meta_col]],
                                               order = order,
                                               fill = fill,
                                               fill_na = fill_na,
                                               color = color,
                                               label_color_inside = label_color_inside,
                                               label_color_outside = label_color_outside,
                                               radius_inside = radius_inside,
                                               label_outside = label_outside,
                                               label_inside = label_inside,
                                               label_rel_cutoff = label_rel_cutoff,
                                               label_size = label_size,
                                               label_radius_inside = label_radius_inside,
                                               label_radius_outside = label_radius_outside,
                                               label_angle_inside = label_angle_inside,
                                               label_angle_outside = label_angle_outside,
                                               label_overlap = label_overlap,
                                               overlap_outside_radius = overlap_outside_radius,
                                               label_rel_pct = label_rel_pct,
                                               label_rel_dec = label_rel_dec,
                                               legend_title = legend_title,
                                               theme_args = theme_args,
                                               theme = theme,
                                               col_pal_args = col_pal_args))

  if (length(pies) == 1) {
    pies <- pies[[1]]
  }
  return(pies)
}
