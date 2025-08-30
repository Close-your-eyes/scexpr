#' Title
#'
#' @param SO Seurat object or data frame
#' @param meta_col column from meta data to use
#' @param order
#' @param fill
#' @param color
#' @param radius_inside
#' @param label_outside
#' @param label_inside
#' @param label_rel_cutoff
#' @param label_size
#' @param label_radius_inside
#' @param label_radius_outside
#' @param label_angle_inside
#' @param label_angle_outside
#' @param label_overlap
#' @param overlap_outside_radius
#' @param label_rel_pct
#' @param label_rel_dec
#' @param legend_title
#' @param theme_args
#'
#' @return
#' @export
#'
#' @examples
freq_pie_chart <- function(SO,
                           meta_col,
                           order = T,
                           fill = colrr::col_pal("custom"),
                           fill_na = "grey50",
                           color = "white",
                           radius_inside = 0.3,
                           label_outside = c("none", "abs", "rel"),
                           label_inside = c("rel", "abs", "none"),
                           label_rel_cutoff = 0,
                           label_size = 5,
                           label_radius_inside = 0.75,
                           label_radius_outside = 1.1,
                           label_angle_inside = NULL, # circle or numeric
                           label_angle_outside = NULL, # circle or numeric
                           label_overlap = c("ignore", "alternate", "outside"),
                           overlap_outside_radius = 1.1,
                           label_rel_pct = F,
                           label_rel_dec = 2,
                           legend_title = NULL,
                           theme_args = list(panel.grid = ggplot2::element_blank(),
                                             axis.title = ggplot2::element_blank(),
                                             axis.text = ggplot2::element_blank(),
                                             axis.ticks = ggplot2::element_blank()),
                           col_pal_args = list(missing_fct_to_na = T)) {
  if (methods::is(SO, "Seurat")) {
    SO <- SO@meta.data
  }

  brathering::piechart(x = SO[[meta_col]],
                       order = order,
                       fill = fill,
                       fill_na = fill_na,
                       color = color,
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
                       col_pal_args = col_pal_args)
}
