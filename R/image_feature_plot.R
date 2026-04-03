#' Plot Spatial Image with Optional Feature Overlay
#'
#' Generates a ggplot2-based visualization of spatial imaging data, optionally
#' overlaying cell segmentation polygons or cell-level points colored by a
#' selected feature. Supports background image rendering and flexible styling.
#'
#' @param obj seurat object
#' @param feature Character scalar. Name of the feature to visualize. If `NULL`,
#'   only structural elements (image or segmentation) are shown.
#' @param plot_img Logical. Whether to render the background image. Default is `TRUE`.
#' @param img_subset_fct Numeric. Downsampling factor for the image to improve
#'   plotting performance. Default is `4`.
#' @param img_col_conv_fun Function. Function used to convert image colors
#'   (e.g., \code{colrr::hex_to_grey}). Default is \code{colrr::hex_to_grey}.
#' @param plot_cell_seg Logical. Whether to plot cell segmentation polygons.
#'   If `FALSE`, cell-level points are used instead. Default is `TRUE`.
#' @param cell_seg_args List. Additional arguments passed to
#'   \code{ggplot2::geom_polygon()} for segmentation plotting (e.g., `color`,
#'   `linewidth`, `alpha`).
#' @param cell_point_args List. Additional arguments passed to
#'   \code{ggplot2::geom_point()} for point-based plotting (e.g., `size`, `alpha`).
#' @param get_data_args List. Additional arguments passed to the internal
#'   \code{get_data()} function (e.g., `layer = "data"`).
#' @param axes_hide Logical. Whether to hide axis labels, ticks, and text.
#'   Default is `TRUE`.
#' @param col_pal_c_args List. Arguments for continuous color palette generation
#'   using \code{colrr::make_col_pal()}. Must include `name`.
#' @param col_pal_d_args List. Arguments for discrete color palette generation
#'   using \code{colrr::make_col_pal()}. Must include `name`.
#' @param hull_df provide data frame of cell hulls directly
#'
#' @details
#' The function supports three main visualization modes:
#' \itemize{
#'   \item Background image only
#'   \item Cell segmentation polygons (optionally colored by feature)
#'   \item Cell-level points (if segmentation is disabled)
#' }
#'
#' When a feature is provided:
#' \itemize{
#'   \item Numeric features are visualized using a continuous color scale
#'   \item Categorical features are converted to factors and visualized using
#'   a discrete palette
#' }
#'
#' Multiple features are currently not supported.
#'
#' @return A \code{ggplot2} object representing the spatial visualization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with segmentation
#' p <- image_feature_plot(obj, feature = "GeneA")
#' print(p)
#'
#' # Plot without image background
#' p <- image_feature_plot(obj, feature = "GeneA", plot_img = FALSE)
#'
#' # Plot using points instead of segmentation
#' p <- image_feature_plot(obj, feature = "GeneA", plot_cell_seg = FALSE)
#'
#' # Customize segmentation appearance
#' p <- image_feature_plot(
#'   obj,
#'   feature = "GeneA",
#'   cell_seg_args = list(color = "black", linewidth = 0.2, alpha = 0.8)
#' )
#'
#' # Use custom color palette
#' p <- image_feature_plot(
#'   obj,
#'   feature = "GeneA",
#'   col_pal_c_args = list(name = "viridis")
#' )
#' }
image_feature_plot <- function(obj,
                               feature = NULL,
                               plot_img = T,
                               img_subset_fct = 4,
                               img_col_conv_fun = colrr::hex_to_grey,
                               plot_cell_seg = T,
                               hull_df = NULL,
                               cell_seg_args = list(color = NA,
                                                    linewidth = 0.1,
                                                    alpha = 0.7),
                               cell_point_args = list(size = 1,
                                                      alpha = 0.1),
                               get_data_args = list(layer = "data"),
                               axes_hide = T,
                               col_pal_c_args = list(name = "spectral", direction = -1),
                               col_pal_d_args = list(name = "custom"),
                               col_binary = F) {

  if (length(feature) > 1) {
    stop("only one feature at a time, yet.")
  }
  # not fully customizable yet
  # for spatial data as from nikos

  if (plot_img) {
    img_df <- get_img(
      obj,
      subset_factor = img_subset_fct,
      col_conv_fun = img_col_conv_fun
    )
  }

  if (plot_cell_seg && is.null(hull_df)) {
    hull_df <- get_cell_hulls(so)
  }

  if (!is.null(feature)) {
    if (plot_cell_seg) {
      feature_df <- Gmisc::fastDoCall(get_data, args = c(list(SO = obj,
                                                              feature = feature,
                                                              reduction = NULL,
                                                              try_df = T), get_data_args))
      hull_df <- dplyr::left_join(hull_df, feature_df, by = c("z" = "id"))
      if (col_binary && is.numeric(feature_df[[feature]])) {
        hull_df[[feature]] <- factor(ifelse(hull_df[[feature]]>0, "+", "-"), levels = c("-", "+"))
      }
    } else {
      feature_df <- Gmisc::fastDoCall(get_data, args = c(list(SO = obj,
                                                              feature = feature,
                                                              reduction = "spatial",
                                                              try_df = T), get_data_args)) |>
        dplyr::rename("x" = spatial_1, "y" = spatial_2)
      if (col_binary && is.numeric(feature_df[[feature]])) {
        feature_df[[feature]] <- factor(ifelse(feature_df[[feature]]>0, "+", "-"), levels = c("-", "+"))
      }
    }
  }

  if (plot_img) {

    plot <- ggplot2::ggplot(img_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_raster(ggplot2::aes(fill = color)) +
      ggplot2::scale_fill_identity()

    if (plot_cell_seg && !is.null(feature)) {

      plot <- plot +
        ggnewscale::new_scale_fill() +
        Gmisc::fastDoCall(ggplot2::geom_polygon,
                          args = c(cell_seg_args, list(data = hull_df,
                                                       mapping = ggplot2::aes(group = z, fill = !!rlang::sym(feature)))))
    } else if (plot_cell_seg && is.null(feature)) {

      plot <- plot +
        ggnewscale::new_scale_fill() +
        Gmisc::fastDoCall(ggplot2::geom_polygon,
                          args = c(cell_seg_args, list(data = hull_df,
                                                       mapping = ggplot2::aes(group = z))))

    } else if (!plot_cell_seg && !is.null(feature)) {

      plot <- plot +
        Gmisc::fastDoCall(ggplot2::geom_point,
                          args = c(cell_point_args, list(data = feature_df,
                                                         mapping = ggplot2::aes(color = !!rlang::sym(feature)))))

    }

  } else {

    if (plot_cell_seg && !is.null(feature)) {

      plot <- ggplot2::ggplot(hull_df, ggplot2::aes(x = x, y = y)) +
        Gmisc::fastDoCall(ggplot2::geom_polygon,
                          args = c(cell_seg_args, list(mapping = ggplot2::aes(group = z, fill = !!rlang::sym(feature)))))

    } else if (plot_cell_seg && is.null(feature)) {

      plot <- ggplot2::ggplot(hull_df, ggplot2::aes(x = x, y = y)) +
        Gmisc::fastDoCall(ggplot2::geom_polygon,
                          args = c(cell_seg_args, list(mapping = ggplot2::aes(group = z))))

    } else if (!plot_cell_seg && !is.null(feature)) {

      plot <- ggplot2::ggplot(feature_df, ggplot2::aes(x = x, y = y)) +
        Gmisc::fastDoCall(ggplot2::geom_point,
                          args = c(cell_point_args, list(mapping = ggplot2::aes(color = !!rlang::sym(feature)))))

    }
  }


  plot <- plot +
    colrr::theme_material() +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion()) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion())

  if (plot_cell_seg) {
    if (is.numeric(hull_df[[feature]])) {
      col.pal <- colrr::make_col_pal(col_vec = col_pal_c_args[["name"]],
                                     col_pal_args = col_pal_c_args[-which(names(col_pal_c_args) == "name")])
      plot <- plot + ggplot2::scale_fill_gradientn(colors = col.pal)
    } else {
      if (!is.factor(hull_df[[feature]])) {
        hull_df[[feature]] <- as.factor(hull_df[[feature]])
      }
      col.pal <- colrr::make_col_pal(col_vec = col_pal_d_args[["name"]],
                                     fct_lvls = levels(forcats::fct_drop(hull_df[[feature]])),
                                     missing_fct_to_na = ifelse("missing_fct_to_na" %in%
                                                                  names(col_pal_d_args), col_pal_d_args[["missing_fct_to_na"]], T), col_pal_args = col_pal_d_args[-which(names(col_pal_d_args) %in%
                                                                                                                                                                           c("name", "missing_fct_to_na"))])
      plot <- plot + ggplot2::scale_fill_manual(values = col.pal)
    }
  }

  if (axes_hide) {
    plot <- plot +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank())
  }

  return(plot)
}
