#' Plot overview of different nHVF
#'
#' Plots mean-variance plots with different cutoffs. Variable features need
#' to be calculated in obj before calling this function.
#'
#' @param obj seurat object
#' @param n_varfeat number of variable deature to visualize (cutoffs)
#'
#' @returns patchwork of ggplots
#' @export
#'
#' @examples
varfeat_plot <- function(obj, n_varfeat = seq(200, 2000, 200)) {

  hvfdat <- Seurat::HVFInfo(obj)
  names(n_varfeat) <- n_varfeat
  if (all(c("mean", "variance.standardized") %in% names(hvfdat))) {
    xvar <- "mean"
    yvar <- "variance.standardized"
  } else if (all(c("gmean", "residual_variance") %in% names(hvfdat))) {
    xvar <- "gmean"
    yvar <- "residual_variance"
  } else {
    message("varfeat_plot: cannot determine what to plot.")
    return(NULL)
  }
  varfeatplots <- purrr::map(n_varfeat, function(nvarfeat) {
    ggplot2::ggplot(hvfdat, ggplot2::aes(!!rlang::sym(xvar), !!rlang::sym(yvar))) +
      ggplot2::geom_point(color = "white", size = 0.2) +
      colrr::theme_material() +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_text(margin = ggplot2::margin(2,0,2,0, unit = "pt")),
                     plot.margin = grid::unit(c(1,1,1,1), "pt"),
                     panel.spacing = grid::unit(2, "pt")) +
      ggplot2::geom_point(data = dplyr::slice_max(hvfdat, order_by = !!rlang::sym(yvar), n = nvarfeat), color = "tomato2", size = 0.2) +
      ggplot2::facet_wrap(nvarfeat)
  })
  varfeatplots <- c(list(density = ggplot2::ggplot(hvfdat, ggplot2::aes(!!rlang::sym(xvar), !!rlang::sym(yvar))) +
                           ggpointdensity::geom_pointdensity(size = 0.2) +
                           colrr::theme_material() +
                           ggplot2::scale_x_log10() +
                           ggplot2::scale_y_log10() +
                           ggplot2::theme(#axis.text.x = ggplot2::element_blank(),
                                          #axis.title.x = ggplot2::element_blank(),
                                          #axis.text.y = ggplot2::element_blank(),
                                          #axis.title.y = ggplot2::element_blank(),
                                          legend.position = "none",
                                          strip.text.x = ggplot2::element_text(margin = ggplot2::margin(2,0,2,0, unit = "pt")),
                                          plot.margin = grid::unit(c(1,1,1,1), "pt"),
                                          panel.spacing = grid::unit(2, "pt")) +
                           ggplot2::scale_color_viridis_c(direction = -1) +
                           ggplot2::facet_wrap(ggplot2::vars("density"))),
                    varfeatplots)

  return(patchwork::wrap_plots(varfeatplots))
}

