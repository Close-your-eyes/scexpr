#' Plot features of single cell transcriptomes on a dimension reduction map
#'
#' This function wraps around scexpr::get_data, scexpr::feature_plot_data
#' and optional combination of plots with patchwork.
#'
#' @param SO one Seurat object or a list of multiple ones
#' @param features vector of features to fetch (genes or column names in
#' meta data)
#' @param reduction reduction to fetch from reductions slot of SO
#' @param dims dimensions of the selected reduction to extract from SO
#' @param assay which assay to get expression data from
#' @param cells vector of cell names to select for regular plotting (=1),
#' deselected ones are plotted with color col_ex_cells (=0) in feature_plot
#' @param get_data_args further arguments to scexpr::get_data
#' @param combine do combine plots for multiple features with patchwork
#' @param ncol_combine number of columns when combine
#' @param nrow_combine number of rows when combine
#' @param strip_select numeric: which facet strips to plot, NULL: plot all
#' @param ... several feature arguments to scexpr::get_data, for convenience,
#' like feature_ex or split_feature
#' @param pt_size point size
#' @param pt_size_fct factor of size increase for expressers
#' @param col_expr color expressers
#' @param col_non_expr color non expressers
#' @param col_ex_cells color excluded cells
#' @param col_split color for all cells across facets, requires
#' plot_all_across_split = T and a split_feature in get_data or ...
#' @param col_na color for NA values
#' @param col_binary dichotomous coloring of expresser and non expresser
#' @param col_pal_c_args for continuous color: arguments to colrr::col_pal,
#' name can be (i) name of a color palette invoking colrr::col_pal, (ii)
#' a single color or (iii) a vector of colors
#' @param col_pal_d_args for discrete color
#' @param col_steps NULL to have normal colorbar, auto for default colorsteps,
#' a single number or a vector of explicit steps; may not work with any number
#' when col_steps_nice is TRUE; fcexpr:::colorscale_heuristic is used
#' @param col_steps_nice algorithmic determination of pretty steps,
#' see ggplot2::scale_color_stepsn
#' @param col_legend_args arguments to ggplot2::guide_colorsteps,
#' ggplot2::guide_colorbar or ggplot2::guide_legend for binned or continuous
#' continuous or discrete legend
#' @param legendbreaks a single number, a vector of explicit breaks, or "auto"
#' for ggplot default or "minmidmax" for three breaks at minimum, middle and
#' maximum of value range
#' @param legendlabels labels for breaks, e.g. c("min", "mid", "max")
#' @param shape_legend_args arguments to ggplot2::guide_legend for shapes,
#' requires a shape_feature in get_data or ...
#' @param shape_legend_hide do not show shape legend?
#' @param feature_alias named vector of aliases for features, old = new; e.g.
#' c("FAIM3" = "FCMR", "MS4A1" = "CD20")
#' @param freq_plot do plot frequency of expressers? TRUE, FALSE or ..auto..
#' for decision based on number of SOs
#' @param freq_pos where to plot, xy coordinates
#' @param freq_size font size of frequency annotation
#' @param name_anno name of the plot, any string possible, understands glue
#' syntax with three possible variables like '{feature} ({freq}) in {feature_cut_ex}',
#' where freq is global expression frequency and feature_cut_ex is a collapsed
#' string of feature_cut, feature_cut_expr, and feature_ex; '..auto..' for
#' algorithmic decision based on freq_plot; freq is removed for meta features;
#' feature_cut_ex is adjusted if either element is missing
#' @param name_anno_pos where to plot name_anno, NULL to omit plotting, '..auto..'
#' for automatic decision between title, annotation or NULL; title and annotation
#' simultaneously also possible
#' @param name_anno_args arguments to ggplot2::annotate when name_anno_pos
#' contains 'annotation'
#' @param theme ggplot theme to use
#' @param theme_args arguments to ggplot2::theme
#' @param facet_scales passed as scales argument of facet_wrap, relevant when
#' number of SO > 1 and/or feature_split is given
#' @param facet_grid_row_var forces using facet_grid instead of facet_wrap,
#' facet_grid_row_var is a string that appears as row label on facets
#' @param nrow_inner row number for facets
#' @param ncol_inner col number for facets
#' @param axes_lim_set set axes limits like
#' list(UMAP_1 = c(-10,10), UMAP_2 = c(-12,13))
#' @param axes_lim_expand extend axes limits, syntax like axes_lim_set
#' @param label_filter_cells when feature_label is provided; only consider
#' cells which are not excluded for label position calculation
#' @param label_center_fun how to calculate label position
#' @param label_nudge named list of xy-values how to shift labels, names must
#' be the values of feature_label
#' @param label_repel do repel labels to avoid overlap with ggrepel?
#' @param label_multi_try try to label split clusters with one label each?
#' @param label_args arguments to ggtext::geom_richtext
#' @param contour_filter_cells when feature_contour is provided; only consider
#' cells which are not excluded for contour calculation
#' @param contour_rm_outlier per value (group, cluster) in feature_contour: remove
#' outlier points to avoid ugly contours? uses brathering::outlier
#' @param contour_rm_lowfreq_subcluster groups in feature_contour are checked for
#' split clusters with diptest and mclust, do remove low frequency subclusters
#' to avoid ugly contours?
#' @param contour_col_pal_args discrete colors of contours, see col_pal_c_args
#' @param contour_args arguments to geom_density_2d
#' @param contour_expr_freq do plot expression frequency per group/cluster of
#' feature_contour
#' @param contour_label_nudge one xy value pair to nudge contour_expr_freq
#' labels
#' @param contour_label_args arguments to ggtext::geom_richtext for
#' contour_expr_freq
#' @param contour_same_across_split have the same contours across feature_split
#' facets
#' @param contour_ggnewscale use ggnewscale for colored contours
#' @param contour_fun function for contour drawing, either one of:
#' ggplot2::geom_density2d, geomtextpath::geom_textpath,
#' geomtextpath::geom_labelpath; geomtexpath functions require a
#' contour_path_label
#' @param contour_path_label if geomtextpath used as contour_fun,
#' what label to draw, could be an attached meta feature that is unique to each
#' group or expression frequency via: 'pct'
#' @param order_discr_explicit explicit plotting order for discrete meta features,
#' can be an incomplete vector of feature values, use leading '^' to plot a
#' group first and trailing '$' to plot last; e.g. c("^healthy", "disease$")
#' @param plot_all_across_split do plot all cells across feature_split in
#' col_split?
#' @param feature_ex
#' @param feature_cut
#' @param feature_cut_expr
#'
#' @return
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
#'
feature_plot2 <- function(SO,
                          features,
                          reduction = "umap",
                          dims = c(1,2),
                          assay = "RNA",
                          cells = NULL,
                          get_data_args = list(qmin = 0,
                                               qmax = 1,
                                               label_feature = NULL),
                          combine = T,
                          ncol_combine = NULL,
                          nrow_combine = NULL,
                          strip_select = NULL,
                          ...,

                          pt_size = 0.3,
                          pt_size_fct = 1,
                          col_expr = "tomato2",
                          col_non_expr = "grey85",
                          col_ex_cells = "grey95",
                          col_split = "grey30",
                          col_na = "grey50",
                          col_binary = F,
                          col_pal_c_args = list(name = "spectral", direction = -1),
                          col_pal_d_args = list(name = "custom"),
                          col_steps = "auto",
                          col_steps_nice = T,
                          col_legend_args = list(barwidth = 1,
                                                 barheight = 8,
                                                 override.aes = list(size = 4),
                                                 title.theme = ggtext::element_markdown(),
                                                 title = "..auto..",
                                                 order = 1),
                          legendbreaks = "minmidmax",
                          legendlabels = "auto",
                          shape_legend_args = list(override.aes = list(size = 4),
                                                   order = 2),
                          shape_legend_hide = F,
                          feature_alias = NULL,
                          feature_ex = NULL,
                          feature_cut = NULL,
                          feature_cut_expr = 0,

                          freq_plot = "..auto..",
                          freq_pos = c(-Inf, Inf),
                          freq_size = 4,
                          name_anno = "..auto..", # "{feature} ({freq}) in {feature_cut_ex}"
                          name_anno_pos = c("..auto..", "title", "annotation"), # or NULL
                          name_anno_args = list(x = Inf,
                                                y = Inf,
                                                hjust = 1.1,
                                                vjust = 1.25,
                                                color = NA,
                                                size = 4,
                                                text.color = "black"),

                          theme = ggplot2::theme_bw(),
                          theme_args = list(axis.ticks = ggplot2::element_blank(),
                                            axis.text.x = ggplot2::element_blank(),
                                            axis.text.y = ggplot2::element_blank(),
                                            axis.title.x = ggplot2::element_blank(),
                                            axis.title.y = ggplot2::element_blank(),
                                            panel.grid = ggplot2::element_blank(),
                                            legend.background = ggplot2::element_blank(),
                                            legend.key.size = ggplot2::unit(0.3, "cm"),
                                            legend.key = ggplot2::element_blank()),

                          facet_scales = c("fixed", "free", "free_x", "free_y"),
                          facet_grid_row_var = NULL,

                          nrow_inner = NULL,
                          ncol_inner = NULL,
                          axes_lim_set = list(),
                          axes_lim_expand = list(),
                          label_filter_cells = T,
                          label_center_fun = c("median", "mean"),
                          label_nudge = list(),
                          label_repel = F,
                          label_multi_try = F,
                          label_multi_max = 3,
                          label_args = list(
                            label.colour = NA,
                            fill = "white",
                            size = 4,
                            color = "black",
                            label.padding = ggplot2::unit(rep(1/20,4), "lines")),

                          contour_filter_cells = T,
                          contour_rm_outlier = F,
                          contour_rm_lowfreq_subcluster = F,
                          contour_multi_try = F,
                          contour_multi_max = 3,
                          contour_col_pal_args = list(name = "custom"),
                          contour_args = list(
                            contour_var = "ndensity",
                            breaks = 0.3,
                            linewidth = 0.5
                          ), # arguments to geom_density_2d
                          contour_label_nudge = list(),
                          contour_label_args = list(
                            label.colour = NA,
                            fill = "white",
                            size = 4,
                            color = "black",
                            label.padding = ggplot2::unit(rep(1/20,4), "lines")),
                          contour_same_across_split = T,
                          contour_expr_freq = F,
                          contour_ggnewscale = F,
                          contour_fun = ggplot2::geom_density2d,
                          contour_path_label = NULL,
                          order_discr_explicit = NULL,
                          plot_all_across_split = F) {

 ## ggnewscale breaks the legend of dot colors; setting to F will avoid that but also does not allow to have a legend for contour lines

  ## add axis arrows, shortened:
  #p_blood <- p_blood + guides(x = ggh4x::guide_axis_truncated(trunc_lower = ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1], trunc_upper = ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1] + abs(min(c(ggplot_build(p)$layout$panel_params[[1]]$x.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1])) - max(c(ggplot_build(p)$layout$panel_params[[1]]$x.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1])))/4 ),
  #                          y = ggh4x::guide_axis_truncated(trunc_lower = ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1], trunc_upper = ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1] + abs(min(c(ggplot_build(p)$layout$panel_params[[1]]$y.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1])) - max(c(ggplot_build(p)$layout$panel_params[[1]]$y.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1])))/4 ))


  ## label position calculation is not facetted!!

  dotlist <- list(...)
  for (i in c("contour_feature", "label_feature", "split_feature", "shape_feature", "feature_ex", "feature_cut", "feature_cut_expr")) {
    if (!is.null(dotlist[[i]])) {get_data_args[[i]] <- dotlist[[i]]}
  }
  label_center_fun <- rlang::arg_match(label_center_fun)


  data <- Gmisc::fastDoCall(what = get_data, args = c(list(SO = SO,
                                                           feature = features,
                                                           assay = assay,
                                                           cells = cells,
                                                           reduction = reduction,
                                                           dims = dims,
                                                           feature_cut = feature_cut,
                                                           feature_cut_expr = feature_cut_expr,
                                                           feature_ex = feature_ex),
                                                      get_data_args))

  plots <- purrr::map(data, ~Gmisc::fastDoCall(what = feature_plot_data,
                                               args = list(data = .x,
                                                           col_binary = col_binary,
                                                           col_expr = col_expr,
                                                           col_ex_cells = col_ex_cells,
                                                           col_na = col_na,
                                                           col_non_expr = col_non_expr,
                                                           col_pal_c_args = col_pal_c_args,
                                                           col_pal_d_args = col_pal_d_args,
                                                           col_legend_args = col_legend_args,
                                                           col_steps = col_steps,
                                                           facet_grid_row_var = facet_grid_row_var,
                                                           feature_alias = feature_alias,
                                                           freq_size = freq_size,
                                                           freq_pos = freq_pos,
                                                           freq_plot = freq_plot,
                                                           name_anno = name_anno,
                                                           name_anno_pos = name_anno_pos,
                                                           name_anno_args = name_anno_args,
                                                           legendbreaks = legendbreaks,
                                                           legendlabels = legendlabels,
                                                           col_steps_nice = col_steps_nice,
                                                           ncol_inner = ncol_inner,
                                                           nrow_inner = nrow_inner,
                                                           pt_size = pt_size,
                                                           pt_size_fct = pt_size_fct,
                                                           shape_legend_args = shape_legend_args,
                                                           shape_legend_hide = shape_legend_hide,
                                                           facet_scales = facet_scales,
                                                           theme = theme,
                                                           theme_args = theme_args,
                                                           axes_lim_set = axes_lim_set,
                                                           axes_lim_expand = axes_lim_expand,
                                                           label_filter_cells = label_filter_cells,
                                                           label_center_fun = label_center_fun,
                                                           label_nudge = label_nudge,
                                                           label_repel = label_repel,
                                                           label_multi_try = label_multi_try,
                                                           label_multi_max = label_multi_max,
                                                           label_args = label_args,
                                                           order_discr_explicit = order_discr_explicit,
                                                           contour_filter_cells = contour_filter_cells,
                                                           contour_rm_outlier = contour_rm_outlier,
                                                           contour_rm_lowfreq_subcluster = contour_rm_lowfreq_subcluster,
                                                           contour_multi_try = contour_multi_try,
                                                           contour_multi_max = contour_multi_max,
                                                           contour_col_pal_args = contour_col_pal_args,
                                                           contour_args = contour_args,
                                                           contour_label_nudge = contour_label_nudge,
                                                           contour_label_args = contour_label_args,
                                                           contour_same_across_split = contour_same_across_split,
                                                           contour_expr_freq = contour_expr_freq,
                                                           contour_ggnewscale = contour_ggnewscale,
                                                           contour_fun = contour_fun,
                                                           contour_path_label = contour_path_label,
                                                           col_split = col_split,
                                                           plot_all_across_split = plot_all_across_split)
  ))



  if (!is.null(strip_select)) {
    for (i in 1:length(plots)) {
      if (!i %in% strip_select) {
        plots[[i]] <- plots[[i]] + ggplot2::theme(strip.text.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank())
      }
    }
  }

  if (combine && length(features) > 1) {
    if (!is.null(ncol_combine) && !is.null(nrow_combine)) {
      message("ncol_combine set NULL as ncol_combine and nrow_combine are provided.")
      ncol_combine <- NULL
    }
    plots <- patchwork::wrap_plots(plots, ncol = ncol_combine, nrow = nrow_combine)
  }
  if (length(plots) == 1 && !combine) {
    plots <- plots[[1]]
  }
  return(plots)
}

