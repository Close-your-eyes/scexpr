#' Feature plot based on a data frame
#'
#' @param data data frame from scexpr::get_data
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
#' @param plot_all_across_split do plot all cells across feature_split in
#' col_split?
#' @param order_discr_explicit
#'
#' @return
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
feature_plot_data <- function(data,
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

                              freq_plot = "..auto..",
                              freq_pos = c(-Inf, Inf),
                              freq_size = 4,
                              name_anno = "..auto..",
                              name_anno_pos = c("..auto..", "title", "annotation"),
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
                              ),
                              contour_expr_freq = F,
                              contour_label_nudge = list(),
                              contour_label_args = list(
                                label.colour = NA,
                                fill = "white",
                                size = 4,
                                color = "black",
                                label.padding = ggplot2::unit(rep(1/20,4), "lines")),
                              contour_same_across_split = T,
                              contour_ggnewscale = F,
                              contour_fun = ggplot2::geom_density2d,
                              contour_path_label = NULL,
                              order_discr_explicit = NULL,
                              plot_all_across_split = F) {

  if (!is.data.frame(data)) {
    stop("data must be data frame from scexpr::get_data")
  }

  if (freq_plot == "..auto..") {
    freq_plot <- nlevels(data[["SO.split"]]) > 1
  } else if (!is.logical(freq_plot)) {
    stop("freq_plot has to be ..auto.., TRUE or FALSE.")
  }
  if (freq_plot && name_anno == "..auto..") {
    name_anno <- "{feature} in {feature_cut_ex}"
  } else if (!freq_plot && name_anno == "..auto..") {
    name_anno <- "{feature} ({freq}) in {feature_cut_ex}"
  }
  if (!is.null(name_anno_pos)) {
    name_anno_pos <- rlang::arg_match(name_anno_pos, multiple = T)
    if (any(name_anno_pos == "..auto..")) {
      if (attr(data, "feature_type") == "gene") {
        if (nlevels(data[["SO.split"]]) > 1) {
          name_anno_pos <- "title"
        } else if (nlevels(data[["SO.split"]]) == 1) {
          name_anno_pos <- "annotation"
        }
      } else if (attr(data, "feature_type") == "meta") {
        # legend title has it
        name_anno_pos <- NULL
      }
    }
  }

  # if (!is.data.frame(data)) {
  #   ## when !is.null(trajectory.slot) and the slot has been found
  #   data_traj <- data[["data_traj"]]
  #   data <- data[["data"]]
  #   plot_traj <- T
  # } else {
  #   plot_traj <- F
  # }

  # get color palette
  col.pal <- scexpr:::get_col_pal(data = data,
                                  col_pal_c_args = col_pal_c_args,
                                  col_pal_d_args = col_pal_d_args)

  data <- scexpr:::check.aliases(feature = attr(data, "feature"), feature_alias, data)

  # excluded cells
  shapeby <- tryCatch(rlang::sym(attr(data, "shape_feature")), error = function(e) NULL)
  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(attr(data, "dim1")), y = !!rlang::sym(attr(data, "dim2")))) +
    ggplot2::geom_point(data = ~dplyr::filter(., cells == 0), ggplot2::aes(shape = !!shapeby), size = pt_size, color = col_ex_cells)

  if (attr(data, "feature_type") == "gene") {
    freqs <- scexpr:::get.freqs2(data = data)
    plot <- do.call(scexpr:::feature_plot_gene, args = list(plot = plot,
                                                            freqs = freqs,
                                                            pt_size = pt_size,
                                                            pt_size_fct = pt_size_fct,
                                                            col_expr = col_expr,
                                                            col_non_expr = col_non_expr,
                                                            col_binary = col_binary,
                                                            freq_plot = freq_plot,
                                                            freq_pos = freq_pos,
                                                            freq_size = freq_size,
                                                            col_split = col_split,
                                                            plot_all_across_split = plot_all_across_split))
  } else {
    freqs <- NULL
    plot <- do.call(scexpr:::feature_plot_meta, args = list(plot = plot,
                                                            order_discr_explicit = order_discr_explicit,
                                                            pt_size = pt_size,
                                                            col_split = col_split,
                                                            plot_all_across_split = plot_all_across_split))
  }

  plot <- plot + theme
  plot <- plot + Gmisc::fastDoCall(ggplot2::theme, args = theme_args)

  plot <-
    plot +
    ggplot2::guides(shape = if (shape_legend_hide) {
      "none"
    } else {
      Gmisc::fastDoCall(ggplot2::guide_legend, args = shape_legend_args)
    })

  plot <- scexpr:::add_facet(plot = plot,
                             facet_grid_row_var = facet_grid_row_var,
                             facet_scales = facet_scales,
                             nrow_inner = nrow_inner,
                             ncol_inner = ncol_inner)

  if (!col_binary) {
    plot <- scexpr:::add_color_scale(plot = plot,
                                     col.pal = col.pal,
                                     col_legend_args = col_legend_args,
                                     col_steps = col_steps,
                                     legendbreaks = legendbreaks,
                                     legendlabels = legendlabels,
                                     col_steps_nice = col_steps_nice,
                                     col_na = col_na)
  }

  if (!is.null(name_anno_pos)) {
    title <- scexpr:::get_title(feature_ex = attr(data, "feature_ex", exact = T),
                                feature_cut = attr(data, "feature_cut", exact = T),
                                feature_cut_expr = attr(data, "feature_cut_expr"),
                                freq = freqs[["freq.expr"]][["freq2"]],
                                feature_italic = attr(data, "feature_type") == "gene",
                                feature = attr(data, "feature"),
                                name_anno = name_anno,
                                #markdown = "element_markdown" %in% class(plot[["theme"]][["plot.title"]])
                                markdown = T)
    if ("title" %in% name_anno_pos) {
      plot <- plot + ggplot2::labs(title = title)
      #plot <- plot + ggplot2::labs(title = substitute(paste(x, sep = ""), list(x = title)))
    }
    if ("annotation" %in% name_anno_pos) {
      library(ggtext)
      plot <- plot + Gmisc::fastDoCall(ggplot2::annotate, args = c(list(geom = "richtext",
                                                                        label = title),
                                                                   name_anno_args))
    }
  }

  plot <- scexpr:::add_axes_expansion(plot = plot,
                                      axes_lim_set = axes_lim_set,
                                      axes_lim_expand = axes_lim_expand)

  label_label <- NULL
  if ("label_feature" %in% names(attributes(plot[["data"]]))) {
    plot <- add_labels(plot = plot,
                       label_filter_cells = label_filter_cells,
                       label_center_fun = label_center_fun,
                       label_nudge = label_nudge,
                       label_repel = label_repel,
                       label_multi_try = label_multi_try,
                       label_multi_max = label_multi_max,
                       label_args = label_args,
                       finalize_plotting = F)
  }

  label_contour <- NULL
  if ("contour_feature" %in% names(attributes(plot[["data"]]))) {
    plot <- add_contour(plot = plot,
                        label_center_fun = label_center_fun,
                        contour_rm_outlier = contour_rm_outlier,
                        contour_rm_lowfreq_subcluster = contour_rm_lowfreq_subcluster,
                        contour_multi_try = contour_multi_try,
                        contour_multi_max = contour_multi_max,
                        contour_filter_cells = contour_filter_cells,
                        contour_col_pal_args = contour_col_pal_args,
                        contour_args = contour_args,
                        contour_label_nudge = contour_label_nudge,
                        contour_label_args = contour_label_args,
                        contour_same_across_split = contour_same_across_split,
                        contour_expr_freq = contour_expr_freq,
                        contour_ggnewscale = contour_ggnewscale,
                        contour_fun = contour_fun,
                        contour_path_label = contour_path_label,
                        finalize_plotting_expr_freq_labels = F)
  }


  if (!is.null(label_label) || !is.null(label_contour)) {
    # label_label and label_contour are assigned in add_labels and add_contour respectively
    # plotting them is combined here to allow repelling
    # contour label have no separate repel argument
    #
    plot <- co_add_feature_and_contour_labels(plot = plot,
                                              label_label = label_label,
                                              label_contour = label_contour,
                                              label_nudge = label_nudge,
                                              label_repel = label_repel,
                                              label_args = label_args,
                                              contour_label_nudge = contour_label_nudge,
                                              contour_label_args = contour_label_args)

  }


  # repel labels from label_feature and contour_expr_freq

  # if (plot_traj) {
  #   plot <- plot + ggplot2::geom_segment(data = data_traj, ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
  #                                        linewidth = trajectory.size,
  #                                        color = trajectory.color,
  #                                        linetype = trajectory.linetype,
  #                                        na.rm = TRUE)
  # }

  # plot cutoff feature inset plot (https://stackoverflow.com/questions/5219671/it-is-possible-to-create-inset-graphs)
  # if (!is.null(feature_cut) && plot.feature_cut.inset) {
  #   if (length(SO) > 1 || !is.null(split_by)) {
  #     message("The inset plot does not appear on every facet.")
  #   }
  #   inset.data <- do.call(rbind, lapply(names(SO), function(y) {
  #     data.frame(cbind(as.matrix(t(Seurat::GetAssayData(SO[[y]], slot = "data", assay = assay)[feature_cut,,drop = F]))))
  #   }))
  #   if (cutoff_expr == 0) {
  #     cutoff_expr.plot <- brathering::floor2(min(inset.data[,feature_cut][which(inset.data[,feature_cut] > 0)]), 0.1)
  #   } else {
  #     cutoff_expr.plot <- cutoff_expr
  #   }
  #   inset <- ggplot2::ggplot(inset.data, ggplot2::aes(!!rlang::sym(feature_cut))) +
  #     ggplot2::geom_density(adjust = 1) +
  #     ggplot2::geom_vline(xintercept = cutoff_expr.plot, color = "red") +
  #     ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), strip.text.x = ggplot2::element_text(size = 9, face = "italic", family = font.family), plot.margin = ggplot2::unit(c(0,0,0,0), "cm"), plot.background = ggplot2::element_rect(fill = "transparent")) +
  #     ggplot2::facet_wrap(ggplot2::vars(!!feature_cut))
  #   plot <- cowplot::ggdraw() + cowplot::draw_plot(plot) + cowplot::draw_plot(inset, x = inset.position[1], y = inset.position[2], width = inset.size[1], height = inset.size[2])
  # }
  return (plot)
}
