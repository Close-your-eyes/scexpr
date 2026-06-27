#' Plot one or more features on a single-cell dimension reduction map
#'
#' Creates feature plots from one Seurat object or a list of Seurat objects.
#' This is a high-level wrapper around `scexpr::get_data()` and
#' `scexpr::feature_plot_data()`. For multiple features, plots can optionally
#' be combined with `patchwork`.
#'
#' @param SO A Seurat object or a list of Seurat objects.
#' @param features Character vector of features to plot. Features may be gene
#'   names or metadata column names. If missing, the function attempts to use
#'   `Idents`, `seurat_clusters`, `orig.ident`, or the first shared metadata
#'   column.
#' @param reduction Name of the dimensionality reduction to use.
#' @param dims Numeric vector of length 2 specifying the dimensions to plot.
#' @param assay Assay used to fetch expression values.
#' @param cells Optional vector of cell names to highlight. Non-selected cells
#'   are drawn as excluded cells.
#' @param get_data_args Additional arguments passed to `scexpr::get_data()`.
#' @param combine Logical; combine plots for multiple features with
#'   `patchwork::wrap_plots()`.
#' @param combine_guides Passed to `patchwork::wrap_plots(guides = )`.
#'   Usually `NULL` or `"collect"`.
#' @param ncol_combine Number of columns in the combined patchwork layout.
#' @param nrow_combine Number of rows in the combined patchwork layout.
#' @param strip_select Numeric vector selecting which plot indices keep facet
#'   strip labels. Other plots have facet strips hidden.
#' @param ... Convenience arguments forwarded to `scexpr::get_data()`, including
#'   `feature_ex`, `feature_cut`, `feature_cut_expr`, `label_feature`,
#'   `contour_feature`, `split_feature`, and `shape_feature`.
#'
#' @section Feature selection and filtering:
#' @param feature_ex Optional feature used to exclude cells.
#' @param feature_cut Optional feature used for expression/value cutoff-based
#'   filtering.
#' @param feature_cut_expr Cutoff value for `feature_cut`.
#' @param label_feature Optional metadata feature used for label placement.
#' @param contour_feature Optional metadata feature used for contour groups.
#' @param split_feature Optional metadata feature used to split plots into
#'   facets.
#' @param shape_feature Optional metadata feature mapped to point shape.
#'
#' @section Point appearance:
#' @param pt_size Base point size.
#' @param pt_size_fct Multiplicative size factor for expressing cells.
#' @param col_expr Colour used for expressing cells.
#' @param col_non_expr Colour used for non-expressing cells. `"..auto.."` chooses
#'   a suitable colour based on the plot background.
#' @param col_ex_cells Colour used for excluded cells.
#' @param col_split Colour used for background cells across split facets when
#'   `plot_all_across_split = TRUE`.
#' @param col_na Colour used for missing values.
#' @param col_binary Logical; show gene expression as binary expressing vs.
#'   non-expressing values.
#'
#' @section Colour scales:
#' @param col_pal_c_args Arguments for continuous colour palettes.
#' @param col_pal_d_args Arguments for discrete colour palettes. If
#'   `name = "..auto.."`, colours may be read from `SO@misc$metacolors`.
#' @param col_steps Colour step specification: `NULL`, `"..auto.."`, a number of
#'   bins, or explicit break values.
#' @param col_steps_nice Logical; use pretty step breaks.
#' @param col_trans_log Logical; apply logarithmic colour transformation.
#' @param col_legend_c_args Arguments passed to the continuous colour guide.
#' @param col_legend_d_args Arguments passed to the discrete colour guide.
#' @param legendbreaks Legend break specification.
#' @param legendlabels Optional labels for legend breaks.
#'
#' @section Legends:
#' @param shape_legend_args Arguments passed to `ggplot2::guide_legend()` for
#'   shape legends.
#' @param shape_legend_hide Logical; hide the shape legend.
#'
#' @section Titles and annotations:
#' @param feature_alias Named character vector mapping original feature names to
#'   display aliases.
#' @param freq_plot Logical or `"..auto.."`; draw expression frequency labels.
#' @param freq_pos Position of frequency labels.
#' @param freq_size Font size of frequency labels.
#' @param freq_col Colour of frequency labels.
#' @param name Plot title template. Supports glue-style variables
#'   `{feature}`, `{freq}`, and `{feature_cut_ex}`.
#' @param anno Plot annotation template. Supports the same variables as `name`.
#' @param name_anno_pos Where to draw name/annotation text: `"title"`,
#'   `"annotation"`, both, `NULL`, or `"..auto.."`.
#' @param name_anno_args Arguments passed to annotation layers.
#' @param title Optional superordinate title. Added as a patchwork annotation
#'   for combined plots, otherwise as a ggplot title.
#'
#' @section Themes and facets:
#' @param theme Base ggplot theme.
#' @param theme_args Additional arguments passed to `ggplot2::theme()`.
#' @param facet_scales Facet scale behaviour.
#' @param facet_grid_row_var Optional variable used as row variable in
#'   `facet_grid()`.
#' @param nrow_inner Number of rows for inner facets.
#' @param ncol_inner Number of columns for inner facets.
#'
#' @section Axes:
#' @param axes_lim_set Named list of fixed axis limits.
#' @param axes_lim_expand Named list of axis expansion values.
#' @param axes_arrows Logical; draw coordinate axes as arrows.
#'
#' @section Labels:
#' @param label_filter_cells Restrict label position calculation to included
#'   cells.
#' @param label_center_fun Method for label centres: `"median"` or `"mean"`.
#' @param label_nudge Named list of x/y nudges for labels.
#' @param label_repel Logical; repel labels to avoid overlap.
#' @param label_multi_try Try to label split clusters separately.
#' @param label_multi_max Maximum number of labels per split cluster attempt.
#' @param label_args Arguments passed to label geoms.
#'
#' @section Contours:
#' @param contour_filter_cells Restrict contour calculation to included cells.
#' @param contour_rm_outlier Remove outlier cells before drawing contours.
#' @param contour_rm_lowfreq_subcluster Remove low-frequency subclusters before
#'   drawing contours.
#' @param contour_multi_try Try to draw contours for split clusters separately.
#' @param contour_multi_max Maximum number of contour groups per split attempt.
#' @param contour_col_pal_args Colour palette arguments for contours.
#' @param contour_args Arguments passed to the contour geom.
#' @param contour_label_nudge Nudges for contour labels.
#' @param contour_label_args Arguments for contour labels.
#' @param contour_same_across_split Use identical contours across split facets.
#' @param contour_expr_freq Logical; label contour groups with expression
#'   frequency.
#' @param contour_fun Contour drawing function.
#' @param contour_path_label Label used when `contour_fun` is a text-path geom.
#'
#' @section Ordering and split plotting:
#' @param order_discr_explicit Explicit drawing order for discrete metadata
#'   values. Use leading `^` to draw a group first and trailing `$` to draw it
#'   last, e.g. `c("^healthy", "disease$")`.
#' @param plot_all_across_split Logical; draw all cells in each split facet
#'   using `col_split`.
#'
#' @section Spatial plotting:
#' @param img_plot Logical; draw the spatial image underneath the points.
#'   Supported only for a single Seurat object.
#' @param img_subset_fct Subsampling factor passed to `get_img()`.
#' @param img_col_conv_fun Function used to convert image colours.
#' @param cell_hull_plot Logical; draw cell hull polygons. Supported only for a
#'   single Seurat object.
#' @param cell_hull_args Arguments passed to the hull polygon layer.
#'
#' @return
#' A `ggplot2` object, a list of `ggplot2` objects, or a `patchwork` object,
#' depending on `features` and `combine`.
#'
#' @details
#' `feature_plot2()` first retrieves plotting data with `scexpr::get_data()`,
#' then calls `scexpr::feature_plot_data()` for each requested feature. When
#' more than one feature is supplied and `combine = TRUE`, the individual plots
#' are combined using `patchwork::wrap_plots()`.
#'
#' If `features` is omitted, the function tries to choose a sensible metadata
#' feature automatically, preferring active Seurat identities, then
#' `seurat_clusters`, then `orig.ident`, and finally the first shared metadata
#' column across all provided Seurat objects.
#'
#' @seealso
#' `scexpr::get_data()`, `scexpr::feature_plot_data()`,
#' `patchwork::wrap_plots()`
#'
#' @examples
#' feature_plot2(
#'   SO,
#'   features = c("MS4A1", "CD3D"),
#'   reduction = "umap",
#'   assay = "RNA",
#'   combine = TRUE,
#'   ncol_combine = 2
#' )
#'
#' @export
feature_plot2 <- function(
    SO,
    features,
    reduction = "umap",
    dims = c(1,2),
    assay = "RNA",
    cells = NULL,
    get_data_args = list(
      qmin = 0,
      qmax = 1,
      order_discr = T,
      shuffle = F
    ),
    combine = T,
    combine_guides = NULL,
    ncol_combine = NULL,
    nrow_combine = NULL,
    strip_select = NULL,
    ...,
    pt_size = 0.3,
    pt_size_fct = 1,
    col_expr = "tomato2",
    col_non_expr = "..auto..",
    col_ex_cells = "grey95",
    col_split = "grey30",
    col_na = "grey50",
    col_binary = F,
    col_pal_c_args = list(
      name = "spectral",
      direction = -1
    ),
    col_pal_d_args = list(
      name = "..auto..",
      missing_fct_to_na = T
    ),
    col_steps = "..auto..",
    col_steps_nice = T,
    col_trans_log = F,
    col_legend_c_args = list(
      barwidth = 0.5,
      barheight = 7,
      title = "..auto..",
      order = 1
    ),
    col_legend_d_args = list(
      nrow = 12,
      override.aes = list(size = 4),
      title = "..auto..",
      order = 1
    ),
    legendbreaks = "minmidmax",
    legendlabels = "..auto..",
    shape_legend_args = list(
      override.aes = list(size = 4),
      order = 2
    ),
    shape_legend_hide = F,
    feature_alias = NULL,
    feature_ex = NULL,
    feature_cut = NULL,
    feature_cut_expr = 0,
    freq_plot = "..auto..",
    freq_pos = c(Inf, Inf),
    freq_size = 4,
    freq_col = "..auto..",
    name = "..auto..",
    # "{feature} ({freq}) in {feature_cut_ex}"
    anno = "..auto..",
    # "{feature} ({freq}) in {feature_cut_ex}"
    name_anno_pos = c("..auto..", "title", "annotation"),
    # or NULL
    name_anno_args = list(
      x = "..auto..",
      y = "..auto..",
      hjust = "..auto..",
      vjust = "..auto..",
      color = NA,
      fill = NA,
      label.color = NA,
      size = 4,
      text.color = "..auto.."
    ),
    theme = colrr::theme_material(
      white = T,
      legend_tight = T,
      text_fun = ggplot2::element_text
    ),
    theme_args = list(#legend.direction = "vertical",
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key.size = grid::unit(0.3, "cm"),
      legend.key = ggplot2::element_blank(),
      title = ggtext::element_markdown(),
      plot.title = ggtext::element_markdown(margin = ggplot2::margin(b = 0, t = 2, unit = "pt")),
      strip.text.x = ggplot2::element_text(margin = ggplot2::margin(2,0,2,0, unit = "pt")),
      #strip.background = ggplot2::element_rect(color = "white", fill = "grey95"),
      plot.margin = grid::unit(c(1,1,1,1), "pt"),
      panel.spacing = grid::unit(2, "pt")),
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
      #contour_var = "ndensity",
      breaks = 0.3,
      linewidth = 0.5,
      linetype = "solid"
    ),
    contour_label_nudge = list(),
    contour_label_args = list(
      label.colour = NA,
      fill = "white",
      size = 4,
      color = "black",
      label.padding = ggplot2::unit(rep(1/20,4), "lines")),
    contour_same_across_split = T,
    contour_expr_freq = F,
    contour_fun = ggplot2::geom_density2d,
    contour_path_label = NULL,
    order_discr_explicit = NULL,
    plot_all_across_split = F,
    label_feature = NULL,
    contour_feature = NULL,
    split_feature = NULL,
    shape_feature = NULL,
    title = NULL,
    axes_arrows = F,
    img_plot = F,
    img_subset_fct = 4,
    img_col_conv_fun = colrr::hex_to_grey,
    cell_hull_plot = F,
    cell_hull_args = list(color = "grey30",
                          linewidth = 0.1,
                          alpha = 1)
) {

  ## ggnewscale breaks the legend of dot colors; setting to F will avoid that but also does not allow to have a legend for contour lines

  ## label position calculation is not facetted!!

  contour_ggnewscale <- F # T not tested yet

  if (!requireNamespace("colrr", quietly = T)) {
    pak::pak("Close-your-eyes/colrr")
  }

  SO <- scexpr:::check.SO(SO = SO)

  if (missing(features)) {
    if (all(purrr::map_lgl(SO, ~!is.null(Seurat::Idents(.x))))) {
      SO <- purrr::map(SO, ~Seurat::AddMetaData(.x, Seurat::Idents(.x), col.name = "Idents"))
      features <- "Idents"
    } else if (all(purrr::map_lgl(SO, ~"seurat_clusters" %in% names(.x@meta.data)))) {
      features <- "seurat_clusters"
    } else if (all(purrr::map_lgl(SO, ~"orig.ident" %in% names(.x@meta.data)))) {
      features <- "orig.ident"
    } else {
      sharedfeat <- Reduce(intersect, purrr::map(SO, ~names(.x@meta.data)))
      if (!length(sharedfeat)) {
        stop("no common feature found.")
      }
      features <- sharedfeat[1]
    }
  }

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
                                                           label_feature = label_feature,
                                                           contour_feature = contour_feature,
                                                           split_feature = split_feature,
                                                           shape_feature = shape_feature,
                                                           feature_cut = feature_cut,
                                                           feature_cut_expr = feature_cut_expr,
                                                           feature_ex = feature_ex),
                                                      get_data_args))

  if (img_plot) {
    if (length(SO)>1) {
      stop("img_plot only with one SO.")
    }
    img_df <- get_img(
      obj = SO[[1]],
      subset_factor = img_subset_fct,
      col_conv_fun = img_col_conv_fun)
  } else {
    img_df <- NULL
  }

  if (cell_hull_plot) {
    if (length(SO)>1) {
      stop("cell_hull_plot only with one SO.")
    }
    hull_df <- get_cell_hulls(obj = SO[[1]])
  } else {
    hull_df <- NULL
  }


  if (is.null(col_pal_d_args[["name"]])) {
    col_pal_d_args[["name"]] <- "custom"
  }
  col_pal_d_args_lst <- purrr::map(data, function(x) {
    y <- col_pal_d_args
    if (y[["name"]][1] == "..auto..") {
      y[["name"]] <- "custom"
    }
    return(y)
  })

#browser()
  if (col_pal_d_args[["name"]][1] == "..auto..") {
    if (all(purrr::map_lgl(SO, ~"metacolors" %in% names(.x@misc))) && all(purrr::map_lgl(SO, ~is.list(.x@misc[["metacolors"]])))) {
      for (i in names(col_pal_d_args_lst)) {
        if (all(purrr::map_lgl(SO, ~i %in% names(.x@misc[["metacolors"]])))) {
          y <- purrr::list_c(purrr::map(SO, ~.x@misc[["metacolors"]][[i]]))
          col_pal_d_args_lst[[i]][["name"]] <- y[which(!duplicated(names(y)))]
        }
      }
    }
  }


  plots <- purrr::map2(data,
                       col_pal_d_args_lst,
                       ~Gmisc::fastDoCall(what = feature_plot_data,
                                          args = list(data = .x,
                                                      col_binary = col_binary,
                                                      col_expr = col_expr,
                                                      col_ex_cells = col_ex_cells,
                                                      col_na = col_na,
                                                      col_non_expr = col_non_expr,
                                                      col_pal_c_args = col_pal_c_args,
                                                      col_pal_d_args = .y,
                                                      col_legend_c_args = col_legend_c_args,
                                                      col_legend_d_args = col_legend_d_args,
                                                      col_steps = col_steps,
                                                      col_trans_log = col_trans_log,
                                                      facet_grid_row_var = facet_grid_row_var,
                                                      feature_alias = feature_alias,
                                                      freq_size = freq_size,
                                                      freq_pos = freq_pos,
                                                      freq_plot = freq_plot,
                                                      freq_col = freq_col,
                                                      name = name,
                                                      anno = anno,
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
                                                      plot_all_across_split = plot_all_across_split,
                                                      axes_arrows = axes_arrows,
                                                      img_df = img_df,
                                                      hull_df = hull_df,
                                                      cell_hull_args = cell_hull_args)
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
    plots <- patchwork::wrap_plots(plots, ncol = ncol_combine, nrow = nrow_combine, guides = combine_guides)
  }
  if ((length(plots) == 1 && !combine) || length(features) == 1) {
    plots <- plots[[1]]
  }

  if (!is.null(title)) {
    if ("patchwork" %in% class(plots)) {
      plots <- plots +
        patchwork::plot_annotation(
          title = title
        )
    } else {
      plots <- plots + ggplot2::labs(title = title)
    }
  }


  return(plots)
}

