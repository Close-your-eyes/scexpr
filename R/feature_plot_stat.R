#' Plot Feature Expression with Optional Statistics
#'
#' Creates faceted expression plots for one or more features (e.g., genes)
#' from a Seurat object, grouped by a metadata column. Supports jitter/point
#' overlays, violin/boxplots, and optional annotation of the percentage of
#' expressing cells per group.
#'
#' @param SO A Seurat object or a list of Seurat objects.
#' @param features Character vector of feature names (e.g., genes) to plot.
#' @param meta_col Character scalar specifying the metadata column used for grouping.
#' @param get_data_args List of arguments passed to \code{get_data()}, including
#'   \code{assay}, \code{layer}, and \code{cells}.
#' @param geom1 Character. Type of point layer to use: \code{"jitter"} or \code{"point"}.
#' @param geom2 Character. Secondary geometry: \code{"boxplot"}, \code{"violin"}, or \code{"none"}.
#' @param plot_first Character. Which geometry to plot first: \code{"geom1"} or \code{"geom2"}.
#' @param jitterwidth Numeric. Width of jitter for point placement.
#' @param dodgewidth Numeric. Width for position dodging.
#' @param expr_freq_size Numeric. Text size for expression frequency labels; set to 0 to disable.
#' @param expr_freq_angle Numeric. Angle of expression frequency text labels.
#' @param expr_freq_decimals Integer. Number of decimal places for expression frequency.
#' @param expr_freq_pct Logical. Whether to display expression frequency as percentage.
#' @param expr_freq_hjust Numeric. Vertical offset for expression frequency labels.
#' @param plot_non_expr Logical. Whether to include cells with zero expression.
#' @param plot_strip Logical. Whether to display facet strip labels.
#' @param legend_title Character. Title for the legend (used when multiple objects are provided).
#' @param pt_size Numeric. Size of points in the plot.
#' @param col_pal Color palette (vector or function) used for coloring groups.
#' @param col_na Color used for missing values.
#' @param col_pal_args List of additional arguments passed to \code{colrr::make_col_pal()}.
#' @param feature_cut Optional filtering threshold for feature inclusion.
#' @param feature_cut_expr Numeric. Expression threshold used with \code{feature_cut}.
#' @param feature_ex Optional vector of features to exclude.
#' @param theme A ggplot2 theme object applied to the plot.
#' @param geom2_args List of additional arguments passed to the secondary geometry.
#' @param facetting_args List of arguments passed to \code{ggplot2::facet_wrap()}.
#' @param axis_expansion_y_mult Numeric vector of length 2 controlling y-axis expansion.
#' @param geom3 geom only plotted on expressing cells (feature>0). when ..auto..
#' a boxplot is plotted only on those groups where median across all cells is zero.
#' @param expr_freq_y_equal_max equal max y across all features to plot expr freq annotation?
#' @param theme_args args to theme
#' @param geom3_args args to geom3
#' @param add_pwc add default ggpubr::geom_pwc?
#'
#' @details
#' The function extracts expression data using \code{get_data()} and visualizes
#' it across groups defined by \code{meta_col}. Multiple features are displayed
#' as facets. Optionally, the percentage (or fraction) of expressing cells
#' (expression > 0) is computed per group and displayed as text annotations.
#'
#' When multiple Seurat objects are provided, they are combined and colored
#' according to their origin. Custom color palettes are generated using the
#' \pkg{colrr} package.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' feature_plot_stat(
#'   SO = seurat_obj,
#'   features = c("CD3D", "MS4A1"),
#'   meta_col = "cell_type",
#'   geom1 = "jitter",
#'   geom2 = "violin",
#'   expr_freq_pct = TRUE
#' )
#' }
#'
#' @seealso \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{geom_violin}},
#'   \code{\link[ggplot2]{geom_boxplot}}
#'
#' @export
feature_plot_stat <- function(SO,
                              features,
                              meta_col = "orig.ident",
                              get_data_args = list(assay = "RNA",
                                                   layer = "data",
                                                   cells = NULL),
                              geom1 = c("jitter", "point", "sina", "none"),
                              geom2 = c("boxplot", "violin", "none"),
                              geom3 = c("none", "boxplot", "violin", "..auto.."),
                              plot_first = c("geom1", "geom2"),
                              jitterwidth = 0.2,
                              dodgewidth = 0.9,
                              expr_freq_size = 4,
                              expr_freq_angle = 90,
                              expr_freq_decimals = 2,
                              expr_freq_pct = F,
                              expr_freq_hjust = 0.3,
                              expr_freq_y_equal_max = F,
                              plot_non_expr = T,
                              plot_strip = T,
                              legend_title = "object",
                              pt_size = 0.5,
                              col_pal = "..auto..",
                              col_na = "grey50",
                              col_pal_args = list(missing_fct_to_na = T),
                              feature_cut = NULL,
                              feature_cut_expr = 0,
                              feature_ex = NULL,
                              theme = colrr::theme_material(white = T, style = "prismy"),
                              theme_args = list(strip.text.x = ggtext::element_markdown(margin = ggplot2::margin(2,0,2,0, unit = "pt")),
                                                axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)),
                              geom2_args = list(outlier.shape = NA, alpha = 0.7),
                              geom3_args = list(outlier.shape = NA, alpha = 0.5, width = 0.4),
                              facetting_args = list(scales = "free_y",
                                                    axes = "all",
                                                    axis.labels = "margins"),
                              axis_expansion_y_mult = c(0.02,0.125),
                              add_pwc = F) {


  # handle facetting
  # check fun from promotion



  if (!requireNamespace("colrr", quietly = T)) {
    pak::pak("Close-your-eyes/colrr")
  }

  if (missing(SO)) {
    stop("Seurat object list or feature vector is missing.")
  }
  if (is.na(meta_col) || missing(meta_col) || meta_col == "" || is.null(meta_col)) {
    stop("meta_col required.")
  }
  if (length(meta_col) > 1) {
    meta_col <- meta_col[1]
    message("Please provide only one meta_col. Using first element of meta_col.")
  }

  geom1 <- rlang::arg_match(geom1)
  geom2 <- rlang::arg_match(geom2)
  geom3 <- rlang::arg_match(geom3)
  plot_first <- rlang::arg_match(plot_first)

  SO <- scexpr:::check.SO(SO = SO, assay = get_data_args[["assay"]])
  assay <- Seurat::DefaultAssay(SO[[1]])
  features <- scexpr:::check.features(SO = SO, features = features, meta.data = T)
  meta_col <- scexpr:::check.features(SO = SO, features = meta_col, rownames = F)

  get_data_args[["cells"]] <- scexpr:::check.and.get.cells(SO = SO,
                                                           assay = get_data_args[["assay"]],
                                                           cells = get_data_args[["cells"]],
                                                           feature_cut = feature_cut,
                                                           feature_cut_expr = feature_cut_expr,
                                                           feature_ex = feature_ex,
                                                           included_only = T)

  data <- Gmisc::fastDoCall(get_data, args = c(list(SO = SO,
                                                    reduction = NULL,
                                                    meta_col = meta_col,
                                                    feature = features), get_data_args))

  # data <- get_data(SO,
  #                  feature = features,
  #                  reduction = NULL,
  #                  assay = assay,
  #                  layer = layer,
  #                  cells = cells,
  #                  meta_col = meta_col,
  #                  #try_df = T,
  #                  split_feature = split_feature)
  all_gene_feat <- unique(unlist(purrr::map(SO, scexpr:::get_gene_features)))
  data <- dplyr::bind_rows(data, .id = "feature_split") ## different features

  # split_feature requires testing - include in pivoting etc and geom_text
  #data <- tidyr::pivot_longer(data, dplyr::all_of(features), names_to = "feature_split", values_to = "expr")

  ylab <- "feature"
  if (all(features %in% all_gene_feat)) {
    if (get_data_args[["layer"]] == "data") {
      ylab <- "norm UMI"
    } else if (get_data_args[["layer"]] == "counts") {
      ylab <- "UMI"
    }
  }

  if (expr_freq_size>0) {
    stat <- data |>
      dplyr::mutate(max.feat.expr = max(feature), .by = feature_split) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c("feature_split", meta_col, "SO.split", "max.feat.expr")))) |>
      dplyr::summarise(pct.expr = sum(feature > 0)/dplyr::n(), .groups = "drop") |>
      dplyr::mutate(pct.expr.adjust.pct = dplyr::case_when(pct.expr == 0 ~ "0 %",
                                                           pct.expr > 0 & pct.expr < 0.01 ~ "> 1 %",
                                                           pct.expr >= 0.01 ~ paste0(round(pct.expr*100, expr_freq_decimals), " %"))) |>
      dplyr::mutate(pct.expr.adjust = dplyr::case_when(pct.expr == 0 ~ "0",
                                                       pct.expr > 0 & pct.expr < 0.01 ~ "> 0.01",
                                                       pct.expr >= 0.01 ~ as.character(round(pct.expr, expr_freq_decimals)))) |>
      dplyr::filter(feature_split %in% all_gene_feat) |>
      dplyr::mutate(feature_split = paste0("*", feature_split, "*"))

    names(stat)[which(names(stat) == "SO.split")] <- legend_title
  }
  names(data)[which(names(data) == "SO.split")] <- legend_title
  data <- data |>
    dplyr::mutate(feature_split = ifelse(feature_split %in% all_gene_feat,
                                         paste0("*", feature_split, "*"),
                                         feature_split))
  features <- ifelse(features %in% all_gene_feat,
                     paste0("*", features, "*"),
                     features)

  if (expr_freq_y_equal_max) {
    stat$max.feat.expr <- max(stat$max.feat.expr)
  }

  if (!plot_non_expr) {
    data <- dplyr::filter(data, feature > 0)
  }

  geom1pos <- if (geom1 == "jitter" && jitterwidth > 0) ggplot2::position_jitterdodge(jitter.width = jitterwidth, dodge.width = dodgewidth) else ggplot2::position_dodge(width = dodgewidth)

  if (geom3 != "none") {
    if (geom3 == "..auto..") {
      datasummary <- dplyr::summarise(data, median = median(feature), .by = c(!!rlang::sym(meta_col), feature_split))
      median0 <- datasummary |> dplyr::filter(median == 0)
      datageom3 <- data |>
        dplyr::filter(feature > 0) |>
        dplyr::filter(!!rlang::sym(meta_col) %in% median0[[meta_col]] & feature_split %in% median0[["feature_split"]])
    } else {
      datageom3 <- dplyr::filter(data, feature > 0)
    }
  } else {
    datageom3 <- data
  }

  geom1_layer <- switch(geom1,
                  "jitter" = ggplot2::geom_point(size = pt_size, position = geom1pos),
                  "point" = ggplot2::geom_point(size = pt_size, position = geom1pos),
                  "sina" = ggforce::geom_sina(size = pt_size),
                  "none" = ggplot2::geom_blank())
  geom2 <- switch(geom2,
                  "violin" = ggplot2::geom_violin,
                  "boxplot" = ggplot2::geom_boxplot,
                  "none" = ggplot2::geom_blank)
  geom3 <- switch(geom3,
                  "violin" = ggplot2::geom_violin,
                  "boxplot" = ggplot2::geom_boxplot,
                  "..auto.." = ggplot2::geom_boxplot,
                  "none" = ggplot2::geom_blank)

  color_aes <- ifelse(length(SO) > 1, legend_title, meta_col)


  if (col_pal[1] == "..auto..") {
    if ("metacolors" %in% names(SO[[1]]@misc) && is.list(SO[[1]]@misc[["metacolors"]])) {
      if (meta_col %in% names(SO[[1]]@misc[["metacolors"]])) {
        col_pal <- SO[[1]]@misc[["metacolors"]][[meta_col]]
      } else {
        col_pal <- colrr::col_pal("custom")
      }
    } else {
      col_pal <- colrr::col_pal("custom")
    }
  }

  if (length(SO) == 1 && !"legend.position" %in% names(theme_args)) {
    theme_args <- c(theme_args, list(legend.position = "none"))
  }

  if (!plot_strip) {
    theme_args <- c(theme_args, list(strip.background = ggplot2::element_blank(),
                                     strip.text.x = ggplot2::element_blank()))
  }

  col_pal <- colrr::make_col_pal(col_vec = col_pal,
                                 fct_lvls = if (is.factor(data[[color_aes]])) levels(data[[color_aes]]) else sort(unique(data[[color_aes]])),
                                 missing_fct_to_na = ifelse("missing_fct_to_na" %in% names(col_pal_args), col_pal_args[["missing_fct_to_na"]], T),
                                 col_pal_args = col_pal_args[-which(names(col_pal_args) %in% c("name", "missing_fct_to_na"))])

  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(meta_col), y = feature, color = !!rlang::sym(color_aes))) +
    theme +
    ggplot2::scale_color_manual(values = col_pal, na.value = col_na) +
    Gmisc::fastDoCall(ggplot2::facet_wrap, args = c(list(facets = ggplot2::vars(factor(feature_split, features))), facetting_args)) +
    ggplot2::labs(y = ylab) +
    Gmisc::fastDoCall(ggplot2::theme, args = theme_args) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = axis_expansion_y_mult))

  bckgr_col <- brathering::gg_get_theme_element(plot, element = "plot.background")@fill

  # geom1_layer <- ggplot2::geom_point(size = pt_size, position = geom1pos)

  geom2_layer <- suppressWarnings(Gmisc::fastDoCall(
    geom2,
    args = c(list(
      position = ggplot2::position_dodge(width = dodgewidth),
      fill = bckgr_col
    ), geom2_args)
  ))

  geom3_layer <- suppressWarnings(Gmisc::fastDoCall(
    geom3,
    args = c(list(
      position = ggplot2::position_dodge(width = dodgewidth),
      fill = bckgr_col,
      data = datageom3
    ), geom3_args)
  ))

  plot <- plot + switch(
    plot_first,
    geom1 = list(geom1_layer, geom2_layer, geom3_layer),
    geom2 = list(geom2_layer, geom3_layer, geom1_layer)
  )

  if (expr_freq_size>0) {
    plot <- plot +
      ggplot2::geom_text(data = stat,
                         color = brathering:::bw_txt(bckgr_col),
                         ggplot2::aes(label = !!rlang::sym(ifelse(expr_freq_pct,
                                                                  "pct.expr.adjust.pct",
                                                                  "pct.expr.adjust")),
                                      y = max.feat.expr + expr_freq_hjust),
                         position = ggplot2::position_dodge(width = 0.75),
                         size = expr_freq_size, show.legend = F,
                         angle = expr_freq_angle)
  }

  if (add_pwc) {
    plot <- plot + ggpubr::geom_pwc(label = "p.signif", tip.length = 0)
  }


  return(plot)
}

