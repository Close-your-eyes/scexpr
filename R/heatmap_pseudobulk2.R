#' Plot a pseudobulk expression heatmap across groups
#'
#' Generates a heatmap (or dot plot) of average gene expression across groups
#' (e.g. clusters or cell types) from one or more Seurat objects. Expression
#' values are calculated on pseudobulked groups, optionally scaled, and can be
#' filtered to display marker genes or user-defined features.
#'
#' @param SO A Seurat object or a list of Seurat objects.
#'
#' @param meta_col Character vector specifying the metadata column(s) used to
#' define groups. If `NULL`, the current `Idents(SO)` are used.
#'
#' @param levels_calc Character vector (or list of vectors when multiple Seurat
#' objects are supplied) defining which group levels are included in expression
#' calculations. If `NULL`, all levels are used. The order determines the group
#' order during calculations and clustering.
#'
#' @param levels_plot Character vector (or list of vectors) defining which
#' groups are displayed. Must be a subset of `levels_calc`. This affects only
#' plotting and ordering, not expression calculations.
#'
#' @param assay Assay from which expression values are obtained.
#'
#' @param features Character vector of features (genes) to plot. Alternatively,
#' a named list can be supplied to define feature groups. If `NULL`,
#' `features_topn` determines which genes are selected.
#'
#' @param features_topn Integer specifying the number of top marker genes to
#' display per group when `features` is `NULL`. If both `features` and
#' `features_topn` are `NULL`, all detected features are plotted.
#'
#' @param topn_metric Ranking metric used for marker selection. One or more of
#' `"padj"`, `"auc"`, or `"logFC"` may be supplied.
#'
#' @param min_pct Minimum percentage of expressing cells required for a gene to
#' be considered during marker selection.
#'
#' @param max_padj Maximum adjusted p-value for a gene to be considered a
#' marker during feature selection.
#'
#' @param min_pct_force Logical. If `TRUE`, genes below `min_pct` are removed
#' even when features are provided manually.
#'
#' @param featuregroup_style How feature groups should be displayed. Can include
#' `"facet"` and/or `"color"`.
#'
#' @param featuregroup_col_name Legend title for feature group colours.
#'
#' @param featuregroup_col_pal Colour palette passed to
#' `colrr::col_pal()`.
#'
#' @param feature_order Method used to order features:
#' * `"custom"`: preserve supplied order,
#' * `"hclust"`: hierarchical clustering,
#' * `"none"`: no explicit ordering.
#'
#' @param group_order Method used to order groups:
#' * `"custom"`: preserve supplied order,
#' * `"hclust"`: hierarchical clustering,
#' * `"none"`: no explicit ordering.
#'
#' @param topn_ties Logical; if `TRUE`, tied genes are retained, potentially
#' displaying more than `features_topn` genes per group.
#'
#' @param dotplot Logical. If `TRUE`, plot a dot plot where dot size indicates
#' the percentage of expressing cells instead of a heatmap.
#'
#' @param dotsize_range Numeric vector of length two defining the minimum and
#' maximum dot sizes.
#'
#' @param fill Colour palette used for expression values. `"auto"` selects a
#' default diverging palette.
#'
#' @param color Border colour for heatmap tiles or dots. `"NA"` disables
#' borders. `"auto"` chooses a suitable default.
#'
#' @param scale Expression scaling method:
#' * `"zscore"`: gene-wise z-score,
#' * `"1"`: scale each gene to [-1, 1],
#' * `"none"`: no scaling.
#'
#' @param featurelabels Feature labels to display. `NULL` plots all labels,
#' `""` plots none, `"auto"` suppresses labels for large heatmaps. A named
#' vector may be used to provide aliases.
#'
#' @param featurelabels_repel Logical; repel overlapping feature labels.
#'
#' @param featuresitalic Logical; display feature labels in italic font.
#'
#' @param color_linewidth Line width of tile or dot borders.
#'
#' @param legendbreaks Legend break positions. May be `"auto"`,
#' `"minmidmax"`, a numeric vector, or a single number.
#'
#' @param legendlabels Labels corresponding to `legendbreaks`.
#'
#' @param colorsteps Controls discretization of the colour scale. Can be
#' `NULL`, `"auto"`, a numeric vector, or a single integer.
#'
#' @param colorsteps_nice Logical; use aesthetically pleasing colour breaks.
#'
#' @param axes_flip Logical; transpose the heatmap axes.
#'
#' @param group_seplines Logical; draw separator lines between feature groups.
#'
#' @param seplines_args Named list of additional arguments passed to
#' `geom_hline()`.
#'
#' @param legend_fill_args Named list of arguments forwarded to
#' `guide_colorbar()` or `guide_colorsteps()`.
#'
#' @param legend_size_args Named list of arguments forwarded to
#' `guide_legend()` for the dot-size legend.
#'
#' @param theme ggplot2 theme applied to the plot.
#'
#' @param theme_args Named list of additional arguments passed to
#' `ggplot2::theme()`.
#'
#' @param repel_args Named list of arguments controlling feature label
#' repulsion.
#'
#' @param sec_axis Logical; if `TRUE`, add a secondary y-axis showing gene
#' names converted with `scexpr::convert_gene_identifier()`.
#'
#' @param convert_gene_identifier_args Named list of arguments passed to
#' `scexpr::convert_gene_identifier()`.
#'
#' @param pvals Optional data frame containing p-values for annotation.
#'
#' @param pval_features Features for which p-value annotations should be shown.
#'
#' @param pval_max Maximum p-value displayed.
#'
#' @param pval_symnum_args Arguments passed to `symnum()` for significance
#' symbol generation.
#'
#' @param pval_filter Strategy for filtering p-value annotations.
#'
#' @param pval_logfc Name of the log-fold change column in `pvals`.
#'
#' @param pval_text_args Named list of arguments passed to the text annotation
#' layer.
#'
#' @return A named list containing:
#' \describe{
#'   \item{plot}{The ggplot object.}
#'   \item{data}{Processed plotting data after filtering and scaling.}
#'   \item{complete_data}{Complete marker statistics prior to filtering.}
#' }
#'
#' @importFrom zeallot %<-%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' hm <- heatmap_pseudobulk2(
#'   SO,
#'   meta_col = "celltype",
#'   features_topn = 10
#' )
#'
#' hm$plot
#'
#' hm$plot +
#'   scale_size_continuous(
#'     breaks = c(1, 5, 10, 20, 50),
#'     range = c(1, 10)
#'   )
#' }
heatmap_pseudobulk2 <- function(SO,
                                meta_col = NULL,
                                levels_calc = NULL,
                                levels_plot = NULL,
                                assay = "RNA",
                                min_pct = 20,
                                max_padj = 0.05,
                                min_pct_force = T,
                                features = NULL,
                                featuregroup_style = c("facet", "color"),
                                featuregroup_col_name = "",
                                featuregroup_col_pal = "custom_light",
                                feature_order = c("custom", "hclust", "none"),
                                group_order = c("hclust", "none", "custom"),
                                topn_metric = c("padj", "auc", "logFC"),
                                topn_ties = F,
                                dotplot = F,
                                dotsize_range = c(2,7),
                                fill = "..auto..",
                                color = "..auto..",
                                scale = c("zscore", "1", "none"),
                                features_topn = NULL,
                                featurelabels = "..auto..",
                                featurelabels_repel = F,
                                featuresitalic = T,
                                color_linewidth = 0.2,
                                legendbreaks = "..auto..",
                                legendlabels = "..auto..",
                                colorsteps = "..auto..",
                                colorsteps_nice = T,
                                axes_flip = F,
                                group_seplines = F,
                                seplines_args = list(),
                                theme = colrr::theme_material(white = T),
                                legend_fill_args = list(
                                  label.theme = ggplot2::element_text(size = 10),
                                  title.theme = ggplot2::element_text(size = 10),
                                  title.position = "top",
                                  title = "..auto..",
                                  title.hjust = 0.5,
                                  barwidth = 1,
                                  barheight = 8,
                                  order = 1
                                ),
                                legend_size_args = list(
                                  label.theme = ggplot2::element_text(size = 10),
                                  title.theme = ggplot2::element_text(size = 10),
                                  title.position = "top",
                                  title = "transcription\nfrequency\n[%]",
                                  title.hjust = 0.5,
                                  label.position = "right",
                                  order = 2,
                                  ncol = NULL,
                                  nrow = NULL,
                                  override.aes = list(color = "..auto..", fill = "..auto..")
                                  #size = c(2,4, 10))
                                ),
                                theme_args = list(
                                  axis.title = ggplot2::element_blank(),
                                  panel.grid = ggplot2::element_blank(),
                                  axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                                  legend.position = "right",
                                  legend.direction = "vertical",
                                  legend.box = "horizontal"
                                ),
                                repel_args = list(featurelabels_width = 0.2,
                                                  featurelabels_nudhe_x = -1),
                                sec_axis = F,
                                convert_gene_identifier_args = list(ident_in = "SYMBOL", ident_out = "GENENAME"),
                                pvals = NULL,
                                pval_features = NULL,
                                pval_max = 0.01,
                                pval_symnum_args = list(cutpoints = c(0, 0.0001, Inf),
                                                        symbols = c("*", "ns")),
                                pval_filter = c("top", "pos_fc"),
                                pval_logfc = "logFC",
                                pval_text_args = list(size = 5, vjust = 0.75)) {

  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (!requireNamespace("presto", quietly = T)) {
    pak::pak("immunogenomics/presto")
  }
  if (!requireNamespace("fcexpr", quietly = T)) {
    pak::pak("Close-your-eyes/fcexpr")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }

  feature_order <- rlang::arg_match(feature_order)
  group_order <- rlang::arg_match(group_order)
  topn_metric <- rlang::arg_match(topn_metric, multiple = T)
  featuregroup_style <- rlang::arg_match(featuregroup_style, multiple = T)
  pval_filter <- rlang::arg_match(pval_filter)

  if (!is.null(features_topn)) {
    features_topn <- max(1, features_topn)
  }
  scale <- rlang::arg_match(scale)

  if (!is.null(features) && !is.null(features_topn)) {
    message("Providing features and features_topn may not be meaningful.")
  }

  if ("title" %in% names(legend_fill_args) && legend_fill_args[["title"]] == "..auto..") {
    if (scale == "zscore") {
      legend_fill_args[["title"]] <- "transcription\nlevel\n[z-score]"
    } else if (scale == "none") {
      legend_fill_args[["title"]] <- "transcription\nlevel\n[log (UMI)]"
    } else if (scale == "1") {
      legend_fill_args[["title"]] <- "transcription\nlevel"
    }
  }
  if ("legend.direction" %in% names(theme_args)) {
    theme_args$legend.direction <- match.arg(theme_args$legend.direction, c("horizontal", "vertical"))
    if (theme_args$legend.direction == "horizontal" &&
        "barwidth" %in% names(legend_fill_args) &&
        "barheight" %in% names(legend_fill_args)) {
      temp <- legend_fill_args$barwidth
      legend_fill_args$barwidth <- legend_fill_args$barheight
      legend_fill_args$barheight <- temp
    }
  }

  SO <- scexpr:::check.SO(SO = SO, assay = assay)
  assay <- Seurat::DefaultAssay(SO[[1]])
  feature_groups <- NULL
  if (!is.null(features) && is.list(features)) {
    feature_groups <- features
    features <- unlist(features)
  }

  c(SO,
    meta_col,
    levels_calc,
    levels_plot) %<-% check_meta_and_levels(
      SO = SO,
      meta_col = meta_col,
      levels_calc = levels_calc,
      levels_plot = levels_plot)

  c(wil_auc_raw,
    wil_auc,
    features,
    feature_order) %<-% determine_features(
      SO = SO,
      assay = assay,
      features = features,
      levels_plot = levels_plot,
      meta_col = meta_col,
      feature_order = feature_order)
  # groups not present in levels_plot become NA in wil_auc

  # optional scaling
  # issue: scaling before vs. after selection for topn features
  # do here before feature selection to not disturb zscoring result e.g.
  # but what if features are selected manually with features argument?
  # for scaling of each feature it is only relevant that all groups are still included
  # determine_features does not remove groups but features only
  if (scale != "none") {
    wil_auc <- dplyr::mutate(wil_auc, avgExprScale = dplyr::case_when(
      scale == "zscore" ~ as.vector(scale(avgExpr)),
      scale == "1" ~ brathering::scale2(avgExpr, min = -1, max = 1),
      .default = avgExpr  # fallback (optional), not needed here but anyway
    ), .by = feature)
    values <- "avgExprScale"
  } else {
    values <- "avgExpr"
  }

  # do this here but not in heatmap_long_df to allow pre-filter for pct_in and padj
  if (!is.null(features_topn)) {
    select <-
      wil_auc |>
      dplyr::mutate(padj = -padj) |> # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc")),
      dplyr::filter(pct_in >= min_pct) |>
      dplyr::filter(padj <= max_padj) |>
      dplyr::slice_max(order_by = avgExpr, n = 1, by = feature) |>
      dplyr::slice_max(order_by = tibble::tibble(!!!rlang::syms(as.list(topn_metric))), n = features_topn,
                       by = group, with_ties = topn_ties) |>
      # filter here but not in determine_features to not bias feature selection by levels_plot
      dplyr::filter(group %in% unlist(levels_plot)) |>
      dplyr::pull(feature)
    wil_auc <- wil_auc[which(wil_auc[["feature"]] %in% select),,drop = F]
  } else if (is.null(features_topn) && min_pct_force) {
    select <-
      wil_auc |>
      dplyr::filter(pct_in >= min_pct) |>
      dplyr::filter(group %in% unlist(levels_plot)) |>
      dplyr::pull(feature)
    rm_feat <- setdiff(as.character(unique(wil_auc$feature)), select)
    if (length(rm_feat)) {
      message("features rm due to min_pct and min_pct_force (", length(rm_feat), ")")
      #message(paste(rm_feat, collapse = ", "))
    }

    wil_auc <- wil_auc[which(wil_auc[["feature"]] %in% select),,drop = F]
  }
  # filter here after feature selection
  wil_auc <- dplyr::filter(wil_auc, group %in% unlist(levels_plot))

  featuregroup <- NULL
  if (!is.null(feature_groups)) {
    feat_group <- stats::setNames(utils::stack(feature_groups), c("feature", "featgroup"))
    wil_auc <- dplyr::left_join(wil_auc, feat_group, by = "feature")
    featuregroup <- "featgroup"
  }



  dotsizes <- if (dotplot) "pct_in" else NULL
  # fcexpr::

  plot <- fcexpr::heatmap_long_df(df = wil_auc,
                                  groups = "group",
                                  features = "feature",
                                  featuregroup = featuregroup,
                                  featuregroup_style = featuregroup_style,
                                  featuregroup_col_name = featuregroup_col_name,
                                  featuregroup_col_pal = featuregroup_col_pal,
                                  values = values,
                                  dotsizes = dotsizes,
                                  dotsize_range = dotsize_range,
                                  features_topn = NULL,
                                  # irrelevant as features_topn is handled above
                                  #topn_metric
                                  #topn_ties
                                  fill = fill,
                                  color = color,
                                  scale = "none",
                                  featurelabels = featurelabels,
                                  featurelabels_repel = featurelabels_repel,
                                  featuresitalic = featuresitalic,
                                  color_linewidth = color_linewidth,
                                  legendbreaks = legendbreaks,
                                  legendlabels = legendlabels,
                                  colorsteps = colorsteps,
                                  colorsteps_nice = colorsteps_nice,
                                  axes_flip = axes_flip,
                                  group_seplines = group_seplines,
                                  seplines_args = seplines_args,
                                  legend_fill_args = legend_fill_args,
                                  legend_size_args = legend_size_args,
                                  theme_args = theme_args,
                                  repel_args = repel_args,
                                  theme = theme,
                                  heatmap_ordering_args = list(feature_order = feature_order,
                                                               group_order = group_order),
                                  pvals = pvals,
                                  pval_features = pval_features,
                                  pval_max = pval_max,
                                  pval_symnum_args = pval_symnum_args,
                                  pval_filter = pval_filter,
                                  pval_logfc = pval_logfc,
                                  pval_text_args = pval_text_args)


  # if (!is.null(feature_groups)) {
  #   marker_df <- stats::setNames(utils::stack(feature_groups), c("gene", "cell_type"))
  #   color_conv <- stats::setNames(as.character(colrr::col_pal("custom_light", n = length(unique(marker_df$cell_type)))),
  #                                 unique(marker_df$cell_type))
  #   marker_df$color <- color_conv[marker_df$cell_type]
  #   colman <- setNames(marker_df$color, marker_df$cell_type)
  #
  #   plot <- plot +
  #     ggplot2::scale_y_discrete(labels = function(y) color_labels(y, stats::setNames(marker_df$color, marker_df$gene))) +
  #     ggplot2::theme(axis.text.y = ggtext::element_markdown()) +
  #     ggplot2::geom_point(
  #       data = data.frame(
  #         x = 0,
  #         y = 0,
  #         yaxis = marker_df$cell_type
  #       ),
  #       ggplot2::aes(x = x, y = y, color = yaxis),
  #       inherit.aes = F
  #     ) +
  #     ggplot2::scale_color_manual(
  #       name = "",  #"Cell type",
  #       values = colman[!duplicated(colman)]
  #     )
  # }

  if (sec_axis) {
    plot <- add_sec_axis(
      plot = plot,
      convert_gene_identifier_args = convert_gene_identifier_args
    )
  }

  return(list(
    plot = plot,
    data = dplyr::mutate(wil_auc, padj = -padj),
    complete_data = wil_auc_raw
  ))
}

check.levels <- function(SO, meta_col, levels = NULL, append_by_missing = F) {

  if (any(is.null(levels)) || any(is.na(levels))) {
    levels <- as.character(unique(SO@meta.data[,meta_col,drop=T]))
    # problem: converts 01, 02, 03 to 1,2,3
    # if (suppressWarnings(!any(is.na(as.numeric(levels))))) {
    #   levels <- as.character(sort(as.numeric(unique(levels))))
    # }
  } else {
    if (any(!levels %in% unique(SO@meta.data[,meta_col,drop=T]))) {
      print(paste0("levels not found in meta_col: ", paste(levels[which(!levels %in% unique(SO@meta.data[,meta_col,drop=T]))], collapse = ", ")))
    }
    levels <- as.character(unique(levels[which(levels %in% unique(SO@meta.data[,meta_col,drop=T]))]))
  }
  if (append_by_missing) {
    levels <- c(levels, as.character(unique(SO@meta.data[,meta_col,drop=T])[which(!as.character(unique(SO@meta.data[,meta_col,drop=T])) %in% levels)]))
  }
  return(levels)
}

# feature_order_fun <- function(wil_auc_raw,
#                               features_topn,
#                               break.ties,
#                               topn_metric,
#                               feature_selection_strategy,
#                               features,
#                               min_pct,
#                               max_padj,
#                               levels_plot) {
#   # select multiple columns for ordering, in order to break ties; especially if padj is used as primary element to order
#   if (break.ties) {
#     all_metrics <- c("padj", "logFC", "auc")
#     second_metric <- all_metrics[which(!all_metrics %in% topn_metric)][1]
#     third_metric <- all_metrics[which(!all_metrics %in% topn_metric)][2]
#   }
#
#   if (feature_selection_strategy == "2") {
#     # this is better for when features are provided; yields better order on y-axis
#     features2 <-
#       wil_auc_raw |>
#       dplyr::filter(feature %in% features) |>
#       dplyr::filter(pct_in >= min_pct) |>
#       dplyr::filter(padj <= max_padj) |>
#       dplyr::mutate(padj = -padj) |> # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc"))
#       dplyr::group_by(feature) |>
#       dplyr::slice_max(order_by = avgExpr, n = 1) |>
#       dplyr::ungroup() |>
#       dplyr::mutate(group = factor(group, levels = unlist(levels_plot))) |>
#       dplyr::group_by(group)
#     if (break.ties) {
#       features3 <-
#         features2 |>
#         dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn_metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = features_topn)
#     } else {
#       features3 <-
#         features2 |>
#         dplyr::slice_max(order_by = !!rlang::sym(topn_metric), n = features_topn, with_ties = T)
#     }
#     features3 <-
#       features3 |>
#       dplyr::ungroup() |>
#       dplyr::arrange(group, avgExpr) #!!rlang::sym(topn_metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)
#
#     ## and maybe add additional ones, but how?
#     features_check <-
#       features3 |>
#       dplyr::group_by(feature) |>
#       dplyr::slice_max(order_by = avgExpr, n = 1) |>
#       dplyr::ungroup() |>
#       dplyr::group_by(group) |>
#       dplyr::count()
#
#   } else if (feature_selection_strategy == "1") {
#     ## alternative feature selection; this procedure avoid that a diff gene is "stolen" by another group (due to max avgExpr) where it does not lie within features_topn
#     ## this would cause that a feature actually relevant (diff expressed by 2 group does not appear on the heatmap)
#     ## it may happen that a few more features per group are plotted
#     ## for very large heatmaps (overview with 100 of features) is may look ugly; then feature_selection_strategy 2 may be better suited
#     features2 <-
#       wil_auc_raw |>
#       dplyr::filter(feature %in% features) |>
#       dplyr::filter(pct_in >= min_pct) |>
#       dplyr::filter(padj <= max_padj) |>
#       dplyr::mutate(padj = -padj) |> # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc"))
#       dplyr::mutate(group = factor(group, levels = unlist(levels_plot))) |>
#       dplyr::group_by(group)
#
#     if (break.ties) {
#       features3 <-
#         features2 |>
#         dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn_metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = features_topn) |>
#         dplyr::ungroup()
#
#       ## check for number of features per group
#       ## and maybe add additional ones, but how?
#       features_check <-
#         features3 |>
#         dplyr::group_by(feature, .drop = F) |>
#         dplyr::slice_max(order_by = avgExpr, n = 1) |>
#         dplyr::ungroup() |>
#         dplyr::group_by(group) |>
#         dplyr::count()
#
#       # try to make sure features_topn is reach for every group
#       max_round <- min(features_topn, 20)
#       n <- 1
#       while(any(features_check$n < features_topn) && n <= max_round) {
#         for (i in features_check[which(features_check$n<features_topn),"group",drop=T]) {
#           nplus <- features_topn - features_check[which(features_check$group == i),"n",drop=T]
#           # in case factor level i is missing in features_check
#           if (length(nplus) == 0) {
#             nplus <- features_topn
#           }
#           features_plus <-
#             dplyr::anti_join(features2, features3, by = dplyr::join_by(feature, group, avgExpr, logFC, statistic, auc, pval, padj, pct_in, pct_out)) |>
#             dplyr::filter(group == i) |>
#             dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn_metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = nplus)
#
#           features3 <-
#             rbind(features3, features_plus) |>
#             dplyr::distinct()
#
#           features_check <-
#             features3 |>
#             dplyr::group_by(feature) |>
#             dplyr::slice_max(order_by = avgExpr, n = 1) |>
#             dplyr::ungroup() |>
#             dplyr::group_by(group) |>
#             dplyr::count()
#         }
#         n<-n+1
#       }
#       features3 <-
#         features3 |>
#         dplyr::arrange(group, avgExpr) # avgExpr instead of metric?! # -!!rlang::sym(topn_metric), -!!rlang::sym(second_metric), -!!rlang::sym(third_metric)
#
#     } else {
#
#       features3 <-
#         features2 |>
#         dplyr::slice_max(order_by = !!rlang::sym(topn_metric), n = features_topn, with_ties = T) |>
#         dplyr::ungroup() |>
#         dplyr::arrange(group, avgExpr) # avgExpr instead of metric?! #-!!rlang::sym(topn_metric)
#
#       features_check <-
#         features3 |>
#         dplyr::group_by(feature) |>
#         dplyr::slice_max(order_by = avgExpr, n = 1) |>
#         dplyr::ungroup() |>
#         dplyr::group_by(group) |>
#         dplyr::count()
#     }
#   }
#
#   print(as.data.frame(features_check))
#   return(features3)
# }


check_meta_and_levels <- function(SO, meta_col, levels_calc, levels_plot) {

  ## checking meta_col, levels_plot and levels_calc
  if (is.null(meta_col)) {
    SO <- lapply(SO, function(x) {
      x@meta.data$idents <- Seurat::Idents(x)
      return(x)
    })
    meta_col <- rep("idents", length(SO))
  } else {
    if (length(meta_col) != length(SO)) {
      meta_col <- brathering::recycle(meta_col, seq_along(SO))
    }
    if (any(unlist(mapply(x = SO, y = meta_col, function(x,y) !y %in% names(x@meta.data), SIMPLIFY = F)))) {
      stop("One of meta_col not found in respective SO.")
    }
    if (any(unlist(mapply(x = SO, y = meta_col, function(x,y) is.numeric(x@meta.data[,y,drop=T]), SIMPLIFY = F)))) {
      # maybe handle internally
      stop("One of meta_col across SO is numeric. Please make it factor or character.")
    }
  }



  # this will only apply if length(SO)==1 and levels_plot and levels_calc are vectors
  if (!is.null(levels_plot) && !is.list(levels_plot)) {
    levels_plot <- list(levels_plot)
  }
  if (!is.null(levels_calc) && !is.list(levels_calc)) {
    levels_calc <- list(levels_calc)
  }


  if (!is.null(levels_plot) && !is.null(levels_calc)) {
    levels_plot <- mapply(x = levels_plot, y = levels_calc, function(x,y) {
      if (any(!x %in% y)) {
        print("levels_plot has more levels than levels_calc. levels_plot is reduced to levels_calc.")
        x <- x[which(x %in% y)]
      }
      return(x)
    }, SIMPLIFY = F)
  }

  if (!is.null(levels_calc)) {
    if (length(levels_calc) != length(SO)) {
      stop("levels_calc has to have same length as SO.")
    }
  } else {
    levels_calc <- rep(NA, length(SO))
  }
  if (!is.null(levels_plot)) {
    if (length(levels_plot) != length(SO)) {
      stop("levels_plot has to have same length as SO.")
    }
  } else {
    levels_plot <- rep(NA, length(SO))
  }

  # this will replace NA with all available levels
  levels_calc <- purrr::pmap(list(x = SO, y = levels_calc, z = meta_col), function(x,y,z) check.levels(SO = x, meta_col = z, levels = y, append_by_missing = F))
  levels_plot <- purrr::pmap(list(x = SO, y = levels_plot, z = meta_col), function(x,y,z) check.levels(SO = x, meta_col = z, levels = y, append_by_missing = F))

  # subset levels for calculation
  SO <- purrr::pmap(list(x = SO, y = levels_calc, z = meta_col), function(x,y,z) {
    if (length(y) < length(unique(x@meta.data[,z,drop=T]))) {
      x <- subset(x, cells = rownames(x@meta.data[,z,drop=F][which(x@meta.data[,z] %in% y),,drop=F]))
    }
    return(x)
  })

  # this problem was discovered by chance (":" are replaced by Seurat::AverageExpression by "_"); maybe other symbols will also cause problems.
  ## use make portable?
  for (i in length(SO)) {
    if (any(grepl(":", unique(SO[[i]]@meta.data[,meta_col[[i]]])))) {
      message("Found colon (:) in factor levels of meta_col. Will replace those with underscores (_).")
      SO[[i]]@meta.data[,meta_col[[i]]] <- gsub(":", "_", SO[[i]]@meta.data[,meta_col[[i]]])
      levels_plot[[i]] <- gsub(":", "_", levels_plot[[i]])
      levels_calc[[i]] <- gsub(":", "_", levels_calc[[i]])
    }
  }

  # make cluster names unique if there is intersection
  ## this will prohibit to define cluster order
  if (length(SO) > 1) {
    if (length(Reduce(intersect, levels_plot)) > 1 || length(Reduce(intersect, levels_calc)) > 1) {
      SO <- purrr::pmap(list(x = SO, y = meta_col, z = names(SO)), function(x,y,z) {
        x@meta.data[,y] <- paste0(z, "___", x@meta.data[,y,drop=T])
        return(x)
      })
    }
    if (length(Reduce(intersect, levels_plot)) > 1) {
      levels_plot <- mapply(x = names(SO), y = levels_plot, function(x,y) paste0(x, "___", y), SIMPLIFY = F)
    }
    if (length(Reduce(intersect, levels_calc)) > 1) {
      levels_calc <- mapply(x = names(SO), y = levels_calc, function(x,y) paste0(x, "___", y), SIMPLIFY = F)
    }
  }
  return(list(SO, meta_col, levels_calc, levels_plot))
}


determine_features <- function(SO,
                               assay,
                               features,
                               levels_plot,
                               meta_col,
                               feature_order) {

  wil_auc_raw <- find_all_marker(obj = SO,
                                 meta_col = meta_col,
                                 assay = assay,
                                 calc_avglog2fc = F, # for speed? or change fun to avg_log2FC
                                 eps = 1e-6,
                                 features = features)

  if (is.null(features)) {
    if (feature_order == "none") {
      feature_order <- "custom"
    }
    features <- wil_auc_raw$feature
  }

  wil_auc <-
    wil_auc_raw |>
    dplyr::filter(feature %in% features) |>
    # pre-filtering for levels_plot before selecting features_topn may cause deviating topn features depend upon the groups included
    # so somehow the selection procedure still has relativeness in it
    #dplyr::filter(group %in% unlist(levels_plot)) |>
    dplyr::mutate(group = factor(group, levels = unlist(levels_plot))) |>
    dplyr::mutate(feature = factor(feature, levels = unique(features)))

  return(list(wil_auc_raw, wil_auc, features, feature_order))
}


add_sec_axis <- function(plot, convert_gene_identifier_args) {
  axis.df <-
    Gmisc::fastDoCall(scexpr::convert_gene_identifier,
                      args = c(list(idents = levels(plot[["data"]]$feature)), convert_gene_identifier_args)) |>
    dplyr::mutate(y = dplyr::row_number())
  axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = GENENAME)) +
    ggplot2::geom_text(hjust = 0) +
    ggplot2::scale_x_continuous(
      limits = c(0, 2),
      expand = c(0, 0),
      breaks = NULL,
      labels = NULL,
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.4, nrow(axis.df) + 0.4),
      expand = c(0, 0),
      breaks = NULL,
      labels = NULL,
      name = NULL
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  plot <- cowplot::plot_grid(
    plot + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")),
    axis,
    align = "h",
    axis = "tb",
    nrow = 1,
    rel_widths = c(0.6,0.4)
  )
  return(plot)
}

