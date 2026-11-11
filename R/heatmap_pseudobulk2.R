#' Plot a pseudobulk expression heatmap across groups
#'
#' Generates a heatmap (or dot plot) of average gene expression across groups
#' (e.g. clusters or cell types) from one or more Seurat objects. The function
#' averages expression across cells in each group and can scale those averages
#' and filter them to marker genes or user-defined features.
#'
#' @param SO A Seurat object or a list of Seurat objects.
#'
#' @param meta_col Name of the metadata column used to define groups. For
#' multiple objects, supply a shared column name or one name per object. If
#' `NULL`, use each object's active identities.
#'
#' @param levels_calc Character vector (or list of vectors when multiple Seurat
#' objects are supplied) defining which groups are included when calculating
#' group means and marker statistics. If `NULL`, all levels are used.
#' Objects with no selected cells in the assay's `data` layer are skipped.
#' Marker statistics require at least two calculated groups.
#'
#' @param levels_plot Character vector (or list of vectors) defining which
#' groups are displayed. `NULL` displays all calculated groups with a
#' nonmissing metadata label. Entries outside `levels_calc` are omitted.
#' With `group_order = "none"`, the plot follows the order in `levels_plot`.
#' Group means, marker statistics, and scaling are calculated before filtering
#' to the displayed groups. The displayed groups can still affect which
#' features are retained.
#'
#' @param assay Assay whose `data` layer supplies expression values. The assay
#' must be present in every object.
#'
#' @param features Character vector of features (genes) to plot. Alternatively,
#' supply a named list of feature vectors to define feature groups. If `NULL`,
#' start with detected features shared across the supplied objects.
#'
#' @param featuregroup_by_peak Logical. When `TRUE`, build a named feature
#' list from the calculation groups where `find_all_marker()` reports positive
#' `logFC` and `padj` below 1e-4. Plot features with at least one qualifying
#' group. Label color follows the highest mean expression among qualifying
#' groups. A supplied named list takes precedence.
#'
#' @param featuregroup_membership When using automatic feature groups,
#' `"all"` puts each feature in all qualifying groups. `"top"` keeps the
#' qualifying group with the highest unscaled mean expression. Ties follow
#' the `levels_calc` order.
#'
#' @param features_topn Number of top marker genes to select per group. Each
#' eligible gene is assigned to its group or tied groups with the highest mean
#' expression among eligible rows, then ranked within each group by
#' `topn_metric`. If `features` is supplied, selection is limited to those
#' genes. If `NULL`, retain candidate features subject to `min_pct_force`.
#'
#' @param topn_metric One or more columns used in the supplied order to rank
#' candidates for `features_topn`. `"padj"` ranks smaller values first,
#' while `"auc"` and `"logFC"` rank larger values first.
#'
#' @param min_pct Minimum percentage of expressing cells within a group. Used
#' to filter candidates when `features_topn` is set, or to filter by expression
#' in the displayed groups when `min_pct_force` applies.
#'
#' @param max_padj Maximum adjusted p-value for a feature to qualify in a
#' group when `features_topn` is set.
#'
#' @param min_pct_force Logical. When `features_topn` is `NULL` and this is
#' `TRUE`, retain a feature only if it reaches `min_pct` in at least one
#' displayed group. This also applies to manually supplied features.
#'
#' @param featuregroup_style How feature groups should be displayed. Can include
#' `"facet"` and/or `"color"`. The default uses both. For genes assigned to
#' multiple feature groups, use the color of the matching cluster with the
#' highest mean expression, including calculated clusters omitted from the
#' plot. Ties and assignments without matching clusters use the list order.
#'
#' @param featuregroup_col_name Legend title for feature group colours.
#'
#' @param featuregroup_col_pal Colour palette passed to
#' `colrr::col_pal()`.
#'
#' @param feature_order Method used to order features.
#' * `"custom"` sorts features by the group where their plotted expression is
#'   highest, following the group order, then by that expression in ascending
#'   order.
#' * `"hclust"` uses hierarchical clustering of feature expression profiles.
#' * `"none"` keeps the order supplied in `features`. If `features` is `NULL`,
#'   it keeps the order of features in the marker results.
#'
#' @param group_order Method used to order groups.
#' * `"hclust"` uses hierarchical clustering of group expression profiles.
#' * `"none"` keeps the order in `levels_plot`.
#' * `"custom"` sorts groups by the feature where their plotted expression is
#'   highest, following the feature order, then by that expression in
#'   ascending order. With feature facets, it follows the visible feature
#'   order from the bottom panel upward, or from the left panel rightward
#'   when `axes_flip = TRUE`.
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
#' @param fill Vector of colours or a palette name for expression values.
#' `"..auto.."` uses a reversed 11-colour `"RdBu"` palette.
#'
#' @param color Border colour for heatmap tiles or dots. `"NA"` disables
#' borders. `"..auto.."` uses grey borders for tile heatmaps with at most 100
#' feature levels and omits borders for larger tile heatmaps and dot plots.
#'
#' @param scale Expression scaling method applied to each gene's group means.
#' * `"zscore"` standardizes the values within each gene.
#' * `"1"` rescales each gene's values to the range from minus one to one.
#' * `"none"` uses the unscaled group means.
#'
#' @param featurelabels Feature labels to display. `NULL` plots all labels,
#' `""` plots none, and `"..auto.."` omits labels when more than 200 features
#' are plotted. Supply a character vector of feature IDs to label a subset.
#' The names of this vector, when supplied, are the display labels.
#'
#' @param featurelabels_repel Logical. Draw feature labels in a separate panel
#' with repulsion to reduce overlap.
#'
#' @param featuresitalic Logical; display feature labels in italic font.
#'
#' @param color_linewidth Line width of tile borders or stroke width of dot
#' borders.
#'
#' @param legendbreaks Fill legend breaks. Use `"..auto.."` for automatic
#' breaks, `"minmidmax"` for breaks at the minimum, midpoint, and maximum,
#' a single number for the break count, or a numeric vector of positions.
#'
#' @param legendlabels Labels corresponding to `legendbreaks`, or
#' `"..auto.."` for automatic labels.
#'
#' @param colorsteps Use `NULL` for a continuous colour bar, `"..auto.."` for
#' automatic steps, a single number for the requested number of steps, or a
#' numeric vector of step boundaries.
#'
#' @param colorsteps_nice Logical; use aesthetically pleasing colour breaks.
#'
#' @param axes_flip Logical; transpose the heatmap axes.
#'
#' @param group_seplines Logical. Draw separator lines based on runs of
#' features assigned to the same group by their highest plotted expression.
#' This option is designed for `feature_order = "custom"`.
#'
#' @param seplines_args Named list of additional arguments passed to
#' `geom_hline()`.
#'
#' @param legend_fill_args Named list of arguments forwarded to
#' `guide_colorbar()` or `guide_colorsteps()`. Use `title = NULL` to suppress
#' the title.
#'
#' @param legend_size_args Named list of arguments forwarded to
#' `guide_legend()` for the dot-size legend.
#'
#' @param theme ggplot2 theme applied to the plot.
#'
#' @param theme_args Named list of additional arguments passed to
#' `ggplot2::theme()`.
#'
#' @param repel_args Named list controlling the label panel when
#' `featurelabels_repel` is `TRUE`. Use `featurelabels_width` for the relative
#' panel width and `featurelabels_nudge_x` for the horizontal nudge.
#'
#' @param sec_axis Logical. Add a panel on the right with feature IDs
#' converted by `scexpr::convert_gene_identifier()`.
#'
#' @param convert_gene_identifier_args Named list of arguments passed to
#' `scexpr::convert_gene_identifier()`.
#'
#' @param pvals Optional name of a p-value column in the plotting data, such
#' as `"pval"` or `"padj"`, used to annotate cells with significance symbols.
#'
#' @param pval_features Feature IDs eligible for p-value annotations. `NULL`
#' makes all plotted features eligible.
#'
#' @param pval_max Maximum p-value eligible for annotation.
#'
#' @param pval_symnum_args Arguments passed to `symnum()` for significance
#' symbol generation.
#'
#' @param pval_filter Strategy for filtering p-value annotations. `"top"`
#' selects the cells with the highest plotted value for each feature.
#' `"pos_fc"` selects cells with a positive value in the column named by
#' `pval_logfc`. In either case, `pval_features` and `pval_max` further
#' restrict which cells are annotated.
#'
#' @param pval_logfc Name of the effect-size column in the plotting data used
#' by `pval_filter = "pos_fc"`. The default `"logFC"` column from `presto`
#' contains differences in mean expression on the assay data scale.
#'
#' @param pval_text_args Named list of arguments passed to the text annotation
#' layer.
#' @param featuregroup_color_on_x_axis_text Logical. When features are grouped
#' and `featuregroup_style` includes `"color"`, color matching group labels on
#' the horizontal axis, or on the vertical axis if `axes_flip = TRUE`.
#'
#' @details
#' The function passes `feature_order` and `group_order` in the
#' `heatmap_ordering_args` list supplied to `fcexpr::heatmap_long_df()`.
#' Ordering takes place after feature selection and scaling. For a standard
#' heatmap, the feature and group factor levels in `plot$data` contain the plot
#' order. Hierarchical clustering requires at least two entries on the axis
#' being clustered. A single feature or displayed group keeps its position.
#'
#' When calculated group labels overlap across objects, the result prefixes
#' all group labels with their object names and adds suffixes to resolve any
#' remaining collisions.
#'
#' @return A named list with the following components.
#' \describe{
#'   \item{plot}{The ggplot object, with an additional panel when repelled
#'   labels or converted identifiers are requested.}
#'   \item{data}{Marker statistics after feature and group filtering and
#'   optional scaling, before plot ordering. Scaled group means are in
#'   \code{avgExprScale} when scaling is enabled.}
#'   \item{complete_data}{Marker statistics for the calculated groups and
#'   tested features before plot filtering and scaling.}
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
#'   features_topn = 10,
#'   dotplot = TRUE
#' )
#'
#' hm$plot
#'
#' hm$plot +
#'   ggplot2::scale_size_continuous(
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
                                                  featurelabels_nudge_x = -0.05),
                                sec_axis = F,
                                convert_gene_identifier_args = list(ident_in = "SYMBOL", ident_out = "GENENAME"),
                                pvals = NULL,
                                pval_features = NULL,
                                pval_max = 0.01,
                                pval_symnum_args = list(cutpoints = c(0, 0.0001, Inf),
                                                        symbols = c("*", "ns")),
                                pval_filter = c("top", "pos_fc"),
                                pval_logfc = "logFC",
                                pval_text_args = list(size = 5, vjust = 0.75),
                                featuregroup_color_on_x_axis_text = F,
                                featuregroup_by_peak = FALSE,
                                featuregroup_membership = c( "top", "all")) {

  scexpr:::.ensure_packages(c("presto", "fcexpr", "brathering", "cowplot", "colrr"))

  feature_order <- rlang::arg_match(feature_order)
  group_order <- rlang::arg_match(group_order)
  topn_metric <- rlang::arg_match(topn_metric, multiple = T)
  featuregroup_style <- rlang::arg_match(featuregroup_style, multiple = T)
  featuregroup_membership <- rlang::arg_match(featuregroup_membership)
  pval_filter <- rlang::arg_match(pval_filter)

  if (!is.null(features_topn)) {
    features_topn <- max(1, features_topn)
  }
  scale <- rlang::arg_match(scale)

  if (!is.null(features) && !is.null(features_topn)) {
    message("Providing features and features_topn may not be meaningful.")
  }

  if (identical(legend_fill_args[["title"]], "..auto..")) {
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
      assay = assay,
      meta_col = meta_col,
      levels_calc = levels_calc,
      levels_plot = levels_plot)
  plot_groups <- unlist(levels_plot, use.names = FALSE)
  plot_groups <- plot_groups[!is.na(plot_groups)]

  c(wil_auc_raw,
    wil_auc) %<-% determine_features(
      SO = SO,
      assay = assay,
      features = features,
      meta_col = meta_col)

  # optional scaling
  # issue: scaling before vs. after selection for topn features
  # do here before feature selection to not disturb zscoring result e.g.
  # but what if features are selected manually with features argument?
  # for scaling of each feature it is only relevant that all groups are still included
  # determine_features does not remove groups but features only

  # Inside mutate(), dplyr looks for bare names in the data first. .env$scale explicitly reads the scale argument passed to heatmap_pseudobulk2(), even if the data later gains a column named scale.
  # The old .default = avgExpr fallback was redundant. arg_match() allows only "zscore", "1", and "none". The outer branch handles "none" by plotting avgExpr directly, so the inner if/`else` only needs to handle the two scaling modes.
  if (scale != "none") {
    wil_auc <- dplyr::mutate(
      wil_auc,
      avgExprScale = if (.env$scale == "zscore") {
        as.vector(base::scale(avgExpr))
      } else {
        as.vector(brathering::scale2(avgExpr, min = -1, max = 1))
      },
      .by = feature
    )
    values <- "avgExprScale"
  } else {
    values <- "avgExpr"
  }

  # do this here but not in heatmap_long_df to allow pre-filter for pct_in and padj
  if (!is.null(features_topn)) {
    select <-
      wil_auc |>
      dplyr::filter(pct_in >= min_pct) |>
      dplyr::filter(padj <= max_padj) |>
      dplyr::mutate(padj = -padj) |> # rank smaller adjusted p-values first
      dplyr::slice_max(order_by = avgExpr, n = 1, by = feature) |>
      dplyr::slice_max(order_by = tibble::tibble(!!!rlang::syms(as.list(topn_metric))),
                       n = features_topn,
                       by = group, with_ties = topn_ties) |>
      dplyr::filter(group %in% .env$plot_groups) |>
      dplyr::pull(feature)
    wil_auc <- wil_auc[which(wil_auc[["feature"]] %in% select),,drop = F]
  } else if (is.null(features_topn) && min_pct_force) {
    select <-
      wil_auc |>
      dplyr::filter(pct_in >= min_pct) |>
      dplyr::filter(group %in% .env$plot_groups) |>
      dplyr::pull(feature)
    rm_feat <- setdiff(as.character(unique(wil_auc$feature)), select)
    if (length(rm_feat)) {
      message("features rm due to min_pct and min_pct_force (", length(rm_feat), ")")
      #message(paste(rm_feat, collapse = ", "))
    }

    wil_auc <- wil_auc[which(wil_auc[["feature"]] %in% select),,drop = F]
  }

  # Apply the display order after filtering by the original group labels.
  wil_auc <- wil_auc |>
    dplyr::filter(group %in% .env$plot_groups) |>
    dplyr::mutate(group = factor(group, levels = .env$plot_groups))
  if (!any(!is.na(wil_auc[[values]]))) {
    stop("No expression values remain for plotting. Check levels_plot, feature filters, and scale.")
  }

  if (featuregroup_by_peak && is.null(feature_groups)) {
    calc_groups <- unlist(levels_calc, use.names = FALSE)
    calc_groups <- calc_groups[!is.na(calc_groups)]
    retained_features <- levels(droplevels(wil_auc$feature))

    markers <- wil_auc_raw |>
      dplyr::filter(
        feature %in% .env$retained_features,
        group %in% .env$calc_groups,
        logFC > 0,
        padj < 1e-4
      )
    if (featuregroup_membership == "top") {
      markers <- markers |>
        dplyr::arrange(match(group, .env$calc_groups)) |>
        dplyr::slice_max(order_by = avgExpr, n = 1L, by = feature, with_ties = FALSE)
    }
    markers <- dplyr::arrange(markers, match(feature, .env$retained_features))
    marker_groups <- factor(markers$group, levels = calc_groups)
    feature_groups <- split(markers$feature, marker_groups, drop = TRUE)
    marker_features <- unique(unlist(feature_groups, use.names = FALSE))
    if (!length(marker_features)) {
      stop("No retained features have positive logFC and padj below 1e-4 in any named group.",
           call. = FALSE)
    }
    excluded <- length(setdiff(retained_features, marker_features))
    if (excluded > 0L) {
      message("Omitting ", excluded, " features without a qualifying marker group.")
    }
    wil_auc <- wil_auc |>
      dplyr::filter(feature %in% .env$marker_features) |>
      dplyr::mutate(feature = factor(feature, levels = .env$marker_features))
  }

  featuregroup <- NULL
  feature_label_groups <- NULL
  if (!is.null(feature_groups)) {
    feat_group <- stats::setNames(utils::stack(feature_groups), c("feature", "featgroup"))
    if ("color" %in% featuregroup_style) {
      marker_peaks <- feat_group |>
        dplyr::inner_join(
          wil_auc_raw[c("feature", "group", "avgExpr")],
          by = c("feature", "featgroup" = "group"),
          relationship = "many-to-one"
        ) |>
        dplyr::slice_max(
          order_by = avgExpr, n = 1L, by = feature,
          with_ties = FALSE, na_rm = TRUE
        )
      feature_label_groups <- stats::setNames(as.character(marker_peaks$featgroup), marker_peaks$feature)
    }
    # many-to-many is fine??
    # same genes can occur in different feat_group
    # features in wil_auc exist for every group
    feature_levels <- levels(droplevels(wil_auc$feature))
    wil_auc <- dplyr::left_join(wil_auc, feat_group, by = "feature", relationship = "many-to-many")
    wil_auc$feature <- factor(wil_auc$feature, levels = feature_levels)
    featuregroup <- "featgroup"
  }



  dotsizes <- if (dotplot) "pct_in" else NULL
  color_group_labels <- featuregroup_color_on_x_axis_text &&
    !is.null(featuregroup) && "color" %in% featuregroup_style &&
    any(levels(wil_auc$group) %in% wil_auc[[featuregroup]])
  defer_repel <- featurelabels_repel && color_group_labels

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
                                  featurelabels_repel = featurelabels_repel && !defer_repel,
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
                                  pval_text_args = pval_text_args,
                                  feature_label_groups = feature_label_groups)


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

  if (color_group_labels) {
    plot <- color_heatmap_group_labels(plot, axes_flip)
  }
  if (sec_axis) {
    sec_axis_features <- levels(as.factor(plot$data$feature))
  }
  if (defer_repel) {
    # Reuse the resolved feature labels before the repeller adds its own scale.
    labels <- plot$scales$get_scales("y")$breaks
    plot$scales$scales <- plot$scales$scales[!plot$scales$find("y")]
    plot_data <- plot$data
    plot_data$feature <- as.factor(plot_data$feature)
    plot <- fcexpr:::repel_features(
      df = plot_data,
      plot = plot,
      repel_args = repel_args,
      featurelabels = labels,
      featuresitalic = featuresitalic
    )
  }
  if (sec_axis) {
    plot <- add_sec_axis(
      plot = plot,
      features = sec_axis_features,
      convert_gene_identifier_args = convert_gene_identifier_args
    )
  }

  return(list(
    plot = plot,
    data = wil_auc,
    complete_data = wil_auc_raw
  ))
}

color_heatmap_group_labels <- function(plot, axes_flip) {
  group_colors <- brathering::gg_get_scale_mapping(plot, aesthetic = "colour")
  axis_element <- if (axes_flip) "axis.text.y.left" else "axis.text.x.bottom"
  axis_text <- ggplot2::calc_element(axis_element, plot$theme)
  if (inherits(axis_text, "element_blank") ||
      !any(levels(plot$data$group) %in% names(group_colors))) {
    return(plot)
  }

  # The group scale still maps x when coord_flip() draws it vertically.
  group_scale <- plot$scales$get_scales("x")
  group_scale$labels <- function(groups) {
    colors <- unname(group_colors[groups])
    matched <- !is.na(colors)
    labels <- gsub("([[:punct:]])", "\\\\\\1", groups, perl = TRUE)
    labels[matched] <- sprintf('<span style="color:%s">%s</span>', colors[matched], labels[matched])
    labels
  }

  if (!inherits(axis_text, "element_markdown")) {
    fields <- c("family", "face", "size", "colour", "hjust", "vjust",
                "angle", "lineheight", "margin", "debug")
    text_args <- stats::setNames(lapply(fields, function(field) axis_text[[field]]), fields)
    axis_text <- do.call(ggtext::element_markdown, text_args)
  }
  plot$theme[[axis_element]] <- axis_text
  plot
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


check_meta_and_levels <- function(SO, assay, meta_col, levels_calc, levels_plot) {

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


  levels_plot_supplied <- !is.null(levels_plot)
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
  if (levels_plot_supplied && any(purrr::map2_lgl(levels_plot, levels_calc, ~ any(!.x %in% .y)))) {
    message("Some levels_plot are absent from levels_calc and will be omitted.")
  }
  levels_plot <- purrr::map2(levels_plot, levels_calc, intersect)

  selected_cells <- purrr::pmap(list(SO, meta_col, levels_calc), function(obj, column, groups) {
    metadata <- obj@meta.data
    cells <- rownames(metadata)[metadata[[column]] %in% groups]
    if (length(cells) && any(cells %in% scexpr:::get_layer_cells(obj[[assay]], "data"))) {
      return(cells)
    }
    character()
  })
  keep <- lengths(selected_cells) > 0L
  if (!any(keep)) {
    stop("No cells in the assay data layers match levels_calc.", call. = FALSE)
  }
  if (any(!keep)) {
    message("Skipping objects with no selected cells in the data layer: ",
            paste(names(SO)[!keep], collapse = ", "), ".")
  }

  SO <- purrr::map2(SO[keep], selected_cells[keep], function(obj, cells) {
    if (length(cells) < ncol(obj)) {
      obj <- subset2(obj, cells = cells)
    }
    obj
  })
  meta_col <- meta_col[keep]
  levels_calc <- levels_calc[keep]
  levels_plot <- levels_plot[keep]

  groups <- unlist(levels_calc, use.names = FALSE)
  labeled <- !is.na(groups)
  if (length(SO) > 1 && anyDuplicated(groups[labeled])) {
    prefixes <- rep(make.unique(names(SO)), lengths(levels_calc))
    labels <- rep(NA_character_, length(groups))
    labels[labeled] <- make.unique(paste0(prefixes[labeled], "___", groups[labeled]))
    labels_by_object <- utils::relist(labels, levels_calc)

    SO <- purrr::pmap(
      list(SO, meta_col, levels_calc, labels_by_object),
      function(obj, col, original, labels) {
        obj@meta.data[[col]] <- labels[match(obj@meta.data[[col]], original)]
        obj
      }
    )
    levels_plot <- purrr::pmap(
      list(levels_plot, levels_calc, labels_by_object),
      function(plot_groups, original, labels) labels[match(plot_groups, original)]
    )
    levels_calc <- labels_by_object
  }
  return(list(SO, meta_col, levels_calc, levels_plot))
}


determine_features <- function(SO,
                               assay,
                               features,
                               meta_col) {

  wil_auc_raw <- find_all_marker(obj = SO,
                                 meta_col = meta_col,
                                 assay = assay,
                                 calc_avglog2fc = F, # for speed? or change fun to avg_log2FC
                                 eps = 1e-6,
                                 features = features)

  if (is.null(features)) {
    features <- wil_auc_raw$feature
  }

  wil_auc <- wil_auc_raw |>
    dplyr::filter(feature %in% features) |>
    dplyr::mutate(feature = factor(feature, levels = unique(features)))

  return(list(wil_auc_raw, wil_auc))
}


add_sec_axis <- function(plot, features, convert_gene_identifier_args) {
  axis.df <-
    Gmisc::fastDoCall(scexpr::convert_gene_identifier,
                      args = c(list(idents = features), convert_gene_identifier_args)) |>
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
