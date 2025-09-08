#' Plot a heatmap of average gene transcription of transcriptomes by
#' different groups (clusters)
#'
#' @param SO Seurat object
#' @param meta_col which column from meta.data of SO to use as x-axis;
#' if NULL current Idents(SO) are used
#' @param levels_calc which levels in meta_col to include in calculation;
#' all levels if NULL; the order provided defines order on x-axis; thisn level
#' selection will affect scaling calculations
#' @param levels_plot which levels in meta_col to include in platting;
#' if NULL this equals to levels_calc; must be subset of levels_calc;
#' the order provided defines order on x-axis; this will not affect scaling
#' calculation but only select levels for plotting;
#' ordering will not work perfectly when more than 1 SO is provided and factor levels
#' in meta_cols are not unique
#' @param assay which assay to obtain expression values from
#' @param features choose which features to plot; if NULL features_topn becomes
#' relevant
#' @param features_topn if no features are selected, this will select how many
#' features to plot per level in meta_col; selection is done based on the metric
#' selected in topn_metric; respective features with greatest difference between
#' meta_col levels are selected (best DE features so to say);
#' if NULL and features is NULL, all features are plotted
#' @param topn_metric which differential expression metric to apply for feature
#' selection; only relevant if features_topn is not NULL
#' @param min_pct for marker feature selection per only: min fraction of expressing
#' cells for a gene to become marker gene for a group
#' @param max_padj for marker feature selection per only: max adjusted
#' p-value for a gene to become marker gene
#' @param dotplot plot dots instead of tiles with percent of expressing cells
#' as dot size
#' @param fill color palette vector for fill of tiles or dots. when auto,
#' RColorBrewer::RdBu is used
#' @param color color of stroke (border) around tiles or dots; "NA" means no
#' stroke is plotted; NA has
#' to be put in quotation mark ("NA"), such that geom_point accepts it.
#' other choices may be black, white or any other color code; when "auto",
#' grey70 is used by default when is.null(dotsizes) and a the number of
#' features is below 100.
#' @param scale how to scale values: not, zscore or from -1 to 1
#' @param features_topn only plot the top n features (ordered by value) per
#' group
#' @param color_linewidth linewidth of borders (stroke) around tiles or dots
#' @param legendbreaks a single number, a vector of explicit breaks, or "auto"
#' for ggplot default or "minmidmax" for three breaks at minimum, middle and
#' maximum of value range
#' @param legendlabels labels for breaks, e.g. c("min", "mid", "max")
#' @param colorsteps NULL to have normal colorbar, auto for default colorsteps,
#' a single number or a vector of explicit steps; may not work with any number
#' when colorsteps_nice is TRUE
#' @param colorsteps_nice heuristic for pretty steps
#' @param axes_flip do flip them?
#' @param group_seplines plot lines that separate features belonging to
#' different groups
#' @param seplines_args arguments to geom_hline for seplines
#' @param legend_fill_args arguments to ggplot2::guide_colorsteps or
#' ggplot2::guide_colorbar, depend upon the color scale
#' @param legend_size_args arguments to ggplot2::guide_legend to modify the
#' size legend; e.g. use override.aes = list(size = c(1,3,5)) to adjust dot size
#' in legend in contrast to dotsize_range, one number for each dot in legend
#' needed
#' @param featurelabels subset of feature labels to plot; NULL to plot all,
#' "" to plot none; can be a named vector with names being aliases to use for
#' plotting e.g. c("CD20" = "MS4A1", "CD3", "KLRG1") to only alter MS4A1;
#' auto to omit labels by default when more than 100 features are there
#' @param featurelabels_repel do repel feature axis labels?
#' @param featuresitalic shorthand to plot feature labels in italic
#' @param theme_args arguments to ggplot2::theme
#' @param repel_args fine tuning for featurelabel repelling
#' @param theme ggplot2 theme
#' @param dotsize_range range for dot size
#' @param topn_ties break ties for features_topn? if TRUE, more then
#' features_topn may be plotted
#' @param feature_order if and how to order features?
#' @param group_order if and how to order groups?
#' @param sec_axis plot a secondary feature axis with full gene names
#' @param convert_gene_identifier_args arguments to
#' scexpr::convert_gene_identifier
#'
#' @importFrom zeallot %<-%
#'
#' @return
#' @export
#'
#' @examples
heatmap_pseudobulk2 <- function(SO,
                                meta_col = NULL,
                                levels_calc = NULL,
                                levels_plot = NULL,
                                assay = "RNA",
                                min_pct = 0.1,
                                max_padj = 0.05,
                                features = NULL,
                                feature_order = c("custom", "hclust", "none"),
                                group_order = c("hclust", "none", "custom"),
                                topn_metric = c("padj", "auc", "logFC"),
                                topn_ties = F,
                                dotplot = F,
                                dotsize_range = c(2,7),
                                fill = "auto",
                                color = "auto",
                                scale = c("zscore", "1", "none"),
                                features_topn = NULL,
                                featurelabels = "auto",
                                featurelabels_repel = F,
                                featuresitalic = T,
                                color_linewidth = 0.2,
                                legendbreaks = "auto",
                                legendlabels = "auto",
                                colorsteps = "auto",
                                colorsteps_nice = T,
                                axes_flip = F,
                                group_seplines = F,
                                seplines_args = list(),
                                theme = ggplot2::theme_classic(),
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
                                  title = "transcription\nfrequency [%]",
                                  title.hjust = 0.5,
                                  label.position = "bottom",
                                  order = 2,
                                  ncol = NULL,
                                  nrow = NULL,
                                  override.aes = list(color = "black")
                                  #size = c(2,4, 10))
                                ),
                                theme_args = list(
                                  axis.title = ggplot2::element_blank(),
                                  panel.grid = ggplot2::element_blank(),
                                  axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                                  legend.position = "right",
                                  legend.direction = "vertical"
                                ),
                                repel_args = list(featurelabels_width = 0.2,
                                                  featurelabels_nudhe_x = -1),
                                sec_axis = F,
                                convert_gene_identifier_args = list(ident_in = "SYMBOL", ident_out = "GENENAME")) {

  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github("immunogenomics/presto")
  }

  feature_order <- rlang::arg_match(feature_order)
  group_order <- rlang::arg_match(group_order)
  topn_metric <- rlang::arg_match(topn_metric, multiple = T)
  if (!is.null(features_topn)) {
    features_topn <- max(1, features_topn)
  }
  scale <- rlang::arg_match(scale)

  if (!is.null(features) && !is.null(features_topn)) {
    message("Providing features and features_topn may not be meaningful.")
  }

  if ("title" %in% names(legend_fill_args) && legend_fill_args[["title"]] == "..auto..") {
    if (scale == "zscore") {
      legend_fill_args[["title"]] <- "transcription\nlevel [z-score]"
    } else if (scale == "none") {
      legend_fill_args[["title"]] <- "transcription\nlevel [log (UMI)]"
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

  SO <- check.SO(SO = SO, assay = assay)
  assay <- Seurat::DefaultAssay(SO[[1]])

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
      feature_order)
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
      dplyr::slice_max(order_by = tibble::tibble(!!!rlang::syms(as.list(topn_metric))), n = features_topn, by = group, with_ties = topn_ties) |>
      # filter here but not in determine_features to not bias feature selection by levels_plot
      dplyr::filter(group %in% unlist(levels_plot)) |>
      dplyr::pull(feature)
    wil_auc <- wil_auc[which(wil_auc[["feature"]] %in% select),,drop = F]
  }
  # filter here after feature selection
  wil_auc <- dplyr::filter(wil_auc, group %in% unlist(levels_plot))

  plot <- fcexpr::heatmap_long_df(df = wil_auc,
                                  groups = "group",
                                  features = "feature",
                                  values = values,
                                  dotsizes = if (dotplot) "pct_in" else NULL,
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
                                  feature_order = feature_order,
                                  group_order = group_order)

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
    if (suppressWarnings(!any(is.na(as.numeric(levels))))) {
      levels <- as.character(sort(as.numeric(unique(levels))))
    }
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


determine_features <- function(SO, assay, features, levels_plot, meta_col, feature_order) {



  presto_feat <- unique(unlist(purrr::map(SO, ~names(which(Matrix::rowSums(get_layer(obj = .x, layer = "data", assay = assay)) > 0)))))
  if (!is.null(features)) {
    features <- check.features(SO = SO, features = unique(features), meta.data = F)
    if (any(!features %in% presto_feat)) {
      message("No expressers found for: ", paste(features[which(!features %in% presto_feat)], collapse = ","), ". Will not be plotted.")
    }
    presto_feat <- intersect(presto_feat, features)
  }

  # by default, presto gives deviating results with respect to avgExpr and logFC
  # use expm1 to match results Seurats procedure, see example below (comparison of presto and Seurat)
  wil_auc_raw <- presto::wilcoxauc(X = Gmisc::fastDoCall(cbind, purrr::map(SO, ~expm1(get_layer(obj = .x, layer = "data", assay = assay, features = presto_feat)))),
                                   y = unlist(purrr::pmap(list(x = SO, y = meta_col), function(x,y) x@meta.data[,y,drop=T]), use.names = F))

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

