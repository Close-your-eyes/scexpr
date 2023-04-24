#' Title
#'
#' @param SO Seurat object
#' @param meta.col which column from meta.data of SO to use as x-axis; if NULL current Idents(SO) are used
#' @param levels.calc which levels in meta.col to include in calculation;
#' all levels if NULL; the order provided defines order on x-axis; level selection
#' will affect scaling calculations
#' @param levels.plot which levels in meta.col to include in platting;
#' if NULL this equals to levels.calc; must be subset of levels.calc;
#' the order provided defines order on x-axis; this will not affect scaling
#' calculation but only select levels for plotting
#' @param assay which assay to obtain expression values from
#' @param features optionally choose which features to plot (supervised)
#' @param normalization how to scale expression values; may be "scale" to
#' use base::scale and transform average expression value to a standardized
#' normal distribution or may be a numeric vector of length 2 for scaling
#' each feature from min (first value) to max (second value); in the latter case
#' c(-1,1) is most meaningful
#' @param topn.features if no features are selected, this will select how many
#' features to plot per level in meta.col; selection is done based on the metric
#' selected in topn.metric; respective features with greatest difference between
#' meta.col levels are selected (best DE features so to say)
#' @param topn.metric which differential expression metric to apply for feature
#' selection
#' @param min.pct filter features across levels in meta.col for a minimal fraction
#' of expressing cells; 0.1 stands for min 10 % expressing cells
#' @param max.padj filter features across levels in meta.col for a significance
#' level of differential expression (one level vs. all others)
#' @param title which title to add to plot; if NULL no title is plotted
#' @param title.font.size font size of title
#' @param font.family which font type (family) to use for plotting,
#' e.g. mono or sans
#' @param y.font.size font size of features names on y-axis
#' @param tile.borders plot black tile border (TRUE) or not (FALSE);
#' better set to FALSE for may features and optionally set to TRUE for few
#' features only
#' @param dotplot plot heatmap as a dotplot which will also reveal the fraction
#' of expressing cells by dot sizes
#' @param plot.feature.breaks
#' @param plot.sec.axis
#' @param legend.position
#' @param legend.direction
#' @param legend.barheight
#' @param legend.barwidth
#' @param legend.text.size
#' @param legend.title.text.size
#' @param legend.title.fill
#' @param legend.title.size
#' @param legend.title.position
#' @param legend.labels
#' @param feature.labels
#' @param feature.labels.nudge_x
#' @param feature.labels.axis.width
#' @param ...
#' @param order_features when features are provided should they be ordered automatically to generate a pretty heatmap;
#' not relevant when is.null(features)
#' @param plot_hlines_between_groups plot horizontal lines to separate groups of genes that are
#' differentially expressed by a group on the x-axis; if hlines is NULL the breaks are
#' determined automatically
#' @param break.ties strategy to break ties of topn.metric; if 'second_metric' then a secondary metric from
#' c("padj", "logFC", "auc") is used (if topn.metric is 'padj', then 'logFC' is secondary which is the most relevant case
#' and auc is tertiary; if topn.metric is 'logFC' then 'padj' is second followed by auc); if break.ties is logical
#' (TRUE or FALSE) it is passed to 'with_ties' to dplyr::slice_max and only topn.metric is considered
#' @param hlines provide manual breaks of where to plot horizontal lines;
#' even though the y-axis is discrete, a numeric vector works here:
#' hlines = c(1.5, 3.5) will plot lines between the first and second, and, third and fourth
#' entry on the y-axis, respectively
#' @param hlines_args arguments to geom_hline like linewidth or linetype
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
heatmap_pseudobulk <- function(SO,
                               meta.col = NULL,
                               levels.calc = NULL,
                               levels.plot = NULL,
                               assay = c("RNA", "SCT"),
                               features = NULL,
                               order_features = F,
                               normalization = "scale", # scale
                               topn.features = 10,
                               break.ties = c("second_metric", T, F),
                               topn.metric = c("padj", "logFC", "auc"),
                               min.pct = 0.1,
                               max.padj = 0.05,
                               title = NULL,
                               title.font.size = 14,
                               font.family = "sans",

                               y.font.size = 10,
                               tile.borders = T,
                               dotplot = F,
                               plot.feature.breaks = T,
                               plot.sec.axis = F,

                               legend.position = "right",
                               legend.direction = "vertical",
                               legend.barheight = 8,
                               legend.barwidth = 1,
                               legend.text.size = 10,
                               legend.title.text.size = 10,
                               legend.title.fill = NULL,
                               legend.title.size = NULL,
                               legend.title.position = "top",
                               legend.labels = c("min", "", "max"),
                               feature.labels = NULL,
                               feature.labels.nudge_x = -0.1,
                               feature.labels.axis.width = 0.2,
                               plot_hlines_between_groups = F,
                               hlines = NULL,
                               hlines_args = list(),
                               ...) {

  # ... arguments to ggrepel, like nudge_y and scexpr::convert_gene_identifier

  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github("immunogenomics/presto")
  }

  assay <- match.arg(assay, c("RNA", "SCT"))
  topn.metric <- match.arg(topn.metric, c("padj", "logFC", "auc"))
  break.ties <- match.arg(break.ties, c("second_metric", T, F))

  if (topn.features < 1) {
    message("topn.features has to be minimum 1.")
    topn.features <- 1
  }

  if (plot_hlines_between_groups && !is.null(features) && !order_features && is.null(hlines)) {
    message("horizontal lines between group (hlines) can only be inferred automatically if order_features is set to TRUE.
            Will set plot_hlines_between_groups to FALSE now.
            Alternatively provide them with the hlines argument manually")
    plot_hlines_between_groups <- F
  }

  if (!plot_hlines_between_groups && !is.null(hlines)) {
    message("hlines provided but plot_hlines_between_groups is FALSE.")
  }

  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  if (methods::is(normalization, "numeric")) {
    if (length(normalization) != 2) {
      stop("Please set normalization to 'scale' or a numeric vector of length 2 (e.g. c(-1,1)).")
    }
    if (length(legend.labels) > 3) {
      stop("legend.labels can have max length of 3 only")
    }
    if (length(legend.labels) == 2) {
      legend.labels <- c(legend.labels[1], "", legend.labels[2])
    }
  } else {
    normalization <- match.arg(normalization, choices = c("scale"))
  }

  if (is.null(meta.col)) {
    SO@meta.data$idents <- Seurat::Idents(SO)
    meta.col <- "idents"
  } else {
    meta.col <- match.arg(meta.col, names(SO@meta.data))
    if (is.numeric(SO@meta.data[,meta.col])) {
      stop("meta.col is numeric. Please make it factor or character.")
    }
  }

  legend.direction <- match.arg(legend.direction, c("horizontal", "vertical"))
  if (legend.direction == "horizontal") {
    temp <- legend.barwidth
    legend.barwidth <- legend.barheight
    legend.barheight <- temp
  }

  if (!is.null(levels.plot) && !is.null(levels.calc)) {
    if (any(!levels.plot %in% levels.calc)) {
      print("levels.plot has more levels than levels.calc. levels.plot is reduced to levels.calc.")
      levels.plot <- levels.plot[which(levels.plot %in% levels.calc)]
    }
  }

  # this will replace NULL with all available levels
  levels.calc <- .check.levels(SO = SO, meta.col = meta.col, levels = levels.calc, append_by_missing = F)
  levels.plot <- .check.levels(SO = SO, meta.col = meta.col, levels = levels.plot, append_by_missing = F)

  # subset levels for calculation
  if (length(levels.calc) < length(unique(SO@meta.data[,meta.col,drop=T]))) {
    SO <- subset(SO, cells = rownames(SO@meta.data[,meta.col,drop=F][which(SO@meta.data[,meta.col] %in% levels.calc),,drop=F]))
  }

  # this problem was discovered by chance (":" are replaced by Seurat::AverageExpression by "_"); maybe other symbols will also cause problems.
  if (any(grepl(":", unique(SO@meta.data[,meta.col])))) {
    message("Found colon (:) in factor levels of meta.col. Will replace those with underscores (_).")
    SO@meta.data[,meta.col] <- gsub(":", "_", SO@meta.data[,meta.col])
    levels.plot <- gsub(":", "_", levels.plot)
    levels.calc <- gsub(":", "_", levels.calc)
  }

  ## presto gives deviating results with respect to avgExpr and logFC (maybe due to approximation which makes calculation faster)
  ## nevertheless, in order to select marker genes by logFC or other statistics fast perfromance of presto is useful
  ## filter all features with zero expression across the data set
  wil_auc <-
    presto::wilcoxauc(SO, meta.col, seurat_assay = assay) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(sum(pct_in) > 0) %>% # filter non-expressed features
    dplyr::ungroup()

  # prep features variable here to have one common pipeline below
  if (!is.null(features)) {
    features <- .check.features(SO = SO, features = unique(features), meta.data = F)
    if (any(!features %in% wil_auc$feature)) {
      message("No expressers found for: ", paste(features[which(!features %in% wil_auc$feature)], collapse = ","), ". Will not be plotted.")
    }
  } else {
    order_features <- T
    features <- wil_auc$feature
  }

  if (order_features) {

'    if (topn.metric == "padj") {
      slice_fun <- dplyr::slice_min
    } else {
      slice_fun <- dplyr::slice_max
    }'

    # select multiple columns for ordering, in order to break ties; especially if padj is used as primary element to order
    if (is.character(break.ties) && break.ties == "second_metric") {
      all_metrics <- c("padj", "logFC", "auc")
      second_metric <- all_metrics[which(!all_metrics %in% topn.metric)][1]
      third_metric <- all_metrics[which(!all_metrics %in% topn.metric)][2]
    }

    features <-
      wil_auc %>%
      dplyr::filter(feature %in% features) %>%
      dplyr::filter(pct_in >= min.pct) %>%
      dplyr::filter(padj <= max.padj) %>%
      dplyr::mutate(padj = -padj) %>% # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc"))
      dplyr::group_by(feature) %>%
      dplyr::slice_max(order_by = avgExpr, n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(group = factor(group, levels = levels.plot)) %>%
      dplyr::group_by(group)
    if (is.character(break.ties) && break.ties == "second_metric") {
      features <-
        features %>%
        dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = topn.features)
    } else {
      features <-
        features %>%
        dplyr::slice_max(order_by = !!rlang::sym(topn.metric), n = topn.features, with_ties = break.ties)
    }
    features <-
      features %>%
      dplyr::ungroup() %>%
      dplyr::arrange(group, !!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric))

'
    ## alternative feature selection:
    features <-
      wil_auc %>%
      dplyr::filter(feature %in% features) %>%
      dplyr::filter(pct_in >= min.pct) %>%
      dplyr::filter(padj <= max.padj) %>%
      dplyr::mutate(padj = -padj) %>% # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc"))
      dplyr::mutate(group = factor(group, levels = levels.plot)) %>%
      dplyr::group_by(group)
    if (is.character(break.ties) && break.ties == "second_metric") {
      features <-
        features %>%
        dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = topn.features)
    } else {
      features <-
        features %>%
        dplyr::slice_max(order_by = !!rlang::sym(topn.metric), n = topn.features, with_ties = break.ties)
    }

    features <-
      features %>%
      dplyr::ungroup() %>%
      dplyr::arrange(group, !!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric))


    ## check for number of features per group
    ## and maybe add additional ones, but how?
    features_check <-
      features %>%
      dplyr::group_by(group) %>%
      dplyr::slice_max(order_by = avgExpr, n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(padj <= max.padj) %>%
      dplyr::mutate(group = factor(group, levels = levels.plot)) %>%
      dplyr::group_by(group)
'



    ## continue
    if (is.null(hlines)) {
      hlines <- cumsum(rle(as.character(features$group))[["lengths"]]) + 0.5
      hlines <- hlines[-length(hlines)]
    }
    features <-
      features %>%
      dplyr::pull(feature)
  }

  raw_tab <- Seurat::AverageExpression(SO, assays = assay, group.by = meta.col, slot = "data", verbose = F)[[1]][features,,drop=F]

  if (methods::is(normalization, "character") && normalization == "scale") {
    tab <- row_scale(raw_tab, add_attr = F)
  } else if (methods::is(normalization, "numeric")) {
    tab <- scale_min_max(raw_tab, min = min(normalization), max = max(normalization), margin = 1)
  }

  htp <-
    as.data.frame(tab) %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::filter(Feature %in% features) %>%
    tidyr::pivot_longer(cols = -Feature, names_to = "cluster", values_to = "norm_avgexpr") %>%
    dplyr::filter(cluster %in% levels.plot) %>%
    dplyr::left_join(wil_auc[,c(which(names(wil_auc) %in% c("feature", "group", "pct_in")))], by = c("Feature" = "feature", "cluster" = "group")) %>%
    dplyr::mutate(cluster = factor(cluster, levels = levels.plot)) %>%
    dplyr::mutate(Feature = factor(Feature, levels = features))

  scale.max <- as.numeric(format(floor_any(max(htp$norm_avgexpr), 0.1), nsmall = 1))
  scale.min <- as.numeric(format(ceiling_any(min(htp$norm_avgexpr), 0.1), nsmall = 1))
  #scale.mid <- as.numeric(format(round(scale.min + ((scale.max - scale.min) / 2), 1), nsmall = 1))
  scale.mid <- 0

  if (methods::is(normalization, "numeric")) {
    labels <- c("min", "int", "max")
  } else {
    labels <- c(scale.min, scale.mid, scale.max)
  }

  heatmap.plot <-
    ggplot2::ggplot(htp, ggplot2::aes(x = cluster, y = Feature, fill = norm_avgexpr)) +
    ggplot2::scale_fill_gradientn(values = scales::rescale(c(scale.min, scale.mid, scale.max)), colours = col_pal(name = "RdBu", nbrew = 9, reverse = T), breaks = c(scale.min, scale.mid, scale.max), labels = labels) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(title) +
    ggplot2::theme(title = ggplot2::element_text(size = title.font.size, family = font.family), axis.title = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(family = font.family), axis.text.y = ggplot2::element_text(size = y.font.size, face = "italic", family = font.family), legend.position = legend.position, legend.direction = legend.direction)


  if (dotplot) {
    heatmap.plot <- heatmap.plot + ggplot2::geom_point(ggplot2::aes(size = pct_in), shape = 21)
  } else {
    if (tile.borders) {
      heatmap.plot <- heatmap.plot + ggplot2::geom_tile(color = "black")
    } else {
      heatmap.plot <- heatmap.plot + ggplot2::geom_tile()
    }
  }
  heatmap.plot <-
    heatmap.plot +
    ggplot2::guides(fill = ggplot2::guide_colourbar(label.theme = ggplot2::element_text(size = legend.text.size, family = font.family), title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family), title.position = legend.title.position, title = legend.title.fill, barwidth = legend.barwidth, barheight = legend.barheight,  draw.ulim = FALSE, draw.llim = FALSE),
                    size = ggplot2::guide_legend(label.theme = ggplot2::element_text(size = legend.text.size, family = font.family), title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family), title.position = legend.title.position, title = legend.title.size))


  if (plot.feature.breaks & !is.null(feature.labels)) {
    axis.df <- data.frame(y = 1:length(levels(htp$Feature)), Feature = levels(htp$Feature))
    axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = Feature)) +
      ggrepel::geom_text_repel(fontface = "italic", family = font.family, data = axis.df[which(axis.df$Feature %in% feature.labels),], ggplot2::aes(label = Feature), nudge_x = feature.labels.nudge_x, direction = "y", ...) +
      ggplot2::scale_x_continuous(limits = c(-0.1, 0), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0, length(levels(htp$Feature)) + 0.5), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::theme_void()
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
    heatmap.plot <- cowplot::plot_grid(axis, heatmap.plot, align = "h", axis = "tb", nrow = 1, rel_widths = c(feature.labels.axis.width,1)) ## check how to replace with patchwork
  }

  if (plot.sec.axis) {
    out <- convert_gene_identifier(idents = levels(htp$Feature), ident_in = "SYMBOL", ident_out = "GENENAME", ...)
    #out <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = levels(htp$Feature), keytype = "SYMBOL", column = i, multiVals = "first")
    axis.df <- data.frame(y = 1:length(levels(htp$Feature)), Feature = levels(htp$Feature)) %>% dplyr::left_join(out, by = c("Feature" = "SYMBOL"))
    axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = GENENAME)) +
      ggplot2::geom_text(hjust = 0, family = font.family) +
      ggplot2::scale_x_continuous(limits = c(0, 2), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0.4, length(levels(htp$Feature)) + 0.4), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::theme(panel.background = ggplot2::element_blank(), plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"), axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    heatmap.plot <- cowplot::plot_grid(heatmap.plot + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")), axis, align = "h", axis = "tb", nrow = 1, rel_widths = c(0.6,0.4))
  }

  if (plot_hlines_between_groups) {
    heatmap.plot <-
      heatmap.plot +
      Gmisc::fastDoCall(ggplot2::geom_hline, args = c(list(yintercept = hlines), hlines_args))
  } else {
    hlines = NULL
  }


  if (missing(features) & is.null(feature.labels) & nlevels(htp$Feature) > 200 & plot.feature.breaks) {
    print("Large number of rows without selecting feature.labels. Consider setting plot.feature.breaks to FALSE or selecting feature.labels to plot.")
  }

  if (!plot.feature.breaks) {
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }

  return(list(plot = heatmap.plot, data = htp, complete_data = wil_auc, hlines = hlines))
}

.check.levels <- function(SO, meta.col, levels = NULL, append_by_missing = F) {
  if (any(is.null(levels)) || any(is.na(levels))) {
    levels <- as.character(unique(SO@meta.data[,meta.col,drop=T]))
    if (suppressWarnings(!any(is.na(as.numeric(levels))))) {
      levels <- as.character(sort(as.numeric(unique(levels))))
    }
  } else {
    if (any(!levels %in% unique(SO@meta.data[,meta.col,drop=T]))) {
      print(paste0("levels not found in meta.col: ", paste(levels[which(!levels %in% unique(SO@meta.data[,meta.col,drop=T]))], collapse = ", ")))
    }
    levels <- as.character(unique(levels[which(levels %in% unique(SO@meta.data[,meta.col,drop=T]))]))
  }
  if (append_by_missing) {
    levels <- c(levels, as.character(unique(SO@meta.data[,meta.col,drop=T])[which(!as.character(unique(SO@meta.data[,meta.col,drop=T])) %in% levels)]))
  }
  return(levels)
}
