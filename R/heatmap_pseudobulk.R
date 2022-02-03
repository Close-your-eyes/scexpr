#' Title
#'
#' @param SO
#' @param assay
#' @param meta.col
#' @param levels.calc which levels in meta.col to include in calculation; all levels if NULL; also defines order on x-axis
#' @param levels.plot which levels in meta.col (and subset of levels.calc) to include in platting; set to levels.calc if NULL; also defines order on x-axis
#' @param features
#' @param normalization
#' @param topn.features
#' @param topn.metric
#' @param min.pct
#' @param max.padj
#' @param title
#' @param title.font.size
#' @param font.family
#' @param y.font.size
#' @param tile.borders
#' @param dotplot
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
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
heatmap_pseudobulk <- function(SO,
                               assay = c("RNA", "SCT"),
                               meta.col = NULL,
                               levels.calc = NULL,
                               levels.plot = NULL,
                               features = NULL,
                               normalization = c(-1,1), # scale
                               topn.features = 10,
                               topn.metric = c("logFC", "auc", "padj"),
                               min.pct = 0.1,
                               max.padj = 0.05,
                               title = NULL,
                               title.font.size = 14,
                               font.family = "Courier",

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
                               ...) {

  # ... arguments to ggrepel, like nudge_y and scexpr::convert_gene_identifier

  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github("immunogenomics/presto")
  }
  assay <- match.arg(assay, c("RNA", "SCT"))


  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  if (class(normalization) == "numeric") {
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
  }

  legend.direction <- match.arg(legend.direction, c("horizontal", "vertical"))
  if (legend.direction == "horizontal") {
    temp <- legend.barwidth
    legend.barwidth <- legend.barheight
    legend.barheight <- temp
  }

  levels.calc <- .check.levels(SO = SO, meta.col = meta.col, levels = levels.calc, append_by_missing = F)

  if (is.null(levels.plot)) {
    levels.plot <- levels.calc
  } else {
    if (any(!levels.plot %in% unique(SO@meta.data[,meta.col,drop=T]))) {
      print(paste0("levels.plot not found in meta.col: ", paste(levels.plot[which(!levels.plot %in% unique(SO@meta.data[,meta.col,drop=T]))], collapse = ", ")))
    }
    levels.plot <- as.character(unique(levels.plot[which(levels.plot %in% unique(SO@meta.data[,meta.col,drop=T]))]))
    if (any(!levels.plot %in% levels.calc)) {
      print("levels.plot has more levels than levels.calc. levels.plot is reduced to levels.calc.")
      levels.plot <- levels.plot[which(levels.plot %in% levels.calc)]
    }
  }

  # subset levels for calculation
  if (length(levels.calc) < length(unique(SO@meta.data[,meta.col,drop=T]))) {
    SO <- subset(SO, cells = rownames(SO@meta.data[,meta.col,drop=F][which(SO@meta.data[,meta.col] %in% levels.calc),,drop=F]))
  }

  if (missing(features)) {
    ## presto gives deviating results with respect to avgExpr and logFC (maybe due to approximation which makes calculation faster)
    ## nevertheless, in order to select marker genes by logFC or other statistics the presto output is useful
    wil_auc <- presto::wilcoxauc(SO, meta.col, seurat_assay = assay)
    topn.metric <- match.arg(topn.metric, c("logFC", "auc", "padj"))

    features <-
      wil_auc %>%
      dplyr::group_by(feature) %>%
      dplyr::slice_max(order_by = avgExpr, n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(pct_in >= min.pct) %>%
      dplyr::filter(padj <= max.padj) %>%
      dplyr::mutate(group = factor(group, levels = levels.calc)) %>%
      dplyr::group_by(group) %>%
      dplyr::slice_max(order_by = !!rlang::sym(topn.metric), n = topn.features) %>%
      dplyr::arrange(group, avgExpr) %>%
      dplyr::pull(feature)
  } else {
    features <- .check.features(SO = SO, features = unique(features), meta.data = F)
    #features <- sapply(features, function(x) grep(x, rownames(Seurat::GetAssayData(SO, assay = assay, slot = "data")), ignore.case = T, value = T))
  }

  raw_tab <- Seurat::AverageExpression(SO, assays = assay, group.by = meta.col, slot = "data", verbose = F)[[1]][features,,drop=F]
  if (class(normalization) == "character" && normalization == "scale") {
    tab <- row_scale(raw_tab, add_attr = F)
  } else if (class(normalization) == "numeric") {
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
  scale.mid <- as.numeric(format(round(scale.min + ((scale.max - scale.min) / 2), 1), nsmall = 1))

  if (class(normalization) == "numeric") {
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
      heatmap.plot <- heatmap.plot + ggplot2::geom_tile(colour = "black")
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
    heatmap.plot <- cowplot::plot_grid(axis, heatmap.plot, align = "h", axis = "tb", nrow = 1, rel_widths = c(feature.labels.axis.width,1))
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


  if (missing(features) & is.null(feature.labels) & nlevels(htp$Feature) > 200 & plot.feature.breaks) {
    print("Large number of rows without selecting feature.labels. Consider setting plot.feature.breaks to FALSE or selecting feature.labels to plot.")
  }

  if (!plot.feature.breaks) {
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }

  return(list(plot = heatmap.plot, data = htp))
}

.check.levels <- function(SO, meta.col, levels = NULL, append_by_missing = F) {
  if (is.null(levels) || is.na(levels)) {
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
