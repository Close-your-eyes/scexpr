heatmap_pseudobulk <- function(SO,
                               assay = c("RNA", "SCT"),
                               meta.col,
                               meta.col.levels,
                               meta.col.levels.select,
                               features = NULL,
                               feature.labels,
                               topn.features = 10e5,
                               topn.metric = c("logFC", "auc", "padj"),
                               min.pct = 0.1,
                               max.padj = 0.05,
                               title = NULL,
                               title.font.size = 14,

                               y.font.size = 10,
                               legend.position = "right",
                               legend.direction = "vertical",
                               tile.borders = T,
                               dotplot = F,
                               plot.feature.breaks = T,
                               plot.secondary.axis.with.full.genenames = F,

                               legend.title.text.size = 10,
                               legend.barheight = 8,
                               legend.barwidth = 1,
                               legend.text.size = 10,
                               legend.title = "z-score of avg expr",
                               legend.title.position = "top",
                               feature.labels.nudge_x = -0.1,
                               feature.labels.axis.width = 0.2,
                               ...) {

  if (!requireNamespace("presto", quietly = T)){
    devtools::install_github("immunogenomics/presto")
  }
  assay <- match.arg(assay, c("RNA", "SCT"))
  meta.col <- match.arg(meta.col, names(SO@meta.data))
  legend.direction <- match.arg(legend.direction, c("horizontal", "vertical"))
  if (legend.direction == "horizontal") {
    temp <- legend.barwidth
    legend.barwidth <- legend.barheight
    legend.barheight <- temp
  }

  if (missing(meta.col.levels)) {
    meta.col.levels <- as.character(unique(SO@meta.data[,meta.col,drop=T]))
    if (suppressWarnings(!NA %in% as.numeric(meta.col.levels))) {
      meta.col.levels <- as.character(sort(as.numeric(unique(meta.col.levels))))
    }
  } else {
    meta.col.levels <- as.character(unique(meta.col.levels[which(meta.col.levels %in% unique(SO@meta.data[,meta.col,drop=T]))]))
  }

  # subset by meta.col.levels
  if (length(meta.col.levels) < length(unique(SO@meta.data[,meta.col,drop=T]))) {
    SO <- subset(SO, cells = rownames(SO@meta.data[,meta.col,drop=F][which(SO@meta.data[,meta.col] %in% meta.col.levels),,drop=F]))
  }

  if (missing(meta.col.levels.select)) {
    meta.col.levels.select <- meta.col.levels
  }

  if (missing(features)) {
    ## presto gives deviating results with respect to avgExpr and logFC (maybe due to approximation which makes calculation faster)
    ## nevertheless, in order to select marker genes by logFC or other statistics the presto output may be useful
    wil_auc <- presto::wilcoxauc(SO, meta.col, seurat_assay = assay)
    topn.metric <- match.arg(topn.metric, c("logFC", "auc", "padj"))

    features <-
      wil_auc %>%
      dplyr::group_by(feature) %>%
      dplyr::slice_max(order_by = avgExpr, n = 1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(pct_in >= min.pct) %>%
      dplyr::filter(padj <= max.padj) %>%
      dplyr::mutate(group = factor(group, levels = meta.col.levels)) %>%
      dplyr::group_by(group) %>%
      dplyr::slice_max(order_by = !!sym(topn.metric), n = topn.features) %>%
      dplyr::arrange(group, avgExpr) %>%
      dplyr::pull(feature)
  } else {
    features <- sapply(features, function(x) {grep(x, rownames(GetAssayData(SO, assay = assay, slot = "data")), ignore.case = T, value = T)})
  }

  htp <-
    as.data.frame(Seurat::AverageExpression(SO, assays = assay, group.by = meta.col, slot = "data", verbose = F)[[1]]) %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::filter(Feature %in% features) %>%
    tidyr::pivot_longer(cols = -Feature, names_to = "cluster", values_to = "avgExpr") %>%
    dplyr::group_by(Feature) %>%
    dplyr::mutate(scaledAvgExpr = scale(avgExpr)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cluster %in% meta.col.levels.select) %>%
    dplyr::left_join(wil_auc[,c(which(names(wil_auc) %in% c("feature", "group", "pct_in")))], by = c("Feature" = "feature", "cluster" = "group")) %>%
    dplyr::mutate(cluster = factor(cluster, levels = meta.col.levels)) %>%
    dplyr::mutate(Feature = factor(Feature, levels = features))

  scale.max <- as.numeric(format(floor_any(max(htp$scaledAvgExpr), 0.1), nsmall = 1))
  scale.min <- as.numeric(format(ceiling_any(min(htp$scaledAvgExpr), 0.1), nsmall = 1))
  scale.mid <- as.numeric(format(round(0, 1), nsmall = 1))

  heatmap.plot <-
    ggplot(htp, aes(x = cluster, y = Feature, fill = scaledAvgExpr)) +
    scale_fill_gradientn(values = scales::rescale(c(scale.min, scale.mid, scale.max)), colours = colour.scale.function(scale.name = "RdBu", n.brewer.colors = 9), breaks = c(scale.min, scale.mid, scale.max)) +
    xlab("Cluster") +
    theme_classic() +
    ggtitle(title) +
    theme(title = element_text(size = title.font.size, family = "Courier"), axis.title = element_blank(), axis.text.x = element_text(family = "Courier"), axis.text.y = element_text(size = y.font.size, face = "italic", family = "Courier"), legend.position = legend.position, legend.direction = legend.direction) +
    guides(fill = guide_colourbar(barwidth = legend.barwidth, barheight = legend.barheight, label.theme = element_text(size = legend.text.size, family = "Courier"), title.theme = element_text(size = legend.title.text.size, family = "Courier"), title.position = legend.title.position, title = legend.title, draw.ulim = FALSE, draw.llim = FALSE))

  if (dotplot) {
    # shape legend not working yet
    heatmap.plot <- heatmap.plot + geom_point(aes(size = pct_in), shape = 21) + guides(shape = guide_legend(override.aes = list(label.theme = element_text(size = legend.text.size, family = "Courier"), title.theme = element_text(size = legend.title.text.size, family = "Courier"))))
  } else {
    if (tile.borders) {
      heatmap.plot <- heatmap.plot + geom_tile(colour = "black")
    } else {
      heatmap.plot <- heatmap.plot + geom_tile()
    }
  }

  if (plot.feature.breaks & !missing(feature.labels)) {
    axis.df <- data.frame(y = 1:length(levels(htp$Feature)), Feature = levels(htp$Feature))
    axis <- ggplot(axis.df, aes(x = 0, y = y, label = Feature)) +
      ggrepel::geom_text_repel(fontface = "italic", family = "Courier", data = axis.df[which(axis.df$Feature %in% feature.labels),], aes(label = Feature), nudge_x = feature.labels.nudge_x, direction = "y", ...) +
      scale_x_continuous(limits = c(-0.1, 0), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      scale_y_continuous(limits = c(0, length(levels(htp$Feature)) + 0.5), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      theme_void()
    heatmap.plot <- heatmap.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
    heatmap.plot <- cowplot::plot_grid(axis, heatmap.plot, align = "h", axis = "tb", nrow = 1, rel_widths = c(feature.labels.axis.width,1))
  }

  if (plot.secondary.axis.with.full.genenames) {
    out <- convert_gene_identifier(unique(htp$Feature), species = "Hs")
    axis.df <- data.frame(y = 1:length(levels(htp$Feature)), Feature = levels(htp$Feature)) %>% dplyr::left_join(out, by = c("Feature" = "SYMBOL"))
    axis <- ggplot(axis.df, aes(x = 0, y = y, label = GENENAME)) +
      geom_text(hjust = 0) +
      scale_x_continuous(limits = c(0, 2), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      scale_y_continuous(limits = c(0.4, length(levels(htp$Feature)) + 0.4), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      theme(panel.background = element_blank(), plot.margin = margin(0, 0, 0, 0, "pt"), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
    heatmap.plot <- cowplot::plot_grid(heatmap.plot + theme(plot.margin = margin(0, 0, 0, 0, "pt")), axis, align = "h", axis = "tb", nrow = 1, rel_widths = c(0.6,0.4))
  }


  if (missing(features) & missing(feature.labels) & nlevels(htp$Feature) > 200 & plot.feature.breaks) {
    print("Large number of rows without selecting feature.labels. Consider setting plot.feature.breaks to FALSE or selecting feature.labels to plot.")
  }

  if (!plot.feature.breaks) {
    heatmap.plot <- heatmap.plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  return(list(plot = heatmap.plot, data = htp))
}
