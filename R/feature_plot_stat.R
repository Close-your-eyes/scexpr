#' Title
#'
#' @param SO
#' @param features
#' @param assay
#' @param geom1
#' @param geom2
#' @param jitterwidth
#' @param label.size
#' @param meta.col
#' @param plot.expr.freq to plot labels of expression frequencies
#' @param plot.non.expr do plot dots of non-expressers (0 UMIs)
#' @param cells
#' @param downsample
#' @param split.by
#' @param pt.size
#' @param feature.aliases
#' @param cutoff.feature
#' @param cutoff.expression
#' @param exclusion.feature
#' @param font.family
#' @param theme
#' @param ...
#' @param expr.freq.hjust
#' @param legend.title
#' @param col.pal
#' @param col.pal.dir
#' @param expr.freq.decimals decimal precision of frequency labels
#' @param expr.freq.pct plot expression frequency in percent?
#'
#' @return
#' @export
#'
#' @examples
feature_plot_stat <- function(SO,
                              features,
                              meta.col,
                              assay = c("RNA", "SCT"),
                              geom1 = c("jitter", "point"),
                              geom2 = c("boxplot", "violin", "none"),
                              jitterwidth = 0.2,
                              label.size = 4,
                              plot.expr.freq = F,
                              expr.freq.decimals = 2,
                              expr.freq.pct = F,
                              expr.freq.hjust = 0.3,
                              plot.non.expr = T,
                              plot.strip = T,
                              cells = NULL,
                              downsample = 1,
                              legend.title = "SO.split",
                              split.by = NULL,
                              pt.size = 0.5,
                              col.pal = "custom",
                              col.pal.dir = 1,
                              feature.aliases = NULL,
                              cutoff.feature = NULL,
                              cutoff.expression = 0,
                              exclusion.feature = NULL,
                              font.family = "sans",
                              theme = ggplot2::theme_bw(),
                              facetting_args = list(scales = "free_y"),
                              ...) {

  # add option to plot facet labels in italics

  if (missing(SO)) {
    stop("Seurat object list or feature vector is missing.")
  }
  if (missing(meta.col)) {
    stop("meta.col reauired.")
  }
  if (!is.null(split.by)) {
    warning("split.by requires testing.")
  }

  assay <- match.arg(assay, c("RNA", "SCT"))
  geom1 <- match.arg(geom1, c("jitter", "point"))
  geom2 <- match.arg(geom2, c("boxplot", "violin", "none"))

  SO <- .check.SO(SO = SO, assay = assay, split.by = split.by) # length = 1 only one SO currently
  features <- .check.features(SO = SO, features = unique(features), meta.data = T)
  if (length(meta.col) > 1) {
    stop("Please provide only one meta.col.")
  }

  meta.col <- .check.features(SO = SO, features = meta.col, rownames = F)
  cells <- .check.and.get.cells(SO = SO,
                                assay = assay,
                                cells = cells,
                                cutoff.feature = cutoff.feature,
                                cutoff.expression = cutoff.expression,
                                exclusion.feature = exclusion.feature,
                                downsample = downsample,
                                return.included.cells.only = T)

  # get data with SO being no list does not work
  data <- .get.data(SO, #stats::setNames(list(SO), "1")
                    feature = features,
                    assay = assay,
                    slot = "data",
                    cells = cells,
                    split.by = split.by,
                    reduction = NULL,
                    meta.col = meta.col)


  # split.by requires testing - include in pivoting etc and geom_text
  data <- tidyr::pivot_longer(data, dplyr::all_of(features), names_to = "Feature", values_to = "expr")

  if (plot.expr.freq) {
    stat <-
      data %>%
      dplyr::group_by(Feature) %>%
      dplyr::mutate(max.feat.expr = max(expr)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("Feature", meta.col, "SO.split", "max.feat.expr")))) %>%
      dplyr::summarise(pct.expr = sum(expr > 0)/dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(pct.expr.adjust.pct = dplyr::case_when(pct.expr == 0 ~ "0 %",
                                                           pct.expr > 0 & pct.expr < 0.01 ~ "> 1 %",
                                                           pct.expr >= 0.01 ~ paste0(round(pct.expr*100, expr.freq.decimals), " %"))) %>%
      dplyr::mutate(pct.expr.adjust = dplyr::case_when(pct.expr == 0 ~ "0",
                                                       pct.expr > 0 & pct.expr < 0.01 ~ "> 0.01",
                                                       pct.expr >= 0.01 ~ as.character(round(pct.expr, expr.freq.decimals))))
    names(stat)[which(names(stat) == "SO.split")] <- legend.title
  }
  names(data)[which(names(data) == "SO.split")] <- legend.title


  if (!plot.non.expr) {
    data <- dplyr::filter(data, expr > 0)
  }

  data <- as.data.frame(data)
  my_geom2 <- switch(geom2,
                     "violin" = ggplot2::geom_violin,
                     "boxplot" = ggplot2::geom_boxplot,
                     "none" = NULL)

  if (length(SO) > 1) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(meta.col), y = expr, color = !!rlang::sym(legend.title)))
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(meta.col), y = expr, color = !!rlang::sym(meta.col)))
  }

  if (!is.null(my_geom2)) {
    plot <- plot + suppressWarnings(my_geom2(outlier.shape = NA, position = ggplot2::position_dodge(width = 0.75)))
  }
  ##ggforce::geom_sina()
  if (geom1 == "jitter" && jitterwidth > 0) {
    plot <- plot + ggplot2::geom_point(size = pt.size, position = ggplot2::position_jitterdodge(jitter.width = jitterwidth, dodge.width = 0.75))
  } else if (geom1 == "point") {
    plot <- plot + ggplot2::geom_point(size = pt.size, position = ggplot2::position_dodge(width = 0.75))
  }

  #plot <- plot + ggplot2::geom_dotplot(binaxis = "y", binwidth = (max(data$expr) - min(data$expr))/40, stackdir = "center", fill = "black", stackratio = 0.7)

  if (plot.expr.freq) {
    if (expr.freq.pct) {
      label_col <- "pct.expr.adjust.pct"
    } else {
      label_col <- "pct.expr.adjust"
    }
    plot <- plot + ggplot2::geom_text(data = stat,
                                      ggplot2::aes(label = !!rlang::sym(label_col),
                                                   y = max.feat.expr + expr.freq.hjust),
                                      position = ggplot2::position_dodge(width = 0.75),
                                      size = label.size,
                                      family = font.family, show.legend = F)
    #expand_limits
  }

  if (length(col.pal) == 1 && !col.pal %in% grDevices::colors()) {
    col.pal <- col_pal(name = col.pal, direction = col.pal.dir)
  }

  plot <- plot + theme
  if (length(SO) == 1) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }

  plot <- plot + ggplot2::theme(...)
  plot <- plot + ggplot2::scale_color_manual(values = col.pal)
  plot <- plot + Gmisc::fastDoCall(ggh4x::facet_wrap2, args = c(list(facets = ggplot2::vars(Feature)), facetting_args))

  if (!plot.strip) {
    plot <- plot + theme(strip.background = ggplot2::element_blank(), strip.text.x = ggplot2::element_blank())
  }

  return(plot)
}

