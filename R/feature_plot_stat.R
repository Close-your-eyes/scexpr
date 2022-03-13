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
#' @param meta.col.levels
#' @param plot.expr.freq
#' @param filter.non.expr
#' @param cells
#' @param downsample
#' @param split.by
#' @param pt.size
#' @param feature.aliases
#' @param cutoff.feature
#' @param cutoff.expression
#' @param exclusion.feature
#' @param plot.panel.grid
#' @param font.family
#' @param make.cells.unique
#' @param theme
#' @param ...
#' @param expr.freq.hjust
#' @param legend.title
#' @param col.pal
#' @param col.pal.rev
#' @param nrow
#' @param ncol
#'
#' @return
#' @export
#'
#' @examples
feature_plot_stat <- function(SO,
                              features,
                              assay = c("RNA", "SCT"),
                              geom1 = c("jitter", "point"),
                              geom2 = c("violin", "boxplot", "none"),
                              jitterwidth = 0.2,
                              label.size = 4,
                              meta.col,
                              meta.col.levels = NULL,
                              plot.expr.freq = F,
                              expr.freq.hjust = 0.3,
                              filter.non.expr = F,
                              cells = NULL,
                              downsample = 1,
                              legend.title = "SO.split",
                              split.by = NULL,
                              pt.size = 0.5,
                              col.pal = "custom",
                              col.pal.rev = F,
                              nrow = NULL,
                              ncol = NULL,
                              feature.aliases = NULL,
                              cutoff.feature = NULL,
                              cutoff.expression = 0,
                              exclusion.feature = NULL,
                              plot.panel.grid = F,
                              font.family = "sans",
                              make.cells.unique = F,
                              theme = ggplot2::theme_bw(),
                              ...) {

  if (missing(SO)) {
    stop("Seurat object list or feature vector is missing.")
  }
  if (missing(meta.col)) {
    stop("meta.col reauired.")
  }
  if (!is.null(ncol) && !is.null(nrow)) {
    stop("Please only select one, ncol or nrow Leave the other NULL.")
  }
  if (!is.null(split.by)) {
    warning("split.by requires testing.")
  }

  assay <- match.arg(assay, c("RNA", "SCT"))
  geom1 <- match.arg(geom1, c("jitter", "point"))
  geom2 <- match.arg(geom2, c("violin", "boxplot", "none"))


  SO <- .check.SO(SO = SO, assay = assay, split.by = split.by) # length = 1 only one SO currently
  features <- .check.features(SO = SO, features = unique(features), meta.data = F)
  if (length(meta.col) > 1) {
    stop("Please provide only one meta.col.")
  }


  ## procedure to allow subsetting by SOs
  if (length(SO) > 1) {
    if (!is.null(meta.col.levels)) {
      if (!is.list(meta.col.levels) || length(meta.col.levels) != length(SO)) {
        stop("meta.col.levels has to be list of same length as SO.")
      }
      if (is.null(names(meta.col.levels)) || length(intersect(names(meta.col.levels), names(SO))) < length(meta.col.levels)) {
        stop("meta.col.levels has to have names which match names of SO.")
      }
      meta.col.levels <- sapply(names(meta.col.levels), function(x) paste0(.check.levels(SO = SO[[x]], meta.col = meta.col, levels = meta.col.levels[[x]], append_by_missing = F), "__", x), simplify = F)
    } else {
      meta.col.levels <- lapply(SO, function(x) unique(x@meta.data[,meta.col]))
      meta.col.levels <- sapply(names(meta.col.levels), function(x) paste0(meta.col.levels[[x]], "__", x), simplify = F)
    }
  } else {
    meta.col.levels <- .check.levels(SO = SO[[1]], meta.col = meta.col, levels = meta.col.levels, append_by_missing = F)
  }
  meta.col <- .check.features(SO = SO, features = meta.col, rownames = F)
  cells <- .check.and.get.cells(SO = SO,
                                assay = assay,
                                cells = cells,
                                make.cells.unique = make.cells.unique,
                                cutoff.feature = cutoff.feature,
                                cutoff.expression = cutoff.expression,
                                exclusion.feature = exclusion.feature,
                                downsample = downsample,
                                make.cells.unique.warning = 1,
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


  if (length(SO) > 1) {
    data[,meta.col] <- paste0(data[,meta.col], "__", data[,"SO.split"])
  }

  data <- data[which(data[,meta.col] %in% unlist(meta.col.levels)),]
  data[,meta.col] <- stringr::str_replace(data[,meta.col], paste0("__",data[,"SO.split"]), "")


  # split.by requires testing - include in pivoting etc and geom_text
  data <- tidyr::pivot_longer(data, dplyr::all_of(features), names_to = "Feature", values_to = "expr")

  if (plot.expr.freq) {
    stat <-
      data %>%
      dplyr::group_by(Feature) %>%
      dplyr::mutate(max.feat.expr = max(expr)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("Feature", meta.col, "SO.split", "max.feat.expr")))) %>%
      dplyr::summarise(pct.expr = sum(expr > 0)/dplyr::n(), .groups = "drop")
    names(stat)[which(names(stat) == "SO.split")] <- legend.title
  }
  names(data)[which(names(data) == "SO.split")] <- legend.title


  if (filter.non.expr) {
    data <- dplyr::filter(data, expr > 0)
  }

  data <- as.data.frame(data)
  data[,meta.col] <- factor(data[,meta.col], levels = meta.col.levels)
  my_geom2 <- switch(geom2,
                     "violin" = ggplot2::geom_violin,
                     "boxplot" = ggplot2::geom_boxplot,
                     "none" = NULL)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(meta.col), y = expr, color = !!rlang::sym(legend.title)))
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
    plot <- plot + ggplot2::geom_text(data = stat, ggplot2::aes(label = round(pct.expr, 2), y = max.feat.expr + expr.freq.hjust), position = ggplot2::position_dodge(width = 0.75), size = label.size, family = font.family, show.legend = F)
    #expand_limits
  }

  if (length(col.pal) == 1 && !col.pal %in% grDevices::colors()) {
    col.pal <- col_pal(name = col.pal, reverse = col.pal.rev)
  }

  plot <- plot + theme
  if (length(SO) == 1) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }
  if (!plot.panel.grid) {plot <- plot + ggplot2::theme(panel.grid = ggplot2::element_blank())}
  plot <- plot + ggplot2::theme(...)
  plot <- plot + ggplot2::scale_color_manual(values = col.pal)
  plot <- plot + ggplot2::facet_wrap(ggplot2::vars(Feature), scales = "free_y", nrow = nrow, ncol = ncol)

  return(plot)
}

