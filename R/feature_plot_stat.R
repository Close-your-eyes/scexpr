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
#' @param combine
#' @param ncol.combine
#' @param nrow.combine
#' @param nrow.inner
#' @param ncol.inner
#' @param feature.aliases
#' @param cutoff.feature
#' @param cutoff.expression
#' @param exclusion.feature
#' @param plot.panel.grid
#' @param font.family
#' @param make.cells.unique
#' @param theme
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
feature_plot_stat <- function(SO,
                              features,
                              assay = c("RNA", "SCT"),
                              geom1 = c("jitter", "point", "dotplot"),
                              geom2 = c("violin", "boxplot", "none"),
                              jitterwidth = 0.2,
                              label.size = 4,
                              meta.col,
                              meta.col.levels = NULL,
                              plot.expr.freq = F,
                              filter.non.expr = F,
                              cells = NULL,
                              downsample = 1,
                              split.by = NULL,
                              pt.size = 1,
                              combine = T,
                              ncol.combine = NULL,
                              nrow.combine = NULL,
                              nrow.inner = NULL,
                              ncol.inner = NULL,
                              feature.aliases = NULL,
                              cutoff.feature = NULL,
                              cutoff.expression = 0,
                              exclusion.feature = NULL,
                              plot.panel.grid = F,
                              font.family = "sans",
                              make.cells.unique = F,
                              theme = ggplot2::theme_bw(),
                              ...) {

  if (missing(SO)) {stop("Seurat object list or feature vector is missing.")}
  if (missing(meta.col)) {stop("meta.col reauired.")}
  if (length(features) == 1) {combine <- F}
  if (combine && (!is.null(ncol.combine) && !is.null(nrow.combine))) {stop("Please only select one, ncol.combine or nrow.combine. Leave the other NULL.")}
  if (!is.null(ncol.inner) && !is.null(nrow.inner)) {stop("Please only select one, ncol.inner or nrow.inner. Leave the other NULL.")}
  if (!is.null(split.by)) {warning("split.by requires testing.")}

  assay <- match.arg(assay, c("RNA", "SCT"))
  geom1 <- match.arg(geom1, c("jitter", "point", "dotplot"))
  geom2 <- match.arg(geom2, c("violin", "boxplot", "none"))

  browser()
  SO <- .check.SO(SO = SO, assay = assay, split.by = split.by) # length = 1 only one SO currently
  features <- .check.features(SO = SO, features = unique(features), meta.data = F)
  if (length(meta.col) > 1) {
    stop("Please provide only one meta.col.")
  }
  if (length(SO) > 1) {
    if (!is.null(meta.col.levels)) {
      if (!is.list(meta.col.levels) || length(meta.col.levels) != length(SO)) {
        stop("meta.col.levels has to be list of same length as SO.")
      }
      if (is.null(names(meta.col.levels)) || length(intersect(names(meta.col.levels), names(SO))) < length(meta.col.levels)) {
        stop("meta.col.levels has to have names which match names of SO.")
      }
    }
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
    if (is.null(meta.col.levels)) {
      ### FIX!!
      meta.col.levels <- stats::setNames(paste0(lapply(names(SO), function(x) unique(SO[[x]]@meta.data[,meta.col]) ), "__", x), nm = names(SO))
    } else {
      meta.col.levels <- lapply(names(meta.col.levels), function(x) paste0(.check.levels(SO = SO[[x]], meta.col = meta.col, levels = meta.col.levels[[x]], append_by_missing = F), "__", x))
    }
    data[,meta.col] <- paste0(data[,meta.col], "__", data[,"SO.split"])
  } else {
    meta.col.levels <- .check.levels(SO = SO, meta.col = meta.col, levels = meta.col.levels, append_by_missing = F)
  }
  data <- data[which(data[,meta.col] %in% meta.col.levels),]

  # split.by requires testing - include in pivoting etc and geom_text

  data <- tidyr::pivot_longer(data, dplyr::all_of(features), names_to = "Feature", values_to = "expr")

  if (plot.expr.freq) {
    stat <-
      data %>%
      dplyr::group_by(Feature) %>%
      dplyr::mutate(max.feat.expr = max(expr)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("Feature", meta.col, "max.feat.expr")))) %>%
      dplyr::summarise(pct.expr = sum(expr > 0)/dplyr::n(), .groups = "drop")
  }

  if (filter.non.expr) {
    data <- dplyr::filter(data, expr > 0)
  }

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!sym(meta.col), y = expr))
  if (geom1 == "jitter") {
    plot <- plot + ggplot2::geom_jitter(width = jitterwidth, size = pt.size, color = "tomato2")
  }
  if (geom1 == "point") {
    plot <- plot + ggplot2:: geom_point(size = pt.size, color = "tomato2")
  }
  if (geom1 == "dotplot") {
    plot <- plot + ggplot2::geom_dotplot(binaxis = "y", binwidth = (max(dd$expr) - min(dd$expr))/40, stackdir = "center", fill = "tomato2", stackratio = 0.7)
  }
  if (geom2 == "violin") {
    plot <- plot + ggplot2::geom_violin(alpha = 0.5)
  }
  if (geom2 == "boxplot") {
    plot <- plot + ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA)
  }
  if (plot.expr.freq) {
    plot <- plot + ggplot2::geom_text(data = stat, ggplot2::aes(label = round(pct.expr, 2), y = max.feat.expr + 0.4), size = label.size, family = font.family)
    #expand_limits
  }

  plot <- plot + theme
  if (!plot.panel.grid) {plot <- plot + ggplot2::theme(panel.grid = ggplot2::element_blank())}
  plot <- plot + ggplot2::theme(...)



  plot <- plot + ggplot2::facet_wrap(ggplot2::vars(Feature), scales = "free_y")

  return(plot)
}
