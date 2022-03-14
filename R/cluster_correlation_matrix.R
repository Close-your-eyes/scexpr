#' Find groups (clusters) of cells (single cell transcriptomes) with highest gene expression correlation in two Seurat objects
#'
#' Sometimes one wants to avoid combining data sets into a single seurat object as this often requires to employ an integration procedure or a
#' batch correction and respective uncertainties. At some point then one may be interested in which clusters correspond to
#' each other in separate seurat object. This function quickly returns a graphic and data frames to judge best matching clusters.
#' Seurat::AverageExpression is used to derive average gene expressions per cluster. Selected gene (features) are then used
#' for correlation calculation.
#'
#' @param SO named list of exactly 2 Seurat objects
#' @param assay which assay to obtain expression values from
#' @param meta.cols character vector of length 2, indicating the column names of meta.data in SO[[1]] and SO[[2]] to use for comparison,
#' e.g. the clustering columns in both SOs
#' @param levels list of length 2 indicating the factor levels in meta.cols to include; this also defines axis orders;
#' can be incomplete e.g. if order matters only for some levels, if complement.levels = T the missing levels are added randomly
#' @param complement.levels add all missing levels in random order
#' @param features features to use for correlation calculation; will always be reduced to intersecting features between SOs;
#' if all, all intersecting features in assay are considered; if pca, intersecting rownames of feature.loadings in pca-reduction;
#' or a vector of features
#' @param avg.expr avg.expr from a previous return list of cluster_correlation_matrix; e.g. provide that to  only change axis
#' orders or similar but avoid repeated calculation
#' @param method which correlation metric to calculate, passed to stats::cor
#' @param corr.in.percent display correlation in percent (TRUE, T) or as a fraction (FALSE, F)
#' @param round.corr decimal places to round correlation values to
#'
#' @return list with (i) ggplot object of correlation matrix plot, (ii) the data frame to that plot and (iii) calculated average expressions from Seurat::AverageExpression
#' @export
#'
#' @examples
cluster_correlation_matrix <- function(SO,
                                       assay = c("RNA", "SCT"),
                                       meta.cols,
                                       levels = NULL, # NA possible
                                       complement.levels = F,
                                       features = c("all", "pca"),
                                       avg.expr,
                                       method = c("pearson", "kendall", "spearman"),
                                       corr.in.percent = FALSE,
                                       round.corr = 2) {

  if (!requireNamespace("reshape2", quietly = T)) {
    utils::install.packages("reshape2")
  }

  if (length(features) == 1) {
    match.arg(features, choices = c("all", "pca"))
  }

  if (missing(meta.cols)) {
    stop("Please provide meta.cols for the first and second SO.")
  } else {
    if (length(meta.cols) != 2) {
      stop("meta.cols has to be of length 2.")
    }
  }
  if (!missing(avg.expr)) {
    if (!is.list(avg.expr)) {
      stop("avg.expr has to be a list of matrices.")
    }
    if (length(avg.expr) != 2) {
      stop("avg.expr has to be of length 2.")
    }
  }

  if (!is.null(levels)) {
    if (!is.list(levels)) {
      stop("levels has to be a list.")
    }
    if (length(levels) != 2) {
      stop("levels has to be of length 2.")
    }
    levels <- sapply(seq_along(SO), function(x) {
      .check.levels(SO[[x]], meta.cols[x], levels = levels[[x]], append_by_missing = complement.levels)
    })

  } else {
    levels <- sapply(seq_along(SO), function(x) {
      .check.levels(SO[[x]], meta.cols[x], levels = levels, append_by_missing = F)
    })
  }

  SO <- .check.SO(SO, assay = assay, length = 2)
  assay <- match.arg(assay, c("RNA", "SCT"))
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  meta.cols[1] <- .check.features(SO[[1]], features = meta.cols[1], rownames = F)
  meta.cols[2] <- .check.features(SO[[2]], features = meta.cols[2], rownames = F)

  if (features == "all") {
    features <- Reduce(intersect, lapply(SO, function(x) rownames(x)))
  } else if (features == "pca") {
    features <- Reduce(intersect, lapply(SO, function(x) rownames(x@reductions[["pca"]]@feature.loadings)))
  } else {
    features <- .check.features(SO, features, meta.data = F)
  }

  if (missing(avg.expr)) {
    avg.expr <- lapply(seq_along(SO), function (x) Seurat::AverageExpression(SO[[x]], assays = assay, features = features, group.by = meta.cols[[x]], slot = "data", verbose = F)[[assay]])
  } else {
    avg.expr <- lapply(avg.expr, function (x) x[features,,drop=F])
  }

  cm <- stats::cor(avg.expr[[1]][,levels[[1]],drop=F], avg.expr[[2]][,levels[[2]],drop=F], method = method)

  '  if (correlation == "pearson") {
    # conduct linear model in case of pearson
    # https://sebastianraschka.com/faq/docs/pearson-r-vs-linear-regr.html

    lms <- lapply(colnames(avg.expr[[1]]), function(x) {
      lms <- lapply(colnames(avg.expr[[2]]), function(y) {
        lm(y~x, data.frame(x = scale(avg.expr[[2]][,y]), y = scale(avg.expr[[1]][,x])))
      })
      names(lms) <- colnames(avg.expr[[2]])
      return(lms)
    })
    names(lms) <- colnames(avg.expr[[1]])


    cm <- sapply(lms, function(x) {
      z <- sapply(x, function(y) {
        y[["coefficients"]][["x"]]
      })
      ## just to avoid dropping
      if (length(z) == 1) {
        z <- c(z, 0)
      }
      return(z)
    })
    cm <- cm[which(!rownames(cm) == ""),,drop = F]
    cm <- t(cm)


  } else {

  }
'

  cm.melt <- reshape2::melt(cm)
  cm.melt$Var1 <- factor(as.character(cm.melt$Var1), levels = levels[[1]])
  cm.melt$Var2 <- factor(as.character(cm.melt$Var2), levels = levels[[2]])

  pp <- c(col_pal("RdBu", n = 11, reverse = T)[2:5], rep(col_pal("RdBu", n = 11, reverse = T)[6], 6), col_pal("RdBu", n = 11, reverse = T)[7:10])

  cm.plot <- ggplot2::ggplot(cm.melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(colour = "black") +
    #ggplot2::scale_fill_gradient2(high = "#BC3F2B", low = "#4C6CA6", mid = "#edf9ff", midpoint = median(cm.melt$value)) +
    ggplot2::scale_fill_gradientn(colors = pp) +
    ggplot2::labs(x = paste0(names(SO)[1], " ", meta.cols[1]), y = paste0(names(SO)[2], " ", meta.cols[2])) +
    ggplot2::theme_classic()

  if (corr.in.percent) {
    cm.plot <- cm.plot + ggplot2::geom_text(size = cm.plot[["theme"]][["text"]][["size"]] *(5/14), ggplot2::aes(label = paste0(round(value*100, 0), " %")))
  } else {
    cm.plot <- cm.plot + ggplot2::geom_text(size = cm.plot[["theme"]][["text"]][["size"]] *(5/14), ggplot2::aes(label = format(round(value, round.corr), nsmall = round.corr)))
  }

  # dot plot of genes
  '  avg.expr <- lapply(avg.expr, function(x) reshape2::melt(as.matrix(x)))
  names(avg.expr[[1]]) <- c("Gene.symbol", "cluster.x", "avg.expr.x")
  names(avg.expr[[2]]) <- c("Gene.symbol", "cluster.y", "avg.expr.y")

  dot.plot.matrix.df <- left_join(avg.expr[[1]], avg.expr[[2]], by = "Gene.symbol")

  dot.plot.matrix <- ggplot(dot.plot.matrix.df, aes(x = scales::squish(avg.expr.x, c(0.001,100000)), y = scales::squish(avg.expr.y, c(0.001,100000)))) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(rows = vars(cluster.x), cols = vars(cluster.y), scales = "free")'

  return(list(cm.plot = cm.plot, cm.melt = cm.melt, avg.expr = avg.expr))
}
