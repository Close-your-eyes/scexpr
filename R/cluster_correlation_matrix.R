#' Title
#'
#' @param SO
#' @param assay
#' @param meta.cols
#' @param levels
#' @param complement.levels
#' @param features
#' @param avg.expr
#' @param method
#' @param corr.in.percent
#' @param round.corr
#'
#' @return
#' @export
#'
#' @examples
cluster_correlation_matrix <- function (SO,
                                        assay = c("RNA", "SCT"),
                                        meta.cols,
                                        levels = NULL, # NA possible
                                        complement.levels = F,
                                        features = NULL, # all, pca
                                        avg.expr,
                                        method = c("pearson", "kendall", "spearman"),
                                        corr.in.percent = FALSE,
                                        round.corr = 2) {

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

  if (is.null(features)) {
    print("Using all intersecting features. Provide features = 'all' to avoid this message.")
    features <- Reduce(intersect, lapply(SO, function(x) rownames(x)))
  } else if (features == "all") {
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
