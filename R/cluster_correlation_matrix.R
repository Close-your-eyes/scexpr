cluster_correlation_matrix <- function (SO,
                                        assay = c("RNA", "SCT"),
                                        meta.cols,
                                        axis.label.orders,
                                        features.to.use = "intersecting.var.features",
                                        avg.expr,
                                        legend.position = "none",
                                        method = c("pearson", "kendall", "spearman"),
                                        corr.in.percent = FALSE,
                                        decimal.places = 2) {

  if (missing(meta.cols)) {
    stop("Please provide meta.cols for the first and second SO.")
  } else {
    if (length(meta.cols) != 2) {
      stop("meta.cols has to be of length 2.")
    }
  }

  SO <- .check.SO(SO, assay = assay, length = 2)
  assay <- match.arg(assay, c("RNA", "SCT"))
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  meta.cols[1] <- .check.features(SO[[1]], features = meta.cols[1], rownames = F)
  meta.cols[2] <- .check.features(SO[[2]], features = meta.cols[2], rownames = F)



  if (!missing(avg.expr) && !is.list(avg.expr)) {
    stop("avg.expr has to be a list.")
  }

  if (!missing(axis.label.orders) && !is.list(axis.label.orders)) {
    stop("axis.label.orders has to be a list.")
  }

  if (missing(axis.label.orders)) {
    axis.label.orders <- list(sort(unique(SO.list[[1]]@meta.data[,resolutions[1]])), sort(unique(SO.list[[2]]@meta.data[,resolutions[2]])))
  }


  if (length(features.to.use) > 1) {
    features.to.use <- unique(features.to.use)
    features.to.use <- features.to.use[which(features.to.use %in% Reduce(intersect, lapply(SO.list, function(x) rownames(x))))]
    if (length(features.to.use) < 3) {
      stop("Less than 3 genes left after filtering to use for correlation analysis: ", paste(features.to.use, collapse = ", "), " Please provide more.")
    }
  } else {
    if (features.to.use == "intersecting.var.features") {
      #VariableFeatures(SO) equals rownames(SO@reductions[["pca"]]@feature.loadings)
      features.to.use <- Reduce(intersect, lapply(SO.list, function(x) VariableFeatures(x)))
    } else if (features.to.use == "intersecting.features") {
      features.to.use <- Reduce(intersect, lapply(SO.list, function(x) rownames(GetAssayData(x, assay = assay))))
    } else if (features.to.use == "var.features") {
      features.to.use <- unique(unlist(lapply(SO.list, function(x) VariableFeatures(x))))
    }
  }

  if (missing(avg.expr)) {
    avg.expr <- lapply(SO.list, function (x) {Seurat::AverageExpression(x, assays = assay)[[assay]]})
    #Seurat changes column name if all cells are within one group only
    avg.expr <- lapply(seq_along(avg.expr), function(x) {
      if (ncol(avg.expr[[x]]) == 1) {
        colnames(avg.expr[[x]]) <- axis.label.orders[[x]]
      }
      return(avg.expr[[x]][features.to.use,which(colnames(avg.expr[[x]]) %in% axis.label.orders[[x]]), drop = F])
    })
  } else {
    avg.expr <- lapply(avg.expr, function (x) {x[features.to.use,,drop=F]})
  }

  if (missing(axis.label.orders)) {
    axis.label.orders <- list(colnames(avg.expr[[1]]), colnames(avg.expr[[2]]))
  }

  if (correlation == "pearson") {
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
    '    lmss <<- lms
    avg.exprs <<- avg.expr'

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

    '    test <- lms[[1]][[1]]
    res <- stack(test[["residuals"]]) %>% dplyr::rename("res" = values) %>%
    fit <- stack(test[["fitted.values"]]) %>% dplyr::rename("fit" = values)
    res.corr <- stack(test[["residuals"]]/test[["fitted.values"]]) %>% dplyr::rename("res.corr" = values)
    df <-
      res %>%
      dplyr::left_join(fit) %>%
      dplyr::left_join(res.corr)'

  } else {
    cm <- cor(avg.expr[[1]], avg.expr[[2]], method = correlation)
  }

  '  rr <- data.frame(x = scale(avg.expr[[2]][,1]), y = scale(avg.expr[[1]][,1]))
  ggplot(rr , aes(x,y)) +
    geom_point() +
    geom_smooth(method = "lm")
  '

  cm.melt <- reshape2::melt(cm)
  cm.melt$Var1 <- factor(as.character(cm.melt$Var1), levels = axis.label.orders[[1]])
  cm.melt$Var2 <- factor(as.character(cm.melt$Var2), levels = axis.label.orders[[2]])

  cm.plot <- ggplot(cm.melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(colour = "black") +
    scale_fill_gradient2(high = "#BC3F2B", low = "#4C6CA6", mid = "#edf9ff", midpoint = mean(cm.melt$value)) +
    xlab(names(SO.list)[1]) +
    ylab(names(SO.list)[2]) +
    theme_classic() +
    theme(text = element_text(size=24), plot.title = element_text(size=16), axis.text.y = element_text(face = "bold"), axis.text.x = element_text(face = "bold"), legend.position = legend.position)

  if (corr.in.percent) {
    cm.plot <- cm.plot + geom_text(size = 8, aes(label = paste0(round(value*100, 0), " %")))
  } else {
    cm.plot <- cm.plot + geom_text(size = 8, aes(label = format(round(value, decimal.places), nsmall = decimal.places)))
  }

  # dot plot of genes
  avg.expr <- lapply(avg.expr, function(x) reshape2::melt(as.matrix(x)))
  names(avg.expr[[1]]) <- c("Gene.symbol", "cluster.x", "avg.expr.x")
  names(avg.expr[[2]]) <- c("Gene.symbol", "cluster.y", "avg.expr.y")

  dot.plot.matrix.df <- left_join(avg.expr[[1]], avg.expr[[2]], by = "Gene.symbol")

  dot.plot.matrix <- ggplot(dot.plot.matrix.df, aes(x = scales::squish(avg.expr.x, c(0.001,100000)), y = scales::squish(avg.expr.y, c(0.001,100000)))) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(rows = vars(cluster.x), cols = vars(cluster.y), scales = "free")

  return(list(cm.plot = cm.plot, dot.plot.matrix = dot.plot.matrix, dot.plot.matrix.df = dot.plot.matrix.df, cm.melt = cm.melt, avg.expr = avg.expr))
}
