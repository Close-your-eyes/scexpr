#' Title
#'
#' @param SO
#' @param assay
#' @param volcano.data
#' @param negative.group.cells
#' @param positive.group.cells
#' @param negative.group.name
#' @param positive.group.name
#' @param x.axis.symmetric
#' @param x.axis.extension
#' @param y.axis.pseudo.log
#' @param pseudo.log.sigma
#' @param inf.fc.shift
#' @param pt.size
#' @param pt.alpha
#' @param font.size
#' @param pval.tick
#' @param min.pct
#' @param max.iter
#' @param max.overlaps
#' @param label.neg.pos.sep
#' @param label.col
#' @param label.face
#' @param font.family
#' @param label.size
#' @param labels.topn
#' @param label.features
#' @param topn.metric
#' @param nudge.x
#' @param nudge.y
#' @param p.plot
#' @param p.adjust
#' @param p.cut
#' @param p.signif
#' @param fc.cut
#' @param features.exclude
#' @param meta.cols
#' @param save.path.interactive
#' @param gsea.param
#' @param interactive.only
#' @param use.limma limma for DE gene detection; intended to use subsequent MetaVolcanoR
#' @param ... arguments to scexpr::feature_plot like  col.pal.d = setNames(c("grey90", scexpr::col_pal()[2:3]), c("other", "name1", "name2")), order.discrete = "^other",plot.title = F,
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
volcano_plot <- function(SO,
                         assay = "RNA",
                         volcano.data = NULL,
                         negative.group.cells,
                         positive.group.cells,
                         negative.group.name = "negative.group",
                         positive.group.name = "positive.group",
                         meta.cols = NULL,
                         x.axis.symmetric = T,
                         x.axis.extension = 0,
                         y.axis.pseudo.log = F,
                         pseudo.log.sigma = 1,
                         inf.fc.shift = 2,
                         pt.size = 1,
                         pt.alpha = 0.8,
                         font.size = 14,
                         pval.tick = 0.01,
                         min.pct = 0.1,
                         max.iter = 10000,
                         max.overlaps = 50,
                         label.neg.pos.sep = T,
                         label.col = "black",
                         label.face = "bold",
                         font.family = "Courier",
                         label.size = 4,
                         labels.topn = 30,
                         label.features = NULL,
                         topn.metric = "p.value",
                         nudge.x = 0,
                         nudge.y = 0,
                         p.plot = "adj.p.val",
                         p.adjust = "bonferroni",
                         p.cut = NA,
                         p.signif = 0.001,
                         fc.cut = NA,
                         features.exclude = NULL,
                         save.path.interactive = NULL,
                         gsea.param = NULL,
                         interactive.only = F,
                         use.limma = F,
                         ...) {

  if (missing(negative.group.cells) || missing(positive.group.cells)) {
    stop("positive.group.cells and negative.group.cells are required.")
  }
  if (!is.null(gsea.param)) {
    if (!is.list(gsea.param)) {
      stop("gsea.param has to be a list.")
    }
    if (length(gsea.param) != 2) {
      stop("gsea.param has to be of length 2, index 1 being log2.fc limits and index 2 being p-value limits (p.plot).")
    }
    if (any(lengths(gsea.param) != 2)) {
      stop("each index of gsea.param should be a numeric vector of length 2 indicating the limits for negative and positive DEGs (log2.fc and p val (p.plot)).")
    }
  }

  assay <- match.arg(assay, c("RNA", "SCT"))
  p.plot <- match.arg(p.plot, c("adj.p.val", "p.val"))
  p.adjust <- match.arg(p.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL)

  dots <- list(...)

  # label.features
  # add option to label all features indicated with p.cut and fc.cut

  ## check cells
  ngc <- do.call(.check.and.get.cells, args = c(list(SO = SO, assay = assay, cells = negative.group.cells, return.included.cells.only = T), dots[which(names(dots) %in% names(formals(.check.and.get.cells)))]))
  pgc <- do.call(.check.and.get.cells, args = c(list(SO = SO, assay = assay, cells = positive.group.cells, return.included.cells.only = T), dots[which(names(dots) %in% names(formals(.check.and.get.cells)))]))
  if (length(intersect(ngc, pgc)) > 0) {
    warning("Equal cell names in negative.group.cells and positive.group.cells found.")
  }

  # make meta data column in SOs to identify ngc and pgc
  colname <- paste0("cellident_", format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d%H%M%S"))
  SO <- lapply(SO, function(x) {
    x@meta.data[,colname] <- ifelse(!rownames(x@meta.data) %in% c(pgc, ngc),
                                    "other",
                                    ifelse(rownames(x@meta.data) %in% pgc,
                                           positive.group.name,
                                           negative.group.name))
    return(x)
  })
  # call here as DietSO happens below
  if (!interactive.only) {
    feat_plots <- tryCatch({
      do.call(feature_plot, args = c(list(SO = SO, assay = assay, features = colname), dots[which(names(dots) %in% names(formals(feature_plot)))]))
    }, error = function(e) {
      NULL
    })
  }

  intersect_features <- Reduce(intersect, lapply(SO, rownames))
  if (length(unique(c(sapply(SO, nrow), length(intersect_features)))) != 1) {
    warning("Different features across SOs detected. Will carry on with common ones only.")
  }

  SO <- lapply(SO, function(x) Seurat::DietSeurat(subset(x, features = intersect(rownames(x), intersect_features), cells = intersect(colnames(x), c(ngc, pgc))), assays = assay, counts = F))
  if (length(SO) > 1) {
    # this restores the counts matrix
    SO <- merge(x = SO[[1]], y = SO[2:length(SO)], merge.data = T)
  } else {
    SO <- SO[[1]]
  }

  # exclude features which are below min.pct.set in both populations
  SO <- Seurat::DietSeurat(SO, assays = assay, features = unique(c(names(which(Matrix::rowSums(Seurat::GetAssayData(SO)[, ngc] != 0)/length(ngc) > min.pct)),
                                                                   names(which(Matrix::rowSums(Seurat::GetAssayData(SO)[, pgc] != 0)/length(pgc) > min.pct)))), counts = F)

  # filter meta.col
  SO@meta.data <- SO@meta.data[,c(colname, meta.cols),drop=F]
  Seurat::Misc(SO, "volcano.group.col") <- colname
  SO@meta.data[,colname] <- factor(SO@meta.data[,colname], levels = c(negative.group.name, positive.group.name))

  # get Assay data, overwrite SO to save memory
  # Diet Seurat
  # subset Seurat, keeping ngc, pgc and additional cols if needed, assing meta.data cols for pgn, ngn
  #SO <- lapply(SO, function(x) Seurat::GetAssayData(x, assay = assay, slot = "data")[,intersect(colnames(Seurat::GetAssayData(x, assay = assay)), c(ngc, pgc))])
  #SO <- lapply(SO, function(x) x[intersect_features,])
  #SO <- do.call(cbind, SO)

  # equal order of intersecting features which are taken into non-log space for wilcox test and FC calculation
  if (is.null(volcano.data)) {
    vd <- .calc_vd(assay_data = Seurat::GetAssayData(SO, slot = "data"),
                   ngc = ngc,
                   pgc = pgc,
                   pgn = positive.group.name,
                   ngn = negative.group.name,
                   inf.fc.shift = inf.fc.shift,
                   n.feat.for.p.adj = length(intersect_features),
                   p.adjust = p.adjust,
                   use.limma = use.limma)
  } else {
    vd <- volcano.data
  }

  if (!interactive.only) {
    vp <- .plot_vp(vd = vd,
                   y = p.plot,
                   x.axis.symmetric = x.axis.symmetric,
                   y.axis.pseudo.log = y.axis.pseudo.log,
                   pseudo.log.sigma = pseudo.log.sigma,
                   pt.size = pt.size,
                   pt.alpha = pt.alpha,
                   font.size = font.size,
                   font.family = font.family,
                   ngn = negative.group.name,
                   pgn = positive.group.name,
                   x.axis.extension = x.axis.extension,
                   pval.tick = pval.tick,
                   features.exclude = features.exclude,
                   min.pct = min.pct,
                   p.cut = p.cut,
                   fc.cut = fc.cut)


    vp <- .label_vp(vp = vp,
                    vd = vd,
                    p.plot = p.plot,
                    label.features = label.features,
                    labels.topn = labels.topn,
                    label.size = label.size,
                    features.exclude = features.exclude,
                    nudge.x = nudge.x,
                    nudge.y = nudge.y,
                    max.iter = max.iter,
                    label.neg.pos.sep = label.neg.pos.sep,
                    label.col = label.col,
                    label.face = label.face,
                    font.family = font.family,
                    max.overlaps = max.overlaps,
                    p.signif = p.signif)


    ## to fgsea_on_msigdbr
    if (!is.null(gsea.param)) {
      # run on positive and negative DEG - which ones? additional arguments?
      # 1 log2.fc
      # 2 adj.p.val
      # 1 negative
      # 2 positive

      ### negative/positive volcano not working yet

      vd.gsea1 <- vd[intersect(which(vd$log2.fc <= gsea.param[[1]][1]), which(vd[,p.plot,drop=T] <= gsea.param[[2]][1])),]
      vd.gsea2 <- vd[intersect(which(vd$log2.fc >= gsea.param[[1]][2]), which(vd[,p.plot,drop=T] <= gsea.param[[2]][2])),]
      ranks.1 <- sort(stats::setNames(vd.gsea1[,"log2.fc",drop=T],vd.gsea1[,"Feature",drop=T]))
      ranks.2 <- sort(stats::setNames(vd.gsea2[,"log2.fc",drop=T],vd.gsea2[,"Feature",drop=T]))

      if (!"maxSize" %in% names(dots)) {
        arg_add <- list(gene.ranks = ranks.1, maxSize = 500)
      } else {
        arg_add <- list(gene.ranks = ranks.1)
      }
      g1 <- do.call(fgsea_on_msigdbr, args = c(arg_add, dots[which(names(dots) %in% names(formals(fgsea_on_msigdbr)))]))
      if (!"maxSize" %in% names(dots)) {
        arg_add <- list(gene.ranks = ranks.2, maxSize = 500)
      } else {
        arg_add <- list(gene.ranks = ranks.2)
      }
      g2 <- do.call(fgsea_on_msigdbr, args = c(arg_add, dots[which(names(dots) %in% names(formals(fgsea_on_msigdbr)))]))
      vp_gsea.1 <- .plot_vp(vd = vd,
                            y = p.plot,
                            x.axis.symmetric = x.axis.symmetric,
                            y.axis.pseudo.log = y.axis.pseudo.log,
                            pseudo.log.sigma = pseudo.log.sigma,
                            pt.size = pt.size,
                            pt.alpha = pt.alpha,
                            font.size = font.size,
                            font.family = font.family,
                            ngn = negative.group.name,
                            pgn = positive.group.name,
                            x.axis.extension = x.axis.extension,
                            pval.tick = pval.tick,
                            features.exclude = features.exclude,
                            min.pct = min.pct,
                            p.cut = gsea.param[[2]][1],
                            fc.cut = gsea.param[[1]][1])
      vp_gsea.1 <- .label_vp(vp = vp_gsea.1,
                             vd = vd,
                             p.plot = p.plot,
                             label.features = names(ranks.1),
                             label.size = label.size,
                             features.exclude = features.exclude,
                             label.neg.pos.sep = label.neg.pos.sep,
                             nudge.x = nudge.x,
                             nudge.y = nudge.y,
                             max.iter = max.iter,
                             label.col = label.col,
                             label.face = label.face,
                             font.family = font.family,
                             max.overlaps = max.overlaps,
                             color.only = T)
      vp_gsea.2 <- .plot_vp(vd = vd,
                            y = p.plot,
                            x.axis.symmetric = x.axis.symmetric,
                            y.axis.pseudo.log = y.axis.pseudo.log,
                            pseudo.log.sigma = pseudo.log.sigma,
                            pt.size = pt.size,
                            pt.alpha = pt.alpha,
                            font.size = font.size,
                            font.family = font.family,
                            ngn = negative.group.name,
                            pgn = positive.group.name,
                            x.axis.extension = x.axis.extension,
                            pval.tick = pval.tick,
                            features.exclude = features.exclude,
                            min.pct = min.pct,
                            p.cut = gsea.param[[2]][2],
                            fc.cut = gsea.param[[1]][2])
      vp_gsea.2 <- .label_vp(vp = vp_gsea.2,
                             vd = vd,
                             p.plot = p.plot,
                             label.features = names(ranks.2),
                             label.size = label.size,
                             features.exclude = features.exclude,
                             label.neg.pos.sep = label.neg.pos.sep,
                             nudge.x = nudge.x,
                             nudge.y = nudge.y,
                             max.iter = max.iter,
                             label.col = label.col,
                             label.face = label.face,
                             font.family = font.family,
                             max.overlaps = max.overlaps,
                             color.only = T)

      gsea <- list(negative.gsea = g1, positive.gsea = g2, negative.volcano = vp_gsea.1, positive.volcano = vp_gsea.2)
    } else {
      gsea <- NULL
    }
  }

  if (interactive.only) {
    #return(Seurat::DietSeurat(SO, features = rownames(vd)))
    return(list(vd = vd, ngn = negative.group.name, pgn = positive.group.name)) #ngc = ngc, pgc = pgc
  }

  # SO for interactive volcano
  Seurat::Misc(SO, "volcano.data") <- vd
  Seurat::Misc(SO, "features.exclude") <- features.exclude
  Seurat::Misc(SO, "ngn") <- negative.group.name
  Seurat::Misc(SO, "pgn") <- positive.group.name

  if (!is.null(save.path.interactive)) {
    dir.create(save.path.interactive, showWarnings = FALSE, recursive = TRUE)
    file.copy(from = system.file("extdata", "app.R", package = "scexpr"), to = file.path(save.path.interactive, "app.R"))
    saveRDS(stats::setNames(list(stats::setNames(list(list(Seurat_object = Seurat::DietSeurat(SO, assays = assay, features = rownames(vd))), reduction_plot = feat_plots), c("1", "reduction_plot"))), "1"), file = file.path(save.path.interactive, "data.rds"))
  }

  return(list(plot = vp,
              data = vd,
              feat.plots = feat_plots,
              gsea = gsea,
              interactive.data = stats::setNames(list(stats::setNames(list(list(Seurat_object = Seurat::DietSeurat(SO, assays = assay, features = rownames(vd))), reduction_plot = feat_plots), c("1", "reduction_plot"))), "1")))
}



.calc_vd <- function(assay_data,
                     ngc,
                     pgc,
                     pgn = "positive",
                     ngn = "negative",
                     n.feat.for.p.adj = NULL,
                     inf.fc.shift = 2,
                     p.adjust = "bonferroni",
                     use.limma = F) {

  if (use.limma) {
    if (!requireNamespace("BiocManager", quietly = T)) {
      utils::install.packages("BiocManager")
    }
    if (!requireNamespace("limma", quietly = T)) {
      BiocManager::install("limma")
    }

    message("(i) limma: p.adjust ignored. (ii) caution: limma calculates log2.fc differently as Seurat, see here: https://support.bioconductor.org/p/82478/ ")
    vd <- limma::topTable(limma::eBayes(limma::lmFit(assay_data[, c(ngc,pgc)], design = stats::model.matrix(~c(colnames(assay_data[, c(ngc,pgc)]) %in% pgc)))), number = nrow(assay_data), confint = T)

    assay_data <- expm1(assay_data) + 1
    apm <- Matrix::rowMeans(assay_data[, pgc])
    anm <- Matrix::rowMeans(assay_data[, ngc])

    vd <- data.frame(log2.fc = vd$logFC,
                     CI.L = vd$CI.L,
                     CI.R = vd$CI.R,
                     p.val = vd$P.Value,
                     adj.p.val = as.numeric(formatC(vd$adj.P.Val, format = "e", digits = 2)),
                     stats::setNames(list(round(log2(anm[rownames(vd)]), 2)), ngn),
                     stats::setNames(list(round(log2(apm[rownames(vd)]), 2)), pgn),
                     stats::setNames(list(round(apply(assay_data[, ngc]-1, 1, function(x) sum(x != 0))/ncol(assay_data[, ngc]), 2)[rownames(vd)]), paste0("pct.", ngn)),
                     stats::setNames(list(round(apply(assay_data[, pgc]-1, 1, function(x) sum(x != 0))/ncol(assay_data[, pgc]), 2)[rownames(vd)]), paste0("pct.", pgn)),
                     infinite.FC = ifelse(is.infinite(log2(apm) - log2(anm)), 1, 0)[rownames(vd)],
                     check.names = F)

  } else {
    if (!requireNamespace("matrixTests", quietly = T)) {
      utils::install.packages("matrixTests")
    }

    p.adjust <- match.arg(p.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    assay_data <- expm1(assay_data) + 1
    p <- matrixTests::row_wilcoxon_twosample(as.matrix(assay_data[, ngc]), as.matrix(assay_data[, pgc]))$pvalue
    apm <- Matrix::rowMeans(assay_data[, pgc])
    anm <- Matrix::rowMeans(assay_data[, ngc])

    if (is.null(n.feat.for.p.adj)) {
      n.feat.for.p.adj <- nrow(assay_data)
    }

    vd <- data.frame(log2.fc = log2(apm) - log2(anm), #log2(apm / anm),
                     p.val = p,
                     adj.p.val = as.numeric(formatC(stats::p.adjust(p, method = p.adjust, n = n.feat.for.p.adj), format = "e", digits = 2)),
                     stats::setNames(list(round(log2(anm), 2)), ngn),
                     stats::setNames(list(round(log2(apm), 2)), pgn),
                     stats::setNames(list(round(apply(assay_data[, ngc]-1, 1, function(x) sum(x != 0))/ncol(assay_data[, ngc]), 2)), paste0("pct.", ngn)),
                     stats::setNames(list(round(apply(assay_data[, pgc]-1, 1, function(x) sum(x != 0))/ncol(assay_data[, pgc]), 2)), paste0("pct.", pgn)),
                     infinite.FC = ifelse(is.infinite(log2(apm) - log2(anm)), 1, 0),
                     check.names = F)
  }



  vd$log2.fc <- scales::oob_squish_infinite(vd$log2.fc, range = c(min(vd$log2.fc[!is.infinite(vd$log2.fc)], na.rm = T) - inf.fc.shift,
                                                                  max(vd$log2.fc[!is.infinite(vd$log2.fc)], na.rm = T) + inf.fc.shift))
  return(as.matrix(vd))
}



.plot_vp <- function (vd,
                      x = "log2.fc",
                      y = "adj.p.val",
                      x_label = NULL,
                      y_label = NULL,
                      x.axis.symmetric = T,
                      y.axis.pseudo.log = F,
                      pseudo.log.sigma = 1,
                      features.to.color = NULL, #which features to plot with color on top
                      features.color.by = NULL, #which column to color them by
                      errorbar.low.col = NULL, # absolute coordinate of lower errorbar
                      errorbar.up.col = NULL, # absolute coordinate of upper errorbar
                      errorbar.size = 0.2,
                      errorbar.width = 0.2,
                      col.pal = "RdBu",
                      col.pal.rev = T,
                      col.type = c("c", "d"), # continuous or discrete
                      pt.size = 1,
                      pt.alpha = 0.8,
                      font.size = 14,
                      font.family = "Courier",
                      ngn = "negative.group",
                      pgn = "positive.group",
                      x.axis.extension = 0,
                      pval.tick = 0.01,
                      features.exclude = NULL,
                      min.pct = 0,
                      p.cut = NA,
                      fc.cut = NA) {

  x <- match.arg(x, colnames(vd))
  y <- match.arg(y, colnames(vd)) #c("adj.p.val", "p.val")
  col.type <- match.arg(col.type, c("c", "d"))

  vd <- as.data.frame(vd)
  if (!"Feature" %in% names(vd)) {
    # also row numbers would become a Feature column, but not relevant
    vd$Feature <- rownames(vd)
  }

  if (!is.null(features.exclude)) {
    print(paste0("The following features are excluded from the volcano plot: ", paste(vd[which(grepl(paste(features.exclude, collapse = "|"), vd$Feature)),"Feature"], collapse = ",")))
    vd <- vd[which(!grepl(paste0(features.exclude, collapse = "|"), vd$Feature)),]
  }

  if (!is.null(features.to.color)) {
    if (any(!features.to.color %in% vd$Feature)) {
      message("features.to.color: ", paste(features.to.color[which(!features.to.color %in% vd$Feature)], collapse = ", "), " not found.")
    }
    features.to.color <- features.to.color[which(features.to.color %in% vd$Feature)]
    if (length(features.to.color) == 0) {
      features.to.color <- NULL
    }
  }

  if (!is.null(features.color.by) && !features.color.by %in% names(vd)) {
    message("features.color.by not found as column.")
    features.color.by <- NULL
  }

  if (!is.null(errorbar.low.col) && !errorbar.low.col %in% names(vd)) {
    message("errorbar.low.col not found as column. Both errorbar limit will be ignored.")
    errorbar.low.col <- NULL
    errorbar.up.col <- NULL
  }

  if (!is.null(errorbar.up.col) && !errorbar.up.col %in% names(vd)) {
    message("errorbar.up.col not found as column. Both errorbar limit will be ignored.")
    errorbar.low.col <- NULL
    errorbar.up.col <- NULL
  }


  if (!is.null(ngn) && !is.null(pgn) && min.pct > 0) {
    vd <- rbind(vd[intersect(which(vd[,paste0("pct.", ngn)] >= min.pct), which(vd[,x] < 0)),],
                vd[intersect(which(vd[,paste0("pct.", pgn)] >= min.pct), which(vd[,x] > 0)),])
  }

  vp <- ggplot2::ggplot(vd, ggplot2::aes(x = !!rlang::sym(x), y = round(scales::oob_squish_infinite(-log10(!!rlang::sym(y)), range = c(0,300)), 2), label = Feature)) +
    ggplot2::geom_point(color = "#999999", alpha = pt.alpha, size = pt.size) +
    ggplot2::geom_point(data = vd[which(vd$infinite.FC == 1),], color = "cornflowerblue", size = pt.size) +
    ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


  if (!is.null(features.color.by) && !is.null(features.to.color)) {
    # select color scale
    if (col.type == "c") {
      if (length(col.pal) == 1 && !col.pal %in% grDevices::colors()) {
        col.pal <- col_pal(name = col.pal, reverse = col.pal.rev)
      }
    } else if (col.type == "d") {
      if (length(col.pal) == 1 && !col.pal %in% grDevices::colors()) {
        col.pal <- col_pal(name = col.pal, reverse = col.pal.rev, n = nlevels(as.factor(vd[which(vd$Feature %in% features.to.color),features.color.by])))
      }
    }
    if (col.type == "c") {
      vp <-
        vp +
        ggplot2::geom_point(data = vd[which(vd$Feature %in% features.to.color),], aes(color = !!rlang::sym(features.color.by)), size = pt.size) +
        ggplot2::scale_color_gradientn(colors = col.pal)
    }
    if (col.type == "d") {
      vd[,features.color.by] <- as.factor(vd[,features.color.by])
      vp <-
        vp +
        ggplot2::geom_point(data = vd[which(vd$Feature %in% features.to.color),], aes(color = !!rlang::sym(features.color.by)), size = pt.size) +
        ggplot2::scale_color_manual(values = col.pal)
    }

    if (!is.null(errorbar.up.col)) {
      # checking one of errorbar.up.col, errorbar.low.col is enough
      vp <-
        vp +
        ggplot2::geom_errorbar(data = vd[which(vd$Feature %in% features.to.color),], aes(color = !!rlang::sym(features.color.by),
                                                                                         xmin = !!rlang::sym(errorbar.low.col),
                                                                                         xmax = !!rlang::sym(errorbar.up.col)),
                               size = errorbar.size, width = errorbar.width)
      ## 95 % conf-interval in case of metavolcanoR
    }
  }

  if (!is.null(pval.tick) && pval.tick > 0) {
    gg.brk <- ggplot2::ggplot_build(vp)[["layout"]][["panel_params"]][[1]][["y"]][["breaks"]]
    ord <- order(c(gg.brk, -log10(pval.tick)))
    brk <- c(gg.brk, -log10(pval.tick))
    lab <- c(gg.brk, paste0("p = ", pval.tick))
    lab <- gsub("^0$", "", lab)
  }


  if (is.null(y_label)) {
    if (y == "adj.p.val") {
      y_label <- "(adj p-val)"
    } else if (y == "p.val") {
      y_label <- "(p-val)"
    } else {
      y_label <- "p"
    }
  }


  if (y.axis.pseudo.log) {
    vp <- vp + ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = pseudo.log.sigma),
                                           breaks = brk[ord], labels = lab[ord], limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  } else {
    vp <- vp + ggplot2::scale_y_continuous(breaks = brk[ord], labels = lab[ord], limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  }

  if (is.null(ngn) && is.null(pgn)) {
    if (is.null(x_label)) {
      vp <- vp + ggplot2::labs(x = bquote("avg" ~ log[2] ~ "FC"), y = bquote(-log[10]~.(rlang::sym(y_label))))
    } else {
      vp <- vp + ggplot2::labs(x = x_label, y = bquote(-log[10]~.(rlang::sym(y_label))))
    }
  } else {
    if (is.null(x_label)) {
      vp <- vp + ggplot2::labs(x = bquote(bold(.(ngn)) ~ "  <====  " ~ log[2] ~ "FC" ~ "  ====>  " ~ bold(.(pgn))), y = bquote(-log[10]~.(rlang::sym(y_label))))
    } else {
      vp <- vp + ggplot2::labs(x = bquote(bold(.(ngn)) ~ "  <====  " ~ .(rlang::sym(x_label)) ~ "  ====>  " ~ bold(.(pgn))), y = bquote(-log[10]~.(rlang::sym(y_label))))
    }
  }

  if (x.axis.symmetric) {
    vp <- vp + ggplot2::xlim(-round(max(abs(vd[,x]))) - 0.5 - x.axis.extension, round(max(abs(vd[,x]))) + 0.5 + x.axis.extension)
  } else {
    vp <- vp + ggplot2::xlim(round(min(vd[,x])) - 0.5 - x.axis.extension, round(max(vd[,x])) + 0.5 + x.axis.extension)
  }

  if (!any(c(is.na(p.cut), is.na(fc.cut)))) {
    vp <- vp + ggplot2::geom_hline(yintercept = -log10(p.cut), linetype = "dashed") + ggplot2::geom_vline(xintercept = c(-fc.cut, fc.cut), linetype = "dashed")

    vd.cut <-
      vd %>%
      dplyr::filter(abs(!!rlang::sym(x)) >= fc.cut) %>%
      dplyr::filter(!!rlang::sym(y) <= p.cut)

    vd.cut.sum <-
      vd.cut %>%
      dplyr::mutate(sign = ifelse(!!rlang::sym(x) > 0, "plus", "minus")) %>%
      dplyr::group_by(sign) %>%
      dplyr::summarise(n.genes = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(xpos = ifelse(sign == "plus", fc.cut + 0.75*(max(abs(vd.cut[,x]))-fc.cut), -(fc.cut + 0.75*(max(abs(vd.cut[,x]))-fc.cut)))) %>%
      dplyr::mutate(ypos = -log10(p.cut)*2)

    vp <-
      vp +
      ggplot2::geom_point(data = vd.cut, color = "black", size = pt.size) +
      ggplot2::geom_text(data = vd.cut.sum, ggplot2::aes(x = xpos, y = ypos, label = n.genes), inherit.aes = F, size = 5/14*font.size, family = font.family)

  }
  return(vp)
}

.label_vp <- function(vp,
                      vd,
                      p.plot = "adj.p.val",
                      label.features = NULL,
                      topn.metric = "p.value",
                      labels.topn = 30,
                      label.size = 4,
                      nudge.x = 0,
                      nudge.y = 0,
                      max.iter = 10000,
                      label.neg.pos.sep = T,
                      label.col = "black",
                      label.face = "bold",
                      font.family = "Courier",
                      max.overlaps = 50,
                      p.signif = 0.001,
                      features.exclude = NULL,
                      color.only = F,
                      label.only = F) {

  vd <- as.data.frame(vd)
  if (!"Feature" %in% names(vd)) {
    # also row numbers would become a Feature column, but not relevant
    vd$Feature <- rownames(vd)
  }

  if (!is.null(features.exclude)) {
    vd <- vd[which(!grepl(paste0(features.exclude, collapse = "|"), vd$Feature)),]
  }

  if (is.null(label.features)) {
    if (topn.metric == "p.value") {
      f_lab <- vd %>% dplyr::top_n(-labels.topn, !!rlang::sym(p.plot))
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    } else if (topn.metric == "both") {
      f_lab.p.val <- vd %>% dplyr::top_n(-labels.topn, !!rlang::sym(p.plot))
      f_lab.logfc <- dplyr::bind_rows(vd %>% dplyr::top_n(labels.topn/2, log2.fc), vd %>% dplyr::top_n(-(labels.topn/2), log2.fc))
      f_lab <- dplyr::bind_rows(f_lab.logfc, f_lab.p.val) %>% dplyr::distinct()
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    } else {
      f_lab <- dplyr::bind_rows(vd %>% dplyr::top_n(labels.topn/2, log2.fc), vd %>% dplyr::top_n(-(labels.topn/2), log2.fc))
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    }
    if (!label.only) {
      vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% f_lab$Feature), colour = "tomato2")
    }
    if (!color.only) {
      if (label.neg.pos.sep) {
        vp <- vp +
          ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% f_lab.pos), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
          ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% f_lab.neg), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
      } else {
        vp <- vp +
          ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% f_lab$Feature), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge.x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
      }
    }
  } else {
    if (length(label.features) == 1) {
      if (label.features == "significant") {
        label.features <- vd[which(as.numeric(vd[,p.plot]) < p.signif), "Feature"]
      } else {
        if (!label.only) {
          vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = "tomato2")
        }
        if (!color.only) {
          if (label.neg.pos.sep) {
            vp <- vp +
              ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc > 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
              ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc < 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
          } else {
            vp <- vp + ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
          }
        }
      }
    } else {
      if (!label.only) {
        vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = "tomato2")
      }
      if (!color.only) {
        if (label.neg.pos.sep) {
          vp <- vp +
            ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc > 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
            ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc < 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
        } else {
          vp <- vp + ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
        }
      }
    }
  }
  return(vp)
}

