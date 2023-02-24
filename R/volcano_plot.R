#' Make a volcano plot comparing the gene expression of selected groups of cells from one or more preprocessed Seurat objects
#'
#'
#'
#' @param SO one or more Seurat object(s); provide multiple objects as a (named) list
#' @param assay which assay to get expression data from
#' @param volcano.data optionally provide a data frame with data that have been calculated before
#' with this function. e.g. if only style element of the returned plot are to be modified, providing this
#' data frame will avoid to re-calculate DE genes
#' @param negative.group.cells vector of cell names for the negative group of cells; genes which
#' are expressed at higher level in this group will have a negative sign (minus); use Seurat::Cells()
#' or Seurat::WhichCells() to select or use filter operations on SO@meta.data to select cells
#' @param positive.group.cells vector of cell names for the positive group of cells; genes which
#' are expressed at higher level in this group will have a positive sign (plus); use Seurat::Cells()
#' or Seurat::WhichCells() to select or use filter operations on SO@meta.data to select cells
#' @param negative.group.name name for the negative group of cells; this will appear on plot axes
#' @param positive.group.name name for the positive group of cells; this will appear on plot axes
#' @param x.axis.symmetric logical; make the x-axis which indicates log-fold changes symmetric (same range
#' in plus and minus direction)
#' @param x.axis.extension numeric; positive or negative value to extend the x-axis by
#' @param y.axis.pseudo.log logical; use a pseudo-log y-axis for -log10(pvalue); may be helpful
#' to better visualize and label dense area of genes in the plot; done with scales::pseudo_log_trans
#' @param pseudo.log.sigma numeric; adjustment to the pseudo-log tranasformation; passed to scales::pseudo_log_trans
#' @param inf.fc.shift numeric; genes with an infinite fold-change (happens when a gene not expressed at all in one group)
#' are not left infinite but will be put at: "largest finite fold-change +/- inf.fc.shift"; such genes are indicated
#' specifically
#' @param pt.size numeric; point size passed to ggplot2::geom_point
#' @param pt.alpha numeric; point opacity (alpha value) passed to ggplot2::geom_point
#' @param font.size numeric; font size of point labels; passed as base_size to ggplot2::theme_bw() and
#' as size = 5/14*font.size to ggplot2::geom_text (5/14 is the conversion between how geom_text handles the
#' size argument and how theme does it)
#' @param pval.tick NULL or numeric, if not NULL: plot an extra y-axis-tick to indicate a significance level
#' in the linear space (not -log10); e.g.: pval.tick = 0.01 will plot an axis-tick at p = 0.01
#' @param min.pct numeric; only consider genes for the volcano plot which are expressed at least by this
#' fraction of cells in both groups (negative.group.cells, positive.group.cells)
#' @param max.iter max.iter passed to ggrepel::geom_text_repel
#' @param max.overlaps max.overlaps passed to ggrepel::geom_text_repel
#' @param label.neg.pos.sep logical whether to label genes with negative and positive fold change separately;
#' this will make sure that labels on the left (negative log fc) point to the left and labels on the right
#' (positive log fc) point to the right
#' @param label.col font color of labels
#' @param label.face font face of labels
#' @param font.family font family (font type) of labels, e.g. Courier or Sans
#' @param label.size size of labels
#' @param labels.topn numeric; how many labels (roughly) to plot based on the selected metric
#' @param label.features character vector; select genes which are to be labeled; OR: "significant"
#' which will let p.signif come into usage
#' @param topn.metric which metric to use for selection of top genes for labeling ()
#' @param nudge.x nudge.x passed to ggrepel::geom_text_repel; shift labels into x-direction;
#' for genes with negative fold change -nudge.x will be used; genes with positive fold change
#' use nudge.x
#' @param nudge.y nudge.x passed to ggrepel::geom_text_repel; shift labels into y-direction;
#' @param p.plot which p-value to use for plotting, adjusted p-value or unadjusted p-value
#' @param p.adjust method for p-value adjustment
#' @param p.cut plot a horizontal line at this p-value; in combination with fc.cut
#' all genes above this cut (and above fc.cut) are counted and the number is plotted
#' @param p.signif a significance level from which on genes are labeled;
#' supply label.features = "significant" to make use of this
#' @param fc.cut plot vertical lines at these fold changes (plus and minus); in combination with p.cut
#' all genes above this cut (and above p.cut) are counted and the number is plotted
#' @param features.exclude character vector of features to exclude from plotting; you may
#' supply regular expressions like "^RPL" and/or "^RPS" to exclude all ribosomal genes
#' @param meta.cols which meta.cols to keep in SO for interactive plotting
#' @param save.path.interactive path on disk where to save data for interactive
#' analysis of the volcano plot; this will initiate a directory with a shiny script
#' and an rds file of with SO and DE genes
#' @param gsea.param list of length 2, each index holding a numeric vector of length 2:
#' fold change (list index 1) and p-value (list index 2) limits for GSEA; first entry of each vector
#' applies for genes with negative fold change, second entry for genes with positive fold change;
#' 2 GSEA will be run separately on genes with negative and positive fold change;
#' if argument is NULL no GSEA is performed; the function called for GSEA is scexpr::fgsea_on_msigdbr which
#' uses all gene sets from \href{https://igordot.github.io/msigdbr/articles/msigdbr-intro.html}{msigdb(r)} by default
#' and runs \href{https://bioconductor.org/packages/release/bioc/html/fgsea.html}{fgsea} for GSEA; alternative
#' arguments to scexpr::fgsea_on_msigdbr can be supplied in ... (except for gene.ranks)
#' @param interactive.only logical; do calculations only for interactive analysis and only
#' return these values
#' @param use.limma logical; use limma for DE gene detection; intended for subsequent use MetaVolcanoR
#' which relies on values returned from limma
#' @param ... arguments to scexpr::feature_plot like  col.pal.d = setNames(c("grey90", scexpr::col_pal()[c(2,3)]), c("other", "name1", "name2")), order.discrete = "^other", plot.title = F, etc
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
                         pval.tick = NULL,
                         min.pct = 0.1,
                         max.iter = 10000,
                         max.overlaps = 50,
                         label.neg.pos.sep = T,
                         label.col = "black",
                         label.face = "italic",
                         font.family = "sans",
                         label.size = 4,
                         labels.topn = 30,
                         label.features = NULL,
                         topn.metric = c("p.value", "fc", "both"),
                         nudge.x = 0,
                         nudge.y = 0,
                         p.plot = c("adj.p.val", "p.val"),
                         p.adjust = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"),
                         p.cut = NA,
                         p.signif = 0.001,
                         fc.cut = NA,
                         features.exclude = NULL,
                         save.path.interactive = NULL,
                         gsea.param = NULL,
                         interactive.only = F,
                         use.limma = F,
                         ...) {

  ### add option to only label positive or negative features
  ### attach attributes to vd; like the ngc, pgc, min.pct, assay, p.adjust, use.limma, and so on, this can then be tested for, when volcano.data as input

  ## fix shiny app for interactive volcano plot

  if (missing(negative.group.cells) || missing(positive.group.cells)) {
    stop("positive.group.cells and negative.group.cells are required.")
  }
  if (length(positive.group.name) > 1 || length(negative.group.name) > 1) {
    stop("Only provide one name for negative.group.cells and positive.group.name, each.")
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
  p.adjust <- match.arg(p.adjust, c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"))
  topn.metric <- match.arg(topn.metric, c("p.value", "fc", "both"))
  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL)

  dots <- list(...)

  # label.features
  # add option to label all features indicated with p.cut and fc.cut

  ## check cells
  if (length(intersect(negative.group.cells, positive.group.cells)) > 0) {
    stop("Same cell names (barcodes) found in positive.group.cells and negative.group.cells. Please change selection or fix SOs with Seurat::RenameCells first.")
  }

  #temp <- .check.and.get.cells(SO = SO, assay = assay, cells = c(negative.group.cells, positive.group.cells), return.included.cells.only = T)
  # done like this as cell names in SOs may be subject to prefixing when there are duplicates and .check.and.get.cells is called once
  temp <- .check.and.get.cells(SO = SO, assay = assay, cells = c(negative.group.cells, positive.group.cells), return.included.cells.only = T)
  #pgc <- do.call(.check.and.get.cells, args = c(list(SO = SO, assay = assay, cells = positive.group.cells, return.included.cells.only = T), dots[which(names(dots) %in% names(formals(.check.and.get.cells)))]))

  ## slow procedure; how to speed up? why needed actually? 2022 11 07
  # https://stackoverflow.com/questions/35726028/understanding-grep-with-fixed-t-in-r
  # https://stackoverflow.com/questions/10128617/test-if-characters-are-in-a-string

  if (grepl("SO_[[:digit:]]{1,}_", temp[1])) {
    pgc <- purrr::map_chr(positive.group.cells, function(x) {
      m <- grep(pattern = paste0("SO_[[:digit:]]{1,}_", x, "$"), x = temp, value = T)
      if (length(m) > 1) {
        stop("positive.group.cells could not be identified unambigously due to duplicate cell names (barcodes). Please change selection or fix SOs with Seurat::RenameCells first.")
      }
      return(m)
    })
    ngc <- purrr::map_chr(negative.group.cells, function(x) {
      m <- grep(pattern = paste0("SO_[[:digit:]]{1,}_", x, "$"), x = temp, value = T)
      if (length(m) > 1) {
        stop("negative.group.cells could not be identified unambigously due to duplicate cell names (barcodes). Please change selection or fix SOs with Seurat::RenameCells first.")
      }
      return(m)
    })
  } else {
    pgc <- purrr::map_chr(positive.group.cells, function(x) {
      #m <- grep(pattern = paste0(x, "$"), x = temp, value = T)
      m <- grep(pattern = x, x = temp, value = T, fixed = T)
      if (length(m) > 1) {
        stop("positive.group.cells could not be identified unambigously. That means one of the barcodes in positive.group.cells is a substring of another cells barcode.")
      }
      return(m)
    })
    ngc <- purrr::map_chr(negative.group.cells, function(x) {
      #m <- grep(pattern = paste0(x, "$"), x = temp, value = T)
      m <- grep(pattern = x, x = temp, value = T, fixed = T)
      if (length(m) > 1) {
        stop("negative.group.cells could not be identified unambigously. That means one of the barcodes in negative.group.cells is a substring of another cells barcode.")
      }
      return(m)
    })
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
    if(!is.matrix(volcano.data)) {
      stop("volcano.data has to be a matrix.")
    }
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
                    p.signif = p.signif,
                    topn.metric = topn.metric)


    ## to fgsea_on_msigdbr
    if (!is.null(gsea.param)) {
      # run on positive and negative DEG - which ones? additional arguments?
      # 1 log2.fc
      # 2 adj.p.val
      # 1 negative
      # 2 positive

      ### negative/positive volcano not working yet
      vd.gsea1 <- vd[intersect(which(vd[,"log2.fc"] <= gsea.param[[1]][1]), which(vd[,p.plot,drop=T] <= gsea.param[[2]][1])),]
      vd.gsea2 <- vd[intersect(which(vd[,"log2.fc"] >= gsea.param[[1]][2]), which(vd[,p.plot,drop=T] <= gsea.param[[2]][2])),]
      ranks.1 <- sort(stats::setNames(vd.gsea1[,"log2.fc",drop=T],rownames(vd.gsea1)))
      ranks.2 <- sort(stats::setNames(vd.gsea2[,"log2.fc",drop=T],rownames(vd.gsea2)))

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
                             topn.metric = topn.metric,
                             plot.label = F)
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
                             topn.metric = topn.metric,
                             plot.label = F)

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

    p.adjust <- match.arg(p.adjust, c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"))
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
                      font.family = "sans", # serif = Times New Roman, sans = Arial, mono = Courier New
                      ngn = "negative.group",
                      pgn = "positive.group",
                      x.axis.extension = 0,
                      pval.tick = NULL,
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
        ggplot2::geom_point(data = vd[which(vd$Feature %in% features.to.color),], ggplot2::aes(color = !!rlang::sym(features.color.by)), size = pt.size) +
        ggplot2::scale_color_gradientn(colors = col.pal)
    }
    if (col.type == "d") {
      vd[,features.color.by] <- as.factor(vd[,features.color.by])
      vp <-
        vp +
        ggplot2::geom_point(data = vd[which(vd$Feature %in% features.to.color),], ggplot2::aes(color = !!rlang::sym(features.color.by)), size = pt.size) +
        ggplot2::scale_color_manual(values = col.pal)
    }

    if (!is.null(errorbar.up.col)) {
      # checking one of errorbar.up.col, errorbar.low.col is enough
      vp <-
        vp +
        ggplot2::geom_errorbar(data = vd[which(vd$Feature %in% features.to.color),], ggplot2::aes(color = !!rlang::sym(features.color.by),
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
    brk <- brk[ord]
    lab <- lab[ord]
  } else {
    brk <- ggplot2::waiver()
    lab <- ggplot2::waiver()
  }


  if (is.null(y_label)) {
    if (grepl("adj", y, ignore.case = T) {
      y_label <- "(adj p-val)"
    } else {
      y_label <- "(p-val)"
    }
  }


  if (y.axis.pseudo.log) {
    vp <- vp + ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = pseudo.log.sigma),
                                           breaks = brk, labels = lab, limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  } else {
    vp <- vp + ggplot2::scale_y_continuous(breaks = brk, labels = lab, limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
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
                      x = "log2.fc",
                      p.plot = "adj.p.val",
                      label.features = NULL,
                      topn.metric = "p.value",
                      labels.topn = 30,
                      label.size = 4,
                      nudge.x = 0,
                      nudge.y = 0,
                      max.iter = 10000,
                      dot.color = "tomato2", # set NA to have no color
                      label.neg.pos.sep = T,
                      label.col = "black",
                      label.face = "bold",
                      font.family = "sans", # serif = Times New Roman, sans = Arial, mono = Courier New
                      max.overlaps = 50,
                      p.signif = 0.001,
                      features.exclude = NULL,
                      plot.label = T) {

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
      f_lab.pos <- f_lab %>% dplyr::filter(!!rlang::sym(x) > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(!!rlang::sym(x) < 0) %>% dplyr::pull(Feature)
    } else if (topn.metric == "both") {
      f_lab.p.val <- vd %>% dplyr::top_n(-labels.topn, !!rlang::sym(p.plot))
      f_lab.logfc <- dplyr::bind_rows(vd %>% dplyr::top_n(labels.topn/2, !!rlang::sym(x)), vd %>% dplyr::top_n(-(labels.topn/2), !!rlang::sym(x)))
      f_lab <- dplyr::bind_rows(f_lab.logfc, f_lab.p.val) %>% dplyr::distinct()
      f_lab.pos <- f_lab %>% dplyr::filter(!!rlang::sym(x) > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(!!rlang::sym(x) < 0) %>% dplyr::pull(Feature)
    } else if (topn.metric == "fc") {
      f_lab <- dplyr::bind_rows(vd %>% dplyr::top_n(labels.topn/2, !!rlang::sym(x)), vd %>% dplyr::top_n(-(labels.topn/2), !!rlang::sym(x)))
      f_lab.pos <- f_lab %>% dplyr::filter(!!rlang::sym(x) > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(!!rlang::sym(x) < 0) %>% dplyr::pull(Feature)
    }
    vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% f_lab$Feature), colour = dot.color)
    if (plot.label) {
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
        vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = dot.color)
        if (plot.label) {
          if (label.neg.pos.sep) {
            vp <- vp +
              ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & !!rlang::sym(x) > 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
              ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & !!rlang::sym(x) < 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
          } else {
            vp <- vp + ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
          }
        }
      }
    } else {
      vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = dot.color)
      if (plot.label) {
        if (label.neg.pos.sep) {
          vp <- vp +
            ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & !!rlang::sym(x) > 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
            ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & !!rlang::sym(x) < 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
        } else {
          vp <- vp + ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
        }
      }
    }
  }
  return(vp)
}

