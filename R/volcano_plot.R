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
#' @param ...
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
  if (any(duplicated(negative.group.cells))) {
    print("Duplicates found in negative.group.cells")
  }
  ngc <- unique(negative.group.cells)
  if (any(duplicated(positive.group.cells))) {
    print("Duplicates found in positive.group.cells")
  }
  pgc <- unique(positive.group.cells)
  if (length(intersect(negative.group.cells, positive.group.cells)) > 0) {
    stop("Equal cell names in negative.group.cells and positive.group.cells. Please fix this.")
  }

  all.cells <- unlist(lapply(SO, function(x) Seurat::Cells(x)))
  if (length(SO) > 1 && any(duplicated(all.cells))) {
    stop("Cell names are not unique across SOs. Please fix that manually with Seurat::RenameCells. Then pass SOs again.")
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
  # call here as SO is overwritten below
  if (!interactive.only) {
    feat_plots <- do.call(feature_plot, args = c(list(SO = SO, assay = assay, features = colname), dots[which(names(dots) %in% names(formals(feature_plot)))]))
  }

  # get Assay data, overwrite SO to save memory
  SO <- lapply(SO, function(x) Seurat::GetAssayData(x, assay = assay, slot = "data")[,intersect(colnames(Seurat::GetAssayData(x, assay = assay)), c(ngc, pgc))])

  intersect_features <- Reduce(intersect, lapply(SO, rownames))
  if (length(unique(c(sapply(SO, nrow), length(intersect_features)))) != 1) {
    print("Different features across SOs detected. Will only carry on with common ones.")
    SO <- lapply(SO, function(x) x[intersect_features,])
  }
  SO <- do.call(cbind, SO)
  # exclude features which are below min.pct.set in both populations
  SO <- SO[unique(c(which(Matrix::rowSums(SO[, ngc] != 0)/length(ngc) > min.pct),
                    which(Matrix::rowSums(SO[, pgc] != 0)/length(pgc) > min.pct))), ]
  # equal order of intersecting features which are taken into non-log space for wilcox test and FC calculation
  if (is.null(volcano.data)) {
    vd <- .calc_vd(assay_data = expm1(SO),
                   ngc = ngc,
                   pgc = pgc,
                   pgn = positive.group.name,
                   ngn = negative.group.name,
                   inf.fc.shift = inf.fc.shift,
                   n.feat.for.p.adj = length(intersect_features),
                   p.adjust = p.adjust)
  } else {
    vd <- volcano.data
  }

  if (!interactive.only) {
    vp <- .plot_vp(vd = vd,
                   p.plot = p.plot,
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
                            p.plot = p.plot,
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
                            p.plot = p.plot,
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

  interactive.data = list(non.aggr.data = SO[vd$Feature,c(ngc, pgc)], data = vd, negative.group.cells = ngc, positive.group.cells = pgc, negative.group.name = negative.group.name, positive.group.name = positive.group.name, features.exclude = features.exclude)
  if (!is.null(save.path.interactive)) {
    dir.create(save.path.interactive, showWarnings = FALSE, recursive = TRUE)
    file.copy(from = system.file("extdata", "app.R", package = "scexpr"), to = file.path(save.path.interactive, "app.R"))
    saveRDS(stats::setNames(list(stats::setNames(list(interactive.data), "1")), "1"), file = file.path(save.path.interactive, "data.rds"))
  }
  if (interactive.only) {
    return(interactive.data)
  }
  return(list(plot = vp,
              data = vd,
              feat.plots = feat_plots,
              gsea = gsea,
              interactive.data = interactive.data))
}



.calc_vd <- function(assay_data,
                     ngc,
                     pgc,
                     pgn = "positive",
                     ngn = "negative",
                     n.feat.for.p.adj = NULL,
                     inf.fc.shift = 2,
                     p.adjust = "bonferroni") {

  p.adjust <- match.arg(p.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))

  p <- matrixTests::row_wilcoxon_twosample(as.matrix(assay_data[, ngc]), as.matrix(assay_data[, pgc]))$pvalue
  apm <- apply(assay_data[, pgc], 1, mean)
  anm <- apply(assay_data[, ngc], 1, mean)
  if (is.null(n.feat.for.p.adj)) {
    n.feat.for.p.adj <- nrow(an)
  }

  vd <- data.frame(Feature = rownames(assay_data),
                   log2.fc = log2(apm / anm),
                   p.val = p,
                   adj.p.val = as.numeric(formatC(stats::p.adjust(p, method = p.adjust, n = n.feat.for.p.adj), format = "e", digits = 2)),
                   stats::setNames(list(round(anm, 2)), ngn),
                   stats::setNames(list(round(apm, 2)), pgn),
                   stats::setNames(list(round(apply(assay_data[, ngc], 1, function(x) sum(x != 0))/ncol(assay_data[, ngc]), 2)), paste0("pct.", ngn)),
                   stats::setNames(list(round(apply(assay_data[, pgc], 1, function(x) sum(x != 0))/ncol(assay_data[, pgc]), 2)), paste0("pct.", pgn)),
                   infinite.FC = ifelse(is.infinite(log2(apm / anm)), TRUE, FALSE),
                   check.names = F)

  vd$log2.fc <- scales::oob_squish_infinite(vd$log2.fc, range = c(min(vd$log2.fc[!is.infinite(vd$log2.fc)], na.rm = T) - inf.fc.shift,
                                                                  max(vd$log2.fc[!is.infinite(vd$log2.fc)], na.rm = T) + inf.fc.shift))
  return(vd)
}



.plot_vp <- function (vd,
                      p.plot = "adj.p.val",
                      x.axis.symmetric = T,
                      y.axis.pseudo.log = F,
                      pseudo.log.sigma = 1,
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

  p.plot <- match.arg(p.plot, c("adj.p.val", "p.val"))

  if (!is.null(features.exclude)) {
    print(paste0("The following features are excluded from the volcano plot: ", paste(vd[which(grepl(paste(features.exclude, collapse = "|"), vd$Feature)),"Feature"], collapse = ",")))
    vd <- vd[which(!grepl(paste0(features.exclude, collapse = "|"), vd$Feature)),]
  }

  vd <- rbind(vd[intersect(which(vd[,paste0("pct.", ngn)] >= min.pct), which(vd[,"log2.fc"] < 0)),],
              vd[intersect(which(vd[,paste0("pct.", pgn)] >= min.pct), which(vd[,"log2.fc"] > 0)),])

  vp <- ggplot2::ggplot(vd, ggplot2::aes(x = log2.fc, y = round(scales::oob_squish_infinite(-log10(!!rlang::sym(p.plot)), range = c(0,300)), 2), label = Feature)) +
    ggplot2::geom_point(color = "#999999", alpha = pt.alpha, size = pt.size) +
    ggplot2::geom_point(data = vd[which(vd$infinite.FC),], color = "cornflowerblue", size = pt.size) +
    ggplot2::theme_bw(base_size = font.size, base_family = font.family) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (pval.tick > 0) {
    gg.brk <- ggplot2::ggplot_build(vp)[["layout"]][["panel_params"]][[1]][["y"]][["breaks"]]
    ord <- order(c(gg.brk, -log10(pval.tick)))
    brk <- c(gg.brk, -log10(pval.tick))
    lab <- c(gg.brk, paste0("p = ", pval.tick))
    lab <- gsub("^0$", "", lab)
  }


  if (p.plot == "adj.p.val") {
    p.plot_label <- "(adj p-val)"
  }
  if (p.plot == "p.val") {
    p.plot_label <- "(p-val)"
  }

  if (y.axis.pseudo.log) {
    vp <- vp + ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = pseudo.log.sigma),
                                           breaks = brk[ord], labels = lab[ord], limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  } else {
    vp <- vp + ggplot2::scale_y_continuous(breaks = brk[ord], labels = lab[ord], limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  }

  if (is.null(ngn) && is.null(pgn)) {
    vp <- vp + ggplot2::labs(x = bquote("avg" ~ log[2] ~ "FC"), y = bquote(-log[10]~.(rlang::sym(p.plot_label))))
  } else {
    vp <- vp + ggplot2::labs(x = bquote(bold(.(ngn)) ~ "  <====  " ~ log[2] ~ "FC" ~ "  ====>  " ~ bold(.(pgn))), y = bquote(-log[10]~.(rlang::sym(p.plot_label))))
  }

  if (x.axis.symmetric) {
    vp <- vp + ggplot2::xlim(-round(max(abs(vd$log2.fc))) - 0.5 - x.axis.extension, round(max(abs(vd$log2.fc))) + 0.5 + x.axis.extension)
  } else {
    vp <- vp + ggplot2::xlim(round(min(vd$log2.fc)) - 0.5 - x.axis.extension, round(max(vd$log2.fc)) + 0.5 + x.axis.extension)
  }

  if (!any(c(is.na(p.cut), is.na(fc.cut)))) {
    vp <- vp + ggplot2::geom_hline(yintercept = -log10(p.cut), linetype = "dashed") + ggplot2::geom_vline(xintercept = c(-fc.cut, fc.cut), linetype = "dashed")

    vd.cut <-
      vd %>%
      dplyr::filter(abs(log2.fc) >= fc.cut) %>%
      dplyr::filter(!!rlang::sym(p.plot) <= p.cut)

    vd.cut.sum <-
      vd.cut %>%
      dplyr::mutate(sign = ifelse(log2.fc > 0, "plus", "minus")) %>%
      dplyr::group_by(sign) %>%
      dplyr::summarise(n.genes = n()) %>%
      dplyr::mutate(xpos = ifelse(sign == "plus", fc.cut + 0.75*(max(abs(vd.cut$log2.fc))-fc.cut), -(fc.cut + 0.75*(max(abs(vd.cut$log2.fc))-fc.cut)))) %>%
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
                      color.only = F) {

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
      f_lab.logfc <- bind_rows(vd %>% dplyr::top_n(labels.topn/2, log2.fc), vd %>% dplyr::top_n(-(labels.topn/2), log2.fc))
      f_lab <- dplyr::bind_rows(f_lab.logfc, f_lab.p.val) %>% dplyr::distinct()
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    } else {
      f_lab <- dplyr::bind_rows(vd %>% dplyr::top_n(labels.topn/2, log2.fc), vd %>% dplyr::top_n(-(labels.topn/2), log2.fc))
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    }
    vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% f_lab$Feature), colour = "tomato2")

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
        label.features <- vd[which(as.numeric(vd[,p]) < p.signif), "Feature"]
      } else {
        vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = "tomato2")
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
      vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = "tomato2")
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


