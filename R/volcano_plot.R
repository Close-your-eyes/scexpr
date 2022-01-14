volcano_plot <- function(SO,
                         assay = "RNA",
                         volcano.data = NULL,
                         negative.group.cells,
                         positive.group.cells,
                         negative.group.name = "negative.group",
                         positive.group.name = "positive.group",
                         plot.symmetric = T,
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
                         pos.neg.label.separate = T,
                         label.col = "black",
                         label.face = "bold",
                         font.family = "Courier",
                         label.size = 4,
                         lables.topn = 30,
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
                         ...) {

  assay <- match.arg(assay, c("RNA", "SCT"))
  p.plot <- match.arg(p.plot, c("adj.p.val", "p.val"))
  p.adjust <- match.arg(p.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL)


  # label.features
  # add option to label all features indicated with p.cut and fc.cut

  ## check cells
  if (any(duplicated(negative.group.cells))) {
    print("Duplicates found in negative.group.cells")
    negative.group.cells <- unique(negative.group.cells)
  }
  if (any(duplicated(positive.group.cells))) {
    print("Duplicates found in positive.group.cells")
    positive.group.cells <- unique(positive.group.cells)
  }

  if (length(intersect(negative.group.cells, positive.group.cells)) > 0) {
    stop("Equal cell names in negative.group.cells and positive.group.cells. Please fix this.")
  }

  all.cells <- unlist(lapply(SO, function(x) Seurat::Cells(x)))
  if (length(SO) > 1 && any(duplicated(all.cells)) && !make.cells.unique) {
    stop("Cell names are not unique across SOs. Please fix that manually with Seurat::RenameCells. Then pass SOs again.")
  }


  # make meta data column in SOs to identify negative.group.cells and positive.group.cells
  colname <- paste0("cellident_", format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d%H%M%S"))
  SO <- lapply(SO, function(x) {
    x@meta.data[,colname] <- ifelse(!rownames(x@meta.data) %in% c(positive.group.cells, negative.group.cells),
                                    "other",
                                    ifelse(rownames(x@meta.data) %in% positive.group.cells,
                                           positive.group.name,
                                           negative.group.name))
    return(x)
  })

  feat_plots <- feature_plot(SO,
                             features = colname,
                             order.discrete = c("other", negative.group.name, positive.group.name),
                             ...)


  ## check if data has been scaled (compare to counts)
  check <- unlist(lapply(SO, function(x) identical(Seurat::GetAssayData(x, assay = assay, slot = "data"), Seurat::GetAssayData(x, assay = assay, slot = "counts"))))
  if (any(check)) {
    stop("data slot in at least one SO does not seem to contain normalized data since it is equal to the counts slot. You may want to normalize before using volcano_plot.")
  }

  # get Assay data (check features first)
  assay_data <- lapply(SO, function(x) as.matrix(Seurat::GetAssayData(x, assay = assay, slot = "data"))[,intersect(colnames(Seurat::GetAssayData(x, assay = assay, slot = "data")), c(negative.group.cells, positive.group.cells))])

  all_features <- lapply(assay_data, rownames)
  n_features <- sapply(assay_data, nrow)
  if (length(unique(n_features)) != 1) {
    print("Different number of features detected in SOs.")
  }

  intersect_features <- Reduce(intersect, all_features)
  if (length(unique(c(n_features, length(intersect_features)))) != 1) {
    print("Different features across SOs detected. Will only carry on with common ones.")
    assay_data <- lapply(assay_data, function(x) x[intersect_features,])
  }
  assay_data <- do.call(cbind, assay_data)

  #get Assay data
  an <- assay_data[, negative.group.cells]
  ap <- assay_data[, positive.group.cells]

  # exclude features which are below min.pct.set in both populations
  # # equal order of intersecting features which are taken into non-log space for wilcox test and FC calculation
  neg.pct <- which(Matrix::rowSums(an != 0)/ncol(an) > min.pct)
  pos.pct <- which(Matrix::rowSums(ap != 0)/ncol(ap) > min.pct)
  an <- expm1(an[unique(c(neg.pct, pos.pct)),])
  ap <- expm1(ap[unique(c(neg.pct, pos.pct)),])

  # for non.aggr.data in interactive volcano
  assay_data <- assay_data[,c(negative.group.cells, positive.group.cells)]

  if (is.null(volcano.data)) {
    vd <- .calc_vd(an = an,
                   ap = ap,
                   pgn = positive.group.name,
                   ngn = negative.group.name,
                   inf.fc.shift = inf.fc.shift,
                   n.feat.for.p.adj = length(intersect_features),
                   p.adjust = p.adjust)
  } else {
    vd <- volcano.data
  }


  vp <- .plot_vp(vd = vd,
                 p.plot = p.plot,
                 plot.symmetric = plot.symmetric,
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
                  lables.topn = lables.topn,
                  label.size = label.size,
                  features.exclude = features.exclude,
                  nudge.x = nudge.x,
                  nudge.y = nudge.y,
                  max.iter = max.iter,
                  pos.neg.label.separate = pos.neg.label.separate,
                  label.col = label.col,
                  label.face = label.face,
                  font.family = font.family,
                  max.overlaps = max.overlaps,
                  p.signif = p.signif)



  interactive.data = list(non.aggr.data = assay_data, data = vd, negative.group.cells = negative.group.cells, positive.group.cells = positive.group.cells, negative.group.name = negative.group.name, positive.group.name = positive.group.name, features.exclude = features.exclude)

  return(list(plot = vp,
              data = vd,
              feat.plots = feat_plots,
              interactive.data = interactive.data))

  '  if (!missing(path.interactive.volcano.folder) && file.exists(file.path(dirname(wd), "shiny_apps", "single_cell_volcano", "app.R"))) {
    dir.create(path.interactive.volcano.folder, showWarnings = F, recursive = T)
    saveRDS(interactive.plot.data, file.path(path.interactive.volcano.folder, "data.rds"))
    file.copy(file.path(dirname(wd), "shiny_apps", "single_cell_volcano", "app.R"), file.path(path.interactive.volcano.folder, "app.R"), overwrite = T)
    print(paste0("Interactive volcano folder saved to ", path.interactive.volcano.folder))
  } else {
    print("For interactive volcano: Save interactive.plot.data as data.rds and make a folder with the app.R file from the scRNAseq base folder.")
  }'

}



.calc_vd <- function(an,
                     ap,
                     pgn = "positive",
                     ngn = "negative",
                     n.feat.for.p.adj = NULL,
                     inf.fc.shift = 2,
                     p.adjust = "bonferroni") {


  if (nrow(an) != nrow(ap)) {
    stop("an and ap need to have the same number of rows!")
  }
  p.adjust <- match.arg(p.adjust, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))

  p <- matrixTests::row_wilcoxon_twosample(as.matrix(an), as.matrix(ap))$pvalue
  apm <- apply(ap, 1, mean)
  anm <- apply(an, 1, mean)
  if (is.null(n.feat.for.p.adj)) {
    n.feat.for.p.adj <- nrow(an)
  }

  vd <- data.frame(Feature = rownames(an),
                   log2.fc = log2(apm / anm),
                   p.val = p,
                   adj.p.val = as.numeric(formatC(stats::p.adjust(p, method = p.adjust, n = n.feat.for.p.adj), format = "e", digits = 2)),
                   stats::setNames(list(round(anm, 2)), ngn),
                   stats::setNames(list(round(apm, 2)), pgn),
                   stats::setNames(list(round(apply(an, 1, function(x) sum(x != 0))/ncol(an), 2)), paste0("pct.", ngn)),
                   stats::setNames(list(round(apply(ap, 1, function(x) sum(x != 0))/ncol(ap), 2)), paste0("pct.", pgn)),
                   infinite.FC = ifelse(is.infinite(log2(apm / anm)), TRUE, FALSE))

  vd$log2.fc <- scales::oob_squish_infinite(vd$log2.fc, range = c(min(vd$log2.fc[!is.infinite(vd$log2.fc)], na.rm = T) - inf.fc.shift,
                                                                  max(vd$log2.fc[!is.infinite(vd$log2.fc)], na.rm = T) + inf.fc.shift))
  return(vd)
}



.plot_vp <- function (vd,
                      p.plot = "adj.p.val",
                      plot.symmetric = T,
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

  vp <- ggplot2::ggplot(vd, ggplot2::aes(x = log2.fc, y = round(scales::oob_squish_infinite(-log10(!!sym(p.plot)), range = c(0,300)), 2), label = Feature)) +
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
    vp <- vp + ggplot2::labs(x = bquote("avg" ~ log[2] ~ "FC"), y = bquote(-log[10]~.(sym(p.plot_label))))
  } else {
    vp <- vp + ggplot2::labs(x = bquote(bold(.(ngn)) ~ "  <====  " ~ log[2] ~ "FC" ~ "  ====>  " ~ bold(.(pgn))), y = bquote(-log[10]~.(sym(p.plot_label))))
  }

  if (plot.symmetric) {
    vp <- vp + ggplot2::xlim(-round(max(abs(vd$log2.fc))) - 0.5 - x.axis.extension, round(max(abs(vd$log2.fc))) + 0.5 + x.axis.extension)
  } else {
    vp <- vp + ggplot2::xlim(round(min(vd$log2.fc)) - 0.5 - x.axis.extension, round(max(vd$log2.fc)) + 0.5 + x.axis.extension)
  }

  if (!any(c(is.na(p.cut), is.na(fc.cut)))) {
    vp <- vp + ggplot2::geom_hline(yintercept = -log10(p.cut), linetype = "dashed") + ggplot2::geom_vline(xintercept = c(-fc.cut, fc.cut), linetype = "dashed")

    vd.cut <-
      vd %>%
      dplyr::filter(abs(log2.fc) >= fc.cut) %>%
      dplyr::filter(!!sym(p.plot) <= p.cut)

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
                      lables.topn = 30,
                      label.size = 4,
                      nudge.x = 0,
                      nudge.y = 0,
                      max.iter = 10000,
                      pos.neg.label.separate = T,
                      label.col = "black",
                      label.face = "bold",
                      font.family = "Courier",
                      max.overlaps = 50,
                      p.signif = 0.001,
                      features.exclude = NULL) {

  if (!is.null(features.exclude)) {
    vd <- vd[which(!grepl(paste0(features.exclude, collapse = "|"), vd$Feature)),]
  }

  if (is.null(label.features)) {
    if (topn.metric == "p.value") {
      f_lab <- vd %>% dplyr::top_n(-lables.topn, !!sym(p.plot))
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    } else if (topn.metric == "both") {
      f_lab.p.val <- vd %>% dplyr::top_n(-lables.topn, !!sym(p.plot))
      f_lab.logfc <- bind_rows(vd %>% dplyr::top_n(lables.topn/2, log2.fc), vd %>% dplyr::top_n(-(lables.topn/2), log2.fc))
      f_lab <- dplyr::bind_rows(f_lab.logfc, f_lab.p.val) %>% dplyr::distinct()
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    } else {
      f_lab <- dplyr::bind_rows(vd %>% dplyr::top_n(lables.topn/2, log2.fc), vd %>% dplyr::top_n(-(lables.topn/2), log2.fc))
      f_lab.pos <- f_lab %>% dplyr::filter(log2.fc > 0) %>% dplyr::pull(Feature)
      f_lab.neg <- f_lab %>% dplyr::filter(log2.fc < 0) %>% dplyr::pull(Feature)
    }
    vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% f_lab$Feature), colour = "tomato2")


    if (pos.neg.label.separate) {
      vp <- vp +
        ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% f_lab.pos), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
        ggrepel::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% f_lab.neg), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
    } else {
      vp <- vp +
        ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% f_lab$Feature), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
    }

  } else {
    if (length(label.features) == 1) {
      if (label.features == "significant") {
        label.features <- vd[which(as.numeric(vd[,p]) < p.signif), "Feature"]
      } else {
        vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = "tomato2")
        if (pos.neg.label.separate) {
          vp <- vp +
            ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc > 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
            ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc < 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
        } else {
          vp <- vp + ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
        }
      }
    } else {
      vp <- vp + ggplot2::geom_point(data = vd %>% dplyr::filter(Feature %in% label.features), colour = "tomato2")
      if (pos.neg.label.separate) {
        vp <- vp +
          ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc > 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps) +
          ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features & log2.fc < 0), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
      } else {
        vp <- vp + ggplot2::geom_text_repel(data = vd %>% dplyr::filter(Feature %in% label.features), family = font.family, color = label.col, fontface = label.face, size = label.size, max.iter = max.iter, nudge_x = -nudge.x, nudge_y = nudge.y, max.overlaps = max.overlaps)
      }
    }
  }
  return(vp)
}


