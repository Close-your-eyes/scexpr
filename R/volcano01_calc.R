#' Make a volcano plot comparing the gene expression of selected groups of cells from one or more preprocessed Seurat objects
#'
#'
#'
#' @param SO one or more Seurat object(s); provide multiple objects as a (named) list
#' @param assay which assay to get expression data from
#' @param neg_cells vector of cell names for the negative group of cells; genes which
#' are expressed at higher level in this group will have a negative sign (minus); use Seurat::Cells()
#' or Seurat::WhichCells() to select or use filter operations on SO@meta.data to select cells
#' @param pos_cells vector of cell names for the positive group of cells; genes which
#' are expressed at higher level in this group will have a positive sign (plus); use Seurat::Cells()
#' or Seurat::WhichCells() to select or use filter operations on SO@meta.data to select cells
#' @param neg_name name for the negative group of cells; this will appear on plot axes
#' @param pos_name name for the positive group of cells; this will appear on plot axes
#' @param fc_inf_squish numeric; genes with an infinite fold-change (happens when a gene not expressed at all in one group)
#' are not left infinite but will be put at: "largest finite fold-change +/- fc_inf_squish"; such genes are indicated
#' specifically
#' @param min_pct numeric; only consider genes for the volcano plot which are expressed at least by this
#' fraction of cells at least one group (neg_cells or pos_cells)
#' @param p_adj method for p-value adjustment
#' @param method backend for DE calculation;limma:intended for subsequent use MetaVolcanoR
#' which relies on values returned from limma; MAST: Seurat::FindMarkers;
#' with MAST as method; default: matrixTests::row_wilcoxon_twosample
#'
#' @return
#' @export
#'
#' @examples
volcano01_calc <- function(SO,
                           assay = "RNA",
                           layer = "data",
                           neg_cells,
                           pos_cells,
                           neg_name = "negative.group",
                           pos_name = "positive.group",
                           method = c("MAST", "wilcox", "wilcox_limma", "bimod", "roc", "t", "negbinom", "poisson", "LR", "DESeq2", "limma", "default"),
                           fc_inf_squish = 2,
                           p_zero_squish = 1,
                           min_pct = 0.1,
                           p_adj = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"),
                           fc_thresh = 0,
                           feature_plot_args = list(col_pal_d_args = list(name = stats::setNames(colrr::col_pal("custom")[2:3],
                                                                                                 c(neg_name, pos_name)),
                                                                          missing_fct_to_na = F))) {

  if (!requireNamespace("colrr", quietly = T)) {
    devtools::install_github("Close-your-eyes/colrr")
  }

  default_warn <- getOption("warn")
  options(warn = 1)

  if (missing(neg_cells) || missing(pos_cells)) {
    stop("pos_cells and neg_cells are required.")
  }
  if (!length(neg_cells) || !length(pos_cells)) {
    stop("Length of neg_cells or pos_cells is zero.")
  }
  if (length(intersect(neg_cells, pos_cells)) > 0) {
    stop("Some neg_cells and pos_cells intersect. Fix this.")
  }


  SO <- check.SO(SO = SO, assay = assay)
  assay <- Seurat::DefaultAssay(SO[[1]])
  min_pct <- max(0, min_pct)
  p_adj <- rlang::arg_match(p_adj)
  method <- rlang::arg_match(method)

  pgc <- pos_cells
  ngc <- neg_cells

  intersect_features <- Reduce(intersect, lapply(SO, rownames))
  nonintersect_features <- lapply(SO, function(x) setdiff(rownames(x), intersect_features))
  if (length(unique(c(sapply(SO, nrow), length(intersect_features)))) != 1) {
    message("Different features across SOs detected. Will use intersecting ones only.")
  }

  colname <- "volcano_groups"
  SO <- lapply(SO, function(x) {
    # addmetadata?
    x@meta.data[,colname] <- ifelse(!rownames(x@meta.data) %in% c(pgc, ngc), "other",
                                    ifelse(rownames(x@meta.data) %in% pgc,
                                           pos_name[1],
                                           neg_name[1]))
    return(x)
  })

  feat_plots <- tryCatch(
    expr = {
      do.call(feature_plot2,
              args = c(list(SO = SO, features = colname),
                       feature_plot_args))
    }, error = function(e) NULL)



  # DietSeurat is slow ?!
  # Seurat::DietSeurat(, assays = assay, counts = F)
  SO <- lapply(SO, function(x) Seurat::DietSeurat(subset(
    x,
    features = intersect(rownames(x), intersect_features),
    cells = intersect(colnames(x), c(ngc, pgc))
  ), assays = assay, counts = F))

  if (length(SO) > 1) {
    # this restores the counts matrix
    ## fix this like in so_prep
    SO <- merge(x = SO[[1]], y = SO[2:length(SO)], merge.data = T)
  } else {
    SO <- SO[[1]]
  }

  # exclude features which are below min_pct.set in both populations
  # = only continue with those which are ">= min_pct" in at least one of them
  if (min_pct > 0) {
    comparefun <- `>=`
  } else if (min_pct == 0) {
    comparefun <- `>`
  }

  min_pct_features <- unique(c(names(which(comparefun(Matrix::rowSums(get_layer(obj = SO, assay = assay, layer = layer, cells = ngc) != 0)/length(ngc), min_pct))),
                               names(which(comparefun(Matrix::rowSums(get_layer(obj = SO, assay = assay, layer = layer, cells = pgc) != 0)/length(pgc), min_pct)))))

  min_pct_features_removed <- setdiff(intersect_features, min_pct_features)
  SO <- Seurat::DietSeurat(SO, assays = assay, features = min_pct_features, counts = F)

  # equal order of intersecting features which are taken into non-log space for wilcox test and FC calculation
  # DefaultAssay set above
  vd <- calculate_DEG(SO = SO,
                      layer = layer,
                      ngc = ngc,
                      pgc = pgc,
                      pgn = pos_name[1],
                      ngn = neg_name[1],
                      fc_inf_squish = fc_inf_squish,
                      n.feat.for.p.adj = length(intersect_features),
                      p_adj = p_adj,
                      method = method,
                      fc_thresh = fc_thresh)

  df <- tibble::rownames_to_column(as.data.frame(vd), "feature")
  attr(df, "pos_name") <- pos_name[1]
  attr(df, "neg_name") <- neg_name[1]
  attr(df, "method") <- method
  attr(df, "p_adj") <- p_adj
  attr(df, "min_pct") <- min_pct
  attr(df, "assay") <- assay

  options(warn = default_warn)
  return(list(df = df,
              mat = vd,
              feat.plots = feat_plots,
              intersect_features = intersect_features,
              nonintersect_features = nonintersect_features,
              min_pct_features = min_pct_features,
              min_pct_features_removed = min_pct_features_removed))
}



calculate_DEG <- function(SO,
                          ngc,
                          pgc,
                          layer = "data",
                          pgn = "positive",
                          ngn = "negative",
                          n.feat.for.p.adj = NULL,
                          fc_inf_squish = 2,
                          p_zero_squish = 1,
                          p_adj = "bonferroni",
                          method = "MAST",
                          fc_thresh = 0) {


  if (method %in% c("MAST", "wilcox", "wilcox_limma", "bimod", "roc", "t", "negbinom", "poisson", "LR", "DESeq2")) {
    if (!requireNamespace("MAST", quietly = T)) {
      BiocManager::install("MAST")
    }
    Seurat::Idents(SO) <- SO@meta.data[,1,drop=T]
    df <- Seurat::FindMarkers(Seurat::GetAssay(SO, Seurat::DefaultAssay(SO)),
                              test.use = method,
                              min.pct = 0, # filtered before
                              logfc.threshold = fc_thresh,
                              cells.1 = pgc, cells.2 = ngc)

    names(df)[which(names(df) == "avg_log2FC")] <- "log2.fc"
    names(df)[which(names(df) == "p_val")] <- "p.val"
    names(df)[which(names(df) == "p_val_adj")] <- "adj.p.val"
    names(df)[which(names(df) == "pct.1")] <- paste0("pct.", pgn)
    names(df)[which(names(df) == "pct.2")] <- paste0("pct.", ngn)

    SO <- get_layer(obj = SO, layer = layer)
    SO <- expm1(SO) + 1 # actually works only for data slot
    if (layer != "data") {
      message("calculate_DEG: layer != 'data'. log expression values in returned df may be imprecise.")
    }
    apm <- Matrix::rowMeans(SO[, pgc])
    anm <- Matrix::rowMeans(SO[, ngc])

    df[[ngn]] <- round(log2(anm[rownames(df)]), 2)
    df[[pgn]] <- round(log2(apm[rownames(df)]), 2)

  }



  if (method == "limma") {
    if (!requireNamespace("BiocManager", quietly = T)) {
      utils::install.packages("BiocManager")
    }
    if (!requireNamespace("limma", quietly = T)) {
      BiocManager::install("limma")
    }

    SO <- get_layer(obj = SO, layer = layer)

    message("(i) limma: p_adj ignored. (ii) caution: limma calculates log2.fc differently as Seurat, see here: https://support.bioconductor.org/p/82478/ ")
    df <- limma::topTable(limma::eBayes(limma::lmFit(
      SO[, c(ngc,pgc)],
      design = stats::model.matrix(~c(colnames(SO[, c(ngc,pgc)]) %in% pgc))
    )), number = nrow(SO), confint = T)

    SO <- expm1(SO) + 1
    apm <- Matrix::rowMeans(SO[, pgc])
    anm <- Matrix::rowMeans(SO[, ngc])

    df <- data.frame(log2.fc = df$logFC,
                     CI.L = df$CI.L,
                     CI.R = df$CI.R,
                     p.val = df$P.Value,
                     adj.p.val = as.numeric(formatC(df$adj.P.Val, format = "e", digits = 2)),
                     stats::setNames(list(round(log2(anm[rownames(df)]), 2)), ngn),
                     stats::setNames(list(round(log2(apm[rownames(df)]), 2)), pgn),
                     stats::setNames(list(round(apply(SO[, ngc]-1, 1, function(x) sum(x != 0))/ncol(SO[, ngc]), 2)[rownames(df)]), paste0("pct.", ngn)),
                     stats::setNames(list(round(apply(SO[, pgc]-1, 1, function(x) sum(x != 0))/ncol(SO[, pgc]), 2)[rownames(df)]), paste0("pct.", pgn)),
                     #infinite.FC = ifelse(is.infinite(log2(apm) - log2(anm)), 1, 0)[rownames(df)],
                     check.names = F)

  }

  if (method == "default") {
    if (!requireNamespace("matrixTests", quietly = T)) {
      utils::install.packages("matrixTests")
    }

    SO <- get_layer(obj = SO, layer = layer)

    SO <- expm1(SO) + 1
    p <- matrixTests::row_wilcoxon_twosample(as.matrix(SO[, ngc]), as.matrix(SO[, pgc]))$pvalue
    apm <- Matrix::rowMeans(SO[, pgc])
    anm <- Matrix::rowMeans(SO[, ngc])

    if (is.null(n.feat.for.p.adj)) {
      n.feat.for.p.adj <- nrow(SO)
    }

    df <- data.frame(log2.fc = log2(apm) - log2(anm), #log2(apm / anm),
                     p.val = p,
                     adj.p.val = as.numeric(formatC(stats::p.adjust(p, method = p_adj, n = n.feat.for.p.adj), format = "e", digits = 2)),
                     stats::setNames(list(round(log2(anm), 2)), ngn),
                     stats::setNames(list(round(log2(apm), 2)), pgn),
                     stats::setNames(list(round(apply(SO[, ngc]-1, 1, function(x) sum(x != 0))/ncol(SO[, ngc]), 2)), paste0("pct.", ngn)),
                     stats::setNames(list(round(apply(SO[, pgc]-1, 1, function(x) sum(x != 0))/ncol(SO[, pgc]), 2)), paste0("pct.", pgn)),
                     #infinite.FC = ifelse(is.infinite(log2(apm) - log2(anm)), 1, 0),
                     check.names = F)
  }

  df$infinite.FC <- ifelse(is.infinite(df$log2.fc), 1, 0)
  if (!is.null(fc_inf_squish)) {
    df$log2.fc <- scales::oob_squish_infinite(df$log2.fc, range = c(min(df$log2.fc[!is.infinite(df$log2.fc)], na.rm = T) - fc_inf_squish,
                                                                    max(df$log2.fc[!is.infinite(df$log2.fc)], na.rm = T) + fc_inf_squish))
  }

  df$zero.p <- ifelse(df$adj.p.val == 0, 1, 0)
  if (!is.null(p_zero_squish)) {
    # p_zero_squish log units
    df$adj.p.val <- scales::oob_squish(df$adj.p.val, range = c(min(df$adj.p.val)/10^p_zero_squish, 1))
  }


  return(as.matrix(df))
}
