#' Get features whose expression is correlated or anti-correlated across all cells or groups of cells
#'
#' Find out the expression of which features is correlated or anti-correlated using simple correlation metrics like pearson or spearman.
#' The analysis may be applied to a subset of cells or subset of features (see arguments). Due to dropouts in some scRNAseq technologies
#' this analysis is not super-clean but may still provide a valid, relative, comparison of feature correlations.
#' Other methods considering the dropout-effect exist: \href{"https://academic.oup.com/nargab/article/3/3/lqab072/6348150"}{COTAN}.
#'
#'
#' @param SO Seurat object
#' @param features which features to calculate correlations for (must be rownames in the selected assay)
#' @param assay which assay to obtain expression values from; the data slot will be used in any case
#' @param cells vector of cell names to consider for correlation anaylsis; if NULL (default) all cells are used
#' @param min.pct minimum percentage of expressing cells (> 0 UMIs) to include a feature in correlation analysis
#' @param limit_p p-value which p-values of 0 will be set to; this avoids obtaining INF when deriving -log10(p-val)
#' @param method which metric of correlation to calculate
#' @param bar.fill which bar fill to apply
#' @param theme which ggplot theme to set as basis
#' @param split.by groups for correlation analysis; must be a categorical column in meta.data of SO
#' @param min.group.size required number of cell in one group to be considered in analysis
#' @param topn numeric vector of length two; number of top anti-correlated and correlated features, respectively
#' @param ... additional arguments to psych::corr.test and to ggplot2::theme
#'
#' @return a default plot (ggplot2 object) and the underlying data frame with correlation values
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' }
feature_correlation <- function(SO,
                                features,
                                assay = c("RNA", "SCT"),
                                method = c("pearson", "spearman", "kendall"),
                                cells = NULL,
                                min.pct = 0.1,
                                limit_p = 1e-303,
                                bar.fill = c("correlation_sign", "ref_feature_pct", "none"),
                                theme = ggplot2::theme_bw(),
                                split.by = NULL,
                                min.group.size = 20,
                                topn = c(10,10), # n for min and max
                                ...) {

  if (!requireNamespace("psych", quietly = T)) {
    utils::install.packages("psych")
  }
  if (!requireNamespace("reshape2", quietly = T)) {
    utils::install.packages("reshape2")
  }
  if (!requireNamespace("Matrix", quietly = T)) {
    utils::install.packages("Matrix")
  }

  if (missing(features)) {
    stop("Please provide features.")
  }
  dots <- list(...)
  assay <- match.arg(assay, c("RNA", "SCT"))
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  bar.fill <- match.arg(bar.fill, c("correlation_sign", "ref_feature_pct", "none"))

  if (length(topn) != 2 || !is.numeric(topn)) {
    warning("topn should be a numeric vector of length 2, providing the numbers of lowest ranked features (index 1) and highest ranked features (index 2) with respect to correlation. Will be set to default c(10,10)")
    topn <- c(10,10)
  }

  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  features <- .check.features(SO = SO, features = unique(features), meta.data = F, meta.data.numeric = T)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells, return.included.cells.only = T)

  ref_mat <- as.matrix(Seurat::GetAssayData(SO, assay = assay)[filter_feature(SO = SO, assay = assay, min.pct = min.pct, cells = cells), cells, drop=F])

  # putative dichotomous meta features which are TRUE / FALSE: they are coerced to 1 / 0
  # applying pearson correlation of a 0/1 dichotomous variable and a numeric one is called point-biseral correlation
  mat <- rbind(as.matrix(Seurat::GetAssayData(SO, assay = assay)[features[which(features %in% rownames(SO))], cells, drop=F]),
                         as.matrix(t(SO@meta.data[cells,features[which(features %in% names(SO@meta.data))], drop = F])))

  if (!"ci" %in% names(dots)) {
    ci <- F
  } else {
    ci <- dots[["ci"]]
    dots <- dots[-which(names(dots) == "ci")]
  }

  ## do all the checking
  if (!is.null(split.by)) {
    groups <-
      SO@meta.data[ ,split.by, drop=F] %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::filter(ID %in% cells)
    groups <- split(groups$ID, groups[,split.by,drop=F])
  } else {
    groups <- stats::setNames(list(cells), "all")
  }

  if (any(lengths(groups) < min.group.size)) {
    print(paste0("These groups are not considered due to having less cells than min.group.size: ", paste(names(which(lengths(groups) < min.group.size)), collapse = ", ")))
    groups <- groups[which(lengths(groups) >= min.group.size)]
  }


  # NAs are removed by psych::corr.test
  out <- lapply(names(groups), function(x) {
    corr_obj <- do.call(psych::corr.test, args = c(list(x = t(mat[,groups[[x]],drop=F]),
                                                        y = t(ref_mat[,groups[[x]],drop=F]),
                                                        ci = ci,
                                                        method = method),
                                                   dots[which(names(dots) %in% names(formals(psych::corr.test)))]))

    corr_df <- merge(merge(reshape2::melt(t(corr_obj[["r"]]), value.name = "r"),
                           reshape2::melt(t(corr_obj[["p"]]), value.name = "p")),
                     reshape2::melt(t(corr_obj[["p.adj"]]), value.name = "p.adj"))

    if (is.numeric(limit_p)) {
      corr_df$p.adj[which(corr_df$p.adj == 0)] <- limit_p
      corr_df$p[which(corr_df$p == 0)] <- limit_p
    }

    corr_df$minus.log10.p <- -log10(corr_df$p)
    corr_df$minus.log10.p.adj <- -log10(corr_df$p.adj)
    names(corr_df)[1:2] <- c("ref_feature", "feature")
    ## add pcts
    corr_df <- dplyr::left_join(corr_df, stats::setNames(utils::stack(pct_feature(SO, assay = assay, features = unique(corr_df$ref_feature))), c("ref_feature_pct", "ref_feature")), by = "ref_feature")
    ## handle numeric meta col which do not have pct expression
    pctx <- pct_feature(SO, assay = assay, features = unique(corr_df$feature))
    if (length(pctx) > 0) {
      corr_df <- dplyr::left_join(corr_df, stats::setNames(utils::stack(pctx), c("feature_pct", "feature")), by = "feature")
    }

    corr_df <- dplyr::mutate(corr_df, feature = as.character(feature), ref_feature = as.character(ref_feature))
    corr_df_plot <- dplyr::filter(corr_df, feature != ref_feature)

    corr_df_plot <- rbind(corr_df_plot %>%
                            dplyr::group_by(feature) %>%
                            dplyr::slice_min(n = topn[1], order_by = r) %>%
                            dplyr::ungroup(),
                          corr_df_plot %>%
                            dplyr::group_by(feature) %>%
                            dplyr::slice_max(n = topn[2], order_by = r) %>%
                            dplyr::ungroup())
    corr_df_plot <- dplyr::mutate(corr_df_plot, correlation_sign = factor(ifelse(r > 0, "+", "-"), levels = c("+", "-")))
    corr_df_plot[,"group"] <- x
    corr_df_plot[,"feature_group"] <- paste0(corr_df_plot[,"feature", drop = T], "_", x)
    corr_df[,"group"] <- x
    corr_df[,"feature_group"] <- paste0(corr_df[,"feature"], "_", x)

    return(list(corr_df, corr_df_plot))
  })

  corr_df <- do.call(rbind, sapply(out, "[", 1))
  corr_df_plot <- do.call(rbind, sapply(out, "[", 2))

  if (is.null(split.by) && length(features) == 1) {
    plot <- ggplot2::ggplot(corr_df_plot, ggplot2::aes(x = r, y = stats::reorder(ref_feature, r), fill = !!rlang::sym(bar.fill))) +
      ggplot2::geom_bar(stat = "identity", color = "black") +
      theme +
      do.call(ggplot2::theme, args = dots[which(names(dots) %in% names(formals(ggplot2::theme)))])
  } else {
    plot <- ggplot2::ggplot(corr_df_plot, ggplot2::aes(x = r, y = .reorder_within(ref_feature, r, feature_group), fill = !!rlang::sym(bar.fill))) +
      ggplot2::geom_bar(stat = "identity", color = "black") +
      theme +
      .scale_y_reordered() +
      do.call(ggplot2::theme, args = dots[which(names(dots) %in% names(formals(ggplot2::theme)))])
  }

  wrap_by <- function(...) {ggplot2::facet_wrap(ggplot2::vars(...), scales = "free")}
  if (is.null(split.by)) {
    plot <- plot + wrap_by(feature)
  } else {
    plot <- plot + wrap_by(feature, group)
  }

  #https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
  #https://stackoverflow.com/questions/5293715/how-to-use-greek-symbols-in-ggplot2
  xlabel <- switch(method,
                   "pearson" = "pearson r",
                   "spearman" = paste0("spearman ", "\u03C1"),
                   "kendall" = paste0("kendall ", "\u03C4 "))

  if (bar.fill == "correlation_sign") {
    plot <- plot + ggplot2::scale_fill_manual(values = c("forestgreen", "tomato2")) + ggplot2::labs(y = "Feature", x = xlabel, fill = "direction")
  } else if (bar.fill == "ref_feature_pct") {
    plot <- plot + ggplot2::scale_fill_viridis_c() + ggplot2::labs(y = "Feature", x = xlabel, fill = "pct")
  }

  return(list(plot = plot, data = corr_df))

  # lm dist or dist from corr-line (which goes through 0 by definition, so only has slope)
  # only for pearson
  # https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line#35209406

  # for spearman:
  # get rank diff by observation
  # https://www.simplilearn.com/tutorials/statistics-tutorial/spearmans-rank-correlation

  '  ggplot(corr_df, aes(x = r, y = minus.log10.p.adj)) +
    geom_point() +
    theme_bw() +
    facet_wrap(vars(feature))'



  # many models
  # https://r4ds.had.co.nz/many-models.html
  '  mat2 <- reshape2::melt(mat)
  ref_mat2 <- reshape2::melt(ref_mat)
  groups <-
    corr_df %>%
    dplyr::split.by(feature, ref_feature) %>%
    tidyr::nest()'


  '  df <- as.data.frame(cbind(t(mat[1,,drop=F]),t(ref_mat[which(rownames(ref_mat) == "B2M"),,drop=F])))
  model <- lm(CD8B ~ B2M, data = df)
  broom::tidy(model)
  plot(df)
  model <- lm(B2M ~ CD8B, data = df)
  df <- modelr::add_residuals(df, model = model, var = "resid")
  ggplot(df, aes(CD8B,B2M, color = abs(resid))) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_color_gradientn(colors = colrr::col_pal("spectral"))

'

}



filter_feature <- function(SO,
                           assay = c("RNA", "SCT"),
                           cells = NULL,
                           min.pct = 0.1) {

  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells, return.included.cells.only = T)
  assay <- match.arg(assay, c("RNA", "SCT"))

  # filter non expressed features first
  f1 <- names(which(Matrix::rowSums(Seurat::GetAssayData(SO, assay = assay)[,cells]) > 0))
  # get those above min.pct
  f2 <- names(which(pct_feature(SO = SO, assay = assay, features = f1, cells = cells) >= min.pct))

  return(f2)
}


pct_feature <- function(SO,
                        features,
                        cells = NULL,
                        assay = c("RNA", "SCT")) {

  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells, return.included.cells.only = T)
  #features <- .check.features(SO = SO, features = unique(features), meta.data = F) # need speed up first
  assay <- match.arg(assay, c("RNA", "SCT"))

  ## case of dichtomous meta.col
  features <- features[which(features %in% rownames(Seurat::GetAssayData(SO, assay = assay)))]

  return(Matrix::rowSums(Seurat::GetAssayData(SO, assay = assay)[features[which(features %in% rownames(Seurat::GetAssayData(SO, assay = assay)))], cells, drop=F] > 0)/length(cells))
}

.reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

.scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

.scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

# try this:
# http://www.bioconductor.org/packages/release/bioc/html/scHOT.html

# or this:
# https://github.com/edvanburen/twosigma
# https://academic.oup.com/bib/article-abstract/23/3/bbac084/6553609?login=true
# https://www.biorxiv.org/content/10.1101/2021.01.24.427979v1

# from Twitter:
# https://divingintogeneticsandgenomics.rbind.io/post/how-to-do-gene-correlation-for-single-cell-rnaseq-data-part-1/


