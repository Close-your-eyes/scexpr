#' Title
#'
#' @param SO
#' @param assay
#' @param features
#' @param cells
#' @param min_pct
#' @param limit_p
#' @param lm_resid
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
feature_correlation <- function(SO,
                                assay = c("RNA", "SCT"),
                                features,
                                cells = NULL,
                                min_pct = 0.1,
                                limit_p = 1e-303,
                                lm_resid = F,
                                ...) {

  if (missing(features)) {
    stop("Please provide features.")
  }
  assay <- match.arg(assay, c("RNA", "SCT"))
  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  features <- .check.features(SO = SO, features = unique(features), meta.data = F)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells)
  cells <- names(cells[which(cells == 1)])

  ref_mat <- as.matrix(Seurat::GetAssayData(SO, assay = assay)[filter_feature(SO = SO, assay = assay, min_pct = min_pct, cells = cells), cells, drop=F])
  mat <- as.matrix(Seurat::GetAssayData(SO, assay = assay)[features, cells, drop=F])


  corr_obj <- psych::corr.test(t(mat), t(ref_mat), ci = F) #...
  corr_df <- merge(merge(reshape2::melt(t(corr_obj[["r"]]), value.name = "r"), reshape2::melt(t(corr_obj[["p"]]), value.name = "p")), reshape2::melt(t(corr_obj[["p.adj"]]), value.name = "p.adj"))

  if (is.numeric(limit_p)) {
    corr_df$p.adj[which(corr_df$p.adj == 0)] <- limit_p
    corr_df$p[which(corr_df$p == 0)] <- limit_p
  }

  corr_df$minus.log10.p <- -log10(corr_df$p)
  corr_df$minus.log10.p.adj <- -log10(corr_df$p.adj)
  names(corr_df)[1:2] <- c("ref_feature", "feature")
  ## add pcts
  corr_df <- corr_df %>% dplyr::left_join(stats::setNames(utils::stack(pct_feature(SO, assay = assay, features = unique(corr_df$ref_feature))), c("ref_feature_pct", "ref_feature")), by = "ref_feature")
  corr_df <- corr_df %>% dplyr::left_join(stats::setNames(utils::stack(pct_feature(SO, assay = assay, features = unique(corr_df$feature))), c("feature_pct", "feature")), by = "feature")

  return(corr_df)


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
    dplyr::group_by(feature, ref_feature) %>%
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
    scale_color_gradientn(colors = scexpr::col_pal("spectral"))

'

}



filter_feature <- function(SO,
                           assay = c("RNA", "SCT"),
                           cells = NULL,
                           min_pct = 0.1) {

  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells)
  cells <- names(cells[which(cells == 1)])
  assay <- match.arg(assay, c("RNA", "SCT"))

  # filter non expressed features first
  f1 <- names(which(Matrix::rowSums(Seurat::GetAssayData(SO, assay = assay)[,cells]) > 0))
  # get those above min_pct
  f2 <- names(which(pct_feature(SO = SO, assay = assay, features = f1, cells = cells) >= min_pct))

  return(f2)
}


pct_feature <- function(SO,
                        features,
                        cells = NULL,
                        assay = c("RNA", "SCT")) {

  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, length = 1)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells)
  cells <- names(cells[which(cells == 1)])
  #features <- .check.features(SO = SO, features = unique(features), meta.data = F) # need speed up first
  assay <- match.arg(assay, c("RNA", "SCT"))

  return(Matrix::rowSums(Seurat::GetAssayData(SO, assay = assay)[features,cells,drop=F] > 0)/length(cells))
}
