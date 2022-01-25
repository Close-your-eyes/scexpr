feature_correlation <- function(SO,
                                assay = c("RNA", "SCT"),
                                features,
                                cells = NULL,
                                min_pct = 0.1,
                                ...) {


  SO <- .check.SO(SO = SO, assay = assay, split.by = NULL, shape.by = NULL, max.length = 1)
  features <- .check.features(SO = SO, features = unique(features), meta.data = F)
  cells <- .check.and.get.cells(SO = SO, assay = assay, cells = cells)

  ref_mat <- as.matrix(Seurat::GetAssayData(SO, assay = assay)[which(Matrix::rowSums(Seurat::GetAssayData(SO, assay = assay)) > 0), names(cells[which(cells == 1)])])
  ref_mat <- ref_mat[Matrix::rowSums(ref_mat > 0)/ncol(ref_mat) >= min_pct, ]


  corr_obj <- psych::corr.test(t(as.matrix(Seurat::GetAssayData(SO, assay = assay)[features,])), t(ref_mat), ci = F) #...
  corr_df <- merge(reshape2::melt(t(corr_obj[["r"]]), value.name = "r"), reshape2::melt(t(corr_obj[["p"]]), value.name = "p"))
  corr_df$m.log10.p <- -log10(corr_df$p)
  names(corr_df)[1:2] <- c("ref_feature", "features")

  xx <- reshape2::melt(t(corr_obj[["r"]]), value.name = "r")

  ggplot(corr_df, aes(x = r, y = m.log10.p)) +
    geom_point() +
    theme_bw() +
    facet_wrap(vars(features))

}
