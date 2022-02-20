#' Title
#'
#' @param data.dir
#' @param nhvf
#' @param npcs
#' @param min_nCount_RNA
#' @param resolutions
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
qc_diagnostic <- function(data.dir, nhvf = 500, npcs = 10, min_nCount_RNA = 300, resolutions = seq(0.6,0.9,0.1),  ...) {
  if (length(data.dir) > 1) {
    message("Advice: Better provide only one data.dir, e.g. one sample only.")
  }
  counts <- as.matrix(Seurat::Read10X(data.dir = data.dir, strip.suffix = T))
  SO <-
    Seurat::CreateSeuratObject(counts = counts) %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures(nfeatures = nhvf) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = npcs) %>%
    Seurat::RunUMAP(dims = 1:npcs) %>%
    Seurat::FindNeighbors(dims = 1:npcs) %>%
    Seurat::FindClusters(algorithm = 1, resolution = resolutions)

  if (any(matrixStats::colSums2(counts) < min_nCount_RNA)) {
    warning(paste0("Transcriptomes (cells) with less than ", min_nCount_RNA, " total transcripts found. These will be excluded from finding doublets. They will have NA as dbl_score."))
  }
  dbl_score <- setNames(scDblFinder::computeDoubletDensity(x = counts[,which(matrixStats::colSums2(counts) >= min_nCount_RNA)],
                                                           subset.row = Seurat::VariableFeatures(SO),
                                                           dims = npcs),
                        nm = colnames(counts[,which(matrixStats::colSums2(counts) >= min_nCount_RNA)]))

  SO <- Seurat::AddMetaData(SO, dbl_score, "dbl_score")
  SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^MT-"), "pct_mt")
  for (i in c(c("nCount_RNA", "nFeature_RNA", "pct_mt", "dbl_score"))) {
    SO@meta.data[,paste0(i, "_log")] <- log1p(SO@meta.data[,i])
  }
  lm <- stats::lm(nFeature_RNA_log~nCount_RNA_log, data = SO@meta.data)
  SO@meta.data$residuals <- stats::residuals(lm)
  SO@meta.data$residuals_norm <- SO@meta.data$residuals/SO@meta.data$nCount_RNA_log/SO@meta.data$nFeature_RNA_log


  qc_p1 <- scexpr::feature_plot(SO, features = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log", "residuals_norm", paste0("RNA_snn_res.", resolutions[length(resolutions)])),
                                reduction = "UMAP", legend.position = c(0,1), plot.labels = "text", label.size = 10)

  qc_p2 <- cowplot::plot_grid(scexpr::qc_params_meta_cols(SO,
                                                          qc.cols = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log", "residuals_norm"),
                                                          meta.cols = paste0("RNA_snn_res.", resolutions[length(resolutions)])),
                              ggplot2::ggplot(SO@meta.data, ggplot2::aes(nCount_RNA_log, nFeature_RNA_log, color = !!rlang::sym(paste0("RNA_snn_res.", resolutions[length(resolutions)])))) +
                                ggplot2::geom_point() +
                                ggplot2::theme_bw() +
                                ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank()) +
                                ggplot2::scale_color_manual(values = scexpr::col_pal()) +
                                ggplot2::labs(color = ""),
                              ggplot2::ggplot(dplyr::arrange(SO@meta.data, dbl_score_log), ggplot2::aes(nCount_RNA_log, nFeature_RNA_log, color = dbl_score_log)) +
                                ggplot2::geom_point() +
                                ggplot2::theme_bw() +
                                ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank()) +
                                ggplot2::scale_color_gradientn(colors = scexpr::col_pal("spectral")),
                              ggplot2::ggplot(dplyr::arrange(SO@meta.data, pct_mt_log), ggplot2::aes(nCount_RNA_log, nFeature_RNA_log, color = pct_mt_log)) +
                                ggplot2::geom_point() +
                                ggplot2::theme_bw() +
                                ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank()) +
                                ggplot2::scale_color_gradientn(colors = scexpr::col_pal("spectral")),
                              align = "hv", axis = "tblr")

  qc_p3 <- cowplot::plot_grid(scexpr::feature_plot(SO,
                                                   features = paste0("RNA_snn_res.", resolutions),
                                                   reduction = "UMAP", legend.position = "none", plot.labels = "text", label.size = 10, nrow.combine = 1, plot.title = F),
                              scexpr::qc_params_meta_cols(SO,
                                                          qc.cols = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log"),
                                                          meta.cols = paste0("RNA_snn_res.", resolutions)),
                              nrow = 2,
                              rel_heights = c(0.3,0.7),
                              align = "hv", axis = "tblr")

  return(list(SO = SO, qc_p1 = qc_p1, qc_p2 = qc_p2))

}
