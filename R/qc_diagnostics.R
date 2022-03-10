#' Title
#'
#' @param data.dir list or vector of full paths to filtered_feature_bc_matrix and raw_feature_bc_matrix, or filtered_feature_bc_matrix only (prohibiting the use of SoupX);
#' if not named, then the name of the common parent folder is uses as name(s)
#' @param nhvf
#' @param npcs
#' @param min_nCount_RNA
#' @param resolutions
#' @param SoupX logical whether to run SoupX. If T, raw_feature_bc_matrix is needed.
#' @param SoupX.resolution
#' @param ...
#'
#' On error in scDblFinder: Increase nhvf or try to change any other parameter
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
qc_diagnostic <- function(data.dir, nhvf = 1000, npcs = 10, min_nCount_RNA = 300, resolutions = seq(0.4,0.9,0.1), SoupX = F, SoupX.resolution = 0.6, ...) {

  # set ... for computeDoubletDensity
  # set ... for autoEstCont
  # decontX as alternative if raw_feature_bc_matrix is missing (https://github.com/campbio/celda)

  # check data.dir
  if (class(data.dir) == "list") {
    if (length(data.dir) > 1) {
      stop("Only provide one data.dir which may contain two folders, filtered_feature_bc_matrix and raw_feature_bc_matrix.")
    }
    if (!is.null(names(data.dir))) {
      name <- names(data.dir)
    } else {
      if (length(unique(dirname(data.dir[[1]]))) != 1) {
        stop("data.dirs should have common parent dir, otherwise provide a named list of filtered_feature_bc_matrix and raw_feature_bc_matrix.")
      }
      name <- basename(unique(dirname(data.dir[[1]])))
    }
    data.dir <- stats::setNames(data.dir[[1]], nm = rep(name, length(data.dir[[1]])))
  } else {
    if (is.null(names(data.dir))) {
      if (length(unique(dirname(data.dir))) != 1) {
        stop("data.dirs should have common parent dir, otherwise provide a named list of filtered_feature_bc_matrix and raw_feature_bc_matrix.")
      }
      names(data.dir) <- rep(basename(unique(dirname(data.dir))), 2)
    }
  }

  if (length(data.dir) == 1) {
    if (basename(data.dir) != "filtered_feature_bc_matrix") {
      stop("data.dir is expected to contain full paths of two dirs, filtered_feature_bc_matrix and raw_feature_bc_matrix, or filtered_feature_bc_matrix only.")
    }
  } else if (length(data.dir) == 2) {
    if (length(intersect(basename(data.dir), c("filtered_feature_bc_matrix", "raw_feature_bc_matrix"))) < 2) {
      stop("data.dir is expected to contain full paths of two dirs, filtered_feature_bc_matrix and raw_feature_bc_matrix, or filtered_feature_bc_matrix only.")
    }
  }

  if (any(basename(data.dir) == data.dir)) {
    stop("Please provide full path to filtered_feature_bc_matrix and/or raw_feature_bc_matrix.")
  }


  if (SoupX && length(data.dir) != 2) {
    warning("SoupX requires filtered_feature_bc_matrix and raw_feature_bc_matrix. SoupX will not run.")
    SoupX <- F
  }

  if (SoupX && !any(grepl(SoupX.resolution, resolutions))) {
    stop("SoupX.resolution not found in resolutions.")
  }


  filt_data <- Seurat::Read10X(data.dir = grep("filtered_feature_bc_matrix", data.dir, value = T))
  if (is.list(filt_data)) {
    filt_data <- filt_data[["Gene Expression"]]
  }

  SO <-
    Seurat::CreateSeuratObject(counts = as.matrix(filt_data)) %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures(nfeatures = nhvf) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = npcs) %>%
    Seurat::RunUMAP(dims = 1:npcs) %>%
    Seurat::FindNeighbors(dims = 1:npcs) %>%
    Seurat::FindClusters(algorithm = 1, resolution = resolutions)

  if (SoupX) {
    raw_data <- Seurat::Read10X(data.dir = grep("raw_feature_bc_matrix", data.dir, value = T))
    if (is.list(raw_data)) {
      raw_data <- raw_data[["Gene Expression"]]
    }
    sc <- SoupX::SoupChannel(tod = raw_data,
                             toc = filt_data)
    sc = SoupX::setClusters(sc, SO@meta.data[rownames(sc$metaData), paste0("RNA_snn_res.", SoupX.resolution)])
    sc = SoupX::autoEstCont(sc)
    out = SoupX::adjustCounts(sc)

    SO_sx <-
      Seurat::CreateSeuratObject(counts = out) %>%
      Seurat::NormalizeData() %>%
      Seurat::FindVariableFeatures(nfeatures = nhvf) %>%
      Seurat::ScaleData() %>%
      Seurat::RunPCA(npcs = npcs) %>%
      Seurat::RunUMAP(dims = 1:npcs) %>%
      Seurat::FindNeighbors(dims = 1:npcs) %>%
      Seurat::FindClusters(algorithm = 1, resolution = resolutions)

    sc <- SoupX::setDR(sc, SO_sx@reductions$umap@cell.embeddings)

    sc_info_df <- data.frame(n_expr_uncorrected = Matrix::rowSums(sc$toc > 0),
                             n_expr_corrected = Matrix::rowSums(out > 0)) %>%
      dplyr::mutate(abs_diff = n_expr_uncorrected-n_expr_corrected) %>%
      dplyr::mutate(rel_diff = abs_diff/n_expr_uncorrected) %>%
      dplyr::filter(abs_diff > 0) %>%
      tibble::rownames_to_column("Feature")

    message("Create a SoupX RNA assay as follows: SO[['SoupXRNA']] <- Seurat::CreateAssayObject(counts = soupx_matrix).")
    SOs <- list(SO, SO_sx)
    names(SOs) <- c("original", "SoupX")
  } else {
    SOs <- list(SO)
  }

  results <- lapply(SOs, function(SO) {

    counts <- as.matrix(Seurat::GetAssayData(SO, slot = "counts"))
    if (any(matrixStats::colSums2(counts) < min_nCount_RNA)) {
      warning(paste0("Transcriptomes (cells) with less than ", min_nCount_RNA, " total transcripts found. These will be excluded from finding doublets. They will have NA as dbl_score."))
    }


    dbl_score <- setNames(scDblFinder::computeDoubletDensity(x = counts[,which(matrixStats::colSums2(counts) >= min_nCount_RNA)],
                                                             subset.row = Seurat::VariableFeatures(SO),
                                                             dims = npcs, ...),
                          nm = colnames(counts[,which(matrixStats::colSums2(counts) >= min_nCount_RNA)]))

    SO <- Seurat::AddMetaData(SO, dbl_score, "dbl_score")
    SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^MT-"), "pct_mt")
    for (i in c("nCount_RNA", "nFeature_RNA", "pct_mt", "dbl_score")) {
      SO@meta.data[,paste0(i, "_log")] <- log1p(SO@meta.data[,i])
    }
    lm <- stats::lm(nFeature_RNA_log~nCount_RNA_log, data = SO@meta.data)
    SO@meta.data$residuals <- stats::residuals(lm)
    SO@meta.data$residuals_norm <- SO@meta.data$residuals/SO@meta.data$nCount_RNA_log/SO@meta.data$nFeature_RNA_log

    qc_p1 <- scexpr::feature_plot(SO, features = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log", "residuals_norm", paste0("RNA_snn_res.", resolutions[length(resolutions)])),
                                  reduction = "UMAP", legend.position = c(0,1), plot.labels = "text", label.size = 6)

    qc_p2 <- cowplot::plot_grid(scexpr::qc_params_meta_cols(SO,
                                                            qc.cols = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log", "residuals_norm"),
                                                            meta.cols = paste0("RNA_snn_res.", resolutions[length(resolutions)]),
                                                            strip.text = element_text(size = 6), panel.grid.minor = ggplot2::element_blank()),
                                ggplot2::ggplot(SO@meta.data, ggplot2::aes(nCount_RNA_log, nFeature_RNA_log, color = !!rlang::sym(paste0("RNA_snn_res.", resolutions[length(resolutions)])))) +
                                  ggplot2::geom_point() +
                                  ggplot2::theme_bw() +
                                  ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                                                 legend.key = element_blank(), legend.key.size = ggplot2::unit(0.3, "cm")) +
                                  ggplot2::scale_color_manual(values = scexpr::col_pal()) +
                                  ggplot2::labs(color = ""),
                                ggplot2::ggplot(dplyr::arrange(SO@meta.data, dbl_score_log), ggplot2::aes(nCount_RNA_log, nFeature_RNA_log, color = dbl_score_log)) +
                                  ggplot2::geom_point() +
                                  ggplot2::theme_bw() +
                                  ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
                                  ggplot2::scale_color_gradientn(colors = scexpr::col_pal("spectral")),
                                ggplot2::ggplot(dplyr::arrange(SO@meta.data, pct_mt_log), ggplot2::aes(nCount_RNA_log, nFeature_RNA_log, color = pct_mt_log)) +
                                  ggplot2::geom_point() +
                                  ggplot2::theme_bw() +
                                  ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
                                  ggplot2::scale_color_gradientn(colors = scexpr::col_pal("spectral")),
                                align = "hv", axis = "tblr")

    qc_p3 <- cowplot::plot_grid(scexpr::feature_plot(SO,
                                                     features = paste0("RNA_snn_res.", resolutions),
                                                     reduction = "UMAP", legend.position = "none", plot.labels = "text", label.size = 3, nrow.combine = 1, plot.title = F),
                                scexpr::qc_params_meta_cols(SO,
                                                            qc.cols = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log"),
                                                            meta.cols = paste0("RNA_snn_res.", resolutions),
                                                            panel.grid.minor = ggplot2::element_blank()),
                                nrow = 2,
                                rel_heights = c(0.3,0.7),
                                align = "hv", axis = "tblr")

    #remove count slot to save memory
    return(list(SO = Seurat::DietSeurat(SO, assays = names(SO@assays), counts = F, dimreducs = names(SO@reductions))
                qc_p1 = qc_p1, qc_p2 = qc_p2, qc_p3 = qc_p3))
  })

  if (length(results) == 1) {
    return(results[[1]])
  } else {
    message("Use e.g. SoupX::plotChangeMap(x[['sc']], cleanedMatrix = SoupX::adjustCounts(x[['sc']]), geneSet = 'GNLY') + theme_bw() + theme(panel.grid = element_blank()) + scale_color_gradientn(colors = scexpr::col_pal('spectral'), na.value = 'grey95') to plot changes in counts.")
    results <- c(results, list(sc), list(sc_info_df))
    names(results) <- c(names(SOs), "sc", "sc_info")
    return(results)
  }
}


