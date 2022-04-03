#' Title
#'
#' On error in scDblFinder: Increase nhvf or try to change any other parameter
#' SoupX and/or decontX to filter cells or for correction
#'
#' @param data.dir list or vector of full paths to filtered_feature_bc_matrix and raw_feature_bc_matrix, or filtered_feature_bc_matrix only (prohibiting the use of SoupX);
#' if not named, then the name of the common parent folder is uses as name(s)
#' @param nhvf
#' @param npcs
#' @param min_nCount_RNA
#' @param resolution
#' @param SoupX logical whether to run SoupX. If T, raw_feature_bc_matrix is needed.
#' @param SoupX.resolution
#' @param cells select to include
#' @param invert_cells invert cell selection, if TRUE cells are excluded
#' @param ...
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
qc_diagnostic <- function(data.dir,
                          nhvf = 1000,
                          npcs = 10,
                          min_nCount_RNA = 300,
                          resolution = 0.9,
                          SoupX = F,
                          decontX = F,
                          SoupX.resolution = 0.6,
                          cells = NULL,
                          invert_cells = F,
                          qc_meta_resolution = 0.6,
                          ...) {

  if (!requireNamespace("matrixStats", quietly = T)) {
    utils::install.packages("matrixStats")
  }
  if (!requireNamespace("uwot", quietly = T)) {
    utils::install.packages("uwot")
  }
  #if (!requireNamespace("ggpointdensity", quietly = T)) {utils::install.packages("ggpointdensity")}
  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github("immunogenomics/presto")
  }
  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("scDblFinder", quietly = T)) {
    BiocManager::install("scDblFinder")
  }
  if (decontX && !requireNamespace("celda", quietly = T)) {
    BiocManager::install("celda")
  }
  if (SoupX && !requireNamespace("SoupX", quietly = T)) {
    install.packages("SoupX")
  }

  resolution <- as.numeric(gsub("1.0", "1", resolution))

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
      names(data.dir) <- basename(unique(dirname(data.dir)))
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

  if (!is.numeric(qc_meta_resolution) | length(qc_meta_resolution) != 1) {
    stop("qc_meta_resolution has to a numeric of length one. It is passed to the resolution parameter of Seurat::FindClusters.")
  }


  filt_data <- Seurat::Read10X(data.dir = grep("filtered_feature_bc_matrix", data.dir, value = T))
  if (is.list(filt_data)) {
    filt_data <- filt_data[["Gene Expression"]]
  }

  if (!is.null(cells)) {
    if (any(cells) %in% colnames(filt_data)) {
      message(length(which((any(cells) %in% colnames(filt_data)))), " names from cells found in column names of data.")
      cells_sub <- cells[which(cells %in% colnames(filt_data))]
      if (invert_cells) {
        filt_data <- filt_data[,which(!colnames(filt_data) %in% cells_sub)]
      } else {
        filt_data <- filt_data[,cells_sub]
      }
    } else {
      message("Non of cells found in column names of data.")
    }

    cells_sub <- cells[which(cells %in% colnames(filt_data))]
    if (invert_cells) {
      filt_data <- filt_data[,which(!colnames(filt_data) %in% cells_sub)]
    } else {
      filt_data <- filt_data[,cells_sub]
    }
  }

  SO <-
    Seurat::CreateSeuratObject(counts = as.matrix(filt_data)) %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures(nfeatures = nhvf) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = npcs) %>%
    Seurat::RunUMAP(dims = 1:npcs) %>%
    Seurat::FindNeighbors(dims = 1:npcs) %>%
    Seurat::FindClusters(algorithm = 1, resolution = resolution)

  if (SoupX) {
    # use filt_data which may have been reduced 'cells' selection; raw_feature_bc_matrix will provide the picture of the soup
    raw_data <- Seurat::Read10X(data.dir = grep("raw_feature_bc_matrix", data.dir, value = T))
    if (is.list(raw_data)) {
      raw_data <- raw_data[["Gene Expression"]]
    }
    sc <- SoupX::SoupChannel(tod = raw_data, toc = filt_data)
    if (SoupX.resolution != resolution) {
      SO <- Seurat::FindClusters(SO, algorithm = 1, resolution = SoupX.resolution)
    }
    sc = SoupX::setClusters(sc, SO@meta.data[rownames(sc$metaData), paste0("RNA_snn_res.", SoupX.resolution)])
    sc = SoupX::autoEstCont(sc)
    out = SoupX::adjustCounts(sc)
    # add percentage of soup as meta data, similar to decontX
    SO <- Seurat::AddMetaData(SO, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(out))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

    SO_sx <-
      Seurat::CreateSeuratObject(counts = out) %>%
      Seurat::NormalizeData() %>%
      Seurat::FindVariableFeatures(nfeatures = nhvf) %>%
      Seurat::ScaleData() %>%
      Seurat::RunPCA(npcs = npcs) %>%
      Seurat::RunUMAP(dims = 1:npcs) %>%
      Seurat::FindNeighbors(dims = 1:npcs) %>%
      Seurat::FindClusters(algorithm = 1, resolution = resolution)
    SO_sx <- Seurat::AddMetaData(SO_sx, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(out))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

    sc <- SoupX::setDR(sc, SO_sx@reductions$umap@cell.embeddings)

    sc_info_df <- data.frame(n_expr_uncorrected = Matrix::rowSums(sc$toc > 0),
                             n_expr_corrected = Matrix::rowSums(out > 0)) %>%
      dplyr::mutate(abs_diff = n_expr_uncorrected-n_expr_corrected) %>%
      dplyr::mutate(rel_diff = abs_diff/n_expr_uncorrected) %>%
      dplyr::filter(abs_diff > 0) %>%
      tibble::rownames_to_column("Feature")

    message("Create a SoupX RNA assay as follows: SO[['SoupXRNA']] <- Seurat::CreateAssayObject(counts = soupx_matrix).")
    SO <- list(SO, SO_sx)
    names(SO) <- c("original", "SoupX")
  } else {
    SO <- list(SO)
  }

  results <- lapply(SO, function(SOx) {

    counts <- as.matrix(Seurat::GetAssayData(SOx, slot = "counts"))
    if (any(matrixStats::colSums2(counts) < min_nCount_RNA)) {
      warning(paste0("Transcriptomes (cells) with less than ", min_nCount_RNA, " total transcripts found. These will be excluded from finding doublets. They will have NA as dbl_score."))
    }


    dbl_score <- log1p(stats::setNames(scDblFinder::computeDoubletDensity(x = counts[,which(matrixStats::colSums2(counts) >= min_nCount_RNA)],
                                                                          subset.row = Seurat::VariableFeatures(SOx),
                                                                          dims = npcs, ...),
                                       nm = colnames(counts[,which(matrixStats::colSums2(counts) >= min_nCount_RNA)])))

    SOx <- Seurat::AddMetaData(SOx, dbl_score, "dbl_score_log")
    # differentiate mouse, human or no MT-genes at all
    qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt", "dbl_score")
    if (any(grepl("^MT-", rownames(SOx))) && !any(grepl("^mt-", rownames(SOx)))) {
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^MT-"), "pct_mt")
    } else if (!any(grepl("^MT-", rownames(SOx))) && any(grepl("^mt-", rownames(SOx)))) {
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^mt-"), "pct_mt")
    } else {
      message("No mitochondrial genes could be identified from gene names - none starting with MT- (human) or mt- (mouse).")
      qc_cols <- c("nCount_RNA", "nFeature_RNA", "dbl_score")
    }
    qc_cols <- paste0(qc_cols, "_log")

    if (decontX) {
      SOx <- Seurat::AddMetaData(SOx, celda::decontX(x = Seurat::GetAssayData(SOx),
                                                     z = SOx@meta.data[[paste0("RNA_snn_res.", resolution)]],
                                                     varGenes = nhvf)[["contamination"]], "pct_soup_decontX")
      qc_cols <- c(qc_cols, "pct_soup_decontX")
    }
    if (SoupX) {
      qc_cols <- c(qc_cols, "pct_soup_SoupX")
    }

    SOx@meta.data$nFeature_RNA_log <- log1p(SOx@meta.data$nFeature_RNA)
    SOx@meta.data$nCount_RNA_log <- log1p(SOx@meta.data$nCount_RNA)
    SOx@meta.data$pct_mt_log <- log1p(SOx@meta.data$pct_mt)
    SOx@meta.data$residuals <- stats::residuals(stats::lm(nFeature_RNA_log~nCount_RNA_log, data = SOx@meta.data))
    ## clustering on meta data (quality metrics)
    meta <- dplyr::select(SOx@meta.data, dplyr::all_of(qc_cols), residuals)
    #meta <- scale_min_max(meta) # leave unscaled, only scale for umap and finding neighbors
    # scale(cbind(meta, SO@reductions[["pca"]]@cell.embeddings))
    umap_dims <- uwot::umap(scale_min_max(meta), metric = "cosine")
    colnames(umap_dims) <- c("meta_UMAP_1", "meta_UMAP_2")
    SOx[["umapmeta"]] <- Seurat::CreateDimReducObject(embeddings = umap_dims, key = "UMAPMETA_", assay = "RNA")

    clusters <- Seurat::FindClusters(Seurat::FindNeighbors(scale_min_max(meta), annoy.metric = "cosine")$snn, resolution = qc_meta_resolution)
    colnames(clusters) <- paste0("meta_", colnames(clusters))
    SOx <- Seurat::AddMetaData(SOx, cbind(umap_dims, clusters))

    meta <- cbind(meta, cbind(umap_dims, clusters))
    # dependencies below!
    meta_cols <- c(paste0("RNA_snn_res.", resolution), colnames(clusters))

    qc_p1 <- suppressMessages(feature_plot(SOx,
                                           features = c(qc_cols, meta_cols[1]),
                                           reduction = "UMAP", legend.position = "none",
                                           plot.labels = "text", label.size = 6))

    qc_p2 <- cowplot::plot_grid(ggplot2::ggplot(tidyr::pivot_longer(SOx@meta.data[,c(qc_cols, meta_cols)], cols = dplyr::all_of(qc_cols), names_to = "qc_param", values_to = "value"),
                                                ggplot2::aes(x = !!rlang::sym(meta_cols[1]), y = value, color = !!rlang::sym(meta_cols[2]))) +
                                  ggplot2::geom_jitter(width = 0.2, size = 0.3) +
                                  ggplot2::theme_bw() +
                                  ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                                                 panel.grid.minor.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                                                 legend.key.size = ggplot2::unit(0.3, "cm"), legend.key = ggplot2::element_blank()) +
                                  ggplot2::scale_color_manual(values = col_pal()) +
                                  ggplot2::scale_y_continuous(sec.axis = sec_axis(~ expm1(.))) +
                                  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
                                  ggplot2::facet_wrap(ggplot2::vars(qc_param), scales = "free_y", ncol = 1),
                                align = "hv", axis = "tblr")



    qc_p3 <- cowplot::plot_grid(
      cowplot::plot_grid(feature_plot(SOx, features = meta_cols[2], reduction = "umapmeta", pt.size = 0.5, legend.position = "none",
                                      label.size = 6, plot.labels = "text"),
                         suppressMessages(scexpr:::freq_pie_chart(SO = SOx, meta.col = meta_cols[2])),
                         ncol = 1, align = "hv", axis = "tblr"),
      feature_plot_stat(SOx,
                        features = c(qc_cols, "residuals"),
                        meta.col = meta_cols[2],
                        geom2 = "none",
                        jitterwidth = 0.9,
                        panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                        axis.title.y = ggplot2::element_blank()),
      align = "hv", axis = "tblr", rel_widths = c(1/3,2/3))


    cluster_marker_list <- lapply(paste0("RNA_snn_res.", resolution), function(x) {
      presto::wilcoxauc(SOx, group_by = x, seurat_assay = "RNA", assay = "data")  %>%
        dplyr::filter(padj < 0.0001) %>%
        dplyr::filter(abs(logFC) > 0.3) %>%
        dplyr::select(-c(statistic, pval)) %>%
        dplyr::mutate(pct_in = round(pct_in,2), pct_out = round(pct_out,2))
    })
    names(cluster_marker_list) <- paste0("RNA_snn_res.", resolution)

    #remove count slot to save memory
    return(list(SO = Seurat::DietSeurat(SOx, assays = names(SOx@assays), counts = F, dimreducs = names(SOx@reductions)),
                qc_p1 = qc_p1, qc_p2 = qc_p2, qc_p3 = qc_p3,
                cluster_markers = cluster_marker_list))
  })

  if (length(results) == 1) {
    return(results[[1]])
  } else {
    message("Use e.g. SoupX::plotChangeMap(x[['sc']], cleanedMatrix = SoupX::adjustCounts(x[['sc']]), geneSet = 'GNLY') + theme_bw() + theme(panel.grid = element_blank()) + scale_color_gradientn(colors = col_pal('spectral'), na.value = 'grey95') to plot changes in counts.")
    results <- c(results, list(sc), list(sc_info_df))
    names(results) <- c(names(SO), "sc", "sc_info")
    return(results)
  }
}



'cowplot::plot_grid(
                                  ggplot2::ggplot(SOx@meta.data, ggplot2::aes(nCount_RNA, nFeature_RNA, color = !!rlang::sym(meta_cols[1]))) +
                                    ggplot2::geom_point(size = 0.4) +
                                    ggplot2::theme_bw() +
                                    ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0),
                                                   legend.background = ggplot2::element_blank(),
                                                   panel.grid = ggplot2::element_blank(), legend.key = ggplot2::element_blank(),
                                                   legend.key.size = ggplot2::unit(0.4, "cm")) +
                                    ggplot2::scale_color_manual(values = col_pal()) +
                                    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3), ncol = 3)) +
                                    ggplot2::labs(color = ""),
                                  ggplot2::ggplot(dplyr::arrange(SOx@meta.data, dbl_score), ggplot2::aes(nCount_RNA, nFeature_RNA, color = dbl_score)) +
                                    ggplot2::geom_point(size = 0.4) +
                                    ggplot2::theme_bw() +
                                    ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
                                    ggplot2::scale_color_gradientn(colors = col_pal("spectral"))  +
                                    ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 0.5, barheight = 3, label.theme = ggplot2::element_text(size = 6), title.theme = ggplot2::element_text(size = 8))),
                                  ggplot2::ggplot(dplyr::arrange(SOx@meta.data, pct_mt), ggplot2::aes(nCount_RNA, nFeature_RNA, color = pct_mt)) +
                                    ggplot2::geom_point(size = 0.4) +
                                    ggplot2::theme_bw() +
                                    ggplot2::theme(legend.justification = c(1,0), legend.position = c(1,0), legend.background = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
                                    ggplot2::scale_color_gradientn(colors = col_pal("spectral"))  +
                                    ggplot2::guides(color = ggplot2::guide_colorbar(barwidth = 0.5, barheight = 3, label.theme = ggplot2::element_text(size = 6), title.theme = ggplot2::element_text(size = 8))),
                                  ncol = 1, align = "hv", axis = "tblr")'


'ggplot2::ggplot(meta, ggplot2::aes(x = meta_UMAP_1, y = meta_UMAP_2)) +
                           ggpointdensity::geom_pointdensity(size = 0.3) +
                           ggplot2::theme_bw() +
                           ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank()) +
                           ggplot2::scale_color_gradientn(colors = col_pal("spectral"))
'

'    feature_plot_stat(SOx, features = qc_cols, meta.col = meta_cols[2], col.pal = "custom", geom2 = "none") +
      scale_color_manual(values = col_pal())
'
