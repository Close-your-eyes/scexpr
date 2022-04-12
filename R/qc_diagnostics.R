#' Run a bunch of default quality checks on scRNAseq data from 10X genomics (Cell Ranger)
#'
#' Taking Cell Rangers output, namely the (i) feature, (ii) barcode, (iii) matrix files
#' within respective folder(s), filtered_feature_bc_matrix and optionally additionally raw_feature_bc_matrix,
#' a very default Seurat Object is generated and quality metrics are computed. Moreover, these
#' metrics are used to cluster the cells. Groups of cells with low quality scores (e.g. high
#' mt-fraction plus low number of detected features) will likely clusters together. This allows to filter them easily.
#' If raw_feature_bc_matrix is provided, \href{https://github.com/constantAmateur/SoupX}{SoupX} may be run.
#' In that case the whole pipeline in run twice.
#' Namely once with the original count matrix and once with a corrected count matrix after running SoupX with default
#' settings. If raw_feature_bc_matrix is not available, only \href{https://github.com/campbio/celda}{DecontX}
#' is available to check for ambient RNA contamination. Both of these metrics (SoupX and DecontX estimation of ambient RNA)
#' become part of clustering by qc metrics if set to TRUE. \href{https://github.com/plger/scDblFinder}{scDblFinder} is run
#' to detect doublets which is another quality metric.
#' In addition to qc metrics, a number of principle components (PCs) from feature expression may be added to clustering and
#' dimension reduction as also the feature composition of low quality transcriptomes may be skewed.
#' This will likely cause these cells to cluster separately even more. Apart from clustering with meta data,
#' also a pure analysis with feature expression values only is run. That may allow subset-wise application
#' of additional filters for qc metrics after very definite low quality transcriptomes have been eliminated in
#' a first round.
#'
#' On error in scDblFinder: Increase nhvf or try to change any other parameter.
#'
#' @param data.dir list or vector of full paths to filtered_feature_bc_matrix and raw_feature_bc_matrix; or filtered_feature_bc_matrix only (prohibiting the use of SoupX);
#' if not named, then the name of the common parent folder is used
#' @param nhvf number of highly variable features for every of the procedures; may be subject to lower numbers (e.g. 500 or 2000)
#' at the risk of scDblFinder returning an error 'size factors should be positive'
#' @param npcs number or principle components to calculate, e.g. 12 for diverse data sets and 8 for isolated subsets
#' @param min_nCount_RNA minimum number of transcripts per cells to be considered for scDblFinder::computeDoubletDensity
#' @param resolution resolution (louvain algorithm) for clustering based on feature expression (UMI count matrix)
#' @param SoupX logical whether to run SoupX. If TRUE, raw_feature_bc_matrix is needed.
#' @param SoupX.resolution resolution for (louvain algorithm) SoupX analysis
#' @param cells vector of cell names to include, consider the trailing '-1' in cell names
#' @param invert_cells invert cell selection, if TRUE cell names provides in 'cells' are excluded
#' @param ... arguments to SoupX::autoEstCont prefixed with 'SoupX__'
#' @param decontX logical whether to run celda::decontX to estimate RNA soup (contaminating ambient RNA molecules)
#' @param qc_meta_resolution resolution (louvain algorithm) for clustering based on qc meta data and optionally additional PC dimensions (n_PCs_to_meta_clustering)
#' @param n_PCs_to_meta_clustering how many principle compoments (PCs) from phenotypic clustering to add to qc meta data;
#' this will generate a mixed clustering (PCs from phenotypes (RNA) and qc meta data like pct mt and nCount_RNA); the more PCs are added the greater the
#' phenotypic influence becomes
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
qc_diagnostic <- function(data.dir,
                          nhvf = 5000,
                          npcs = 10,
                          min_nCount_RNA = 100,
                          resolution = 0.8,
                          scDblFinder = T,
                          SoupX = F,
                          decontX = F,
                          SoupX.resolution = 0.6,
                          cells = NULL,
                          invert_cells = F,
                          qc_meta_resolution = 0.8,
                          n_PCs_to_meta_clustering = 2,
                          ...) {

  if (!requireNamespace("matrixStats", quietly = T)) {
    utils::install.packages("matrixStats")
  }
  if (!requireNamespace("uwot", quietly = T)) {
    utils::install.packages("uwot")
  }
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

  resolution <- as.numeric(gsub("^1.0$", "1", resolution))

  dots <- list(...)
  ## check dots for first letters

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

  message("Reading filtered_feature_bc_matrix data.")
  filt_data <- Seurat::Read10X(data.dir = grep("filtered_feature_bc_matrix", data.dir, value = T))
  if (is.list(filt_data)) {
    message("filtered_feature_bc_matrix is a list. Using 'Gene Expression' index")
    filt_data <- filt_data[["Gene Expression"]]
  }
  # filter cells with very low RNA_count
  nCount_RNA <- Matrix::colSums(filt_data) ## fix
  if (any(nCount_RNA < min_nCount_RNA)) {
    message(length(which(nCount_RNA < min_nCount_RNA)), " cells removed due to min_nCount_RNA.")
    filt_data <- filt_data[,which(nCount_RNA >= min_nCount_RNA)]
    if (ncol(filt_data) == 0) {
      stop("No cells left after filtering for min_nCount_RNA.")
    }
  }

  if (!is.null(cells)) {
    if (any(cells %in% colnames(filt_data))) {
      message(length(cells), " cells provided.")
      message(length(which(cells %in% colnames(filt_data))), " of cells from a total of ", ncol(filt_data), " cells found in data (", length(which(cells %in% colnames(filt_data)))/ncol(filt_data), " %).")
      cells <- cells[which(cells %in% colnames(filt_data))]
      if (invert_cells) {
        filt_data <- filt_data[,which(!colnames(filt_data) %in% cells)]
      } else {
        filt_data <- filt_data[,cells]
      }
    } else {
      message("Non of cells found in data.")
    }
  }

  message("Creating initial Seurat object with ", ncol(filt_data), " cells.")
  SO <-
    Seurat::CreateSeuratObject(counts = as.matrix(filt_data)) %>%
    Seurat::NormalizeData(verbose = F) %>%
    Seurat::FindVariableFeatures(nfeatures = nhvf, verbose = F) %>%
    Seurat::ScaleData(verbose = F) %>%
    Seurat::RunPCA(npcs = npcs, verbose = F) %>%
    Seurat::RunUMAP(dims = 1:npcs, verbose = F) %>%
    Seurat::FindNeighbors(dims = 1:npcs, verbose = F) %>%
    Seurat::FindClusters(algorithm = 1, resolution = resolution, verbose = F)

  if (SoupX) {
    # use filt_data which may have been reduced 'cells' selection; raw_feature_bc_matrix will provide the picture of the soup
    message("Reading raw_feature_bc_matrix data.")
    raw_data <- Seurat::Read10X(data.dir = grep("raw_feature_bc_matrix", data.dir, value = T))
    if (is.list(raw_data)) {
      raw_data <- raw_data[["Gene Expression"]]
    }
    message("Running SoupX.")
    sc <- SoupX::SoupChannel(tod = raw_data, toc = filt_data)
    if (SoupX.resolution != resolution) {
      SO <- Seurat::FindClusters(SO, algorithm = 1, resolution = SoupX.resolution, verbose = F)
    }
    sc = SoupX::setClusters(sc, SO@meta.data[rownames(sc$metaData), paste0("RNA_snn_res.", SoupX.resolution)])

    temp_dots <- dots[which(grepl("^SoupX__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^SoupX__", "", names(temp_dots), ignore.case = T)
    temp_dots <- temp_dots[which(names(temp_dots) %in% formals(SoupX::autoEstCont))]
    sc <- do.call(SoupX::autoEstCont, args = c(list(sc = sc), temp_dots))

    out = SoupX::adjustCounts(sc)
    # add percentage of soup as meta data, similar to decontX
    SO <- Seurat::AddMetaData(SO, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(out))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

    message("Creating Seurat object on SoupX-corrected count matrix with ", ncol(filt_data), " cells.")
    SO_sx <-
      Seurat::CreateSeuratObject(counts = out, verbose = F) %>%
      Seurat::NormalizeData(verbose = F) %>%
      Seurat::FindVariableFeatures(nfeatures = nhvf, verbose = F) %>%
      Seurat::ScaleData(verbose = F) %>%
      Seurat::RunPCA(npcs = npcs, verbose = F) %>%
      Seurat::RunUMAP(dims = 1:npcs, verbose = F) %>%
      Seurat::FindNeighbors(dims = 1:npcs, verbose = F) %>%
      Seurat::FindClusters(algorithm = 1, resolution = resolution, verbose = F)

    SO_sx <- Seurat::AddMetaData(SO_sx, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(out))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

    sc <- SoupX::setDR(sc, SO_sx@reductions$umap@cell.embeddings)

    sc_info_df <- data.frame(n_expr_uncorrected = Matrix::rowSums(sc$toc > 0),
                             n_expr_corrected = Matrix::rowSums(out > 0)) %>%
      dplyr::mutate(abs_diff = n_expr_uncorrected-n_expr_corrected) %>%
      dplyr::mutate(rel_diff = abs_diff/n_expr_uncorrected) %>%
      dplyr::filter(abs_diff > 0) %>%
      tibble::rownames_to_column("Feature")

    message("Optionally: Create a SoupX RNA assay as follows: SO[['SoupXRNA']] <- Seurat::CreateAssayObject(counts = soupx_matrix).")
    SO <- list(SO, SO_sx)
    names(SO) <- c("original", "SoupX")
  } else {
    SO <- list(SO)
  }

  results <- lapply(SO, function(SOx) {

    message("Running scDblFinder.")
    dbl_score <- NULL
    if (scDblFinder) {
      dbl_score <- tryCatch({
        log1p(scDblFinder::computeDoubletDensity(x = as.matrix(Seurat::GetAssayData(SOx, slot = "counts")),
                                                 subset.row = Seurat::VariableFeatures(SOx),
                                                 dims = npcs))
      }, error = function(error_condition) {
        message("doublet calculation failed. Ty to increase nhvf.")
        message(error_condition)
        return(NULL)
      })
    }

    if (is.null(dbl_score)) {
      qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt")
    } else {
      dbl_score <- stats::setNames(dbl_score, nm = colnames(Seurat::GetAssayData(SOx, slot = "counts")))
      SOx <- Seurat::AddMetaData(SOx, dbl_score, "dbl_score_log")
      qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt", "dbl_score")
    }


    # slow !!
    'temp_dots <- dots[which(grepl("^scDbl__", names(dots), ignore.case = T))]
    names(temp_dots) <- gsub("^scDbl__", "", names(temp_dots), ignore.case = T)
    dbl_score <- do.call(scDblFinder::computeDoubletDensity, args = c(list(x = counts[,which(matrixStats::colSums2(counts) >= 2000)],
                                                                           subset.row = Seurat::VariableFeatures(SOx),
                                                                           dims = npcs),
                                                                      temp_dots))'


    # impute missing dbl score values with medians
    #SOx@meta.data$dbl_score_log[which(is.na(SOx@meta.data$dbl_score_log))] <- median(SOx@meta.data$dbl_score_log, na.rm = T)

    # differentiate mouse, human or no MT-genes at all
    if (any(grepl("^MT-", rownames(SOx))) && !any(grepl("^mt-", rownames(SOx)))) {
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^MT-"), "pct_mt")
    } else if (!any(grepl("^MT-", rownames(SOx))) && any(grepl("^mt-", rownames(SOx)))) {
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^mt-"), "pct_mt")
    } else {
      message("No mitochondrial genes could be identified from gene names - none starting with MT- (human) or mt- (mouse).")
      qc_cols <- qc_cols[-which(qc_cols == "pct_mt")]
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
    message("Running dimension reduction and clustering on qc meta data.")
    meta <- dplyr::select(SOx@meta.data, dplyr::all_of(qc_cols), residuals)
    if (n_PCs_to_meta_clustering > 0) {
      meta <- cbind(meta, SOx@reductions[["pca"]]@cell.embeddings[,1:n_PCs_to_meta_clustering])
    }
    umap_dims <- uwot::umap(scale_min_max(meta), metric = "cosine")
    colnames(umap_dims) <- c("meta_UMAP_1", "meta_UMAP_2")
    SOx[["umapmeta"]] <- Seurat::CreateDimReducObject(embeddings = umap_dims, key = "UMAPMETA_", assay = "RNA")

    clusters <- Seurat::FindClusters(Seurat::FindNeighbors(scale_min_max(meta), annoy.metric = "cosine", verbose = F)$snn, resolution = qc_meta_resolution, verbose = F)
    colnames(clusters) <- paste0("meta_", colnames(clusters))
    SOx <- Seurat::AddMetaData(SOx, cbind(umap_dims, clusters))

    meta <- cbind(meta, cbind(umap_dims, clusters))
    # dependencies below!
    meta_cols <- c(paste0("RNA_snn_res.", resolution), colnames(clusters))

    qc_p1 <- suppressMessages(feature_plot(SOx,
                                           features = c(qc_cols, meta_cols[1]),
                                           reduction = "UMAP", legend.position = "none",
                                           plot.labels = "text", label.size = 6))

    qc_p2 <- ggplot2::ggplot(tidyr::pivot_longer(SOx@meta.data[,c(qc_cols, meta_cols)], cols = dplyr::all_of(qc_cols), names_to = "qc_param", values_to = "value"),
                             ggplot2::aes(x = !!rlang::sym(meta_cols[1]), y = value, color = !!rlang::sym(meta_cols[2]))) +
      ggplot2::geom_jitter(width = 0.2, size = 0.3) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                     legend.key.size = ggplot2::unit(0.3, "cm"), legend.key = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = col_pal()) +
      ggplot2::scale_y_continuous(sec.axis = sec_axis(~ expm1(.))) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
      ggplot2::facet_wrap(ggplot2::vars(qc_param), scales = "free_y", ncol = 1)

    p3_1 <- patchwork::wrap_plots(feature_plot(SOx, features = meta_cols[2], reduction = "umapmeta", pt.size = 0.5, legend.position = "none",
                                               label.size = 6, plot.labels = "text", plot.title = F),
                                  suppressMessages(freq_pie_chart(SO = SOx, meta.col = meta_cols[2])),
                                  ncol = 1)

    rrr <- c(seq(0, 1e1, 2e0),
             seq(0, 1e2, 2e1),
             seq(0, 1e3, 2e2),
             seq(0, 1e4, 2e3),
             seq(0, 1e5, 2e4))
    rrr <- rrr[which(rrr != 0)]
    p3_2 <- feature_plot_stat(SOx,
                              features = "nCount_RNA_log",
                              meta.col = meta_cols[2],
                              geom2 = "none",
                              jitterwidth = 0.9,
                              panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                              axis.title.y = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_blank(),
                              axis.title.x = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(sec.axis = sec_axis(~ expm1(.), breaks = rrr[intersect(which(rrr > min(expm1(SOx@meta.data$nCount_RNA_log))),
                                                                                         which(rrr < max(expm1(SOx@meta.data$nCount_RNA_log))))]))

    p3_3 <- feature_plot_stat(SOx,
                              features = "nFeature_RNA_log",
                              meta.col = meta_cols[2],
                              geom2 = "none",
                              jitterwidth = 0.9,
                              panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                              axis.title.y = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_blank(),
                              axis.title.x = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(sec.axis = sec_axis(~ expm1(.), breaks = rrr[intersect(which(rrr > min(expm1(SOx@meta.data$nFeature_RNA_log))),
                                                                                         which(rrr < max(expm1(SOx@meta.data$nFeature_RNA_log))))]))

    p3_4 <- feature_plot_stat(SOx,
                              features = "pct_mt_log",
                              meta.col = meta_cols[2],
                              geom2 = "none",
                              jitterwidth = 0.9,
                              panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                              axis.title.y = ggplot2::element_blank()) +
      ggplot2::scale_y_continuous(sec.axis = sec_axis(~ expm1(.), breaks = rrr[intersect(which(rrr > min(expm1(SOx@meta.data$pct_mt_log))),
                                                                                         which(rrr < max(expm1(SOx@meta.data$pct_mt_log))))]))

    p3_x <- lapply(c(qc_cols[which(!grepl("nCount_RNA_log|nFeature_RNA_log|pct_mt_log", qc_cols))], "residuals"), function(qcf) {
      p3_x <- feature_plot_stat(SOx,
                                features = qcf,
                                meta.col = meta_cols[2],
                                geom2 = "none",
                                jitterwidth = 0.9,
                                panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                                axis.title.y = ggplot2::element_blank())
      if (qcf != "residuals") {
        p3_x <-
          p3_x +
          ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())
      }
      return(p3_x)
    })


    p3_2 <- patchwork::wrap_plots(p3_2, p3_3, p3_4, ncol = 1)
    p3_3 <- patchwork::wrap_plots(p3_x, ncol = 1)
    qc_p3 <- patchwork::wrap_plots(p3_1, p3_2, p3_3, nrow = 1)


    message("Calculating phenotype cluster markers.")
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
                phenotype_clusters_plot = qc_p1, meta_vs_phenotype_clusters_plot = qc_p2, meta_clusters_plot = qc_p3,
                phenotype_cluster_markers = cluster_marker_list))
  })

  if (length(results) == 1) {
    return(results)
  } else {
    message("Use e.g. SoupX::plotChangeMap(x[['sc']], cleanedMatrix = SoupX::adjustCounts(x[['sc']]), geneSet = 'GNLY') + theme_bw() + theme(panel.grid = element_blank()) + scale_color_gradientn(colors = col_pal('spectral'), na.value = 'grey95') to plot changes in counts.")
    results <- c(results, list(sc), list(sc_info_df))
    names(results) <- c(names(SO), "sc", "sc_info")
    return(results)
  }
}

ceiling_any = function(x, accuracy, f = ceiling) {
  f(x/ accuracy) * accuracy
}

floor_any = function(x, accuracy, f = floor) {
  f(x/ accuracy) * accuracy
}
