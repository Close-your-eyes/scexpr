#' Run a bunch of default quality checks on scRNAseq data from 10X genomics (Cell Ranger)
#'
#' Taking Cell Rangers output, namely the (i) feature, (ii) barcode, (iii) matrix files
#' within respective folder(s), filtered_feature_bc_matrix and optionally additionally raw_feature_bc_matrix,
#' a very default Seurat Object is generated and quality metrics are computed. Moreover, these
#' metrics are used to cluster the cells. Groups of cells with low quality scores (e.g. high
#' mt-fraction plus low number of detected features) will likely clusters together. This allows to filter them easily.
#' If raw_feature_bc_matrix is provided, \href{https://github.com/constantAmateur/SoupX}{SoupX} may be run.
#' In that case the whole pipeline in run twice.
#' Namely, once with the original count matrix and once with a corrected count matrix after running SoupX with default
#' settings. If raw_feature_bc_matrix is not available, only \href{https://github.com/campbio/celda}{DecontX}
#' is available to check for ambient RNA contamination. Both of these metrics (SoupX and DecontX estimation of ambient RNA)
#' become part of clustering by qc metrics if set to TRUE. \href{https://github.com/plger/scDblFinder}{scDblFinder} is run
#' to detect doublets which is another quality metric.
#' In addition to qc metrics, a number of principle components (PCs) from feature expression may be added to clustering and
#' dimension reduction as also the feature composition of low quality transcriptomes may be skewed.
#' This will likely cause these cells to cluster separately even more. Apart from clustering with meta data,
#' also a pure analysis with feature expression values only is run. That may allow cluster-wise application
#' of additional filters for qc metrics after very definite low quality transcriptomes have been eliminated in
#' a first round based on qc metric clustering.
#' If data_dirs contains multiple samples then integration of samples is done with \href{"https://github.com/immunogenomics/harmony"}{harmony}.
#' (i) Detection of soup (ambient RNA) by SoupX and decontX, (ii) detection of doublets and (iii) calculation of residuals from the linear model
#' of nCount_RNA_log vs nFeature_RNA_log is done sample-wise when multiple data_dirs are detected/provided. Results are written into
#' the common Seurat object though, the merged and harmonzized PCA space of which is subject for clustering the cells based on feature expression (phenotypes)
#'
#'
#' @param data_dirs list or vector of parent direction(s) which will be search for folders called "filtered_feature_bc_matrix";
#' on the same level where each of these folders is found, a raw_feature_bc_matrix folder may exist to enable SoupX; if one "raw_feature_bc_matrix"
#' is missing, SoupX is disable for all others
#' @param nhvf number of highly variable features for every of the procedures
#' @param npcs number or principle components to calculate, e.g. 12 for diverse data sets and 8 for isolated subsets
#' @param resolution resolution (louvain algorithm) for clustering based on feature expression
#' @param SoupX logical whether to run SoupX. If TRUE, raw_feature_bc_matrix is needed.
#' @param resolution_SoupX resolution (louvain algorithm) for SoupX analysis
#' @param cells vector of cell names to include, consider the trailing '-1' in cell names
#' @param invert_cells invert cell selection, if TRUE cell names provides in 'cells' are excluded
#' @param decontX logical whether to run celda::decontX to estimate RNA soup (contaminating ambient RNA molecules)
#' @param resolution_meta resolution(s) (louvain algorithm) for clustering based on qc meta data and optionally additional PC dimensions (n_PCs_to_meta_clustering)
#' @param n_PCs_to_meta_clustering how many principle components (PCs) from phenotypic clustering to add to qc meta data;
#' this will generate a mixed clustering (PCs from phenotypes (RNA) and qc meta data like pct mt and nCount_RNA); the more PCs are added the greater the
#' phenotypic influence becomes; one or more integers can be supplied to explore the effect; pass 0, to have no PCs included in meta clustering; e.g. when
#' n_PCs_to_meta_clustering = 3 PCs 1-3 are used, when n_PCs_to_meta_clustering = 1 only PC 1 is used.
#' @param ... additional arguments to SoupX::autoEstCont
#' @param scDblFinder logical, whether to run doublet detection algorithm from scDblFinder
#' @param return_SoupX logical whether to return a full Seurat-object and diagnostics from SoupX (TRUE) or whether to run SoupX without these returns and just
#' have the Soup-metric included as an additional quality-control metric along with pct_mt and nCount_RNA etc. Will be set to FALSE if more than one data_dir with
#' filtered_feature_bc_matrix is supplied. So, only possible when data set are provided one by one.
#' @param feature_rm character vector of features to remove from count matrices;
#' removal is done after aggregation (if feature_aggr is provided)
#' @param feature_aggr named list of character vectors of features to aggregate;
#' names of of list entries are names of the aggregated feature; aggregation of counts
#' is simply done by addition; aggregation is done before feature removal (if feature_rm is provided)
#' @param min_UMI min number of UMI per cell (colSums)
#' @param min_UMI_var_feat min number of UMI in variable featured per cell (colSums); only relevant if scDblFinder = T
#' @param ffbms named paths to folders with raw/filtered feature matrix files or .h5 files; this will skip any procedure with data_dirs
#' @param rfbms named paths to folders with filtered/raw feature matrix files or .h5 files; this will skip any procedure with data_dirs
#' @param batch_corr which batch correction to apply for integration of different samples
#'
#' @return a list of Seurat object and data frame with marker genes for clusters based on feature expression
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
qc_diagnostic <- function(data_dirs,
                          nhvf = 2000,
                          npcs = 10,
                          resolution = 0.8,
                          resolution_SoupX = 0.6,
                          resolution_meta = 0.8,
                          n_PCs_to_meta_clustering = 2,
                          scDblFinder = T,
                          min_UMI = 30,
                          min_UMI_var_feat = 5,
                          SoupX = F,
                          decontX = F,
                          return_SoupX = T,
                          cells = NULL,
                          invert_cells = F,
                          feature_rm = NULL,
                          feature_aggr = NULL,
                          ffbms = NULL,
                          rfbms = NULL,
                          batch_corr = c("harmony", "none"),
                          ...) {

  if (!requireNamespace("matrixStats", quietly = T)) {
    utils::install.packages("matrixStats")
  }
  if (!requireNamespace("hdf5r", quietly = T)) {
    utils::install.packages("hdf5r")
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
    utils::install.packages("SoupX")
  }
  if (!requireNamespace("patchwork", quietly = T)) {
    utils::install.packages("patchwork")
  }
  if (!requireNamespace("scuttle", quietly = T)) {
    BiocManager::install("scuttle")
  }


  resolution <- as.numeric(gsub("^1.0$", "1", resolution))

  if (!is.numeric(resolution_meta)) {
    stop("resolution_meta has to be numeric.")
  }
  if (!is.numeric(resolution_SoupX) || length(resolution_SoupX) != 1) {
    stop("resolution_SoupX has to be numeric and of length 1.")
  }
  if (!is.numeric(resolution) || length(resolution) != 1) {
    stop("resolution has to be numeric and of length 1.")
  }
  if (!is.numeric(n_PCs_to_meta_clustering)) {
    stop("n_PCs_to_meta_clustering should be numeric.")
  }

  batch_corr <- match.arg(batch_corr, c("harmony", "none"))

  #dots <- list(...)

  if (is.null(ffbms) && is.null(rfbms)) {
    checked_dirs <- check_dir(data_dirs = data_dirs, SoupX = SoupX)
    data_dirs <- checked_dirs[[1]]

    if (SoupX && !checked_dirs[[2]]) {
      message("raw_feature_bc_matrix not found in every data_dir. SoupX set to FALSE.")
    }
    SoupX <- checked_dirs[[2]]
    if (return_SoupX && !checked_dirs[[3]] && SoupX) {
      message("More than one data_dir provided. return_SoupX set to FALSE.")
    }
    return_SoupX <- checked_dirs[[3]]
    if (!SoupX) {
      return_SoupX <- F
    }
    ffbms <- unlist(lapply(data_dirs, function(x) x[which(grepl("filtered_feature_bc_matrix|filtered_gene_bc_matrices", x))]))
    rfbms <- unlist(lapply(data_dirs, function(x) x[which(grepl("raw_feature_bc_matrix|raw_gene_bc_matrices", x))]))
  } else {
    if (SoupX) {
      message("ffbms and/or rfbms provided directly. SoupX and returnSoupX set to FALSE.")
      SoupX <- F
    }
    return_SoupX <- F
    if (decontX) {
      message("ffbms and/or rfbms provided directly. decontX set to FALSE.")
      decontX <- F
    }
    if (!is.null(ffbms)) {
      if (is.null(names(ffbms))) {
        stop("ffbms need to have names.")
      }
    }
    if (!is.null(rfbms)) {
      if (is.null(names(rfbms))) {
        stop("rfbms need to have names.")
      }
    }
  }

  if (length(ffbms) == 1 && batch_corr == "harmony") {
    message("Only one sample provided. Setting batch_corr to 'none'.")
    batch_corr <- "none"
  }


  if (any(duplicated(names(ffbms)))) {
    print(ffbms)
    stop("Duplicate names for paths not allowed.")
  }

  message("Reading filtered_feature_bc_matrix data.")
  SO <- lapply(names(ffbms), function(x) {

    message(x)
    # this adds folder name prefix to cell names
    #names(x) <- basename(dirname(x))


    if (any(grepl("\\.h5$", list.files(ffbms[x])))) {
      if (length(grepl("\\.h5$", list.files(ffbms[x]))) > 1) {
        message("Found more than one .h5 file in ", ffbms[x], ". Will use the first: ", grep("\\.h5$", list.files(ffbms[x], full.names = T), value = T)[1])
      }
      filt_data <- Seurat::Read10X_h5(filename = grep("\\.h5$", list.files(ffbms[x], full.names = T)[1], value = T))
    } else {
      filt_data <- Seurat::Read10X(data.dir = ffbms[x])
    }

    if (is.list(filt_data)) {
      message("filtered_feature_bc_matrix is a list. Using 'Gene Expression' index")
      filt_data <- filt_data[["Gene Expression"]]
    }
    rownames(filt_data) <- gsub("_", "-", rownames(filt_data)) # handle error with scDblFinder below by manual correction of matrix, rather than have it corrected in SO only

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

    if (ncol(filt_data) == 0) {
      message("No cells left after filtering. Return NULL for this sample.")
      return(NULL)
      ## check that (giving names after loop.)
    }

    message("Creating initial Seurat object with ", ncol(filt_data), " cells.")

    ## feature_rm, feature_aggr
    if (!is.null(feature_rm) || !is.null(feature_aggr)) {
      if (!is.null(feature_aggr)) {
        if (!is.list(feature_aggr) || is.null(names(feature_aggr)) || anyDuplicated(names(feature_aggr))) {
          stop("feature_aggr has to be a named list. Each list entry should contain features to aggregate, names are new feature names and should be unique")
        }
        aggr_rows <- lapply(names(feature_aggr), function(x) {
          y <- SeuratObject::as.sparse(matrix(Matrix::colSums(filt_data[which(rownames(filt_data) %in% feature_aggr[[x]]),,drop=F]), nrow = 1))
          rownames(y) <- x
          return(y)
        })
        aggr_rows <- Reduce(rbind, aggr_rows)
        filt_data <- rbind(filt_data, aggr_rows)
      }

      if (!is.null(feature_rm)) {
        if (!is.character(feature_rm)) {
          stop("feature_rm has to be a character vector of features to remove.")
        }
        filt_data <- filt_data[which(!rownames(filt_data) %in% feature_rm),]
      }
    }

    # this is generally not a bad idea and it was necessary to get scDblFinder running once: https://github.com/LTLA/BiocNeighbors/issues/24
    if (!is.null(min_UMI)) {
      UMI_sum <- colSums(filt_data)
      if (any(UMI_sum < min_UMI)) {
        message(sum(UMI_sum < min_UMI), " cells removed for having less UMI then min_UMI.")
      }
      filt_data <- filt_data[,names(UMI_sum[which(UMI_sum >= min_UMI)])]
    }

    SO <- Seurat::CreateSeuratObject(counts = filt_data)
    SO@meta.data$orig.ident <- x
    Seurat::RenameCells(SO, paste0(Seurat::Cells(SO), "__", x))

    # filter cells which have library size of zero to enable scDblFinder without error
    # https://github.com/plger/scDblFinder/issues/55
    if (scDblFinder) {
      var_feat <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(SO,
                                                                        selection.method = "vst",
                                                                        nfeatures = nhvf,
                                                                        verbose = F,
                                                                        assay = "RNA"))
      factors <- scuttle::librarySizeFactors(filt_data[var_feat,])
      zero_libsize_cells <- names(which(factors == 0))
      if (length(zero_libsize_cells) > 0) {
        message(length(zero_libsize_cells), " cell(s) found which have zero library size based on hvf. These are removed to allow running scDblFinder. See https://github.com/plger/scDblFinder/issues/55.")
        SO <- subset(SO, cells = setdiff(names(factors), zero_libsize_cells))
      }

      if (!is.null(min_UMI_var_feat)) {
        UMI_sum <- colSums(Seurat::GetAssayData(SO, slot = "counts", assay = "RNA")[var_feat,])
        if (any(UMI_sum < min_UMI_var_feat)) {
          message(sum(UMI_sum < min_UMI_var_feat), " cells removed for having less UMI in variable features as min_UMI_var_feat. See https://github.com/LTLA/BiocNeighbors/issues/24,")
        }
        SO <- subset(SO, cells = names(UMI_sum[which(UMI_sum >= min_UMI_var_feat)]))
      }

      # scDblFinder
      SO@meta.data$dbl_score <- scDblFinder::computeDoubletDensity(x = Seurat::GetAssayData(SO, slot = "counts", assay = "RNA"),
                                                                   subset.row = var_feat,
                                                                   dims = npcs)
      SO@meta.data$dbl_score_log <- log1p(SO@meta.data$dbl_score)
    }
    return(SO)
  })
  names(SO) <- names(ffbms)

  if (length(SO) > 1) {
    message("Preparing merged and harmonized Seurat object with ", length(unlist(lapply(SO, Seurat::Cells))), " cells.")
  } else {
    message("Preparing Seurat object.")
  }

  SO <- prep_SO(SO_unprocessed = SO,
                reductions = "umap",
                nhvf = nhvf,
                npcs = npcs,
                batch_corr = batch_corr,
                RunHarmony_args = list(group.by.vars = "orig.ident"),
                FindClusters_args = list(resolution = resolution),
                normalization = "LogNormalize",
                diet_seurat = F)

  if (SoupX) {
    # use filt_data which may have been reduced 'cells' selection; raw_feature_bc_matrix will provide the whole picture of the soup
    message("Reading filtered and raw_feature_bc_matrix data.")

    SoupX_results <- lapply(names(rfbms), ..., FUN = function(x) {
      message(x)
      filt_data <- Seurat::Read10X(data.dir = ffbms[x])
      if (is.list(filt_data)) {
        filt_data <- filt_data[["Gene Expression"]]
      }
      ## filter for existing cells in SO (potential scDblFinder filtering, above)
      filt_data <- filt_data[,intersect(Seurat::Cells(SO), colnames(filt_data))]
      raw_data <- Seurat::Read10X(data.dir = rfbms[x])
      if (is.list(raw_data)) {
        raw_data <- raw_data[["Gene Expression"]]
      }

      sc <- SoupX::SoupChannel(tod = raw_data, toc = filt_data)

      ## https://github.com/constantAmateur/SoupX/issues/93
      ## run clustering on the specified subset only, otherwise an error may occur
      ## use intersect here as some cells may have been excluded above for scDblFinder
      SO_sub <-
        subset(SO, cells = intersect(Seurat::Cells(SO), rownames(sc$metaData))) %>%
        Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nhvf, verbose = F, assay = "RNA") %>%
        Seurat::ScaleData(verbose = F) %>%
        Seurat::RunPCA(verbose = F) %>%
        Seurat::FindNeighbors(verbose = F) %>%
        Seurat::FindClusters(verbose = F, algorithm = 1, resolution = resolution_SoupX)
      #Seurat::FetchData(vars = paste0("RNA_snn_res.",resolution_SoupX))

      #clusts <- Seurat::FindClusters(Seurat::FindNeighbors(SO@reductions[["pca"]]@cell.embeddings[rownames(sc$metaData),], verbose = F)$snn , algorithm = 1, resolution = resolution_SoupX, verbose = F)
      sc <- SoupX::setClusters(sc, stats::setNames(as.character(SO_sub@meta.data[,paste0("RNA_snn_res.",resolution_SoupX)]), rownames(SO_sub@meta.data)))

      ## old
      #temp_dots <- dots[which(grepl("^SoupX__", names(dots), ignore.case = T))]
      #names(temp_dots) <- gsub("^SoupX__", "", names(temp_dots), ignore.case = T)
      #temp_dots <- temp_dots[which(names(temp_dots) %in% formals(SoupX::autoEstCont))]
      #sc <- do.call(SoupX::autoEstCont, args = c(list(sc = sc, verbose = F), temp_dots))

      message("Running SoupX.")

      sc <- SoupX::autoEstCont(sc, verbose = F, ...)
      sx_counts <- suppressWarnings(SoupX::adjustCounts(sc, verbose = 0))

      if (return_SoupX) {
        message("Creating Seurat object on SoupX-corrected count matrix with ", ncol(filt_data), " cells.")
        SO_sx <- Seurat::CreateSeuratObject(counts = sx_counts)
        SO_sx@meta.data$orig.ident <- x

        SO_sx <- prep_SO(SO_unprocessed = SO_sx,
                         reductions = "umap",
                         nhvf = nhvf,
                         npcs = npcs,
                         batch_corr = batch_corr,
                         RunHarmony_args = list(group.by.vars = "orig.ident"),
                         FindClusters_args = list(resolution = resolution),
                         normalization = "LogNormalize",
                         diet_seurat = F)

        SO_sx <- Seurat::AddMetaData(SO_sx, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

        sc <- SoupX::setDR(sc, SO_sx@reductions$umap@cell.embeddings)

        sc_info_df <- data.frame(n_expr_uncorrected = Matrix::rowSums(sc$toc > 0),
                                 n_expr_corrected = Matrix::rowSums(sx_counts > 0)) %>%
          dplyr::mutate(abs_diff = n_expr_uncorrected-n_expr_corrected) %>%
          dplyr::mutate(rel_diff = abs_diff/n_expr_uncorrected) %>%
          dplyr::filter(abs_diff > 0) %>%
          tibble::rownames_to_column("Feature")

        message("Optionally: Create a SoupX RNA assay as follows: SO[['SoupXRNA']] <- Seurat::CreateAssayObject(counts = soupx_matrix).")
        return(list(SO = SO_sx,
                    sc = sc,
                    sc_info = sc_info_df,
                    pct_soup_SoupX = (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100))
      } else {
        return(list(pct_soup_SoupX = (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100))
      }
    })
    # add percentage of soup as meta data, similar to decontX
    SO <- Seurat::AddMetaData(SO, unlist(sapply(SoupX_results, "[", "pct_soup_SoupX")), "pct_soup_SoupX")
  }

  if (return_SoupX) {
    SO <- list(SO, SO_sx[["SO"]])
    names(SO) <- c("original", "SoupX")
  } else {
    SO <- list(SO)
    names(SO) <- "original"
  }

  SO <- lapply(SO, function(SOx) {

    if (!scDblFinder) {
      qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt")
    } else {
      qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt", "dbl_score")
    }

    # differentiate mouse, human or no MT-genes at all
    # and add freq of RPS / RPL and MRPS / MRPL genes
    if (any(grepl("^MT-", rownames(SOx))) && !any(grepl("^mt-", rownames(SOx)))) {
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^MT-"), "pct_mt")
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^RP[SL]"), "pct_ribo")
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^MRP[SL]"), "pct_mribo")
    } else if (!any(grepl("^MT-", rownames(SOx))) && any(grepl("^mt-", rownames(SOx)))) {
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^mt-"), "pct_mt")
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^Rp[sl]-"), "pct_ribo")
      SOx <- Seurat::AddMetaData(SOx, Seurat::PercentageFeatureSet(SOx, pattern = "^Mrp[sl]-"), "pct_mribo")
    } else {
      message("No mitochondrial genes could be identified from gene names - none starting with MT- (human) or mt- (mouse).")
      qc_cols <- qc_cols[-which(qc_cols == "pct_mt")]
    }
    qc_cols <- paste0(qc_cols, "_log")


    if (decontX) {
      message("Running decontX.")
      ## multi dirs: split matrix
      split_mats <- lapply(split(x = seq_len(ncol(Seurat::GetAssayData(SOx, slot = "counts"))), f = SOx@meta.data$orig.ident), function(ind) Seurat::GetAssayData(SOx, slot = "counts")[,ind , drop = FALSE])
      names(split_mats) <- basename(dirname(ffbms))
      split_idents <- split(SOx@meta.data[[paste0("RNA_snn_res.", resolution)]], SOx@meta.data$orig.ident)
      names(split_idents) <- basename(dirname(ffbms))

      dx <- lapply(names(split_mats), function(x) {
        message(x)
        celda::decontX(x = split_mats[[x]],
                       z = split_idents[[x]],
                       varGenes = nhvf)[["contamination"]]
      })
      SOx <- Seurat::AddMetaData(SOx, unlist(dx), "pct_soup_decontX")
      qc_cols <- c(qc_cols, "pct_soup_decontX")
    }
    if (SoupX) {
      qc_cols <- c(qc_cols, "pct_soup_SoupX")
    }

    SOx@meta.data$nFeature_RNA_log <- log1p(SOx@meta.data$nFeature_RNA)
    SOx@meta.data$nCount_RNA_log <- log1p(SOx@meta.data$nCount_RNA)
    SOx@meta.data$pct_mt_log <- log1p(SOx@meta.data$pct_mt)
    # this could be done above with SO
    SOx@meta.data$pct_mt <- ifelse(is.na(SOx@meta.data$pct_mt), 0,SOx@meta.data$pct_mt)
    SOx@meta.data$pct_mt_log <- ifelse(is.na(SOx@meta.data$pct_mt_log ), 0, SOx@meta.data$pct_mt_log )

    ## multi-dirs: split matrix!
    SOx@meta.data$residuals <- unlist(lapply(unique(SOx@meta.data$orig.ident), function(x) {
      stats::residuals(stats::lm(nCount_RNA_log~nFeature_RNA_log, data = SOx@meta.data[which(SOx@meta.data$orig.ident == x),]))
    }))

    ## clustering on meta data (quality metrics)
    message("Running dimension reduction and clustering on qc meta data.")
    meta <- dplyr::select(SOx@meta.data, dplyr::all_of(qc_cols), residuals)

    for (nn in n_PCs_to_meta_clustering) {
      if (nn > 0) {
        if (batch_corr == "harmony" && length(ffbms) > 1) {
          temp_slot <- "harmony"
        } else {
          temp_slot <- "pca"
        }
        meta2 <- scale_min_max(cbind(meta, SOx@reductions[[temp_slot]]@cell.embeddings[,1:nn])) # length(ffbms) or length(data_dirs)
      } else {
        meta2 <- scale_min_max(meta)
      }


      #https://datascience.stackexchange.com/questions/27726/when-to-use-cosine-simlarity-over-euclidean-similarity
      # with cosine metric: relative composition is more important
      umap_dims <- uwot::umap(X = meta2, metric = "cosine")
      colnames(umap_dims) <- c(paste0("meta_UMAP_1_PC", nn), paste0("meta_UMAP_2_PC", nn))
      SOx[[paste0("umapmetaPC", nn)]] <- Seurat::CreateDimReducObject(embeddings = umap_dims, key = paste0("UMAPMETAPC", nn, "_"), assay = "RNA")

      # fix colnames manually as Seurat::CreateDimReducObject makes a mistake in taking the trailing number of umap_dims-colnames as trailing dim-number
      colnames(SOx@reductions[[paste0("umapmetaPC", nn)]]@cell.embeddings) <- paste0(toupper(paste0("umapmetaPC", nn, "_")), c(1,2))

      # matrix has to be supplied to FindNeighbors
      clusters <- Seurat::FindClusters(Seurat::FindNeighbors(meta2, annoy.metric = "cosine", verbose = F)$snn,
                                       resolution = resolution_meta,
                                       verbose = F)
      colnames(clusters) <- paste0("meta_", colnames(clusters), "_PC", nn)
      SOx <- Seurat::AddMetaData(SOx, cbind(umap_dims, clusters))
    }

    return(SOx)
  })

  message("Run scexpr:::qc_plots() on the Seurat object for visualization of clustering on phenotype and qc data.")

  if (return_SoupX) {
    message("Use e.g. SoupX::plotChangeMap(x[['sc']], cleanedMatrix = SoupX::adjustCounts(x[['sc']]), geneSet = 'GNLY') + theme_bw() + theme(panel.grid = element_blank()) + scale_color_gradientn(colors = col_pal('spectral'), na.value = 'grey95') to plot changes in counts.")
    SO <- c(SO, list(SO_sx[["sc"]]), list(SO_sx[["sc_info"]]))
    names(SO) <- c(names(SO), "sc", "sc_info")
    return(SO)
  } else {
    return(SO)
  }
}

check_dir <- function(data_dirs, SoupX = F) {

  dir_roots <- unlist(lapply(data_dirs, function(x) {
    dd <- list.dirs(x)
    dirname(dd[which(grepl("filtered_feature_bc_matrix|filtered_gene_bc_matrices", basename(dd)))])
  }))

  if (is.null(names(dir_roots))) {
    names(dir_roots) <- basename(dir_roots)
  }

  if (any(duplicated(names(dir_roots)))) {
    dup_inds <- which(duplicated(basename(dir_roots)))
    dup_names <- basename(dir_roots)[dup_inds]
    stop("Duplicated names found: ", paste(dup_names, collapse = ", "), ". Please fix.")
  }
  message(length(dir_roots), " folders with filtered_feature_bc_matrix or filtered_gene_bc_matrices found in data_dirs.")

  raw_and_filt <- unlist(lapply(dir_roots, function(x) {
    sub_folders <- basename(list.dirs(x, recursive = F))
    all(c("filtered_feature_bc_matrix", "raw_feature_bc_matrix") %in% sub_folders) | all(c("filtered_gene_bc_matrices", "raw_gene_bc_matrices") %in% sub_folders)
  }))

  if (any(!raw_and_filt) && SoupX) {
    message(length(raw_and_filt) - sum(raw_and_filt), " of ", length(raw_and_filt), " folder(s) (", paste(names(dir_roots[which(!raw_and_filt)]), collapse = ", "), ") do not contain raw_feature_bc_matrix. Hence, SoupX cannot be run.")
    SoupX <- F
  }

  dir_folders <- lapply(dir_roots, function(x) {
    dir_temp <- list.dirs(x, recursive = F)
    dir_temp[which(grepl("filtered_feature_bc_matrix|filtered_gene_bc_matrices|raw_feature_bc_matrix|raw_gene_bc_matrices", basename(dir_temp)))]
  })

  # length(dir_roots) == 1 sets return_SoupX
  return(list(dir_folders, SoupX, length(dir_roots) == 1))
}



qc_plots <- function(SO,
                     qc_cols = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log", "residuals"),
                     clustering_cols = c("RNA_snn_res.0.8", "meta_res.0.8_PC2"),
                     reduction = "umapmetaPC2",
                     geom2 = "boxplot") {


  # extra to do: nCount_RNA_log, nFeature_RNA_log, pct_mt_log - check and calc if needed

  if (any(!qc_cols %in% names(SO@meta.data))) {
    message(paste(qc_cols[which(!qc_cols %in% names(SO@meta.data))], collapse = ", "), " not found in SO@meta.data. These will be excluded from plotting.")
    qc_cols <- intersect(qc_cols, names(SO@meta.data))
    if (length(qc_cols) == 0) {
      stop("No qc_cols left for plotting.")
    }
  }

  if (length(clustering_cols) != 2) {
    stop("clustering_cols should be of length 2 containing one resolution for clustering that has been conducted on feature expression (index 1) and one clustering conducted on qc meta data (index 2).")
  }
  if (any(!clustering_cols %in% names(SO@meta.data))) {
    stop("At least one of clustering_cols has not been found in SO@meta.data. Check and fix that.")
  }

  breaks <- c(seq(0, 1e1, 2e0),
              seq(0, 1e2, 2e1),
              seq(0, 1e3, 2e2),
              seq(0, 1e4, 2e3),
              seq(0, 1e5, 2e4))
  breaks <- breaks[which(breaks != 0)]

  qc_p1 <- suppressMessages(feature_plot(SO,
                                         features = c(qc_cols, clustering_cols[1]),
                                         reduction = "UMAP", legend.position = "none",
                                         plot.labels = T))

  qc_p2 <- ggplot2::ggplot(tidyr::pivot_longer(SO@meta.data[,c(qc_cols, clustering_cols)], cols = dplyr::all_of(qc_cols), names_to = "qc_param", values_to = "value"),
                           ggplot2::aes(x = !!rlang::sym(clustering_cols[1]), y = value, color = !!rlang::sym(clustering_cols[2]))) +
    ggplot2::geom_boxplot(color = "grey30", outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, size = 0.3) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(0.3, "cm"), legend.key = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = col_pal("custom")) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ expm1(.), breaks = breaks)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
    ggplot2::facet_wrap(ggplot2::vars(qc_param), scales = "free_y", ncol = 1)

  p3_1 <- patchwork::wrap_plots(feature_plot(SO,
                                             features = clustering_cols[2],
                                             reduction = reduction,
                                             pt.size = 0.5,
                                             col.pal.dir = 1,
                                             legend.position = "none",
                                             plot.labels = T,
                                             plot.title = F),
                                suppressMessages(freq_pie_chart(SO = SO, meta.col = clustering_cols[2])[["plot"]]),
                                ncol = 1)

  p3_2 <- feature_plot_stat(SO,
                            features = "nCount_RNA_log",
                            meta.col = clustering_cols[2],
                            geom2 = geom2,
                            jitterwidth = 0.9,
                            panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                            axis.title.y = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ expm1(.), breaks = breaks[intersect(which(breaks > min(expm1(SO@meta.data$nCount_RNA_log))),
                                                                                                   which(breaks < max(expm1(SO@meta.data$nCount_RNA_log))))]))

  p3_3 <- feature_plot_stat(SO,
                            features = "nFeature_RNA_log",
                            meta.col = clustering_cols[2],
                            geom2 = geom2,
                            jitterwidth = 0.9,
                            panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                            axis.title.y = ggplot2::element_blank(),
                            axis.text.x = ggplot2::element_blank(),
                            axis.title.x = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ expm1(.), breaks = breaks[intersect(which(breaks > min(expm1(SO@meta.data$nFeature_RNA_log))),
                                                                                                   which(breaks < max(expm1(SO@meta.data$nFeature_RNA_log))))]))

  p3_4 <- feature_plot_stat(SO,
                            features = "pct_mt_log",
                            meta.col = clustering_cols[2],
                            geom2 = geom2,
                            jitterwidth = 0.9,
                            panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                            axis.title.y = ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ expm1(.), breaks = breaks[intersect(which(breaks > min(expm1(SO@meta.data$pct_mt_log))),
                                                                                                   which(breaks < max(expm1(SO@meta.data$pct_mt_log))))]))

  p3_x <- lapply(qc_cols[which(!grepl("nCount_RNA_log|nFeature_RNA_log|pct_mt_log", qc_cols))], function(qcf) {
    p3_x <- feature_plot_stat(SO,
                              features = qcf,
                              meta.col = clustering_cols[2],
                              geom2 = geom2,
                              jitterwidth = 0.9,
                              panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                              axis.title.y = ggplot2::element_blank())
    if (qcf != rev(qc_cols[which(!grepl("nCount_RNA_log|nFeature_RNA_log|pct_mt_log", qc_cols))])[1]) {
      p3_x <-
        p3_x +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())
    }
    return(p3_x)
  })


  p3_2 <- patchwork::wrap_plots(p3_2, p3_3, p3_4, ncol = 1)
  p3_3 <- patchwork::wrap_plots(p3_x, ncol = 1)
  qc_p3 <- patchwork::wrap_plots(p3_1, p3_2, p3_3, nrow = 1)

  return(list(qc_p1 = qc_p1, qc_p2 = qc_p2, qc_p3 = qc_p3))

}

ceiling_any = function(x, accuracy, f = ceiling) {
  f(x/ accuracy) * accuracy
}

floor_any = function(x, accuracy, f = floor) {
  f(x/ accuracy) * accuracy
}


'
  # check data_dirs
  if (class(data_dirs) == "list") {
    if (length(data_dirs) > 1) {
      stop("Only provide one data_dirs which may contain two folders, filtered_feature_bc_matrix and raw_feature_bc_matrix.")
    }
    if (!is.null(names(data_dirs))) {
      name <- names(data_dirs)
    } else {
      if (length(unique(dirname(data_dirs[[1]]))) != 1) {
        stop("data_dirss should have common parent dir, otherwise provide a named list of filtered_feature_bc_matrix and raw_feature_bc_matrix.")
      }
      name <- basename(unique(dirname(data_dirs[[1]])))
    }
    data_dirs <- stats::setNames(data_dirs[[1]], nm = rep(name, length(data_dirs[[1]])))
  } else {
    if (is.null(names(data_dirs))) {
      if (length(unique(dirname(data_dirs))) != 1) {
        stop("data_dirss should have common parent dir, otherwise provide a named list of filtered_feature_bc_matrix and raw_feature_bc_matrix.")
      }
      names(data_dirs) <- basename(unique(dirname(data_dirs)))
    }
  }

  if (length(data_dirs) == 1) {
    if (basename(data_dirs) != "filtered_feature_bc_matrix") {
      stop("data_dirs is expected to contain full paths of two dirs, filtered_feature_bc_matrix and raw_feature_bc_matrix, or filtered_feature_bc_matrix only.")
    }
  } else if (length(data_dirs) == 2) {
    if (length(intersect(basename(data_dirs), c("filtered_feature_bc_matrix", "raw_feature_bc_matrix"))) < 2) {
      stop("data_dirs is expected to contain full paths of two dirs, filtered_feature_bc_matrix and raw_feature_bc_matrix, or filtered_feature_bc_matrix only.")
    }
  }

  if (any(basename(data_dirs) == data_dirs)) {
    stop("Please provide full path to filtered_feature_bc_matrix and/or raw_feature_bc_matrix.")
  }

  if (SoupX && length(data_dirs) != 2) {
    warning("SoupX requires filtered_feature_bc_matrix and raw_feature_bc_matrix. SoupX will not run.")
    SoupX <- F
  }
'

'
    ### old soupX
    raw_data <- Seurat::Read10X(data.dir = grep("raw_feature_bc_matrix", data_dirs, value = T))
    if (is.list(raw_data)) {
      raw_data <- raw_data[["Gene Expression"]]
    }
    message("Running SoupX.")
    sc <- SoupX::SoupChannel(tod = raw_data, toc = filt_data)
    if (resolution_SoupX != resolution) {
      SO <- Seurat::FindClusters(SO, algorithm = 1, resolution = resolution_SoupX, verbose = F)
    }
    sc = SoupX::setClusters(sc, SO@meta.data[rownames(sc$metaData), paste0("RNA_snn_res.", resolution_SoupX)])

    #temp_dots <- dots[which(grepl("^SoupX__", names(dots), ignore.case = T))]
    #names(temp_dots) <- gsub("^SoupX__", "", names(temp_dots), ignore.case = T)
    #temp_dots <- temp_dots[which(names(temp_dots) %in% formals(SoupX::autoEstCont))]
    #sc <- do.call(SoupX::autoEstCont, args = c(list(sc = sc, verbose = F), temp_dots))

    sc <- SoupX::autoEstCont(sc, verbose = F, ...)
    sx_counts = SoupX::adjustCounts(sc, verbose = 0)

    # add percentage of soup as meta data, similar to decontX
    SO <- Seurat::AddMetaData(SO, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

    if (return_SoupX) {
      message("Creating Seurat object on SoupX-corrected count matrix with ", ncol(filt_data), " cells.")
      SO_sx <-
        Seurat::CreateSeuratObject(counts = sx_counts, verbose = F) %>%
        Seurat::NormalizeData(verbose = F) %>%
        Seurat::FindVariableFeatures(nfeatures = nhvf, verbose = F) %>%
        Seurat::ScaleData(verbose = F) %>%
        Seurat::RunPCA(npcs = npcs, verbose = F) %>%
        Seurat::RunUMAP(dims = 1:npcs, umap.method = "uwot", metric = "cosine", verbose = F) %>%
        Seurat::FindNeighbors(dims = 1:npcs, verbose = F) %>%
        Seurat::FindClusters(algorithm = 1, resolution = resolution, verbose = F)
      SO_sx@meta.data$orig.ident <- names(data_dirs)[1]

      SO_sx <- Seurat::AddMetaData(SO_sx, (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100, "pct_soup_SoupX")

      sc <- SoupX::setDR(sc, SO_sx@reductions$umap@cell.embeddings)

      sc_info_df <- data.frame(n_expr_uncorrected = Matrix::rowSums(sc$toc > 0),
                               n_expr_corrected = Matrix::rowSums(sx_counts > 0)) %>%
        dplyr::mutate(abs_diff = n_expr_uncorrected-n_expr_corrected) %>%
        dplyr::mutate(rel_diff = abs_diff/n_expr_uncorrected) %>%
        dplyr::filter(abs_diff > 0) %>%
        tibble::rownames_to_column("Feature")

      message("Optionally: Create a SoupX RNA assay as follows: SO[["SoupXRNA"]] <- Seurat::CreateAssayObject(counts = soupx_matrix).")
      SO <- list(SO, SO_sx)
      names(SO) <- c("original", "SoupX")
    } else {
      SO <- list(SO)
      names(SO) <- "original"
    }'


'
        SO_sx <-
          Seurat::CreateSeuratObject(counts = sx_counts, verbose = F) %>%
          Seurat::NormalizeData(verbose = F) %>%
          Seurat::FindVariableFeatures(nfeatures = nhvf, verbose = F) %>%
          Seurat::ScaleData(verbose = F) %>%
          Seurat::RunPCA(npcs = npcs, verbose = F) %>%
          Seurat::RunUMAP(dims = 1:npcs, umap.method = "uwot", metric = "cosine", verbose = F) %>%
          Seurat::FindNeighbors(dims = 1:npcs, verbose = F) %>%
          Seurat::FindClusters(algorithm = 1, resolution = resolution, verbose = F)
        SO_sx@meta.data$orig.ident <- names(data_dirs)[1]
        '
## old scDblFinder
'    if (scDblFinder) {
      message("Running scDblFinder.")

      ## https://stackoverflow.com/questions/62161916/is-there-a-function-in-r-that-splits-a-matrix-along-a-margin-using-a-factor-or-c
      ## modified from base function split.data.frame (to avoid 2 x transposation)
      ## multi dirs: split count matrix by orig.idents
      split_mats <- lapply(split(x = seq_len(ncol(Seurat::GetAssayData(SOx, slot = "counts"))), f = SOx@meta.data$orig.ident), function(ind) Seurat::GetAssayData(SOx, slot = "counts")[,ind , drop = FALSE])
      names(split_mats) <- basename(dirname(ffbms))

      ## use for-loop to allow break-statement
      for (x in names(split_mats)) {
        message(x)
        dbl_score_temp <- tryCatch({
          log1p(computeDoubletDensity(x = split_mats[[x]],
                                      subset.row = Seurat::VariableFeatures(Seurat::FindVariableFeatures(subset(SOx, cells = colnames(split_mats[[x]])),
                                                                                                         selection.method = "vst",
                                                                                                         nfeatures = nhvf,
                                                                                                         verbose = F,
                                                                                                         assay = "RNA")),
                                      dims = npcs))
        }, error = function(error_condition) {
          message(error_condition, " ... doublet calculation failed. Try to increase nhvf. Seurat object is exported as global variable: SO_qc_export_rescue")
          message("")
          SO_qc_export_rescue <<- SO
          return(NULL)
        })

        if (is.null(dbl_score_temp)) {
          ## in case of error an any data set: dbl_score set NULL and hence excluded (see below)
          dbl_score <- NULL
          break
        } else {
          dbl_score <- c(dbl_score, dbl_score_temp)
        }
      }
    }'
