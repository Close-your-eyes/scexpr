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
#' the common Seurat object though, the merged and harmonized PCA space of which is subject for clustering the cells based on feature expression (phenotypes)
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
#' @param resolution_meta resolution(s) (louvain algorithm) for clustering based on qc meta data and optionally additional PC dimensions (PCs_to_meta_clustering)
#' @param PCs_to_meta_clustering how many principle components (PCs) from phenotypic clustering to add to qc meta data;
#' this will generate a mixed clustering (PCs from phenotypes (RNA) and qc meta data like pct mt and nCount_RNA); the more PCs are added the greater the
#' phenotypic influence becomes; one or more integers can be supplied to explore the effect; pass 0, to have no PCs included in meta clustering; e.g. when
#' PCs_to_meta_clustering = 3 PCs 1-3 are used, when PCs_to_meta_clustering = 1 only PC 1 is used.
#' @param scDblFinder logical, whether to run doublet detection algorithm from scDblFinder
#' @param SoupX_return logical whether to return a full Seurat-object and diagnostics from SoupX (TRUE) or whether to run SoupX without these returns and just
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
#' @param diet_seurat
#'
#' @return a list of Seurat object and data frame with marker genes for clusters based on feature expression
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
SO_prep01 <- function(data_dirs,
                      nhvf = 2000,
                      npcs = 10,
                      resolution = 0.8,
                      resolution_SoupX = 0.6,
                      resolution_meta = seq(0.1,0.8,0.1),
                      PCs_to_meta_clustering = 2,
                      scDblFinder = F,
                      min_UMI = 30,
                      min_UMI_var_feat = 10,
                      SoupX = F,
                      decontX = F,
                      SoupX_return = F,
                      SoupX_autoEstCont_args = list(),
                      cells = NULL,
                      invert_cells = F,
                      feature_rm = NULL,
                      feature_aggr = NULL,
                      ffbms = NULL,
                      rfbms = NULL,
                      batch_corr = c("harmony", "none"),
                      diet_seurat = T) {


  install_pkgs(SoupX, scDblFinder, decontX)
  resolution <- checks(resolution_meta, resolution_SoupX, resolution, PCs_to_meta_clustering)
  batch_corr <- rlang::arg_match(batch_corr)

  c(ffbms,
    rfbms,
    SoupX,
    SoupX_return,
    decontX,
    batch_corr) %<-% check_inputs(ffbms = ffbms,
                                  rfbms = rfbms,
                                  data_dirs = data_dirs,
                                  SoupX = SoupX,
                                  SoupX_return = SoupX_return,
                                  decontX = decontX,
                                  batch_corr = batch_corr)
  #SoupX_return <- T

  message("Reading filtered_feature_bc_matrix data.")
  SO <- purrr::map(stats::setNames(names(ffbms), names(ffbms)), function(x) {

    message(x)
    # this adds folder name prefix to cell names
    #names(x) <- basename(dirname(x))

    # not SO yet but matrix only
    SO <- read_10X_data(path = ffbms[x],
                        cells = cells,
                        min_UMI = min_UMI,
                        name = x)

    if (!is.null(feature_rm) || !is.null(feature_aggr)) {
      SO <- aggregate_or_remove_features(filt_data = SO,
                                         feature_rm = feature_rm,
                                         feature_aggr = feature_aggr)
    }

    message("Creating initial Seurat object with ", ncol(SO), " cells.")
    SO <- Seurat::CreateSeuratObject(counts = SO)
    SO@meta.data$orig.ident <- x

    # doublet score calculation with not yet merged data
    if (scDblFinder) {
      SO <- add_dbl_score_to_metadata(SO = SO,
                                      nhvf = nhvf,
                                      min_UMI_var_feat = min_UMI_var_feat,
                                      npcs = npcs)
    }
    return(SO)
  })

  if (length(SO) > 1) {
    message("Preparing merged and harmonized Seurat object with ", sum(lengths(lapply(SO, Seurat::Cells))), " cells.")
  } else {
    message("Preparing Seurat object.")
  }

  SO <- SO_prep02(SO_unprocessed = SO,
                  reductions = "umap",
                  nhvf = nhvf,
                  npcs = npcs,
                  min_cells = 1,
                  batch_corr = batch_corr,
                  RunHarmony_args = list(group.by.vars = "orig.ident"),
                  FindClusters_args = list(resolution = resolution),
                  normalization = "LogNormalize",
                  diet_seurat = F)
  suppressWarnings(SeuratObject::Misc(SO, "clusterings") <- paste0("RNA_snn_res.", resolution))

  if (SoupX) {
    # ffmbs and rfbms are paired by name
    SoupX_results <- run_soupx(ffbms = ffbms,
                               rfbms = rfbms,
                               SO = SO,
                               nhvf = nhvf,
                               min_UMI = min_UMI,
                               resolution_SoupX = resolution_SoupX,
                               npcs = npcs,
                               batch_corr = batch_corr,
                               resolution = resolution,
                               SoupX_return = SoupX_return,
                               SoupX_autoEstCont_args = SoupX_autoEstCont_args,
                               feature_aggr = feature_aggr,
                               feature_rm = feature_rm)

    # add percentage of soup as meta data, similar to decontX

    if (SoupX_return) {
      pct_soup_SoupX <- unlist(unname(sapply(SoupX_results[["SoupX_results"]], "[", "pct_soup_SoupX")))
    } else {
      pct_soup_SoupX <- unlist(unname(sapply(SoupX_results, "[", "pct_soup_SoupX")))
    }
    SO <- Seurat::AddMetaData(object = SO,
                              metadata = pct_soup_SoupX,
                              col.name = "pct_soup_SoupX")

  }

  # two SO if SoupX is to be returned
  if (SoupX && SoupX_return) {
    SO <- list(SO, SoupX_results[[1]])
    names(SO) <- c("original", "SoupX")
    SoupX_results <- SoupX_results[[2]]
  } else {
    SO <- list(SO)
    names(SO) <- "original"
  }

  if (decontX) {
    SO <- run_decontx(SO = SO,
                      resolution = resolution,
                      nhvf = nhvf)
  }

  SO <- cluster_on_metadata(SO = SO,
                            batch_corr = batch_corr,
                            PCs_to_meta_clustering = PCs_to_meta_clustering,
                            ffbms = ffbms,
                            resolution_meta = resolution_meta)


  if (diet_seurat) {
    SO <- lapply(SO, function(SO) Seurat::DietSeurat(SO, dimreducs = c(names(SO@reductions))))
  }

  message("Run scexpr::qc_plot2() on the Seurat object for visualization of clustering on phenotype and qc data.")
  if (SoupX && SoupX_return) {
    message("Use e.g. SoupX::plotChangeMap(x[['sc']], cleanedMatrix = SoupX::adjustCounts(x[['sc']]), geneSet = 'GNLY')")
    return(c(SO, list(SoupX_results = SoupX_results)))
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
  message(length(dir_roots), " folders with filtered_feature_bc_matrix and/or filtered_gene_bc_matrices found in data_dirs.")

  message("searching for partial matches with filtered_feature_bc_matrix|raw_feature_bc_matrix or filtered_gene_bc_matrices|raw_gene_bc_matrices.")
  raw_and_filt <- unlist(lapply(dir_roots, function(x) {
    sub_folders <- basename(list.dirs(x, recursive = F))
    return(sum(grepl("filtered_feature_bc_matrix|raw_feature_bc_matrix", sub_folders)) == 2 ||
             sum(grepl("filtered_gene_bc_matrices|raw_gene_bc_matrices", sub_folders)) == 2)
  }))

  if (any(!raw_and_filt) && SoupX) {
    message(length(raw_and_filt) - sum(raw_and_filt), " of ", length(raw_and_filt), " folder(s) (", paste(names(dir_roots[which(!raw_and_filt)]), collapse = ", "), ") do not contain raw_feature_bc_matrix. Hence, SoupX cannot be run and is set to FALSE.")
    SoupX <- F
  }

  dir_folders <- lapply(dir_roots, function(x) {
    dir_temp <- list.dirs(x, recursive = F)
    dir_temp[which(grepl("filtered_feature_bc_matrix|filtered_gene_bc_matrices|raw_feature_bc_matrix|raw_gene_bc_matrices", basename(dir_temp)))]
  })

  # length(dir_roots) == 1 sets SoupX_return
  return(list(dir_folders, SoupX, length(dir_roots) == 1))
}




install_pkgs <- function(SoupX, scDblFinder, decontX) {
  for (p in c("matrixStats", "hdf5r", "uwot", "devtools", "patchwork", "BiocManager")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      utils::install.packages(p)
    }
  }

  if (SoupX && !requireNamespace("SoupX", quietly = T)) {
    utils::install.packages("SoupX")
  }
  if (scDblFinder && !requireNamespace("scDblFinder", quietly = T)) {
    BiocManager::install("scDblFinder")
  }
  if (decontX && !requireNamespace("celda", quietly = T)) {
    BiocManager::install("celda")
  }
  if (!requireNamespace("scuttle", quietly = T)) {
    BiocManager::install("scuttle")
  }

  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github("immunogenomics/presto")
  }
}

checks <- function(resolution_meta, resolution_SoupX, resolution, PCs_to_meta_clustering) {
  if (!is.numeric(resolution_meta)) {
    stop("resolution_meta has to be numeric.")
  }
  if (!is.numeric(resolution_SoupX) || length(resolution_SoupX) != 1) {
    stop("resolution_SoupX has to be numeric and of length 1.")
  }
  if (!is.numeric(resolution)) {
    stop("resolution has to be numeric.")
  }
  resolution <- as.numeric(gsub("^1.0$", "1", resolution))
  # if (!is.integer(resolution) && dplyr::near(resolution, 1)) {
  #   resolution <- as.integer(resolution)
  # }
  if (!is.numeric(PCs_to_meta_clustering)) {
    stop("PCs_to_meta_clustering should be numeric.")
  }
  return(resolution)
}

add_dbl_score_to_metadata <- function(SO, nhvf, min_UMI_var_feat, npcs) {
  var_feat <- Seurat::VariableFeatures(Seurat::FindVariableFeatures(SO,
                                                                    selection.method = "vst",
                                                                    nfeatures = nhvf,
                                                                    verbose = F,
                                                                    assay = "RNA"))

  # filter out cells which have library size of zero to enable scDblFinder without error
  # https://github.com/plger/scDblFinder/issues/55
  factors <- scuttle::librarySizeFactors(SeuratObject::LayerData(SO, layer = "counts", assay = "RNA")[var_feat,])
  zero_libsize_cells <- names(which(factors == 0))
  if (length(zero_libsize_cells) > 0) {
    message(length(zero_libsize_cells), " cell(s) found which have zero library size based on hvf. These are removed to allow running scDblFinder. See https://github.com/plger/scDblFinder/issues/55.")
    SO <- subset(SO, cells = setdiff(names(factors), zero_libsize_cells))
  }

  if (!is.null(min_UMI_var_feat)) {
    UMI_sum <- Matrix::colSums(SeuratObject::LayerData(SO, layer = "counts", assay = "RNA")[var_feat,])
    if (any(UMI_sum < min_UMI_var_feat)) {
      message(sum(UMI_sum < min_UMI_var_feat), " cells removed for having less UMI in variable features as min_UMI_var_feat. See https://github.com/LTLA/BiocNeighbors/issues/24.")
    }
    SO <- subset(SO, cells = names(UMI_sum[which(UMI_sum >= min_UMI_var_feat)]))
  }
  # scDblFinder
  SO@meta.data$dbl_score <- scDblFinder::computeDoubletDensity(x = SeuratObject::LayerData(SO, layer = "counts", assay = "RNA"),
                                                               subset.row = var_feat,
                                                               dims = npcs)
  SO@meta.data$dbl_score_log <- log1p(SO@meta.data$dbl_score)
  return(SO)
}

read_10X_data <- function(path, cells, min_UMI, verbose = T, name) {

  h5files <- list.files(path, pattern = "\\.h5$", full.names = T)
  if (length(h5files)) {
    if (length(h5files) > 1 && verbose) {
      message("Found more than one .h5 file in ", path, ". Will use the first: ", h5files[1])
    }
    filt_data <- Seurat::Read10X_h5(filename = h5files[1])
  } else {
    filt_data <- Seurat::Read10X(data.dir = path)
  }

  if (is.list(filt_data)) {
    if (verbose) message("filtered_feature_bc_matrix is a list. Using 'Gene Expression' index")
    filt_data <- filt_data[["Gene Expression"]]
  }
  rownames(filt_data) <- gsub("_", "-", rownames(filt_data)) # handle error with scDblFinder below by manual correction of matrix, rather than have it corrected in SO only


  if (!is.null(cells)) {
    if (any(cells %in% colnames(filt_data))) {
      if (verbose) {
        message(length(cells), " cells provided.")
        message(length(which(cells %in% colnames(filt_data))), " of cells from a total of ", ncol(filt_data), " cells found in data (", round(length(which(cells %in% colnames(filt_data)))/ncol(filt_data)*100, 1), " %).")
      }
      cells <- cells[which(cells %in% colnames(filt_data))]
      if (invert_cells) {
        filt_data <- filt_data[,which(!colnames(filt_data) %in% cells)]
      } else {
        filt_data <- filt_data[,cells]
      }
    } else {
      if (verbose) message("Non of cells found in data.")
    }
  }
  if (ncol(filt_data) == 0) {
    if (verbose) message("No cells left after filtering for cells. Return NULL for this sample.")
    return(NULL)
    ## check that (giving names after loop.)
  }

  # this is generally not a bad idea and it was necessary to get scDblFinder running once: https://github.com/LTLA/BiocNeighbors/issues/24
  if (!is.null(min_UMI)) {
    UMI_sum <- Matrix::colSums(filt_data)
    if (any(UMI_sum < min_UMI) && verbose) {
      message(sum(UMI_sum < min_UMI), " cells removed for having less UMI then min_UMI.")
    }
    filt_data <- filt_data[,names(UMI_sum[which(UMI_sum >= min_UMI)])]
  }

  if (ncol(filt_data) == 0) {
    if (verbose) message("No cells left after filtering for min_UMI. Return NULL for this sample.")
    return(NULL)
    ## check that (giving names after loop.)
  }

  # change cell names here to avoid duplicate names from multiple samples
  # but only if not already there
  if (!all(grepl(paste0("^", name), colnames(filt_data)))) {
    colnames(filt_data) <- paste0(name, "__", colnames(filt_data))
  }


  return(filt_data)
}

aggregate_or_remove_features <- function(filt_data,
                                         feature_rm,
                                         feature_aggr) {
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
  return(filt_data)
}

check_inputs <- function(ffbms, rfbms, data_dirs, SoupX, SoupX_return, decontX, batch_corr) {

  if (is.null(ffbms) && is.null(rfbms)) {

    checked_dirs <- check_dir(data_dirs = data_dirs, SoupX = SoupX)
    data_dirs <- checked_dirs[[1]]

    if (SoupX && !checked_dirs[[2]]) {
      message("raw_feature_bc_matrix not found in every data_dir. SoupX set to FALSE.")
    }
    SoupX <- checked_dirs[[2]]
    # more than one input dir works now
    # if (SoupX_return && !checked_dirs[[3]] && SoupX) {
    #   message("More than one data_dir provided. SoupX_return set to FALSE.")
    # }
    # SoupX_return <- checked_dirs[[3]]

    if (!SoupX) {
      SoupX_return <- F
    }
    ffbms <- unlist(lapply(data_dirs, function(x) x[which(grepl("filtered_feature_bc_matrix|filtered_gene_bc_matrices", x))]))
    rfbms <- unlist(lapply(data_dirs, function(x) x[which(grepl("raw_feature_bc_matrix|raw_gene_bc_matrices", x))]))


  } else {


    if (SoupX) {
      message("ffbms and/or rfbms provided directly. SoupX and returnSoupX set to FALSE.")
      SoupX <- F
    }
    SoupX_return <- F
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

  if (length(ffbms) == 1 && batch_corr != "none") {
    message("Only one sample provided. Setting batch_corr to 'none'.")
    batch_corr <- "none"
  }

  if (anyDuplicated(names(ffbms))) {
    print(ffbms)
    stop("Duplicate names for paths not allowed.")
  }

  return(list(ffbms, rfbms, SoupX, SoupX_return, decontX, batch_corr))
}

run_soupx <- function(ffbms,
                      rfbms,
                      SO,
                      nhvf,
                      min_UMI,
                      resolution_SoupX,
                      npcs,
                      batch_corr,
                      resolution,
                      SoupX_return,
                      SoupX_autoEstCont_args,
                      feature_aggr,
                      feature_rm) {
  # ffmbs and rfbms are paired by name

  # use filt_data which may have been reduced 'cells' selection; raw_feature_bc_matrix will provide the whole picture of the soup
  message("SoupX: Reading filtered and raw_feature_bc_matrix data.")

  SoupX_results <- lapply(stats::setNames(names(rfbms), names(rfbms)), function(x) {
    message(x)
    # just read again
    filt_data <- read_10X_data(path = ffbms[x],
                               cells = Seurat::Cells(SO), ## filter for existing cells in SO (potential scDblFinder filtering, above)
                               min_UMI = min_UMI,
                               verbose = F,
                               name = x)
    raw_data <- read_10X_data(path = rfbms[x],
                              cells = NULL,
                              min_UMI = NULL,
                              verbose = F,
                              name = x)

    # will never be many features, so irrelevant to remove them for soup estimation
    if (!is.null(feature_rm) || !is.null(feature_aggr)) {
      filt_data <- aggregate_or_remove_features(filt_data = filt_data,
                                                feature_rm = feature_rm,
                                                feature_aggr = feature_aggr)
      raw_data <- aggregate_or_remove_features(filt_data = raw_data,
                                               feature_rm = feature_rm,
                                               feature_aggr = feature_aggr)
    }

    filt_data <- filt_data[order(rownames(filt_data)),]
    raw_data <- raw_data[order(rownames(raw_data)),]

    # maybe this is a property of CellRanger8 or of Fixed RNA profiling
    if (!identical(rownames(filt_data), rownames(raw_data))) {
      message("SoupX: feature names of raw_data and filt_data are not identical.")
      if (all(rownames(filt_data) %in% rownames(raw_data))) {
        message("However, all features from filt_data are in raw_data. Running SoupX on features from filt_data only.")
      }
      raw_data <- raw_data[rownames(filt_data),]
    }

    sc <- SoupX::SoupChannel(tod = raw_data, toc = filt_data)

    ## https://github.com/constantAmateur/SoupX/issues/93
    ## run clustering on the specified subset only, otherwise an error may occur
    ## use intersect here as some cells may have been excluded above for scDblFinder
    clusters <-
      subset(Seurat::DietSeurat(SO, layers = "counts"), cells = intersect(Seurat::Cells(SO), rownames(sc$metaData))) |>
      Seurat::NormalizeData(verbose = F) |>
      Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nhvf, verbose = F, assay = "RNA") |>
      Seurat::ScaleData(verbose = F) |>
      Seurat::RunPCA(verbose = F) |>
      Seurat::FindNeighbors(verbose = F) |>
      Seurat::FindClusters(verbose = F, algorithm = 1, resolution = resolution_SoupX) |>
      SeuratObject::FetchData(vars = paste0("RNA_snn_res.", resolution_SoupX))

    sc <- SoupX::setClusters(sc = sc,
                             clusters = stats::setNames(as.character(clusters[,1,drop = T]),
                                                        rownames(clusters)))

    message("Running SoupX.")

    sc <- Gmisc::fastDoCall(what = SoupX::autoEstCont,
                            args = c(list(sc = sc), SoupX_autoEstCont_args))
    sx_counts <- SoupX::adjustCounts(sc = sc, verbose = 0)

    if (SoupX_return) {
      message("Creating Seurat object on SoupX-corrected count matrix with ", ncol(filt_data), " cells.")
      SO_sx <- Seurat::CreateSeuratObject(counts = sx_counts)
      SO_sx@meta.data$orig.ident <- x

      SO_sx <- SO_prep02(SO_unprocessed = stats::setNames(list(SO_sx), x),
                         reductions = "umap",
                         nhvf = nhvf,
                         npcs = npcs,
                         batch_corr = batch_corr,
                         RunHarmony_args = list(group.by.vars = "orig.ident"),
                         FindClusters_args = list(resolution = resolution),
                         normalization = "LogNormalize",
                         diet_seurat = F)

      SO_sx <- SeuratObject::AddMetaData(object = SO_sx,
                                         metadata = (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100,
                                         col.name = "pct_soup_SoupX")

      sc <- SoupX::setDR(sc = sc,
                         DR = SO_sx@reductions$umap@cell.embeddings)

      sc_info_df <- data.frame(n_expr_uncorrected = Matrix::rowSums(sc$toc > 0),
                               n_expr_corrected = Matrix::rowSums(sx_counts > 0)) |>
        dplyr::mutate(abs_diff = n_expr_uncorrected-n_expr_corrected) |>
        dplyr::mutate(rel_diff = abs_diff/n_expr_uncorrected) |>
        dplyr::filter(abs_diff > 0) |>
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

  if (SoupX_return) {
    SOx <- SO_prep02(SO_unprocessed = purrr::map(SoupX_results, ~purrr::pluck(.x, "SO")),
                     reductions = "umap",
                     nhvf = nhvf,
                     npcs = npcs,
                     min_cells = 1,
                     batch_corr = batch_corr,
                     RunHarmony_args = list(group.by.vars = "orig.ident"),
                     FindClusters_args = list(resolution = resolution),
                     normalization = "LogNormalize",
                     diet_seurat = F)
    # rm seurat to save ram
    SoupX_results <- purrr::map(SoupX_results, ~.x[-1])
    SoupX_results <- list(SoupX = SOx, SoupX_results = SoupX_results)
  }

  return(SoupX_results)
}

run_decontx <- function(SO, resolution, nhvf) {
  message("Running decontX.")
  SO <- purrr::map(SO, function(SO) {
    ## multi dirs: split matrix
    split_mats <- brathering::split_mat(x = SeuratObject::LayerData(SO, assay = "RNA", layer = "counts"),
                                        f = SO@meta.data$orig.ident,
                                        byrow = F)
    split_idents <- split(x = SO@meta.data[[paste0("RNA_snn_res.", resolution)]],
                          f = SO@meta.data$orig.ident)

    dx <- purrr::map2(.x = split_mats,
                      .y = split_idents,
                      ~celda::decontX(x = .x,
                                      z = .y,
                                      varGenes = nhvf)[["contamination"]],
                      nhvf = nhvf)
    dx <- unname(unlist(dx))
    names(dx) <- unname(unlist(purrr::map(split_mats, colnames)))
    SO <- Seurat::AddMetaData(object = SO,
                              metadata = dx,
                              col.name = "pct_soup_decontX")
    return(SO)
  })
  return(SO)
}

cluster_on_metadata <- function(SO,
                                batch_corr,
                                PCs_to_meta_clustering,
                                ffbms,
                                resolution_meta) {
  SO <- purrr::map(SO, function(SO) {

    # differentiate mouse, human or no MT-genes at all
    # and add freq of RPS / RPL and MRPS / MRPL genes
    #grep("^MT-|mt-", c("MT-iii", "mt-zzz"), value = T)
    # regex for or: |
    qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt")
    if (any(grepl("^MT-", rownames(SO))) && !any(grepl("^mt-", rownames(SO)))) {
      SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^MT-"), "pct_mt")
      SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^RP[SL]"), "pct_ribo")
      SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^MRP[SL]"), "pct_mribo")
    } else if (!any(grepl("^MT-", rownames(SO))) && any(grepl("^mt-", rownames(SO)))) {
      SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^mt-"), "pct_mt")
      SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^Rp[sl]-"), "pct_ribo")
      SO <- Seurat::AddMetaData(SO, Seurat::PercentageFeatureSet(SO, pattern = "^Mrp[sl]-"), "pct_mribo")
    } else {
      message("No mitochondrial genes could be identified from gene names - none starting with MT- (human) or mt- (mouse).")
      qc_cols <- qc_cols[-which(qc_cols == "pct_mt")]
    }
    if ("dbl_score" %in% names(SO@meta.data)) {
      qc_cols <- c(qc_cols, "dbl_score")
    }
    qc_cols <- paste0(qc_cols, "_log")

    if ("pct_soup_decontX" %in% names(SO@meta.data)) {
      qc_cols <- c(qc_cols, "pct_soup_decontX")
    }
    if ("pct_soup_SoupX" %in% names(SO@meta.data)) {
      qc_cols <- c(qc_cols, "pct_soup_SoupX")
    }

    SO@meta.data$nFeature_RNA_log <- log1p(SO@meta.data$nFeature_RNA)
    SO@meta.data$nCount_RNA_log <- log1p(SO@meta.data$nCount_RNA)
    if ("pct_mt" %in% names(SO@meta.data)) {
      SO@meta.data$pct_mt_log <- log1p(SO@meta.data$pct_mt)
      # this could be done above with SO
      # replace NA wiht 0 to avoid error in umap calculation
      SO@meta.data$pct_mt <- ifelse(is.na(SO@meta.data$pct_mt), 0, SO@meta.data$pct_mt)
      SO@meta.data$pct_mt_log <- ifelse(is.na(SO@meta.data$pct_mt_log), 0, SO@meta.data$pct_mt_log)
    }

    ## multi-dirs: split matrix!
    SO@meta.data$residuals <- unlist(lapply(unique(SO@meta.data$orig.ident),
                                            function(x) stats::residuals(stats::lm(nCount_RNA_log~nFeature_RNA_log,
                                                                                   data = SO@meta.data[which(SO@meta.data$orig.ident == x),]))))

    ## clustering on meta data (quality metrics)
    message("Running dimension reduction and clustering on qc meta data.")

    for (nn in PCs_to_meta_clustering) {
      meta2 <- dplyr::select(SO@meta.data, dplyr::all_of(qc_cols), residuals)
      if (nn > 0) {
        meta2 <- cbind(meta2, SO@reductions[[ifelse(batch_corr == "harmony" && length(ffbms) > 1, "harmony", "pca")]]@cell.embeddings[,1:nn]) # length(ffbms) or length(data_dirs)
      }

      meta2 <- brathering::scale2(meta2)
      #https://datascience.stackexchange.com/questions/27726/when-to-use-cosine-simlarity-over-euclidean-similarity
      # with cosine metric: relative composition is more important
      umap_dims <- suppressWarnings(uwot::umap(X = meta2, metric = "cosine"))
      colnames(umap_dims) <- paste0("meta_PC", nn, "_UMAP_", c(1,2))

      reductioname <- paste0("umapmetaPC", nn)
      SO[[reductioname]] <- Seurat::CreateDimReducObject(
        embeddings = umap_dims,
        key = paste0(toupper(reductioname), "_"),
        assay = "RNA"
      )

      # fix colnames manually as Seurat::CreateDimReducObject makes a mistake in taking the trailing number of umap_dims-colnames as trailing dim-number
      colnames(SO@reductions[[reductioname]]@cell.embeddings) <- paste0(toupper(paste0(reductioname, "_")), c(1,2))

      # matrix has to be supplied to FindNeighbors
      clusters <- Seurat::FindClusters(Seurat::FindNeighbors(meta2, annoy.metric = "cosine", verbose = F)$snn,
                                       resolution = resolution_meta,
                                       verbose = F)
      nclust <- apply(clusters, 2, function(x) length(unique(x)))
      candidates <- names(nclust)[which(dplyr::between(nclust, 1,12))]
      choice <- ifelse(!length(candidates), names(nlust)[1], candidates[length(candidates)])
      clusters <- clusters[,choice,drop = F]

      # add reduction belonging to the meta clustering as name
      metaclustname <- stats::setNames(paste0("meta_PC", nn, "_", colnames(clusters)), reductioname)
      colnames(clusters) <- metaclustname
      SO <- Seurat::AddMetaData(SO, cbind(umap_dims, clusters))

      suppressWarnings(SeuratObject::Misc(SO, "meta_clusterings") <- c(SeuratObject::Misc(SO, "meta_clusterings"), metaclustname))
    }

    return(SO)
  })

}
