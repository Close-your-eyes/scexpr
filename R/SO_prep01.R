#' Prepare and quality-control 10x single-cell RNA-seq data
#'
#' @description
#' Reads one or more 10x Genomics gene-expression count matrices, creates
#' per-sample Seurat objects, and runs a standard preprocessing and
#' quality-control workflow. Samples are merged before expression-based
#' dimensional reduction and clustering. When more than one sample is present,
#' Harmony can be used for batch correction.
#'
#' @details
#' The function supports two input modes. By default, `data_dirs` is searched
#' recursively for directories named `filtered_feature_bc_matrix` or
#' `filtered_gene_bc_matrices` and, when requested for SoupX, matching raw
#' matrix directories. Alternatively, named filtered-matrix directories can be
#' supplied through `ffbms`. Direct input currently disables SoupX and DecontX.
#'
#' Quality-control metadata include log-transformed UMI and detected-feature
#' counts, species-dependent mitochondrial and ribosomal feature percentages,
#' and, when requested and successfully calculated, scDblFinder and ambient-RNA
#' estimates. These variables are rescaled and used for a separate UMAP and
#' graph clustering. Setting `PCs_to_meta_clustering` above zero adds the
#' indicated number of leading expression PCs to this QC representation.
#'
#' SoupX is run only when a raw matrix is found for every discovered sample.
#' If any raw matrix is missing, SoupX is disabled for the complete run. DecontX
#' and scDblFinder are run separately by sample. SoupX- and DecontX-derived
#' values are estimates and are not automatic cell-removal thresholds.
#'
#' With `early_exit = TRUE`, the function stops after reading the matrices,
#' applying cell and feature filters, creating per-sample Seurat objects, and
#' adding the available basic QC or doublet metadata. Expression preprocessing,
#' ambient-RNA estimation, merged clustering, and QC-metadata clustering are
#' skipped.
#'
#' @param data_dirs Character vector or list of directories to search
#'   recursively for 10x filtered-matrix folders. This argument is ignored when
#'   `ffbms` or `rfbms` is supplied.
#' @param nhvf Positive integer giving the number of highly variable features
#'   used during expression preprocessing, doublet detection, SoupX processing,
#'   and DecontX estimation.
#' @param npcs Positive integer giving the number of expression principal
#'   components to calculate and use in downstream analyses.
#' @param resolution Numeric vector of graph-clustering resolutions for the
#'   expression-based analysis. One resulting clustering is selected for
#'   downstream QC calculations.
#' @param resolution_SoupX Single numeric graph-clustering resolution used to
#'   define clusters for SoupX contamination estimation.
#' @param resolution_meta Numeric vector of graph-clustering resolutions for
#'   clustering the rescaled QC metadata.
#' @param PCs_to_meta_clustering Numeric vector. For each value `n`, a separate
#'   QC representation is created using the QC variables and the first `n`
#'   expression PCs. Use `0` for QC variables alone. Positive values should be
#'   whole numbers no greater than `npcs`.
#' @param scDblFinder Logical; run scDblFinder independently on each sample and
#'   add `dbl_score` and `dbl_class` to cell metadata. A failed sample-level run
#'   is reported and processing continues without those fields for that sample.
#' @param min_UMI Numeric scalar or `NULL`. Cells with fewer total counts than
#'   this value are removed while reading each filtered matrix. Use `NULL` to
#'   disable this filter.
#' @param min_UMI_var_feat Numeric scalar or `NULL`. When `scDblFinder = TRUE`,
#'   cells with fewer counts across the selected variable features are removed
#'   before doublet detection. Use `NULL` to disable this filter.
#' @param SoupX Logical; estimate ambient-RNA contamination with SoupX. Requires
#'   matching filtered and raw matrix directories for every sample discovered
#'   through `data_dirs`.
#' @param decontX Logical; estimate ambient-RNA contamination with
#'   `celda::decontX()` separately for each sample and add the result as
#'   `pct_soup_decontX`. Direct `ffbms`/`rfbms` input currently disables this
#'   option.
#' @param SoupX_return Logical; if SoupX runs, also preprocess and return a
#'   Seurat object made from the SoupX-adjusted counts together with per-sample
#'   SoupX diagnostics. Otherwise, only the estimated contamination percentage
#'   is added to the original object.
#' @param SoupX_autoEstCont_args Named list of additional arguments passed to
#'   `SoupX::autoEstCont()`. Do not supply `sc` or `tfidfMin`, which are set by
#'   this function.
#' @param cells Optional character vector of unprefixed 10x cell barcodes to
#'   retain in each sample. Selection occurs before sample prefixes are added.
#' @param invert_cells Logical; if `TRUE`, exclude the matching barcodes in
#'   `cells` instead of retaining them.
#' @param feature_rm Optional character vector of feature names to remove after
#'   any feature aggregation.
#' @param feature_aggr Optional named list of character vectors. Each element
#'   identifies features whose counts are summed into a new feature named after
#'   that list element. Aggregation precedes `feature_rm`.
#' @param ffbms Optional named character vector of filtered 10x matrix
#'   directories. Each directory must contain a matrix HDF5 file or the standard
#'   compressed Matrix Market files. Supplying this argument bypasses discovery
#'   through `data_dirs`.
#' @param rfbms Optional named character vector of raw 10x matrix directories,
#'   paired to `ffbms` by name. In the current direct-input branch, SoupX and
#'   DecontX are disabled, so these paths are not used for contamination
#'   estimation.
#' @param diet_seurat Logical; after the complete workflow, reduce returned
#'   Seurat objects with `Seurat::DietSeurat()` and remove the RNA `data` and
#'   `scale.data` layers. Ignored when `early_exit = TRUE`.
#' @param batch_corr Character scalar, either `"harmony"` or `"none"`. Harmony
#'   is used only for multiple samples; for a single sample this is changed to
#'   `"none"`.
#' @param reduction Character scalar, either `"tsne"` or `"umap"`, specifying
#'   the expression-based nonlinear reduction requested from `SO_prep02()`.
#' @param sample_prefix_to_cell_id Logical; prepend each sample name and `"__"`
#'   to cell barcodes to keep cell names unique across samples.
#' @param early_exit Logical; return per-sample objects before merged expression
#'   processing and contamination or QC-metadata analyses.
#' @param equalize_feature_order Logical; with `early_exit = TRUE`, restrict and
#'   reorder features consistently across samples using
#'   `make_equal_feature_order()`.
#' @param common_cells Logical; with `early_exit = TRUE`, restrict samples to
#'   barcodes shared across all inputs after removing sample prefixes.
#' @param mc.cores Positive integer giving the number of forked workers supplied
#'   to `parallel::mclapply()` when reading and creating sample objects.
#'
#' @return A named list. With `early_exit = TRUE`, each element is a per-sample
#'   Seurat object. Otherwise, the list contains `original`, the processed
#'   Seurat object based on the original counts. When both `SoupX = TRUE` and
#'   `SoupX_return = TRUE` remain enabled after input validation, it additionally
#'   contains `SoupX`, a processed Seurat object based on adjusted counts, and
#'   `SoupX_results`, a named list of per-sample SoupX diagnostic objects and
#'   summaries.
#'
#' @export
#' @importFrom zeallot %<-%
#'
#' @examples
#' \dontrun{
#' # Search sample directories for 10x matrix folders.
#' samples <- list.dirs("path/to/cellranger", recursive = FALSE)
#' out <- SO_prep01(samples)
#'
#' # Read named filtered-matrix directories and return before integration.
#' filtered <- c(sample_A = "sample_A/filtered_feature_bc_matrix",
#'               sample_B = "sample_B/filtered_feature_bc_matrix")
#' out_early <- SO_prep01(ffbms = filtered, early_exit = TRUE)
#' }
SO_prep01 <- function(data_dirs,
                      nhvf = 2000,
                      npcs = 10,
                      resolution = seq(0.1,0.8,0.1),
                      resolution_SoupX = 0.6,
                      resolution_meta = seq(0.1,0.8,0.1),
                      PCs_to_meta_clustering = c(0,2),
                      scDblFinder = F,
                      min_UMI = 200,
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
                      diet_seurat = F,
                      batch_corr = c("harmony", "none"),
                      reduction = c("umap", "tsne"),
                      sample_prefix_to_cell_id = T,
                      early_exit = F,
                      equalize_feature_order = F,
                      common_cells = F,
                      mc.cores = 4) {


  install_pkgs(SoupX, scDblFinder, decontX)
  resolution <- resolution_checks(resolution_meta, resolution_SoupX, resolution, PCs_to_meta_clustering)
  batch_corr <- rlang::arg_match(batch_corr)
  reduction <- rlang::arg_match(reduction)

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
  check_inputs2(sample_folders = ffbms)


  message("Reading filtered_feature_bc_matrix data.")
  SO <- parallel::mclapply(purrr::set_names(names(ffbms)), function(x) {
    read_10X_data(path = ffbms[x],
                  cells = cells,
                  min_UMI = min_UMI,
                  name = x,
                  sample_prefix_to_cell_id = sample_prefix_to_cell_id,
                  invert_cells = invert_cells)
  }, mc.cores = mc.cores)

  SO <- purrr::discard(SO, is.null)

  if (!length(SO)) {
    stop("No filtered_feature_bc_matrix successfully read.")
  }

  if (early_exit) {
    ## check it here, otherwise it is done in SO_prep02
    if (equalize_feature_order) {
      SO <- scexpr:::make_equal_feature_order(SO)
    }
    if (common_cells) {
      SO <- make_equal_cells(SO)
    }
  }


  SO <- parallel::mclapply(SO, function(x) {
    if (!is.null(feature_rm) || !is.null(feature_aggr)) {
      x <- aggregate_or_remove_features(filt_data = x,
                                        feature_rm = feature_rm,
                                        feature_aggr = feature_aggr)
    }

    message("Creating initial Seurat object with ", ncol(x), " cells.")
    x <- Seurat::CreateSeuratObject(counts = x)

    if (early_exit) {
      x <- add_pct_featset_and_cc(x)
    }

    # doublet score calculation with not yet merged data
    if (scDblFinder) {
      tryCatch(expr = {
        x <- add_dbl_score_to_metadata(SO = x,
                                       nhvf = nhvf,
                                       min_UMI_var_feat = min_UMI_var_feat,
                                       npcs = npcs)
      }, error = function(err){
        print(err)
      })

    }
    return(x)
  }, mc.cores = mc.cores)

  if (scDblFinder && any(purrr::map_lgl(SO, ~!"dbl_score" %in% names(.x@meta.data)))) {
    message("doublet calculation failed in at least one object.")
  }

  for (i in names(SO)) {
    SO[[i]]@meta.data$orig.ident <- i
  }

  if (early_exit) {
    return(SO)
  }

  if (length(SO) > 1) {
    message("Preparing merged and harmonized Seurat object with ", sum(lengths(purrr::map(SO, Seurat::Cells))), " cells.")
  } else {
    message("Preparing Seurat object.")
  }

  SO <- SO_prep02(SO_unprocessed = SO,
                  reductions = reduction,
                  nhvf = nhvf,
                  npcs = npcs,
                  min_cells = 1,
                  batch_corr = batch_corr,
                  RunHarmony_args = list(group.by.vars = "orig.ident"),
                  FindClusters_args = list(resolution = resolution),
                  normalization = "LogNormalize",
                  interactive_varfeat_selection = F,
                  interactive_pc_selection = F)

  # only save a useful cluster resolution
  allclust <- SO@misc$clusterings[-length(SO@misc$clusterings)]
  nclust <- apply(SO@meta.data[,allclust], 2, function(x) length(unique(x)))
  candidates <- names(nclust)[which(dplyr::between(nclust, 1,12))]
  choice <- ifelse(!length(candidates), names(nclust)[1], candidates[length(candidates)])
  # for qc_plot2; in analogy to meta_clustering
  names(choice) <- reduction
  resolution <- as.numeric(brathering::strsplit2(choice, "_", -1)[[1]][2])
  suppressWarnings(SeuratObject::Misc(SO, "clusterings") <- choice)

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
                               feature_rm = feature_rm,
                               reduction = reduction)

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
    SO <- purrr::map(SO, ~Seurat::DietSeurat(.x, dimreducs = names(.x@reductions)))

    SO <- purrr::map(SO, function(x) {
      x@assays[["RNA"]]@layers[["data"]] <- NULL
      x@assays[["RNA"]]@layers[["scale.data"]] <- NULL
      return(x)
    })
  }

  message("Run scexpr::qc_plot2() on the Seurat object for visualization of clustering on phenotype and qc data.")
  if (SoupX && SoupX_return) {
    message("Use e.g. SoupX::plotChangeMap(x[['sc']], cleanedMatrix = SoupX::adjustCounts(x[['sc']]), geneSet = 'GNLY')")
    return(c(SO, list(SoupX_results = SoupX_results)))
  } else {
    return(SO)
  }
}

# alldir <- list.dirs(datadir, recursive = TRUE, full.names = TRUE)
# seldir <- alldir[basename(alldir) %in% c("filtered_feature_bc_matrix", "filtered_gene_bc_matrices", "raw_feature_bc_matrix", "raw_gene_bc_matrices")]
# allfile <- list.files(datadir, recursive = TRUE, full.names = TRUE)
# selfile <- allfile[basename(allfile) %in% c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")]
#
# find_divergence_folder <- function(paths) {
#   # Split paths into components
#   split_paths <- strsplit(normalizePath(paths), .Platform$file.sep)
#
#   # Convert to matrix (pad with NA if needed)
#   max_len <- max(lengths(split_paths))
#   mat <- t(sapply(split_paths, function(x) {
#     c(x, rep(NA, max_len - length(x)))
#   }))
#
#   # Find last column where all values are identical
#   same_cols <- apply(mat, 2, function(col) length(unique(col[!is.na(col)])) == 1)
#   last_common_idx <- max(which(same_cols))
#
#   list(
#     lowest_common_folder = file.path(mat[1, seq_len(last_common_idx)]),
#     divergence_level = mat[, last_common_idx + 1]
#   )
# }
#
# prefer_h5 <- function(x) {
#   h5 <- x[grepl("\\.h5$", x)]
#   if (length(h5) > 0) {
#     unique(h5)
#   } else {
#     unique(x)
#   }
# }
# prefer_h5_per_type <- function(x) {
#   split_types <- split(
#     x,
#     ifelse(grepl("filtered", x), "filtered", "raw")
#   )
#
#   out <- lapply(split_types, function(paths) {
#     h5 <- paths[grepl("\\.h5$", paths)]
#     if (length(h5) > 0) unique(h5) else unique(paths)
#   })
#   unlist(out)
# }
#
#
# divfolder1 <- find_divergence_folder(seldir)$divergence_level
# divfolder2 <- find_divergence_folder(selfile)$divergence_level
# names(seldir) <- divfolder1
# names(selfile) <- divfolder2
# out <- split(c(seldir, selfile), names(c(seldir, selfile)))
# out <- lapply(out, prefer_h5_per_type)
#

check_dir <- function(data_dirs, SoupX = F) {

  dir_roots <- unlist(lapply(data_dirs, function(x) {
    dd <- list.dirs(x)
    dirname(dd[which(grepl("filtered_feature_bc_matrix|filtered_gene_bc_matrices", basename(dd)))])
  }))

  if (is.null(names(dir_roots))) {
    names(dir_roots) <- basename(dir_roots)
    names(dir_roots)[which(names(dir_roots) == "outs")] <- basename(dirname(dir_roots))[which(names(dir_roots) == "outs")]
    names(dir_roots)[which(names(dir_roots) == "count")] <- basename(dirname(dir_roots))[which(names(dir_roots) == "count")]
    names(dir_roots)[which(names(dir_roots) == "SoupX")] <- basename(dirname(dir_roots))[which(names(dir_roots) == "SoupX")]
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
  if (scDblFinder && !requireNamespace("scDblFinder", quietly = T)) {
    pak::pak("plger/scDblFinder")
  }
  if (decontX && !requireNamespace("celda", quietly = T)) {
    BiocManager::install("celda")
  }
  if (!requireNamespace("scuttle", quietly = T)) {
    BiocManager::install("scuttle")
  }
  if (!requireNamespace("presto", quietly = T)) {
    pak::pak("immunogenomics/presto")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }
}


resolution_checks <- function(resolution_meta, resolution_SoupX, resolution, PCs_to_meta_clustering) {
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

  factors <- scuttle::librarySizeFactors(get_layer(obj = SO, assay = "RNA", layer = "counts", features = var_feat))
  zero_libsize_cells <- names(which(factors == 0))
  if (length(zero_libsize_cells) > 0) {
    message(length(zero_libsize_cells), " cell(s) found which have zero library size based on hvf. These are removed to allow running scDblFinder. See https://github.com/plger/scDblFinder/issues/55.")
    SO <- subset2(SO, cells = setdiff(names(factors), zero_libsize_cells))
  }

  if (!is.null(min_UMI_var_feat)) {
    UMI_sum <- Matrix::colSums(get_layer(obj = SO, assay = "RNA", layer = "counts", features = var_feat))
    if (any(UMI_sum < min_UMI_var_feat)) {
      message(sum(UMI_sum < min_UMI_var_feat), " cells removed for having less UMI in variable features as min_UMI_var_feat. See https://github.com/LTLA/BiocNeighbors/issues/24.")
    }
    SO <- subset2(SO, cells = names(UMI_sum[which(UMI_sum >= min_UMI_var_feat)]))
  }
  # scDblFinder

  # SO@meta.data$dbl_score <- scDblFinder::computeDoubletDensity(x = get_layer(obj = SO, assay = "RNA", layer = "counts"),
  #                                                              subset.row = var_feat,
  #                                                              dims = npcs)
  #   SO@meta.data$dbl_score_log <- log1p(SO@meta.data$dbl_score)
  scf <- scDblFinder::scDblFinder(sce = get_layer(obj = SO,
                                                  assay = "RNA",
                                                  layer = "counts"),
                                  nfeatures = var_feat,
                                  dims = npcs)
  SO@meta.data$dbl_score <- scf@colData$scDblFinder.score
  SO@meta.data$dbl_class <- scf@colData$scDblFinder.class


  return(SO)
}

read_10X_data <- function(path,
                          name = NULL,
                          cells = NULL,
                          min_UMI = 1,
                          verbose = T,
                          sample_prefix_to_cell_id = T,
                          invert_cells = F) {


  h5files <- list.files(path, pattern = "\\.h5$", full.names = T)
  if (length(h5files)) {
    if (length(h5files) > 1 && verbose) {
      message("Found more than one .h5 file in ", path, ". Will use the first: ", h5files[1])
    }
    filt_data <- Seurat::Read10X_h5(filename = h5files[1])
  } else {
    filt_data <- Seurat::Read10X(data.dir = path)
  }
  #filt_data <- Seurat::ReadSTARsolo(data.dir = path)
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
      return(NULL)
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
  if (sample_prefix_to_cell_id) {
    if (is.null(name)) {
      # assumption of min prefix length of 3; maybe wrong sometimes
      prefixes <- substr(colnames(filt_data), 1, 3)
      if (length(unique(prefixes)) != 1) {
        message("random prefix added to cell names.")
        name <- stringr::str_pad(sample(1:1000, 1), width = 4, pad = "0")
        colnames(filt_data) <- paste0(name, "__", colnames(filt_data))
      }
    } else {
      if (!all(grepl(paste0("^", name), colnames(filt_data)))) {
        colnames(filt_data) <- paste0(name, "__", colnames(filt_data))
      }
    }

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

check_inputs <- function(data_dirs,
                         ffbms = NULL,
                         rfbms = NULL,
                         SoupX = F,
                         SoupX_return = F,
                         decontX = F,
                         batch_corr = "none") {

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

  return(list(
    ffbms = ffbms,
    rfbms = rfbms,
    SoupX = SoupX,
    SoupX_return = SoupX_return,
    decontX = decontX,
    batch_corr = batch_corr
  ))
}

run_soupx <- function(ffbms,
                      rfbms,
                      SO = NULL,
                      nhvf = 2000,
                      min_UMI = 1,
                      resolution_SoupX = 0.6,
                      npcs = 10,
                      batch_corr = "harmony",
                      resolution = 0.8,
                      SoupX_return = F,
                      SoupX_autoEstCont_args = list(),
                      feature_aggr = NULL,
                      feature_rm = NULL,
                      reduction = "umap") {
  # ffmbs and rfbms are paired by name

  # use filt_data which may have been reduced 'cells' selection; raw_feature_bc_matrix will provide the whole picture of the soup
  message("SoupX: Reading filtered and raw_feature_bc_matrix data.")

  # if (!is.list(ffbms)) {
  #   ffbms <- list(sample = ffbms)
  # }
  # if (!is.list(rfbms)) {
  #   rfbms <- list(sample = rfbms)
  # }

  SoupX_results <- purrr::map(purrr::set_names(names(rfbms)), function(x) {
    message(x)
    # just read again
    cells <- NULL
    if (!is.null(SO)) {
      cells <- Seurat::Cells(SO)
    }
    filt_data <- read_10X_data(path = ffbms[x],
                               cells = cells, ## filter for existing cells in SO (potential scDblFinder filtering, above)
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
    if (is.null(SO)) {
      clusters <- Seurat::CreateSeuratObject(filt_data)
    } else {
      clusters <- subset2(Seurat::DietSeurat(SO, layers = "counts"), cells = intersect(Seurat::Cells(SO), rownames(sc$metaData)))
    }
    clusters <- clusters |>
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

    sc <- tryCatch(expr = {
      Gmisc::fastDoCall(what = SoupX::autoEstCont,
                        args = c(list(sc = sc, tfidfMin = 1), SoupX_autoEstCont_args))
    },
    error = function(err) {
      Gmisc::fastDoCall(what = SoupX::autoEstCont,
                        args = c(list(sc = sc, tfidfMin = 0.5), SoupX_autoEstCont_args))
    })

    sx_counts <- SoupX::adjustCounts(sc = sc, verbose = 0)

    if (SoupX_return) {
      message("Creating Seurat object on SoupX-corrected count matrix with ", ncol(filt_data), " cells.")
      SO_sx <- Seurat::CreateSeuratObject(counts = sx_counts)
      SO_sx@meta.data$orig.ident <- x

      SO_sx <- SO_prep02(SO_unprocessed = stats::setNames(list(SO_sx), x),
                         reductions = reduction,
                         nhvf = nhvf,
                         npcs = npcs,
                         batch_corr = batch_corr,
                         RunHarmony_args = list(group.by.vars = "orig.ident"),
                         FindClusters_args = list(resolution = resolution),
                         normalization = "LogNormalize",
                         interactive_varfeat_selection = F,
                         interactive_pc_selection = F)

      SO_sx <- SeuratObject::AddMetaData(object = SO_sx,
                                         metadata = (Matrix::colSums(sc[["toc"]]) - Matrix::colSums(sx_counts))/Matrix::colSums(sc[["toc"]])*100,
                                         col.name = "pct_soup_SoupX")

      sc <- SoupX::setDR(sc = sc,
                         DR = SO_sx@reductions[[length(SO_sx@reductions)]]@cell.embeddings)

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
                     reductions = reduction,
                     nhvf = nhvf,
                     npcs = npcs,
                     min_cells = 1,
                     batch_corr = batch_corr,
                     RunHarmony_args = list(group.by.vars = "orig.ident"),
                     FindClusters_args = list(resolution = resolution),
                     normalization = "LogNormalize",
                     interactive_varfeat_selection = F,
                     interactive_pc_selection = F)
    # rm seurat to save ram
    SoupX_results <- purrr::map(SoupX_results, ~.x[-1])
    SoupX_results <- list(SoupX = SOx, SoupX_results = SoupX_results)
  }

  return(SoupX_results)
}



run_decontx <- function(SO, resolution, nhvf) {
  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }
  message("Running decontX.")
  SO <- purrr::map(SO, function(SO) {
    ## multi dirs: split matrix
    split_mats <- brathering::split_mat(x = get_layer(obj = SO, assay = "RNA", layer = "counts"),
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
    #tt <- grep("^MT-", rownames(SO), value=T)

    qc_cols <- c("nCount_RNA", "nFeature_RNA", "pct_mt")
    SO <- add_pct_featset_and_cc(SO)
    if (!"pct_mt" %in% names(SO@meta.data)) {
      qc_cols <- qc_cols[-which(qc_cols == "pct_mt")]
    }
    qc_cols <- paste0(qc_cols, "_log")

    if ("dbl_score" %in% names(SO@meta.data)) {
      qc_cols <- c(qc_cols, "dbl_score")
    }
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
      # replace NA with 0 to avoid error in umap calculation
      SO@meta.data$pct_mt <- ifelse(is.na(SO@meta.data$pct_mt), 0, SO@meta.data$pct_mt)
      SO@meta.data$pct_mt_log <- ifelse(is.na(SO@meta.data$pct_mt_log), 0, SO@meta.data$pct_mt_log)
    }

    ## multi-dirs: split matrix!
    # SO@meta.data$residuals <- unlist(lapply(unique(SO@meta.data$orig.ident),
    #                                         function(x) stats::residuals(stats::lm(nCount_RNA_log~nFeature_RNA_log,
    #                                                                                data = SO@meta.data[which(SO@meta.data$orig.ident == x),]))))

    ## clustering on meta data (quality metrics)
    message("Running dimension reduction and clustering on qc meta data.")

    for (nn in PCs_to_meta_clustering) {
      meta2 <- dplyr::select(SO@meta.data, dplyr::all_of(qc_cols))
      if (nn > 0) {
        meta2 <- cbind(meta2, SO@reductions[[ifelse(batch_corr == "harmony" && length(ffbms) > 1,
                                                    grep("^harmony", names(SO@reductions), value = T),
                                                    grep("^pca", names(SO@reductions), value = T))]]@cell.embeddings[,1:nn]) # length(ffbms) or length(data_dirs)
      }

      meta2 <- apply(meta2, 2, function(x) scales::rescale(x, to = c(0,1)))
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

      clusters <- scexpr:::pad_default_cluster_numbers(x = clusters)
      clusters <- choose_top_ncluster_in_range(x = clusters)


      # add reduction belonging to the meta clustering as name
      metaclustname <- stats::setNames(paste0("meta_PC", nn, "_", colnames(clusters)), reductioname)
      colnames(clusters) <- metaclustname
      SO <- Seurat::AddMetaData(SO, cbind(umap_dims, clusters))

      suppressWarnings(SeuratObject::Misc(SO, "metacolors") <- c(SeuratObject::Misc(SO, "metacolors"),
                                                                 stats::setNames(list(colrr::col_pal("custom", n = sort(unique(clusters[,1])), return = "c")), nm = unname(metaclustname))))
      suppressWarnings(SeuratObject::Misc(SO, "meta_clusterings") <- c(SeuratObject::Misc(SO, "meta_clusterings"), metaclustname))
    }

    return(SO)
  })

}

choose_top_ncluster_in_range <- function(x, range = c(1,12)) {

  nclust <- apply(x, 2, function(y) length(unique(y)))
  nclust <- sort(nclust)
  candidates <- names(nclust)[which(dplyr::between(nclust, range[1], range[2]))]
  choice <- ifelse(!length(candidates), names(nclust)[1], candidates[length(candidates)])
  x <- x[,choice,drop = F]
  return(x)

}


make_equal_cells <- function(x) {
  # x: list of vectors
  cells <- purrr::map(x, colnames)

  # remove prefix
  cells_rmprefix <- purrr::map(purrr::set_names(names(cells)), function(x) {
    cells[[x]] <- gsub(x, "", cells[[x]])
    return(cells[[x]])
  })

  common_cells <- purrr::reduce(cells_rmprefix, intersect)
  if (any(length(common_cells) < lengths(cells_rmprefix))) {
    message("different cells across input samples. reducing to common: n = ", length(common_cells))
    match_inds <- purrr::map(cells_rmprefix, ~which(.x %in% common_cells))
    x <- purrr::map2(x, match_inds, function(x,y) x[,y])
  }

  return(x)
}



#' Add typical feature set pct and cell cycle score
#'
#' @param obj seurat object
#' @param species human or mouse or ..auto.. for guessing
#'
#' @returns obj with updated meta data
#' @export
#'
#' @examples
#' \dontrun{
#' so <- add_pct_featset_and_cc(so)
#' }
add_pct_featset_and_cc <- function(obj, species = "..auto..") {

  if (!requireNamespace("UCell", quietly = T)) {
    BiocManager::install("UCell")
  }

  if (species == "..auto..") {
    species <- guess_species(get_gene_features(obj))
  } else {
    species <- rlang::arg_match(species, values = c("human", "mouse"))
  }

  # RPS* (cytosolic small subunit)
  # RPL* (cytosolic large subunit)
  # MRPS* (mitochondrial small subunit)
  # MRPL* (mitochondrial large subunit)
  if (species == "human") {
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^MT-"), "pct_mt")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^HSP"), "pct_hsp")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^(FOS|FOSB|JUN|JUNB|ATF3|EGR1|DUSP1|HSPA1A|HSPA1B)$"), "pct_iestress")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^HB[ABDEGMQZ]"), "pct_hb")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^RP[SL][0-9]+[A-Z]?$"), "pct_ribo")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^MRP[SL]"), "pct_mribo")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^MT[1-2][A-Z]$"), "pct_metallothionein")

    cc_genes <- get_cell_cycle_genesets()[["cc_lst"]]
    cc_genes_seu <- cc_genes[c("seurat_2019_g2m", "seurat_2019_s")]
    names(cc_genes_seu) <- c("G2M_score", "S_score")
    obj <- UCell::AddModuleScore_UCell(obj, features = cc_genes_seu, ncores = 4, name = "")
  } else if (species == "mouse") {
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^mt-"), "pct_mt")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^Rp[sl]"), "pct_ribo")
    obj <- Seurat::AddMetaData(obj, Seurat::PercentageFeatureSet(obj, pattern = "^Mrp[sl]"), "pct_mribo")
  }

  return(obj)
}

guess_species <- function(genes) {

  human_like <- mean(grepl("^(MT-|RPL|RPS|HLA-)", genes))
  mouse_like <- mean(grepl("^(mt-|Rpl|Rps|H2-)", genes))

  species <- if (human_like > mouse_like) "human" else "mouse"

  return(species)
}

check_inputs2 <- function(sample_folders) {
  files_missing <- purrr::map_lgl(sample_folders, function(x) {
    files <- list.files(x)
    if (any(tools::file_ext(files) == "h5")) {
      return(FALSE)
    }
    return(any(!"barcodes.tsv.gz" %in% files | !"matrix.mtx.gz" %in% files | !"features.tsv.gz" %in% files))
  })

  if (any(files_missing)) {
    files_missing2 <- paste(names(files_missing)[which(files_missing)], collapse = ", ")
    stop("these folders do not have all required files: ", files_missing2, "\nrequired are barcodes.tsv.gz, matrix.mtx.gz, features.tsv.gz.")
  }
}
