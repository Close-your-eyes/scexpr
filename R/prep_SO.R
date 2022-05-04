#' Prepare a somewhat standardized Seurat Object for further analyses
#'
#'
#'
#' @param SO_unprocessed named list of split Seurat objects
#' @param samples names of SO_unprocessed that are to include; if missing all are used
#' @param cells vector of cell names to include; if missint all cells are used
#' @param min_cells min number of cells per object in SO_unprocessed; objects with lower cell number are excluded
#' @param downsample downsample a fraction of cells from each object in SO_unprocessed (downsample < 0);
#' or an absolute number of cells per object (downsample > 1 & downsample >= min_cells)
#' @param export_prefix prefix to the output rds file
#' @param cluster_resolutions vector of resolutions to compute by Seurat::FindClusters
#' @param reductions which reduction to calculate, tsne --> Seurat::RunTSNE, umap --> Seurat::RunUMAP
#' @param nhvf number high variables features, passed to Seurat::SCTransform
#' or Seurat::FindVariableFeatures
#' @param npcs number of principle components in calculate in pca and to
#' consider for downstream functions like tSNE or UMAP
#' @param nintdims number of integration dimension, when batch_corr = "integration"
#' @param normalization algorithm for normalization of UMIs
#' @param batch_corr algorithm for batch correction between samples in SO_unprocessed;
#' only relevant if more than 1 sample is passed
#' @param ref_sample optional names of the reference sample in SO_unprocessed for
#' Seurat::FindIntegrationAnchors or RunHarmony; reference; if NULL no reference is used
#' @param integr_reduction reduction method for Seurat::FindIntegrationAnchors; rpca is quicker than cca
#' @param vars.to.regress passed to Seurat::SCTransform or Seurat::ScaleData, applied if batch_corr == "regression"
#' @param seeed seed passed to several methods
#' @param save_path folder to save the resulting Seurat object as rds-file to
#' @param ... additional arguments passed to harmony::RunHarmony, Seurat::RunUMAP,
#' Seurat::RunTSNE, Seurat::FindNeighbors, Seurat::FindClusters,
#' Seurat::PrepSCTIntegration, Seurat::FindIntegrationAnchors, Seurat::IntegrateData;
#' prefix arguments for RunHarmony by "RunHarmony__"; so e.g. RunHarmony__theta, prefix arguments
#' for EmbedSOM::SOM with "SOM__" (e.g. SOM__batch = T or SOM__rlen = 20,), for EmbedSOM::GQTSOM with "GQTSOM__" (e.g. GQTSOM__distf = 4)
#' and for EmbedSOM::EmbedSOM with "EmbedSOM__"
#' @param celltype_refs list(prim_cell_atlas = celldex::HumanPrimaryCellAtlasData())
#' @param celltype_label label name to use from the reference data set
#' @param celltype_ref_clusters use a clustering as basis for grouped celltype annotation with SingleR. This will
#' fundamentally speed up the calculation!
#' e.g. if SCT assay is used and cluster_resolutions includes 0.8 then pass: SCT_snn_res.0.8.
#' when integration procedure assay is used: integrated_snn_res.0.8. For RNA assay: RNA_snn_res.0.8.
#' @param diet_seurat logical whether to run Seurat::DietSeurat
#' @param verbose print messages and progress bars from functions
#'
#' @return Seurat Object, as R object and saved to disk as rds file
#' @export
#'
#' @examples
#' \dontrun{
#' }
prep_SO <- function(SO_unprocessed,
                    samples = NULL,
                    cells = NULL,
                    min_cells = 50,
                    downsample = 1,
                    export_prefix = NULL,
                    cluster_resolutions = 0.8,
                    reductions = c("umap", "som", "gqtsom", "tsne"),
                    nhvf = 800,
                    npcs = 20,
                    nintdims = 30,
                    normalization = c("SCT", "LogNormalize"),
                    batch_corr = c("harmony", "integration", "regression", "none"),
                    ref_sample = NULL,
                    integr_reduction = c("rpca", "cca", "rlsi"),
                    vars.to.regress = NULL,
                    seeed = 42,
                    save_path = NULL,
                    celltype_refs = NULL, # list of celldex::objects
                    celltype_label = "label.main",
                    celltype_ref_clusters = NULL,
                    diet_seurat = T,
                    verbose = F,
                    ...) {


  mydots <- list(...)

  options(warn = 1)
  if (is.null(save_path)) {
    message("No save_path to save Seurat objects to provided.")
  } else if (!is.character(save_path) || length(save_path) != 1) {
    stop("save_path has to be a character; a path to a folder where to save Seurat objects to.")
  }

  if (!is.null(celltype_refs)) {
    if (!is.list(celltype_refs)) {
      stop("celltype_refs should be a named list of reference data sets from celldex or similar.")
    }
    if (is.null(names(celltype_refs))) {
      stop("Provide names for celltype_refs.")
    }
    if (any(duplicated(names(celltype_refs)))) {
      stop("names for celltype_refs must be unique.")
    }
    if (length(celltype_label) == 1) {
      celltype_label <- rep(celltype_label, length(celltype_refs))
    }
    if (length(celltype_label) != length(celltype_refs)) {
      stop("celltype_label and celltype_refs must have the same lengths.")
    }
    for (i in seq_along(celltype_refs)) {
      if (!celltype_label[i] %in% names(celltype_refs[[i]]@colData@listData)) {
        stop(paste0(celltype_label[i], " not found in ", names(celltype_refs)[i], "."))
      }
    }

    if (!requireNamespace("SingleR", quietly = T)) {
      BiocManager::install("SingleR")
    }
  }

  reductions <- match.arg(tolower(reductions), c("tsne", "umap", "som", "gqtsom"), several.ok = T)
  normalization <- match.arg(normalization, c("SCT", "LogNormalize"))
  batch_corr <- match.arg(batch_corr, c("harmony", "integration", "regression", "none"))
  integr_reduction <- match.arg(integr_reduction, c("rpca", "cca", "rlsi"))

  if (!is.null(celltype_ref_clusters)) {
    #check if celltype_ref_clusters can exist
    pref <- ifelse(batch_corr == "integration", "integrated", ifelse(normalization == "SCT", "SCT", "RNA"))
    mid <- "_snn_res."
    suf <- gsub("\\.0$", "", cluster_resolutions)
    candidates <- paste0(pref,mid,suf)
    if (!celltype_ref_clusters %in% candidates) {
      stop("celltype_ref_clusters not found in candidates based on batch_corr, normalization and cluster_resolutions: ", paste(candidates, collapse = ", "))
    }
  }

  if (class(SO_unprocessed) == "list") {
    SO.list <- SO_unprocessed
  } else if (class(SO_unprocessed) == "character") {
    if (!file.exists(SO_unprocessed)) {
      stop(paste0(SO_unprocessed, "not found."))
    } else {
      if (!grepl("\\.rds$", SO_unprocessed, ignore.case = T)) {
        stop("SO_unprocessed has to be an .rds file.")
      }
      SO.list <- readRDS(SO_unprocessed)
    }
  } else {
    stop("SO_unprocessed has to be named list of splitted Seurat objects or a path (character) to an .rds file of those. If it is only one Seurat object (one sample)
         make it a list of length 1.")
  }

  if (class(SO.list) == "Seurat") {
    # if only one Seurat object is provided
    SO.list <- list(SO.list)
    names(SO.list) <- "sample"
  }

  if (is.null(names(SO.list))) {
    stop("SO_unprocessed has no names.")
  }

  if (downsample > 1 && downsample < min_cells) {
    message("downsample set to min_cells.")
    downsample <- min_cells
  }

  if (is.null(samples)) {
    samples <- names(SO.list)
  } else {
    if (any(!samples %in% names(SO.list))) {
      message("samples not found in SO_unprocessed: ", samples[which(!samples %in% names(SO.list))])
    }

  }
  samples <- names(SO.list)[which(grepl(paste(samples, collapse = "|"), names(SO.list)))]
  SO.list <- SO.list[which(names(SO.list) %in% samples)]


  if (length(SO.list) > 1 && batch_corr == "harmony" && !any(grepl("RunHarmony__group.by.vars", names(mydots)))) {
    stop("Please provide one or more group.by.vars from meta.data (with prefix: RunHarmony__group.by.vars) for RunHarmony: ", paste(names(SO.list[[1]]@meta.data), collapse = ", "), ".")
  }

  if (!is.null(cells)) {
    SO.list <- lapply(SO.list, function (x) {
      inds <- which(Seurat::Cells(x) %in% cells)
      if (length(inds) == 0) {
        return(NULL)
      }
      return(subset(x, cells = Seurat::Cells(x)[which(Seurat::Cells(x) %in% cells)]))
    })
  }
  if (any(sapply(SO.list, is.null))) {
    message("No cells found for: ", paste(names(SO.list)[which(sapply(SO.list, is.null))], collapse = ", "))
    SO.list <- SO.list[which(!sapply(SO.list, is.null))]
    if (length(SO.list) == 0) {
      stop("No Seurat objects left after filtering for cells.")
    }
  }

  if (downsample < 1) {
    SO.list <- lapply(SO.list, function (x) subset(x, cells = sample(Seurat::Cells(x), as.integer(downsample*length(Seurat::Cells(x))), replace = FALSE)))
  } else if (downsample > 1) {
    SO.list <- lapply(SO.list, function (x) subset(x, cells = sample(Seurat::Cells(x), downsample, replace = FALSE)))
  }

  # remove samples with insufficient number of cells
  rm.nm <- names(SO.list[which(sapply(SO.list, function(x){length(Seurat::Cells(x)) < min_cells}))])
  if (length(rm.nm) > 0) {
    SO.list <- SO.list[which(!names(SO.list) %in% rm.nm)]
    message("samples removed due to min.cells: ", paste(rm.nm, collapse = ","))
  }

  message("Samples included (", length(SO.list), "): ", paste(names(SO.list), collapse=", "))

  if (length(SO.list) == 1) {
    batch_corr <- "none"
  }

  # cases:
  if (length(SO.list) == 1) {
    SO <- SO.list[[1]]
    if (normalization == "SCT") {
      SO <- suppressWarnings(Seurat::SCTransform(SO, variable.features.n = nhvf, vars.to.regress = vars.to.regress, verbose = verbose, assay = "RNA"))
      SO <- Seurat::NormalizeData(SO, verbose = verbose, assay = "RNA")
      SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
    } else if (normalization == "LogNormalize") {
      SO <- Seurat::NormalizeData(SO, verbose = verbose, assay = "RNA")
      SO <- Seurat::FindVariableFeatures(SO, selection.method = "vst", nfeatures = nhvf, verbose = verbose, assay = "RNA")
      SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
    }
    SO <- suppressWarnings(Seurat::ProjectDim(Seurat::RunPCA(object = SO, npcs = npcs, verbose = verbose, seed.use = seeed), reduction = "pca", do.center = T, overwrite = F, verbose = verbose))
  }

  if (length(SO.list) > 1) {
    ## 3 cases
    if (batch_corr %in% c("regression", "none", "harmony")) {
      SO <- merge(x = SO.list[[1]], y = SO.list[2:length(SO.list)], merge.data = T)
      if (batch_corr %in% c("none", "harmony") && !is.null(vars.to.regress)) {
        vars.to.regress <- NULL
        message("vars.to.regress set to NULL as batch_corr %in% c('none', 'harmony').")
      }
      if (normalization == "SCT") {
        SO <- suppressWarnings(Seurat::SCTransform(SO, variable.features.n = nhvf, vars.to.regress = vars.to.regress, seed.use = seeed, verbose = verbose))
        SO <- Seurat::NormalizeData(SO, verbose = verbose, assay = "RNA")
        SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
      } else if (normalization == "LogNormalize") {
        SO <- Seurat::ScaleData(Seurat::FindVariableFeatures(Seurat::NormalizeData(SO, verbose = verbose), selection.method = "vst", nfeatures = nhvf, verbose = verbose), vars.to.regress = vars.to.regress, verbose = verbose)
      }
      SO <- suppressWarnings(Seurat::ProjectDim(Seurat::RunPCA(object = SO, npcs = npcs, verbose = verbose, seed.use = seeed), reduction = "pca", do.center = T, overwrite = F, verbose = verbose))
    }

    if (batch_corr == "harmony") {
      if (!requireNamespace("harmony", quietly = T)) {
        utils::install.packages("harmony")
      }
      dots <- mydots[which(grepl("^RunHarmony__", names(mydots), ignore.case = T))]
      names(dots) <- gsub("^RunHarmony__", "", names(dots), ignore.case = T)
      SO <- suppressWarnings(do.call(harmony::RunHarmony, args = c(list(object = SO, assay.use = switch(normalization, SCT = "SCT", LogNormalize = "RNA"), reference_values = ref_sample, verbose = verbose), dots)))
      SO <- Seurat::NormalizeData(SO, verbose = verbose, assay = "RNA")
      SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
    }

    if (batch_corr == "integration") {
      k.filter <- as.integer(min(200, min(sapply(SO.list, ncol))/2))
      k.score <- as.integer(min(30, min(sapply(SO.list, ncol))/6))

      if (normalization == "SCT") {
        SO.list <- lapply(SO.list, FUN = Seurat::SCTransform, verbose = verbose, ...)
        anchor_features <- Seurat::SelectIntegrationFeatures(SO.list, verbose = verbose, ...)
        SO.list <- Seurat::PrepSCTIntegration(object.list = SO.list, anchor.features = anchor_features, verbose = verbose, ...) # error by wrong argument
      } else {
        anchor_features <- Seurat::SelectIntegrationFeatures(SO.list, verbose = verbose, ...)
      }

      if (integr_reduction == "rpca") {
        if (normalization == "SCT") {
          SO.list <- lapply(SO.list, function(x) {
            x <- Seurat::RunPCA(x, features = anchor_features, verbose = verbose)
            return(x)
          })
        } else {
          SO.list <- lapply(SO.list, function(x) {
            x <- Seurat::ScaleData(x, features = anchor_features, verbose = verbose)
            x <- Seurat::RunPCA(x, features = anchor_features, verbose = verbose)
            return(x)
          })
        }
      }

      anchorset <- Seurat::FindIntegrationAnchors(object.list = SO.list, anchor.features = anchor_features, dims = 1:nintdims, normalization.method = normalization, reference = ref_sample, k.filter = k.filter, k.score = k.score, reduction = integr_reduction, ...)
      SO <- Seurat::IntegrateData(anchorset = anchorset, dims = 1:nintdims, normalization.method = normalization, features.to.integrate = rownames(Seurat::GetAssayData(SO.list[[1]], assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA"))), k.weight = k.filter, ...)
      Seurat::DefaultAssay(SO) <- "integrated"

      if (normalization == "SCT") {
        SO <- suppressWarnings(Seurat::ProjectDim(Seurat::RunPCA(object = SO, npcs = npcs, verbose = verbose, seed.use = seeed), reduction = "pca", do.center = T, overwrite = F, verbose = verbose))
        SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = verbose)
      }
      if (normalization == "LogNormalize") {
        SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = verbose)
        SO <- Seurat::ScaleData(SO, verbose = verbose)
        SO <- suppressWarnings(Seurat::ProjectDim(Seurat::RunPCA(object = SO, npcs = npcs, verbose = verbose, seed.use = seeed), reduction = "pca", do.center = T, overwrite = F, verbose = verbose))
      }
    }
  }

  #dots <- mydots[which(grepl("^PrepSCTIntegration__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^PrepSCTIntegration__", "", names(dots), ignore.case = T)
  #SO.list <- do.call(Seurat::PrepSCTIntegration, args = c(list(object.list = SO.list, verbose = verbose), dots))

  #dots <- mydots[which(grepl("^FindIntegrationAnchors__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^FindIntegrationAnchors__", "", names(dots), ignore.case = T)
  #anchorset <- do.call(Seurat::FindIntegrationAnchors, args = c(list(object.list = SO.list, dims = 1:nintdims, normalization.method = normalization, reference = ref_sample, k.filter = k.filter, k.score = k.score, reduction = integr_reduction), dots))

  #dots <- mydots[which(grepl("^IntegrateData__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^IntegrateData__", "", names(dots), ignore.case = T)
  #SO <- do.call(Seurat::IntegrateData, args = c(list(anchorset = anchorset, dims = 1:nintdims, normalization.method = normalization, features.to.integrate = rownames(Seurat::GetAssayData(SO.list[[1]], assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA"))), k.weight = k.filter), dots))


  ### do.call on large SeuratObject became super slow, not practicable!

  red <- switch(batch_corr, harmony = "harmony", integration = "pca", regression = "pca", none = "pca")
  if (any(grepl("umap", reductions, ignore.case = T))) {
    SO <- Seurat::RunUMAP(object = SO, umap.method = "uwot", metric = "cosine", dims = 1:npcs, seed.use = seeed, reduction = red, verbose = verbose, ...)
  }

  #dots <- mydots[which(grepl("^RunUMAP__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^RunUMAP__", "", names(dots), ignore.case = T)
  #SO <- do.call(Seurat::RunUMAP, args = c(list(object = SO, umap.method = "uwot", dims = 1:npcs, seed.use = seeed, reduction = red, verbose = T), dots))

  #dots <- mydots[which(grepl("^RunTSNE__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^RunTSNE__", "", names(dots), ignore.case = T)
  #SO <- do.call(Seurat::RunTSNE, args = c(list(object = SO, dims = 1:npcs, seed.use = seeed, reduction = red, verbose = T, num_threads = 0), dots))

  if (any(grepl("tsne", reductions, ignore.case = T))) {
    SO <- Seurat::RunTSNE(object = SO, dims = 1:npcs, seed.use = seeed, reduction = red, verbose = verbose, num_threads = 0, ...)
  }

  if (any(grepl("^som$", reductions, ignore.case = T))) {
    if (!requireNamespace("devtools", quietly = T)) {
      utils::install.packages("devtools")
    }
    if (!requireNamespace("EmbedSOM", quietly = T)) {
      devtools::install_github("exaexa/EmbedSOM")
    }

    dots <- mydots[which(grepl("^SOM__", names(mydots), ignore.case = T))]
    names(dots) <- gsub("^SOM__", "", names(dots), ignore.case = T)
    map <- do.call(EmbedSOM::SOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings), dots))

    dots <- mydots[which(grepl("^EmbedSOM__", names(mydots), ignore.case = T))]
    names(dots) <- gsub("^EmbedSOM__", "", names(dots), ignore.case = T)
    ES <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings, map = map), dots))

    SO[["SOM"]] <- Seurat::CreateDimReducObject(embeddings = ES, key = "SOM_", assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA"), misc = mydots[which(grepl("^SOM__|^EmbedSOM__", names(mydots), ignore.case = T))])
  }

  if (any(grepl("^gqtsom$", reductions, ignore.case = T))) {
    if (!requireNamespace("devtools", quietly = T)) {
      utils::install.packages("devtools")
    }
    if (!requireNamespace("EmbedSOM", quietly = T)) {
      devtools::install_github("exaexa/EmbedSOM")
    }

    dots <- mydots[which(grepl("^GQTSOM__", names(mydots), ignore.case = T))]
    names(dots) <- gsub("^GQTSOM__", "", names(dots), ignore.case = T)
    map <- do.call(EmbedSOM::GQTSOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings), dots))

    dots <- mydots[which(grepl("^EmbedSOM__", names(mydots), ignore.case = T))]
    names(dots) <- gsub("^EmbedSOM__", "", names(dots), ignore.case = T)
    ES <- do.call(EmbedSOM::EmbedSOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings, map = map), dots))

    SO[["GQTSOM"]] <- Seurat::CreateDimReducObject(embeddings = ES, key = "GQTSOM_", assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA"), misc = mydots[which(grepl("^GQTSOM__|^EmbedSOM__", names(mydots), ignore.case = T))])
  }

  #dots <- mydots[which(grepl("^FindNeighbors__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^FindNeighbors__", "", names(dots), ignore.case = T)
  #SO <- do.call(Seurat::FindNeighbors, args = c(list(object = SO, reduction = red, dims = 1:npcs), dots))

  #dots <- mydots[which(grepl("^FindClusters__", names(mydots), ignore.case = T))]
  #names(dots) <- gsub("^FindClusters__", "", names(dots), ignore.case = T)
  #SO <- do.call(Seurat::FindClusters, args = c(list(object = SO, resolution = cluster_resolutions), dots))


  SO <- Seurat::FindNeighbors(object = SO, reduction = red, dims = 1:npcs, verbose = verbose, ...)
  SO <- Seurat::FindClusters(object = SO, resolution = cluster_resolutions, verbose = verbose, ...)

  if (!is.null(celltype_refs)) {
    if (!celltype_ref_clusters %in% names(SO@meta.data) && !is.null(celltype_ref_clusters)) {
      message("celltype_ref_clusters not found in SO meta data. Will be set to NULL and SingleR will operate on single cell level.")
      celltype_ref_clusters <- NULL
    }
    for (i in seq_along(celltype_refs)) {
      celltypes <- SingleR::SingleR(test = Seurat::GetAssayData(SO, slot = "data", assay = "RNA"),
                                    ref = celltype_refs[[i]],
                                    labels = celltype_refs[[i]]@colData@listData[[celltype_label[i]]],
                                    clusters = SO@meta.data[,celltype_ref_clusters])
      if (is.null(celltype_ref_clusters)) {
        SO@meta.data[,paste0(names(celltype_refs)[i], "_labels")] <- celltypes$labels
      } else {
        celltypes_df <- utils::stack(stats::setNames(celltypes$labels, levels(SO@meta.data[,celltype_ref_clusters])))
        names(celltypes_df) <- c(paste0(names(celltype_refs)[i], "_labels"), celltype_ref_clusters)
        celltypes_df <- tibble::column_to_rownames(dplyr::left_join(tibble::rownames_to_column(SO@meta.data[,celltype_ref_clusters,drop=F], "ID"), celltypes_df, by = celltype_ref_clusters), "ID")
        SO <- Seurat::AddMetaData(SO, celltypes_df)
      }
      Seurat::Misc(SO, paste0(names(celltype_refs)[i], "_object")) <- celltypes
    }
  }

  # remove counts as they can be recalculated with rev_lognorm
  if (diet_seurat) {
    Seurat::Misc(SO, slot = "RNA_count_colSums") <- unname(Matrix::colSums(Seurat::GetAssayData(SO, slot = "counts", assay = "RNA")))
    SO <- Seurat::DietSeurat(SO, assays = names(SO@assays), counts = F, dimreducs = names(SO@reductions))
  }

  if (!is.null(save_path)) {
    dir.create(save_path, showWarnings = F, recursive = T)
    save.time <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%y%m%d-%H%M%S")
    saveRDS(SO, compress = T, file = file.path(save_path, paste("SO", export_prefix, normalization, batch_corr, downsample, nhvf, npcs, paste0(save.time, ".rds"), sep = "_")))
    message(paste0("SO saved to: ", save_path))
  }

  return(SO)
}

