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
#' @param reductions which reduction to calculate, one of multiple of
#' c("umap", "som", "gqtsom", "tsne")
#' @param nhvf number high variables features, passed to Seurat::SCTransform
#' or Seurat::FindVariableFeatures or Seurat::SelectIntegrationFeatures
#' @param npcs number of principle components in calculate in pca and to
#' consider for downstream functions like tSNE or UMAP
#' @param normalization algorithm for normalization of UMIs
#' @param batch_corr procedure for batch correction between samples in SO_unprocessed;
#' only relevant if more than 1 sample is passed
#' @param vars.to.regress passed to Seurat::SCTransform or Seurat::ScaleData; will be applied
#' independent of what is set for batch_corr; if batch_corr is set to 'none', then only vars.to.regress is
#' used to regress out a variable in meta.data which may be sample-specific; other then that
#' it may not be meaningful to regress a sample-specific variable and perform batch_corr;
#' rather vars.to.regress may be used in combination with batch_corr to regress percent_mito or so
#' @param seeed seed passed to several methods
#' @param save_path folder to save the resulting Seurat object as rds-file to
#' @param celltype_refs list(prim_cell_atlas = celldex::HumanPrimaryCellAtlasData(), MonacoImmune = celldex::MonacoImmuneData())
#' @param celltype_label list of label names to use from the reference data set (one list entry per celltype_refs; each entry may contain a vector of labels)
#' @param celltype_ref_clusters use a clustering as basis for grouped celltype annotation with SingleR. This will
#' fundamentally speed up the calculation!
#' e.g. if SCT assay is used and cluster_resolutions includes 0.8 then pass: SCT_snn_res.0.8.
#' when integration procedure assay is used: integrated_snn_res.0.8. For RNA assay: RNA_snn_res.0.8.
#' @param diet_seurat logical whether to run Seurat::DietSeurat
#' @param verbose print messages and progress bars from functions
#' @param var_feature_filter character vector of features which are to exclude from variable features;
#' this will affect downstream PCA, dimension reduction and clustering;
#' e.g.: var_feature_filter = grep("^TR[ABGD]V", rownames(SOqc_split[[1]]), value = T) to
#' exclude T cell receptor gene segments
#' @param hvf_determination_before_merge determine variable features before or after merging multiple samples;
#' either by the standard LogNormalize workflow or by SCTransform function; this is most important for
#' deciding whether SCtransform is run on mulitple samples separately before merging (set to TRUE) or
#' after merging (set to FALSE); in my experience and when batch_corr = harmony, setting
#' it to FALSE yields better results; see: https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06a_integration_harmony.md and subsequent links
#' @param FindVariableFeatures_args
#' @param SCtransform_args
#' @param RunUMAP_args
#' @param RunTSNE_args
#' @param FindNeighbors_args
#' @param FindClusters_args
#' @param RunHarmony_args
#' @param SOM_args
#' @param GQTSOM_args
#' @param EmbedSOM_args
#' @param FindIntegrationAnchors_args
#' @param IntegrateData_args by default features.to.integrate = rownames(Seurat::GetAssayData(SO_unprocessed[[1]], assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA")))
#' @param RunPCA_args
#' @param ...
#'
#' @return Seurat Object, as R object and saved to disk as rds file
#' @export
#' @importFrom zeallot %<-%
#'
#' @examples
#' \dontrun{
#' }
SO_prep02 <- function(SO_unprocessed,
                    samples = NULL,
                    cells = NULL,
                    min_cells = 50,
                    downsample = 1,
                    export_prefix = NULL,
                    reductions = c("umap"),
                    nhvf = 800,
                    npcs = 20,
                    normalization = c("SCT", "LogNormalize"),
                    hvf_determination_before_merge = F,
                    batch_corr = c("harmony", "integration", "none"),
                    vars.to.regress = NULL,
                    seeed = 42,
                    save_path = NULL,
                    celltype_refs = NULL, # list of celldex::objects
                    celltype_label = "label.main",
                    celltype_ref_clusters = NULL,
                    diet_seurat = F,
                    var_feature_filter = NULL,
                    verbose = F,
                    FindVariableFeatures_args = list(),
                    SCtransform_args = list(vst.flavor = "v2", method = "glmGamPoi"),
                    RunPCA_args = list(),
                    RunUMAP_args = list(),
                    RunTSNE_args = list(theta = 0),
                    FindNeighbors_args = list(),
                    FindClusters_args = list(resolution = seq(0.1,1,0.1)),
                    RunHarmony_args = list(group.by.vars = "orig.ident"),
                    SOM_args = list(),
                    GQTSOM_args = list(),
                    EmbedSOM_args = list(),
                    FindIntegrationAnchors_args = list(reduction = "rpca"),
                    IntegrateData_args = list(),
                    scale_RNA_assay_when_SCT = T,
                    ...) {

  mydots <- list(...)
  options(warn = 1)

  reductions <- match.arg(tolower(reductions), c("umap", "som", "gqtsom", "tsne"), several.ok = T)
  normalization <- rlang::arg_match(normalization)
  batch_corr <- rlang::arg_match(batch_corr)

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
  if (!requireNamespace("Gmisc", quietly = TRUE)) {
    utils::install.packages("Gmisc")
  }
  if (batch_corr == "harmony") {
    if (!requireNamespace("harmony", quietly = T)) {
      utils::install.packages("harmony")
    }
  }
  if (!requireNamespace("zap", quietly = T)) {
    devtools::install_github("coolbutuseless/zap")
  }

  if (is.null(save_path)) {
    message("No save_path to save Seurat objects to provided.")
  } else if (!is.character(save_path) || length(save_path) != 1) {
    stop("save_path has to be a character; a path to a folder where to save Seurat objects to.")
  }

  celltype_label <- check_celltype_refs(celltype_refs, celltype_label)
  check_RunPCA_args(RunPCA_args)
  check_celltype_ref_clusters(celltype_ref_clusters,
                              batch_corr,
                              normalization,
                              FindClusters_args)
  c(SO_unprocessed, samples) %<-% check_SO_unprocessed_and_samples(SO_unprocessed,
                                                                   samples,
                                                                   batch_corr,
                                                                   RunHarmony_args)
  SO_unprocessed <- subset_SO_unprocessed(SO_unprocessed, cells, downsample, min_cells)
  if (length(SO_unprocessed) == 1) {batch_corr <- "none"}

  SCtransform_args <- check_SCtransform_args(SCtransform_args,
                                             nhvf,
                                             vars.to.regress)
  FindVariableFeatures_args <- check_FindVariableFeatures_args(FindVariableFeatures_args,
                                                               nhvf)



  # cases:
  if (length(SO_unprocessed) == 1) {
    SO <- make_so_single(SO_unprocessed,
                         normalization,
                         SCtransform_args,
                         seeed,
                         verbose,
                         var_feature_filter,
                         nhvf,
                         vars.to.regress,
                         FindVariableFeatures_args,
                         scale_RNA_assay_when_SCT,
                         RunPCA_args,
                         npcs)
  }


  if (length(SO_unprocessed) > 1) {
    if (batch_corr %in% c("none", "harmony")) {
      SO <- make_so_multi_harmony(SO_unprocessed,
                                  normalization,
                                  SCtransform_args,
                                  seeed,
                                  verbose,
                                  var_feature_filter,
                                  nhvf,
                                  vars.to.regress,
                                  FindVariableFeatures_args,
                                  scale_RNA_assay_when_SCT,
                                  RunPCA_args,
                                  npcs,
                                  RunHarmony_args,
                                  batch_corr,
                                  hvf_determination_before_merge)
    } else if (batch_corr == "integration") {
      SO <- make_so_multi_integrate(SO_unprocessed,
                                    normalization,
                                    SCtransform_args,
                                    seeed,
                                    verbose,
                                    var_feature_filter,
                                    nhvf,
                                    vars.to.regress,
                                    FindVariableFeatures_args,
                                    IntegrateData_args,
                                    hvf_determination_before_merge,
                                    scale_RNA_assay_when_SCT,
                                    npcs)
    }
  }

  ### do.call on large SeuratObject became super slow, not practicable!
  # https://stackoverflow.com/questions/28198103/alternative-to-do-call-for-large-datasets

  # alternatives:
  # https://rlang.r-lib.org/reference/exec.html
  # gmisc::fastDoCall

  red <- switch(batch_corr, harmony = "harmony", integration = "pca", none = "pca")
  if (any(grepl("umap", reductions, ignore.case = T))) {
    RunUMAP_args <- RunUMAP_args[which(!names(RunUMAP_args) %in% c("object", "seed.use", "reduction", "verbose"))]
    if (!"dims" %in% names(RunUMAP_args)) {
      RunUMAP_args <- c(list(dims = 1:npcs), RunUMAP_args)
    }
    SO <- Gmisc::fastDoCall(Seurat::RunUMAP, args = c(list(object = SO,
                                                           reduction = red,
                                                           seed.use = seeed,
                                                           verbose = verbose),
                                                      RunUMAP_args))
    #SO <- Seurat::RunUMAP(object = SO, umap.method = "uwot", metric = "cosine", dims = 1:npcs, seed.use = seeed, reduction = red, verbose = verbose, ...)
  }

  if (any(grepl("tsne", reductions, ignore.case = T))) {
    RunTSNE_args <- RunTSNE_args[which(!names(RunTSNE_args) %in% c("object", "seed.use", "reduction", "verbose"))]
    if (!"dims" %in% names(RunTSNE_args)) {
      RunTSNE_args <- c(list(dims = 1:npcs), RunTSNE_args)
    }
    if (!"num_threads" %in% names(RunTSNE_args)) {
      RunTSNE_args <- c(list(num_threads = 0), RunTSNE_args)
    }
    SO <- Gmisc::fastDoCall(Seurat::RunTSNE, args = c(list(object = SO,
                                                           reduction = red,
                                                           seed.use = seeed,
                                                           verbose = verbose),
                                                      RunTSNE_args))
    #SO <- Seurat::RunTSNE(object = SO, dims = 1:npcs, seed.use = seeed, reduction = red, verbose = verbose, num_threads = 0, ...)
  }

  if (any(grepl("^som$", reductions, ignore.case = T))) {
    if (!requireNamespace("devtools", quietly = T)) {
      utils::install.packages("devtools")
    }
    if (!requireNamespace("EmbedSOM", quietly = T)) {
      devtools::install_github("exaexa/EmbedSOM")
    }

    SOM_args <- SOM_args[which(!names(SOM_args) %in% c("data"))]
    map <- Gmisc::fastDoCall(EmbedSOM::SOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings)), SOM_args)
    EmbedSOM_args <- EmbedSOM_args[which(!names(EmbedSOM_args) %in% c("data", "map"))]
    ES <- Gmisc::fastDoCall(EmbedSOM::EmbedSOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings, map = map)), EmbedSOM_args)
    SO[["SOM"]] <- Seurat::CreateDimReducObject(embeddings = ES, key = "SOM_", assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA"), misc = list(SOM_args, EmbedSOM_args))
  }

  if (any(grepl("^gqtsom$", reductions, ignore.case = T))) {
    if (!requireNamespace("devtools", quietly = T)) {
      utils::install.packages("devtools")
    }
    if (!requireNamespace("EmbedSOM", quietly = T)) {
      devtools::install_github("exaexa/EmbedSOM")
    }

    GQTSOM_args <- GQTSOM_args[which(!names(GQTSOM_args) %in% c("data"))]
    map <- Gmisc::fastDoCall(EmbedSOM::GQTSOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings)), GQTSOM_args)
    EmbedSOM_args <- EmbedSOM_args[which(!names(EmbedSOM_args) %in% c("data", "map"))]
    ES <- Gmisc::fastDoCall(EmbedSOM::EmbedSOM, args = c(list(data = SO@reductions[[red]]@cell.embeddings, map = map)), EmbedSOM_args)
    SO[["GQTSOM"]] <- Seurat::CreateDimReducObject(embeddings = ES, key = "GQTSOM_", assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA"), misc = list(GQTSOM_args, EmbedSOM_args))
  }

  FindNeighbors_args <- FindNeighbors_args[which(!names(FindNeighbors_args) %in% c("object", "reduction", "verbose"))]
  if (!"dims" %in% names(FindNeighbors_args)) {
    FindNeighbors_args <- c(list(dims = 1:npcs), FindNeighbors_args)
  }
  SO <- Gmisc::fastDoCall(Seurat::FindNeighbors, args = c(list(object = SO,
                                                               reduction = red,
                                                               verbose = verbose),
                                                          FindNeighbors_args))
  FindClusters_args <- FindClusters_args[which(!names(FindClusters_args) %in% c("object", "verbose"))]
  SO <- Gmisc::fastDoCall(Seurat::FindClusters, args = c(list(object = SO,
                                                              verbose = verbose),
                                                         FindClusters_args))

  SO <- run_celltyping(SO,
                       celltype_refs,
                       celltype_ref_clusters,
                       celltype_label)


  # remove counts as they can be recalculated with rev_lognorm
  if (diet_seurat) {
    Seurat::Misc(SO, slot = "RNA_count_colSums") <- unname(Matrix::colSums(Gmisc::fastDoCall(what = Seurat::GetAssayData, args = GetAssayData_args)))
    SO <- Seurat::DietSeurat(SO, assays = names(SO@assays), counts = F, dimreducs = names(SO@reductions))
  }

  save.time <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%y%m%d_%H%M%S")
  ext <- ifelse(!requireNamespace("zap", quietly = T), "rds", "zap")
  save.name <- paste("SO", export_prefix, normalization, batch_corr, downsample, nhvf, npcs, paste0(save.time, ".", ext), sep = "_")
  Seurat::Misc(SO, "object_name") <- save.name

  if (!is.null(save_path)) {
    dir.create(save_path, showWarnings = F, recursive = T)
    if (ext == "rds") {
      saveRDS(SO, compress = F, file = file.path(save_path, save.name))
    } else if (ext == "zap") {
      zap::zap_write(SO, dst = file.path(save_path, save.name), compress = "zstd")
    }
    message(paste0("SO saved to: ", save_path, " as ", save.name, "."))
  }

  return(SO)
}

.var_feature_filter_removal <- function(SO,
                                        var_feature_filter,
                                        max_rounds = 5,
                                        normalization,
                                        nhvf,
                                        vars.to.regress,
                                        seeed,
                                        verbose,
                                        SCtransform_args,
                                        FindVariableFeatures_args) {

  ## rerun SCTransform (or FindVariableFeatures) until nhvf is met while filtering for var_feature_filter
  ## max iterations is 5; maybe print some messages about progress
  ## filtering for var_feature_filter will affect PCA and all other downstream calculation which are based on PCA

  ## intention: filter TCR/BCR chain genes
  ## alternatively one could also regress them out?!

  if (any(!var_feature_filter %in% rownames(SO))) {
    ## print message?!
    var_feature_filter <- var_feature_filter[which(var_feature_filter %in% rownames(SO))]
    if (length(var_feature_filter) == 0) {
      return(SO)
    }
  }

  n <- 1
  while (any(Seurat::VariableFeatures(SO) %in% var_feature_filter) &&
         length(Seurat::VariableFeatures(SO)[which(!Seurat::VariableFeatures(SO) %in% var_feature_filter)]) < nhvf &&
         n <= max_rounds) {
    message(sum(Seurat::VariableFeatures(SO) %in% var_feature_filter), " of var_feature_filter found in Variable Features. Round ", n, " of increasing nhvf to ", nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), " so that var_feature_filter can be removed while nhvf = " , nhvf, " is met.")

    if (normalization == "SCT") {
      SCtransform_args[which(names(SCtransform_args) == "variable.features.n")] <- SCtransform_args[["variable.features.n"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter))
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 seed.use = seeed,
                                                                 verbose = verbose),
                                                            SCtransform_args))

      #SO <- Seurat::SCTransform(SO, assay = "RNA", vst.flavor = "v2", method = "glmGamPoi", variable.features.n = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), vars.to.regress = vars.to.regress, seed.use = seeed, verbose = verbose)
    }
    if (normalization == "LogNormalize") {
      FindVariableFeatures_args[which(names(FindVariableFeatures_args) == "nfeatures")] <- FindVariableFeatures_args[["nfeatures"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter))
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = "RNA",
                                                                          verbose = verbose),
                                                                     FindVariableFeatures_args))
      #SO <- Seurat::FindVariableFeatures(SO, selection.method = "vst", nfeatures = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), verbose = verbose, assay = "RNA")
    }
    n <- n + 1
  }
  Seurat::VariableFeatures(SO) <- Seurat::VariableFeatures(SO)[which(!Seurat::VariableFeatures(SO) %in% var_feature_filter)]

  return(SO)
}


check_RunPCA_args <- function(RunPCA_args) {
  # actually only very few arguments are allowed to be passed by RunPCA_args. Otherwise the function would break.
  if ("npcs" %in% names(RunPCA_args)) {
    stop("Please do not pass npcs in RunPCA_args. Use the npcs argument.")
  }
  if ("seed.use" %in% names(RunPCA_args)) {
    stop("Please do not pass seed.use in RunPCA_args. Use the seeed argument.")
  }
  if ("verbose" %in% names(RunPCA_args)) {
    stop("Please do not pass verbose in RunPCA_args. Use the verbose argument.")
  }
  if ("reduction.name" %in% names(RunPCA_args)) {
    stop("Please do not pass reduction.name in RunPCA_args.")
  }
  if ("reduction.key" %in% names(RunPCA_args)) {
    stop("Please do not pass reduction.key in RunPCA_args.")
  }
  if ("assay" %in% names(RunPCA_args)) {
    stop("Please do not pass assay in RunPCA_args.")
  }
}


check_celltype_refs <- function(celltype_refs, celltype_label) {
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
      stop("celltype_label and celltype_refs must have the same lengths. Please choose one label for each ref.")
    }
    for (i in seq_along(celltype_refs)) {
      if (any(!celltype_label[[i]] %in% names(celltype_refs[[i]]@colData@listData))) {
        stop(paste0(paste(celltype_label[[i]][which(!celltype_label[[i]] %in% names(celltype_refs[[i]]@colData@listData))], collapse = ", "), " not found in ", names(celltype_refs)[i], "."))
      }
    }

    if (!requireNamespace("SingleR", quietly = T)) {
      BiocManager::install("SingleR")
    }
    return(celltype_label)
  }
}

check_celltype_ref_clusters <- function(celltype_ref_clusters,
                                        batch_corr,
                                        normalization,
                                        FindClusters_args) {
  if (!is.null(celltype_ref_clusters)) {
    #check if celltype_ref_clusters can exist
    pref <- ifelse(batch_corr == "integration", "integrated", ifelse(normalization == "SCT", "SCT", "RNA"))
    mid <- "_snn_res."
    if (length(FindClusters_args[["resolution"]])) {
      suf <- gsub("\\.0$", "", FindClusters_args[["resolution"]])
    } else {
      suf <- "0.8"
    }
    candidates <- paste0(pref,mid,suf)
    if (!celltype_ref_clusters %in% candidates) {
      stop("celltype_ref_clusters not found in candidates based on batch_corr, normalization and cluster_resolutions: ", paste(candidates, collapse = ", "))
    }
  }
}

check_SO_unprocessed_and_samples <- function(SO_unprocessed,
                                             samples,
                                             batch_corr,
                                             RunHarmony_args) {
  if (methods::is(SO_unprocessed, "list")) {
    SO_unprocessed <- SO_unprocessed
  } else if (methods::is(SO_unprocessed, "character")) {
    if (!file.exists(SO_unprocessed)) {
      stop(paste0(SO_unprocessed, "not found."))
    } else {
      if (!grepl("\\.rds$", SO_unprocessed, ignore.case = T)) {
        stop("SO_unprocessed has to be an .rds file.")
      }
      SO_unprocessed <- readRDS(SO_unprocessed)
    }
  } else {
    stop("SO_unprocessed has to be named list of splitted Seurat objects or a path (character) to an .rds file of those. If it is only one Seurat object (one sample)
         make it a list of length 1.")
  }

  if (methods::is(SO_unprocessed, "Seurat")) {
    # if only one Seurat object is provided
    SO_unprocessed <- list(SO_unprocessed)
    names(SO_unprocessed) <- "sample"
  }

  if (is.null(names(SO_unprocessed))) {
    stop("SO_unprocessed has no names.")
  }


  if (is.null(samples)) {
    samples <- names(SO_unprocessed)
  } else {
    if (any(!samples %in% names(SO_unprocessed))) {
      message("samples not found in SO_unprocessed: ", samples[which(!samples %in% names(SO_unprocessed))])
    }
  }
  samples <- names(SO_unprocessed)[which(grepl(paste(samples, collapse = "|"), names(SO_unprocessed)))]
  SO_unprocessed <- SO_unprocessed[which(names(SO_unprocessed) %in% samples)]

  if (length(SO_unprocessed) > 1 && batch_corr == "harmony" && !"group.by.vars" %in% names(RunHarmony_args)) {
    stop("Please provide one or more group.by.vars from meta.data in RunHarmony_args as a list: ", paste(names(SO_unprocessed[[1]]@meta.data), collapse = ", "), ".")
  } else if (length(SO_unprocessed) > 1 && batch_corr == "harmony" && "group.by.vars" %in% names(RunHarmony_args)) {
    #if (any(RunHarmony_args[["group.by.vars"]] %in% names()))
    # check all SO
  }

  return(list(SO_unprocessed, samples))
}

subset_SO_unprocessed <- function(SO_unprocessed,
                                  cells,
                                  downsample,
                                  min_cells) {

  if (!is.null(cells)) {
    SO_unprocessed <- lapply(SO_unprocessed, function(x) {
      inds <- which(Seurat::Cells(x) %in% cells)
      if (length(inds) == 0) {
        return(NULL)
      }
      return(subset(x, cells = Seurat::Cells(x)[which(Seurat::Cells(x) %in% cells)]))
    })
  }

  if (any(sapply(SO_unprocessed, is.null))) {
    message("No cells found for: ", paste(names(SO_unprocessed)[which(sapply(SO_unprocessed, is.null))], collapse = ", "))
    SO_unprocessed <- SO_unprocessed[which(!sapply(SO_unprocessed, is.null))]
    if (length(SO_unprocessed) == 0) {
      stop("No Seurat objects left after filtering for cells.")
    }
  }

  ## use downsample method from Seurat dev package?!
  # ?Seurat::SketchData()
  if (downsample > 1 && downsample < min_cells) {
    message("downsample set to min_cells.")
    downsample <- min_cells
  } else if (downsample < 1) {
    # check if downsampling forces samples in SO_unprocessed below min_cells
    if (any(sapply(SO_unprocessed, function(x) as.integer(downsample*length(Seurat::Cells(x)))) < min_cells)) {
      message("downsampling leads some samples to drop below min_cells.")
    }
  }

  if (downsample < 1) {
    SO_unprocessed <- lapply(SO_unprocessed, function(x) subset(x, cells = sample(Seurat::Cells(x), as.integer(downsample*length(Seurat::Cells(x))), replace = FALSE)))
  } else if (downsample > 1) {
    SO_unprocessed <- lapply(SO_unprocessed, function(x) subset(x, cells = sample(Seurat::Cells(x), downsample, replace = FALSE)))
  }

  # remove samples with insufficient number of cells
  rm.nm <- names(SO_unprocessed[which(sapply(SO_unprocessed, function(x){length(Seurat::Cells(x)) < min_cells}))])
  if (length(rm.nm) > 0) {
    SO_unprocessed <- SO_unprocessed[which(!names(SO_unprocessed) %in% rm.nm)]
    message("samples removed due to min.cells: ", paste(rm.nm, collapse = ","))
  }

  nsamples <- length(SO_unprocessed)
  message("Samples included (", nsamples, "): ", paste(names(SO_unprocessed), collapse=", "))

  return(SO_unprocessed)
}


check_SCtransform_args <- function(SCtransform_args,
                                   nhvf,
                                   vars.to.regress) {
  # prep SCtransform args once for all
  SCtransform_args <- SCtransform_args[which(!names(SCtransform_args) %in% c("object", "assay", "new.assay.name", "seed.use", "verbose"))]
  if (!"variable.features.n" %in% names(SCtransform_args)) {
    SCtransform_args <- c(SCtransform_args, list(variable.features.n = nhvf))
  }
  if (!"vst.flavor" %in% names(SCtransform_args)) {
    SCtransform_args <- c(SCtransform_args, list(vst.flavor = "v2"))
  }
  if (!"method" %in% names(SCtransform_args)) {
    SCtransform_args <- c(SCtransform_args, list(method = "glmGamPoi"))
  }
  if (!"vars.to.regress" %in% names(SCtransform_args)) {
    SCtransform_args <- c(SCtransform_args, list(vars.to.regress = vars.to.regress))
  }
  return(SCtransform_args)
}

check_FindVariableFeatures_args <- function(FindVariableFeatures_args,
                                            nhvf) {
  # prep FindVariableFeatures args once for all
  FindVariableFeatures_args <- FindVariableFeatures_args[which(!names(FindVariableFeatures_args) %in% c("object", "assay", "verbose"))]
  if (!"nfeatures" %in% names(FindVariableFeatures_args)) {
    FindVariableFeatures_args <- c(FindVariableFeatures_args, list(nfeatures = nhvf))
  }
}

make_so_single <- function(SO_unprocessed,
                           normalization,
                           SCtransform_args,
                           seeed,
                           verbose,
                           var_feature_filter,
                           nhvf,
                           vars.to.regress,
                           FindVariableFeatures_args,
                           scale_RNA_assay_when_SCT,
                           RunPCA_args,
                           npcs) {
  SO <- SO_unprocessed[[1]]
  rm(SO_unprocessed)
  if (normalization == "SCT") {
    SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                               seed.use = seeed,
                                                               verbose = verbose),
                                                          SCtransform_args))

    # remove var features which are to filter
    if (!is.null(var_feature_filter)) {
      SO <- .var_feature_filter_removal(SO = SO,
                                        var_feature_filter = var_feature_filter,
                                        normalization = normalization,
                                        nhvf = nhvf,
                                        vars.to.regress = vars.to.regress,
                                        seeed = seeed,
                                        verbose = verbose,
                                        SCtransform_args = SCtransform_args,
                                        FindVariableFeatures_args = FindVariableFeatures_args)
    }

    SO <- Seurat::NormalizeData(SO, verbose = verbose, assay = "RNA")
    if (scale_RNA_assay_when_SCT) {
      SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
    }
  } else if (normalization == "LogNormalize") {
    SO <- Seurat::NormalizeData(SO, verbose = verbose, assay = "RNA")
    SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                        assay = "RNA",
                                                                        verbose = verbose),
                                                                   FindVariableFeatures_args))
    #SO <- Seurat::FindVariableFeatures(SO, selection.method = "vst", nfeatures = nhvf, verbose = verbose, assay = "RNA")
    if (!is.null(var_feature_filter)) {
      SO <- .var_feature_filter_removal(SO = SO,
                                        var_feature_filter = var_feature_filter,
                                        normalization = normalization,
                                        nhvf = nhvf,
                                        vars.to.regress = vars.to.regress,
                                        seeed = seeed,
                                        verbose = verbose,
                                        SCtransform_args = SCtransform_args,
                                        FindVariableFeatures_args = FindVariableFeatures_args)
    }
    SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
  }
  SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO,
                                                        npcs = npcs,
                                                        seed.use = seeed,
                                                        verbose = verbose),
                                                   RunPCA_args))
  SO <- Seurat::ProjectDim(SO, reduction = "pca", do.center = T, overwrite = F, verbose = verbose)
}

check_FindIntegrationAnchors_args <- function(FindIntegrationAnchors_args,
                                              anchor_features,
                                              k.filter,
                                              k.score) {
  FindIntegrationAnchors_args <- FindIntegrationAnchors_args[which(!names(FindIntegrationAnchors_args) %in% c("object.list", "normalization.method", "verbose"))]
  if (!"anchor.features" %in% names(FindIntegrationAnchors_args)) {
    FindIntegrationAnchors_args <- c(list(anchor_features = anchor_features), FindIntegrationAnchors_args)
  }
  if (!"k.filter" %in% names(FindIntegrationAnchors_args)) {
    FindIntegrationAnchors_args <- c(list(k.filter = k.filter), FindIntegrationAnchors_args)
  }
  if (!"k.score" %in% names(FindIntegrationAnchors_args)) {
    FindIntegrationAnchors_args <- c(list(k.score = k.score), FindIntegrationAnchors_args)
  }
  if (!"reduction" %in% names(FindIntegrationAnchors_args)) {
    FindIntegrationAnchors_args <- c(list(reduction = "rpca"), FindIntegrationAnchors_args)
  }
  return(FindIntegrationAnchors_args)
}

make_so_multi_integrate <- function(SO_unprocessed,
                                    normalization,
                                    SCtransform_args,
                                    seeed,
                                    verbose,
                                    var_feature_filter,
                                    nhvf,
                                    vars.to.regress,
                                    FindVariableFeatures_args,
                                    IntegrateData_args,
                                    hvf_determination_before_merge,
                                    scale_RNA_assay_when_SCT,
                                    npcs) {

  k.filter <- as.integer(min(200, min(sapply(SO_unprocessed, ncol))/2))
  k.score <- as.integer(min(30, min(sapply(SO_unprocessed, ncol))/6))

  if (normalization == "SCT") {
    SO_unprocessed <- lapply(SO_unprocessed, function(x) {
      Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = x,
                                                           seed.use = seeed,
                                                           verbose = verbose),
                                                      SCtransform_args))
    })
    # SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::SCTransform, assay = "RNA", variable.features.n = nhvf, vars.to.regress = vars.to.regress, seed.use = seeed, vst.flavor = "v2", method = "glmGamPoi", verbose = verbose)

    # remove var features which are to filter
    if (!is.null(var_feature_filter)) {
      SO_unprocessed <- lapply(SO_unprocessed,
                               FUN = .var_feature_filter_removal,
                               var_feature_filter = var_feature_filter,
                               normalization = normalization,
                               nhvf = nhvf,
                               vars.to.regress = vars.to.regress,
                               seeed = seeed,
                               verbose = verbose,
                               SCtransform_args = SCtransform_args,
                               FindVariableFeatures_args = FindVariableFeatures_args)
    }

    anchor_features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                         verbose = verbose,
                                                         nfeatures = nhvf)
    SO_unprocessed <- Seurat::PrepSCTIntegration(object.list = SO_unprocessed,
                                                 anchor.features = anchor_features,
                                                 verbose = verbose)
  } else {
    SO_unprocessed <- lapply(SO_unprocessed,
                             FUN = Seurat::NormalizeData,
                             assay = "RNA",
                             verbose = verbose)
    SO_unprocessed <- lapply(SO_unprocessed, function(x) {
      Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = x,
                                                                    assay = "RNA",
                                                                    verbose = verbose),
                                                               FindVariableFeatures_args))
    })
    #SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::FindVariableFeatures, assay = "RNA", selection.method = "vst", nfeatures = nhvf, verbose = verbose)
    anchor_features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                         nfeatures = nhvf,
                                                         verbose = verbose)
  }
  FindIntegrationAnchors_args <- check_FindIntegrationAnchors_args(FindIntegrationAnchors_args,
                                                                   anchor_features,
                                                                   k.filter,
                                                                   k.score)

  if (FindIntegrationAnchors_args[["reduction"]] == "rpca") {
    if (normalization == "SCT") {
      ## does passing RunPCA_args like this cause an error?
      SO_unprocessed <- lapply(SO_unprocessed,
                               FUN = Seurat::RunPCA,
                               assay = "SCT",
                               npcs = npcs,
                               seed.use = seeed,
                               features = anchor_features,
                               verbose = verbose,
                               unlist(RunPCA_args))
    }
  } else {
    SO_unprocessed <- lapply(SO_unprocessed, function(x) {
      x <- Seurat::ScaleData(x, features = anchor_features, verbose = verbose)
      x <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = x,
                                                           npcs = npcs,
                                                           seed.use = seeed,
                                                           verbose = verbose),
                                                      RunPCA_args))
      return(x)
    })
  }

  anchorset <- Gmisc::fastDoCall(Seurat::FindIntegrationAnchors, args = c(list(object.list = SO_unprocessed,
                                                                               normalization.method = normalization,
                                                                               verbose = verbose),
                                                                          FindIntegrationAnchors_args))

  IntegrateData_args <- IntegrateData_args[which(!names(IntegrateData_args) %in% c("anchorset", "new.assay.name", "normalization.method", "verbose"))]
  if (!"features.to.integrate" %in% names(IntegrateData_args)) {
    IntegrateData_args <- c(list(features.to.integrate = rownames(Seurat::GetAssayData(SO_unprocessed[[1]],
                                                                                       assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA")))), IntegrateData_args)
  }
  if (!"k.weight" %in% names(IntegrateData_args)) {
    IntegrateData_args <- c(list(k.weight = k.filter), IntegrateData_args)
  }

  SO <- Gmisc::fastDoCall(Seurat::IntegrateData, args = c(list(anchorset = anchorset,
                                                               normalization.method = normalization,
                                                               verbose = verbose),
                                                          IntegrateData_args))
  Seurat::DefaultAssay(SO) <- "integrated"
  Seurat::VariableFeatures(SO) <- anchor_features

  if (normalization == "SCT" && (hvf_determination_before_merge || batch_corr == "integration")) {
    ## see https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
    SO <- Seurat::PrepSCTFindMarkers(SO,
                                     assay = "SCT",
                                     verbose = verbose)
    message("If running on a subset of the original object after running PrepSCTFindMarkers(), FindMarkers() should be invoked with recorrect_umi = FALSE.")
  }

  ## independent of SCT or LogNormalize:
  if (scale_RNA_assay_when_SCT && normalization == "SCT") {
    SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose) # just for completeness
  }
  SO <- Seurat::ScaleData(SO, assay = "integrated", verbose = verbose) # see https://satijalab.org/seurat/articles/integration_introduction.html

  SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO,
                                                        npcs = npcs,
                                                        seed.use = seeed,
                                                        verbose = verbose),
                                                   RunPCA_args))
  SO <- Seurat::ProjectDim(SO, reduction = "pca", do.center = T, overwrite = F, verbose = verbose)

  return(SO)
}

make_so_multi_harmony <- function(SO_unprocessed,
                                  normalization,
                                  SCtransform_args,
                                  seeed,
                                  verbose,
                                  var_feature_filter,
                                  nhvf,
                                  vars.to.regress,
                                  FindVariableFeatures_args,
                                  scale_RNA_assay_when_SCT,
                                  RunPCA_args,
                                  npcs,
                                  RunHarmony_args,
                                  batch_corr,
                                  hvf_determination_before_merge) {


  ### run SCT separately on unmerged SOs?
  # https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
  # https://github.com/immunogenomics/harmony/issues/41
  # https://github.com/satijalab/sctransform/issues/55#issuecomment-633843730
  ## also see: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

  if (normalization == "SCT") {
    if (hvf_determination_before_merge) {

      SO_unprocessed <- lapply(SO_unprocessed, function(x) {
        Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = x,
                                                             seed.use = seeed,
                                                             verbose = verbose),
                                                        SCtransform_args))
      })
      # SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::SCTransform, assay = "RNA", variable.features.n = nhvf, vars.to.regress = vars.to.regress, seed.use = seeed, vst.flavor = "v2", method = "glmGamPoi", verbose = verbose)

      if (!is.null(var_feature_filter)) {
        SO_unprocessed <- lapply(SO_unprocessed,
                                 FUN = .var_feature_filter_removal,
                                 var_feature_filter = var_feature_filter,
                                 normalization = normalization,
                                 nhvf = nhvf,
                                 vars.to.regress = vars.to.regress,
                                 seeed = seeed,
                                 verbose = verbose,
                                 SCtransform_args = SCtransform_args,
                                 FindVariableFeatures_args = FindVariableFeatures_args)
      }

      anchor_features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                           assay = rep("SCT", length(SO_unprocessed)),
                                                           nfeatures = nhvf,
                                                           fvf.nfeatures = nhvf)

      SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, Seurat::GetAssayData, layer = "counts")),
                                       meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))
      #SO <- merge(x = SO_unprocessed[[1]], y = SO_unprocessed[2:length(SO_unprocessed)], merge.data = T, collapse = T)
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 seed.use = seeed,
                                                                 verbose = verbose),
                                                            SCtransform_args))
      Seurat::VariableFeatures(SO) <- anchor_features
    } else {
      SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, Seurat::GetAssayData, layer = "counts")),
                                       meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))
      #SO <- merge(x = SO_unprocessed[[1]], y = SO_unprocessed[2:length(SO_unprocessed)], merge.data = T, collapse = T)
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 seed.use = seeed,
                                                                 verbose = verbose),
                                                            SCtransform_args))
      #SO <- Seurat::SCTransform(SO, assay = "RNA", variable.features.n = nhvf, vars.to.regress = vars.to.regress, seed.use = seeed, vst.flavor = "v2", method = "glmGamPoi", verbose = verbose)

      # remove var features which are to filter
      if (!is.null(var_feature_filter)) {
        SO <- .var_feature_filter_removal(SO = SO,
                                          var_feature_filter = var_feature_filter,
                                          normalization = normalization,
                                          nhvf = nhvf,
                                          vars.to.regress = vars.to.regress,
                                          seeed = seeed,
                                          verbose = verbose,
                                          SCtransform_args = SCtransform_args,
                                          FindVariableFeatures_args = FindVariableFeatures_args)
      }
    }
    SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = verbose)
    if (scale_RNA_assay_when_SCT) {
      SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = verbose)
    }

  } else if (normalization == "LogNormalize") {

    if (hvf_determination_before_merge) {
      SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::NormalizeData, assay = "RNA", verbose = verbose)

      SO_unprocessed <- lapply(SO_unprocessed, function(x) {
        Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = x,
                                                                      assay = "RNA",
                                                                      verbose = verbose),
                                                                 FindVariableFeatures_args))
      })
      #SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::FindVariableFeatures, assay = "RNA", selection.method = "vst", nfeatures = nhvf, verbose = verbose)

      if (!is.null(var_feature_filter)) {
        SO_unprocessed <- lapply(SO_unprocessed,
                                 FUN = .var_feature_filter_removal,
                                 var_feature_filter = var_feature_filter,
                                 normalization = normalization,
                                 nhvf = nhvf,
                                 vars.to.regress = vars.to.regress,
                                 seeed = seeed,
                                 verbose = verbose,
                                 SCtransform_args = SCtransform_args,
                                 FindVariableFeatures_args = FindVariableFeatures_args)
      }
      anchor_features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                           assay = rep("RNA", length(SO_unprocessed)),
                                                           nfeatures = nhvf,
                                                           verbose = verbose)
      SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, Seurat::GetAssayData, layer = "counts")),
                                       meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))
      SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = verbose)
      Seurat::VariableFeatures(SO) <- anchor_features
    } else {
      SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, Seurat::GetAssayData, layer = "counts")),
                                       meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))
      #SO <- merge(x = SO_unprocessed[[1]], y = SO_unprocessed[2:length(SO_unprocessed)], merge.data = T, collapse = T)
      SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = verbose)
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = "RNA",
                                                                          verbose = verbose),
                                                                     FindVariableFeatures_args))
      #SO <- Seurat::FindVariableFeatures(SO, assay = "RNA", selection.method = "vst", nfeatures = nhvf, verbose = verbose)
      if (!is.null(var_feature_filter)) {
        SO <- .var_feature_filter_removal(SO = SO,
                                          var_feature_filter = var_feature_filter,
                                          normalization = normalization,
                                          nhvf = nhvf,
                                          vars.to.regress = vars.to.regress,
                                          seeed = seeed,
                                          verbose = verbose,
                                          SCtransform_args = SCtransform_args,
                                          FindVariableFeatures_args = FindVariableFeatures_args)
      }
    }
    SO <- Seurat::ScaleData(SO, vars.to.regress = vars.to.regress, verbose = verbose)
  }

  SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO,
                                                        npcs = npcs,
                                                        seed.use = seeed,
                                                        verbose = verbose),
                                                   RunPCA_args))
  SO <- Seurat::ProjectDim(SO, reduction = "pca", do.center = T, overwrite = F, verbose = verbose)

  if (batch_corr == "harmony") {
    RunHarmony_args <- RunHarmony_args[which(!names(RunHarmony_args) %in% c("object", "assay.use", "verbose"))]
    SO <- Gmisc::fastDoCall(harmony::RunHarmony, args = c(list(object = SO,
                                                               assay.use = switch(normalization, SCT = "SCT", LogNormalize = "RNA")),
                                                          RunHarmony_args))
  }
  return(SO)
}

run_celltyping <- function(SO,
                           celltype_refs,
                           celltype_ref_clusters,
                           celltype_label) {
  if (!is.null(celltype_refs) && !is.null(celltype_ref_clusters) && !celltype_ref_clusters %in% names(SO@meta.data)) {
    message("celltype_ref_clusters not found in SO meta data. Will be set to NULL and SingleR will operate on single cell level.")
    celltype_ref_clusters <- NULL
  }

  if (is.null(celltype_ref_clusters)) {
    refs <- NULL
  } else {
    refs <- SO@meta.data[,celltype_ref_clusters]
  }

  layer <- "data"
  if (utils::compareVersion(as.character(SO@version), "4.9.9") == 1) {
    GetAssayData_args <- list(object = SO,
                              layer = layer,
                              assay = "RNA")
  } else {
    GetAssayData_args <- list(object = SO,
                              slot = layer,
                              assay = "RNA")
  }

  for (i in seq_along(celltype_refs)) {
    for (j in seq_along(celltype_label[[i]])) {
      celltypes <- SingleR::SingleR(test = Gmisc::fastDoCall(what = Seurat::GetAssayData, args = GetAssayData_args),
                                    ref = celltype_refs[[i]],
                                    labels = celltype_refs[[i]]@colData@listData[[celltype_label[[i]][j]]],
                                    clusters = refs)

      if (is.null(celltype_ref_clusters)) {
        SO@meta.data[,paste0(names(celltype_refs)[i], "__", celltype_label[[i]][j])] <- celltypes$labels
      } else {
        celltypes_df <- utils::stack(stats::setNames(celltypes$labels, levels(SO@meta.data[,celltype_ref_clusters])))
        names(celltypes_df) <- c(paste0(names(celltype_refs)[i], "__", celltype_label[[i]][j]), celltype_ref_clusters)
        celltypes_df <- tibble::column_to_rownames(dplyr::left_join(tibble::rownames_to_column(SO@meta.data[,celltype_ref_clusters,drop=F], "ID"), celltypes_df, by = celltype_ref_clusters), "ID")
        SO <- Seurat::AddMetaData(SO, celltypes_df[-which(names(celltypes_df) == celltype_ref_clusters)])
      }
      Seurat::Misc(SO, paste0(names(celltype_refs)[i], "__", celltype_label[[i]][j])) <- celltypes
    }
  }
  return(SO)
}
