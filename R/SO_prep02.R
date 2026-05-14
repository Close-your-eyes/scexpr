#' Prepare a somewhat standardized Seurat Object for further analyses
#'
#'
#'
#' @param SO_unprocessed named list of split Seurat objects
#' @param samples names of SO_unprocessed that are to include; if missing all are used
#' @param cells vector of cell names to include; if missing all cells are used
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
#' @param normalization algorithm for normalization of UMIs, RNA and LogNormalize
#' are the same
#' @param batch_corr procedure for batch correction between samples in SO_unprocessed;
#' only relevant if more than 1 sample is passed
#' @param vars.to.regress passed to Seurat::SCTransform or Seurat::ScaleData; will be applied
#' independent of what is set for batch_corr; if batch_corr is set to 'none', then only vars.to.regress is
#' used to regress out a variable in meta.data which may be sample-specific; other then that
#' it may not be meaningful to regress a sample-specific variable and perform batch_corr;
#' rather vars.to.regress may be used in combination with batch_corr to regress percent_mito or so
#' @param seed seed passed to several methods
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
#' deciding whether SCtransform is run on multiple samples separately before merging (set to TRUE) or
#' after merging (set to FALSE); in my experience and when batch_corr = harmony, setting
#' it to FALSE yields better results; see: https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/06a_integration_harmony.md and subsequent links
#' @param FindVariableFeatures_args arguments to Seurat::FindVariableFeatures
#' @param SCtransform_args arguments to Seurat::SCTransform but c("object", "assay", "new.assay.name", "seed.use", "verbose")
#' @param RunUMAP_args arguments to Seurat::RunUMAP
#' @param RunTSNE_args arguments to scexpr::run_fft_tsne
#' @param FindNeighbors_args arguments to Seurat::FindNeighbors
#' @param FindClusters_args arguments to Seurat::FindClusters
#' @param RunHarmony_args arguments to harmony::RunHarmony
#' @param SOM_args arguments to EmbedSOM::SOM
#' @param GQTSOM_args arguments to EmbedSOM::GQTSOM
#' @param EmbedSOM_args arguments to EmbedSOM::EmbedSOM
#' @param FindIntegrationAnchors_args arguments to Seurat::FindIntegrationAnchors
#' @param IntegrateData_args by default features.to.integrate = rownames(Seurat::GetAssayData(SO_unprocessed[[1]], assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA")))
#' @param RunPCA_args arguments to Seurat::RunUMAP
#' @param ... not used yet
#' @param downsample_method how to downsample?
#' @param save_ext save as rds or zap file?
#' @param join_layers run SeuratObject::JoinLayers in the end?
#' @param interactive_varfeat_selection only applies when hvf_determination_before_merge = F
#' @param interactive_varfeat_selection_inds only applies when hvf_determination_before_merge = F
#' @param interactive_pc_selection do conduct interactive PC selection?
#' @param var_feature_set set hvf manually
#' @param use_nn_for_umap
#'
#' @return Seurat Object, as R object and saved to disk as rds file
#' @export
#' @importFrom zeallot %<-%
#'
#' @examples
#' \dontrun{
#' # for future issues: options(future.globals.maxSize = 1000 * 1024^2)
#' }
SO_prep02 <- function(SO_unprocessed,
                      samples = NULL,
                      cells = NULL,
                      min_cells = 50,
                      downsample = 1,
                      downsample_method = c("uniform", "leverage"),
                      export_prefix = "",
                      reductions = c("umap"),
                      nhvf = 800,
                      npcs = 20,
                      normalization = c("SCT", "LogNormalize", "RNA"),
                      hvf_determination_before_merge = F,
                      batch_corr = c("harmony", "integration", "none"),
                      vars.to.regress = NULL,
                      seed = 42,
                      save_path = NULL,
                      save_ext = c("rds", "zap"),
                      celltype_refs = NULL, # list of celldex::objects of seurat or matrix
                      celltype_label = c("label.main", "label.fine"),
                      celltype_ref_clusters = NULL,
                      diet_seurat = F,
                      var_feature_filter = NULL,
                      var_feature_set = NULL,
                      verbose = F,
                      FindVariableFeatures_args = list(),
                      SCtransform_args = list(
                        vst.flavor = "v2",
                        method = "glmGamPoi",
                        conserve.memory = T
                      ),
                      RunPCA_args = list(weight.by.var = T),
                      RunUMAP_args = list(metric = "cosine"),
                      RunTSNE_args = list(theta = 0.01),
                      FindNeighbors_args = list(dims = 1:npcs),
                      FindClusters_args = list(resolution = seq(0.1,0.8,0.1)),
                      RunHarmony_args = list(group.by.vars = "orig.ident"),
                      SOM_args = list(),
                      GQTSOM_args = list(),
                      EmbedSOM_args = list(),
                      FindIntegrationAnchors_args = list(reduction = "rpca"),
                      IntegrateData_args = list(),
                      join_layers = T,
                      interactive_varfeat_selection = T,
                      interactive_varfeat_selection_inds = seq(max(nhvf/10,50),
                                                               min(3*nhvf, nrow(SO_unprocessed[[1]])),
                                                               length.out = 11),
                      interactive_pc_selection = T,
                      use_nn_for_umap = F,
                      ...) {

  pkg_checks()
  mydots <- list(...)
  options(warn = 1)
  options(future.globals.maxSize = 20 * 1024^3)

  # options(parallelly.availableCores.custom = function() {
  #   ncores <- max(parallel::detectCores(), 1L, na.rm = TRUE)
  #   ncores <- min(as.integer(0.75 * ncores), 4L)
  #   max(1L, ncores)
  # })

  reductions <- match.arg(tolower(reductions), c("umap", "som", "gqtsom", "tsne"), several.ok = T)
  normalization <- rlang::arg_match(normalization)
  normalization <- ifelse(normalization == "RNA", "LogNormalize", normalization)
  batch_corr <- rlang::arg_match(batch_corr)
  downsample_method <- rlang::arg_match(downsample_method)
  save_ext <- rlang::arg_match(save_ext)
  var_feature_set <- unique(var_feature_set)

  if (is.null(save_path)) {
    message("No save_path provided.")
  } else if (!is.character(save_path) || length(save_path) != 1) {
    stop("save_path has to be a character; a path to a folder where to save Seurat objects to.")
  }

  if (interactive_pc_selection) {
    npcs <- max(50, npcs)
  }

  celltype_label <- check_celltype_refs(celltype_refs = celltype_refs,
                                        celltype_label = celltype_label)

  RunPCA_args <- check_RunPCA_args(#obj = SO,
                                   RunPCA_args = RunPCA_args,
                                   normalization = normalization,
                                   npcs = npcs,
                                   seed = seed,
                                   verbose = verbose)

  RunHarmony_args <- check_RunHarmony_args(RunHarmony_args = RunHarmony_args,
                                           RunPCA_args = RunPCA_args)

  check_celltype_ref_clusters(celltype_ref_clusters,
                              batch_corr,
                              normalization,
                              FindClusters_args)

  c(SO_unprocessed, samples) %<-% check_SO_unprocessed_and_samples(SO_unprocessed = SO_unprocessed,
                                                                   samples = samples,
                                                                   batch_corr = batch_corr,
                                                                   RunHarmony_args = RunHarmony_args,
                                                                   verbose = verbose,
                                                                   cells = cells,
                                                                   downsample = downsample,
                                                                   min_cells = min_cells,
                                                                   downsample_method = downsample_method)

  if (length(SO_unprocessed) == 1) {batch_corr <- "none"}

  SCtransform_args <- check_SCtransform_args(SCtransform_args = SCtransform_args,
                                             nhvf = nhvf,
                                             vars.to.regress = vars.to.regress,
                                             verbose = verbose,
                                             seed = seed)

  FindVariableFeatures_args <- check_FindVariableFeatures_args(FindVariableFeatures_args = FindVariableFeatures_args,
                                                               nhvf = nhvf,
                                                               verbose = verbose)

  # cases:
  if (length(SO_unprocessed) == 1) {
    # %<-%
    c(SO,
      RunPCA_args) %<-% make_so_single(SO_unprocessed = SO_unprocessed,
                                       normalization = normalization,
                                       SCtransform_args = SCtransform_args,
                                       var_feature_filter = var_feature_filter,
                                       vars.to.regress = vars.to.regress,
                                       FindVariableFeatures_args = FindVariableFeatures_args,
                                       RunPCA_args = RunPCA_args,
                                       interactive_varfeat_selection = interactive_varfeat_selection,
                                       interactive_varfeat_selection_inds = interactive_varfeat_selection_inds,
                                       interactive_pc_selection = interactive_pc_selection,
                                       var_feature_set = var_feature_set)
  } else if (length(SO_unprocessed) > 1) {
    if (batch_corr %in% c("none", "harmony")) {
      c(SO,
        RunPCA_args,
        RunHarmony_args) %<-% make_so_multi_harmony(SO_unprocessed = SO_unprocessed,
                                                    normalization = normalization,
                                                    SCtransform_args = SCtransform_args,
                                                    var_feature_filter = var_feature_filter,
                                                    vars.to.regress = vars.to.regress,
                                                    FindVariableFeatures_args = FindVariableFeatures_args,
                                                    RunPCA_args = RunPCA_args,
                                                    RunHarmony_args = RunHarmony_args,
                                                    batch_corr = batch_corr,
                                                    hvf_determination_before_merge = hvf_determination_before_merge,
                                                    interactive_varfeat_selection = interactive_varfeat_selection,
                                                    interactive_varfeat_selection_inds = interactive_varfeat_selection_inds,
                                                    interactive_pc_selection = interactive_pc_selection,
                                                    var_feature_set = var_feature_set)
    } else if (batch_corr == "integration") {
      c(SO,
        RunPCA_args) %<-% make_so_multi_integrate(SO_unprocessed = SO_unprocessed,
                                                  normalization = normalization,
                                                  SCtransform_args = SCtransform_args,
                                                  var_feature_filter = var_feature_filter,
                                                  vars.to.regress = vars.to.regress,
                                                  FindVariableFeatures_args = FindVariableFeatures_args,
                                                  IntegrateData_args = IntegrateData_args,
                                                  hvf_determination_before_merge = hvf_determination_before_merge,
                                                  FindIntegrationAnchors_args = FindIntegrationAnchors_args,
                                                  RunPCA_args = RunPCA_args,
                                                  interactive_varfeat_selection = interactive_varfeat_selection,
                                                  interactive_varfeat_selection_inds = interactive_varfeat_selection_inds,
                                                  interactive_pc_selection = interactive_pc_selection,
                                                  var_feature_set = var_feature_set)
    }
  }

  rm(SO_unprocessed)

  ### do.call on large SeuratObject became super slow, not practicable!
  # https://stackoverflow.com/questions/28198103/alternative-to-do-call-for-large-datasets

  # alternatives:
  # https://rlang.r-lib.org/reference/exec.html
  # gmisc::fastDoCall


  red <- switch(
    batch_corr,
    harmony = RunHarmony_args[["reduction.save"]],
    integration = RunPCA_args[["reduction.name"]],
    none = RunPCA_args[["reduction.name"]]
  )
  names_wo_clust <- names(SO@meta.data)

  SO <- calc_neighbor_and_cluster(SO = SO,
                                  red = red,
                                  npcs = RunPCA_args[["npcs"]],
                                  normalization = normalization,
                                  FindNeighbors_args = FindNeighbors_args,
                                  FindClusters_args = FindClusters_args,
                                  verbose = verbose,
                                  mc.cores = 8)

  ## pick clustering with decent cluster number
  ## derive cluster markers
  ## re-define hvf and re-run pca / harmony
  ## no, do it outside


  if (any(grepl("umap", reductions, ignore.case = T))) {
    RunUMAP_args <- RunUMAP_args[which(!names(RunUMAP_args) %in% c("object", "seed.use", "reduction", "verbose"))]
    if (!"dims" %in% names(RunUMAP_args)) {
      RunUMAP_args <- c(list(dims = 1:RunPCA_args[["npcs"]]), RunUMAP_args)
    }

    tryCatch(expr = {
      if (!use_nn_for_umap) {
        SO <- Gmisc::fastDoCall(Seurat::RunUMAP, args = c(list(object = SO,
                                                               reduction = red,
                                                               reduction.name = paste0("umap_", red),
                                                               seed.use = seed,
                                                               verbose = verbose),
                                                          RunUMAP_args))
      } else {
        RunUMAP_args <- RunUMAP_args[which(names(RunUMAP_args) %in% names(formals(uwot::umap)))]
        RunUMAP_args <- c(list(X = NULL,
                               nn_method = list(idx = SO@misc$nn.ranked@nn.idx, dist = SO@misc$nn.ranked@nn.dist),
                               seed = seed,
                               verbose = verbose),
                          RunUMAP_args)
        if (!"metric" %in% names(RunUMAP_args)) {
          RunUMAP_args <- c(list(metric = "cosine"), RunUMAP_args)
        }
        if (!"n_neighbors" %in% names(RunUMAP_args)) {
          RunUMAP_args <- c(list(n_neighbors = 30), RunUMAP_args)
        }
        if (!"min_dist" %in% names(RunUMAP_args)) {
          RunUMAP_args <- c(list(min_dist = 0.3), RunUMAP_args)
        }

        RunUMAP_args <- RunUMAP_args[which(!duplicated(names(RunUMAP_args)))]
        um <- Gmisc::fastDoCall(uwot::umap, RunUMAP_args)
        rownames(um) <- Seurat::Cells(SO)
        colnames(um) <- paste0("umap", gsub("[^A-Za-z1-9]", "", red), "_", c(1,2))
        SO@reductions[[paste0("umap_", red)]] <- SeuratObject::CreateDimReducObject(embeddings = um, assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA", RNA = "RNA"))
      }

    }, error = function(err) {
      message("umap failed.")
    })
  }

  if (any(grepl("tsne", reductions, ignore.case = T))) {
    RunTSNE_args <- RunTSNE_args[which(!names(RunTSNE_args) %in% c("object", "seed.use", "reduction", "verbose"))]

    tryCatch(expr = {
      SO <- Gmisc::fastDoCall(scexpr::run_fft_tsne, args = c(list(SO = SO,
                                                                  reduction = red,
                                                                  reduction.name = paste0("tsne_", red),
                                                                  rand_seed = seed),
                                                             RunTSNE_args))
    }, error = function(err) {
      message("fallback to barnes-hut tsne; default seurat method.")
      if (!"dims" %in% names(RunTSNE_args)) {
        RunTSNE_args <- c(list(dims = 1:RunPCA_args[["npcs"]]), RunTSNE_args)
      }
      if (!"num_threads" %in% names(RunTSNE_args)) {
        RunTSNE_args <- c(list(num_threads = 0), RunTSNE_args)
      }
      #RunTSNE_args[["tsne.method"]] <- "FIt-SNE"
      SO <- Gmisc::fastDoCall(Seurat::RunTSNE, args = c(list(object = SO,
                                                             reduction = red,
                                                             reduction.name = paste0("tsne_", red),
                                                             seed.use = seed,
                                                             verbose = verbose),
                                                        RunTSNE_args))
    })

    #SO <- Seurat::RunTSNE(object = SO, dims = 1:npcs, seed.use = seed, reduction = red, verbose = verbose, num_threads = 0, ...)
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
    SO[["SOM"]] <- Seurat::CreateDimReducObject(embeddings = ES, key = "SOM_", assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA", RNA = "RNA"), misc = list(SOM_args, EmbedSOM_args))
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
    SO[["GQTSOM"]] <- Seurat::CreateDimReducObject(embeddings = ES, key = "GQTSOM_", assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA", RNA = "RNA"), misc = list(GQTSOM_args, EmbedSOM_args))
  }



  # add clusterings and their colors to Misc
  names_w_clust <- names(SO@meta.data)
  SeuratObject::Misc(SO, "clusterings") <- setdiff(names_w_clust, names_wo_clust)
  try(expr = {
    all_cluster <- unique(unname(unlist(SO@meta.data[,SeuratObject::Misc(SO, "clusterings")])))
    # all_cluster <- all_cluster[order(as.numeric(all_cluster))]
    all_cluster <- sort(all_cluster)
    SeuratObject::Misc(SO, "clustering_colors") <- stats::setNames(as.character(colrr::col_pal("custom", n = length(all_cluster))), nm = all_cluster)
  }, silent = T)
  origids <- sort(unique(SO@meta.data$orig.ident))
  SeuratObject::Misc(SO, "orig.ident_colors") <- stats::setNames(as.character(colrr::col_pal("custom", n = length(origids))), nm = origids)

  # add meta cols
  SO@meta.data$id <- rownames(SO@meta.data)
  try(expr = {
    # when metacols exist from SO_prep01 rownames_to_col (tibble) throws error
    newmeta <- SO@meta.data[,"id",drop = F]
    rownames(newmeta) <- newmeta$id
    newmeta$barcode <- stringr::str_extract(newmeta$id, "[ATCG]{1,}-1$")
    newmeta$barcode <- dplyr::coalesce(newmeta$barcode, stringr::str_extract(newmeta$id, "[ATCG]{1,}$"))
    newmeta$prefix <- stringr::str_replace(newmeta$id, newmeta$barcode, "")
    newmeta$prefix <- gsub("_{1,}$", "", newmeta$prefix)
    newmeta$prefix <- gsub("\\.{1,}$", "", newmeta$prefix)
    newmeta$prefix <- gsub("-{1,}$", "", newmeta$prefix)
    newmeta <- newmeta[,-which(names(newmeta) == "id")]
    SO <- SeuratObject::AddMetaData(SO, newmeta)
  }, silent = T)




  SO <- run_celltyping(SO,
                       celltype_refs,
                       celltype_ref_clusters,
                       celltype_label)

  if (join_layers) {
    for (i in names(SO@assays)) {
      try({
        SO <- SeuratObject::JoinLayers(SO, assay = i)
      }, silent = T)
    }
  }

  # remove counts as they can be recalculated with rev_lognorm
  Seurat::Misc(SO, slot = "RNA_count_colSums") <- Matrix::colSums(get_layer(obj = SO, layer = "counts", assay = "RNA"))
  if (diet_seurat) {
    SO <- Seurat::DietSeurat(SO, assays = names(SO@assays), counts = F, dimreducs = names(SO@reductions))
  }

  try(expr = {
    # quick to calculate
    # save disk space
    SO@assays[["RNA"]]@layers[["scale.data"]] <- NULL
    # SO@assays[["RNA"]]@layers[["scale.data"]] <- matrix(NA)
  }, silent = T)



  save.time <- format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%y%m%d_%H%M%S")
  save.name <- paste(
    "SO",
    export_prefix,
    normalization,
    batch_corr,
    downsample,
    length(Seurat::VariableFeatures(SO)),
    RunPCA_args[["npcs"]],
    paste0(save.time, ".", save_ext),
    sep = "_"
  )
  save.name <- gsub("__", "_", save.name)
  Seurat::Misc(SO, "object_name") <- save.name

  if (!is.null(save_path)) {
    if (verbose) {
      message("saving to disk.")
    }
    dir.create(save_path, showWarnings = F, recursive = T)
    if (save_ext == "rds") {
      saveRDS(SO, compress = T, file = file.path(save_path, save.name))
    } else if (save_ext == "zap") {
      zap::zap_write(SO, dst = file.path(save_path, save.name), compress = "zstd")
    }
    message("SO saved to: ")
    message(file.path(save_path, save.name))
  }

  return(SO)
}

make_so_single <- function(SO_unprocessed,
                           normalization,
                           SCtransform_args,
                           var_feature_filter,
                           vars.to.regress,
                           FindVariableFeatures_args,
                           RunPCA_args,
                           interactive_varfeat_selection,
                           interactive_varfeat_selection_inds,
                           interactive_pc_selection,
                           var_feature_set) {

  SO <- SO_unprocessed[[1]]
  rm(SO_unprocessed)
  if (normalization == "SCT") {

    if (!is.null(var_feature_set)) {
      Seurat::VariableFeatures(SO) <- var_feature_set
    }

    if (interactive_varfeat_selection && is.null(var_feature_set)) {
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO),
                                                                     FindVariableFeatures_args))
      vfplot <- varfeat_plot(SO, n_varfeat = interactive_varfeat_selection_inds)
      print(vfplot)
      SCtransform_args[["variable.features.n"]] <- get_numeric_input("select number of variable features.")
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 assay = "RNA"),
                                                            SCtransform_args))
    } else {
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 assay = "RNA"),
                                                            SCtransform_args))
      varfeat_plot(obj = SO, n_varfeat = SCtransform_args[["variable.features.n"]])
    }

    # remove var features which are to filter
    if (!is.null(var_feature_filter)) {
      SO <- var_feature_filter_removal(SO = SO,
                                       var_feature_filter = var_feature_filter,
                                       normalization = normalization,
                                       vars.to.regress = vars.to.regress,
                                       SCtransform_args = SCtransform_args,
                                       FindVariableFeatures_args = FindVariableFeatures_args)
    }

    SO <- Seurat::NormalizeData(SO, verbose = SCtransform_args[["verbose"]], assay = "RNA")

  } else if (normalization %in% c("LogNormalize", "RNA")) {

    SO <- Seurat::NormalizeData(SO, verbose = FindVariableFeatures_args[["verbose"]], assay = "RNA")

    SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                        assay = "RNA"),
                                                                   FindVariableFeatures_args))
    if (!is.null(var_feature_set)) {
      Seurat::VariableFeatures(SO) <- var_feature_set
    }

    if (interactive_varfeat_selection && is.null(var_feature_set)) {
      vfplot <- varfeat_plot(SO, n_varfeat = interactive_varfeat_selection_inds)
      print(vfplot)
      FindVariableFeatures_args[["nfeatures"]] <- get_numeric_input("select number of variable features.")
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = "RNA"),
                                                                     FindVariableFeatures_args))
    } else {
      varfeat_plot(obj = SO, n_varfeat = FindVariableFeatures_args[["nfeatures"]])
    }

    if (!is.null(var_feature_filter)) {
      SO <- var_feature_filter_removal(SO = SO,
                                       var_feature_filter = var_feature_filter,
                                       normalization = normalization,
                                       vars.to.regress = vars.to.regress,
                                       SCtransform_args = SCtransform_args,
                                       FindVariableFeatures_args = FindVariableFeatures_args)
    }
    SO <- Seurat::ScaleData(SO, assay = "RNA", verbose = FindVariableFeatures_args[["verbose"]])
  }

  SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))

  if (interactive_pc_selection) {
    print(scexpr::elbowplot2(SO, npcs = RunPCA_args[["npcs"]])[["plot"]])
    RunPCA_args[["npcs"]] <- get_numeric_input("select number of pc.")

    RunPCA_args <- check_RunPCA_args(#obj = SO,
      RunPCA_args = RunPCA_args,
      normalization = normalization,
      npcs = RunPCA_args[["npcs"]],
      seed = RunPCA_args[["seed"]],
      verbose = RunPCA_args[["verbose"]])

    SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))
  }
  # SO <- Seurat::ProjectDim(SO, reduction = "pca", do.center = T, overwrite = F, verbose = verbose)

  return(list(SO = SO, RunPCA_args = RunPCA_args))
}



make_so_multi_integrate <- function(SO_unprocessed,
                                    normalization,
                                    SCtransform_args,
                                    var_feature_filter,
                                    vars.to.regress,
                                    FindVariableFeatures_args,
                                    IntegrateData_args,
                                    hvf_determination_before_merge,
                                    FindIntegrationAnchors_args,
                                    RunPCA_args,
                                    interactive_varfeat_selection,
                                    interactive_varfeat_selection_inds,
                                    interactive_pc_selection,
                                    var_feature_set) {

  k.filter <- as.integer(min(200, min(sapply(SO_unprocessed, ncol))/2))
  k.score <- as.integer(min(30, min(sapply(SO_unprocessed, ncol))/6))

  if (normalization == "SCT") {
    SO_unprocessed <- lapply(SO_unprocessed, function(x) {
      x <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = x,
                                                                assay = "RNA"),
                                                           SCtransform_args))
      if (!is.null(var_feature_set)) {
        Seurat::VariableFeatures(x) <- var_feature_set
      }
      return(x)
    })

    if (interactive_varfeat_selection && is.null(var_feature_set)) {
      for (i in names(SO_unprocessed)) {
        vfplot <- varfeat_plot(SO_unprocessed[[i]], n_varfeat = interactive_varfeat_selection_inds) +
          patchwork::plot_annotation(title = i)
        print(vfplot)
      }
      SCtransform_args[["variable.features.n"]] <- get_numeric_input("select number of variable features.")
      SO_unprocessed <- lapply(SO_unprocessed, function(x) {
        Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = x,
                                                             assay = "RNA"),
                                                        SCtransform_args))
      })
    } else {
      for (i in names(SO_unprocessed)) {
        varfeat_plot(SO_unprocessed[[i]], n_varfeat = interactive_varfeat_selection_inds) +
          patchwork::plot_annotation(title = i)
      }
    }

    # remove var features which are to filter
    if (!is.null(var_feature_filter)) {
      SO_unprocessed <- lapply(SO_unprocessed,
                               FUN = var_feature_filter_removal,
                               var_feature_filter = var_feature_filter,
                               normalization = normalization,
                               vars.to.regress = vars.to.regress,
                               SCtransform_args = SCtransform_args,
                               FindVariableFeatures_args = FindVariableFeatures_args)
    }

    anchor.features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                         verbose = SCtransform_args[["verbose"]],
                                                         nfeatures = SCtransform_args[["variable.features.n"]])
    SO_unprocessed <- Seurat::PrepSCTIntegration(object.list = SO_unprocessed,
                                                 anchor.features = anchor.features,
                                                 verbose = SCtransform_args[["verbose"]])
  } else {

    SO_unprocessed <- lapply(SO_unprocessed,
                             FUN = Seurat::NormalizeData,
                             assay = "RNA",
                             verbose = FindVariableFeatures_args[["verbose"]])
    SO_unprocessed <- lapply(SO_unprocessed, function(x) {

      x <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = x,
                                                                         assay = "RNA"),
                                                                    FindVariableFeatures_args))
      if (!is.null(var_feature_set)) {
        Seurat::VariableFeatures(x) <- var_feature_set
      }
      return(x)
    })

    if (interactive_varfeat_selection && is.null(var_feature_set)) {
      for (i in names(SO_unprocessed)) {
        vfplot <- varfeat_plot(SO_unprocessed[[i]], n_varfeat = interactive_varfeat_selection_inds) +
          patchwork::plot_annotation(title = i)
        print(vfplot)
      }
      FindVariableFeatures_args[["nfeatures"]] <- get_numeric_input("select number of variable features.")
      SO_unprocessed <- lapply(SO_unprocessed, function(x) {
        Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = x,
                                                                      assay = "RNA"),
                                                                 FindVariableFeatures_args))
      })
    } else {
      for (i in names(SO_unprocessed)) {
        varfeat_plot(SO_unprocessed[[i]], n_varfeat = interactive_varfeat_selection_inds) +
          patchwork::plot_annotation(title = i)
      }
    }

    anchor.features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                         nfeatures = FindVariableFeatures_args[["nfeatures"]],
                                                         verbose = FindVariableFeatures_args[["verbose"]])
  }

  FindIntegrationAnchors_args <- check_FindIntegrationAnchors_args(FindIntegrationAnchors_args,
                                                                   anchor.features,
                                                                   k.filter,
                                                                   k.score,
                                                                   npcs)



  # if (FindIntegrationAnchors_args[["reduction"]] == "rpca") {
  #   if (normalization == "SCT") {
  #     SO_unprocessed <- lapply(SO_unprocessed, function(x) {
  #       x <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = x,
  #                                                            npcs = npcs,
  #                                                            seed.use = seed,
  #                                                            assay = "SCT",
  #                                                            features = anchor.features,
  #                                                            verbose = verbose),
  #                                                       RunPCA_args))
  #       return(x)
  #     })
  #   }
  # } else {
  #   SO_unprocessed <- lapply(SO_unprocessed, function(x) {
  #     x <- Seurat::ScaleData(x, features = anchor.features, verbose = verbose)
  #     x <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = x,
  #                                                          npcs = npcs,
  #                                                          seed.use = seed,
  #                                                          verbose = verbose),
  #                                                     RunPCA_args))
  #     return(x)
  #   })
  # }

  SO_unprocessed <- lapply(SO_unprocessed, function(x) {
    x <- Seurat::ScaleData(
      x,
      features = anchor.features,
      assay = "RNA",
      verbose = FindVariableFeatures_args[["verbose"]]
    )
    x <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = x,
                                                         features = anchor.features,
                                                         assay = ifelse(normalization == "SCT", "SCT", "RNA")),
                                                    RunPCA_args))
    return(x)
  })

  anchorset <- Gmisc::fastDoCall(Seurat::FindIntegrationAnchors, args = c(list(object.list = SO_unprocessed,
                                                                               normalization.method = normalization,
                                                                               verbose = FindVariableFeatures_args[["verbose"]]),
                                                                          FindIntegrationAnchors_args))

  IntegrateData_args <- IntegrateData_args[which(!names(IntegrateData_args) %in% c("anchorset", "new.assay.name", "normalization.method", "verbose"))]
  if (!"features.to.integrate" %in% names(IntegrateData_args)) {
    IntegrateData_args <- c(list(features.to.integrate = rownames(get_layer(obj = SO_unprocessed[[1]],
                                                                            assay = switch(normalization, SCT = "SCT", LogNormalize = "RNA", RNA = "RNA")))),
                            IntegrateData_args)
  }
  if (!"k.weight" %in% names(IntegrateData_args)) {
    IntegrateData_args <- c(list(k.weight = k.filter), IntegrateData_args)
  }

  SO <- Gmisc::fastDoCall(Seurat::IntegrateData, args = c(list(anchorset = anchorset,
                                                               normalization.method = normalization,
                                                               verbose = FindVariableFeatures_args[["verbose"]]),
                                                          IntegrateData_args))
  Seurat::DefaultAssay(SO) <- "integrated"
  Seurat::VariableFeatures(SO) <- anchor.features

  if (normalization == "SCT" && hvf_determination_before_merge) { # || batch_corr == "integration"
    ## see https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
    SO <- Seurat::PrepSCTFindMarkers(SO,
                                     assay = "SCT",
                                     verbose = SCtransform_args[["verbose"]])
    message("If running on a subset of the original object after running PrepSCTFindMarkers(), FindMarkers() should be invoked with recorrect_umi = FALSE.")
  }

  SO <- Seurat::ScaleData(SO, assay = "integrated", verbose = SCtransform_args[["verbose"]]) # see https://satijalab.org/seurat/articles/integration_introduction.html

  SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))

  if (interactive_pc_selection) {
    print(scexpr::elbowplot2(SO, npcs = RunPCA_args[["npcs"]])[["plot"]])
    RunPCA_args[["npcs"]] <- get_numeric_input("select number of pc.")

    RunPCA_args <- check_RunPCA_args(#obj = SO,
      RunPCA_args = RunPCA_args,
      normalization = normalization,
      npcs = RunPCA_args[["npcs"]],
      seed = RunPCA_args[["seed"]],
      verbose = RunPCA_args[["verbose"]])

    SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))
  }
  # SO <- Seurat::ProjectDim(SO, reduction = "pca", do.center = T, overwrite = F, verbose = verbose)

  return(list(SO = SO, RunPCA_args = RunPCA_args))
}

make_so_multi_harmony <- function(SO_unprocessed,
                                  normalization,
                                  SCtransform_args,
                                  var_feature_filter,
                                  vars.to.regress,
                                  FindVariableFeatures_args,
                                  RunPCA_args,
                                  RunHarmony_args,
                                  batch_corr,
                                  hvf_determination_before_merge,
                                  interactive_varfeat_selection,
                                  interactive_varfeat_selection_inds,
                                  interactive_pc_selection,
                                  var_feature_set) {


  ### run SCT separately on unmerged SOs?
  # https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
  # https://github.com/immunogenomics/harmony/issues/41
  # https://github.com/satijalab/sctransform/issues/55#issuecomment-633843730
  ## also see: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

  if (normalization == "SCT") {
    if (hvf_determination_before_merge) {

      SO_unprocessed <- lapply(SO_unprocessed, function(x) {
        x <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = x,
                                                                  assay = "RNA"),
                                                             SCtransform_args))

        if (!is.null(var_feature_set)) {
          Seurat::VariableFeatures(x) <- var_feature_set
        }
        return(x)
      })

      if (!is.null(var_feature_filter)) {
        SO_unprocessed <- lapply(SO_unprocessed,
                                 FUN = var_feature_filter_removal,
                                 var_feature_filter = var_feature_filter,
                                 normalization = normalization,
                                 vars.to.regress = vars.to.regress,
                                 SCtransform_args = SCtransform_args,
                                 FindVariableFeatures_args = FindVariableFeatures_args)
      }

      anchor.features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                           assay = rep("SCT", length(SO_unprocessed)),
                                                           nfeatures = SCtransform_args[["variable.features.n"]],
                                                           fvf.nfeatures = SCtransform_args[["variable.features.n"]])

      # SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, get_layer, layer = "counts")),
      #                                  meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))

      ## better for memory
      ## and much faster
      meta.data <- purrr::map_dfr(SO_unprocessed, ~.x@meta.data)
      SO_unprocessed <- purrr::map(SO_unprocessed,
                                   get_layer,
                                   layer = "counts")
      SO <- Seurat::CreateSeuratObject(counts = chunk_wise_cbind(x = SO_unprocessed),
                                       meta.data = meta.data)
      rm(SO_unprocessed)



      #SO <- merge(x = SO_unprocessed[[1]], y = SO_unprocessed[2:length(SO_unprocessed)], merge.data = T, collapse = T)
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 assay = "RNA"),
                                                            SCtransform_args))
      Seurat::VariableFeatures(SO) <- anchor.features
    } else {

      intersect_features <- purrr::reduce(purrr::map(SO_unprocessed, rownames), intersect)
      feature_numbers <- purrr::map_int(SO_unprocessed, nrow)
      feature_number_common <- length(intersect_features)
      if (any(feature_number_common != feature_numbers)) {
        message("unequal features across inputs:")
        print(feature_numbers)
        message("common features: ", feature_number_common)
      }

      # SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed,
      #                                                                     get_layer,
      #                                                                     layer = "counts",
      #                                                                     features = intersect_features)),
      #                                  meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))

      ## better for memory
      ## and much faster
      meta.data <- purrr::map_dfr(SO_unprocessed, ~.x@meta.data)
      SO_unprocessed <- purrr::map(SO_unprocessed,
                                   get_layer,
                                   layer = "counts",
                                   features = intersect_features)
      SO <- Seurat::CreateSeuratObject(counts = chunk_wise_cbind(x = SO_unprocessed),
                                       meta.data = meta.data)
      rm(SO_unprocessed)





      ## dont run SCTransform twice.
      # future::plan(strategy = "multisession", workers = 2)
      # future::plan(strategy = "sequential") # was faster
      # SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
      #                                                            assay = "RNA"),
      #                                                       SCtransform_args))

      if (!is.null(var_feature_set)) {
        Seurat::VariableFeatures(SO) <- var_feature_set
      }

      if (interactive_varfeat_selection && is.null(var_feature_set)) {
        SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO),
                                                                       FindVariableFeatures_args))
        vfplot <- varfeat_plot(SO, n_varfeat = interactive_varfeat_selection_inds)
        print(vfplot)
        SCtransform_args[["variable.features.n"]] <- get_numeric_input("select number of variable features.")
        SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                   assay = "RNA"),
                                                              SCtransform_args))
      } else {
        SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                   assay = "RNA"),
                                                              SCtransform_args))
        varfeat_plot(obj = SO, n_varfeat = SCtransform_args[["variable.features.n"]])
      }

      # remove var features which are to filter
      if (!is.null(var_feature_filter)) {
        SO <- var_feature_filter_removal(SO = SO,
                                         var_feature_filter = var_feature_filter,
                                         normalization = normalization,
                                         vars.to.regress = vars.to.regress,
                                         SCtransform_args = SCtransform_args,
                                         FindVariableFeatures_args = FindVariableFeatures_args)
      }
    }
    SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = FindVariableFeatures_args[["verbose"]])

  } else if (normalization == "LogNormalize") {

    if (hvf_determination_before_merge) {
      SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::NormalizeData, assay = "RNA", verbose = FindVariableFeatures_args[["verbose"]])

      SO_unprocessed <- lapply(SO_unprocessed, function(x) {

        x <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = x,
                                                                           assay = "RNA"),
                                                                      FindVariableFeatures_args))
        if (!is.null(var_feature_set)) {
          Seurat::VariableFeatures(x) <- var_feature_set
        }
        return(x)
      })
      #SO_unprocessed <- lapply(SO_unprocessed, FUN = Seurat::FindVariableFeatures, assay = "RNA", selection.method = "vst", nfeatures = nhvf, verbose = verbose)

      if (!is.null(var_feature_filter)) {
        SO_unprocessed <- lapply(SO_unprocessed,
                                 FUN = var_feature_filter_removal,
                                 var_feature_filter = var_feature_filter,
                                 normalization = normalization,
                                 vars.to.regress = vars.to.regress,
                                 SCtransform_args = SCtransform_args,
                                 FindVariableFeatures_args = FindVariableFeatures_args)
      }
      anchor.features <- Seurat::SelectIntegrationFeatures(SO_unprocessed,
                                                           assay = rep("RNA", length(SO_unprocessed)),
                                                           nfeatures = FindVariableFeatures_args[["nfeatures"]],
                                                           verbose = FindVariableFeatures_args[["verbose"]])



      # SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, get_layer, layer = "counts")),
      #                                  meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))
      meta.data <- purrr::map_dfr(SO_unprocessed, ~.x@meta.data)
      SO_unprocessed <- purrr::map(SO_unprocessed,
                                   get_layer,
                                   layer = "counts")
      SO <- Seurat::CreateSeuratObject(counts = chunk_wise_cbind(x = SO_unprocessed),
                                       meta.data = meta.data)
      rm(SO_unprocessed)

      SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = FindVariableFeatures_args[["verbose"]])
      Seurat::VariableFeatures(SO) <- anchor.features
    } else {
      # SO <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(SO_unprocessed, get_layer, layer = "counts")),
      #                                  meta.data = dplyr::bind_rows(purrr::map(SO_unprocessed, ~.x@meta.data)))
      meta.data <- purrr::map_dfr(SO_unprocessed, ~.x@meta.data)
      SO_unprocessed <- purrr::map(SO_unprocessed,
                                   get_layer,
                                   layer = "counts")

      SO <- Seurat::CreateSeuratObject(counts = chunk_wise_cbind(x = SO_unprocessed),
                                       meta.data = meta.data)
      rm(SO_unprocessed)

      #SO <- merge(x = SO_unprocessed[[1]], y = SO_unprocessed[2:length(SO_unprocessed)], merge.data = T, collapse = T)
      SO <- Seurat::NormalizeData(SO, assay = "RNA", verbose = FindVariableFeatures_args[["verbose"]])

      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = "RNA"),
                                                                     FindVariableFeatures_args))

      if (!is.null(var_feature_set)) {
        Seurat::VariableFeatures(SO) <- var_feature_set
      }

      if (interactive_varfeat_selection && is.null(var_feature_set)) {
        vfplot <- varfeat_plot(SO, n_varfeat = interactive_varfeat_selection_inds)
        print(vfplot)
        FindVariableFeatures_args[["nfeatures"]] <- get_numeric_input("select number of variable features.")
        SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                            assay = "RNA"),
                                                                       FindVariableFeatures_args))
      } else {
        varfeat_plot(obj = SO, n_varfeat = FindVariableFeatures_args[["nfeatures"]])
      }

      if (!is.null(var_feature_filter)) {
        SO <- var_feature_filter_removal(SO = SO,
                                         var_feature_filter = var_feature_filter,
                                         normalization = normalization,
                                         vars.to.regress = vars.to.regress,
                                         SCtransform_args = SCtransform_args,
                                         FindVariableFeatures_args = FindVariableFeatures_args)
      }
    }
    SO <- Seurat::ScaleData(SO, vars.to.regress = vars.to.regress, verbose = FindVariableFeatures_args[["verbose"]])
  }

  SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))

  if (interactive_pc_selection) {
    print(scexpr::elbowplot2(SO, npcs = RunPCA_args[["npcs"]])[["plot"]])

    RunPCA_args[["npcs"]] <- get_numeric_input("select number of pc.")
    SO@reductions[[RunPCA_args[["reduction.name"]]]] <- NULL

    RunPCA_args <- check_RunPCA_args(#obj = SO,
      RunPCA_args = RunPCA_args,
      normalization = normalization,
      npcs = RunPCA_args[["npcs"]],
      seed = RunPCA_args[["seed"]],
      verbose = RunPCA_args[["verbose"]])

    SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))

    RunHarmony_args <- check_RunHarmony_args(RunHarmony_args = RunHarmony_args,
                                             RunPCA_args = RunPCA_args)
  }
  # SO <- Seurat::ProjectDim(SO, reduction = "pca", do.center = T, overwrite = F, verbose = verbose)


  if (batch_corr == "harmony") {
    SO <- Gmisc::fastDoCall(harmony:::RunHarmony.Seurat, args = c(list(object = SO), RunHarmony_args))
    # assay.use = switch(normalization, SCT = "SCT", LogNormalize = "RNA", RNA = "RNA")
  }


  return(list(SO = SO, RunPCA_args = RunPCA_args, RunHarmony_args = RunHarmony_args))
}

var_feature_filter_removal <- function(SO,
                                       var_feature_filter,
                                       max_rounds = 5,
                                       normalization,
                                       vars.to.regress,
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
         length(Seurat::VariableFeatures(SO)[which(!Seurat::VariableFeatures(SO) %in% var_feature_filter)]) < SCtransform_args[["variable.features.n"]] &&
         n <= max_rounds) {
    message(sum(Seurat::VariableFeatures(SO) %in% var_feature_filter), " of var_feature_filter found in Variable Features. Round ", n, " of increasing nhvf to ", SCtransform_args[["variable.features.n"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), " so that var_feature_filter can be removed while nhvf = " , SCtransform_args[["variable.features.n"]], " is met.")

    if (normalization == "SCT") {
      SCtransform_args[which(names(SCtransform_args) == "variable.features.n")] <- SCtransform_args[["variable.features.n"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter))
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 assay = "RNA"),
                                                            SCtransform_args))

      #SO <- Seurat::SCTransform(SO, assay = "RNA", vst.flavor = "v2", method = "glmGamPoi", variable.features.n = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), vars.to.regress = vars.to.regress, seed.use = seed, verbose = verbose)
    }
    if (normalization %in% c("LogNormalize", "RNA")) {
      FindVariableFeatures_args[which(names(FindVariableFeatures_args) == "nfeatures")] <- FindVariableFeatures_args[["nfeatures"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter))
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = "RNA"),
                                                                     FindVariableFeatures_args))
      #SO <- Seurat::FindVariableFeatures(SO, selection.method = "vst", nfeatures = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), verbose = verbose, assay = "RNA")
    }
    n <- n + 1
  }
  Seurat::VariableFeatures(SO) <- Seurat::VariableFeatures(SO)[which(!Seurat::VariableFeatures(SO) %in% var_feature_filter)]

  return(SO)
}

check_RunHarmony_args <- function(RunHarmony_args,
                                  RunPCA_args) {
  RunHarmony_args <- RunHarmony_args[which(!names(RunHarmony_args) %in% c("object", "assay.use", "verbose"))]
  RunHarmony_args[["reduction.use"]] <- RunPCA_args[["reduction.name"]]
  RunHarmony_args[["reduction.save"]] <- paste0("harmony_", RunPCA_args[["reduction.name"]])
  RunHarmony_args[["dims.use"]] <- 1:RunPCA_args[["npcs"]]
  return(RunHarmony_args)
}

check_RunPCA_args <- function(RunPCA_args,
                              normalization,
                              npcs,
                              seed = 42,
                              verbose = T) {
  # actually only very few arguments are allowed to be passed by RunPCA_args. Otherwise the function would break.

  pca_ext <- ifelse(normalization == "LogNormalize", "RNA", normalization)
  pca_ext <- tolower(pca_ext)

  RunPCA_args[["npcs"]] <- npcs
  RunPCA_args[["seed.use"]] <- seed
  RunPCA_args[["verbose"]] <- verbose
  redname <- paste0("pca", npcs, "_", pca_ext)
  RunPCA_args[["reduction.name"]] <- redname

  # i <- 0
  # while (redname %in% names(obj@reductions)) {
  #   i <- i+1
  #   redname <- paste0("pca", npcs, "_", i, "_", pca_ext)
  # }
  #RunPCA_args <- c(list(reduction.name = redname), RunPCA_args)


  # if (!"ndims.print" %in% names(RunPCA_args)) {
  #   RunPCA_args <- c(list(ndims.print = 0), RunPCA_args)
  # }
  # if (!"nfeatures.print" %in% names(RunPCA_args)) {
  #   RunPCA_args <- c(list(nfeatures.print = 0), RunPCA_args)
  # }

  if ("reduction.key" %in% names(RunPCA_args)) {
    message("reduction.key in RunPCA_args not used.")
    RunPCA_args <- RunPCA_args[-which(names(RunPCA_args) == "reduction.key")]
  }
  if ("assay" %in% names(RunPCA_args)) {
    message("assay in RunPCA_args not used.")
    RunPCA_args <- RunPCA_args[-which(names(RunPCA_args) == "assay")]
  }
  return(RunPCA_args)
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

    if (!is.list(celltype_label)) {
      celltype_label <- list(celltype_label)
    }
    if (length(celltype_label) == 1) {
      celltype_label <- rep(celltype_label, length(celltype_refs))
    }
    if (length(celltype_label) != length(celltype_refs)) {
      stop("celltype_label and celltype_refs must have the same lengths. Please choose one label for each ref.")
    }
    ## adjust to seurat, summexp or matrix
    for (i in seq_along(celltype_refs)) {
      if (methods::is(celltype_refs[[i]], "SummarizedExperiment")) {
        celltype_label[[i]] <- intersect(celltype_label[[i]], names(celltype_refs[[i]]@colData@listData))
      } else if (methods::is(celltype_refs[[i]], "Seurat")) {
        celltype_label[[i]] <- intersect(celltype_label[[i]], names(celltype_refs[[i]]@meta.data))
      } else if (methods::is(celltype_refs[[i]], "matrix") || methods::is(celltype_refs[[i]], "sparseMatrix")) {
        celltype_label[[i]] <- celltype_label[[i]][which(lengths(celltype_label[[i]]) == ncol(celltype_refs[[i]]))]
      } else {
        stop("celltype_refs must be Seurat, SummarizedExperiment or gene x cell matrix.")
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
                                             RunHarmony_args,
                                             verbose,
                                             cells,
                                             downsample,
                                             min_cells,
                                             downsample_method) {

  if (methods::is(SO_unprocessed, "list")) {
    SO_unprocessed <- SO_unprocessed
  } else if (methods::is(SO_unprocessed, "character")) {
    if (!any(file.exists(SO_unprocessed))) {
      stop(paste0(SO_unprocessed, "not found."))
    } else {
      if (!any(grepl("\\.rds$", SO_unprocessed, ignore.case = T))) {
        stop("all SO_unprocessed have to be .rds files.")
      }
      SO_unprocessed <- unlist(purrr::map(SO_unprocessed, readRDS))
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

  SO_unprocessed <- subset_SO_unprocessed(
    SO_unprocessed = SO_unprocessed,
    cells = cells,
    downsample = downsample,
    min_cells = min_cells,
    downsample_method = downsample_method)

  # create objects from scratch to rm all previous traces like commands or so
  SO_unprocessed <- parallel::mclapply(SO_unprocessed, function(x) {
    get_layer(obj = x, assay = "RNA", layer = "counts") |>
      SeuratObject::CreateSeuratObject() |>
      SeuratObject::AddMetaData(x@meta.data) |>
      Seurat::NormalizeData(verbose = F, assay = "RNA")
  }, mc.cores = min(1, parallel::detectCores()-4))

  if (batch_corr == "harmony" && length(SO_unprocessed) > 1) {
    if (!"group.by.vars" %in% names(RunHarmony_args)) {
      stop("Please provide one or more group.by.vars from meta.data in RunHarmony_args as a list: ", paste(names(SO_unprocessed[[1]]@meta.data), collapse = ", "), ".")
    }
    ids <- unlist(purrr::map(SO_unprocessed, ~unique(.x@meta.data[[RunHarmony_args[["group.by.vars"]]]])))
    if (anyDuplicated(ids)) {
      message("duplicate harmony group.by.vars found across SO. this may be unexpected.")
    }
    if (length(unique(ids)) == 1) {
      message("only one harmony group.by.vars found. setting to names of SO.")
      SO_unprocessed <- purrr::map(names(SO_unprocessed), function(x) {
        SO_unprocessed[[x]]@meta.data[[RunHarmony_args[["group.by.vars"]]]] <- x
        return(SO_unprocessed[[x]])
      })
    }
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


  # if (length(SO_unprocessed) > 1 && batch_corr == "harmony" && !"group.by.vars" %in% names(RunHarmony_args)) {
  #   stop("Please provide one or more group.by.vars from meta.data in RunHarmony_args as a list: ", paste(names(SO_unprocessed[[1]]@meta.data), collapse = ", "), ".")
  # } else if (length(SO_unprocessed) > 1 && batch_corr == "harmony" && "group.by.vars" %in% names(RunHarmony_args)) {
  #   #if (any(RunHarmony_args[["group.by.vars"]] %in% names()))
  #   # check all SO
  # }

  return(list(SO_unprocessed, samples))
}

subset_SO_unprocessed <- function(SO_unprocessed,
                                  cells,
                                  downsample,
                                  min_cells,
                                  downsample_method) {

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

  if (downsample > 1 && downsample < min_cells) {
    message("downsample set to min_cells.")
    downsample <- min_cells
  } else if (downsample < 1) {
    # check if downsampling forces samples in SO_unprocessed below min_cells
    if (any(sapply(SO_unprocessed, function(x) as.integer(downsample*length(Seurat::Cells(x)))) < min_cells)) {
      message("downsampling leads some samples to drop below min_cells.")
    }
  }

  if (downsample_method == "uniform") {
    if (downsample < 1) {
      SO_unprocessed <- lapply(SO_unprocessed, function(x) subset(
        x,
        cells = sample(Seurat::Cells(x), size = as.integer(downsample*length(Seurat::Cells(x))), replace = FALSE)
      ))
    } else if (downsample > 1) {
      SO_unprocessed <- lapply(SO_unprocessed, function(x) subset(
        x,
        cells = sample(Seurat::Cells(x), size = min(downsample, length(Seurat::Cells(x))), replace = FALSE)
      ))
    }
  }
  if (downsample_method == "leverage") {
    if (downsample < 1) {
      ncells <- purrr::map(SO_unprocessed, ~as.integer(downsample*length(Seurat::Cells(.x))))
    } else if (downsample > 1) {
      ncells <- rep(downsample, length(SO_unprocessed))
    }
    # reduce var feature to avoid error in leverage calc; SketchData caused errors, so do it manually
    # no zero variance columns!
    lvs <- purrr::map(SO_unprocessed, function(x) {
      tryCatch(
        expr = {
          x <- Seurat::FindVariableFeatures(x, nfeatures = 500, verbose = F)
          lvs <- Seurat::LeverageScore(get_layer(
            obj = x,
            assay = "RNA",
            layer = "counts",
            features = Seurat::VariableFeatures(x)
          ), verbose = F)
        },
        error = function(err) {
          x <- Seurat::FindVariableFeatures(x, nfeatures = 200, verbose = F)
          lvs <- Seurat::LeverageScore(get_layer(
            obj = x,
            assay = "RNA",
            layer = "counts",
            features = Seurat::VariableFeatures(x)
          ), verbose = F)
        }
      )
      return(lvs)
    })
    SO_unprocessed <- purrr::pmap(list(SO_unprocessed, ncells, lvs),
                                  function(x, y, z) subset(
                                    x,
                                    cells = sample(Seurat::Cells(x), size = y, replace = F, prob = z)
                                  ))
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
                                   vars.to.regress,
                                   verbose,
                                   seed) {
  # prep SCtransform args once for all
  SCtransform_args <- SCtransform_args[which(!names(SCtransform_args) %in% c("object", "assay", "new.assay.name", "seed.use", "verbose"))]

  SCtransform_args[["variable.features.n"]] <- nhvf
  SCtransform_args[["verbose"]] <- verbose
  SCtransform_args[["seed.use"]] <- seed

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
                                            nhvf,
                                            verbose) {

  # prep FindVariableFeatures args once for all
  FindVariableFeatures_args <- FindVariableFeatures_args[which(!names(FindVariableFeatures_args) %in% c("object", "assay", "verbose"))]
  FindVariableFeatures_args[["nfeatures"]] <- nhvf
  FindVariableFeatures_args[["verbose"]] <- verbose

  return(FindVariableFeatures_args)
}

check_FindIntegrationAnchors_args <- function(FindIntegrationAnchors_args,
                                              anchor.features,
                                              k.filter,
                                              k.score,
                                              npcs) {
  FindIntegrationAnchors_args <- FindIntegrationAnchors_args[which(!names(FindIntegrationAnchors_args) %in% c("object.list", "normalization.method", "verbose"))]
  if (!"anchor.features" %in% names(FindIntegrationAnchors_args)) {
    FindIntegrationAnchors_args <- c(list(anchor.features = anchor.features), FindIntegrationAnchors_args)
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
  # https://github.com/satijalab/seurat/issues/4262
  if (!"dims" %in% names(FindIntegrationAnchors_args)) {
    FindIntegrationAnchors_args <- c(list(dims = 1:npcs), FindIntegrationAnchors_args)
  }
  return(FindIntegrationAnchors_args)
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

  for (i in seq_along(celltype_refs)) {
    for (j in celltype_label[[i]]) {
      out <- labeltransfer_singler(test_obj = SO,
                                   test_clusters = celltype_ref_clusters,
                                   ref_obj = celltype_refs[[i]],
                                   ref_labels = j,
                                   name_prefix = names(celltype_refs)[i])
      SO <- Seurat::AddMetaData(SO, out$metadata)
      Seurat::Misc(SO, paste0(names(celltype_refs)[i], "__", j)) <- out[1:4]
    }
  }
  return(SO)
}


get_numeric_input <- function(prompt = "Enter a number: ") {
  promt <- paste0(prompt, "\n")
  repeat {
    input <- readline(prompt)      # read as string
    num <- suppressWarnings(as.numeric(input))  # try to convert to numeric

    if (!is.na(num)) {
      return(num)                  # return if conversion successful
    } else {
      cat("Invalid input. Please enter a numeric value.\n")
    }
  }
}



chunk_wise_cbind <- function(x, nchunk = 0.2) {
  ## better for memory
  ## and much faster

  if (nchunk < 1) {
    nchunk <- length(x)*nchunk
  }

  x <- make_equal_feature_order(x)

  chunks <- split(x, ceiling(seq_along(x)/nchunk))
  partials <- purrr::map(chunks, ~do.call(cbind, .x))
  result <- do.call(cbind, partials)
  return(result)
}


calc_neighbor_and_cluster <- function(SO,
                                      red = "pca",
                                      npcs = 10,
                                      normalization = "RNA",
                                      FindNeighbors_args = list(),
                                      FindClusters_args = list(),
                                      verbose = TRUE,
                                      mc.cores = 8) {

  # see /Volumes/CMS_SSD_2TB/R_scRNAseq/R_scripts/distance_matrices_and_umap.R
  message("FindNeighbors and FindClusters")
  FindNeighbors_args <- FindNeighbors_args[which(!names(FindNeighbors_args) %in% c("object", "reduction", "verbose"))]
  # if (!"dims" %in% names(FindNeighbors_args)) {
  #   FindNeighbors_args <- c(list(dims = 1:npcs), FindNeighbors_args)
  # }
  fun <- FindNeighbors2 # scexpr:::
  # if (use_nn_from_findneighbors_for_umap) {
  #   fun <- Seurat::FindNeighbors
  # } else {
  #   fun <- scexpr::FindNeighbors2
  # }


  nn_list <- Gmisc::fastDoCall(fun, args = c(list(object = SO@reductions[[red]]@cell.embeddings[,1:npcs],
                                                  verbose = verbose),
                                             FindNeighbors_args))

  #  "_", npcs,
  # names(nn_list[["graphs"]]) <- paste0(gsub("_nn$", "", gsub("_snn$", "", names(nn_list[["graphs"]]))), "_", red, "_", c("nn", "snn"))
  # norm <- ifelse(normalization %in% c("RNA", "LogNormalize"), "RNA", "SCT")
  names(nn_list[["graphs"]]) <- paste0(red, "_", c("nn", "snn")) # norm, "_",


  SO@graphs <- nn_list[["graphs"]]
  SO@misc[["nn.ranked"]] <- nn_list[["nn.ranked"]]

  FindClusters_args <- FindClusters_args[which(!names(FindClusters_args) %in% c("object", "verbose"))]
  if (!"resolution" %in% names(FindClusters_args)) {
    FindClusters_args[["resolution"]] <- 0.8
  }

  cl <- parallel::mclapply(X = FindClusters_args[["resolution"]],
                           snn = SO@graphs[[names(SO@graphs)[2]]],
                           args = FindClusters_args,
                           FUN = function(x, snn, args) {
                             args$resolution <- x
                             Gmisc::fastDoCall(Seurat::FindClusters,
                                               args = c(list(object = snn,
                                                             verbose = F),
                                                        args))
                           }, mc.cores = 8)

  cl <- dplyr::bind_cols(cl)
  cl <- pad_default_cluster_numbers(cl)
  names(cl) <- paste0(names(SO@graphs)[2], "_", names(cl))
  #names(cl) <- unname(sapply(brathering::strsplit2(names(cl), "_"), "[", 2))

  SO <- Seurat::AddMetaData(SO, cl)

  return(SO)
}

pad_default_cluster_numbers <- function(x, minlen = 2) {

  padlen <- max(minlen, max(nchar(as.character(unname(unlist(x))))))
  x <- x |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~.x+1)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~brathering::pad_num_in_str(x = .x, len = padlen)))
  return(x)

}


make_equal_feature_order <-  function(x) {
  # x: list of vectors
  feats <- purrr::map(x, rownames)
  lens <- purrr::map_dbl(feats, length)
  if (length(unique(lens)) > 1) {
    # different features
    common_feats <- purrr::reduce(feats, intersect)
    message("different features across input samples. reducing to common: n = ", length(common_feats))
    message("original: ")
    print(lens)
    message("check new global variable: feature_overlap_df")

    featsdf <- purrr::map(names(feats), function(x) data.frame(x = feats[[x]], y = feats[[x]]))
    featsdf <- purrr::reduce(featsdf, dplyr::full_join, by = "x")
    names(featsdf)[-1] <- names(feats)
    feature_overlap_df <<- featsdf

  } else {
    # same lens but also same order?
    if (!all(purrr::map_lgl(feats, ~all(.x == feats[[1]])))) {
      message("different features order across input samples. will equalize order of features.")
    }
    common_feats <- feats[[1]]
  }

  x <- purrr::map(x, ~.x[common_feats,])
  return(x)
}

pkg_checks <- function() {
  if (!requireNamespace("colrr", quietly = T)) {
    devtools::install_github("Close-your-eyes/colrr")
  }
  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("glmGamPoi", quietly = T)) {
    BiocManager::install("glmGamPoi")
  }
  if (!requireNamespace("Gmisc", quietly = T)) {
    utils::install.packages("Gmisc")
  }
  if (!requireNamespace("harmony", quietly = T)) {
    utils::install.packages("harmony")
  }
  if (!requireNamespace("zap", quietly = T)) {
    devtools::install_github("coolbutuseless/zap")
  }
}

