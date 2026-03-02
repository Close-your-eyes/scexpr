#' Process assay of Seurat object
#'
#' normalize, hvf, pca, batch corr, dim red, clustering
#'
#' @param SO seurat object
#' @param assay which assay to process
#' @param reductions which reduction to calculate, tsne and/or umap
#' @param nhvf number high variables features, passed to Seurat::SCTransform
#' or Seurat::FindVariableFeatures or Seurat::SelectIntegrationFeatures
#' @param npcs number of principle components in calculate in pca and to
#' consider for downstream functions like tSNE or UMAP
#' @param normalization algorithm for normalization of UMIs
#' @param batch_corr procedure for batch correction
#' @param vars.to.regress passed to Seurat::SCTransform or Seurat::ScaleData; will be applied
#' independent of what is set for batch_corr; if batch_corr is set to 'none', then only vars.to.regress is
#' used to regress out a variable in meta.data which may be sample-specific; other then that
#' it may not be meaningful to regress a sample-specific variable and perform batch_corr;
#' rather vars.to.regress may be used in combination with batch_corr to regress percent_mito or so
#' @param seeed random seed
#' @param var_feature_filter character vector of features which are to exclude from variable features;
#' this will affect downstream PCA, dimension reduction and clustering;
#' e.g.: var_feature_filter = grep("^TR[ABGD]V", rownames(SOqc_split[[1]]), value = T) to
#' exclude T cell receptor gene segments
#' @param verbose print messages and progress bars from functions
#' @param FindVariableFeatures_args arguments to Seurat::FindVariableFeatures
#' @param SCtransform_args arguments to Seurat::SCTransform but c("object", "assay", "new.assay.name", "seed.use", "verbose")
#' @param RunPCA_args arguments to Seurat::RunPCA
#' @param RunUMAP_args arguments to Seurat::RunUMAP
#' @param RunTSNE_args arguments to scexpr::run_fft_tsne
#' @param FindNeighbors_args arguments to Seurat::FindNeighbors
#' @param FindClusters_args arguments to Seurat::FindClusters
#' @param RunHarmony_args arguments to harmony::RunHarmony
#' @param interactive_varfeat_selection get asked?
#' @param interactive_varfeat_selection_inds get asked?
#' @param interactive_pc_selection do conduct interactive PC selection?
#' @param ... not used
#'
#' @returns seurat object
#' @export
#'
#' @examples
SO_prep04 <- function(SO,
                      assay = "RNA",
                      reductions = c("umap"),
                      nhvf = 800,
                      npcs = 20,
                      normalization = c("SCT", "LogNormalize"),
                      batch_corr = c("harmony", "none"),
                      vars.to.regress = NULL,
                      seeed = 42,
                      var_feature_filter = NULL,
                      verbose = F,
                      FindVariableFeatures_args = list(),
                      SCtransform_args = list(
                        vst.flavor = "v2",
                        method = "glmGamPoi",
                        conserve.memory = F
                      ),
                      RunPCA_args = list(),
                      RunUMAP_args = list(),
                      RunTSNE_args = list(theta = 0.01),
                      FindNeighbors_args = list(dims = 1:npcs),
                      FindClusters_args = list(resolution = seq(0.1,0.8,0.1)),
                      RunHarmony_args = list(group.by.vars = "orig.ident"),
                      interactive_varfeat_selection = F,
                      interactive_varfeat_selection_inds = seq(max(nhvf/10,50),
                                                               min(3*nhvf, nrow(SO)),
                                                               length.out = 11),
                      interactive_pc_selection = F,
                      recalculate = F,
                      ...) {

  reductions <- match.arg(tolower(reductions), c("umap", "tsne"), several.ok = T)
  normalization <- rlang::arg_match(normalization)
  batch_corr <- rlang::arg_match(batch_corr)

  if (!assay %in% names(SO@assays)) {
    stop("assay not found.")
  }

  check_RunPCA_args(RunPCA_args)

  SCtransform_args <- check_SCtransform_args(SCtransform_args,
                                             nhvf,
                                             vars.to.regress)
  FindVariableFeatures_args <- check_FindVariableFeatures_args(FindVariableFeatures_args,
                                                               nhvf)

  if (!"group.by.vars" %in% names(RunHarmony_args)) {
    batch_corr <- "none"
  } else if (batch_corr == "harmony") {
    if (any(!RunHarmony_args[["group.by.vars"]] %in% names(SO@meta.data))) {
      stop("some group.by.vars from RunHarmony_args not in names(SO@meta.data).")
    }
    if (nrow(unique(SO@meta.data[,RunHarmony_args[["group.by.vars"]],drop = F])) == 1) {
      message("only one level in group.by.vars: batch_corr set to none.")
      batch_corr <- "none"
    }
  }




  SO <- make_so_simple(SO = SO,
                       assay = assay,
                       normalization = normalization,
                       SCtransform_args = SCtransform_args,
                       seeed = seeed,
                       verbose = verbose,
                       var_feature_filter = var_feature_filter,
                       nhvf = nhvf,
                       vars.to.regress = vars.to.regress,
                       FindVariableFeatures_args = FindVariableFeatures_args,
                       RunPCA_args = RunPCA_args,
                       npcs = npcs,
                       RunHarmony_args = RunHarmony_args,
                       batch_corr = batch_corr,
                       interactive_varfeat_selection = interactive_varfeat_selection,
                       interactive_varfeat_selection_inds = interactive_varfeat_selection_inds,
                       interactive_pc_selection = interactive_pc_selection,
                       recalculate = recalculate)


  nhvf <- length(Seurat::VariableFeatures(SO))
  red <- SO@misc$latest_reduction



  if (any(grepl("umap", reductions, ignore.case = T))) {
    RunUMAP_args <- RunUMAP_args[which(!names(RunUMAP_args) %in% c("object", "seed.use", "reduction", "verbose"))]
    if (!"dims" %in% names(RunUMAP_args)) {
      RunUMAP_args <- c(list(dims = 1:npcs), RunUMAP_args)
    }

    tryCatch(expr = {
      reduction.name <- paste0("umap_", red)
      if (!reduction.name %in% names(SO@reductions) || recalculate) {
        SO <- Gmisc::fastDoCall(Seurat::RunUMAP, args = c(list(object = SO,
                                                               reduction = red,
                                                               reduction.name = reduction.name,
                                                               seed.use = seeed,
                                                               verbose = verbose),
                                                          RunUMAP_args))
      }
      #SO <- rename_reduction(SO, reduction = "umap", new_name = paste0("umap_", red))
    }, error = function(err) {
      message("umap failed.")
    })
    #SO <- Seurat::RunUMAP(object = SO, umap.method = "uwot", metric = "cosine", dims = 1:npcs, seed.use = seeed, reduction = red, verbose = verbose, ...)
  }

  if (any(grepl("tsne", reductions, ignore.case = T))) {
    RunTSNE_args <- RunTSNE_args[which(!names(RunTSNE_args) %in% c("object", "seed.use", "reduction", "verbose"))]

    tryCatch(expr = {
      reduction.name <- paste0("tsne_", red)
      if (!reduction.name %in% names(SO@reductions) || recalculate) {
        SO <- Gmisc::fastDoCall(scexpr::run_fft_tsne, args = c(list(SO = SO,
                                                                    reduction = red,
                                                                    reduction.name = reduction.name,
                                                                    rand_seed = seeed),
                                                               RunTSNE_args))
      }
    }, error = function(err) {
      message("fallback to barnes-hut tsne; default seurat method.")
      if (!"dims" %in% names(RunTSNE_args)) {
        RunTSNE_args <- c(list(dims = 1:npcs), RunTSNE_args)
      }
      if (!"num_threads" %in% names(RunTSNE_args)) {
        RunTSNE_args <- c(list(num_threads = 0), RunTSNE_args)
      }
      #RunTSNE_args[["tsne.method"]] <- "FIt-SNE"
      SO <- Gmisc::fastDoCall(Seurat::RunTSNE, args = c(list(object = SO,
                                                             reduction = red,
                                                             reduction.name = paste0("tsne_", red),
                                                             seed.use = seeed,
                                                             verbose = verbose),
                                                        RunTSNE_args))
    })

    #SO <- Seurat::RunTSNE(object = SO, dims = 1:npcs, seed.use = seeed, reduction = red, verbose = verbose, num_threads = 0, ...)
  }


  FindNeighbors_args <- FindNeighbors_args[which(!names(FindNeighbors_args) %in% c("object", "reduction", "verbose"))]
  if (!"dims" %in% names(FindNeighbors_args)) {
    FindNeighbors_args <- c(list(dims = 1:npcs), FindNeighbors_args)
  }
#  "_", npcs,
  graph_names <- paste0(SO@reductions[[red]]@assay.used, "_", ifelse(grepl("pca", red), "pca", "harmony"), "_", c("nn", "snn"))
  SO <- Gmisc::fastDoCall(Seurat::FindNeighbors, args = c(list(object = SO,
                                                               reduction = red,
                                                               graph.name = graph_names,
                                                               verbose = verbose),
                                                          FindNeighbors_args))

  FindClusters_args <- FindClusters_args[which(!names(FindClusters_args) %in% c("object", "verbose"))]
  SO <- Gmisc::fastDoCall(Seurat::FindClusters, args = c(list(object = SO,
                                                              graph.name = graph_names[2],
                                                              verbose = verbose),
                                                         FindClusters_args))


  try(expr = {
    # quick to calculate
    # save disk space
    SO@assays[["RNA"]]@layers[["scale.data"]] <- matrix(NA)
  }, silent = T)
  try(expr = {
    # quick to calculate
    # save disk space
    SO@assays[[assay]]@layers[["scale.data"]] <- NULL
  }, silent = T)


  return(SO)
}


make_so_simple <- function(SO,
                           assay,
                           normalization,
                           SCtransform_args,
                           seeed,
                           verbose,
                           var_feature_filter,
                           nhvf,
                           vars.to.regress,
                           FindVariableFeatures_args,
                           RunPCA_args,
                           npcs,
                           RunHarmony_args,
                           batch_corr,
                           interactive_varfeat_selection,
                           interactive_varfeat_selection_inds,
                           interactive_pc_selection,
                           recalculate) {

  if (normalization == "SCT") {
    new.assay.name = ifelse(assay == "RNA", "SCT", paste0("SCT_", assay))
    assay_out <- new.assay.name

    if (!new.assay.name %in% names(SO@assays) || recalculate) {
      SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                 assay = assay,
                                                                 new.assay.name = new.assay.name,
                                                                 seed.use = seeed,
                                                                 verbose = verbose),
                                                            SCtransform_args))

      if (interactive_varfeat_selection) {
        vfplot <- varfeat_plot(SO, n_varfeat = interactive_varfeat_selection_inds)
        print(vfplot)
        new_nhvf <- get_numeric_input("select number of variable features.")
        nhvf <- new_nhvf
        SCtransform_args[["variable.features.n"]] <- nhvf
        SO <- Gmisc::fastDoCall(Seurat::SCTransform, args = c(list(object = SO,
                                                                   assay = assay,
                                                                   new.assay.name = new.assay.name,
                                                                   seed.use = seeed,
                                                                   verbose = verbose),
                                                              SCtransform_args))
      } else {
        varfeat_plot(SO, n_varfeat = nhvf)
      }

      # remove var features which are to filter
      if (!is.null(var_feature_filter)) {
        SO <- .var_feature_filter_removal2(SO = SO,
                                           assay = assay,
                                           new.assay.name = new.assay.name,
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

  } else if (normalization == "LogNormalize") {

    ## always recalculate for now

    SO <- Seurat::NormalizeData(SO, assay = assay, verbose = verbose)
    SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                        assay = assay,
                                                                        verbose = verbose),
                                                                   FindVariableFeatures_args))

    if (interactive_varfeat_selection) {
      vfplot <- varfeat_plot(SO, n_varfeat = interactive_varfeat_selection_inds)
      print(vfplot)
      new_nhvf <- get_numeric_input("select number of variable features.")
      nhvf <- new_nhvf
      FindVariableFeatures_args[["nfeatures"]] <- nhvf
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = assay,
                                                                          verbose = verbose),
                                                                     FindVariableFeatures_args))
    } else {
      varfeat_plot(SO, n_varfeat = nhvf)
    }

    if (!is.null(var_feature_filter)) {
      SO <- .var_feature_filter_removal2(SO = SO,
                                         assay = assay,
                                         var_feature_filter = var_feature_filter,
                                         normalization = normalization,
                                         nhvf = nhvf,
                                         vars.to.regress = vars.to.regress,
                                         seeed = seeed,
                                         verbose = verbose,
                                         SCtransform_args = SCtransform_args,
                                         FindVariableFeatures_args = FindVariableFeatures_args)
    }
    SO <- Seurat::ScaleData(SO, vars.to.regress = vars.to.regress, verbose = verbose)
    assay_out <- assay
  }

  SeuratObject::DefaultAssay(SO) <- assay_out
  message("DefaultAssay: ", SeuratObject::DefaultAssay(SO))

  reduction.name <- ifelse(assay %in% c("RNA", "SCT"), "pca", paste0("pca_", assay_out))
  if (!reduction.name %in% names(SO@reductions) || npcs != ncol(SO@reductions[[reduction.name]]@cell.embeddings) || recalculate) {
    SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO,
                                                          npcs = npcs,
                                                          reduction.name = reduction.name,
                                                          seed.use = seeed,
                                                          verbose = verbose),
                                                     RunPCA_args))
    if (interactive_pc_selection) {
      print(Seurat::ElbowPlot(SO, ndims = npcs))
      npcs <- get_numeric_input("select number of pc.")
      assign("npcs", npcs, envir = parent.frame())
      SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO,
                                                            npcs = npcs,
                                                            reduction.name = reduction.name,
                                                            seed.use = seeed,
                                                            verbose = verbose),
                                                       RunPCA_args))
    }
    SO <- Seurat::ProjectDim(SO, reduction = reduction.name, do.center = T, overwrite = F, verbose = verbose)
  }


  if (batch_corr == "harmony") {
    #reduction.save <- ifelse(reduction.name == "pca", "harmony", paste0("harmony_", reduction.name))
    reduction.save <- gsub("pca", "harmony", reduction.name)
    RunHarmony_args <- RunHarmony_args[which(!names(RunHarmony_args) %in% c("object", "assay.use", "verbose"))]
    if (!reduction.save %in% names(SO@reductions) || npcs != ncol(SO@reductions[[reduction.save]]@cell.embeddings) || recalculate) {
      SO <- Gmisc::fastDoCall(harmony::RunHarmony, args = c(list(object = SO),
                                                            reduction.use = reduction.name,
                                                            reduction.save = reduction.save,
                                                            RunHarmony_args))
    }
    SeuratObject::Misc(SO, "latest_reduction") <- reduction.save
  } else {
    SeuratObject::Misc(SO, "latest_reduction") <- reduction.name
  }


  return(SO)
}


.var_feature_filter_removal2 <- function(SO,
                                         assay,
                                         new.assay.name = NULL,
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
                                                                 assay = assay,
                                                                 new.assay.name = new.assay.name,
                                                                 seed.use = seeed,
                                                                 verbose = verbose),
                                                            SCtransform_args))

      #SO <- Seurat::SCTransform(SO, assay = "RNA", vst.flavor = "v2", method = "glmGamPoi", variable.features.n = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), vars.to.regress = vars.to.regress, seed.use = seeed, verbose = verbose)
    }
    if (normalization == "LogNormalize") {
      FindVariableFeatures_args[which(names(FindVariableFeatures_args) == "nfeatures")] <- FindVariableFeatures_args[["nfeatures"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter))
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = assay,
                                                                          verbose = verbose),
                                                                     FindVariableFeatures_args))
      #SO <- Seurat::FindVariableFeatures(SO, selection.method = "vst", nfeatures = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), verbose = verbose, assay = "RNA")
    }
    n <- n + 1
  }
  Seurat::VariableFeatures(SO) <- Seurat::VariableFeatures(SO)[which(!Seurat::VariableFeatures(SO) %in% var_feature_filter)]

  return(SO)
}

