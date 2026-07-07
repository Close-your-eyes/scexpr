#' Process assay of Seurat object
#'
#' normalize, hvf, pca, batch corr, dim red, clustering
#'
#' @param SO seurat object
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
#' @param seed random seed
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
#' @param var_feature_set set hvf manually
#' @param use_nn_for_umap re-use nearest neighbor graph from cluster calc for umap?
#' @param recalculate force recalculation of SCT assay even when nhvf match.
#' do so, when setting var_feature_set or vars.to.regress
#'
#' @returns seurat object
#' @export
#' @importFrom zeallot %<-%
#'
#' @examples
#' \dontrun{
#' # scan combinations of nhvf and npcs for nice umap and clustering results
#' for (i in seq(500, 4000, 500)) {
#'   print(i)
#'   for (j in seq(10,50,10)) {
#'     print(j)
#'     so <- scexpr::SO_prep04(so,
#'                             reductions = "umap",
#'                             nhvf = i,
#'                             npcs = j,
#'                             SCtransform_args = list(
#'                               vst.flavor = "v2",
#'                               method = "glmGamPoi",
#'                               conserve.memory = F),
#'                             vars.to.regress = "pct_mt")
#'   }
#' }
#' # find a clusterings tha has as few splits as possible compared to the other
#' clustres <- purrr::map(purrr::set_names(so@misc$clusterings), function(x) {
#'   brathering::compare_labels(so@meta.data[[x]], so@meta.data$renalcellgroup)
#' })
#' # find best matching clusterings with simpson concentration index for each row
#' rsimp <- purrr::map(clustres, ~rowSums(.x[["x_props"]][["mat"]]^2))
#' df <- data.frame(value = unlist(rsimp, use.names = F), name = rep(names(rsimp), lengths(rsimp)))
#'
#' ggplot(df, aes(x = reorder(name, value), y = value)) +
#'   geom_boxplot() +
#'   geom_point()
#'
#' dfsumm <- dplyr::summarise(df, meanval = median(value), minval = min(value), .by = name)
#' #' }
SO_prep04 <- function(SO,
                      reductions = c("umap"),
                      nhvf = 800,
                      npcs = 20,
                      normalization = c("SCT", "RNA", "LogNormalize"),
                      batch_corr = c("harmony", "none"),
                      vars.to.regress = NULL,
                      seed = 42,
                      #var_feature_filter = NULL,
                      var_feature_set = NULL,
                      verbose = T,
                      FindVariableFeatures_args = list(),
                      SCtransform_args = list(
                        vst.flavor = "v2",
                        method = "glmGamPoi",
                        conserve.memory = T,
                        residual_type = "pearson"),
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
                      use_nn_for_umap = F,
                      ...) {

  options(future.globals.maxSize = 20 * 1024^3)

  reductions <- match.arg(tolower(reductions), c("umap", "tsne"), several.ok = T)
  normalization <- rlang::arg_match(normalization)
  batch_corr <- rlang::arg_match(batch_corr)

  RunPCA_args <- scexpr:::check_RunPCA_args(
    RunPCA_args = RunPCA_args,
    normalization = normalization,
    npcs = npcs,
    nhvf = nhvf,
    seed = seed,
    verbose = verbose)

  SCtransform_args <- scexpr:::check_SCtransform_args(
    SCtransform_args = SCtransform_args,
    nhvf = nhvf,
    vars.to.regress = vars.to.regress,
    verbose = verbose,
    seed = seed)

  FindVariableFeatures_args <- scexpr:::check_FindVariableFeatures_args(
    FindVariableFeatures_args = FindVariableFeatures_args,
    nhvf = nhvf,
    verbose = verbose
  )

  c(RunHarmony_args, batch_corr) %<-% scexpr:::check_RunHarmony_args(RunHarmony_args = RunHarmony_args,
                                                                     RunPCA_args = RunPCA_args,
                                                                     batch_corr = batch_corr,
                                                                     SO = SO)



  c(SO,
    RunPCA_args,
    RunHarmony_args) %<-% make_so_simple(SO = SO,
                                         normalization = normalization,
                                         SCtransform_args = SCtransform_args,
                                         verbose = verbose,
                                         var_feature_filter = var_feature_filter,
                                         vars.to.regress = vars.to.regress,
                                         FindVariableFeatures_args = FindVariableFeatures_args,
                                         RunPCA_args = RunPCA_args,
                                         RunHarmony_args = RunHarmony_args,
                                         batch_corr = batch_corr,
                                         interactive_varfeat_selection = interactive_varfeat_selection,
                                         interactive_varfeat_selection_inds = interactive_varfeat_selection_inds,
                                         interactive_pc_selection = interactive_pc_selection,
                                         recalculate = recalculate,
                                         var_feature_set = var_feature_set)

  red <- switch(
    batch_corr,
    harmony = RunHarmony_args[["reduction.save"]],
    integration = RunPCA_args[["reduction.name"]],
    none = RunPCA_args[["reduction.name"]]
  )
  names_wo_clust <- names(SO@meta.data)

  SO <- calc_neighbor_and_cluster(obj = SO,
                                  red = red,
                                  npcs = RunPCA_args[["npcs"]],
                                  FindNeighbors_args = FindNeighbors_args,
                                  FindClusters_args = FindClusters_args,
                                  verbose = verbose,
                                  mc.cores = 10)

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

  # add clusterings and their colors to Misc
  names_w_clust <- names(SO@meta.data)
  SeuratObject::Misc(SO, "clusterings") <- c(SeuratObject::Misc(SO, "clusterings"),
                                             setdiff(names_w_clust, names_wo_clust))

  for (i in c(SeuratObject::Misc(SO, "clusterings"), "orig.ident")) {
    SO <- add_group_color_to_misc(obj = SO, meta_col = i)
  }

  try(expr = {
    # quick to calculate
    # save disk space
    SO@assays[["RNA"]]@layers[["scale.data"]] <- NULL
  }, silent = T)


  return(SO)
}


make_so_simple <- function(SO,
                           normalization,
                           SCtransform_args,
                           verbose,
                           var_feature_filter,
                           vars.to.regress,
                           FindVariableFeatures_args,
                           RunPCA_args,
                           RunHarmony_args,
                           batch_corr,
                           interactive_varfeat_selection,
                           interactive_varfeat_selection_inds,
                           interactive_pc_selection,
                           recalculate,
                           var_feature_set) {

  if (normalization == "SCT") {
    if (!normalization %in% names(SO@assays) ||
        (normalization %in% names(SO@assays) && SCtransform_args[["variable.features.n"]] != length(SO@assays[["SCT"]]@var.features)) ||
        recalculate) {

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
      # if (!is.null(var_feature_filter)) {
      #   SO <- .var_feature_filter_removal2(SO = SO,
      #                                      assay = assay,
      #                                      new.assay.name = new.assay.name,
      #                                      var_feature_filter = var_feature_filter,
      #                                      normalization = normalization,
      #                                      nhvf = nhvf,
      #                                      vars.to.regress = vars.to.regress,
      #                                      seed = seed,
      #                                      verbose = verbose,
      #                                      SCtransform_args = SCtransform_args,
      #                                      FindVariableFeatures_args = FindVariableFeatures_args)
      # }

    }

  } else if (normalization %in% c("LogNormalize", "RNA")) {

    ## always recalculate for now

    SO <- Seurat::NormalizeData(SO, verbose = verbose)
    SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO),
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

    # if (!is.null(var_feature_filter)) {
    #   SO <- .var_feature_filter_removal2(SO = SO,
    #                                      assay = assay,
    #                                      var_feature_filter = var_feature_filter,
    #                                      normalization = normalization,
    #                                      nhvf = nhvf,
    #                                      vars.to.regress = vars.to.regress,
    #                                      seed = seed,
    #                                      verbose = verbose,
    #                                      SCtransform_args = SCtransform_args,
    #                                      FindVariableFeatures_args = FindVariableFeatures_args)
    # }

    SO <- Seurat::ScaleData(SO, vars.to.regress = vars.to.regress, verbose = verbose)
    normalization <- "RNA"
  }

  SeuratObject::DefaultAssay(SO) <- normalization
  message("DefaultAssay: ", SeuratObject::DefaultAssay(SO))

  if (!RunPCA_args[["reduction.name"]] %in% names(SO@reductions) || recalculate) {
    SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO),
                                                     RunPCA_args))
    if (interactive_pc_selection) {
      print(scexpr::elbowplot2(SO, npcs = RunPCA_args[["npcs"]])[["plot"]])
      RunPCA_args[["npcs"]] <- get_numeric_input("select number of pc.")

      RunPCA_args <- check_RunPCA_args(
        RunPCA_args = RunPCA_args,
        normalization = normalization,
        npcs = RunPCA_args[["npcs"]],
        seed = RunPCA_args[["seed"]],
        nhvf = SCtransform_args[["variable.features.n"]],
        verbose = RunPCA_args[["verbose"]])

      SO <- Gmisc::fastDoCall(Seurat::RunPCA, args = c(list(object = SO), RunPCA_args))

      c(RunHarmony_args, batch_corr) %<-% check_RunHarmony_args(RunHarmony_args = RunHarmony_args,
                                                                RunPCA_args = RunPCA_args,
                                                                batch_corr = batch_corr)
    }
    # SO <- Seurat::ProjectDim(SO, reduction = reduction.name, do.center = T, overwrite = F, verbose = verbose)
  }

  if (batch_corr == "harmony") {
    if (!RunHarmony_args[["reduction.save"]] %in% names(SO@reductions) || recalculate) {
      SO <- Gmisc::fastDoCall(harmony:::RunHarmony.Seurat, args = c(list(object = SO), RunHarmony_args))
    }
  }


  return(list(SO = SO, RunPCA_args = RunPCA_args, RunHarmony_args = RunHarmony_args))
}


.var_feature_filter_removal2 <- function(SO,
                                         assay,
                                         new.assay.name = NULL,
                                         var_feature_filter,
                                         max_rounds = 5,
                                         normalization,
                                         nhvf,
                                         vars.to.regress,
                                         seed,
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
                                                                 new.assay.name = new.assay.name),
                                                            SCtransform_args))

      #SO <- Seurat::SCTransform(SO, assay = "RNA", vst.flavor = "v2", method = "glmGamPoi", variable.features.n = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), vars.to.regress = vars.to.regress, seed.use = seed, verbose = verbose)
    }
    if (normalization == "LogNormalize") {
      FindVariableFeatures_args[which(names(FindVariableFeatures_args) == "nfeatures")] <- FindVariableFeatures_args[["nfeatures"]] + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter))
      SO <- Gmisc::fastDoCall(Seurat::FindVariableFeatures, args = c(list(object = SO,
                                                                          assay = assay),
                                                                     FindVariableFeatures_args))
      #SO <- Seurat::FindVariableFeatures(SO, selection.method = "vst", nfeatures = nhvf + length(intersect(Seurat::VariableFeatures(SO), var_feature_filter)), verbose = verbose, assay = "RNA")
    }
    n <- n + 1
  }
  Seurat::VariableFeatures(SO) <- Seurat::VariableFeatures(SO)[which(!Seurat::VariableFeatures(SO) %in% var_feature_filter)]

  return(SO)
}

