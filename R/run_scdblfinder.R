#' Add artificial doublets to Seurat object
#'
#' Using scDblFinder::addDoublets with accessory lines of code.
#'
#' @param obj seurat object
#' @param samples vector of cluster idents, e.g. orig.ident from obj;
#' obj is spit by that and doublets added separately; can be one sample only
#' @param clusters vector of cluster idents, e.g. a meta.col from obj
#' @param prep_new_obj make a new seurat to fit doublets into dim red
#' @param ... args to SO_prep02
#'
#' @returns seurat
#'
#' @examples
add_doublets <- function(obj,
                         samples = obj@meta.data$orig.ident,
                         clusters,
                         prep_new_obj = T,
                         ...) {

  if (!methods::is(obj, "Seurat")) {
    stop("obj should be a Seurat.")
  }
  if (missing(clusters)) {
    stop("please provide clusters.")
  }
  if (is.null(clusters)) {
    stop("clusters not found.")
  }
  if (suppressWarnings(!anyNA(as.numeric(clusters)))) {
    message("numric cluster names: in case of sample.int error, avoid numeric cluster names.")
  }
  if (is.null(samples)) {
    stop("samples not found.")
  }
  if (length(clusters) != length(samples)) {
    stop("clusters and samples must be equal in length. Provide the actual meta columns. Not just their names.")
  }

  if (is.matrix(obj)) {
    mat <- obj
  } else {
    mat <- get_layer(obj, layer = "counts")
  }

  mat <- brathering::split_mat(mat, f = samples, byrow = F)
  clusters <- split(clusters, f = samples)

  sodbl <- purrr::map2(purrr::set_names(names(mat)), clusters, function(x,y) {

    cds <- scDblFinder::addDoublets(x = SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat[[x]])), clusters = as.character(y))
    #cds <- scDblFinder::addDoublets(x = mat[[x]], clusters = as.character(y))
    # if (is.nullcds@colData@listData$type) {
    #   return(NULL)
    # }
    cds@colData@listData$type <- ifelse(cds@colData@listData$type == "singlet", "single_transcriptome_as_is", "doublet_simul")

    counts <- SummarizedExperiment::assays(cds)[["counts"]]
    dbl_ind <- which(grepl("^doublet", colnames(counts)))
    colnames(counts)[dbl_ind] <- paste0(colnames(counts)[dbl_ind], "_", x)

    so <- Seurat::CreateSeuratObject(counts, project = x)
    so@meta.data$dbl_type <- cds@colData@listData[["type"]]
    so@meta.data$dbl_cluster <- cds@colData@listData[["cluster"]]

    if (methods::is(obj, "Seurat")) {
      so <- Seurat::AddMetaData(so, obj@meta.data[,setdiff(names(obj@meta.data), names(so@meta.data))])
      # replace NA with unique (enrich doublet meta data)
      so@meta.data <- so@meta.data |>
        dplyr::mutate(dplyr::across(dplyr::everything(), ~{
          if (dplyr::n_distinct(.x, na.rm = TRUE) == 1) {
            tidyr::replace_na(.x, unique(stats::na.omit(.x)))
          } else {
            .x
          }
        }))
    }
    print(table(cds@colData@listData$type))
    return(so)
  })


  # new obj to have added doublets appear in dim reduction
  # other option would be to just run pca and umap and clustering but thats what SO_prep02 does actually
  if (prep_new_obj) {
    nhvf <- 2000
    if (!is.null(Seurat::VariableFeatures(obj))) {
      nhvf <- length(Seurat::VariableFeatures(obj))
    }
    npcs <- 20
    if (!is.null(obj@reductions$pca)) {
      npcs <- ncol(obj@reductions$pca@cell.embeddings)
    }


    sodbl <- SO_prep02(SO_unprocessed = sodbl,
                       normalization = ifelse("SCT" %in% names(obj@assays), "SCT", "LogNormalize"),
                       interactive_varfeat_selection = F,
                       interactive_pc_selection = F,
                       nhvf = nhvf,
                       npcs = npcs,
                       SCtransform_args = list(
                         vst.flavor = "v2",
                         method = "glmGamPoi",
                         conserve.memory = T #just in case
                       ))
  } else {
    sodbl <- scexpr::merge_objects(sodbl)
  }

  return(sodbl)
}



#' Run scDblFinder with accessory code lines
#'
#' @param obj seurat object
#' @param assay which assay
#' @param features gene features to use
#' @param samples samples to treat separately, vector of idents
#' @param clusters cluster dents to help the algorithm and optionally to
#' add_doublets
#' @param min_umi minimum total UMI sum across features per cell; cells below
#' are not considered for doublet inference
#' @param dims arg to scDblFinder::scDblFinder
#' @param add_doublets add artificial doublets
#' @param add_doublets_tell tell scDblFinder::scDblFinder about artificial doublets
#'
#' @returns seurat
#' @export
#'
#' @examples
run_scdblfinder <- function(obj,
                            assay = "RNA",
                            features = scexpr::get_var_features(so),
                            samples_col = "orig.ident",
                            clusters_col = "seurat_clusters",
                            min_umi = 300,
                            min_cells_per_group = 20,
                            dims = ncol(obj@reductions$pca@cell.embeddings),
                            add_doublets = F,
                            add_doublets_tell = F) {

  if (any(!c(samples_col, clusters_col) %in% names(obj@meta.data))) {
    stop("samples_col or clusters_col not found.")
  }
  samples <- obj@meta.data[[samples_col]]
  clusters <- obj@meta.data[[clusters_col]]

  knownDoublets <- NULL
  if (add_doublets) {
    obj <- add_doublets(obj = obj,
                        samples = samples,
                        clusters = clusters,
                        prep_new_obj = T)
    knownDoublets <- obj@meta.data$dbl_type == "doublet_simul"
    if (!add_doublets_tell) {
      knownDoublets <- NULL
    }
  }
  if (!length(features)) {
    message("no hvf found. try to rescue.")
    Seurat::DefaultAssay(obj) <- "RNA"
    obj <- Seurat::FindVariableFeatures(obj)
    features <- Seurat::VariableFeatures(obj)
    if (!length(features)) {
      stop("no hvf.")
    }
  }

  lowumi <- which(Matrix::colSums(get_layer(obj = obj,
                                            assay = assay,
                                            features = features,
                                            layer = "counts")) < min_umi)

  if (length(lowumi)) {
    message("number of cells below min_umi: ", length(lowumi))
  }

  highumi <- setdiff(Seurat::Cells(obj), names(lowumi))
  Seurat::DefaultAssay(obj) <- assay

  nn <- obj@meta.data |>
    dplyr::summarise(n = dplyr::n(), .by = c(clusters_col, samples_col)) |>
    dplyr::filter(n < min_cells_per_group)
  cells_keep <- obj@meta.data |>
    dplyr::filter(!(!!rlang::sym(samples_col) %in% nn[[samples_col]]))

  obj2 <- subset(obj, cells = rownames(cells_keep))
  obj2 <- subset(obj2, cells = highumi)

  ## subset is defects
  scf <- scDblFinder::scDblFinder(sce = get_layer(obj = obj2,
                                                  assay = assay,
                                                  layer = "counts"),
                                  samples = obj2@meta.data[[samples_col]],
                                  nfeatures = features,
                                  dims = dims,
                                  knownDoublets = knownDoublets,
                                  knownUse = "positive")
  scfdf <- as.data.frame(scf@colData)
  scfdf <- scfdf[,-1]
  names(scfdf) <- gsub("scDblFinder\\.", "dbl_", names(scfdf))

  scfdf$dbl_class <- ifelse(scfdf$dbl_class == "singlet", "singlet_infer", "doublet_calc")
  obj <- Seurat::AddMetaData(obj, scfdf)

  return(obj)
}
