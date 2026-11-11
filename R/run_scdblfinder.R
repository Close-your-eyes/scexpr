#' Detect doublets in a Seurat object with scDblFinder
#'
#' Runs scDblFinder on assay counts, treating samples separately, and adds
#' doublet scores and classifications to the object's cell metadata.
#' Optionally adds artificial doublets before inference.
#'
#' @param obj A Seurat object containing a counts layer and the metadata
#'   columns named by \code{samples_col} and \code{clusters_col}.
#' @param assay Name of the assay whose counts are used for doublet inference
#'   and the minimum UMI filter. Defaults to \code{"RNA"}.
#' @param features Character vector of feature names used for doublet
#'   inference and for calculating the UMI total per cell. Passed to
#'   \code{scDblFinder::scDblFinder()} as \code{nfeatures}. If empty, the
#'   function attempts to find variable features in the RNA assay.
#' @param samples_col Name of the cell metadata column identifying samples
#'   to process separately, typically independent cell captures. Defaults to
#'   \code{"orig.ident"}. Supply a column name, not a vector of sample labels.
#' @param clusters_col Name of the cell metadata column containing cluster
#'   labels. Used to identify small sample-cluster groups and, when
#'   \code{add_doublets = TRUE}, to generate artificial doublets. These labels
#'   are not passed to the final \code{scDblFinder::scDblFinder()} call.
#' @param min_umi Minimum total counts across \code{features} required for
#'   a cell to enter doublet inference. Cells below this threshold remain in
#'   the returned object but are excluded from inference. Defaults to 300.
#' @param min_cells_per_group Minimum number of cells in an observed
#'   sample-cluster group. If any such group is smaller than this threshold,
#'   the entire corresponding sample is excluded from inference. Group sizes
#'   are calculated before the minimum UMI filter. Defaults to 20.
#' @param dims Number of dimensions passed to
#'   \code{scDblFinder::scDblFinder()}. Defaults to the number of components
#'   in \code{obj}'s PCA cell embeddings; provide this argument explicitly
#'   when the input has no PCA reduction.
#' @param add_doublets Logical; whether to add artificial doublets separately
#'   within each sample using \code{clusters_col}, then rerun preprocessing
#'   with \code{SO_prep02()}. Defaults to \code{FALSE}.
#' @param add_doublets_tell Logical; whether to pass the artificial-doublet
#'   indicators to scDblFinder as \code{knownDoublets}, with
#'   \code{knownUse = "positive"}. Only used when
#'   \code{add_doublets = TRUE}. Defaults to \code{FALSE}.
#'
#' @details
#' Both metadata columns must exist, even when artificial doublets are not
#' requested. Samples excluded by the group-size filter and cells excluded
#' by the UMI filter remain in the returned object; inference is performed
#' only on the retained subset.
#'
#' Artificial doublets are generated from RNA counts before the
#' assay-specific inference filters are applied. Samples with no expected
#' artificial doublets retain their original cells.
#'
#' @returns A Seurat object with scDblFinder metadata added for cells used
#'   in inference. The \code{scDblFinder.} prefix is replaced by
#'   \code{dbl_}, including \code{dbl_score} and \code{dbl_class}.
#'   Class labels are \code{"singlet_infer"} and \code{"doublet_calc"}.
#'   Cells excluded from inference receive no new inference results.
#'   When \code{add_doublets = TRUE}, the returned object also includes
#'   any artificial doublets and the \code{dbl_type} and \code{dbl_cluster}
#'   metadata created by the helper.
#'
#' @seealso \code{scDblFinder::scDblFinder()},
#'   \code{scDblFinder::addDoublets()}
#' @export
#'
#' @examples
#' \dontrun{
#' # Use variable features from the input object and process each sample
#' # separately. The input must contain both metadata columns.
#' so <- run_scdblfinder(
#'   obj = so,
#'   features = scexpr::get_var_features(so),
#'   samples_col = "orig.ident",
#'   clusters_col = "seurat_clusters",
#'   dims = 20
#' )
#'
#' table(so$dbl_class, useNA = "ifany")
#' }
run_scdblfinder <- function(obj,
                            assay = "RNA",
                            features = NULL,
                            samples_col = "orig.ident",
                            clusters_col = "seurat_clusters",
                            min_umi = 300,
                            min_cells_per_group = 20,
                            dims = ncol(obj@reductions$pca@cell.embeddings),
                            add_doublets = F,
                            add_doublets_tell = F) {

  scexpr:::.ensure_packages(c("scDblFinder", "brathering"))


  if (any(!c(samples_col, clusters_col) %in% names(obj@meta.data))) {
    stop("samples_col or clusters_col not found.")
  }
  samples <- obj@meta.data[[samples_col]]
  clusters <- obj@meta.data[[clusters_col]]

  if (is.null(features)) {
    features <- scexpr::get_var_features(obj)
  }

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

  obj2 <- subset2(obj, cells = rownames(cells_keep))
  obj2 <- subset2(obj2, cells = highumi)

  scf <- scDblFinder::scDblFinder(sce = get_layer(obj = obj2,
                                                  assay = assay,
                                                  layer = "counts"),
                                  samples = obj2@meta.data[[samples_col]],
                                  clusters = obj2@meta.data[[clusters_col]],
                                  nfeatures = features,
                                  dims = dims,
                                  knownDoublets = knownDoublets,
                                  knownUse = "positive")

  scfdf <- as.data.frame(scf@colData)
  scfdf <- scfdf[,-1]
  names(scfdf) <- gsub("scDblFinder\\.", "dbl_", names(scfdf))

  # scfdf$dbl_class <- ifelse(scfdf$dbl_class == "singlet", "singlet_infer", "doublet_calc")
  obj <- Seurat::AddMetaData(obj, scfdf)

  return(obj)
}


#' Add artificial doublets to Seurat object
#'
#' Used within run_scdblfinder function.
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

  scexpr:::.ensure_package("scDblFinder")

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
    cds <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = mat[[x]])
    )
    cds$type <- rep("singlet", ncol(cds))
    cds$cluster <- as.character(y)

    cds <- scDblFinder::addDoublets(
      x = cds,
      clusters = cds$cluster
    )

    cds$type <- ifelse(
      cds$type == "singlet",
      "single_transcriptome_as_is",
      "doublet_simul"
    )

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
