#' Run linear discriminant analysis (LDA) for optimal separation of groups
#'
#' Alternative to PCA. E.g. run clustering first, then LDA, then tsne or umap on
#' LDA. Or separate celltypes as labelled by SingleR in optimal way.
#' What does Seurat::RunLDA do ?
#'
#' @param obj seurat object
#' @param groups column with grouping in obj@meta.data
#' @param assay which assay to pull data from
#' @param layer which layer
#' @param overwrite overwrite existing dim red with sme name?
#' @param features which genes; less of them inreases calculation speed
#' @param reduction.name name of reduction to be created
#' @param reduction.key key of reduction to be created, needs trailing underscore
#'
#' @returns
#' @export
#'
#' @examples
run_lda <- function(obj,
                    groups,
                    assay = "RNA",
                    layer = "data",
                    overwrite = F,
                    features = Seurat::VariableFeatures(obj),
                    reduction.name = "LDA",
                    reduction.key = paste0(reduction.name, "_")) {

  if (missing(groups) || !groups %in% names(obj@meta.data)) {
    stop("groups is missing or not found in obj@meta.data.")
  }

  if (reduction.name %in% names(obj@reductions) && !overwrite) {
    stop(reduction.name, " exists in obj@reductions. set overwrite = T to replace.")
  }

  if (!grepl("_$", reduction.key)) {
    reduction.key <- paste0(reduction.key, "_")
    message("added underscore to reduction.key.")
  }

  lda_model <- MASS::lda(x = Matrix::t(get_layer(
    obj,
    assay = assay,
    layer = layer,
    features = features
  )),
  grouping = obj@meta.data[[groups]])
  lda_coords <- stats::predict(lda_model)
  obj@reductions[[reduction.name]] <- SeuratObject::CreateDimReducObject(
    embeddings = lda_coords$x,
    assay = assay,
    key = reduction.key
  )
  obj@misc[[reduction.name]] <- lda_model

  return(obj)
}
