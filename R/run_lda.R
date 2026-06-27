#' Run linear discriminant analysis on a Seurat object
#'
#' Performs linear discriminant analysis (LDA) to find axes that optimally
#' separate predefined groups of cells. This can be useful as a supervised
#' alternative to PCA, for example after clustering cells or assigning labels
#' with tools such as SingleR. The resulting LDA coordinates are stored as a
#' Seurat dimensional reduction, and the fitted LDA model is stored in
#' `obj@misc`.
#'
#' @param obj A Seurat object.
#' @param groups Character scalar. Name of the column in `obj@meta.data`
#'   containing group labels used for LDA.
#' @param assay Character scalar. Assay from which expression data should be
#'   extracted. Default is `"RNA"`.
#' @param layer Character scalar. Layer to use from the selected assay.
#'   Default is `"data"`.
#' @param overwrite Logical. Whether to overwrite an existing dimensional
#'   reduction with the same `reduction.name`. Default is `FALSE`.
#' @param features Character vector of features to use for LDA. Fewer features
#'   can substantially reduce runtime. Defaults to
#'   `Seurat::VariableFeatures(obj)`.
#' @param reduction.name Character scalar. Name of the dimensional reduction to
#'   create. Default is `"LDA"`.
#' @param reduction.key Character scalar. Key prefix for the dimensional
#'   reduction. A trailing underscore is added automatically if missing.
#'
#' @return A Seurat object with an added dimensional reduction named
#'   `reduction.name`. The fitted `MASS::lda` model is stored in
#'   `obj@misc[[reduction.name]]`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- run_lda(
#'   obj,
#'   groups = "seurat_clusters",
#'   assay = "RNA",
#'   layer = "data",
#'   reduction.name = "LDA"
#' )
#' }
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
