#' Get layer data with a few fallback mechanisms
#'
#' Works with v4 and v5 objects. For v5 assays, an exact layer name takes
#' precedence. If it is absent, use layers whose names start with the
#' requested name followed by a period, joining them when needed.
#' Returns `NULL` for a zero-by-zero matrix.
#'
#' @param obj Seurat object
#' @param assay assay
#' @param layer layer or slot
#' @param features filter for features
#' @param cells filter for cells
#' @param as what to return, sparse or dense matrix or data frame
#' @param transpose apply Matrix::t() ?
#' @param features_invert invert feature selection. exclude the ones provided.
#' @param cells_invert invert cell selection. exclude the ones provided.
#'
#' @returns matrix (sparse/dense) or data frame
#' @export
#'
#' @importFrom rlang %||%
#'
#' @examples
#' \dontrun{
#' get_layer(SO, assay = "SCT", layer = "counts", features = c("CD4", "CD8A"),
#'           as = "dense")
#' }
get_layer <- function(obj,
                      assay = "RNA",
                      layer = "data",
                      features = NULL,
                      cells = NULL,
                      features_invert = F,
                      cells_invert = F,
                      as = c("sparse", "dense", "df"),
                      transpose = F) {

  # split = NULL
  as <- rlang::arg_match(as)
  assay <- rlang::arg_match(assay, names(obj@assays))



  if (!is.null(features) && !length(features)) {
    # empty gene_features argument in get_data, e.g.
    mat <- matrix(nrow = 0, ncol = ncol(obj)) # n cells
    if (transpose) {
      mat <- t(mat)
    }
    if (as == "df") {
      return(data.frame(mat))
    }
    if (as == "dense") {
      return(mat)
    }
    if (as == "sparse") {
      requireNamespace("Matrix")
      return(methods::as(mat, "sparseMatrix"))
    }
  }

  assay_obj <- obj[[assay]]
  if (inherits(assay_obj, "StdAssay")) {
    x <- get_layer_v5(assay_obj, layer)
  } else {
    x <- get_layer_v4(assay_obj, layer)
  }

  if (is.null(x) || all(dim(x) == c(0, 0))) {
    message("get_layer: layer ", layer, " is empty.")
    return(NULL)
  }

  features <- features %||% rownames(x)
  cells <- cells %||% colnames(x)

  if (!all(cells %in% colnames(x))) {
    message("get_layer: some cells not found.")
    cells <- cells[which(cells %in% colnames(x))]
  }
  if (!all(features %in% rownames(x))) {
    message("get_layer: some features not found.")
    features <- features[which(features %in% rownames(x))]
  }

  if (features_invert) {
    features <- setdiff(rownames(x), features)
  }

  if (cells_invert) {
    cells <- setdiff(colnames(x), cells)
  }

  if (!identical(features, rownames(x)) || !identical(cells, colnames(x))) {
    x <- x[features, cells, drop = FALSE]
  }

  if (transpose) {
    x <- Matrix::t(x)
  }
  if (as == "dense") {
    x <- as.matrix(x)
  }
  if (as == "df") {
    x <- data.frame(as.matrix(x), check.names = F)
  }

  return(x)
}

get_layer_v4 <- function(assay_obj, layer) {
  layer <- rlang::arg_match(layer, setdiff(methods::slotNames(assay_obj),
                                           c("assay.orig", "var.featuress", "meta.featuress", "misc", "key")))
  methods::slot(assay_obj, layer)
}

resolve_assay_layers <- function(assay_obj, layer) {
  available <- SeuratObject::Layers(assay_obj)
  if (layer %in% available) {
    return(layer)
  }
  split_layers <- available[startsWith(available, paste0(layer, "."))]
  if (length(split_layers)) {
    return(split_layers)
  }
  rlang::arg_match(layer, available)
}

get_layer_v5 <- function(assay_obj, layer) {
  layers <- resolve_assay_layers(assay_obj, layer)
  if (length(layers) > 1L) {
    message("get_layer: joining layers to pull data.")
    pattern <- paste0("^", stringr::str_escape(layer), "\\.")
    assay_obj <- SeuratObject::JoinLayers(assay_obj, layers = pattern, new = layer)
    layers <- layer
  }
  SeuratObject::LayerData(assay_obj, layer = layers)
}

get_layer_cells <- function(assay_obj, layer) {
  if (inherits(assay_obj, "StdAssay")) {
    layers <- resolve_assay_layers(assay_obj, layer)
    return(SeuratObject::Cells(assay_obj, layer = layers))
  }
  colnames(get_layer_v4(assay_obj, layer))
}
