#' Get layer data with a few fallback mechanisms
#'
#' Works with v4 and v5 objects. Join layers in v5 objects if possible
#' and necessary. Returns NULL if layer is 0x0 matrix.
#'
#' @param obj Seurat object
#' @param assay assay
#' @param layer layer or slot
#' @param features filter for features
#' @param cells filter for cells
#' @param as what to return, sparse or dense matrix or data frame
#' @param transpose apply Matrix::t() ?
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
                      as = c("sparse", "dense", "df"),
                      transpose = F) {

  as <- rlang::arg_match(as)
  assay <- rlang::arg_match(assay, names(obj@assays))

  if (utils::compareVersion(as.character(obj@version), "4.9.9") == 1) {

    x <- tryCatch(expr = {
      get_layer_v5(obj = obj, assay = assay, layer = layer)
    }, error = function(err) {
      # fallback to v4 method, e.g. for integrated assay which looks in v5 object
      # like assay from v4
      return(get_layer_v4(obj = obj, assay = assay, layer = layer))
    })

  } else {

    x <- get_layer_v4(obj = obj, assay = assay, layer = layer)

  }

  if (all(dim(x) == c(0,0))) {
    message("get_layer: layer ", layer, " is empty.")
    return(NULL)
  }

  Seurat::DefaultAssay(obj) <- assay
  rownames(x) <- rownames(x) %||% rownames(obj)
  colnames(x) <- colnames(x) %||% colnames(obj)
  features <- features %||% rownames(obj)
  cells <- cells %||% colnames(obj)

  data <- x[features, cells, drop = F]

  if (transpose) {
    data <- Matrix::t(data)
  }
  if (as == "dense") {
    data <- as.matrix(data)
  }
  if (as == "df") {
    data <- data.frame(as.matrix(data), check.names = F)
  }

  return(data)
}

get_layer_v4 <- function(obj, assay, layer) {
  layer <- rlang::arg_match(layer, setdiff(slotNames(obj@assays[[assay]]),
                                           c("assay.orig", "var.featuress", "meta.featuress", "misc", "key")))
  x <- slot(obj@assays[[assay]], layer)
  return(x)
}

get_layer_v5 <- function(obj, assay, layer) {
  if (!layer %in% names(obj@assays[[assay]]@layers)) {
    # fallback 1: check for unjoined layers
    newlayer <- gsub("\\.[[:digit:]]{1,}$", "", names(obj@assays[[assay]]@layers))
    if (sum(newlayer == layer) > 1) {
      message("get_layer: joining layers to pull data.")
      obj <- SeuratObject::JoinLayers(obj, assay = assay)
    } else {
      # get the error
      layer <- rlang::arg_match(layer, names(obj@assays[[assay]]@layers))
    }
  } else {
    # get the error
    layer <- rlang::arg_match(layer, names(obj@assays[[assay]]@layers))
  }

  x <- obj@assays[[assay]]@layers[[layer]]
  return(x)
}
