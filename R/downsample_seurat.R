#' Downsample a Seurat object
#'
#' Make a compact function with one subset command only because subset still
#' leads to errors for unknown reason in some subsequent functions
#' with seurat object, sometimes.
#'
#' @param obj A Seurat object.
#' @param downsample Either:
#'   - a proportion in (0, 1], or
#'   - the number of cells to retain (> 1).
#' @param group Optional metadata column used for stratified sampling.
#' @param nhvf downsample to n variable features (RNA assay), var features
#'  are calculated before any downsampling or cell removal.
#' @param recreate_RNA_assay sometimes LogMap in `seurat@assays[["RNA"]]@cells`
#'  is not subsetted. Force this by recreating RNA assay.
#' @param features features retain when subsetting, joined with nhvf if provided
#' @param cells_in cells to keep (before downsampling)
#' @param cells_ex cells to exclude (before downsampling)
#'
#' @returns A downsampled Seurat object.
#'
#' @export
#'
#' @examples
#'\dontrun{
#' so02 <- downsample_seurat(so2, downsample = 300,
#'                           recreate_RNA_assay = T,
#'                           nhvf = 200, features = "CD3E",
#'                           cells_in = cells2(subset(so2, !pca14_rna900_snn_res_0.1 %in% c("06", "07"))))
#' }
downsample_seurat <- function(obj,
                              downsample = 1,
                              group = NULL,
                              nhvf = NULL,
                              features = NULL,
                              cells_in = NULL,
                              cells_ex = NULL,
                              recreate_RNA_assay = F) {

  if (downsample == 1) {
    return(obj)
  }
  stopifnot(inherits(obj, "Seurat"))

  if (!is.null(group)) {
    group <- rlang::sym(group)
  }

  if (!is.null(nhvf)) {
    # what if SCT is default assay?
    obj <- Seurat::FindVariableFeatures(obj, nfeature = nhvf)
    nhvf <- Seurat::VariableFeatures(obj)
  }

  meta <- obj@meta.data

  if (!is.null(cells_in)) {
    meta <- meta[cells_in,]
  }
  if (!is.null(cells_ex)) {
    meta <- meta[!which(rownames(cells) %in% cells_ex),]
  }

  if (downsample < 1) {
    cells <- rownames(dplyr::slice_sample(meta, prop = downsample, by = group))
  } else if (downsample > 1) {
    cells <- rownames(dplyr::slice_sample(meta, n = downsample, by = group))
  }

  obj <- subset(obj, cells = cells, features = unique(c(nhvf, features)))

  if (recreate_RNA_assay) {
    assaymeta <- obj@assays[["RNA"]]@meta.data
    layers <- names(obj@assays[["RNA"]]@layers)
    if ("counts" %in%  layers && "data" %in% layers) {
      obj@assays[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = get_layer(obj, layer = "counts"),
                                                              data = get_layer(obj))
    } else if ("counts" %in% layers) {
      obj@assays[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = get_layer(obj, layer = "counts"),
                                                              data = NULL)
    } else if ("data" %in% layers) {
      obj@assays[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = NULL,
                                                              data = get_layer(obj))
    } else {
      message("RNA assay not recreated.")
    }
    obj@assays[["RNA"]]@meta.data <- assaymeta
  }

  return(obj)
}
