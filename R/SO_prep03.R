#' Prepare Seurat from gene x cell matrix with quick defaults
#'
#' SO_prep02 is called.
#'
#' @param matrix_list list of gene x cell matrices, e.g. from read_10X_data
#' @param normalization which normalization; LogNormalize is quicker
#' @param nhvf n var features
#' @param npcs n principle components
#' @param interactive_varfeat_selection see SO_prep02
#' @param interactive_pc_selection see SO_prep02
#' @param ... args to SO_prep02
#'
#' @returns seurat object
#' @export
#'
#' @examples
SO_prep03 <- function(matrix_list,
                      normalization = "LogNormalize",
                      nhvf = 800,
                      npcs = 20,
                      interactive_varfeat_selection = F,
                      interactive_pc_selection = F,
                      ...) {

  if (!is.list(matrix_list)) {
    matrix_list <- list(matrix_list)
  }
  if (is.null(names(matrix_list))) {
    names(matrix_list) <- as.character(seq_along(matrix_list))
  }

  so_unpr <- purrr::map(matrix_list, Seurat::CreateSeuratObject)
  obj <- SO_prep02(
    SO_unprocessed = so_unpr,
    normalization = normalization,
    nhvf = nhvf,
    npcs = npcs,
    interactive_varfeat_selection = interactive_varfeat_selection,
    interactive_pc_selection = interactive_pc_selection,
    ...
  )
  return(obj)

  # clusters <- Seurat::CreateSeuratObject(matrix_list[[1]]) |>
  #   Seurat::NormalizeData(verbose = F) |>
  #   Seurat::FindVariableFeatures(selection.method = "vst", nfeatures = nhvf, verbose = F, assay = "RNA") |>
  #   Seurat::ScaleData(verbose = F) |>
  #   Seurat::RunPCA(verbose = F) |>
  #   Seurat::RunUMAP(verbose = F) |>
  #   Seurat::FindNeighbors(verbose = F) |>
  #   Seurat::FindClusters(verbose = F, algorithm = 1, resolution = resolution_SoupX)
}
