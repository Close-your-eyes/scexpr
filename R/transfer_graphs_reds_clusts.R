#' Transfer reductions, graphs, and metadata between Seurat objects
#'
#' Copies PCA-based reductions, associated graphs, and metadata columns
#' from one Seurat object (`obj2`) to another (`obj1`). If naming conflicts
#' occur in `obj1`, new elements are automatically renamed to avoid overwriting.
#'
#' Both objects must contain identical cells (same names and order),
#' otherwise the function will stop with an error.
#'
#' @param obj1 A \code{Seurat} object that will receive reductions, graphs,
#'   and metadata.
#' @param obj2 A \code{Seurat} object providing reductions, graphs,
#'   and metadata to transfer.
#'
#' @return A \code{Seurat} object (`obj1`) with additional reductions,
#'   graphs, and metadata copied from `obj2`.
#'
#' @details
#' The function searches for reductions in `obj2` with names matching
#' the pattern `"^pca[0-9]+_"` (e.g., \code{pca1_}, \code{pca2_}).
#' Corresponding reductions, graphs, and metadata entries are copied
#' into `obj1`. If a name already exists in `obj1`, a numeric suffix
#' is appended to ensure uniqueness.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # run SO_prep02 with different parameters. then copy
#' # graphs, reduction, clustering from one to another.
#' obj <- readRDS("")
#' deg <- compare_cluster_deg_to_hvf(obj, meta_col = "celltype")
#' deg1 <- compare_cluster_deg_to_hvf(obj, meta_col = "celltype", topn_per_group = 5)
#' raw <- readRDS("")
#' hinze3 <- SO_prep02(raw,
#'                     cells = Seurat::Cells(hinze),
#'                     npcs = 12,
#'                     interactive_pc_selection = F,
#'                     var_feature_set = deg$marker$feature)
#' hinze4 <- SO_prep02(raw,
#'                     cells = Seurat::Cells(hinze),
#'                     npcs = 12,
#'                     interactive_pc_selection = F,
#'                     var_feature_set = deg1$marker$feature)
#' hinze5 <- SO_prep02(raw,
#'                     cells = Seurat::Cells(hinze),
#'                     npcs = 12,
#'                     interactive_pc_selection = F,
#'                     var_feature_set = deg$marker$feature,
#'                     use_nn_for_umap = T)
#' hinze6 <- SO_prep02(raw,
#'                     cells = Seurat::Cells(hinze),
#'                     npcs = 12,
#'                     interactive_pc_selection = F,
#'                     var_feature_set = deg1$marker$feature,
#'                     use_nn_for_umap = T)
#'
#' obj3 <- Seurat::AddMetaData(obj3, metadata = obj@meta.data[,"celltype",drop = F])
#' obj4 <- Seurat::AddMetaData(obj4, metadata = obj@meta.data[,"celltype",drop = F])
#' obj5 <- Seurat::AddMetaData(obj5, metadata = obj@meta.data[,"celltype",drop = F])
#' obj6 <- Seurat::AddMetaData(obj6, metadata = obj@meta.data[,"celltype",drop = F])
#' feature_plot2(list(obj3, obj4, obj5, obj6), "celltype", facet_scales = "free")
#' obj3 <- transfer_graphs_reds_clusts(obj3, obj4)
#' }
transfer_graphs_reds_clusts <- function(obj1, obj2) {

  ## copy elements from obj2 to obj1
  cells1 <- SeuratObject::Cells(obj1)
  cells2 <- SeuratObject::Cells(obj2)

  if (!identical(cells1, cells2)) {
    stop("cells not identical.")
  }

  base_red <- grep("^pca[[:digit:]]{1,}_", names(obj2@reductions), value = T)

  if (!length(base_red)) {
    message("not matching reduction found in obj2.")
  }

  for (i in base_red) {
    i1 <- i
    j <- 1
    while (i %in% names(obj1@reductions)) {
      j <- j + 1
      i <- paste0(i, j)
    }
    names(obj2@reductions) <- gsub(i1, i , names(obj2@reductions))
    names(obj2@graphs) <- gsub(i1, i , names(obj2@graphs))
    names(obj2@meta.data) <- gsub(i1, i , names(obj2@meta.data))

    # copy
    for (z in grep(i, names(obj2@reductions), value = T)) {
      obj1@reductions[[z]] <- obj2@reductions[[z]]
    }
    for (z in grep(i, names(obj2@graphs), value = T)) {
      obj1@graphs[[z]] <- obj2@graphs[[z]]
    }
    for (z in grep(i, names(obj2@meta.data), value = T)) {
      obj1@meta.data[[z]] <- obj2@meta.data[[z]]
    }
  }

  return(obj1)
}
