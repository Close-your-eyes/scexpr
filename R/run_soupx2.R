#' Run default SoupX pipeline
#'
#' @param filt filtered gene x cell matrix or list thereof; when list: everything
#' is treated as one
#' @param raw raw gene x cell matrix or list thereof
#' @param obj provide seurat object instead of filt
#' @param so_prep_args args to SO_prep03/SO_prep02
#' @param autoEstCont_args args to SoupX::autoEstCont
#' @param adjustCounts_method to method in SoupX::adjustCounts; one of
#' c("subtraction", "soupOnly", "multinomial")
#'
#' @returns Seurat object with SOUPX corrected counts
#' @export
#'
#' @examples
#'\dontrun{
#' filt <- Seurat::Read10X_h5("filtered_feature_bc_matrix.h5")
#' raw <- Seurat::Read10X_h5("raw_feature_bc_matrix.h5")
#' so <- run_soupx2(filt, raw)
#' }
run_soupx2 <- function(filt,
                       raw,
                       obj = NULL,
                       so_prep_args = list(FindClusters_args = list(resolution = 0.8),
                                           nhvf = 2000,
                                           npcs = 20),
                       autoEstCont_args = list(),
                       adjustCounts_method = "subtraction") {

  # c("subtraction", "soupOnly", "multinomial")
  adjustCounts_method <- rlang::arg_match(adjustCounts_method, multiple = T)

  if (is.null(obj)) {
    obj <- Gmisc::fastDoCall(what = SO_prep03,
                             args = c(list(matrix_list = filt),
                                      so_prep_args)) ## get cluster and dr; filt list handled internally
  }


  if (is.list(raw)) {
    raw <- purrr::reduce(raw, cbind)
  }
  sc <- SoupX::SoupChannel(tod = raw, toc = get_layer(obj, assay = "RNA", layer = "counts"))
  sc = SoupX::setClusters(sc, stats::setNames(obj@meta.data[[obj@misc[["clusterings"]][length(obj@misc[["clusterings"]])]]], colnames(obj)))
  # last reduction should be umap or tsne; first are pca and maybe harmony
  sc = SoupX::setDR(sc, obj@reductions[[length(obj@reductions)]]@cell.embeddings[colnames(sc$toc),])

  # level of background contamination (represented as rho)
  sc <- Gmisc::fastDoCall(what = SoupX::autoEstCont,
                          args = c(list(sc = sc), autoEstCont_args))

  # run one or multiple methods
  if (length(adjustCounts_method) == 1) {
    sx_counts <- SoupX::adjustCounts(sc = sc, verbose = 0, method = adjustCounts_method)
    obj[["SOUPX"]] <- SeuratObject::CreateAssay5Object(sx_counts)
    obj <- add_soup_pct_meta(obj)
  } else {
    ## naming by method
    for (i in adjustCounts_method) {
      sx_counts <- SoupX::adjustCounts(sc = sc, verbose = 0, method = i)
      assay_adj <- paste0("SOUPX_", i)
      obj[[assay_adj]] <- SeuratObject::CreateAssay5Object(sx_counts)
      obj <- add_soup_pct_meta(obj, assay_adj = assay_adj,
                               name = paste0("pct_soup_SoupX_", i))
    }
  }
  return(list(so = obj, sc = sc))
}

add_soup_pct_meta <- function(so,
                              assay_orig = "RNA",
                              assay_adj = "SOUPX",
                              name = "pct_soup_SoupX") {

  cs1 <- Matrix::colSums(get_layer(obj = so, assay = assay_orig, layer = "counts"))
  cs2 <- Matrix::colSums(get_layer(obj = so, assay = assay_adj, layer = "counts"))
  pct_soup = (cs1 - cs2)/cs1*100
  so <- SeuratObject::AddMetaData(object = so, metadata = pct_soup, col.name = name)
  return(so)
}


