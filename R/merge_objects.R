#' Merge list of Seurat objects
#'
#' Currently counts of RNA assay and meta data only.
#'
#' @param obj_list list of Seurat objects.
#' @param obj_ident_col new column name to identify origin of cells by object
#'
#' @returns merged seurat
#' @export
#'
#' @examples
#' \dontrun{
#' so <- merge_objects(obj_list = list(so1,so2,so3))
#' }
merge_objects <- function(obj_list,
                          obj_ident_col = "obj_ident") {
  ## currently focused on RNA assay
  # no messages or warnings
  featlist <- purrr::map(obj_list, ~rownames(get_layer(.x, layer = "counts", assay = "RNA")))
  common_feat <- purrr::reduce(featlist, intersect)
  obj_merge <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(obj_list, get_layer, layer = "counts", features = common_feat)),
                                          meta.data = dplyr::bind_rows(purrr::map(obj_list, ~.x@meta.data)))
  obj_merge <- Seurat::NormalizeData(obj_merge)
  if (!is.null(names(obj_list))) {
    obj_names <- names(obj_list)
  } else {
    obj_names <- as.character(seq_along(obj_list))
  }
  obj_merge@meta.data[[obj_ident_col]] <- rep(obj_names, purrr::map_int(obj_list, ~nrow(.x@meta.data)))
  return(obj_merge)
}
