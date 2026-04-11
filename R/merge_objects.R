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

  if (!is.list(obj_list)) {
    stop("obj_list must be a list of seurat objects.")
  }

  if (is.null(names(obj_list))) {
    nme <- as.character(deparse(substitute(obj_list)))
    names(obj_list) <- strsplit(gsub("list\\(|\\)", "", nme), ", ")[[1]]
  }

  featlist <- purrr::map(obj_list, ~rownames(get_layer(.x, layer = "counts", assay = "RNA")))
  common_feat <- purrr::reduce(featlist, intersect)

  shared_feat_freq <- purrr::map_dbl(featlist, ~length(intersect(.x, common_feat)))/length(common_feat)
  message("shared feature frequencies:")
  print(shared_feat_freq)

  obj_merge <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(obj_list, get_layer, assay = "RNA", layer = "counts", features = common_feat)),
                                          meta.data = purrr::map_dfr(obj_list, ~.x@meta.data))
  obj_merge <- Seurat::NormalizeData(obj_merge)
  if (!is.null(names(obj_list))) {
    obj_names <- names(obj_list)
  } else {
    obj_names <- as.character(seq_along(obj_list))
  }
  obj_merge@meta.data[[obj_ident_col]] <- rep(obj_names, purrr::map_int(obj_list, ~nrow(.x@meta.data)))
  return(obj_merge)
}
