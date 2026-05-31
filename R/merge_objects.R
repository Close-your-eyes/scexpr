#' Merge list of Seurat objects
#'
#' Currently counts of RNA assay and meta data only.
#'
#' @param obj_list list of Seurat objects.
#' @param obj_ident_col new column name to identify origin of cells by object
#' @param merge_reductions try to merge reductions?
#' @param misc_from_first use misc slot from first obj in list?
#'
#' @returns merged seurat
#' @export
#'
#' @examples
#' \dontrun{
#' so <- merge_objects(obj_list = list(so1,so2,so3))
#' }
merge_objects <- function(obj_list,
                          obj_ident_col = "obj_ident",
                          merge_reductions = F,
                          misc_from_first = F) {

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

  if (merge_reductions) {
    ## no check for equal columns of reduction
    # no check for duplicate cellname
    redlst <- purrr::map(obj_list, ~names(.x@reductions))
    if (any(is.null(redlst))) {
      message("some obj w/o reductions.")
    } else {
      redcommon <- purrr::reduce(redlst, intersect)
      if (!length(redcommon)) {
        message("no common reductions")
      } else {
        for (x in redcommon) {
          obj_merge@reductions[[x]] <- Seurat::CreateDimReducObject(
            embeddings = do.call(rbind, purrr::map(obj_list, function(y) y@reductions[[x]]@cell.embeddings)),
            assay = obj_list[[1]]@reductions[[x]]@assay.used,
            key = obj_list[[1]]@reductions[[x]]@key
          )
        }
      }
    }
  }

  if (misc_from_first) {
    for (x in names(obj_list[[1]]@misc)) {
      obj_merge@misc[[x]] <- obj_list[[1]]@misc[[x]]
    }
  }

  if (!is.null(names(obj_list))) {
    obj_names <- names(obj_list)
  } else {
    obj_names <- as.character(seq_along(obj_list))
  }
  obj_merge@meta.data[[obj_ident_col]] <- rep(obj_names, purrr::map_int(obj_list, ~nrow(.x@meta.data)))
  return(obj_merge)
}
