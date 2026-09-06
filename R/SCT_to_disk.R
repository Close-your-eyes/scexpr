#' Save an SCT assay to disk
#'
#' Saves the `"SCT"` assay from a Seurat object as an RDS file, removes the
#' assay, and returns the modified object. Before removing the SCT assay, the
#' default assay is changed to the first remaining assay.
#'
#' If `obj@misc$object_name` exists, the SCT filename is constructed by adding
#' the suffix `"_SCT"` to the object name using [brathering::suffix_add()].
#' The modified object, without its SCT assay, is also saved using
#' `obj@misc$object_name` as its filename.
#'
#' If `obj@misc$object_name` does not exist, the SCT assay is saved using a
#' timestamp-based filename of the form `"SCT_<timestamp>.rds"`, and the
#' modified object is not saved.
#'
#' If the object has no SCT assay, it is returned unchanged and no files are
#' written.
#'
#' @param obj A Seurat object. When an SCT assay is present, the object must
#'   contain at least one other assay to use as the default assay.
#' @param save_dir Character string specifying an existing directory in which
#' @param update_on_disk update object on disk?
#' @param compression gz compression rate
#'   to save the RDS file or files. Defaults to the current working directory.
#' @return The input Seurat object with the `"SCT"` assay removed and its
#'   default assay set to the first remaining assay. If no SCT assay is
#'   present, the original object is returned unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- SCT_to_disk(obj, save_dir = "saved_assays")
#' }
SCT_to_disk <- function(obj,
                        save_dir = getwd(),
                        update_on_disk = T,
                        compression = 1L) {

  .ensure_package("brathering")

  if (!"SCT" %in% names(obj@assays)) {
    message("no SCT assay.")
    return(obj)
  }

  other_assays <- setdiff(names(obj@assays), "SCT")
  Seurat::DefaultAssay(obj) <- other_assays[1]

  message(save_dir)

  savefile <- paste0("SCT_", gsub("[^0-9]", "", Sys.time()), ".rds")
  if ("object_name" %in% names(obj@misc)) {
    savefile <- brathering::suffix_add(obj@misc[["object_name"]], suffix = "_SCT")
  }


  message("writing SCT assay to disk.")
  brathering::saverds2(x = obj@assays[["SCT"]],
                       file = file.path(save_dir, savefile),
                       compression = compression)
  obj@assays[["SCT"]] <- NULL

  if (update_on_disk) {
    if ("object_name" %in% names(obj@misc)) {
      message("writing obj without SCT assay to disk.")
      brathering::saverds2(x = obj,
                           file = file.path(save_dir, obj@misc[["object_name"]]),
                           compression = compression)
    } else {
      message("object_name not found in misc-slot. NOT saving object w/o SCT to disk.")
    }
  }

  return(obj)
}
