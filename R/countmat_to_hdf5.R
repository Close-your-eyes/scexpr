#' Save UMI counts from matrix to disk in initial 10X format
#'
#'
#' @param mat complete UMI count matrix
#' @param colsplit optionally: vector of factor levels to split mat by
#' (column wise), becomes separate folders in out_dir
#' @param out_dir output directory
#' @param sub_dir sub directory in out_dir or in folders of colsplit
#' @param types format(s) to save
#' @param ... arguments to DropletUtils::write10xCounts
#'
#' @return nothing in R but files on disk
#' @export
#'
#' @examples
#' \dontrun{
#' # obtain count_matrix from Seurat object or any other single cell format
#' # tcell_sub_split is a named list of data frames with one row per cell
#' # count_matrix is a combined mat of all cells in tcell_sub_split
#' countmat_to_hdf5(mat = count_matrix,
#'                  colsplit = rep(names(tcell_sub_split), purrr::map_int(tcell_sub_split, nrow)),
#'                  out_dir = "outdir",
#'                  types = "HDF5")
#' }
countmat_to_hdf5 <- function(mat,
                             colsplit = NULL,
                             out_dir,
                             sub_dir = "filtered_feature_bc_matrix",
                             types = c("sparse", "HDF5"),
                             ...) {
  types <- rlang::arg_match(types, multiple = T)
  types <- sort(types, decreasing = T) # write sparse first

  if (!is.null(colsplit)) {
    mat <- split_mat(
      x = count_matrix,
      f = colsplit,
      byrow = F
    )
  } else {
    mat <- list(mat)
  }

  for (i in seq_along(mat)) {
    dir_to_create <- ifelse(is.null(names(mat)),
                            ifelse(is.null(sub_dir),
                                   out_dir,
                                   file.path(out_dir, sub_dir)),
                            ifelse(is.null(sub_dir),
                                   file.path(out_dir, names(mat)[i]),
                                   file.path(out_dir, names(mat)[i], sub_dir)))
    dir.create(dir_to_create, showWarnings = T, recursive = T)

    for (j in types) {
      path <- ifelse(j == "HDF5", file.path(dir_to_create, "filtered_feature_bc_matrix.h5"), dir_to_create)
      DropletUtils::write10xCounts(
        path = path,
        x = mat[[i]],
        version = "3",
        type = j,
        overwrite = T,
        ...
      )
      message(path)
    }
  }
}
