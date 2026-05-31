#' Add meta_col-specific colors to a Seurat object's misc slot
#'
#' Assigns a named vector of colors to groups defined in a metadata column
#' and stores the result in the \code{misc} slot of a Seurat object.
#'
#' @param obj A Seurat object containing a \code{meta.data} slot.
#' @param meta_col Character scalar. Name of the metadata column used to define groups.
#'   Defaults to \code{"orig.ident"}.
#' @param colors A character vector of colors. Can be named (matching meta_col levels)
#'   or unnamed. If unnamed or partially named, colors are assigned in order to
#'   sorted unique meta_col values.
#'   Defaults to \code{colrr::col_pal("custom")}.
#' @returns A modified Seurat object.
#'
#' @details
#' If \code{colors} is a named vector and all names match the meta_col levels,
#' it is used directly. Otherwise, colors are assigned sequentially to the
#' sorted unique values of the grouping variable.
#'
#' @export
#'
#' @examples
add_group_color_to_misc <- function(obj,
                                    meta_col = "orig.ident",
                                    colors = colrr::col_pal("custom", return = "c")) {

  ## update_metacolor_slot
  ## no checking yet

  un_group <- sort(unique(obj@meta.data[[meta_col]]))

  if (!is.null(names(colors)) && all(names(colors) %in% un_group)) {
    color_vec <- colors
  } else {
    if (length(colors) < length(un_group)) {
      message("more colors required. using hue pal.")
      colors <- colrr::col_pal("hue", n = length(un_group), return = "c")
    }
    color_vec <- stats::setNames(colors[1:length(un_group)], un_group)
  }

  obj@misc[["metacolors"]][[meta_col]] <- color_vec

  return(obj)
}
