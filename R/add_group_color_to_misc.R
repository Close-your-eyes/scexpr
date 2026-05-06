#' Add group-specific colors to a Seurat object's misc slot
#'
#' Assigns a named vector of colors to groups defined in a metadata column
#' and stores the result in the \code{misc} slot of a Seurat object.
#'
#' @param obj A Seurat object containing a \code{meta.data} slot.
#' @param group Character scalar. Name of the metadata column used to define groups.
#'   Defaults to \code{"orig.ident"}.
#' @param colors A character vector of colors. Can be named (matching group levels)
#'   or unnamed. If unnamed or partially named, colors are assigned in order to
#'   sorted unique group values.
#'   Defaults to \code{colrr::col_pal("custom")}.
#' @param name Character scalar. Name under which the color vector will be stored
#'   in \code{obj@misc}. Defaults to \code{paste0(group, "_colors")}.
#'
#' @returns A modified Seurat object with a named color vector stored in
#'   \code{obj@misc[[name]]}.
#'
#' @details
#' If \code{colors} is a named vector and all names match the group levels,
#' it is used directly. Otherwise, colors are assigned sequentially to the
#' sorted unique values of the grouping variable.
#'
#' @export
#'
#' @examples
add_group_color_to_misc <- function(obj,
                                    group = "orig.ident",
                                    colors = colrr::col_pal("custom"),
                                    name = paste0(group, "_colors")) {

  ## no checking yet

  un_group <- sort(unique(obj@meta.data[[group]]))

  if (!is.null(names(colors)) && all(names(colors) %in% un_group)) {
    color_vec <- colors
  } else {
    color_vec <- stats::setNames(colors[1:length(un_group)], un_group)
  }


  obj@misc[[name]] <- color_vec
  return(obj)
}
