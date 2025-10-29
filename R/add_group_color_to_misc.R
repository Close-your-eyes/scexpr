#' Title
#'
#' @param obj
#' @param group
#' @param colors
#' @param name
#'
#' @returns
#' @export
#'
#' @examples
add_group_color_to_misc <- function(obj,
                                    group = "orig.ident",
                                    colors = colrr::col_pal("custom"),
                                    name = paste0(group, "_colors")) {

  ## no checking yet


  un_group <- sort(unique(obj@meta.data[[group]]))
  color_vec <- stats::setNames(colors[1:length(un_group)], un_group)
  obj@misc[[name]] <- color_vec
  return(obj)
}
