#' Find markers for subsets of a Seurat object
#'
#' This function splits a Seurat object by a metadata column or existing identities
#' and computes markers for each group using \code{find_all_marker}. Optionally,
#' specific groups can be included or excluded.
#'
#' @param obj A \code{Seurat} object.
#' @param split Character. Name of a metadata column used to redefine identities
#'   before marker detection. If \code{NULL}, existing identities are used.
#' @param levels Character vector. Subset of identity levels to include or exclude.
#' @param levels_invert Logical. If \code{FALSE} (default), only the specified
#'   \code{levels} are kept. If \code{TRUE}, the specified \code{levels} are excluded.
#' @param ... Deprecated. Additional arguments (not used).
#' @param find_all_marker_args List. Additional arguments passed to
#'   \code{find_all_marker}. Defaults to \code{list(meta_col = NULL)}.
#'
#' @returns A named list where each element corresponds to a group and contains
#'   the marker results returned by \code{find_all_marker}. Returns \code{NULL}
#'   if no groups are available after filtering.
#'
#' @details
#' The function optionally resets identities based on a metadata column and then
#' iterates over all identity groups (or a subset defined via \code{levels}).
#' For each group, the Seurat object is subsetted and passed to
#' \code{find_all_marker}.
#'
#' @export
#'
#' @examples
find_all_marker_split <- function(obj,
                                  split = NULL,
                                  levels = NULL,
                                  levels_invert = F,
                                  find_all_marker_args = list(meta_col = NULL)) {

  if (!is.null(split)) {
    SeuratObject::Idents(obj) <- as.character(obj@meta.data[[split]])
  } else {
    message("using idents: ", paste(as.character(unique(SeuratObject::Idents(obj))), collapse = ", "))
  }
  groups <- as.character(unique(SeuratObject::Idents(obj)))

  if (!is.null(levels)) {
    if (!levels_invert) {
      groups <- intersect(levels, groups)
    } else {
      groups <- groups[which(!groups %in% levels)]
    }
  }
  if (!length(groups)) {
    message("no groups left.")
    return(NULL)
  }
  groups <- purrr::set_names(groups)
  message("finding markers for: ", paste(groups, collapse = ", "))

  #obj <- Seurat::SplitObject(obj, split.by = split)[groups]
  markerlist <- purrr::map(groups, ~Gmisc::fastDoCall(find_all_marker,
                                                      args = c(list(obj = subset(obj, idents = .x)),
                                                               find_all_marker_args)))
  return(markerlist)
}
