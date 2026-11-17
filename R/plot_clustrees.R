#' Plot clustering trees from a Seurat object
#'
#' Visualize how cells move between clusters across clustering resolutions.
#'
#' @param obj A Seurat object whose `misc$clusterings` entry is a character
#'   vector of clustering column names in `meta.data`.
#'
#' @details
#' The function removes the trailing numeric resolution from each listed
#' column name to find its prefix. For example, `RNA_snn_res.0.25` and
#' `RNA_snn_res.1` both have the prefix `RNA_snn_res.`.
#'
#' Each tree uses all metadata columns that start with its prefix. A prefix
#' must have at least two distinct resolutions, with one column per
#' resolution and a cluster assignment for every cell.
#'
#' @return A named list of ggplot objects, one per clustering prefix.
#'
#' @seealso \code{\link[clustree]{clustree}}
#'
#' @examples
#' \dontrun{
#' clustrees <- plot_clustrees(obj)
#' clustrees[["RNA_snn_res."]]
#' }
#' @export
plot_clustrees <- function(obj) {

  scexpr:::.ensure_package("clustree")

  if (!methods::is(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }

  if (!all(c("misc", "meta.data") %in% methods::slotNames(obj))) {
    stop("`obj` must have `misc` and `meta.data` slots.", call. = FALSE)
  }

  clusterings <- obj@misc[["clusterings"]]
  if (!is.character(clusterings) || length(clusterings) == 0L ||
      anyNA(clusterings) || any(!nzchar(trimws(clusterings)))) {
    stop("`obj@misc$clusterings` must contain nonblank column names.",
         call. = FALSE)
  }

  metadata <- obj@meta.data
  if (nrow(metadata) == 0L || anyNA(names(metadata)) ||
      any(!nzchar(names(metadata))) || anyDuplicated(names(metadata))) {
    stop("`obj@meta.data` must have cells and unique, nonblank column names.",
         call. = FALSE)
  }

  missing_columns <- setdiff(clusterings, names(metadata))
  if (length(missing_columns) > 0L) {
    stop("Clustering columns missing from `obj@meta.data`: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  resolution_pattern <- "[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$"
  prefixes <- sub(resolution_pattern, "", clusterings)
  bad_names <- !nzchar(prefixes) | prefixes == clusterings
  if (any(bad_names)) {
    stop("Clustering names need a prefix followed by a numeric resolution: ",
         paste(clusterings[bad_names], collapse = ", "), call. = FALSE)
  }
  prefixes <- unique(prefixes)

  edge_guide <- ggplot2::guides(
    edge_colour = ggraph::guide_edge_colourbar()
  )

  lapply(stats::setNames(prefixes, prefixes), function(prefix) {
    columns <- names(metadata)[startsWith(names(metadata), prefix)]
    resolutions <- suppressWarnings(
      as.numeric(substring(columns, nchar(prefix) + 1L))
    )

    invalid <- !is.finite(resolutions)
    if (any(invalid)) {
      stop("Invalid resolutions in metadata columns: ",
           paste(columns[invalid], collapse = ", "), call. = FALSE)
    }

    if (length(resolutions) < 2L || anyDuplicated(resolutions)) {
      stop("Prefix '", prefix, "' needs at least two distinct resolutions ",
           "and one column per resolution.", call. = FALSE)
    }

    valid_labels <- vapply(metadata[columns], function(labels) {
      is.atomic(labels) && is.null(dim(labels)) && !anyNA(labels)
    }, logical(1))

    if (any(!valid_labels)) {
      stop("Missing or invalid cluster assignments in: ",
           paste(columns[!valid_labels], collapse = ", "), call. = FALSE)
    }

    clustree::clustree(
      metadata,
      prefix = prefix,
      node_colour = "grey50",
      node_text_size = 4,
      show_axis = TRUE
    ) + edge_guide
  })
}
