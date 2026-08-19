#' Print a size tree for a Seurat object
#'
#' Prints the total memory size of a Seurat object and a tree-like summary of
#' its top-level slots using [lobstr::obj_size()].
#'
#' The `assays` slot is expanded to show every assay and each assay's layers.
#' The `graphs`, `reductions`, and `misc` slots are expanded by one level to
#' show their named elements. Other slots, including `meta.data`, are shown
#' only at the top level.
#'
#' Sizes are inclusive. For example, an assay's reported size includes its
#' layers. Because R objects can share memory, the sizes of child objects may
#' not sum exactly to the size of their parent.
#'
#' @param object A Seurat object.
#'
#' @return The input `object`, invisibly. The function is called primarily for
#'   its console output.
#'
#' @details
#' The output hierarchy is:
#'
#' * Total Seurat object
#' * All top-level slots
#' * All assays
#' * All layers within each assay
#' * All graph, reduction, and `misc` elements
#'
#' The function does not expand:
#'
#' * Assay slots other than layers
#' * Graph or reduction internals
#' * Columns of `meta.data`
#' * Internals of `misc` elements
#' * Other top-level slots
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(lobstr)
#'
#' seurat_size_tree(pbmc_small)
#' }
#'
#' @export
seurat_size_tree <- function(object) {
  size_text <- function(x) {
    format(lobstr::obj_size(x), units = "auto")
  }

  print_elements <- function(x, prefix) {
    element_names <- names(x)

    if (is.null(element_names)) {
      element_names <- as.character(seq_along(x))
    }

    empty <- is.na(element_names) | element_names == ""
    element_names[empty] <- which(empty)

    for (i in seq_along(x)) {
      branch <- if (i == length(x)) "└── " else "├── "

      cat(
        prefix, branch,
        element_names[[i]], " — ",
        size_text(x[[i]]), "\n",
        sep = ""
      )
    }
  }

  top_slots <- methods::slotNames(object)

  cat("Total size — ", size_text(object), "\n", sep = "")

  for (i in seq_along(top_slots)) {
    slot_name  <- top_slots[[i]]
    slot_value <- methods::slot(object, slot_name)
    last_slot  <- i == length(top_slots)

    branch <- if (last_slot) "└── " else "├── "
    prefix <- if (last_slot) "    " else "│   "

    cat(
      branch, "@", slot_name,
      " — ", size_text(slot_value), "\n",
      sep = ""
    )

    # Show assays and their layers.
    if (slot_name == "assays" && length(slot_value) > 0L) {
      assay_names <- names(slot_value)

      for (j in seq_along(slot_value)) {
        assay <- slot_value[[j]]
        last_assay <- j == length(slot_value)

        assay_branch <- if (last_assay) "└── " else "├── "
        layer_prefix <- paste0(
          prefix,
          if (last_assay) "    " else "│   "
        )

        cat(
          prefix, assay_branch,
          assay_names[[j]], " — ",
          size_text(assay), "\n",
          sep = ""
        )

        layer_names <- SeuratObject::Layers(assay)

        for (k in seq_along(layer_names)) {
          layer_name <- layer_names[[k]]
          layer <- SeuratObject::LayerData(
            object = assay,
            layer = layer_name
          )

          layer_branch <-
            if (k == length(layer_names)) "└── " else "├── "

          cat(
            layer_prefix, layer_branch,
            layer_name, " — ",
            size_text(layer), "\n",
            sep = ""
          )
        }
      }
    }

    # Show immediate elements, without their internals.
    if (
      slot_name %in% c("graphs", "reductions", "misc") &&
      length(slot_value) > 0L
    ) {
      print_elements(slot_value, prefix)
    }
  }

  invisible(object)
}
