#' Merge RNA counts and metadata from Seurat objects
#'
#' @description
#' Align shared features across Seurat objects, combine their RNA counts and
#' cell metadata, and normalize the merged counts.
#'
#' @param obj_list A nonempty list of Seurat objects with RNA counts. List names
#'   identify the source objects and must be unique, nonmissing, and nonempty.
#'   Unnamed lists receive positional names, such as \code{"1"}, \code{"2"}.
#' @param obj_ident_col Name of the metadata column recording each cell's source
#'   object, using the list names. An existing column with this name is overwritten.
#' @param merge_reductions Logical; whether to stack cell embeddings for reduction
#'   names present in every object. Defaults to \code{FALSE}. See Details for
#'   coordinate-system requirements.
#' @param misc_from_first Logical; whether to copy named entries of the first
#'   object's \code{misc} slot. Defaults to \code{FALSE}.
#' @param align_features_species Species used for gene annotation during feature
#'   alignment: \code{"Hs"} for human (default) or \code{"Mm"} for mouse.
#'
#' @details
#' Features are aligned by \code{align_features_seurat()}, using the first object
#' as the master. Only features shared according to that alignment are retained;
#' output feature symbols and order follow the master. An error is raised if no
#' shared features remain, master symbols are missing, empty, or duplicated, or
#' aligned features cannot be matched uniquely to RNA counts rows.
#'
#' Before alignment, the function prints the fraction of each object's RNA-count
#' feature names present in every object. These fractions use exact feature-name
#' matches before annotation-based alignment, so they may differ from the final
#' retained fractions. An empty feature set is reported as zero.
#'
#' If cell names overlap across objects, all cells are prefixed with their source
#' list name and an underscore. Any remaining collisions are resolved with
#' \code{make.unique()}. Counts, metadata, and reductions are renamed together.
#' Otherwise, original cell names are retained. Input objects are not modified.
#'
#' A new Seurat object is constructed from the retained RNA counts and combined
#' metadata. \code{Seurat::NormalizeData()} is then run with its default
#' normalization settings. Other assays, existing normalized or scaled expression,
#' graphs, and reduction loadings are not copied.
#'
#' With \code{merge_reductions = TRUE}, common reductions are combined by stacking
#' their existing cell embeddings. This is meaningful only when inputs already
#' share the same coordinate system: independently computed PCA or UMAP embeddings
#' are not aligned by this operation. Embedding columns must be compatible; the
#' function does not explicitly validate this. The assay reference and key are
#' taken from the first object. A message explains this limitation when common
#' reductions are merged.
#'
#' @returns A Seurat object containing merged RNA counts, newly normalized data,
#'   combined cell metadata, and the source-object column. Common reduction
#'   embeddings and the first object's named \code{misc} entries are included
#'   only when requested.
#' @seealso \code{\link{align_features_seurat}}, \code{\link{get_layer}}
#' @export
#'
#' @examples
#' \dontrun{
#' merged <- merge_objects(list(sample1 = so1, sample2 = so2))
#'
#' # Mouse data; source labels are stored in the sample column.
#' merged_mouse <- merge_objects(
#'   list(sample1 = mouse1, sample2 = mouse2),
#'   obj_ident_col = "sample",
#'   align_features_species = "Mm"
#' )
#' }
merge_objects <- function(obj_list,
                          obj_ident_col = "obj_ident",
                          merge_reductions = F,
                          misc_from_first = F,
                          align_features_species = "Hs") {

  ## currently focused on RNA assay
  # no messages or warnings

  if (!is.list(obj_list)) {
    stop("obj_list must be a list of seurat objects.")
  }

  if (!length(obj_list)) {
    stop("obj_list must contain at least one Seurat object.")
  }
  if (is.null(names(obj_list))) {
    names(obj_list) <- as.character(seq_along(obj_list))
  }
  if (anyNA(names(obj_list)) || any(!nzchar(trimws(names(obj_list)))) ||
      anyDuplicated(names(obj_list))) {
    stop("obj_list names must be nonmissing, nonempty, and unique.")
  }

  # Rename whole objects so counts, metadata, and reductions stay in sync.
  cell_names <- unlist(lapply(obj_list, colnames), use.names = FALSE)
  if (anyDuplicated(cell_names)) {
    prefixed <- purrr::imap(obj_list, ~paste0(.y, "_", colnames(.x)))
    unique_names <- make.unique(unlist(prefixed, use.names = FALSE))
    offsets <- c(0L, cumsum(lengths(prefixed)))
    obj_list <- lapply(seq_along(obj_list), function(i) {
      SeuratObject::RenameCells(
        obj_list[[i]],
        new.names = unique_names[offsets[i] + seq_along(prefixed[[i]])]
      )
    }) |> stats::setNames(names(prefixed))
  }


  featlist <- purrr::map(obj_list, ~rownames(get_layer(.x, layer = "counts", assay = "RNA")))
  common_feat <- purrr::reduce(featlist, intersect)

  shared_feat_freq <- purrr::map_dbl(featlist, function(x) {
    if (length(x)) {
      length(intersect(x, common_feat))/length(x)
    } else {
      0
    }
  })
  message("shared feature frequencies:")
  print(shared_feat_freq)


  aligned <- align_features_seurat(
    obj_list,
    species = align_features_species,
    master = names(obj_list)[1]
  )


  master_symbols <- unname(aligned[[1]]$symbol_map)

  if (!length(master_symbols)) {
    stop("No shared features found.")
  }
  if (anyNA(master_symbols) ||
      any(!nzchar(master_symbols)) ||
      anyDuplicated(master_symbols)) {
    stop("Master feature symbols must be nonmissing, nonempty, and unique.")
  }

  counts_list <- purrr::imap(obj_list, function(obj, nm) {
    mat <- get_layer(obj, assay = "RNA", layer = "counts")
    mapping <- aligned[[nm]]

    # Match explicitly: get_layer(features = ...) silently drops missing genes.
    idx <- match(mapping$old_symbols, rownames(mat))
    # browser()
    # rownames(mat)[which(is.na(idx))]
    # which(rownames(mat) == "PWWP4")

    if (anyNA(idx)) {
      stop("Mapped symbols are absent from RNA counts in object: ", nm)
    }
    if (anyDuplicated(idx)) {
      stop("Multiple aligned features reference the same counts row in: ", nm)
    }

    # Subset into master order, then assign master symbols.
    mat <- mat[idx, , drop = FALSE]
    rownames(mat) <- unname(mapping$symbol_map)
    mat
  })

  stopifnot(all(vapply(
    counts_list,
    function(mat) identical(rownames(mat), master_symbols),
    logical(1)
  )))

  obj_merge <- Seurat::CreateSeuratObject(
    counts = do.call(cbind, unname(counts_list)),
    meta.data = purrr::map_dfr(obj_list, ~.x@meta.data)
  )

  # obj_merge <- Seurat::CreateSeuratObject(counts = do.call(cbind, purrr::map(obj_list, get_layer, assay = "RNA", layer = "counts", features = common_feat)),
  #                                         meta.data = purrr::map_dfr(obj_list, ~.x@meta.data))

  # dataslotthere <- tryCatch(expr = {
  #   purrr::map_lgl(obj_list, ~"data" %in% names(.x@assays$RNA@layers))
  # },
  # error = function(err) {
  #   # old object type
  #   purrr::map_lgl(obj_list, ~"data" %in% methods::slotNames(obj_list[[2]]@assays$RNA))
  # })
  #
  # if (all(dataslotthere)) {
  #   obj_merge@assays$RNA@layers[["data"]] <- do.call(cbind, purrr::map(obj_list, get_layer, assay = "RNA", layer = "data", features = common_feat))
  # } else {
  #   obj_merge <- Seurat::NormalizeData(obj_merge)
  # }

  # just normalize. safer. and is quick.
  obj_merge <- Seurat::NormalizeData(obj_merge, verbose = F)

  if (merge_reductions) {
    ## no check for equal columns of reduction
    redlst <- purrr::map(obj_list, ~names(.x@reductions))
    if (any(is.null(redlst))) {
      message("some obj w/o reductions.")
    } else {
      redcommon <- purrr::reduce(redlst, intersect)
      if (!length(redcommon)) {
        message("no common reductions")
      } else {
        message(
          "Merging reductions by stacking existing embeddings. ",
          "This is meaningful only when inputs share the same coordinate system; ",
          "independently computed PCA/UMAP embeddings are not aligned by this operation."
        )
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
