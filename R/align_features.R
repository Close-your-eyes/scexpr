#' Align shared features across dataframes using master gene symbols
#'
#' @description
#' Identify features shared by all input dataframes and return each dataframe's
#' original symbols together with a mapping to symbols from a master dataframe.
#'
#' @details
#' Features are matched by Ensembl identifier. Missing or blank identifiers are
#' resolved internally by symbol when that symbol maps to exactly one known
#' identifier across the inputs. An ambiguous symbol-to-identifier mapping
#' raises an error when needed for this fallback.
#'
#' Symbols without a known identifier in any input are matched directly by
#' symbol. Rows missing both an identifier and a symbol are excluded.
#'
#' If an identifier is duplicated in any input after fallback resolution,
#' matching uses both identifier and symbol. Only identifier-symbol combinations
#' present in every input are retained. Repeated identical combinations are
#' paired by occurrence, retaining the minimum number of occurrences across
#' inputs.
#'
#' For identifiers that are not duplicated, symbols may differ between inputs.
#' Returned features follow the master dataframe's order, and updated symbols
#' are taken from the master.
#'
#' Input dataframes and Seurat objects are not modified. Original symbols must
#' correspond to Seurat feature names to be used directly for filtering.
#' Missing or duplicated original symbols cannot uniquely identify features;
#' use actual Seurat feature names in those cases.
#'
#' @param feature_dfs A nonempty list of feature dataframes, typically one per
#'   Seurat object. Each must contain \code{id_col} and \code{symbol_col}.
#'   Supply unique, nonempty list names. If the list is unnamed, names are
#'   assigned from its positions: \code{"1"}, \code{"2"}, and so on.
#' @param master A single character string naming the master dataframe in
#'   \code{feature_dfs}. Defaults to \code{NULL}, which selects the first input.
#' @param id_col A single character string naming the Ensembl identifier
#'   column. Defaults to \code{"ENSEMBL"}.
#' @param symbol_col A single character string naming the gene symbol column.
#'   Defaults to \code{"SYMBOL"}.
#'
#' @return A named list with one element per input dataframe. Each element is
#'   a list containing:
#'   \describe{
#'     \item{old_symbols}{A character vector of original symbols for retained
#'       features, ordered to match the master.}
#'     \item{symbol_map}{A named character vector whose names are
#'       \code{old_symbols} and whose values are the corresponding master
#'       symbols. Names may be duplicated or missing if the original symbols
#'       are duplicated or missing.}
#'   }
#'   If no features are shared, both vectors are empty for every input.
#'
#' @examples
#' feature_dfs <- list(
#'   reference = data.frame(
#'     ENSEMBL = c("ENSG1", "ENSG1", "ENSG2", NA),
#'     SYMBOL = c("A", "B", "C", "D")
#'   ),
#'   sample1 = data.frame(
#'     ENSEMBL = c("ENSG1", "ENSG1", "ENSG2", "ENSG3"),
#'     SYMBOL = c("B", "X", "C_old", "D")
#'   )
#' )
#'
#' aligned <- align_features(feature_dfs, master = "reference")
#'
#' # ENSG1 retains only the shared symbol B.
#' # ENSG2 matches by ID and maps C_old to C.
#' # The missing reference ID for D resolves to ENSG3.
#' aligned$sample1$old_symbols
#' # c("B", "C_old", "D")
#'
#' aligned$sample1$symbol_map
#' # c(B = "B", C_old = "C", D = "D")
#'
#' @export
align_features <- function(feature_dfs,
                           master = NULL,
                           id_col = "ENSEMBL",
                           symbol_col = "SYMBOL") {

  stopifnot(length(feature_dfs) > 0L)

  if (is.null(names(feature_dfs))) {
    names(feature_dfs) <- as.character(seq_along(feature_dfs))
  }

  if (is.null(master)) {
    master <- names(feature_dfs)[1]
  }

  stopifnot(
    master %in% names(feature_dfs),
    all(vapply(feature_dfs, function(df) {
      all(c(id_col, symbol_col) %in% names(df))
    }, logical(1)))
  )

  clean <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & trimws(x) == ""] <- NA_character_
    x
  }

  ids <- lapply(feature_dfs, function(df) clean(df[[id_col]]))
  symbols <- lapply(feature_dfs, function(df) clean(df[[symbol_col]]))

  # Build symbol -> known Ensembl IDs across all objects.
  all_ids <- unlist(ids, use.names = FALSE)
  all_symbols <- unlist(symbols, use.names = FALSE)
  known <- !is.na(all_ids) & !is.na(all_symbols)
  symbol_ids <- lapply(
    split(all_ids[known], all_symbols[known]),
    unique
  )

  # Resolve missing IDs internally; preserve original dataframe columns.
  ids <- Map(function(id, symbol) {
    for (i in which(is.na(id) & !is.na(symbol))) {
      candidates <- symbol_ids[[symbol[i]]]

      if (length(candidates) > 1L) {
        stop("Ambiguous symbol fallback: ", symbol[i],
             " maps to multiple Ensembl IDs.")
      }
      if (length(candidates) == 1L) {
        id[i] <- candidates
      }
    }
    id
  }, ids, symbols)

  duplicated_ids <- unique(unlist(lapply(ids, function(id) {
    id[!is.na(id) & duplicated(id)]
  }), use.names = FALSE))

  # A duplicated ID must have matching symbol sets wherever it occurs.
  # for (id in duplicated_ids) {
  #   symbol_sets <- Map(function(x, s) {
  #     unique(s[which(x == id)])
  #   }, ids, symbols)
  #
  #   symbol_sets <- Filter(function(x) length(x) > 0L, symbol_sets)
  #
  #   if (!all(vapply(symbol_sets, function(x) {
  #     setequal(x, symbol_sets[[1L]])
  #   }, logical(1)))) {
  #     stop("Symbols do not align across objects for duplicated ID: ", id)
  #   }
  # }

  # Encode values so missing symbols remain distinct from literal text.
  encode <- function(x) {
    ifelse(is.na(x), "N", paste0("V", nchar(x), ":", x))
  }

  keys <- Map(function(id, symbol) {
    key <- rep(NA_character_, length(id))

    has_id <- !is.na(id)
    key[has_id] <- paste0("ID:", encode(id[has_id]))

    dup <- has_id & id %in% duplicated_ids
    key[dup] <- paste0(key[dup], "|SYMBOL:", encode(symbol[dup]))

    # Symbols with no known ID anywhere can still match each other.
    fallback <- !has_id & !is.na(symbol)
    key[fallback] <- paste0("SYMBOL:", encode(symbol[fallback]))

    # Pair repeated identical entries by occurrence, preserving order.
    valid <- which(!is.na(key))
    occurrence <- ave(valid, key[valid], FUN = seq_along)
    key[valid] <- paste0(key[valid], "|ROW:", occurrence)
    key
  }, ids, symbols)

  common <- Reduce(intersect, lapply(keys, function(x) x[!is.na(x)]))
  master_rows <- which(keys[[master]] %in% common)
  ordered_keys <- keys[[master]][master_rows]

  master_symbols <- as.character(
    feature_dfs[[master]][[symbol_col]][master_rows]
  )

  result <- lapply(seq_along(feature_dfs), function(i) {
    rows <- match(ordered_keys, keys[[i]])

    old_symbols <- as.character(
      feature_dfs[[i]][[symbol_col]][rows]
    )

    list(
      old_symbols = old_symbols,
      # Names = original symbols; values = master symbols.
      symbol_map = setNames(master_symbols, old_symbols)
    )
  })

  stats::setNames(result, names(feature_dfs))
}


#' Align shared gene features across Seurat objects
#'
#' @description
#' Convert feature names from each Seurat object into gene annotation
#' dataframes, then identify shared features and map their original symbols
#' to symbols from a master object using \code{align_features()}.
#'
#' @param obj_list A nonempty list of Seurat objects. List names identify
#'   objects in the output and can be used to select the master via
#'   \code{master}.
#' @param species species for gene lookup
#' @param ... Additional arguments passed to \code{align_features()},
#'   including \code{master}, \code{id_col}, and \code{symbol_col}.
#'
#' @details
#' Feature names are obtained using \code{rownames()} and annotated with
#' \code{scexpr::convert_gene_identifier()}, requesting \code{GENENAME}
#' and \code{ENSEMBL}.
#'
#' The variable \code{species} must be available in the function's enclosing
#' environment. The converted dataframes must contain the identifier and
#' symbol columns expected by \code{align_features()}; use \code{...} to
#' override its column defaults when necessary.
#'
#' Matching, missing-identifier fallback, duplicate handling, and output
#' ordering follow \code{align_features()}. Seurat objects are not modified.
#' Returned original symbols can be used directly for filtering only when
#' they correspond to the objects' feature names.
#'
#' @return A named list with one element per input object. Each element
#'   contains:
#'   \describe{
#'     \item{old_symbols}{Original symbols from the converted dataframe
#'       for shared features, ordered to match the master.}
#'     \item{symbol_map}{A named character vector with original symbols
#'       as names and corresponding master symbols as values.}
#'   }
#'
#' @seealso \code{\link{align_features}}
#' @export
align_features_seurat <- function(obj_list,
                                  species = "Hs",
                                  ...) {

  obj_list <- purrr::map(obj_list, function(x) {
    if ("RNA" %in% names(x@assays)) {
      SeuratObject::DefaultAssay(x) <- "RNA"
    } else {
      message("align_features_seurat: RNA assay not found.")
    }
    return(x)
  })

  assays <- purrr::map(obj_list, ~SeuratObject::DefaultAssay(.x))
  if (length(unique(assays)) != 1) {
    message("align_features_seurat: different assays used.")  }

  dfs <- purrr::map(obj_list, ~convert_gene_identifier(
    rownames(.x),
    species = species,
    ident_out = c("ENSEMBL")
  ))

  align_features(dfs, ...)
}
