#' Convert gene identifiers
#'
#' @description
#' Map gene identifiers using \pkg{AnnotationDbi} and the organism annotation
#' databases \pkg{org.Hs.eg.db} for human or \pkg{org.Mm.eg.db} for mouse.
#' Optionally collapse alias and Ensembl mappings and retain one row per
#' input identifier.
#'
#' @param idents Vector of gene identifiers. Values are converted to character,
#'   and only unique input values are submitted for annotation.
#' @param ident_in Input identifier type. Must be a valid key type in the
#'   selected organism database. Defaults to \code{"SYMBOL"}.
#' @param ident_out Character vector of requested output identifier types.
#'   Each must be a valid key type in the selected organism database.
#'   Defaults to \code{c("ENTREZID", "ALIAS", "GENENAME")}.
#' @param species Species database: \code{"Hs"} for human or \code{"Mm"}
#'   for mouse. If omitted for \code{ident_in = "SYMBOL"} or \code{"ALIAS"},
#'   all-uppercase identifiers select human; otherwise, identifiers whose
#'   characters after the first are lowercase select mouse. An error is
#'   raised if neither rule applies. For other input types, defaults to
#'   \code{"Hs"}.
#' @param return Output format: \code{"data.frame"} (default) or
#'   \code{"vector"}. If multiple output types are requested,
#'   \code{"vector"} is changed to \code{"data.frame"} with a message.
#' @param collapse_alias How to combine values in the \code{ALIAS} column:
#'   \code{"paste"} (default) joins sorted values with \code{", "};
#'   \code{"list"} stores sorted values in a list-column; \code{"not"}
#'   leaves mappings uncollapsed. Rows are grouped by all other columns.
#'   Applies whenever an \code{ALIAS} column is present.
#' @param collapse_ensembl How to combine values in the \code{ENSEMBL}
#'   column. Accepts \code{"paste"} (default), \code{"list"}, or
#'   \code{"not"}, with the same behavior as \code{collapse_alias}.
#'   Applied after alias collapsing, whenever an \code{ENSEMBL}
#'   column is present.
#' @param make_distinct Logical. If \code{TRUE} (default), retain the first
#'   row per input identifier after filtering and collapsing, and report
#'   duplicated identifiers. If \code{FALSE}, retain remaining multiple
#'   mappings in data-frame output. Vector mode always retains only the
#'   first row per input identifier.
#'
#' @details
#' Supported identifier types can be inspected with
#' \code{AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)} or
#' \code{AnnotationDbi::keytypes(org.Mm.eg.db::org.Mm.eg.db)}.
#'
#' When both \code{SYMBOL} and \code{ENTREZID} are present, the following
#' mappings are excluded: \code{MEMO1}/\code{7795},
#' \code{TEC}/\code{100124696}, \code{MMD2}/\code{100505381}, and
#' \code{HBD}/\code{100187828}. The filtering expressions also drop rows
#' for which the exclusion condition evaluates to \code{NA}.
#'
#' Collapsing sorts values without removing duplicates. Missing values
#' are removed by sorting; an all-missing group therefore produces an
#' empty string in \code{"paste"} mode or an empty vector in
#' \code{"list"} mode.
#'
#' Input identifier strings are passed to annotation lookup without
#' mitochondrial-prefix rewriting. No additional alias-resolution step
#' is performed to fill missing annotations.
#'
#' Output is not guaranteed to have one entry per original input or to
#' preserve input order. Repeated inputs are queried only once, and
#' filtering, collapsing, or multiple mappings may change the number of
#' returned rows. A message is emitted when the resulting size differs
#' from the input length.
#'
#' @return
#' For \code{return = "data.frame"}, a data frame or tibble containing
#' the input identifier column and requested annotation columns.
#' Collapsed columns contain character values or lists, depending on
#' the corresponding collapse option.
#'
#' For \code{return = "vector"}, the requested output column after
#' retaining the first row per input identifier. The current implementation
#' uses single-bracket column extraction, so a tibble result remains a
#' one-column tibble rather than becoming a vector.
#'
#' @examples
#' \dontrun{
#' convert_gene_identifier(
#'   idents = c("MS4A1", "CD3D", "LYZ"),
#'   ident_out = c("ENTREZID", "GENENAME"),
#'   species = "Hs"
#' )
#'
#' # Retain separate Ensembl mappings for feature alignment.
#' convert_gene_identifier(
#'   idents = c("MS4A1", "CD3D", "LYZ"),
#'   ident_out = "ENSEMBL",
#'   species = "Hs",
#'   collapse_ensembl = "not",
#'   make_distinct = FALSE
#' )
#'
#' convert_gene_identifier(
#'   idents = c("Ms4a1", "Cd3d", "Lyz2"),
#'   ident_out = "ENTREZID",
#'   species = "Mm",
#'   return = "vector"
#' )
#' }
#'
#' @export
convert_gene_identifier <- function (idents,
                                     ident_in = "SYMBOL",
                                     ident_out = c("ENTREZID", "ALIAS", "GENENAME"),
                                     species = c("Hs", "Mm"),
                                     return = c("data.frame", "vector"),
                                     collapse_alias = c("paste", "list", "not"),
                                     collapse_ensembl = c("paste", "list", "not"),
                                     make_distinct = TRUE) {

  # https://medium.com/computational-biology/gene-id-mapping-using-r-14ff50eec9ba


  scexpr:::.ensure_package("limma")
  scexpr:::.ensure_package("AnnotationDbi")


  if (missing(species) && ident_in %in% c("SYMBOL", "ALIAS")) {
    # guess the species by case of letters
    if (all(idents == toupper(idents))) {
      species <- "Hs"
    } else if (all(stringr::str_sub(idents, 2) == tolower(stringr::str_sub(idents, 2)))) {
      species <- "Mm"
    } else {
      stop("convert_gene_identifier: Species could not be guessed. Please provide species = Hs or species = Mm.")
    }
  }

  species <- rlang::arg_match(species)
  return <- rlang::arg_match(return)
  collapse_alias <- rlang::arg_match(collapse_alias)
  collapse_ensembl <- rlang::arg_match(collapse_ensembl)

  if (return == "vector" && length(ident_out) > 1) {
    message("ident_out has more than one entry, setting return to 'data.frame'.")
    return <- "data.frame"
  }
  if (species == "Hs") {
    scexpr:::.ensure_package("org.Hs.eg.db")
    my.db <- org.Hs.eg.db::org.Hs.eg.db
    #idents <- gsub("^MT-", "MT", idents, ignore.case = F)
  }
  if (species == "Mm") {
    scexpr:::.ensure_package("org.Mm.eg.db")
    my.db <- org.Mm.eg.db::org.Mm.eg.db
    #idents <- gsub("^mt-", "mt", idents, ignore.case = F)
  }
  ident_in <- match.arg(ident_in, AnnotationDbi::keytypes(my.db))
  ident_out <- match.arg(ident_out, AnnotationDbi::keytypes(my.db), several.ok = T)

  start_len <- length(idents)
  # use select here to get multiple hits which can be filtered below
  idents <- suppressMessages(AnnotationDbi::select(my.db, keys = as.character(unique(idents)), keytype = ident_in, column = ident_out))

  if ("SYMBOL" %in% names(idents) && "ENTREZID" %in% names(idents)) {
    idents <- dplyr::filter(idents, !(SYMBOL == "MEMO1" & ENTREZID == "7795"))
    idents <- dplyr::filter(idents, !(SYMBOL == "TEC" & ENTREZID == "100124696"))
    idents <- dplyr::filter(idents, !(SYMBOL == "MMD2" & ENTREZID == "100505381"))
    idents <- dplyr::filter(idents, !(SYMBOL == "HBD" & ENTREZID == "100187828"))
  }

  if ("ALIAS" %in% names(idents)) {
    if (collapse_alias == "paste") {
      idents <- dplyr::summarise(idents, ALIAS = paste(sort(ALIAS), collapse = ", "), .by = -ALIAS)
    } else if (collapse_alias == "list") {
      idents <- dplyr::summarise(idents, ALIAS = list(sort(ALIAS)), .by = -ALIAS)
    }
  }
  if ("ENSEMBL" %in% names(idents)) {
    if (collapse_ensembl == "paste") {
      idents <- dplyr::summarise(idents, ENSEMBL = paste(sort(ENSEMBL), collapse = ", "), .by = -ENSEMBL)
    } else if (collapse_ensembl == "list") {
      idents <- dplyr::summarise(idents, ENSEMBL = list(sort(ENSEMBL)), .by = -ENSEMBL)
    }
  }

  if (anyDuplicated(idents[[ident_in]]) && make_distinct) {
    print(paste0("Duplicate return by ident_out for: ", paste(idents[[ident_in]][which(duplicated(idents[[ident_in]]))], collapse = ", ")))
    message("Made distinct with dplyr::distinct")
    idents <- dplyr::distinct(idents, !!rlang::sym(ident_in), .keep_all = T)
  }

  # if ("SYMBOL" %in% names(idents)) {
  #   idents$ALIAS <- suppressWarnings(limma::alias2SymbolTable(alias = idents$SYMBOL, species = species))
  #   for (i in ident_out) {
  #     # use ALIAS to find ident_out
  #     rows <- intersect(which(is.na(idents[,i])), which(!is.na(idents[,"ALIAS"])))
  #     if (length(rows) > 0) {
  #       idents[rows, i] <- suppressMessages(AnnotationDbi::mapIds(my.db, keys = as.character(idents[rows, "ALIAS"]), keytype = "SYMBOL", column = i, multiVals = "first"))
  #     }
  #   }
  # }


  if (return == "data.frame") {
    if (nrow(idents) != start_len) {
      message("input length and output length are not identical.")
    }

    return(idents)
  }
  if (return == "vector") {
    idents <- dplyr::distinct(idents, !!rlang::sym(ident_in), .keep_all = T)
    idents <- idents[[ident_out]]
    if (length(idents) != start_len) {
      message("input length and output length are not identical.")
    }
    return(idents)
  }

}
