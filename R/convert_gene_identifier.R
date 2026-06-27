#' Convert gene identifiers
#'
#' Convert gene identifiers between supported identifier types using
#' `AnnotationDbi` organism annotation databases. Human genes are mapped with
#' `org.Hs.eg.db`; mouse genes are mapped with `org.Mm.eg.db`.
#'
#' The function also handles a small set of known duplicate mappings and, when
#' possible, uses `limma::alias2SymbolTable()` to resolve aliases before filling
#' missing output identifiers.
#'
#' @param idents Character vector of gene identifiers to convert.
#' @param ident_in Input identifier type. Must be one of
#'   `AnnotationDbi::keytypes()` for the selected organism database. Defaults
#'   to `"SYMBOL"`.
#' @param ident_out Character vector of output identifier types. Must be one or
#'   more valid key types for the selected organism database. Defaults to
#'   `c("ENTREZID", "ALIAS", "GENENAME")`.
#' @param species Species annotation database to use. Either `"Hs"` for human
#'   or `"Mm"` for mouse. If omitted and `ident_in` is `"SYMBOL"` or `"ALIAS"`,
#'   the function tries to infer species from gene-name capitalisation.
#' @param return Return format. Either `"data.frame"` or `"vector"`. Vector
#'   output is only used when a single `ident_out` is requested.
#'
#' @return
#' If `return = "data.frame"`, a data frame containing the input identifier
#' column and the requested output identifier columns. If `return = "vector"`,
#' a vector containing the requested converted identifier.
#'
#' @details
#' Valid identifier types depend on the selected annotation database and can be
#' inspected with `AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db)` or
#' `AnnotationDbi::keytypes(org.Mm.eg.db::org.Mm.eg.db)`.
#'
#' Common key types include `"SYMBOL"`, `"ALIAS"`, `"ENSEMBL"`, `"ENTREZID"`,
#' `"GENENAME"`, `"REFSEQ"` and `"UNIPROT"`.
#'
#' For human symbols, mitochondrial identifiers beginning with `"MT-"` are
#' normalised to `"MT"`. For mouse symbols, identifiers beginning with `"mt-"`
#' are normalised to `"mt"`.
#'
#' When duplicate mappings are returned by `AnnotationDbi::select()`, the first
#' distinct mapping is retained after reporting the affected identifiers.
#'
#' @examples
#' convert_gene_identifier(
#'   idents = c("MS4A1", "CD3D", "LYZ"),
#'   ident_in = "SYMBOL",
#'   ident_out = c("ENTREZID", "GENENAME"),
#'   species = "Hs"
#' )
#'
#' convert_gene_identifier(
#'   idents = c("Ms4a1", "Cd3d", "Lyz2"),
#'   ident_in = "SYMBOL",
#'   ident_out = "ENTREZID",
#'   species = "Mm",
#'   return = "vector"
#' )
#'
#' @export
convert_gene_identifier <- function (idents,
                                     ident_in = "SYMBOL",
                                     ident_out = c("ENTREZID", "ALIAS", "GENENAME"),
                                     species = c("Hs", "Mm"),
                                     return = c("data.frame", "vector")) {

  # https://medium.com/computational-biology/gene-id-mapping-using-r-14ff50eec9ba

  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("limma", quietly = T)) {
    BiocManager::install("limma")
  }


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

  species <- match.arg(species, c("Hs", "Mm"))
  return <- match.arg(return, c("data.frame", "vector"))

  if (return == "vector" && length(ident_out) > 1) {
    print("ident_out has more than one entry, setting return to 'data.frame'.")
    return <- "data.frame"
  }
  if (species == "Hs") {
    if (!requireNamespace("org.Hs.eg.db", quietly = T)) {
      BiocManager::install("org.Hs.eg.db")
    }
    my.db <- org.Hs.eg.db::org.Hs.eg.db
    idents <- gsub("^MT-", "MT", idents, ignore.case = F)
  }
  if (species == "Mm") {
    if (!requireNamespace("org.Mm.eg.db", quietly = T)) {
      BiocManager::install("org.Mm.eg.db")
    }
    my.db <- org.Mm.eg.db::org.Mm.eg.db
    idents <- gsub("^mt-", "mt", idents, ignore.case = F)
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
  if (any(duplicated(idents[,ident_in]))) {
    print(paste0("Duplicate return by ident_out for: ", paste(idents[,ident_in][which(duplicated(idents[,ident_in]))], collapse = ", ")))
    print("Made distinct with dplyr::distinct")
    idents <- dplyr::distinct(idents, !!rlang::sym(ident_in), .keep_all = T)
  }

  if ("SYMBOL" %in% names(idents)) {
    idents$ALIAS <- suppressWarnings(limma::alias2SymbolTable(alias = idents$SYMBOL, species = species))
    for (i in ident_out) {
      # use ALIAS to find ident_out
      rows <- intersect(which(is.na(idents[,i])), which(!is.na(idents[,"ALIAS"])))
      if (length(rows) > 0) {
        idents[rows, i] <- suppressMessages(AnnotationDbi::mapIds(my.db, keys = as.character(idents[rows, "ALIAS"]), keytype = "SYMBOL", column = i, multiVals = "first"))
      }
    }
  }

  if (nrow(idents) != start_len) {
    print("input length and output length are not identical.")
  }

  if (return == "data.frame") {
    return(idents)
  }
  if (return == "vector") {
    idents <- dplyr::distinct(idents, !!rlang::sym(ident_in), .keep_all = T)
    idents <- idents[,ident_out]
    if (length(idents) != start_len) {
      print("input length and output length are not identical.")
    }
    return(idents)
  }

}

