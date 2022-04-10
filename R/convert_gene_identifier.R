#' Title
#'
#' Convert gene identifiers from ident_in to ident_out with AnnotationDbi while taking care of known duplicates and providing aliases.
#' For human org.Hs.eg.db is used, for mouse org.Mm.eg.db.
#'
#' keytypes(org.Hs.eg.db) or colums(org.Hs.eg.db) are: ACCNUM ALIAS ENSEMBL
#' ENSEMBLPROT ENSEMBLTRANS ENTREZID ENZYME EVIDENCE EVIDENCEALL
#' GENENAME GENETYPE GO GOALL IPI MAP OMIM ONTOLOGY ONTOLOGYALL
#' PATH PFAM PMID PROSITE REFSEQ SYMBOL UCSCKG UNIPROT
#'
#' @param idents character vector of gene indentifiers
#' @param ident_in type of identifier provided in idents
#' @param ident_out type(s) of identifier(s) to return
#' @param species
#' @param return
#'
#' @return
#' @export
#'
#' @examples
convert_gene_identifier <- function (idents,
                                     ident_in = "SYMBOL",
                                     ident_out = c("ENTREZID", "ALIAS", "GENENAME"),
                                     species = c("Hs", "Mm"),
                                     return = c("data.frame", "vector")) {

  # https://medium.com/computational-biology/gene-id-mapping-using-r-14ff50eec9ba

  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
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

