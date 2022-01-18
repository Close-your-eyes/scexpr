#' Title
#'
#' Convert gene identifiers from ident_in to ident_out while taking care
#' of known duplicates and providing aliases.
#'
#' keytypes(org.Hs.eg.db) or colums(org.Hs.eg.db) are: ACCNUM ALIAS ENSEMBL
#' ENSEMBLPROT ENSEMBLTRANS ENTREZID ENZYME EVIDENCE EVIDENCEALL
#' GENENAME GENETYPE GO GOALL IPI MAP OMIM ONTOLOGY ONTOLOGYALL
#' PATH PFAM PMID PROSITE REFSEQ SYMBOL UCSCKG UNIPROT
#'
#' @param idents
#' @param ident_in
#' @param ident_out
#' @param species
#'
#' @return
#' @export
#'
#' @examples
convert_gene_identifier <- function (idents,
                                     ident_in = "SYMBOL",
                                     ident_out = "ENTREZID",
                                     species = c("Hs", "Mm")) {

  #ident_out <-c("ENTREZID", "ENSEMBL")

  species <- match.arg(species, c("Hs", "Mm"))
  ident_in <- match.arg(ident_in, AnnotationDbi::keytypes(my.db))
  ident_out <- match.arg(ident_out, AnnotationDbi::keytypes(my.db), several.ok = T)

  if (species == "Hs") {
    if (!requireNamespace("org.Hs.eg.db", quietly = T)){
      BiocManager::install("org.Hs.eg.db")
    }
    my.db <- org.Hs.eg.db::org.Hs.eg.db
    idents <- gsub("^MT-", "MT", idents, ignore.case = F)
  }
  if (species == "Mm") {
    if (!requireNamespace("org.Mm.eg.db", quietly = T)){
      BiocManager::install("org.Mm.eg.db")
    }
    my.db <- org.Mm.eg.db::org.Mm.eg.db
    idents <- gsub("^mt-", "mt", idents, ignore.case = F)
  }

  start_len <- length(idents)
  idents <- suppressMessages(AnnotationDbi::select(my.db, keys = as.character(unique(idents)), keytype = ident_in, column = ident_out))


  if (ident_in == "SYMBOL" && "ENTREZID" %in% ident_out) {
    # special treatment
    # yet wrong, fix!
    idents <- idents[intersect(which(!idents$SYMBOL == "MEMO1"), which(!idents$ENTREZID == "7795")),]
    idents <- idents[intersect(which(!idents$SYMBOL == "TEC"), which(!idents$ENTREZID == "100124696")),]
    idents <- idents[intersect(which(!idents$SYMBOL == "MMD2"), which(!idents$ENTREZID == "100505381")),]
    idents <- idents[intersect(which(!idents$SYMBOL == "HBD"), which(!idents$ENTREZID == "100187828")),]
    if (any(duplicated(idents[,ident_in]))) {
      print(paste0("Duplicate return for: ", paste(idents[,ident_in][which(duplicated(idents[,ident_in]))], collapse = ", ")))
      print("Made distinct with dplyr::distinct")
      idents <- dplyr::distinct(idents, !!sym(ident_in), .keep_all = T)
    }
    idents$ALIAS <- limma::alias2SymbolTable(alias = idents$SYMBOL, species = species)
    rows <- intersect(which(is.na(idents$ENTREZID)), which(!is.na(idents$ALIAS)))
    if (length(rows) > 0) {
      idents[rows, "ENTREZID"] <- AnnotationDbi::select(my.db, keys = as.character(idents[rows, "ALIAS"]), keytype = ident_in, column = "ENTREZID")$ENTREZID
    }

  } else {
    if (any(duplicated(idents[,ident_in]))) {
      print(paste0("Duplicate return for: ", paste(idents[,ident_in][which(duplicated(idents[,ident_in]))], collapse = ", ")))
      print("Made distinct with dplyr::distinct")
      idents <- dplyr::distinct(idents, !!sym(ident_in), .keep_all = T)
    }
  }
  idents$GENENAME <- dplyr::distinct(suppressMessages(AnnotationDbi::select(my.db, keys = idents[,ident_in], keytype = ident_in, column = "GENENAME")))
  if (nrow(idents) != start_len) {
    print("input length and output length are not identical.")
  }
  return(idents)
}
