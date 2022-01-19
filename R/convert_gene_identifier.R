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
                                     species = c("Hs", "Mm"),
                                     return = "df") {

  #ident_out <-c("ENTREZID", "ENSEMBL")

  return <- match.arg(return, c("df", "vector"))
  if (return == "vector" && length(ident_out) > 1) {
    print("ident_out has more than one entry, setting return to 'df'.")
    return <- "df"
  }
  species <- match.arg(species, c("Hs", "Mm"))
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
  ident_in <- match.arg(ident_in, AnnotationDbi::keytypes(my.db))
  ident_out <- match.arg(ident_out, AnnotationDbi::keytypes(my.db), several.ok = T)

  start_len <- length(idents)
  idents <- suppressMessages(AnnotationDbi::select(my.db, keys = as.character(unique(idents)), keytype = ident_in, column = ifelse(ident_in != "SYMBOL", unique(c(ident_out, "SYMBOL")), ident_out)))

  if ("SYMBOL" %in% names(idents) && "ENTREZID" %in% names(idents)) {
    idents <- dplyr::filter(idents, !(SYMBOL == "MEMO1" & ENTREZID == "7795"))
    idents <- dplyr::filter(idents, !(SYMBOL == "TEC" & ENTREZID == "100124696"))
    idents <- dplyr::filter(idents, !(SYMBOL == "MMD2" & ENTREZID == "100505381"))
    idents <- dplyr::filter(idents, !(SYMBOL == "HBD" & ENTREZID == "100187828"))
  }
  if (any(duplicated(idents[,ident_in]))) {
    print(paste0("Duplicate return for: ", paste(idents[,ident_in][which(duplicated(idents[,ident_in]))], collapse = ", ")))
    print("Made distinct with dplyr::distinct")
    idents <- dplyr::distinct(idents, !!sym(ident_in), .keep_all = T)
  }
  idents$ALIAS <- suppressWarnings(limma::alias2SymbolTable(alias = idents$SYMBOL, species = species))

  for (i in ident_out) {
    # use ALIAS to find ident_out
    rows <- intersect(which(is.na(idents[,i])), which(!is.na(idents[,"ALIAS"])))
    if (length(rows) > 0) {
      idents[rows, i] <- dplyr::distinct(suppressMessages(AnnotationDbi::select(my.db, keys = as.character(idents[rows, "ALIAS"]), keytype = "SYMBOL", column = i)[,i]))
    }
  }

  if (return == "df") {

    #names <- dplyr::distinct(suppressMessages(AnnotationDbi::select(my.db, keys = ifelse(is.na(idents[,ident_out[1]]), idents[,"ALIAS"], idents[,"SYMBOL"]), keytype = "SYMBOL", column = "GENENAME")))
    #names <- dplyr::distinct(suppressMessages(AnnotationDbi::select(my.db, keys = idents[,"SYMBOL"], keytype = "SYMBOL", column = "GENENAME")))
    #idents <- suppressMessages(dplyr::left_join(idents, names))
    #names <- dplyr::distinct(suppressMessages(AnnotationDbi::select(my.db, keys = idents[which(is.na(idents[,"GENENAME"])),"ALIAS"], keytype = "SYMBOL", column = "GENENAME")))
    #for (i in which(is.na(idents[,"GENENAME"]))) {idents[i,"GENENAME"] <- names[which(names$SYMBOL == idents[i,"ALIAS"]),"GENENAME"]}

    names <- dplyr::distinct(suppressMessages(AnnotationDbi::select(my.db, keys = idents[,idents[,ident_out[1]]], keytype = ident_out[1], column = "GENENAME")))
    idents <- suppressMessages(dplyr::left_join(idents, names))
    idents <- dplyr::distinct(idents, !!sym(ident_in), .keep_all = T)
    if (nrow(idents) != start_len) {
      print("input length and output length are not identical.")
    }
    return(idents)
  }
  if (return == "vector") {
    idents <- dplyr::distinct(idents, !!sym(ident_in), .keep_all = T)
    idents <- idents[,ident_out]
    if (length(idents) != start_len) {
      print("input length and output length are not identical.")
    }
    return(idents)
  }

}
