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
#' @param gene_idents
#' @param ident_in
#' @param ident_out
#' @param species
#'
#' @return
#' @export
#'
#' @examples
convert_gene_identifier <- function (gene_idents,
                                     ident_in = "SYMBOL",
                                     ident_out = "ENTREZID",
                                     species = c("Hs", "Mm")) {

  species <- match.arg(species, c("Hs", "Mm"))
  if (species == "Hs") {
    if (!requireNamespace("org.Hs.eg.db", quietly = T)){
      BiocManager::install("org.Hs.eg.db")
    }
    library("org.Hs.eg.db")
    my.db <- org.Hs.eg.db
  }
  if (species == "Mm") {
    if (!requireNamespace("org.Mm.eg.db", quietly = T)){
      BiocManager::install("org.Mm.eg.db")
    }
    library("org.Mm.eg.db")
    my.db <- org.Mm.eg.db
  }

  ident_in <- match.arg(ident_in, keytypes(my.db))
  ident_out <- match.arg(ident_out, keytypes(my.db), several.ok = T)

  gene_idents <- unique(gene_idents)
  gene_idents <- gsub("^MT-", "MT", gene_idents, ignore.case = T)
  gene_idents <- AnnotationDbi::select(my.db, keys = as.character(gene_idents), keytype = ident_in, column = ident_out)

  # manual filtering of previously recognized duplicate matchings
  if ("SYMBOL" %in% names(gene_idents) && "ENTREZID" %in% names(gene_idents)) {
    gene_idents <-
      gene_idents %>%
      dplyr::filter(!(SYMBOL == "MEMO1" & ENTREZID == "7795")) %>%
      dplyr::filter(!(SYMBOL == "TEC" & ENTREZID == "100124696")) %>%
      dplyr::filter(!(SYMBOL == "MMD2" & ENTREZID == "100505381")) %>%
      dplyr::filter(!(SYMBOL == "HBD" & ENTREZID == "100187828"))
    print("Duplicated matches:")
    gene_idents[which(gene_idents$SYMBOL %in% gene_idents[which(duplicated(gene_idents$SYMBOL)), "SYMBOL"]), ]
  }

  if ("SYMBOL" %in% names(gene_idents)) {
    gene_idents$ALIAS <- limma::alias2SymbolTable(alias = gene_idents$SYMBOL, species = "Hs")
  }
  if ("SYMBOL" %in% names(gene_idents) && "ENTREZID" %in% names(gene_idents)) {
    gene_idents[intersect(which(is.na(gene_idents$ENTREZID)), which(!is.na(gene_idents$ALIAS))), "ENTREZID"] <- AnnotationDbi::select(my.db, keys = as.character(gene_idents[intersect(which(is.na(gene_idents$ENTREZID)), which(!is.na(gene_idents$ALIAS))), "ALIAS"]), keytype = ident_in, column = ident_out)$ENTREZID

  }

  print("Duplicated ENTREZIDs are filtered for distinct rows of ENTREZID and ALIAS:")
  gene_idents[which(gene_idents$ENTREZID %in% gene_idents[which(duplicated(gene_idents$ENTREZID, incomparables = NA)), "ENTREZID"]), ]
  print("Rows for which an ALIAS was found but no ENTREZID:")
  gene_idents[intersect(which(is.na(gene_idents$ENTREZID)), which(!is.na(gene_idents$ALIAS))), ]
  print("Rows with ENTREZID == NA will also be dropped")

  gene_idents <- dplyr::distinct(gene_idents, ENTREZID, ALIAS, .keep_all = T)
  gene_idents <- tidyr::drop_na(gene_idents, ENTREZID)
  gene_idents <- dplyr::left_join(gene_idents, AnnotationDbi::select(my.db, keys = as.character(gene_idents[which(!is.na(gene_idents$ALIAS)), "ALIAS"]), keytype = "SYMBOL", column = "GENENAME"), by = c("ALIAS" = "SYMBOL"))
  gene_idents <- dplyr::distinct(gene_idents, SYMBOL, ENTREZID, ALIAS, .keep_all = T)

  return(gene_idents)
}
