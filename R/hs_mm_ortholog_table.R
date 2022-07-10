#' Create an table of ortholog genes between human and mouse
#'
#' A combined table with entries from ProjecTILs::Hs2Mm.convert.table and
#' entries from a bioMart query is generated. features_hs is used
#' to query bioMart. So, it is intended to look for mouse orthologs based
#' on human features.
#'
#' @param features_hs vector of human genes (gene symbols, hgnc)
#' @param features_mm vector of human genes (gene symbols)
#'
#' @importFrom magrittr %>%
#'
#' @return a list with (i) a table of orthologs and (ii, iii) vectors of genes which no orthologs were found for
#' @export
#'
#' @examples
hs_mm_ortholog_table <- function(features_hs,
                                 features_mm) {


  if (!requireNamespace("remotes", quietly = T)) {
    utils::install.packages("remotes")
  }
  if (!requireNamespace("UCell", quietly = T)) {
    remotes::install_github("carmonalab/UCell")
  }
  if (!requireNamespace("scGate", quietly = T)) {
    remotes::install_github("carmonalab/scGate")
  }
  if (!requireNamespace("ProjecTILs", quietly = T)) {
    remotes::install_github("carmonalab/ProjecTILs")
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
  }

  hs_mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  hs_ids <- biomaRt::getBM(attributes = unique(c("hgnc_symbol", "ensembl_gene_id")),
                           filters = "hgnc_symbol",
                           values = features_hs,
                           mart = hs_mart) %>%
    tibble::as_tibble() %>%
    dplyr::rename("Gene.HS" = 1, "Gene.stable.ID.HS" = 2) %>%
    dplyr::distinct() %>%
    dplyr::filter(trimws(Gene.HS) != "", trimws(Gene.stable.ID.HS) != "")

  mm_orth <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name"),
                            filters = "ensembl_gene_id",
                            values = hs_ids[,"Gene.stable.ID.HS",drop=T],
                            mart = hs_mart) %>%
    tibble::as_tibble() %>%
    dplyr::rename("Gene.stable.ID.HS" = 1, "Gene.stable.ID.MM" = 2, "Gene.MM" = 3) %>%
    dplyr::distinct() %>%
    dplyr::filter(trimws(Gene.MM) != "", trimws(Gene.stable.ID.MM) != "")

  # make a long data frame from alt.symbol of ProjecTILs::Hs2Mm.convert.table
  pt_table <-
    ProjecTILs::Hs2Mm.convert.table %>%
    dplyr::mutate(Alt.symbol = stringr::str_split(Alt.symbol, ",")) %>%
    dplyr::mutate(Alt.symbol.HS = stringr::str_split(Alt.symbol.HS, ",")) %>%
    tidyr::unnest(Alt.symbol) %>%
    tidyr::unnest(Alt.symbol.HS) %>%
    dplyr::mutate(Gene.stable.ID.HS = ifelse(Alt.symbol.HS == Gene.HS, Gene.stable.ID.HS, NA)) %>%
    dplyr::select(-Gene.HS, -Gene.MM) %>%
    dplyr::rename("Gene.HS" = Alt.symbol.HS, "Gene.MM" = Alt.symbol) %>%
    tibble::as_tibble()

  # bind rows and make distinct
  ortholog_table <-
    hs_ids %>%
    dplyr::left_join(mm_orth, by = "Gene.stable.ID.HS") %>%
    dplyr::select(-Gene.stable.ID.MM) %>%
    dplyr::bind_rows(pt_table) %>%
    tidyr::drop_na() %>%
    dplyr::distinct(Gene.HS, Gene.MM, .keep_all = T) # keep an arbitrary Gene.stable.ID.HS for each match

  features_hs_wo_ortholog <- features_hs[which(!features_hs %in% ortholog_table$Gene.HS)]
  features_mm_wo_ortholog <- features_mm[which(!features_mm %in% ortholog_table$Gene.MM)]

  return(list(ortholog_table = ortholog_table,
              features_hs_wo_ortholog = features_hs_wo_ortholog,
              features_mm_wo_ortholog = features_mm_wo_ortholog))

}
