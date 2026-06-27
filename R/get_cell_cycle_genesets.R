#' Get cell cycle markers from Seurat and Whitfield and scran
#'
#' @returns list of different cell cycle associated gene sets
#' @export
#'
#' @examples
#' ccgenelsz <- get_cell_cycle_genesets()
get_cell_cycle_genesets <- function() {

  # seur <- c(Seurat::cc.genes, Seurat::cc.genes.updated.2019)
  # names(seur) <- c("seurat_legacy_S", "seurat_legacy_G2M", "seurat_2019_S", "seurat_2019_G2M")
  #
  # msgidb <- scexpr::gsea_get_msigdb(collection = "C2")
  # whits <- grep("WHITFIELD", names(msgidb$sets), value = T)
  # whit_sets <- msgidb$sets[whits]
  # names(whit_sets) <- gsub("G1_S", "G1S", names(whit_sets))
  # names(whit_sets) <- gsub("G2_M", "G2M", names(whit_sets))
  # names(whit_sets) <- gsub("M_G1", "MG1", names(whit_sets))
  #
  # seur_whit <- c(seur, whit_sets)
  # seur_whit_df <- stack(seur_whit) |>
  #   dplyr::rename("feature" = values, "gene_set" = ind)
  # vroom::vroom_write(seur_whit_df, file = "/Users/chris/Documents/R_packages/scexpr/inst/extdata/human_cell_cycle.tsv")
  #

  df <- vroom::vroom(system.file("extdata", "human_cell_cycle.tsv", package = "scexpr"), show_col_types = F, progress = F)
  df$gene_set <- tolower(df$gene_set)
  lst <- split(df$feature, df$gene_set)

  # hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
  # mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  # hsconvg1first <- scexpr::convert_gene_identifier(idents = hs.pairs$G1$first,
  #                                                  ident_in = "ENSEMBL",
  #                                                  ident_out = "SYMBOL") |>
  #   dplyr::select(ENSEMBL, SYMBOL) |>
  #   dplyr::rename("first" = ENSEMBL, "first2" = SYMBOL)
  # hsconvg1sec <- scexpr::convert_gene_identifier(idents = hs.pairs$G1$second,
  #                                                ident_in = "ENSEMBL",
  #                                                ident_out = "SYMBOL") |>
  #   dplyr::select(ENSEMBL, SYMBOL) |>
  #   dplyr::rename("second" = ENSEMBL, "second2" = SYMBOL)
  # hs.pairs$G1 <- hs.pairs$G1 |>
  #   dplyr::left_join(hsconvg1first) |>
  #   dplyr::left_join(hsconvg1sec)
  # hsg1 <- hs.pairs[["G1"]] |>
  #   tidyr::drop_na() |>
  #   dplyr::select(first2, second2) |>
  #   dplyr::rename("first" = first2, "second" = second2)
  #
  # hsconvSfirst <- scexpr::convert_gene_identifier(idents = hs.pairs$S$first,
  #                                                  ident_in = "ENSEMBL",
  #                                                  ident_out = "SYMBOL") |>
  #   dplyr::select(ENSEMBL, SYMBOL) |>
  #   dplyr::rename("first" = ENSEMBL, "first2" = SYMBOL)
  # hsconvSsec <- scexpr::convert_gene_identifier(idents = hs.pairs$S$second,
  #                                                ident_in = "ENSEMBL",
  #                                                ident_out = "SYMBOL") |>
  #   dplyr::select(ENSEMBL, SYMBOL) |>
  #   dplyr::rename("second" = ENSEMBL, "second2" = SYMBOL)
  # hs.pairs$S <- hs.pairs$S |>
  #   dplyr::left_join(hsconvSfirst) |>
  #   dplyr::left_join(hsconvSsec)
  # hsS <- hs.pairs[["S"]] |>
  #   tidyr::drop_na() |>
  #   dplyr::select(first2, second2) |>
  #   dplyr::rename("first" = first2, "second" = second2)
  #
  # hsconvg2mfirst <- scexpr::convert_gene_identifier(idents = hs.pairs$G2M$first,
  #                                                  ident_in = "ENSEMBL",
  #                                                  ident_out = "SYMBOL") |>
  #   dplyr::select(ENSEMBL, SYMBOL) |>
  #   dplyr::rename("first" = ENSEMBL, "first2" = SYMBOL)
  # hsconvg2msec <- scexpr::convert_gene_identifier(idents = hs.pairs$G2M$second,
  #                                                ident_in = "ENSEMBL",
  #                                                ident_out = "SYMBOL") |>
  #   dplyr::select(ENSEMBL, SYMBOL) |>
  #   dplyr::rename("second" = ENSEMBL, "second2" = SYMBOL)
  # hs.pairs$G2M <- hs.pairs$G2M |>
  #   dplyr::left_join(hsconvg2mfirst) |>
  #   dplyr::left_join(hsconvg2msec)
  # hsg2m <- hs.pairs[["G2M"]] |>
  #   tidyr::drop_na() |>
  #   dplyr::select(first2, second2) |>
  #   dplyr::rename("first" = first2, "second" = second2)
  # scran_hs_pairs_cellcycle <- list(G2M = hsg2m,
  #                                  G1 = hsg1,
  #                                  S = hsS)
  # saveRDS(scran_hs_pairs_cellcycle, file = "/Users/chris/Documents/R_packages/scexpr/inst/extdata/scran_hs_pairs_cellcycle.rds")

  scran_hs_pairs_cellcycle <- readRDS("/Users/chris/Documents/R_packages/scexpr/inst/extdata/scran_hs_pairs_cellcycle.rds")

  return(list(cc_lst = lst, scran_hs_pairs_cellcycle = scran_hs_pairs_cellcycle))
}
