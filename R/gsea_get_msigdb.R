#' Convenience function to get gene sets from msigdbr
#'
#' @param collection vector of collections
#' @param msigdbr_args arguments to msigdbr::msigdbr
#' @param return what to return; df for manual inspect, list for passing to gsea
#' functions
#' @param return_df when return df, reduce columns to gs_collection, gs_name,
#' gene_symbol?
#'
#' @returns list of gene sets and data frame with category info
#' @export
#'
#' @examples
#' # H: hallmark gene sets
#' # C1: positional gene sets
#' # C2: curated gene sets
#' # C3: regulatory target gene sets
#' # C4: computational gene sets
#' # C5: ontology gene sets
#' # C6: oncogenic signature gene sets
#' # C7: immunologic signature gene sets
#' # C8: cell type signature gene sets
#' sets <- gsea_get_msigdb(collection = c("C1", "H"))
#' # then pass on
#' gsea_groupwise(obj = SO,
#'                 group = "meta.col",
#'                 fgseaMultilevel_args = list(pathways = list(sets$sets)))
gsea_get_msigdb <- function(collection = c("C1","C2","C3","C4","C5","C6","C7","C8","H"),
                       msigdbr_args = list(
                         db_species = "HS",
                         species = "human"),
                       return = c("list", "df"),
                       return_df = c("reduced", "full")) {

  return <- rlang::arg_match(return)
  return_df <- rlang::arg_match(return_df)
  collection <- rlang::arg_match(collection, multiple = T)

  sets <- purrr::map_dfr(collection, function(x) {
    y <- Gmisc::fastDoCall(what = msigdbr::msigdbr,
                           args = c(list(collection = x), msigdbr_args))
    if (return_df == "reduced") {
      # some gene symbols exist in duplicates, but ENSEMBLE differs; make gene symbols unique
      y <- dplyr::distinct(y, gs_collection, gs_name, gene_symbol)
    }
    return(y)
  })

  cats <- dplyr::distinct(sets, gs_collection, gs_name)
  if (return == "list") {
    sets = split(sets$gene_symbol, sets$gs_name)
  }

  return(list(sets = sets,
              collection = cats))
}

