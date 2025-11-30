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
#'\dontrun{
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
#' msigdb <- scexpr::gsea_get_msigdb(return = "df")
#'
#' # numeric summary
#' summ <- msigdb$sets |>
#'   dplyr::summarise(n_genes = dplyr::n(), .by = c(gs_collection, gs_name))
#'
#' ggplot(summ, aes(x = gs_collection, y = n_genes)) +
#'   geom_boxplot(outlier.shape = NA) +
#'   geom_jitter(width = 0.1, size = 0.2)
#'
#' ggplot(summ, aes(x = gs_collection, y = n_genes)) +
#'   geom_boxplot(outlier.shape = NA) +
#'   geom_jitter(width = 0.1, size = 0.2) +
#'   scale_y_log10()
#'
#' quantile(summ$n_genes, probs = seq(0.1,1,0.1))
#'
#' # some checks
#' all_genes <- unique(msigdb$sets$gene_symbol)
#' conv <- convert_gene_ident(x = all_genes,
#'                            input = c("hgnc_symbol", "external_synonym", "external_gene_name", "ensembl_gene_id"))
#' conv_short <- conv[["short"]]
#' sum(is.na(conv_short$hgnc_symbol)) # 2164
#' sum(is.na(conv_short$ensembl_gene_id)) # 0
#' anyDuplicated(conv_short$hgnc_symbol[which(!is.na(conv_short$hgnc_symbol))])
#' check <- conv_short |>
#'   dplyr::add_count(hgnc_symbol) |>
#'   dplyr::filter(!is.na(hgnc_symbol))
#' any(check$n > 1)
#'
#' sum(conv_short$unequal_input_hgnc, na.rm = T) # 158 - very minor
#'
#' # harmonized:
#' # dplyr::mutate(harmonized = ifelse(is.na(hgnc_symbol), input, hgnc_symbol))
#'
#' # cluster gene sets
#' msigdb <- scexpr::gsea_get_msigdb()
#' mat <- brathering::list_to_binary_matrix(msigdb$sets)
#' k <- 50
#' pca_res <- irlba::irlba(mat, nv = k, nu = k)
#' # scores or embedding
#' pc_embedding <- pca_res$u %*% diag(pca_res$d)
#' loadings <- pca_res$v
#'
#' ures <- uwot::umap(pc_embedding, verbose = T)
#' brathering::plot2(ures, size = 2)
#' }
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





