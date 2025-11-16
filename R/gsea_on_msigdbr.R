#' Convenient wrapper around fgsea function
#'
#' @param gene_ranks this becomes stats in fgsea_fun
#' @param gene_sets this becomes pathways in fgsea_fun
#' @param use_msigdbr
#' @param msigdbr_args
#' @param fgsea_fun
#' @param fgsea_args
#' @param return_gene_sets set to FALSE in order to save memory when the same gene sets are used multiple times
#' @param return_gene_sets_subset return the subset of each gene set that is actually found in gene_ranks
#' @param ... arguments to .plotEnrichment
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #' # H: hallmark gene sets
#' # C1: positional gene sets
#' # C2: curated gene sets
#' # C3: regulatory target gene sets
#' # C4: computational gene sets
#' # C5: ontology gene sets
#' # C6: oncogenic signature gene sets
#' # C7: immunologic signature gene sets
#' # C8: cell type signature gene sets
#' s2ntab <- scexpr::gsea_s2n_groupwise(obj = so, group = "integrated_snn_res.0.7")
#' msiggsea <- gsea_on_msigdbr(gene_ranks = s2ntab[,"20"], use_msigdbr = T)
#' # get gene sets from msigdbr, split to a named list
#' gene_sets <- scexpr:::.get.split.msigdbr()
#' msigdbr::msigdbr_collections()
#' # provide gene_sets manually and leave use_msigdbr = F in case gsea_on_msigdbr is called many times with the same gene_sets
#' }
gsea_on_msigdbr <- function(gene_ranks,
                            gene_sets = NULL,
                            return_gene_sets_subset = T,
                            return_gene_sets = T,
                            use_msigdbr = F,
                            msigdbr_args = list(
                              db_species = "HS",
                              species = "human",
                              collection = c("C1","C2","C3","C4","C5","C6","C7","C8","H")
                            ),
                            fgsea_fun = fgsea::fgseaMultilevel,
                            fgsea_args = list(),
                            ...) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    BiocManager::install("fgsea")
  }
  if (missing(gene_ranks)) {
    stop("gene_ranks has to be provided.")
  }
  if (!is.numeric(gene_ranks)) {
    stop("gene_ranks has to be a numeric metric.")
  }
  if (is.null(names(gene_ranks))) {
    stop("gene_ranks has to have names which should be genes.")
  }
  if (!is.null(gene_sets) && (!is.list(gene_sets) || is.null(names(gene_sets)))) {
      stop("gene_sets has to be a named list of gene sets.")
  }
  fgsea_fun <- match.fun(fgsea_fun)

  if (use_msigdbr) {
    if (!requireNamespace("msigdbr", quietly = T)) {
      utils::install.packages("msigdbr")
    }
    gene_sets.list <- c(gene_sets, .get.split.msigdbr(msigdbr_args = msigdbr_args))
    gene_sets <- gene_sets.list[["sets"]]
  }
  if (length(gene_sets) == 0) {
    stop("No gene_sets provided and use_msigdbr=F. Either manually provide gene_sets or set use_msigdbr=T to select sets from their.")
  }
  names(gene_sets) <- make.unique(names(gene_sets))

  if (!"pathways" %in% names(fgsea_args)) {
    fgsea_args <- c(list(pathways = gene_sets), fgsea_args)
  }
  if (!"stats" %in% names(fgsea_args)) {
    fgsea_args <- c(list(stats = gene_ranks), fgsea_args)
  }

  results <- as.data.frame(Gmisc::fastDoCall(fgsea_fun, args = fgsea_args))
  if (use_msigdbr) {
    results$pathway_cat <- gene_sets.list[["cats"]][results$pathway]
  }

  results$leadingEdge_sorted_chr <- purrr::map_chr(purrr::map(results$leadingEdge, sort), paste, collapse = ",")
  results$leadingEdge_size <- lengths(results$leadingEdge)
  results$leadingEdge_size_rel <- results$leadingEdge_size/results$size
  results$leadingEdge_rank <- purrr::map_int(results$pathway, function(x) ifelse(results[["ES"]][[which(results[["pathway"]] == x)]] > 0,
                                                                                 max(match(results[["leadingEdge"]][[which(results[["pathway"]] == x)]], names(rev(gene_ranks)))),
                                                                                 min(match(results[["leadingEdge"]][[which(results[["pathway"]] == x)]], names(rev(gene_ranks))))))

  return(list(data = results,
              gene_sets = if (return_gene_sets) {gene_sets} else {NULL},
              gene_sets.subset = if(return_gene_sets_subset) {lapply(gene_sets, function(x) intersect(x, names(gene_ranks)))} else {NULL},
              gene_ranks = gene_ranks))
}




.get.split.msigdbr <- function(msigdbr_args = list(
  db_species = "HS",
  species = "human",
  collection = NULL
)) {

  sets <- purrr::map_dfr(msigdbr_args$collection, function(x) {
    msigdbr_args$collection <- x
    Gmisc::fastDoCall(msigdbr::msigdbr, args = msigdbr_args)[,c("gs_collection", "gs_name", "gene_symbol")]
  })

  sets <- unique(sets) # some gene symbols exist in duplicates, but ENSEMBLE differs; make gene symbols unique
  sets_cat <- unique(sets[,which(names(sets) %in% c("gs_collection", "gs_name"))])
  return(list(sets = split(sets$gene_symbol, sets$gs_name), cats = stats::setNames(sets_cat$gs_collection, sets_cat$gs_name)))
}


