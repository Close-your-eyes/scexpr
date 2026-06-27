#' Run gene set enrichment analysis with fgsea and optional MSigDB gene sets
#'
#' Convenient wrapper around an `fgsea` function, with optional retrieval of
#' gene sets from MSigDB via `msigdbr`. The input `gene_ranks` is passed as
#' `stats`, and `gene_sets` is passed as `pathways` to `fgsea_fun`.
#'
#' @param gene_ranks Named numeric vector of gene-level ranking statistics.
#' Names must be gene identifiers matching the identifiers in `gene_sets`.
#'
#' @param gene_sets Named list of gene sets, where each element is a character
#' vector of gene identifiers. Used as `pathways` in `fgsea_fun`.
#'
#' @param return_gene_sets Logical. If `TRUE`, return the full gene set list in
#' the output. Set to `FALSE` to reduce memory use when repeatedly using the
#' same gene sets.
#'
#' @param return_gene_sets_subset Logical. If `TRUE`, return each gene set
#' restricted to genes present in `gene_ranks`.
#'
#' @param use_msigdbr Logical. If `TRUE`, retrieve gene sets from MSigDB using
#' `msigdbr::msigdbr()` and combine them with any manually supplied `gene_sets`.
#'
#' @param msigdbr_args Named list of arguments passed to
#' `msigdbr::msigdbr()`. By default, human Hallmark and C1--C8 collections are
#' used.
#'
#' @param fgsea_fun Function used to run GSEA. Defaults to
#' `fgsea::fgseaMultilevel`.
#'
#' @param fgsea_args Named list of additional arguments passed to `fgsea_fun`.
#' If `pathways` or `stats` are not supplied, they are filled from `gene_sets`
#' and `gene_ranks`, respectively.
#'
#' @param ... Currently unused; reserved for future extensions.
#'
#' @return A named list with:
#' \describe{
#'   \item{data}{A data frame of fgsea results, optionally joined with MSigDB
#'   collection metadata. Additional columns include sorted leading-edge genes,
#'   leading-edge size, relative leading-edge size, and leading-edge rank.}
#'   \item{gene_sets}{The gene sets used for enrichment, or `NULL` when
#'   `return_gene_sets = FALSE`.}
#'   \item{gene_sets.subset}{Gene sets restricted to genes found in
#'   `gene_ranks`, or `NULL` when `return_gene_sets_subset = FALSE`.}
#'   \item{gene_ranks}{The input ranking vector.}
#' }
#'
#' @details
#' `gene_ranks` should usually be sorted in decreasing order before running
#' GSEA, although `fgsea` can operate on named numeric vectors directly.
#'
#' When `use_msigdbr = TRUE`, gene sets are retrieved using `.get.split.msigdbr()`
#' and MSigDB metadata such as collection, subcollection, description, and exact
#' source are merged into the result table.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Compute signal-to-noise gene ranks for one cluster
#' s2n <- scexpr::gsea_s2n_groupwise(
#'   obj = so,
#'   group = "integrated_snn_res.0.7"
#' )
#'
#' gene_ranks <- sort(s2n[, "20"], decreasing = TRUE)
#'
#' ## Use MSigDB gene sets through msigdbr
#' gsea_res <- gsea_on_msigdbr(
#'   gene_ranks = gene_ranks,
#'   use_msigdbr = TRUE
#' )
#'
#' head(gsea_res$data)
#'
#' ## Reuse gene sets across many calls
#' gene_sets <- scexpr:::.get.split.msigdbr()$sets
#'
#' gsea_res <- gsea_on_msigdbr(
#'   gene_ranks = gene_ranks,
#'   gene_sets = gene_sets,
#'   use_msigdbr = FALSE,
#'   return_gene_sets = FALSE
#' )
#'
#' ## Available MSigDB collections
#' msigdbr::msigdbr_collections()
#' # H: hallmark gene sets
#' # C1: positional gene sets
#' # C2: curated gene sets
#' # C3: regulatory target gene sets
#' # C4: computational gene sets
#' # C5: ontology gene sets
#' # C6: oncogenic signature gene sets
#' # C7: immunologic signature gene sets
#' # C8: cell type signature gene sets
#' }
gsea_on_msigdbr <- function(gene_ranks,
                            gene_sets = NULL,
                            return_gene_sets_subset = T,
                            return_gene_sets = T,
                            use_msigdbr = F,
                            msigdbr_args = list(
                              db_species = "HS",
                              species = "human",
                              collection = c(
                                "C1",
                                "C2",
                                "C3",
                                "C4",
                                "C5",
                                "C6",
                                "C7",
                                "C8",
                                "H"
                              )
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
    results <- dplyr::left_join(results, gene_sets.list[["sets_cats"]], by = c("pathway" = "gs_name"))
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
    Gmisc::fastDoCall(msigdbr::msigdbr, args = msigdbr_args)[,c("gs_collection", "gs_subcollection", "gs_collection_name", "gs_description", "gs_name", "gs_exact_source", "gene_symbol")]
  })

  sets <- unique(sets) # some gene symbols exist in duplicates, but ENSEMBLE differs; make gene symbols unique
  sets_cat <- sets |> dplyr::select(-gene_symbol) |> dplyr::distinct()
  return(list(sets = split(sets$gene_symbol, sets$gs_name), sets_cats = sets_cat))
}


