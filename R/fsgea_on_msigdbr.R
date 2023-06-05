#' Convenient wrapper around fgsea function
#'
#' @param gene.ranks
#' @param gene.sets
#' @param min.padj
#' @param use.msigdbr
#' @param msigdbr_args
#' @param fgsea_fun
#' @param fgsea_args
#' @param return.gene.sets set to FALSE in order to save memory when always the same external set of gene sets is used
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # get gene sets from msigdbr, split to a named list
#' gene.sets <- scexpr:::.get.split.msigdbr()
#' # provide gene.sets manually and leave use.msigdbr = F in case fgsea_on_msigdbr is called many times with the same gene.sets
#' }
fgsea_on_msigdbr <- function(gene.ranks = NULL,
                             gene.sets = NULL,
                             return.gene.sets = T,
                             min.padj = 0.001, # which plots to generate
                             use.msigdbr = F,
                             msigdbr_args = list(species = "Homo sapiens", category = NULL, subcategory = NULL),
                             fgsea_fun = fgsea::fgseaMultilevel,
                             fgsea_args = list(stats = gene.ranks, pathways = gene.sets)) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    BiocManager::install("fgsea")
  }

  if (is.null(gene.ranks)) {
    stop("gene.ranks has to be provided.")
  }
  if (!is.numeric(gene.ranks)) {
    stop("gene.ranks has to be a numeric metric.")
  }
  if (is.null(names(gene.ranks))) {
    stop("gene.ranks has to have names which should be genes.")
  }
  if (!is.null(gene.sets)) {
    if (!is.list(gene.sets) || is.null(names(gene.sets))) {
      stop("gene.sets has to be a named list of gene sets.")
    }
  }
  fgsea_fun <- match.fun(fgsea_fun)

  if (use.msigdbr) {
    if (!requireNamespace("msigdbr", quietly = T)) {
      utils::install.packages("msigdbr")
    }
    gene.sets <- c(gene.sets, .get.split.msigdbr(msigdbr_args = msigdbr_args))
  }
  if (length(gene.sets) == 0) {
    stop("No gene.sets provided and use.msigdbr=F. Either manually provide gene.sets or set use.msigdbr=T to select sets from their.")
  }
  names(gene.sets) <- make.unique(names(gene.sets))


  if (!"pathways" %in% names(fgsea_args)) {
    fgsea_args <- c(list(pathways = gene.sets), fgsea_args)
  }
  if (!"stats" %in% names(fgsea_args)) {
    fgsea_args <- c(list(stats = gene.ranks), fgsea_args)
  }
  results <- as.data.frame(Gmisc::fastDoCall(fgsea_fun, args = fgsea_args))

  gsea_plots <- lapply(.self_naming(results[which(results$padj <= min.padj),"pathway",drop=T]), function(x) {
    fgsea::plotEnrichment(pathway = gene.sets[[x]], stats = gene.ranks) +
      labs(title = x, subtitle = paste0("p = ", signif(results[which(results$pathway == x), "padj"], 2)))
  })

  return(list(fgsea_table = results,
              fgsea_plots = gsea_plots,
              gene.sets = if (return.gene.sets) {gene.sets} else {NULL},
              gene.ranks = gene.ranks))
}


.get.split.msigdbr <- function(msigdbr_args = list(species = "Homo sapiens", category = NULL, subcategory = NULL)) {
  sets <- Gmisc::fastDoCall(msigdbr::msigdbr, args = msigdbr_args)[,c("gs_name", "gene_symbol")]
  sets <- split(sets$gene_symbol, sets$gs_name)
  return(sets)
}

.self_naming <- function(x) {
  return(stats::setNames(x, x))
}

