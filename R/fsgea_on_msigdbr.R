#' Title
#'
#' @param gene.ranks
#' @param gene.sets
#' @param min.padj
#' @param use.msigdbr
#' @param msigdbr_args
#' @param fgsea_fun
#' @param fgsea_args
#'
#' @return
#' @export
#'
#' @examples
fgsea_on_msigdbr <- function(gene.ranks = NULL,
                             gene.sets = NULL,
                             min.padj = 0.001, # which plots to generate
                             use.msigdbr = T,
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
              gene.sets = gene.sets,
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

