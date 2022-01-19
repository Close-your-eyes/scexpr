fgsea_on_msigdbr <- function(gene.ranks,
                             use.msigdbr = T,
                             gene.sets = NULL,
                             min.padj = 0.001, # which plots to generate
                             ...) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    BiocManager::install("fgsea")
  }

  dots <- list(...)
  if (any(c("pathways", "stats") %in% names(dots))) {
    stop("Please do not provide pathways or stats in ...")
  }

  if (!is.null(gene.sets)) {
    if (!is.list(gene.sets) || is.null(names(list))) {
      stop("gene.sets has to be a named list of gene sets. See examples as returned from msigdbr::msigdbr splitted to list by gene_symbol and gs_name.")
    }
  }

  if (use.msigdbr) {
    sets <- as.data.frame(do.call(msigdbr::msigdbr, args = dots[which(names(dots) %in% names(formals(msigdbr::msigdbr)))])[,c("gs_name", "gene_symbol")])
    sets <- split(sets$gene_symbol, sets$gs_name)
  }
  if (!is.null(gene.sets)) {
    sets <- c(gene.sets, sets)
    names(sets) <- make.unique(names(sets))
  }

  results <- do.call(fgsea::fgseaMultilevel, args = c(unlist(dots[which(names(dots) %in% names(formals(fgsea::fgseaMultilevel)))]),
                                                      list(pathways = sets,
                                                           stats = gene.ranks)))
  results_signif <- results[which(results$padf <= min.padj)]
  gsea_plots <- lapply(results_signif$pathway, function(x) {
    fgsea::plotEnrichment(pathway = sets[[x]], stats = gene.ranks)
  })
  names(gsea_plots) <- results_signif$pathway

  return(list(fgsea_table = results,
              fgsea_plots = gsea_plots,
              gene.sets = sets,
              gene.ranks = gene.ranks))
}
