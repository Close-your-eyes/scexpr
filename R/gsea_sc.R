#' Run single-cell gene set enrichment analysis
#'
#' Performs gene set enrichment analysis independently for each cell. Gene
#' expression values are standardized across all cells for each gene, and the
#' resulting per-cell gene rankings are analyzed against the supplied gene sets.
#' Cells are processed in chunks to reduce peak memory usage.
#'
#' @param obj A single-cell object supported by [get_layer()]. Its expression
#'   matrix must have genes in rows and cells in columns.
#' @param gene_sets A named list of gene sets. Each element should be a character
#'   vector of gene identifiers matching the row names of the expression matrix.
#' @param mc.cores A positive integer specifying the number of parallel worker
#'   processes used within each chunk. Defaults to `16`. On Windows,
#'   [parallel::mclapply()] does not support multicore processing and this value
#'   should generally be set to `1`.
#'
#' @return A data frame containing the enrichment results returned by
#'   [gsea_on_msigdbr()] for every cell. Identifier columns indicate the
#'   originating chunk and cell.
#'
#' @details
#' For each gene, the reference mean and standard deviation are calculated
#' across all cells. Each cell's expression profile is then converted to
#' gene-wise z-scores:
#'
#' \deqn{z_{gc} = \frac{x_{gc} - \bar{x}_g}{s_g + 10^{-8}}}
#'
#' The small constant prevents division by zero for genes with no variation.
#' The expression matrix is split into 50 column-wise chunks before parallel
#' processing. Warnings produced during enrichment analysis are suppressed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gene_sets <- list(
#'   interferon_response = c("STAT1", "IRF1", "ISG15"),
#'   cell_cycle = c("MKI67", "CDK1", "CCNB1")
#' )
#'
#' results <- gsea_sc(
#'   obj = obj,
#'   gene_sets = gene_sets,
#'   mc.cores = 4
#' )
#' }
gsea_sc <- function(obj,
                    gene_sets,
                    mc.cores = 16) {

  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }

  if (missing(gene_sets)) {
    stop("gene_sets is missing.")
  }

  if (!is.list(gene_sets) || is.data.frame(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list.")
  }

  ref_mean <- Matrix::rowMeans(get_layer(obj))
  ref_sd   <- brathering::rowsds(get_layer(obj))


  # densify whole matrix: too much memory
  #cell_z <- sweep(sweep(get_layer(obj), 1, ref_mean, FUN = "-"), 1, ref_sd + 1e-8, FUN = "/")

  # creates chunks of cells to pass as a subset to mc cores - lower memory?!
  nchunks <- 50
  chunks <- brathering::split_chunks(1:ncol(get_layer(obj)), chunks = nchunks)
  factor <- rep(names(chunks), lengths(chunks))
  matlst <- brathering::split_mat(get_layer(obj), factor, byrow = F)
  matlst <- matlst[as.character(sort(as.numeric(names(matlst))))]

  out <- purrr::imap(matlst, function(y, idy) {
    print(paste0("chunk ", idy, "/", nchunks))
    # parallel compute over cells
    out <- parallel::mclapply(purrr::set_names(colnames(y)), function(x) {
      cell_z <- (as.numeric(y[, x]) - ref_mean) / (ref_sd + 1e-8)
      suppressWarnings(gsea <- gsea_on_msigdbr(gene_ranks = cell_z, gene_sets = gene_sets,
                                               return_gene_sets_subset = F,
                                               return_gene_sets = F,
                                               return_leading_edge_cols = F)[["data"]][,1:6,drop = F])
      return(gsea)
    }, mc.cores = mc.cores)
    dplyr::bind_rows(out, .id = "id")
  })
  out <- dplyr::bind_rows(out)

  return(out)
}
