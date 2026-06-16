#' Groupwise gene set enrichment analysis
#'
#' S2N is calculated with scexpr::gsea_s2n_groupwise, then fgsea::fgseaMultilevel is
#' called on each group. You may use an overclustered grouping, GSEA on groups of cells
#' should cancel out noise in single cells and group highly correlates cells.
#'
#' @param obj Seurat object or gene x cell matrix
#' @param group grouping column in obj meta.data or vector of length ncol(obj)
#' @param fgseaMultilevel_args arguments to fgsea::fgseaMultilevel; named list of
#' pathways is mandatory
#' @param get_layer_args arguments to scexpr::get_layer
#'
#' @returns list
#' @export
#'
#' @examples
gsea_groupwise <- function(obj,
                           group,
                           fgseaMultilevel_args = list(pathways = list(), nproc = 0),
                           get_layer_args = list()) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    BiocManager::install("fgsea")
  }
  if (!"pathways" %in% names(fgseaMultilevel_args)) {
    stop("pathways in fgseaMultilevel_args missing.")
  }
  if (!length(fgseaMultilevel_args$pathways) || !is.list(fgseaMultilevel_args$pathways) || is.null(names(fgseaMultilevel_args$pathways))) {
    stop("pathways in fgseaMultilevel_args requires a named list of gene sets.")
  }

  s2ntab <- gsea_s2n_groupwise(
    obj = obj,
    group = group,
    get_layer_args = get_layer_args
  )

  if (length(group) > 1) {
    # group is vector and obj is matrix, not seurat
    group <- "group"
  }

  # only calc once when 2 columns (groups) only
  gseares <- purrr::map_dfr(stats::setNames(colnames(s2ntab), colnames(s2ntab)),
                            ~as.data.frame(Gmisc::fastDoCall(fgsea::fgseaMultilevel,
                                                             args = c(list(stats = s2ntab[,.x]), fgseaMultilevel_args))),
                            .id = group)
  gseares$leadingEdge_sorted_chr <- purrr::map_chr(purrr::map(gseares$leadingEdge, sort), paste, collapse = ",")
  gseares$leadingEdge_size <- lengths(gseares$leadingEdge)
  gseares$leadingEdge_size_rel <- gseares$leadingEdge_size/gseares$size
  # gseares$leadingEdge_rank <- purrr::map_int(gseares$pathway, function(x) ifelse(gseares[["ES"]][[which(gseares[["pathway"]] == x)]] > 0,
  #                                                                                max(match(gseares[["leadingEdge"]][[which(gseares[["pathway"]] == x)]], names(rev(gene_ranks)))),
  #                                                                                min(match(gseares[["leadingEdge"]][[which(gseares[["pathway"]] == x)]], names(rev(gene_ranks))))))

  # try({
  #   # only works with a single pathway
  #   rownames(gseares) <- gseares[[group]]
  # }, silent = T)

  # gseares2 <- split(gseares, gseares$pathway)
  # gseares2 <- purrr::map_dfc(names(gseares2), function(i) {
  #   obj@meta.data |>
  #     tibble::rownames_to_column("ididid") |>
  #     dplyr::select(!!rlang::sym(group), ididid) |>
  #     dplyr::left_join(stats::setNames(dplyr::select(gseares2[[i]], !!rlang::sym(group), "ES", "NES", "padj"),
  #                                      c(group, paste0(i, c("_ES", "_NES", "_padj")))),
  #                      by = group) |>
  #     tibble::column_to_rownames("ididid") |>
  #     dplyr::select(-!!rlang::sym(group))
  # })

  gseares <- split(gseares, gseares[[1]])
  gseares <- purrr::map(gseares, ~list(data = .x,
                                       gene_sets = fgseaMultilevel_args[["pathways"]],
                                       gene_sets.subset = purrr::map(fgseaMultilevel_args[["pathways"]], ~.x[.x %in% rownames(s2ntab)]),
                                       gene_ranks = s2ntab[,unique(.x[[1]])]))

  return(gseares)
}
