#' Groupwise gene set enrichment analysis
#'
#' S2N is calculated with scexpr::s2n_groupwise, then fgsea::fgseaMultilevel is
#' called on each group. You may use an overclustered grouping, GSEA on groups of cells
#' should cancel out noise in single cells and group highly correlates cells.
#'
#' @param obj Seurat object
#' @param group grouping column in obj meta.data
#' @param fgseaMultilevel_args arguments to fgsea::fgseaMultilevel; named list of
#' pathways is mandatory
#' @param get_layer_args arguments to scexpr::get_layer
#'
#' @returns list
#' @export
#'
#' @examples
fgsea_groupwise <- function(obj,
                            group,
                            fgseaMultilevel_args = list(pathways = list()),
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

  s2ntab <- s2n_groupwise(
    obj = obj,
    group = group,
    get_layer_args = get_layer_args
  )
  gseares <- purrr::map_dfr(stats::setNames(colnames(s2ntab), colnames(s2ntab)),
                            ~as.data.frame(Gmisc::fastDoCall(fgsea::fgseaMultilevel,
                                                             args = c(list(stats = s2ntab[,.x]), fgseaMultilevel_args))),
                            .id = group)
  # try({
  #   # only works with a single pathway
  #   rownames(gseares) <- gseares[[group]]
  # }, silent = T)


  gseares2 <- split(gseares, gseares$pathway)
  gseares2 <- purrr::map_dfc(names(gseares2), function(i) {
    obj@meta.data |>
      tibble::rownames_to_column("ididid") |>
      dplyr::select(!!rlang::sym(group), ididid) |>
      dplyr::left_join(stats::setNames(dplyr::select(gseares2[[i]], !!rlang::sym(group), "ES", "NES", "padj"),
                                       c(group, paste0(i, c("_ES", "_NES", "_padj")))),
                       by = group) |>
      tibble::column_to_rownames("ididid") |>
      dplyr::select(-!!rlang::sym(group))
  })

  return(list(gsea = gseares, meta = gseares2))
}
