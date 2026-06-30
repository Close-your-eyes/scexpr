#' Plot summary statistics of many enriched gene sets, i.e. enrichplot
#'
#' @param gsea_df (subset) data frame returned from gsea_on_msigdbr
#' @param x variable/column on x axis
#' @param y variable/column on y axis
#' @param fill variable/column for fill
#' @param size variable/column for dot size
#' @param geom which geom to plot
#' @param add_size_to_y add size of gene set to y axis text
#' @param minus_log10_padj calculate -log10 of adjusted p-value
#' @param order_y_by which column to order y axis breaks by
#' @param size_range dot size range from to
#' @param x_expand_mult x axis expansion (mult) to accommodate dots fully
#' @param theme ggplot theme
#' @param col_pal color palette for continuous fill
#' @param title plot title
#' @param plot_order plot dots in order of which parameter? avoid burying.
#'
#' @returns ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' gsea <- scexpr::gsea_on_msigdbr(gene_ranks = gene_ranks,
#'                                 use_msigdbr = T,
#'                                 fgsea_args = list(minSize = 10),
#'                                 msigdbr_args = list(db_species = "HS",
#'                                                     species = "human",
#'                                                     collection = "C5"))
#' gsea_df <- gsea[["data"]]
#' gsea_plot_summary(gsea_df = gsea_df |>
#'                     dplyr::filter(NES>0) |>
#'                     dplyr::slice_min(padj, n = 10))
#' }
gsea_plot_summary <- function(gsea_df,
                              x = c("NES", "padj", "leadingEdge_size", "leadingEdge_size_rel", "ES"),
                              y = "pathway",
                              fill = c("padj", "NES", "leadingEdge_size", "leadingEdge_size_rel", "ES"),
                              size = c("leadingEdge_size_rel", "leadingEdge_size", "padj", "NES", "ES"),
                              geom = c("point", "col"),
                              add_size_to_y = T,
                              minus_log10_padj = T,
                              order_y_by = x,
                              size_range = c(3,10),
                              x_expand_mult = c(0.1, 0.1),
                              theme = colrr::theme_material(white = T),
                              col_pal = colrr::col_pal(name = "spectral", direction = -1),
                              title = NULL,
                              plot_order = c("NES", "leadingEdge_size_rel", "ES", "padj")) {

  x <- rlang::arg_match(x)
  y <- rlang::arg_match(y)
  fill <- rlang::arg_match(fill)
  size <- rlang::arg_match(size)
  geom <- rlang::arg_match(geom)
  plot_order <- rlang::arg_match(plot_order)

  gsea_df <- tidyr::drop_na(gsea_df)

  if (add_size_to_y) {
    gsea_df <- dplyr::mutate(gsea_df, !!y := paste0(!!rlang::sym(y), " (", size, ")"))
  }
  if (minus_log10_padj) {
    gsea_df <- dplyr::mutate(gsea_df, padj = -log10(padj))
  }
  if (!is.null(order_y_by)) {
    gsea_df <- dplyr::mutate(gsea_df, !!y := forcats::fct_reorder(!!rlang::sym(y), !!rlang::sym(order_y_by)))
  }

  title_x <- x
  title_size <- size
  title_fill <- fill

  if (x == "padj" && minus_log10_padj) {
    title_x <- "-log10(padj)"
  }
  if (size == "padj" && minus_log10_padj) {
    title_size <- "-log10(padj)"
  }
  if (fill == "padj" && minus_log10_padj) {
    title_fill <- "-log10(padj)"
  }
  if (x %in% c("leadingEdge_size_rel", "leadingEdge_size")) {
    title_x <- gsub("_", " ", x)
  }
  if (size %in% c("leadingEdge_size_rel", "leadingEdge_size")) {
    title_size <- gsub("_", "\n", size)
  }
  if (fill %in% c("leadingEdge_size_rel", "leadingEdge_size")) {
    title_fill <- gsub("_", "\n", fill)
  }

  breaks_size <- ggplot2::waiver()
  if (size == "leadingEdge_size_rel") {
    breaks_size <- sort(unique(round(gsea_df[[size]], 1)))
  }


gsea_df <- dplyr::arrange(gsea_df, -!!rlang::sym(plot_order))
  p <- ggplot2::ggplot(gsea_df, ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y))) +
    theme +
    ggplot2::scale_fill_gradientn(colors = col_pal) +
    ggplot2::scale_size_continuous(range = size_range, breaks = breaks_size) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = x_expand_mult), name = title_x) +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "black"),
                                                 title = title_size,
                                                 order = 2),
                    fill = ggplot2::guide_colorbar(barheight = 6, barwidth = 1,
                                                   title = title_fill,
                                                   order = 1)) +
    ggplot2::labs(subtitle = title)

  if (geom == "point") {
    p <- p + ggplot2::geom_point(ggplot2::aes(
      fill = !!rlang::sym(fill),
      size = !!rlang::sym(size)
    ), shape = 21)
  }
  if (geom == "col") {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = !!rlang::sym(fill)), color = "black")
  }

  return(p)
}


#' Call gsea_plot_summary with pre-defined groups
#'
#' Goal is to get a structured overview.
#'
#' @param gsea_df data frame returned from gsea_on_msigdbr
#' @param topn how many gene sets per group, no ties
#' @param metric what to sort by for topn
#' @param join_NES_plus_minus enriched and de-enriched sets are plotted
#' separately by default; join instead
#' @param split which columns to split by
#' @param sep separator for split labels
#' @param ... args to gsea_plot_summary
#'
#' @returns list of filtered split data frames and ggplots
#' @export
#'
#' @examples
#' \dontrun{
#' gsea <- gsea_on_msigdbr(gene_ranks = gene_ranks,
#'                         use_msigdbr = T,
#'                         fgsea_args = list(minSize = 10),
#'                         msigdbr_args = list(db_species = "HS",
#'                                             species = "human",
#'                                             collection = "C5"))
#' gsea_df <- gsea[["data"]]
#' grouped_plots <- gsea_plot_summary_default_groups(gsea_df)
#' grouped_plots[["plots"]][["C5_GO:BP___plus"]]
#'
#' # plot graph of GO gene sets
#' goid <- grouped_plots[["dfs"]][["C5_GO:BP___plus"]]$gs_exact_source
#' plot_GO_graph(goid = goid)
#' plot_GO_graph(goid = goid,
#'               layout = "dendrogram",
#'               multiline = "long")
#' }
gsea_plot_summary_default_groups <- function(gsea_df,
                                             topn = 10,
                                             metric = c("padj", "NES", "leadingEdge_size", "leadingEdge_size_rel", "ES"),
                                             join_NES_plus_minus = F,
                                             split = c("gs_collection", "gs_subcollection"),
                                             sep = "_",
                                             ...) {

  metric <- rlang::arg_match(metric)

  if (any(!split %in% names(gsea_df))) {
    stop("not all split found in names of gsea_df.")
  }
  gsea_df$split <- do.call(paste, c(as.list(gsea_df[,split]), sep = sep))
  gsea_df2 <- split(gsea_df, gsea_df[["split"]])
  gsea_df2 <- purrr::map(gsea_df2, ~split(.x, ifelse(.x[["NES"]]>0, "__plus", "__minus")))
  gsea_df2 <- purrr::list_flatten(gsea_df2)
  if (metric == "padj") {
    slicefun <- dplyr::slice_min
  } else {
    slicefun <- dplyr::slice_max
  }
  gsea_df2 <- purrr::map(gsea_df2, ~slicefun(.x, order_by = abs(!!rlang::sym(metric)), n = topn, with_ties = F))

  if (join_NES_plus_minus) {
    gsea_df2 <- brathering::list_unflatten(gsea_df2, sep = "___")
    gsea_df2 <- purrr::map(gsea_df2, dplyr::bind_rows)
  }

  gsea_plots <- purrr::map2(gsea_df2, names(gsea_df2), function(x,y) gsea_plot_summary(gsea_df = x, title = y, ...))
  return(list(plots = gsea_plots, dfs = gsea_df2))
}
