#' Plot composition of one category by another
#'
#' @param SO
#' @param x_cat
#' @param fill_cat
#'
#' @return
#' @export
#'
#' @examples
composition_barplot <- function(SO,
                                x_cat,
                                fill_cat,
                                col_pal = scexpr::col_pal("custom")) {

  if (!x_cat %in% names(SO@meta.data)) {
    stop("x_cat not found in SO@meta.data.")
  }
  if (!fill_cat %in% names(SO@meta.data)) {
    stop("fill_cat not found in SO@meta.data.")
  }

  table <-
    SO@meta.data %>%
    dplyr::count(!!rlang::sym(x_cat), !!rlang::sym(fill_cat)) %>%
    dplyr::left_join(dplyr::count(SO@meta.data, !!rlang::sym(x_cat), name = "total"), by = x_cat) %>%
    dplyr::mutate(rel = n/total) %>%
    tibble::as_tibble()

  plot <- ggplot2::ggplot(table, aes(x = !!rlang::sym(x_cat), y = rel, fill = !!rlang::sym(fill_cat))) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = col_pal)

  return(list(table = table, plot = plot))
}
