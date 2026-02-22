#' Check and visualize ties in ranked gene list
#'
#' @param x numeric vector to check for ties
#'
#' @returns list with table, data.frame and ggplots
#' @export
#'
#' @examples
gsea_inspect_ranked_list <- function(x) {
  tab <- table(x)
  tabstack <- utils::stack(tab) |>
    dplyr::mutate(ind_ind = as.numeric(ind), ind_num = as.numeric(as.character(ind)))
  n_tied_elements <- sum(tab[tab > 1])
  rel_tied_elements <- n_tied_elements/length(x)

  p <- ggplot2::ggplot(tabstack, aes(ind_ind, values)) +
    ggplot2::geom_col(color = "white") +
    colrr::theme_material() +
    ggplot2::labs(x = "index", y = "count", subtitle = paste0(round(rel_tied_elements*100, 2), " % tied"))

  p2 <- ggplot2::ggplot(tabstack, aes(ind_num, values)) +
    ggplot2::geom_col(color = "white") +
    colrr::theme_material() +
    ggplot2::labs(y = "count", subtitle = paste0(round(rel_tied_elements*100, 2), " % tied"))


  return(list(tab = tab, df = tabstack, tie_freq = rel_tied_elements, plot1 = p, plot2 = p2))
}
