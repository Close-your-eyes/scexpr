#' Plot composition of one category by another
#'
#' @param SO Seurat object or data frame
#' @param x_cat
#' @param fill_cat
#' @param border_color
#' @param labels
#' @param min_label
#' @param in_percent
#' @param label_decimals
#' @param col_pal
#'
#' @return
#' @export
#'
#' @examples
composition_barplot <- function(SO,
                                x_cat,
                                fill_cat,
                                col_pal = scexpr::col_pal("custom"),
                                border_color = "black",
                                labels = F,
                                min_label = 0.05,
                                in_percent = T,
                                label_decimals = 2) {


  if (methods::is(SO, "Seurat")) {
    SO <- SO@meta.data
  }

  if (!x_cat %in% names(SO)) {
    stop("x_cat not found in SO.")
  }
  if (!fill_cat %in% names(SO)) {
    stop("fill_cat not found in SO.")
  }

  if (in_percent) {
    fctr <- 100
  } else {
    fctr <- 1
  }

  table <-
    SO %>%
    dplyr::count(!!rlang::sym(x_cat), !!rlang::sym(fill_cat)) %>%
    dplyr::left_join(dplyr::count(SO, !!rlang::sym(x_cat), name = "x_total"), by = x_cat) %>%
    dplyr::left_join(dplyr::count(SO, !!rlang::sym(fill_cat), name = "fill_total"), by = fill_cat) %>%
    dplyr::mutate(rel_x = n/x_total) %>%
    dplyr::mutate(rel_fill = n/fill_total)


  ## use cumsum approach to define label position in middle of bars
  ## this allow to remove label by min_label and still have the correct coordinate
  ## coordinates became wrong with position = ggplot2::position_stack(vjust = 0.5) when some values were removed
  label_table <-
    table %>%
    dplyr::group_by(!!rlang::sym(x_cat)) %>%
    dplyr::arrange(desc(!!rlang::sym(fill_cat)), .by_group = T) %>%
    dplyr::mutate(rel_x_cumsum = cumsum(rel_x)) %>%
    dplyr::mutate(rel_x_cumsum_lag = dplyr::lag(rel_x_cumsum, default = 0)) %>%
    dplyr::mutate(label_ypos = rel_x_cumsum_lag + (rel_x_cumsum-rel_x_cumsum_lag)/2) %>%
    dplyr::filter(rel_x >= min_label) %>%
    dplyr::mutate(rel_x = rel_x*fctr) %>%
    dplyr::mutate(label_ypos = label_ypos*fctr) %>%
    tibble::as_tibble()

  table <-
    table %>%
    dplyr::mutate(rel_x = rel_x*fctr) %>%
    tibble::as_tibble()

  if (!is.null(names(col_pal)) && all(!names(col_pal) %in% unique(table[,fill_cat,drop=F]))) {
    col_pal <- unname(col_pal)
  }

  plot <-
    ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(x_cat), y = rel_x, fill = !!rlang::sym(fill_cat))) +
    ggplot2::geom_col(color = border_color) +
    ggplot2::scale_fill_manual(values = col_pal)


  if (labels) {
    plot <- plot + ggplot2::geom_text(data = label_table, ggplot2::aes(label = round(rel_x, label_decimals), x = !!rlang::sym(x_cat), y = label_ypos))
    #position = ggplot2::position_stack(vjust = 0.5))
  }

  return(list(table = table, plot = plot))
}
