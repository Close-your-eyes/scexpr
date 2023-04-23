#' Title
#'
#' @param SO Seurat object or data frame
#' @param meta.col
#' @param inset.text.size
#' @param inset.text.radius
#' @param legend.position
#' @param col_pal
#' @param order_pieces
#' @param avoid_label_overlap strategy to elegantly avoid label overlap
#' @param ... arguments to ggplot2::theme()
#' @param border_color color of borders of pie pieces
#' @param outside_radius when avoid_label_overlap is outside; the radius where to place label;
#' 1 roughly equals the outer radius of the circle
#'
#' @return
#' @export
#'
#' @examples
freq_pie_chart <- function(SO,
                           meta.col,
                           inset.text.size = 5,
                           inset.text.radius = 0.75,
                           legend.position = "right",
                           col_pal = scexpr::col_pal("custom"),
                           border_color = "white",
                           order_pieces = T,
                           avoid_label_overlap = c("alternating_shift", "outside"),
                           outside_radius = 1.1,
                           #avoid_overlap_min_frac = 0.05,
                           ...) {

  if (!requireNamespace("ggforce", quietly = T)) {
    utils::install.packages("ggforce")
  }
  if (!requireNamespace("farver", quietly = T)) {
    utils::install.packages("farver")
  }

  if (methods::is(SO, "Seurat")) {
    SO <- SO@meta.data
  }

  ## ... make that a list and check for default elements (see below)
  ## them call theme with Gmisc::fastdocall

  avoid_label_overlap <- match.arg(avoid_label_overlap, c("alternating_shift", "outside"))

  # https://stackoverflow.com/questions/16184188/ggplot-facet-piechart-placing-text-in-the-middle-of-pie-chart-slices (ggforce)
  tab <- table(SO[,meta.col], exclude = c())
  tab <- data.frame(frac = as.numeric(tab/sum(tab)), cluster = factor(names(tab), levels = names(tab)))
  tab <- cbind(tab, utils::stack(table(SO[,meta.col], exclude = c())))
  tab <- tab[-4]
  names(tab)[3] <- "abs"
  if (order_pieces) {
    tab <- tab[order(tab$frac, decreasing = T), ]
  }
  tab$start_angle <- c(0,cumsum(tab$frac))[-(length(tab$frac) + 1)]*pi*2
  tab$end_angle <- c(cumsum(tab$frac))*pi*2
  tab$mid_angle <-  0.5*(tab$start_angle + tab$end_angle)

  # text angle equal to angle of circle but readable
  tab$text_angle <- ifelse(tab$mid_angle > pi, 270 - tab$mid_angle*180/pi, 90 - tab$mid_angle*180/pi)

  # optional: adjust position of text labels
  tab$text_radius <- inset.text.radius
  tab$frac_lag <- lag(tab$frac, default = 0)
  tab$frac_lag_diff <- tab$frac - tab$frac_lag

  ## not used yet
  tab$frac_lag_diff_series <- cumsum(abs(tab$frac_lag_diff) <= 0.05)

  rel_series <- rle(abs(tab$frac_lag_diff) <= 0.05)
  if (avoid_label_overlap == "alternating_shift") {
    tab$text_radius[lag(cumsum(rel_series$lengths)+1)[-1]:cumsum(rel_series$lengths)[-1]] <-
      ifelse(lag(cumsum(rel_series$lengths)+1)[-1]:cumsum(rel_series$lengths)[-1] %% 2 == 0,
             tab$text_radius[lag(cumsum(rel_series$lengths)+1)[-1]:cumsum(rel_series$lengths)[-1]] - 0.1,
             tab$text_radius[lag(cumsum(rel_series$lengths)+1)[-1]:cumsum(rel_series$lengths)[-1]] + 0.1)
  } else if (avoid_label_overlap == "outside") {
    tab$text_radius[lag(cumsum(rel_series$lengths)+1)[-1]:cumsum(rel_series$lengths)[-1]] <- outside_radius
  }


  if (length(col_pal) != length(unique(tab[,"cluster"]))) {
    if (is.null(names(col_pal))) {
      if (length(col_pal) < length(unique(tab[,"cluster"]))) {
        col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
        warning("Number of colors provided not sufficient for number of factor levels. Falling back to scales::hue_pal().")
      } else {
        col_pal <- col_pal[1:length(unique(tab[,"cluster"]))]
      }
    } else {
      if (length(col_pal) > length(unique(tab[,"cluster"])) && all(names(col_pal) %in% unique(tab[,"cluster"]))) {
        col_pal <- col_pal[unique(tab[,"cluster"])]
      } else {
        warning("Number of colors provided not matching the number of factor levels in meta.col. Falling back to scales::hue_pal().")
        col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
      }
    }
  } else {
    if (!is.null(names(col_pal)) && !all(names(col_pal) %in% unique(tab[,"cluster"]))) {
      warning("Not all names of col_pal found in factor levels of meta.col. Falling back to scales::hue_pal().")
      col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
    }
  }
  plot <- ggplot2::ggplot(tab, ggplot2::aes(x0 = 0, y0 = 0, r0 = 0.3, r = 1, start = start_angle, end = end_angle, fill = cluster)) +
    ggforce::geom_arc_bar(colour = border_color) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = legend.position,
                   panel.grid = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   ...) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::coord_fixed(ratio = 1)


  # in case label are outside of the circle, adjust their color (black or white) to the panel.background
  tab$text_color <- ifelse(farver::decode_colour(col_pal, to = "hcl")[, "l"] > 50, "black", "white")
  tab[which(tab$text_radius >= 1), "text_color"] <- ifelse(farver::decode_colour(plot[["theme"]][["panel.background"]][["fill"]], to = "hcl")[,"l"] > 50, "black", "white")

  plot <-
    plot +
    ggplot2::geom_text(data = tab, ggplot2::aes(color = I(text_color),
                                                x = text_radius*sin(mid_angle),
                                                y = text_radius*cos(mid_angle),
                                                angle = text_angle,
                                                label = format(round(frac, 2), nsmall = 2)),
                       size = inset.text.size)

  return(list(plot = plot, data = tab))
}
