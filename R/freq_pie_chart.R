freq_pie_chart <- function(SO,
                           meta.col,
                           inset.text.size = 5,
                           inset.text.radius = 0.75,
                           legend.position = "none") {

  if (!requireNamespace("ggforce", quietly = T)) {
    utils::install.packages("ggforce")
  }
  if (!requireNamespace("farver", quietly = T)) {
    utils::install.packages("farver")
  }

  # https://stackoverflow.com/questions/16184188/ggplot-facet-piechart-placing-text-in-the-middle-of-pie-chart-slices (ggforce)
  tab <- table(SO@meta.data[,meta.col])
  tab <- data.frame(frac = as.numeric(tab/sum(tab)), cluster = factor(names(tab), levels = names(tab)))
  tab <- tab[order(tab$frac, decreasing = T), ]
  tab$start_angle <- c(0,cumsum(tab$frac))[-(length(tab$frac) + 1)]*pi*2
  tab$end_angle <- c(cumsum(tab$frac))*pi*2
  tab$mid_angle <-  0.5*(tab$start_angle + tab$end_angle)

  cols <- col_pal("custom", n = nlevels(as.factor(tab[,"cluster"])))
  ggplot2::ggplot(tab, ggplot2::aes(x0 = 0, y0 = 0, r0 = 0.3, r = 1, start = start_angle, end = end_angle, fill = cluster)) +
    ggforce::geom_arc_bar(colour = "white") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position, panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), text = ggplot2::element_text(family = "mono"), axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::geom_text(ggplot2::aes(color = farver::decode_colour(cols, to = "hcl")[,"l"] > 50,
                                    x = inset.text.radius*sin(mid_angle),
                                    y = inset.text.radius*cos(mid_angle),
                                    label = format(round(frac, 2), nsmall = 2)),
                       size = inset.text.size) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(guide = FALSE, values = c("white", "black")) +
    ggplot2::coord_fixed(ratio = 1)
}
