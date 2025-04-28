#' Plot results from scexpr::fgsea_on_msigdbr
#'
#' @param data named numeric vector of ranking metric, e.g. S2N; names must
#' be genes; no need to sort, this is done internally; if taken times -1 the
#' gsea plot and ES flips
#' @param ticks_linewidth width of gene ticks
#' @param color_ES_line color of the ES line along the gene ranks
#' @param color_ES_lims color of vertical lines which indicate min and max
#' enrichment score
#' @param color_leadingEdge color of vertical line at the max ES
#' @param label_genes plot gene names next to their ticks, try leadingEdge or
#' a vector gene names to plot
#' @param theme ggplot theme
#' @param plot_leadingEdge_rank plot a vertical line at the rank of max ES
#' @param plot_leadingEdge_size plot the size (n genes) of the leading edge
#' @param annotation_pos vector of x and y value where to plot gsea stats
#' @param annotation_size size of annotation (geom_text size)
#' @param annotation_with_name
#' @param ticks_height length or height of gene ticks
#' @param rect_height height of colored rectangles on top of gene ticks
#' @param rect_alpha alpha level (transparency) of colored rectangles on top of gene ticks
#' @param zscore_lims breaks at which to cut the color scale (which is zscore),
#' if very few ranks fill a zscore range it is pointless to plot such very small
#' rectangle
#'
#' @return
#' @export
#'
#' @examples
plot_gsea <- function(data,
                      ticks_linewidth = 0.2,
                      ticks_height = 1,
                      rect_height = 0.65,
                      rect_alpha = 0.85,
                      color_ES_line = "black",
                      color_ES_lims = "black",
                      color_leadingEdge = "hotpink2",
                      label_genes = NULL, # leave NULL by default because it may take some time in case many gene sets are tested
                      theme = theme_bw(),
                      plot_leadingEdge_rank = T,
                      plot_leadingEdge_size = T,
                      annotation_pos = NULL,
                      annotation_size = 4,
                      zscore_lims = c(-2,2),
                      annotation_with_name = F) {


  data_colorbar <-
    as.data.frame(data$stats) %>%
    dplyr::mutate(zscore = as.vector(round(scale(stat),0))) %>%
    dplyr::mutate(zscore = ifelse(zscore < min(zscore_lims), min(zscore_lims), zscore)) %>%
    dplyr::mutate(zscore = ifelse(zscore > max(zscore_lims), max(zscore_lims), zscore)) %>%
    dplyr::group_by(zscore) %>%
    dplyr::filter(rank %in% c(min(rank), max(rank))) %>%
    dplyr::summarise(min_rank = min(rank), max_rank = max(rank))

  color_scale_limits <- range(round(scale(data$stats$stat),0))
  color_breaks <- min(c(abs(min(data_colorbar$zscore)), max(data_colorbar$zscore)))
  color_breaks <- seq(1,color_breaks,1)

  sES <- data$spreadES/10
  ticks_height <- sES * ticks_height
  rect_height <- sES * rect_height

  p <-
    ggplot2::ggplot(data = data$curve) +
    ggplot2::geom_segment(
      data = data$ticks,
      mapping = ggplot2::aes(
        x = rank,
        xend = rank,
        y = -ticks_height,
        yend = ticks_height,
      ),
      linewidth = ticks_linewidth
    ) +
    ggplot2::geom_rect(
      data = data_colorbar,
      mapping = ggplot2::aes(
        xmin = min_rank,
        xmax = max_rank,
        fill = zscore
      ),
      ymin = -rect_height,
      ymax = rect_height,
      alpha = rect_alpha,
      show.legend = T
    ) +
    ggplot2::scale_fill_stepsn(colors = rev(RColorBrewer::brewer.pal(11,"RdBu")),
                               breaks = data_colorbar$zscore,
                               limits = c(color_scale_limits[1],color_scale_limits[2]),
                               # center color scale at 0 by ensuring equal number of values around 0
                               # using color_breaks as it is now gives darker colors as if it would be e.g. c(min,1,0,-1,max)
                               values = scales::rescale(c(
                                 color_scale_limits[1],
                                 -sort(color_breaks, decreasing = T),
                                 0,
                                 color_breaks,
                                 color_scale_limits[2]
                               )),
                               show.limits = T) +
    ggplot2::geom_line(ggplot2::aes(x=rank, y=ES), color = color_ES_line) +
    ggplot2::geom_hline(yintercept = data$posES, colour = color_ES_lims, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = data$negES, colour = color_ES_lims, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "black") +
    theme +
    ggplot2::labs(x = "gene rank", y = "enrichment score (ES)", fill = "ranking metric\n[z-score]") +
    ggplot2::guides(fill = ggplot2::guide_colorsteps())

  if (!is.null(label_genes)) {
    p <- p + ggrepel::geom_text_repel(data = if (label_genes == "leadingEdge") le_gene_df else all_gene_df[which(all_gene_df$gene %in% label_genes)], ggplot2::aes(label = gene, x = x, y = 0),
                                      max.overlaps = length(le_ranks),
                                      max.time = 5,
                                      nudge_y = data[["spreadES"]]/2.5,
                                      nudge_x = data[["leadingEdge_rank"]]*6,
                                      segment.size = 0.2,
                                      segment.color = "grey80")
  }

  if (plot_leadingEdge_rank) {
    p <- p + ggplot2::geom_segment(x = data[["leadingEdge_rank"]],
                                   xend = data[["leadingEdge_rank"]],
                                   y = 0,
                                   yend = c(data[["posES"]], data[["negES"]])[which.max(abs(c(data[["posES"]], data[["negES"]])))],
                                   color = color_leadingEdge,
                                   linetype = "dashed")
    if (plot_leadingEdge_size) {
      p <- p + ggplot2::annotate("text",
                                 label = paste0("n = ", data[["leadingEdge_size"]]), #  which(data[["ticks"]]$rank == data[["leadingEdge_rank"]])
                                 color = color_leadingEdge,
                                 y = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]),
                                            ticks_height + 0.03,
                                            -ticks_height - 0.03),
                                 x = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]),
                                            data[["leadingEdge_rank"]]/2,
                                            data[["leadingEdge_rank"]] + (nrow(data[["stats"]]) - data[["leadingEdge_rank"]])/2 ))
    }
  }

  p <- p + ggplot2::labs(title = paste0(data[["name"]], " (", nrow(data[["ticks"]]), "/", length(data[["gene.set"]]),  ")"),
                         subtitle = paste0("p = ", data[["pval"]], "\nES = ", signif(data[["ES"]], 2), "\nNES = ", data[["NES"]]))


  if (!is.null(annotation_pos)) {
    # use geom_richtext to avoid annotate("richtext") which would require library(ggtext)
    # in geom_richtext create a data.frame to avoid multi plotting of label for every data point in original data frame of ggplot object

    #p1 <- round(data[["pval"]]/10^floor(log10(abs(data[["pval"]]))), 0)
    #p2 <- floor(log10(abs(data[["pval"]])))
    p_raw <- formatC(data[["pval"]], format = "e", digits = 1)
    p1 <- strsplit(p_raw, "e")[[1]][1]
    p2 <- strsplit(p_raw, "e")[[1]][2]

    if (annotation_pos[1] == "auto") {
      annotation_pos[2] <- data[["ES"]] * 0.8 # y
      if (data[["ES"]] > 0) {
        annotation_pos[1] <- which.min(as.data.frame(data[["stats"]])$stat[which(as.data.frame(data[["stats"]])$stat > 0)]) #x
      } else {
        annotation_pos[1] <- 300
      }
      annotation_pos <- as.numeric(annotation_pos)
    }

    if (annotation_with_name) {
      ann_label <- paste0("**", data[["name"]], " (n = ", nrow(data[["ticks"]]), ")", "**<br>*p* = ", p1, " x 10<sup>", p2, "</sup><br>ES = ", signif(data[["ES"]], 2), "<br>NES = ", data[["NES"]])
    } else {
      ann_label <- paste0("*p* = ", p1, " x 10<sup>", p2, "</sup><br>ES = ", signif(data[["ES"]], 2), "<br>NES = ", data[["NES"]])
    }

    p <- p +
      ggtext::geom_richtext(data = data.frame(label = ann_label),
                            ggplot2::aes(label = label),
                            x = annotation_pos[1],
                            y = annotation_pos[2],
                            size = annotation_size,
                            fill = NA,
                            label.color = NA,
                            hjust = 0)
  }
  p_metric <- ggplot2::ggplot(as.data.frame(data[["stats"]]),
                              ggplot2::aes(x = rank, y = stat)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_vline(xintercept = which.min(as.data.frame(data[["stats"]])$stat[which(as.data.frame(data[["stats"]])$stat > 0)]) + 0.5, linetype = "dashed", linewidth = 0.3) +
    theme +
    ggplot2::labs(x="gene rank", y="ranking metric")
  return(list(plot = p, metric_plot = p_metric))
}
