#' Plot results from scexpr::fgsea_on_msigdbr
#'
#' @param data
#' @param ticksSize
#' @param color_ES_line
#' @param color_min_max_ES_line
#' @param color_leadingEdge
#' @param label_genes
#' @param theme
#' @param plot_leadingEdge_rank
#' @param plot_leadingEdge_size
#' @param annotation_pos
#' @param annotation_size
#' @param annotation_with_name
#'
#' @return
#' @export
#'
#' @examples
plot_gsea <- function(data,
                      ticksSize=0.2,
                      color_ES_line = "black",
                      color_min_max_ES_line = "black",
                      color_leadingEdge = "hotpink2",
                      label_genes = NULL, # leave NULL be default because it may take some time in case many gene sets are tested
                      theme = ggplot2::theme_bw(),
                      plot_leadingEdge_rank = T,
                      plot_leadingEdge_size = T,
                      annotation_pos = NULL,
                      annotation_size = 4,
                      annotation_with_name = F) {

  data$stats$zscore <- round(scale(data$stats$stat),1)
  data2 <-
    as.data.frame(data$stats) %>%
    dplyr::group_by(zscore) %>%
    dplyr::filter(rank %in% c(min(rank), max(rank))) %>%
    dplyr::summarise(min_rank = min(rank), max_rank = max(rank)) %>%
    dplyr::mutate(zscore2 = round(zscore,0)) %>%
    dplyr::mutate(zscore2 = ifelse(zscore2 < -4, -4, zscore2)) %>%
    dplyr::mutate(zscore2 = ifelse(zscore2 > 4, 4, zscore2))
  color_break_limit <- min(c(abs(min(data2$zscore2)), max(data2$zscore2)))-1


  # manually provide bluish and reddish color to scale_fill_stepsn in order to have switch from blue to red at zero
  #cols1 <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))[1:5])(round(sum(data2$zscore<0)/length(data2$zscore)*10))
  #cols2 <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))[7:11])(round(sum(data2$zscore>0)/length(data2$zscore)*10))
  #scales::show_col(c(cols1, cols2))
  #scales::show_col(rev(RColorBrewer::brewer.pal(11,"RdBu")))

  p <-
    ggplot2::ggplot(data = data$curve) +
    ggplot2::geom_segment(data = data$ticks,
                          mapping=ggplot2::aes(x = rank, y = -data$spreadES/10, xend = rank, yend = data$spreadES/10),
                          linewidth = ticksSize) +
    ggplot2::geom_rect(data = data2, aes(xmin = min_rank,
                                         xmax = max_rank,
                                         fill = zscore2),
                       ymin = -data[["spreadES"]]/16,
                       ymax = data[["spreadES"]]/16, alpha = 0.8,
                       show.legend = T) +
    ggplot2::scale_fill_stepsn(colors = rev(RColorBrewer::brewer.pal(11,"RdBu")), breaks = seq(-color_break_limit,color_break_limit,1)) +
    ggplot2::geom_line(ggplot2::aes(x=rank, y=ES), color = color_ES_line) +
    ggplot2::geom_hline(yintercept = data$posES, colour = color_min_max_ES_line, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = data$negES, colour = color_min_max_ES_line, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "black") +
    theme +
    ggplot2::labs(x = "gene rank", y = "enrichment score (ES)", fill = "ranking metric\n[z-score]")

browser()
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
                                   linetype="dashed")
    if (plot_leadingEdge_size) {
      p <- p + ggplot2::annotate("text",
                                 label = paste0("n = ", data[["leadingEdge_size"]]), #  which(data[["ticks"]]$rank == data[["leadingEdge_rank"]])
                                 color = color_leadingEdge,
                                 y = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), data[["spreadES"]]/7, -data[["spreadES"]]/7),
                                 x = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), data[["leadingEdge_rank"]]/2, data[["leadingEdge_rank"]] + (nrow(data[["stats"]]) - data[["leadingEdge_rank"]])/2 ))
    }
  }

  p <- p + ggplot2::labs(title = paste0(data[["name"]], " (n = ", nrow(data[["ticks"]]), ")"),
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
                            aes(label = label),
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
