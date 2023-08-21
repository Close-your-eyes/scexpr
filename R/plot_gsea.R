#' Title
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
                      plot_leadingEdge_size = T) {

  p <- with(data,
            p <- ggplot2::ggplot(data=curve) +
              ggplot2::geom_line(ggplot2::aes(x=rank, y=ES), color = color_ES_line) +
              ggplot2::geom_segment(data=ticks,
                                    mapping=ggplot2::aes(x=rank, y=-spreadES/10,
                                                         xend=rank, yend=spreadES/10),
                                    linewidth=ticksSize) +
              ggplot2::geom_hline(yintercept=posES, colour=color_min_max_ES_line, linetype="dashed") +
              ggplot2::geom_hline(yintercept=negES, colour=color_min_max_ES_line, linetype="dashed") +
              ggplot2::geom_hline(yintercept=0, colour="black") +
              theme +
              ggplot2::labs(x="gene rank", y="enrichment score (ES)")
  )
  p <- p + ggplot2::geom_rect(data = data[["colorbar_df"]], aes(xmin = min_rank,
                                                                xmax = max_rank,
                                                                fill = I(fill_col)),
                              ymin = -data[["spreadES"]]/16,
                              ymax = data[["spreadES"]]/16, alpha = 0.8)


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
                                 y = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), data[["spreadES"]]/8, -data[["spreadES"]]/8),
                                 x = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), data[["leadingEdge_rank"]]/2, data[["leadingEdge_rank"]] + (nrow(data[["stats"]]) - data[["leadingEdge_rank"]])/2 ))
      #hjust = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), 0, 1)
      # data[["leadingEdge_rank"]]*1.1, data[["leadingEdge_rank"]]*0.98
    }
  }
  p <- p + ggplot2::labs(title = paste0(data[["name"]], " (n = ", nrow(data[["ticks"]]), ")"), subtitle = paste0("p = ", data[["pval"]], "\nES = ", signif(data[["ES"]], 2), "\nNES = ", data[["NES"]]))

  p_metric <- ggplot2::ggplot(as.data.frame(data[["stats"]]),
                              ggplot2::aes(x = rank, y = stat)) +
    ggplot2::geom_col() +
    theme +
    ggplot2::labs(x="gene rank", y="ranking metric")
  return(list(plot = p, metric_plot = p_metric))
}
