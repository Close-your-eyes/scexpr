#' Plot results from scexpr::gsea_on_msigdbr
#'
#' @param data list as returned from gsea_on_msigdbr; filter rows of data$data
#' to select gene sets for plotting
#' @param padj_min filter data$data by this p-value? set NULL for no filtering
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
gsea_plot <- function(data,
                      padj_min = 0.01,
                      return_metric_plot = F,
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

  results <- data$data
  if (!is.null(padj_min)) {
    results <- results |>
      dplyr::filter(padj <= padj_min)
  }

  if (nrow(results) > 100) {
    choi <- menu(c("yes", "no"), title = paste0(nrow(results), " gene sets to plot?"))
    if (choi == 2) {
      return(NULL)
    }
  }

  gsea_plots <- lapply(stats::setNames(results$pathway, results$pathway), function(x) {
    gsea_plot_data <- prep_gsea(gene.set = data$gene_sets[[x]],
                                gene_ranks = data$gene_ranks,
                                gseaParam = 1,
                                pval = signif(results[which(results$pathway == x), "padj"], 2),
                                NES = signif(results[which(results$pathway == x), "NES"], 2))
    gsea_plot_data$name <- x
    # leadingEdge_rank differs by 1 between methods ?!
    # leadingEdge_rank = results[which(results$pathway == x), "leadingEdge_rank"]

    return(make_gsea_plot(data = gsea_plot_data,
                          ticks_linewidth = ticks_linewidth,
                          ticks_height = ticks_height,
                          rect_height = rect_height,
                          rect_alpha = rect_alpha,
                          color_ES_line = color_ES_line,
                          color_ES_lims = color_ES_lims,
                          color_leadingEdge = color_leadingEdge,
                          label_genes = label_genes,
                          theme = theme,
                          plot_leadingEdge_rank = plot_leadingEdge_rank,
                          plot_leadingEdge_size = plot_leadingEdge_size,
                          annotation_pos = annotation_pos,
                          annotation_size = annotation_size,
                          zscore_lims = zscore_lims,
                          annotation_with_name = annotation_with_name))
  })

  if (!return_metric_plot) {
    gsea_plots <- purrr::map(gsea_plots, `[[`, 1)
  }

  return(gsea_plots)

}


make_gsea_plot <- function(data,
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

  p <- ggplot2::ggplot(data = data$curve) +
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





prep_gsea <- function(gene.set,
                      gene_ranks,
                      gseaParam = 1,
                      pval,
                      NES) {

  data <- .plotEnrichmentData(pathway = gene.set,
                              stats = gene_ranks,
                              gseaParam = gseaParam)

  #colorbar_df <- .prep_gsea_colorbar(x = as.data.frame(data[["stats"]]))
  rank_df <- data.frame(gene = names(rev(gene_ranks)[data[["ticks"]][["rank"]]]), rank = data[["ticks"]][["rank"]])

  leadingEdge_rank <- as.numeric(data[["curve"]][which(data[["curve"]]$ES == c(data[["posES"]], data[["negES"]])[which.max(abs(c(data[["posES"]], data[["negES"]])))]), "rank"])
  if (abs(data[["posES"]]) > abs(data[["negES"]])) {
    rank_df$leadingEdge <- rank_df$rank <= leadingEdge_rank
  } else if (abs(data[["posES"]]) < abs(data[["negES"]])) {
    rank_df$leadingEdge <- rank_df$rank >= leadingEdge_rank
  } else {
    stop("posES equal to negES. What now?")
  }
  leadingEdge_size = sum(rank_df$leadingEdge)



  data <- c(data, list(gene.set = gene.set,
                       rank_df = rank_df,
                       leadingEdge_rank = leadingEdge_rank,
                       leadingEdge_size = leadingEdge_size,
                       pval = pval,
                       ES = c(data[["posES"]], data[["negES"]])[which.max(abs(c(data[["posES"]], data[["negES"]])))],
                       NES = NES))
  return(data)
}

.prep_gsea_colorbar <- function(x) {
  zscore_cuts <- NULL
  x$stat_scale <- scale(x$stat)
  if (is.null(zscore_cuts)) {
    zscore_cuts <- seq(floor(min(x$stat_scale)), ceiling(max(x$stat_scale)), 1)


    unique(c(seq(min(zscore_cuts), max(zscore_cuts), 4), max(zscore_cuts)))

    x <- unique(c(min(zscore_cuts), seq(min(zscore_cuts), 0, 2), 0))
    y <- unique(c(min(zscore_cuts), seq(min(x), 0, 3), 0))

    z <- unique(c(seq(0, max(zscore_cuts), 2), max(zscore_cuts)))
    z1 <- unique(c(0, seq(0, max(z), 3), max(z)))

    split_func(zscore_cuts[which(zscore_cuts <= 0)],4)
  }
  data <-
    x %>%
    dplyr::mutate(zscore = cut(stat_scale, breaks = zscore_cuts)) %>% # as.numeric(as.factor())
    dplyr::group_by(zscore) %>%
    dplyr::summarise(min_rank = min(rank), max_rank = max(rank))
  data$zscore <- factor(data$zscore, levels = data %>% dplyr::arrange(min_rank) %>% dplyr::pull(zscore))
  data <- dplyr::arrange(data, zscore)
  data$fill_col <- ""

  # manually take care, that positive stats become reddish color and negative bluish
  pos_col <- RColorBrewer::brewer.pal(10, "RdBu")[1:5]
  neg_col <- RColorBrewer::brewer.pal(10, "RdBu")[6:10]
  data[which(grepl("^\\(-", data$zscore)),"fill_col"] <- grDevices::colorRampPalette(neg_col, interpolate = "linear")(length(which(grepl("^\\(-", data$zscore))))
  data[which(grepl("^\\([[:digit:]]", data$zscore)),"fill_col"] <- grDevices::colorRampPalette(pos_col, interpolate = "linear")(length(which(grepl("^\\([[:digit:]]", data$zscore))))


  return(list(data = data, breaks = zscore_cuts))
}

.plotEnrichmentData <- function(pathway,
                                stats,
                                gseaParam=1) {
  # copied from fgsea pkg

  if (any(!is.finite(stats))){
    stop("Not all stats values are finite numbers")
  }

  rnk <- rank(-stats)
  ord <- order(rnk)

  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)

  pathway <- unname(as.vector(stats::na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  pathway <- unique(pathway)

  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.table::data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
  ticks <- data.table::data.table(rank=pathway, stat=statsAdj[pathway])
  stats <- data.table::data.table(rank=seq_along(stats), stat=statsAdj)

  res <- list(
    curve=toPlot,
    ticks=ticks,
    stats=stats,
    posES=max(tops),
    negES=min(bottoms),
    spreadES=max(tops)-min(bottoms),
    maxAbsStat=max(abs(statsAdj)))
}

split_func <- function(x, by) {
  r <- diff(range(x))
  out <- seq(0, r - by - 1, by = by)
  c(round(min(x) + c(0, out - 0.51 + (max(x) - max(out)) / 2), 0), max(x))
}
