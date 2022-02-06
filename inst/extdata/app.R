library(ggraph)


data <- readRDS(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "data.rds"))

server <- function(input, output, session) {

  shiny::observeEvent(input$clustering, {
    cluster.names <- names(data[[input$clustering]])
    cluster.names <- cluster.names[which(cluster.names != "reduction_plot")]
    shiny::updateSelectInput(session, "cc", choices = cluster.names, selected = cluster.names[1])
  })

  ds <- shiny::reactiveValues()
  shiny::observeEvent({
    input$min.pct
    input$clustering
    input$cc}, {
      ds$ds <- as.data.frame(Seurat::Misc(data[[input$clustering]][[input$cc]])[["volcano.data"]])
      ds$ds$Feature <- rownames(ds$ds)
      #df <- as.data.frame(data[[input$clustering]][[input$cc]][["data"]])
      #ds$ds <- rbind(df[intersect(which(df[,paste0("pct.", data[[input$clustering]][[input$cc]][["negative.group.name"]])] >= input$min.pct), which(df[,"log2.fc"] < 0)),],df[intersect(which(df[,paste0("pct.", data[[input$clustering]][[input$cc]][["positive.group.name"]])] >= input$min.pct), which(df[,"log2.fc"] > 0)),])
      #ds$ds$Feature <- rownames(ds$ds)
    })

  feats <- shiny::reactiveValues()

  shiny::observeEvent({
    input$plot_brush
    input$pval
    input$pseudo.log
    input$pseudo.log.sigma}, {
      if (input$pseudo.log) {
        trafo <- scales::pseudo_log_trans(base = 10, sigma = input$pseudo.log.sigma)
        tr.fun <- function(x) {trafo$transform(round(scales::oob_squish_infinite(-log10(x), range = c(0,300)), 2))} # same as in .vp_plot
      } else {
        tr.fun <- function(x) {round(scales::oob_squish_infinite(-log10(x), range = c(0,300)), 2)} # same as in .vp_plot
      }
      feats$bf <- ds$ds[Reduce(intersect, list(which(tr.fun(ds$ds[,input$pval,drop=T]) >= input$plot_brush$ymin),
                                               which(tr.fun(ds$ds[,input$pval,drop=T]) <= input$plot_brush$ymax),
                                               which(ds$ds[,"log2.fc",drop=T] >= input$plot_brush$xmin),
                                               which(ds$ds[,"log2.fc",drop=T] <= input$plot_brush$xmax))), "Feature", drop = T]
    })

  shiny::observeEvent(input$feature.labels, {
    feats$text <- gsub("\"", "", unlist(strsplit(unlist(strsplit(input$feature.labels, ",")), " ")))
  })

  output$volcano_table = DT::renderDataTable(dplyr::mutate(ds$ds[,c("Feature", "log2.fc", "p.val", "adj.p.val")], log2.fc = format(round(log2.fc,2), nsmall = 2), p.val = signif(p.val, 3), adj.p.val = signif(adj.p.val, 3)), server = T)

  output$clustree_plot = shiny::renderPlot({
    return(data[["clustree_plot"]])
  })

  output$reduction_plot = shiny::renderPlot({
    return(data[[input$clustering]][["reduction_plot"]])
  })

  output$volcano_plot = shiny::renderPlot({
    plot <- scexpr:::.plot_vp(vd = ds$ds,
                              p.plot = input$pval,
                              pt.size = input$pt.size2,
                              font.size = input$font.size,
                              pval.tick = input$pval.tick,
                              fc.cut = input$log2.fc.cut,
                              p.cut = input$p.cut,
                              y.axis.pseudo.log = input$pseudo.log,
                              pseudo.log.sigma = input$pseudo.log.sigma,
                              features.exclude = Seurat::Misc(data[[input$clustering]][[input$cc]])[["features.exclude"]],
                              ngn = Seurat::Misc(data[[input$clustering]][[input$cc]])[["ngn"]],
                              pgn = Seurat::Misc(data[[input$clustering]][[input$cc]])[["pgn"]]) +
      ggplot2::labs(title = input$clustering)

    f <- unique(c(ds$ds[input$volcano_table_rows_selected, "Feature", drop = T], feats$bf, feats$text))
    f <- f[which(f %in% ds$ds$Feature)]
    plot <- plot + ggplot2::geom_point(color = "black", data = ds$ds[f,,drop=F])
    if (length(f) <= input$max.volcano.label) {
      plot <- plot + ggrepel::geom_text_repel(color = "red", max.overlaps = 500, size = input$label.size, data = ds$ds[f,,drop=F], fontface = "italic", family = "Courier")
    }
    saveplot <<- plot
    return(plot)
  })

  output$jitter_plot <- shiny::renderPlot({
    f <- unique(c(ds$ds[input$volcano_table_rows_selected, "Feature", drop = T], feats$bf, feats$text))
    f <- f[which(f %in% ds$ds$Feature)]
    if (!is.null(f) && length(f) <= input$n.max && length(f) > 0) {

      plots <- feature_plot_stat(SO = data[[input$clustering]][[input$cc]],
                        pt.size = input$pt.size,
                        geom2 = input$geom2,
                        #font.size = input$font.size,
                        geom1 = input$geom1,
                        filter.non.expr = input$filter,
                        plot.expr.freq = input$freq.expr,
                        label.size = input$label.size,
                        features = f,
                        meta.col = Seurat::Misc(data[[input$clustering]][[input$cc]])[["volcano.group.col"]],
                        assay = Seurat::DefaultAssay(data[[input$clustering]][[input$cc]]),
                        strip.background = element_rect(fill = "white"),
                        axis.title = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text = element_text(face = "italic"))

'      plots <- Seurat::VlnPlot(object = data[[input$clustering]][[input$cc]],
                              features = f,
                              group.by = Seurat::Misc(data[[input$clustering]][[input$cc]])[["volcano.group.col"]],
                              assay = Seurat::DefaultAssay(data[[input$clustering]][[input$cc]]),
                              cols = col_pal("custom"),
                              combine = F)
      plots <- lapply(plots, function(x) x +ggplot2::theme_bw() +
               ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "none"))
      return(cowplot::plot_grid(plotlist = plots))'
      return(plots)
    }


    ' if (!is.null(f) && length(f) <= input$n.max && length(f) > 0) {
      if ("non.aggr.data" %in% names(data)) {
        d <- data[["non.aggr.data"]]
      } else if ("non.aggr.data" %in% names(data[[input$clustering]][[input$cc]])) {
        d <- data[[input$clustering]][[input$cc]][["non.aggr.data"]]
      } else {
        stop("non.aggr.data not found.")
      }
      return(scexpr:::.expr_jitter(d = d,
                                   feat = f,
                                   pt.size = input$pt.size,
                                   geom2 = input$geom2,
                                   font.size = input$font.size,
                                   geom1 = input$geom1,
                                   filter.non.expr = input$filter,
                                   plot.expr.freq = input$freq.expr,
                                   label.size = input$label.size,
                                   ngc = data[[input$clustering]][[input$cc]][["negative.group.cells"]],
                                   pgc = data[[input$clustering]][[input$cc]][["positive.group.cells"]],
                                   ngn = data[[input$clustering]][[input$cc]][["negative.group.name"]],
                                   pgn = data[[input$clustering]][[input$cc]][["positive.group.name"]]))
    }'

  })

  shiny::observeEvent(input$save, {
    ggplot2::ggsave(plot = saveplot, filename = paste0(format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d.%H%M%S"), "_plot.png"), path = dirname(rstudioapi::getActiveDocumentContext()$path), dpi = "retina", width = input$plot.width, height = input$plot.height)
  })

}

# shiny_uis("multi_clustering_volcano_sc)
shiny::shinyApp(ui = scexpr::shiny_uis("single_volcano_sc"), server = server)
shiny::shinyApp(ui = scexpr::shiny_uis("multi_clustering_volcano_sc"), server = server)

