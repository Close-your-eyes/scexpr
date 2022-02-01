#library(shiny)
#library(tidyverse)
library(clustree) # necessary:

# hide functions in package
# fix split violin

data <- readRDS(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "data.rds"))

ui <- shiny::fluidPage(
  shiny::fluidRow(
    shiny::column(2,shiny::selectInput(inputId = "pval", label = "p value:", choices = c("adj.p.val", "p.val"), selected = "adj.p.val")),
    shiny::column(2,shiny::numericInput(inputId = "pval.tick", label = "pval tick", min = 1e-300, max = 0.99, value = 0.01)),
    shiny::column(2,shiny::numericInput(inputId = "font.size", label = "font size", value = 20)),
    shiny::column(2,shiny::numericInput(inputId = "pt.size2", label = "dot size", min = 0.01, max = 100, value = 2)),
    shiny::column(2,shiny::numericInput(inputId = "label.size", label = "label sizes", min = 1, max = 100, value = 6)),
    shiny::column(2,shiny::textInput(inputId = "feature.labels", label = "label genes", value = ""))
  ),
  shiny::fluidRow(
    shiny::column(2,shiny::numericInput(inputId = "log2.fc.cut", label = "log2.fc cut", value = NA)),
    shiny::column(2,shiny::numericInput(inputId = "p.cut", label = "p value cut", value = NA)),
    shiny::column(2,shiny::numericInput(inputId = "max.volcano.label", label = "Max labels:", min = 0, max = 1000, value = 20)),
    shiny::column(2,shiny::selectInput(inputId = "pseudo.log", label = "y axis pseudo log", choices = c(T, F), selected = F)),
    shiny::column(2,shiny::numericInput(inputId = "pseudo.log.sigma", label = "y axis pseudo log", min = 0, max = 100, value = 1)),
    shiny::column(2,shiny::numericInput(inputId = "min.pct", label = "min pct", min = 0, max = 1, value = 0))
  ),
  shiny::fluidRow(
    shiny::column(2,shiny::selectInput(inputId = "clustering", label = "clustering", choices = names(data)[which(!names(data) %in% c("clustree_plot", "non.aggr.data"))], selected = names(data)[which(!names(data) %in% c("clustree_plot", "non.aggr.data"))][1])),
    shiny::column(2,shiny::selectInput(inputId = "cc", label = "Clusters", choices = names(data[[1]])[which(names(data[[1]]) != "reduction_plot")], selected = names(data[[1]])[which(names(data[[1]]) != "reduction_plot")][1])), # ## works??
    shiny::column(2,shiny::numericInput(inputId = "plot.width", label = "plot.width", value = 12, min = 1, max = 50)),
    shiny::column(2,shiny::numericInput(inputId = "plot.height", label = "plot.height", value = 8, min = 1, max = 50)),
    shiny::column(2,shiny::actionButton(inputId = "save", label = "save.plot"))
  ),
  shiny::fluidRow(
    shiny::column(8, shiny::plotOutput("volcano_plot", brush = "plot_brush")),
    shiny::column(4, shiny::plotOutput("reduction_plot"))
  ),
  shiny::fluidRow(
    shiny::column(12, shiny::plotOutput("clustree_plot"))
  ),

  shiny::fluidRow(
    shiny::column(10,shiny::plotOutput("jitter_plot")),
    shiny::column(2,
                  shiny::fluidRow(shiny::numericInput(inputId = "n.max", label = "Max plots:", min = 0, max = 48, value = 9)),
                  shiny::fluidRow(shiny::selectInput(inputId = "geom1", label = "geom1", choices = c("jitter", "point"), selected = "jitter")),
                  shiny::fluidRow(shiny::selectInput(inputId = "geom2", label = "geom2", choices = c("violin", "boxplot", "split violin", "none"), selected = "none")),
                  shiny::fluidRow(shiny::numericInput(inputId = "pt.size", label = "dot size", min = 0.01, max = 100, value = 2)),
                  shiny::fluidRow(shiny::selectInput(inputId = "freq.expr", label = "plot freq expr", choices = c(T, F), selected = F)),
                  shiny::fluidRow(shiny::selectInput(inputId = "filter", label = "filter non expr", choices = c(T, F), selected = F)),
    )
  ),

  DT::dataTableOutput("volcano_table")
)

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
      ds$ds <- rbind(data[[input$clustering]][[input$cc]][["data"]][intersect(which(data[[input$clustering]][[input$cc]][["data"]][,paste0("pct.", data[[input$clustering]][[input$cc]][["negative.group.name"]])] >= input$min.pct), which(data[[input$clustering]][[input$cc]][["data"]][,"log2.fc"] < 0)),],
                     data[[input$clustering]][[input$cc]][["data"]][intersect(which(data[[input$clustering]][[input$cc]][["data"]][,paste0("pct.", data[[input$clustering]][[input$cc]][["positive.group.name"]])] >= input$min.pct), which(data[[input$clustering]][[input$cc]][["data"]][,"log2.fc"] > 0)),])
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
    f <- unique(c(ds$ds[input$volcano_table_rows_selected, "Feature", drop = T], feats$bf, feats$text))
    f <- f[which(f %in% ds$ds$Feature)]
    plot <- scexpr:::.plot_vp(vd = ds$ds,
                              p.plot = input$pval,
                              pt.size = input$pt.size2,
                              pt.alpha = 0.8,
                              font.size = input$font.size,
                              pval.tick = input$pval.tick,
                              features.exclude = data[[input$clustering]][[input$cc]][["features.exclude"]],
                              fc.cut = input$log2.fc.cut,
                              p.cut = input$p.cut,
                              ngn = data[[input$clustering]][[input$cc]][["negative.group.name"]],
                              pgn = data[[input$clustering]][[input$cc]][["positive.group.name"]],
                              y.axis.pseudo.log = input$pseudo.log,
                              pseudo.log.sigma = input$pseudo.log.sigma) +
      ggplot2::labs(title = input$clustering)

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
      return(scexpr:::.expr_jitter(d = data[["non.aggr.data"]],
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
    }
  })

  shiny::observeEvent(input$save, {
    ggplot2::ggsave(plot = saveplot, filename = paste0(format(as.POSIXct(Sys.time(), format = "%d-%b-%Y-%H:%M:%S"), "%Y%m%d.%H%M%S"), "_plot.png"), path = dirname(rstudioapi::getActiveDocumentContext()$path), dpi = "retina", width = input$plot.width, height = input$plot.height)
  })
}

shiny::shinyApp(ui = ui, server = server)
