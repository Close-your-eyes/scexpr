library(shiny)
library(tidyverse)

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
    shiny::column(2,shiny::selectInput(inputId = "clustering", label = "clustering", choices = names(data), selected = names(data)[1])),
    shiny::column(2,shiny::selectInput(inputId = "cc", label = "Clusters", choices = names(data[[1]]), selected = names(data[[1]])[1]))
  ),
  shiny::plotOutput("volcano_plot", brush = "plot_brush"),
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
    cluster.names <- names(data[[input$clustering]][[input$cc]][[input$clustering]])
    shiny::updateSelectInput(session, "cc", choices = cluster.names, selected = cluster.names[1])
  })

  ds <- reactiveValues()
  shiny::observeEvent({
    input$min.pct
    input$clustering
    input$cc}, {
      ds$ds <- rbind(data[[input$clustering]][[input$cc]][["data"]][intersect(which(data[[input$clustering]][[input$cc]][["data"]][,paste0("pct.", data[[input$clustering]][[input$cc]][["negative.group.name"]])] >= input$min.pct), which(data[[input$clustering]][[input$cc]][["data"]][,"log2.fc"] < 0)),],
                     data[[input$clustering]][[input$cc]][["data"]][intersect(which(data[[input$clustering]][[input$cc]][["data"]][,paste0("pct.", data[[input$clustering]][[input$cc]][["positive.group.name"]])] >= input$min.pct), which(data[[input$clustering]][[input$cc]][["data"]][,"log2.fc"] > 0)),])
    })

  feats <- reactiveValues()

  shiny::observeEvent({
    input$plot_brush
    input$pval
    input$pseudo.log
    input$pseudo.log.sigma}, {


      if (input$pseudo.log) {
        trafo <- scales::pseudo_log_trans(base = 10, sigma = input$pseudo.log.sigma)
        tr.fun <- function(x) {trafo$transform(-log10(x))}
      } else {
        tr.fun <- function(x) {-log10(x)}
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


  output$volcano_plot = shiny::renderPlot({
    f <- unique(c(ds$ds[input$volcano_table_rows_selected, "Feature", drop = T], feats$bf, feats$text))
    f <- f[which(f %in% ds$ds$Feature)]
    plot <- .plot_vp(vd = ds$ds,
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
                     pseudo.log.sigma = input$pseudo.log.sigma)


    plot <- plot + geom_point(color = "black", data = ds$ds[f,,drop=F])
    if (length(f) <= input$max.volcano.label) {
      plot <- plot + geom_text_repel(color = "red", max.overlaps = 500, size = input$label.size, data = ds$ds[f,,drop=F])
    }
    return(plot)
  })

  output$jitter_plot <- shiny::renderPlot({
    f <- unique(c(ds$ds[input$volcano_table_rows_selected, "Feature", drop = T], feats$bf, feats$text))
    f <- f[which(f %in% ds$ds$Feature)]
    if (!is.null(f) && length(f) <= input$n.max && length(f) > 0) {
      return(.expr_jitter(d = data[[input$clustering]][[input$cc]][["non.aggr.data"]],
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
}

shiny::shinyApp(ui = ui, server = server)



