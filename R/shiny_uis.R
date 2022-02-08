#' Title
#'
#' @param name
#'
#' @return
#' @export
#'
#' @importFrom ggraph guide_edge_colourbar
#'
#' @examples
shiny_uis <- function(name = c("single_volcano_sc",
                               "multi_clustering_volcano_sc")) {

  name <- match.arg(name, c("single_volcano_sc", "multi_clustering_volcano_sc"))

  if (name == "multi_clustering_volcano_sc") {
    ui <- shiny::fluidPage(
      shiny::fluidRow(
        shiny::column(2,shiny::selectInput(inputId = "pval", label = "p value:", choices = c("adj.p.val", "p.val"), selected = "adj.p.val")),
        shiny::column(2,shiny::numericInput(inputId = "pval.tick", label = "pval tick", min = 1e-300, max = 0.99, value = 0.01)),
        shiny::column(2,shiny::numericInput(inputId = "font.size", label = "font size", value = 20)),
        shiny::column(2,shiny::numericInput(inputId = "pt.size2", label = "dot size", min = 0.01, max = 100, value = 0.5)),
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
        shiny::column(2,shiny::selectInput(inputId = "clustering", label = "clustering", choices = names(data)[which(!names(data) %in% c("clustree_plot", "Seurat_object"))], selected = names(data)[which(!names(data) %in% c("clustree_plot", "Seurat_object"))][1])),
        shiny::column(2,shiny::selectInput(inputId = "cc", label = "Clusters", choices = names(data[[1]])[which(names(data[[1]]) != "reduction_plot")], selected = names(data[[1]])[which(names(data[[1]]) != "reduction_plot")][1])),
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
                      shiny::fluidRow(shiny::selectInput(inputId = "geom1", label = "geom1", choices = as.character(formals(scexpr::feature_plot_stat)$geom1)[-1], selected = as.character(formals(scexpr::feature_plot_stat)$geom1)[-1][1])),
                      shiny::fluidRow(shiny::selectInput(inputId = "geom2", label = "geom2", choices = as.character(formals(scexpr::feature_plot_stat)$geom2)[-1], selected = as.character(formals(scexpr::feature_plot_stat)$geom2)[-1][1])),
                      shiny::fluidRow(shiny::numericInput(inputId = "pt.size", label = "dot size", min = 0.01, max = 100, value = 0.5)),
                      shiny::fluidRow(shiny::selectInput(inputId = "freq.expr", label = "plot freq expr", choices = c(T, F), selected = F)),
                      shiny::fluidRow(shiny::selectInput(inputId = "filter", label = "filter non expr", choices = c(T, F), selected = F)),
        )
      ),
      DT::dataTableOutput("volcano_table")
    )
  }

  if (name == "single_volcano_sc") {
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
        shiny::column(2,shiny::selectInput(inputId = "cc", label = "Clusters", choices = names(data[[1]])[which(names(data[[1]]) != "reduction_plot")], selected = names(data[[1]])[which(names(data[[1]]) != "reduction_plot")][1])),
        shiny::column(2,shiny::numericInput(inputId = "plot.width", label = "plot.width", value = 12, min = 1, max = 50)),
        shiny::column(2,shiny::numericInput(inputId = "plot.height", label = "plot.height", value = 8, min = 1, max = 50)),
        shiny::column(2,shiny::actionButton(inputId = "save", label = "save.plot"))
      ),
      shiny::fluidRow(
        shiny::column(8, shiny::plotOutput("volcano_plot", brush = "plot_brush")),
        shiny::column(4, shiny::plotOutput("reduction_plot"))
      ),
      shiny::fluidRow(
        shiny::column(10,shiny::plotOutput("jitter_plot")),
        shiny::column(2,
                      shiny::fluidRow(shiny::numericInput(inputId = "n.max", label = "Max plots:", min = 0, max = 48, value = 9)),
                      shiny::fluidRow(shiny::selectInput(inputId = "geom1", label = "geom1", choices = as.character(formals(scexpr::feature_plot_stat)$geom1)[-1], selected = as.character(formals(scexpr::feature_plot_stat)$geom1)[-1][1])),
                      shiny::fluidRow(shiny::selectInput(inputId = "geom2", label = "geom2", choices = as.character(formals(scexpr::feature_plot_stat)$geom2)[-1], selected = as.character(formals(scexpr::feature_plot_stat)$geom2)[-1][1])),
                      shiny::fluidRow(shiny::numericInput(inputId = "pt.size", label = "dot size", min = 0.01, max = 100, value = 2)),
                      shiny::fluidRow(shiny::selectInput(inputId = "freq.expr", label = "plot freq expr", choices = c(T, F), selected = F)),
                      shiny::fluidRow(shiny::selectInput(inputId = "filter", label = "filter non expr", choices = c(T, F), selected = F)),
        )
      ),
      DT::dataTableOutput("volcano_table")
    )
  }

  return(ui)
}
