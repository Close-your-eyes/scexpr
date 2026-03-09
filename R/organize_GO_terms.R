#' Get parents of GO.db id
#'
#' @param goid GO id(s)
#'
#' @returns data frame
#' @export
#'
#' @examples
#' \dontrun{
#' get_GO_parents(c("GO:0006084", "GO:0015800"))
#' c5 <- msigdbr::msigdbr(collection = "C5", subcollection = "BP")
#' go <- unique(c5$gs_exact_source)
#' get_GO_parents(go[1])
#' #' }
get_GO_parents <- function(goid) {
  y <- AnnotationDbi::mget(unique(goid), GO.db::GOBPPARENTS, ifnotfound=NA)
  y <- purrr::discard(y, anyNA)

  if (!length(y)) {
    message("no GO parents.")
    return(NULL)
  }
  y <- purrr::map_dfr(y, utils::stack, .id = "child") |>
    dplyr::rename("parent" = values)

  y$child_term <- get_GO_term(y$child)
  y$parent_term <- get_GO_term(y$parent)

  return(y)
}


#' Get parents of GO.db id recursively
#'
#' Build a graph from it.
#'
#' @param goid GO id vector
#' @param seen used by recursion, leave empty
#'
#' @returns from-to-dataframe
#' @export
#'
#' @examples
#' \dontrun{
#' c5 <- gsea_get_msigdb("C5")
#' c5 <- msigdbr::msigdbr(collection = "C5", subcollection = "BP")
#' go <- unique(c5$gs_exact_source)
#' goterm_fromtodf <- get_GO_parents_recursive(go[1:3])
#' plot_GO_graph(fromtodf = goterm_fromtodf)
#' }
get_GO_parents_recursive <- function(goid, seen = character()) {

  res <- get_GO_parents(goid)
  if (is.null(res) || nrow(res) == 0) {
    return(dplyr::tibble(
      child = character(),
      parent = character(),
      ind = character(),
      child_term = character(),
      parent_term = character()
    ))
  }
  # avoid revisiting already processed terms
  new_parents <- setdiff(res$parent, seen)

  if (length(new_parents) == 0) {
    return(res)
  }
  # recursive call
  dplyr::bind_rows(
    res,
    get_GO_parents_recursive(new_parents, seen = c(seen, new_parents))
  )
}



#' Meta data from GO sets
#'
#' @param goid GO id(s)
#'
#' @returns list, one entry for each goid
#' @export
#'
#' @examples
#' \dontrun{
#' get_GO_meta(c("GO:0044281", "GO:0006082", "GO:0006520", "GO:0032787", "GO:0044283",
#'               "GO:0006629", "GO:0016054", "GO:0044282", "GO:0006631", "GO:0016053"))
#' }
get_GO_meta <- function(goid) {
  AnnotationDbi::mget(goid, GO.db::GOTERM, ifnotfound=NA)
}


#' Term or readable name of GO pathway/process/set
#'
#' @param goid GO id(s)
#'
#' @returns named vector
#' @export
#'
#' @examples
#' \dontrun{
#' get_GO_term(c("GO:0044281", "GO:0006082", "GO:0006520", "GO:0032787", "GO:0044283",
#'               "GO:0006629", "GO:0016054", "GO:0044282", "GO:0006631", "GO:0016053"))
#' }
get_GO_term <- function(goid) {
  purrr::map_chr(get_GO_meta(goid), ~.x@Term)
}



#' Plot graph/network of GO terms
#'
#' Filtered for isa-entries
#'
#' @param goid GO id(s)
#' @param fromtodf get_GO_parents_recursive(goid) instead of goid; modify optionally
#' @param layout layout of ggraph,
#' e.g. auto, tree, dendrogram, circlepack, htree, cactustree, treemap, partition,
#' fr, kk, drl, graphopt, lgl, stress, circle, linear, focus, centrality,
#' grid, star, sphere, random, hive, sugiyama
#' @param label_repel repel labels?
#' @param label_multiline wrap labels to multiple lines
#' @param label_goid label GO id instead of term?
#' @param return_graph return graph instead of ggplot flr own plotting
#' @param label_plot plot labels?
#' @param multiline split at every space (long) or at second space only (wide)
#' @param ... args to ggraph::geom_node_label
#' @param expand_x mult x axis expansion
#' @param expand_y mult y axis expansion
#'
#' @returns ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' c5 <- msigdbr::msigdbr(collection = "C5", subcollection = "BP")
#' go <- unique(c5$gs_exact_source)
#' plot_GO_graph(go, label_plot = F) # whole network
#' plot_GO_graph(go[1])
#' plot_GO_graph(go[1:2])
#' }
plot_GO_graph <- function(goid = "",
                          fromtodf = NULL,
                          layout = "tree",
                          label_plot = T,
                          label_repel = F,
                          label_multiline = T,
                          multiline = c("wide", "long"),
                          label_goid = F,
                          return_graph = F,
                          expand_x = c(0.3,0.3),
                          expand_y = c(0.1,0.1),
                          label_split = " |_|-",
                          ...) {

  multiline <- rlang::arg_match(multiline)

  if (is.null(fromtodf)) {
    fromtodf <- get_GO_parents_recursive(goid)
  }
  if (nrow(fromtodf) == 0) {
    message("fromtodf empty. check get_GO_parents_recursive(goid).")
    return(NULL)
  }

  fromtodf_isa <- fromtodf |>
    dplyr::filter(ind == "isa" & parent != "all") |>
    dplyr::distinct() |>
    dplyr::mutate(root = child %in% goid) |>
    dplyr::relocate(child, .after = parent)
  graph <- igraph::graph_from_data_frame(d = fromtodf_isa, directed = T)

  igraph::V(graph)$name2 <- get_GO_term(igraph::V(graph)$name)
  igraph::V(graph)$name3 <- gsub(label_split, "\n", igraph::V(graph)$name2)
  igraph::V(graph)$name4 <- purrr::map_chr(brathering::strsplit2(igraph::V(graph)$name2, " ", 2), ~paste(.x, collapse = "\n"))
  igraph::V(graph)$enriched_set <- igraph::V(graph)$name %in% goid

  if (return_graph) {
    return(graph)
  }

  which_name <- ifelse(multiline == "wide", "name4", "name3")
  label <- ifelse(label_goid, "name", ifelse(label_multiline, which_name, "name2"))

  plot <- custom_ggraph_plot(graph = graph,
                             layout = layout,
                             expand_x = expand_x,
                             expand_y = expand_y,
                             label_plot = label_plot,
                             label_repel = label_repel,
                             label = label,
                             ...)

  return(plot)

}

custom_ggraph_plot <- function(graph,
                               layout = "dendrogram",
                               expand_x = 0.05,
                               expand_y = 0.05,
                               label_plot = T,
                               label_repel = F,
                               label = "name3",
                               ...) {
  plot <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link() +
    ggraph::geom_node_point(size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = expand_x)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = expand_y))

  if (label_plot) {
    plot <- plot +
      ggraph::geom_node_label(ggplot2::aes(label = !!rlang::sym(label),
                                           color = enriched_set),
                              repel = label_repel,
                              label.padding = ggplot2::unit(2, "pt"),
                              lineheight = 0.75,
                              ...)
  }
  return(plot)
}

set_root_and_add_rootdist <- function(graph) {
  root <- which(igraph::degree(graph, mode = "in") == 0)
  igraph::V(graph)$root <- FALSE
  igraph::V(graph)$root[root] <- TRUE
  dist_root <- igraph::distances(graph, v = root)
  igraph::V(graph)$dist_root <- as.numeric(dist_root)
  return(graph)
}

subset_below_root <- function(graph, k = 1) {
  igraph::induced_subgraph(graph, igraph::V(graph)[igraph::V(graph)$dist_root <= k])
}

set_start_vertex_and_subset_all_below <- function(graph, key = "name", value) {
  start <- igraph::V(graph)[which(igraph::vertex_attr(graph)[[key]] == value)]
  desc <- igraph::subcomponent(graph, start, mode = "out")
  return(igraph::induced_subgraph(graph, desc))
}


##rrvgo
#
# simMatrix <- calculateSimMatrix(
#   purrr::set_names(gseadf$gs_exact_source)[1:20],
#   orgdb = "org.Hs.eg.db",
#   ont = "BP",
#   method = "Rel"
# )
# reduced <- reduceSimMatrix(
#   simMatrix,
#   setNames(-log10(gseadf$padj)[1:20], gseadf$gs_exact_source[1:20]),
#   threshold = 0.7,
#   orgdb = "org.Hs.eg.db"
# )

