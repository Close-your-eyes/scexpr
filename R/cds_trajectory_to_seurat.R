#' Copy a trajectory from a cell data set (cds) object to respective the Misc slot of a Seurat object (SO)
#'
#' The code lines are copied from monocle3::plot_cells
#' The aim is to have all data in a Seurat object to allow plot creation with a common style/layout.
#'
#' @param cds a cell data set
#' @param SO an Seurat object
#' @param x column number of reduction_method to use as x-axis
#' @param y column number of reduction_method to use as y-axis
#' @param reduction_method name of reduction, currently only UMAP allows, according to monocle3
#' @param name target slot name in Seurat::Misc; will be overwritten if already existing
#' @param pseudotime_metacol_name name of new column in SO@meta.data for pseudotime information; set to NULL to not add a column to meta.data of SO
#'
#' @return
#' @export
#'
#' @examples
cds_trajectory_to_seurat <- function(cds,
                                     SO,
                                     x = 1,
                                     y = 2,
                                     reduction_method = "UMAP",
                                     name = "trajectory",
                                     pseudotime_metacol_name = "pseudotime") {

  if (!methods::is(cds, "cell_data_set")) {
    stop("cds has to be an cell_data_set object.")
  }
  if (!methods::is(SO, "Seurat")) {
    stop("SO has to be a Seurat object.")
  }

  if (!requireNamespace("igraph", quietly = T)) {
    utils::install.packages("igraph")
  }

  reduction_method <- match.arg(reduction_method, "UMAP")

  if (!max(c(x,y)) >= ncol(t(cds@principal_graph_aux[[reduction_method]]$dp_mst))) {
    stop("cds@principal_graph_aux[[reduction_method]]$dp_mst does not have sufficient columns for x/y selection.")
  }

  ### these lines have been taken (and adjusted) from moncole3::plot_cells
  ## prepare data frame for plotting with geom_segment
  ## plotting with ggraph should also be possible though; or not? maybe ggraph on top of ggplot object is not possible

  ica_space_df <-
    t(monocle3::principal_graph_aux(cds)[[reduction_method]][["dp_mst"]]) %>%
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(c(x,y))) %>%
    dplyr::mutate(sample_name = rownames(.))

  df <-
    monocle3::principal_graph(cds)[[reduction_method]] %>%
    igraph::as_data_frame() %>%
    dplyr::select(from, to) %>%
    dplyr::left_join(ica_space_df %>% dplyr::select(sample_name, V1, V2), by = c("from" = "sample_name")) %>%
    dplyr::rename("from_x" = V1, "from_y" = V2) %>%
    dplyr::left_join(ica_space_df %>% dplyr::select(sample_name, V1, V2), by = c("to" = "sample_name")) %>%
    dplyr::rename("to_x" = V1, "to_y" = V2)


  if (!is.null(pseudotime_metacol_name)) {
    if ("pseudotime" %in% names(monocle3::principal_graph_aux(cds)[[reduction_method]])) {
      pt <- utils::stack(monocle3::principal_graph_aux(cds)[[reduction_method]][["pseudotime"]])
      rownames(pt) <- pt[,2]
      pt <- pt[,-2,drop=F]
      names(pt) <- pseudotime_metacol_name
      SO <- Seurat::AddMetaData(SO, pt)
    }
  }

  mapping_cells_to_vertices <-
    as.data.frame(t(monocle3::principal_graph_aux(cds)[[reduction_method]][["dp_mst"]])) %>%
    dplyr::mutate(vertex = rownames(.)) %>%
    dplyr::left_join(as.data.frame(monocle3::principal_graph_aux(cds)[[reduction_method]][["pr_graph_cell_proj_closest_vertex"]]) %>%
                       dplyr::rename("vertex" = 1) %>%
                       dplyr::mutate(vertex = paste0("Y_", vertex), ID = rownames(.)),
                     by = "vertex",
                     multiple = "all") %>%
    dplyr::left_join(as.data.frame(SingleCellExperiment::reducedDims(cds)[["UMAP"]]) %>% dplyr::mutate(ID = rownames(.)) %>% dplyr::rename("from_x" = V1, "from_y" = V2), by = "ID") %>%
    dplyr::rename("to_x" = V1, "to_y" = V2)

  vertices <-
    mapping_cells_to_vertices %>%
    dplyr::distinct(vertex, to_x, to_y) %>%
    dplyr::rename("UMAP_1" = to_x, "UMAP_2" = to_y)

  Seurat::Misc(SO, slot = name) <- list(df = df,
                                        vertices_df = vertices,
                                        mapping_cells_to_vertices_df = mapping_cells_to_vertices,
                                        principle_graph = monocle3::principal_graph(cds)[[reduction_method]],
                                        principle_graph_aux = monocle3::principal_graph_aux(cds)[[reduction_method]])

  # add trajectory to ggplot with:
  #g <- g + geom_segment(aes(x = from_x, y = from_y, xend = to_x, yend = to_y), size=trajectory_graph_segment_size, color=I(trajectory_graph_color), linetype="solid", na.rm=TRUE, data=as.data.frame(Seurat::Misc(SO, slot = name)))

  return(SO)
}
