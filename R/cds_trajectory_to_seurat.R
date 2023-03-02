#' Copy a trajectory from a cell data set (cds) object to respective the Misc slot of a Seurat object (SO)
#'
#' The code lines are copied from monocle3::plot_cells
#' The aim is to have all data in a Seurat object to allow plot creation with a common style/layout.
#'
#' @param cds a cell data set
#' @param SO an Seurat object
#' @param x
#' @param y
#' @param reduction_method
#' @param name
#'
#' @return
#' @export
#'
#' @examples
cds_trajectory_to_seurat <- function(cds, SO, x = 1, y = 2, reduction_method = "UMAP", name = "trajectory") {


  reduction_method <- match.arg(reduction_method, "UMAP")

  if (!max(c(x,y)) >= ncol(t(cds@principal_graph_aux[[reduction_method]]$dp_mst))) {
    stop("cds@principal_graph_aux[[reduction_method]]$dp_mst does not have sufficient columns for x/y selection.")
  }

  ### these lines are taken from moncole3::plot_cells

  ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
    dplyr::mutate(sample_name = rownames(.),
                  sample_state = rownames(.))

  Seurat::Misc(SO, slot = name) <-
    cds@principal_graph[[reduction_method]] %>%
    igraph::as_data_frame() %>%
    dplyr::select(from, to) %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select(
                         from="sample_name",
                         from_x="prin_graph_dim_1",
                         from_y="prin_graph_dim_2"),
                     by = "from") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select(
                         to="sample_name",
                         to_x="prin_graph_dim_1",
                         to_y="prin_graph_dim_2"),
                     by = "to")

  # add trajectory to ggplot with:
  #g <- g + geom_segment(aes(x = from_x, y = from_y, xend = to_x, yend = to_y), size=trajectory_graph_segment_size, color=I(trajectory_graph_color), linetype="solid", na.rm=TRUE, data=as.data.frame(Seurat::Misc(SO, slot = name)))

  return(SO)
}
