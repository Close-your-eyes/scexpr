#' Calculate monocle trajectory and add to seurat object
#'
#' @param obj seurat object with umap dr
#'
#' @returns seurat object
#' @export
#'
#' @examples
#'\dontrun{
#' scexpr::feature_plot2(obj, "umap_partition", "umap") +
#'   ggnewscale::new_scale_color() +
#'   ggplot2::geom_segment(ggplot2::aes(x=source_prin_graph_dim_1,
#'                                      y=source_prin_graph_dim_2,
#'                                      xend=target_prin_graph_dim_1,
#'                                      yend=target_prin_graph_dim_2,
#'                                      color = source_umap_partition),
#'                         size=0.75,
#'                         linetype="solid",
#'                         na.rm=TRUE,
#'                         data=obj@misc[["umap_edge_df"]])
#'}
run_monocle <- function(obj) {

  # BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
  #                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
  #                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
  #                        'ggrastr'))
  # remotes::install_github("bnprks/BPCells/r")
  # remotes::install_github("https://github.com/cran/grr")
  # #sudo port install hdf5
  # pak::pak('cole-trapnell-lab/monocle3')
  # remotes::install_github('satijalab/seurat-wrappers')

  #obj <- readRDS("/Volumes/CMS_SSD_2TB/R_scRNAseq/2020_10XGenomics_PBMCs/data/SO_processed/full_objects/SO_SC3_v3_NextGem_SI_PBMC_10K_SCT_none_1_500_10_220617-154703.rds")

  if (!requireNamespace("seurat-wrappers", quietly = T)) {
    remotes::install_github("satijalab/seurat-wrappers")
  }
  reduction_method <- "UMAP"
  cds <- SeuratWrappers::as.cell_data_set(obj)
  # use pre-calculated umap
  cds <- monocle3::cluster_cells(
    cds,
    reduction_method = reduction_method,
    cluster_method = "louvain"
  )
  cds <- monocle3::learn_graph(cds)
  cds <- monocle3::order_cells(cds)

  obj <- Seurat::AddMetaData(obj, monocle3::pseudotime(cds, reduction_method = reduction_method), col.name = "umap_pseudotime")
  obj <- Seurat::AddMetaData(obj, monocle3::partitions(cds, reduction_method = reduction_method), col.name = "umap_partition") # cds@clusters@listData[["UMAP"]][["partitions"]]
  obj <- Seurat::AddMetaData(obj, monocle3::clusters(cds, reduction_method = reduction_method), col.name = "umap_cluster") # cds@clusters@listData[["UMAP"]][["clusters"]]
  obj <- Seurat::AddMetaData(obj, cds@principal_graph_aux[[reduction_method]][["pr_graph_cell_proj_closest_vertex"]][,1], col.name = "umap_node")

  # mode per node
  node_avg <- obj@meta.data[,c("umap_pseudotime", "umap_partition", "umap_cluster", "umap_node")] |>
    dplyr::summarise(umap_pseudotime = mean(umap_pseudotime),
                     umap_partition = names(which.max(table(umap_partition))),
                     umap_cluster = names(which.max(table(umap_cluster))),
                     .by = "umap_node")

  x <- 1
  y <- 2
  ica_space_df <- as.data.frame(t(cds@principal_graph_aux[[reduction_method]][["dp_mst"]])) |>
    dplyr::select(prin_graph_dim_1 = {{x}}, prin_graph_dim_2 = {{y}}) |>
    tibble::rownames_to_column("sample_name")
  edge_df <- igraph::as_data_frame(monocle3::principal_graph(cds)[[reduction_method]]) |>
    dplyr::select(source_node = "from", target_node = "to") |>
    dplyr::left_join(dplyr::select(ica_space_df,
                                   source_node="sample_name",
                                   source_prin_graph_dim_1="prin_graph_dim_1",
                                   source_prin_graph_dim_2="prin_graph_dim_2"),
                     by = "source_node") |>
    dplyr::left_join(dplyr::select(ica_space_df,
                                   target_node="sample_name",
                                   target_prin_graph_dim_1="prin_graph_dim_1",
                                   target_prin_graph_dim_2="prin_graph_dim_2"),
                     by = "target_node") |>
    dplyr::mutate(target_node = as.numeric(gsub("Y_", "", target_node)), source_node = as.numeric(gsub("Y_", "", source_node))) |>
    dplyr::left_join(dplyr::select(node_avg,
                                   source_node = "umap_node",
                                   source_umap_pseudotime = "umap_pseudotime",
                                   source_umap_partition = "umap_partition",
                                   source_umap_cluster = "umap_cluster"),
                     by = "source_node") |>
    dplyr::left_join(dplyr::select(node_avg,
                                   target_node = "umap_node",
                                   target_umap_pseudotime = "umap_pseudotime",
                                   target_umap_partition = "umap_partition",
                                   target_umap_cluster = "umap_cluster"),
                     by = "target_node")
  Seurat::Misc(obj, "umap_edge_df") <- edge_df

  return(obj)

  # scexpr::feature_plot2(obj, "umap_partition", "umap") +
  #   ggnewscale::new_scale_color() +
  #   ggplot2::geom_segment(ggplot2::aes(x=source_prin_graph_dim_1,
  #                                      y=source_prin_graph_dim_2,
  #                                      xend=target_prin_graph_dim_1,
  #                                      yend=target_prin_graph_dim_2,
  #                                      color = source_umap_partition),
  #                         size=0.75,
  #                         linetype="solid",
  #                         na.rm=TRUE,
  #                         data=edge_df)
}


# expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
# cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
# gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))
# cds <- new_cell_data_set(expression_matrix,
#                          cell_metadata = cell_metadata,
#                          gene_metadata = gene_annotation)
