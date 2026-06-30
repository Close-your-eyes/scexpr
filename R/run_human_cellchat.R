#' Run CellChat Analysis on a Human Seurat Object
#'
#' Creates and processes a CellChat object from a Seurat object using the human
#' CellChat ligand-receptor database. The function groups cells by the specified
#' metadata column, identifies overexpressed genes and interactions, computes
#' communication probabilities, filters low-support interactions, computes pathway
#' probabilities, and aggregates the inferred communication network.
#'
#' @param obj A Seurat object containing RNA expression data.
#' @param metal_col A character string specifying the metadata column in `obj`
#'   used to group cells for CellChat analysis.
#'
#' @returns A processed `CellChat` object with inferred cell-cell communication
#'   probabilities and aggregated networks.
#' @export
#'
#' @examples
#' \dontrun{
#' cellchat <- run_human_cellchat(seurat_obj, metal_col = "cell_type")
#' }
run_human_cellchat <- function(obj, metal_col) {
  Seurat::DefaultAssay(obj) <- "RNA"
  cellchat <- CellChat::createCellChat(object = obj, group.by = metal_col)
  cellchat@DB <- CellChat::CellChatDB.human
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  return(cellchat)
}
