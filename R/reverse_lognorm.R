#' Recreate counts layer from data layer
#'
#' Remove counts layer to save memory. Recreate it when needed.
#'
#' @param obj Seurat object
#' @param assay assay to use
#' @param nCount_RNA vector of colsums of RNA counts (lib sizes)
#' @param scale.factor factor as in Seurat::LogNormalize
#'
#' @return counts matrix to write to assay slot
#' @export
#'
#' @examples
#' \dontrun{
#'   ## lognorm forward:
#' csums <- Matrix::colSums(counts)
#' xx <- sweep(counts, 2, csums, "/")
#' xx <- xx*1e4
#' data <- log1p(xx)
#'
#' ## SCT reverse:
#' # simply expm1 from data slot, see ?Seurat::SCTransform
#'
#' # compare count matrices
#' rmse <- scexpr:::rmse_sparse(count_re, counts_orig)
#' equal <- all.equal(count_re, counts_orig)
#' }
reverse_lognorm <- function(obj,
                            assay = "RNA",
                            nCount_RNA = Seurat::Misc(obj, slot = "RNA_count_colSums"),
                            scale.factor = 10000) {

  assay <- match.arg(assay, names(obj@assays))

  if (assay == "RNA") {

    if (is.null(nCount_RNA)) {
      stop("nCount_RNA not found.")
      # warning("Misc slot in obj not found. Assuming lowest value > 0 per column to represent UMI = 1 originally.")
      # nCount_RNA <- apply(SeuratObject::LayerData(obj, layer = "data", assay = assay), 2, function(x) min(x[which(x > 0)]))
      # nCount_RNA <- unname(1 / (expm1(nCount_RNA) / 10000))
    }

    if (!is.vector(nCount_RNA)) {
      stop("nCount_RNA must be a vector.")
    }
    if (length(nCount_RNA) != ncol(get_layer(obj = obj, layer = "data", assay = "RNA"))) {
      stop("length of nCount_RNA does not match data layer columns.")
    }

    # keep sparsity
    counts <- expm1(get_layer(obj = obj, layer = "data", assay = "RNA")) %*% Matrix::Diagonal(x = nCount_RNA / scale.factor)
    # make them integers
    counts@x <- round(counts@x)
    counts <- Matrix::drop0(counts)
    return(counts)
    #return(sweep(expm1(get_layer(obj = obj, layer = "data", assay = "RNA"))/scale.factor, 2, nCount_RNA, "*"))
  }



  if (assay == "SCT") {
    counts <- expm1(get_layer(obj = obj, layer = "data", assay = "SCT"))
    counts@x <- round(counts@x)
    counts <- Matrix::drop0(counts)
    return(counts)
  }


}

#' Remove the RNA counts layer while preserving cell totals
#'
#' Removes the `"counts"` layer from the RNA assay to reduce the size of a
#' Seurat object. Before removal, the total count for each cell is stored in
#' the object's miscellaneous data under `"RNA_count_colSums"`.
#'
#' The stored column sums can be used by [reverse_lognorm()] to approximate
#' the original counts from the log-normalized data layer.
#'
#' @param obj A Seurat object containing an RNA assay with a `"counts"` layer.
#'
#' @return A Seurat object with the RNA counts layer removed and its per-cell
#'   column sums stored in `Seurat::Misc()`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- sacrifice_count_slot(obj)
#'
#' # Reconstruct the counts matrix
#' counts <- reverse_lognorm(obj)
#' }
sacrifice_count_slot <- function(obj, rm_count = F) {

  Seurat::Misc(obj, slot = "RNA_count_colSums") <-
    Matrix::colSums(
      get_layer(obj = obj, layer = "counts", assay = "RNA")
    )

  if (rm_count) {
    obj@assays$RNA@layers$counts <- NULL
  }
  return(obj)

  # Use scexpr::reverse_lognorm() to reconstruct the counts matrix
}



rmse_sparse <- function(A, B) {
  diff <- A - B
  sqrt(sum(diff@x^2) / (nrow(diff) * ncol(diff)))
}
