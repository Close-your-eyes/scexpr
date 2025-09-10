#' Recreate counts layer from data layer
#'
#' Remove counts layer to save memory. Recreate it when needed.
#'
#' @param SO Seurat object
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
reverse_lognorm <- function(SO,
                            assay = "RNA",
                            nCount_RNA = Seurat::Misc(SO, slot = "RNA_count_colSums"),
                            scale.factor = 10000) {

  assay <- match.arg(assay, names(SO@assays))



  if (assay == "RNA") {

    if (is.null(nCount_RNA)) {
      stop("nCount_RNA not found.")
      # warning("Misc slot in SO not found. Assuming lowest value > 0 per column to represent UMI = 1 originally.")
      # nCount_RNA <- apply(SeuratObject::LayerData(SO, layer = "data", assay = assay), 2, function(x) min(x[which(x > 0)]))
      # nCount_RNA <- unname(1 / (expm1(nCount_RNA) / 10000))
    }

    if (!is.vector(nCount_RNA)) {
      stop("nCount_RNA must be a vector.")
    }
    if (length(nCount_RNA) != ncol(get_layer(obj = SO, layer = "data", assay = "RNA"))) {
      stop("length of nCount_RNA does not match data layer columns.")
    }

    return(sweep(expm1(get_layer(obj = SO, layer = "data", assay = "RNA"))/scale.factor, 2, nCount_RNA, "*"))
  }


  if (assay == "SCT") {
    return(expm1(get_layer(obj = SO, layer = "data", assay = "SCT")))
  }


}



rmse_sparse <- function(A, B) {
  diff <- A - B
  sqrt(sum(diff@x^2) / (nrow(diff) * ncol(diff)))
}
