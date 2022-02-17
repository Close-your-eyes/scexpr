#' Title
#'
#' @param SO
#' @param colsSums
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' #### for RNA assay:
#' # Default LogNormalize from Seurat does this:
#' nCount_RNA <- Matrix::colSums(GetAssayData(SO, slot = "counts", assay = "RNA"))
#' # devide each column by corresponding index in nCount_RNA
#' norm_UMI <- sweep(Seurat::GetAssayData(SO, slot = "counts", assay = "RNA"), 2, nCount_RNA, FUN = '/')
#' norm_UMI <- norm_UMI*10000
#' lognorm_UMI <- log1p(norm_UMI)
#'
#' # to reverse (data --> Counts) one needs original nCount_RNA (colSums)
#' # this must be saved (e.g. Misc slot) before deleting the Count slot
#' # if not provided one could assume that the lowest norm_UMI
#' # per column equals to 1 UMI originally
#' norm_UMI <- expm1(GetAssayData(SO, slot = "data", assay = "RNA")) #lognorm_UMI
#' norm_UMI <- norm_UMI/10000
#' UMI_count <- sweep(norm_UMI, 2, nCount_RNA, FUN = '*')
#'
#' #### for SCT assay:
#' # simply expm1 from data slot, see ?Seurat::SCTransform
#' }
reverse_lognorm <- function(SO,
                            assay = c("RNA", "SCT"),
                            nCount_RNA = Seurat::Misc(SO, slot = "RNA_count_colSums"),
                            scale.factor = 10000,
                            return.SO = T) {

  assay <- match.arg(assay, c("RNA", "SCT"))

  if (assay == "RNA") {
    if (is.null(nCount_RNA)) {
      warning("Misc slot in SO not found. Assuming lowest value > 0 per column to represent UMI = 1 originally.")
      nCount_RNA <- apply(Seurat::GetAssayData(SO, slot = "data", assay = "RNA"), 2, function(x) min(x[which(x>0)]))
      nCount_RNA <- unname(1/(expm1(nCount_RNA)/10000))
    } else {
      if (!is.vector(nCount_RNA)) {
        stop("nCount_RNA must be a vector.")
      }
      if (length(nCount_RNA) != ncol(Seurat::GetAssayData(SO, slot = "data", assay = "RNA"))) {
        stop("length(nCount_RNA) != ncol(GetAssayData(SO, slot = 'data', assay = 'RNA'))")
      }
    }

    if (return.SO) {
      if (any(dim(SO@assays[["RNA"]]@counts) != c(0,0))) {
        warning("Counts slot of RNA assay has not been empty before and will be overwritten.")
      }
      SO@assays[["RNA"]]@counts <- sweep(expm1(Seurat::GetAssayData(SO, slot = "data", assay = "RNA"))/scale.factor, 2, nCount_RNA, FUN = '*')
      return(SO)
    } else {
      return(sweep(expm1(Seurat::GetAssayData(SO, slot = "data", assay = "RNA"))/scale.factor, 2, nCount_RNA, FUN = '*'))
    }
  }

  if (assay == "SCT") {
    if (return.SO) {
      if (any(dim(SO@assays[["SCT"]]@counts) != c(0,0))) {
        warning("Counts slot of SCT assay has not been empty before and will be overwritten.")
      }
      SO@assays[["SCT"]]@counts <- expm1(Seurat::GetAssayData(SO, slot = "data", assay = "SCT"))
      return(SO)
    } else {
      return(expm1(Seurat::GetAssayData(SO, slot = "data", assay = "SCT")))
    }
  }


  # did not work yet, strange problem (R hell?!)
  #identical(out, GetAssayData(SO, slot = "counts", assay = "RNA"))

  # LogNorm
  #res <- log1p(sweep(GetAssayData(SO, slot = "counts", assay = "RNA"), 2, Matrix::colSums(GetAssayData(SO, slot = "counts", assay = "RNA")), FUN = '/')*10000)
  #identical(res, GetAssayData(SO, slot = "data", assay = "RNA"))
}


