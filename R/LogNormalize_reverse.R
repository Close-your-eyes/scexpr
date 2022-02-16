LogNormalize_reverse <- function(SO, colsSums = "RNA_count_colSums") {
  # LogNorm
  #res <- log1p(sweep(GetAssayData(SO_urine, slot = "counts", assay = "RNA"), 2, colSums(GetAssayData(SO_urine, slot = "counts", assay = "RNA")), FUN = '/')*10000)
  #identical(res, GetAssayData(SO_urine, slot = "data", assay = "RNA"))

  # reverse LogNorm
  # where is it wrong?
  out <- sweep(expm1(GetAssayData(SO, slot = "data", assay = "RNA"))/10000, 2, colSums(GetAssayData(SO, slot = "counts", assay = "RNA")), FUN = '*')
  identical(out, GetAssayData(SO, slot = "counts", assay = "RNA"))

  ## make neat

  ## assume lowest value per column to represent an original count = 1
}


