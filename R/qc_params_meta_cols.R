#' Title
#'
#' @param SO
#' @param qc_cols
#' @param meta_cols
#' @param theme
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
qc_params_meta_cols <- function (SO,
                                 cells = NULL,
                                 meta_cols,
                                 qc_cols = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                                 theme = ggplot2::theme_bw()) {

  if (missing(meta_cols)) {
    stop("Please provide meta_cols.")
  }
  SO <- check.SO(SO = SO, assay = SeuratObject::DefaultAssay(SO), length = 1)
  meta_cols <- check.features(SO = SO, features = unique(meta_cols), rownames = F)
  qc_cols <- check.features(SO = SO, features = unique(qc_cols), rownames = F)
  cells <- check.and.get.cells(SO,  assay = SeuratObject::DefaultAssay(SO), cells = cells)

  if (length(intersect(qc_cols, meta_cols)) > 0) {
    stop("qc_cols and meta_cols cannot contain the same columns.")
  }
  if (any(sapply(SO@meta.data[,meta_cols,drop=F], is.numeric))) {
    stop(paste0("numeric meta_cols cannot be handled yet: ", paste(names(which(sapply(SO@meta.data[,c(meta_cols)], is.numeric))), collapse = ", ")))
  }

  data <- tidyr::pivot_longer(SO@meta.data[names(cells[which(cells == 1)]),c(qc_cols, meta_cols)], cols = dplyr::all_of(qc_cols), names_to = "qc_param", values_to = "value")
  data <- tidyr::pivot_longer(data, cols = dplyr::all_of(meta_cols), names_to = "meta.col", values_to = "level")

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = level, y = value)) +
    ggplot2::geom_jitter(width = 0.2, size = 0.4) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    theme +
    ggplot2::facet_grid(rows = ggplot2::vars(qc_param), cols = ggplot2::vars(meta.col), scales = "free")

 return(plot)
}
