#' Title
#'
#' @param SO
#' @param qc.cols
#' @param meta.cols
#' @param theme
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
qc_params_meta_cols <- function (SO,
                                 qc.cols = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                                 meta.cols,
                                 theme = ggplot2::theme_bw(),
                                 ...) {

  dots <- list(...)

  if (missing(meta.cols)) {
    stop("Please provide meta.cols.")
  }
  SO <- .check.SO(SO = SO, length = 1)
  cells <- do.call(.check.and.get.cells, args = c(list(SO = SO, return.included.cells.only = T), dots[which(names(dots) %in% names(formals(.check.and.get.cells)))]))
  meta.cols <- .check.features(SO = SO, features = unique(meta.cols), rownames = F)
  qc.cols <- .check.features(SO = SO, features = unique(qc.cols), rownames = F)

  if (length(intersect(qc.cols, meta.cols)) > 0) {
    stop("qc.cols and meta.cols cannot contain the same columns.")
  }
  if (any(sapply(SO@meta.data[,meta.cols,drop=F], is.numeric))) {
    stop(paste0("numeric meta.cols cannot be handled yet: ", paste(names(which(sapply(SO@meta.data[,c(meta.cols)], is.numeric))), collapse = ", ")))
  }

  data <- tidyr::pivot_longer(SO@meta.data[cells,c(qc.cols, meta.cols)], cols = dplyr::all_of(qc.cols), names_to = "qc_param", values_to = "value")
  data <- tidyr::pivot_longer(data, cols = dplyr::all_of(meta.cols), names_to = "meta.col", values_to = "level")

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = level, y = value)) +
    ggplot2::geom_jitter(width = 0.2, size = 0.4) +
    theme +
    ggplot2::theme(...) +
    ggplot2::facet_grid(rows = ggplot2::vars(qc_param), cols = ggplot2::vars(meta.col), scales = "free")

 return(plot)
}
