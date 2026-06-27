#' Plot QC Parameters by Metadata Columns
#'
#' Generates jitter and optional boxplots of quality-control metadata columns
#' grouped by one or more categorical metadata columns from a Seurat object.
#'
#' @param SO A Seurat object.
#' @param meta_cols Character vector of non-numeric metadata columns used for grouping on the x-axis.
#' @param qc_cols Character vector of numeric QC metadata columns to plot.
#' @param cells Optional cells to include. If `NULL`, all cells are used.
#' @param theme A ggplot2 theme object.
#' @param color_by Optional non-numeric metadata column used to color points.
#' @param log Logical. If `TRUE`, applies a log10 scale to the y-axis.
#' @param boxplot Logical. If `TRUE`, overlays boxplots on the jittered points.
#' @param width Numeric jitter width passed to `ggplot2::geom_jitter()`.
#' @param size Numeric point size passed to `ggplot2::geom_jitter()`.
#' @param col_pal Color palette passed to `colrr::make_col_pal()` and used for `scale_color_manual()`.
#'
#' @return A `ggplot` object showing QC parameter distributions faceted by QC
#'   parameter and metadata column.
#' @export
#'
#' @examples
#' \dontrun{
#' qc_params_meta_cols(
#'   SO,
#'   meta_cols = "orig.idents",
#'   qc_cols = c("nFeature_RNA", "nCount_RNA", "pct_mt")
#' )
#' }

qc_params_meta_cols <- function (SO,
                                 meta_cols = "orig.idents",
                                 qc_cols = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                                 cells = NULL,
                                 theme = colrr::theme_material(style = "prismy", white = T),
                                 color_by = NULL,
                                 log = F,
                                 boxplot = T,
                                 width = 0.2,
                                 size = 0.4,
                                 col_pal = "custom") {

  if (missing(meta_cols)) {
    stop("Please provide meta_cols.")
  }
  SO <- scexpr:::check.SO(SO = SO, assay = SeuratObject::DefaultAssay(SO), length = 1)
  meta_cols <- scexpr:::check.features(SO = SO, features = unique(meta_cols), rownames = F)
  qc_cols <- scexpr:::check.features(SO = SO, features = unique(qc_cols), rownames = F)
  color_by <- scexpr:::check.features(SO = SO, features = unique(color_by), rownames = F)
  cells <- scexpr:::check.and.get.cells(SO, assay = SeuratObject::DefaultAssay(SO), cells = cells)

  if (length(intersect(qc_cols, meta_cols)) > 0) {
    stop("qc_cols and meta_cols cannot contain the same columns.")
  }
  if (any(sapply(SO@meta.data[,meta_cols,drop=F], is.numeric))) {
    stop(paste0("numeric meta_cols cannot be handled yet: ", paste(names(which(sapply(SO@meta.data[,c(meta_cols)], is.numeric))), collapse = ", ")))
  }

  data <- tidyr::pivot_longer(
    SO@meta.data[names(cells[which(cells == 1)]),c(qc_cols, meta_cols, color_by)],
    cols = dplyr::all_of(qc_cols),
    names_to = "qc_param",
    values_to = "value"
  )
  data <- tidyr::pivot_longer(
    data,
    cols = dplyr::all_of(meta_cols),
    names_to = "meta.col",
    values_to = "level"
  )


  if (!is.null(color_by)) {
    if (is.numeric(data[[color_by]])) {
      stop("color_by cannot be numeric.")
    }
    data[[color_by]] <- forcats::fct_infreq(as.character(data[[color_by]]))
    data <- data |> dplyr::arrange(!!rlang::sym(color_by))
  }

  col_pal <- colrr::make_col_pal(col_pal)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = level, y = value))
  if (!is.null(color_by)) {
    plot <- plot + ggplot2::geom_jitter(width = width, size = size, ggplot2::aes(color = !!rlang::sym(color_by)))
  } else {
    plot <- plot + ggplot2::geom_jitter(width = width, size = size)
  }

  if (boxplot) {
    plot <- plot + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6)
  }
  plot <- plot +
    theme +
    ggplot2::scale_color_manual(values = col_pal) +
    ggplot2::facet_grid(rows = ggplot2::vars(qc_param), cols = ggplot2::vars(meta.col), scales = "free", axes = "all", axis.labels = "margins") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))

  if (log) {
    plot <- plot + ggplot2::scale_y_log10()
  }

  return(plot)
}
