#' Quickly plot common qc metrics
#'
#' Plot qc metric, optional with a linear secondary axis.
#' When SO1 only is provided: Only plot cells from one Seurat object.
#' When SO1 and SO2 are provided:
#' Plot additionally those cells that have been filtered from that object with more cells.
#' Orig.
#'
#' @param SO1 Seurat object 1
#' @param SO2 Seurat object 2, optional
#' @param qc_cols QC columns on y-axes
#' @param sec_axis_lin plot a secondary axis with linear transformation of qc cols
#' @param x_cat meta data column for x-axis
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
qc_plot <- function(SO1,
                    SO2 = NULL,
                    qc_cols = c("nCount_RNA_log", "nFeature_RNA_log", "pct_mt_log", "dbl_score_log"),
                    x_cat = "orig.ident",
                    sec_axis_lin = T) {

  breaks <- c(seq(0, 1e1, 2e0),
              seq(0, 1e2, 2e1),
              seq(0, 1e3, 2e2),
              seq(0, 1e4, 2e3),
              seq(0, 1e5, 2e4))
  breaks <- breaks[which(breaks != 0)]

  qc_cols <- intersect(qc_cols, names(SO1@meta.data))
  if (!is.null(SO2)) {
    qc_cols <- intersect(qc_cols, names(SO2@meta.data))
  }
  if (length(qc_cols) == 0) {
    stop("Non of qc_cols found in SOs.")
  }

  if(!x_cat %in% names(SO1@meta.data)) {
    stop("x_cat not found in SO1.")
  }

  if(!is.null(SO2) && !x_cat %in% names(SO1@meta.data)) {
    stop("x_cat not found in SO2.")
  }

  if (!is.logical(sec_axis_lin)) {
    stop("sec_axis_lin has to be logical: TRUE or FALSE.")
  }

  if (!is.null(SO2)) {
    SO_list <- list(SO1, SO2)
  } else {
    SO_list <- list(SO1)
  }
  SO_list <- SO_list[order(sapply(SO_list, ncol))]

  meta_data <-
    SO_list[[1]]@meta.data %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::select(ID, !!rlang::sym(x_cat), dplyr::all_of(qc_cols)) %>%
    tidyr::pivot_longer(cols = -c(ID, dplyr::all_of(x_cat)))

  ppp <- ggplot2::ggplot(meta_data, ggplot2::aes(x = !!rlang::sym(x_cat), y = value)) +
    ggplot2::geom_jitter(width = 0.1, color = "grey60", size = 0.3) +
    ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::facet_wrap(ggplot2::vars(name), scales = "free_y")


  if (length(SO_list) > 1) {
    meta_data2 <-
      SO_list[[2]]@meta.data %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::select(ID, !!rlang::sym(x_cat), dplyr::all_of(qc_cols)) %>%
      tidyr::pivot_longer(cols = -c(ID, !!rlang::sym(x_cat))) %>%
      dplyr::filter(!ID %in% meta_data$ID)
    diff <- ncol(SO_list[[2]]) - ncol(SO_list[[1]])
    message(diff, " cells were removed.")
    ppp <- ppp + ggplot2::geom_jitter(data = meta_data2, width = 0.1, color = "tomato2", size = 0.3)
  }

  if (sec_axis_lin) {
    ppp <- ppp + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~expm1(.), breaks = breaks))
  }
  return(ppp)
}
