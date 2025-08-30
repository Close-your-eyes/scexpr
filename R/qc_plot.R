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




#' Title
#'
#' @param SO
#' @param qc_cols
#' @param geom2
#' @param save_path
#'
#' @returns
#' @export
#'
#' @examples
qc_plot2 <- function(SO,
                     qc_cols = c(
                       "nCount_RNA_log",
                       "nFeature_RNA_log",
                       "pct_mt_log",
                       "dbl_score_log",
                       "pct_soup_decontX",
                       "pct_soup_SoupX"
                     ),
                     geom2 = "boxplot",
                     save_path = NULL,
                     save_as = c("png", "pdf")) {


  # SO@misc data is used here!
  save_as <- rlang::arg_match(save_as)

  if (any(!qc_cols %in% names(SO@meta.data))) {
    message(paste(qc_cols[which(!qc_cols %in% names(SO@meta.data))], collapse = ", "), " not found in SO@meta.data. These will be excluded from plotting.")
    qc_cols <- intersect(qc_cols, names(SO@meta.data))
    if (length(qc_cols) == 0) {
      stop("No qc_cols left for plotting.")
    }
  }

  breaks <- c(seq(0, 1e1, 2e0),
              seq(0, 1e2, 2e1),
              seq(0, 1e3, 2e2),
              seq(0, 1e4, 2e3),
              seq(0, 1e5, 2e4))
  breaks <- breaks[which(breaks != 0)]


  qc_p1 <- suppressMessages(feature_plot2(SO,
                                          features = c(gsub("_log$", "", qc_cols), "orig.ident", SO@misc$clusterings),
                                          reduction = "UMAP"))

  # qc_p2 <- ggplot2::ggplot(tidyr::pivot_longer(SO@meta.data[,c(qc_cols, clustering_cols)], cols = dplyr::all_of(qc_cols), names_to = "qc_param", values_to = "value"),
  #                          ggplot2::aes(x = !!rlang::sym(clustering_cols[1]), y = value, color = !!rlang::sym(clust_name))) +
  #   ggplot2::geom_boxplot(color = "grey30", outlier.shape = NA) +
  #   ggplot2::geom_jitter(width = 0.1, size = 0.3) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
  #                  panel.grid.minor.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
  #                  legend.key.size = ggplot2::unit(0.3, "cm"), legend.key = ggplot2::element_blank()) +
  #   ggplot2::scale_color_manual(values = colrr::col_pal("custom")) +
  #   ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ expm1(.), breaks = breaks)) +
  #   ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3))) +
  #   ggplot2::facet_wrap(ggplot2::vars(qc_param), scales = "free_y", ncol = 1)

  mc <- names(SO@misc$meta_clusterings)
  qc_p2 <- purrr::map(stats::setNames(mc, mc), function(red_name) {
    clust_name <- SO@misc$meta_clusterings[red_name]
    p1 <- patchwork::wrap_plots(feature_plot2(SO,
                                              features = clust_name,
                                              reduction = red_name,
                                              pt.size = 0.5,
                                              col_pal_d_args = list(name = SO@misc$clustering_colors)),
                                feature_plot2(SO,
                                              features = "orig.ident",
                                              reduction = red_name,
                                              pt.size = 0.5,
                                              col_pal_d_args = list(name = SO@misc$orig.ident_colors)),
                                freq_pie_chart(
                                  SO = SO,
                                  meta_col = clust_name,
                                  fill = SO@misc$clustering_colors,
                                  label_rel_cutoff = 0.05)[["plot"]],
                                ncol = 1,
                                heights = c(0.35,0.35,0.3))
    p2 <- purrr::map(stats::setNames(qc_cols, qc_cols), function(feature) {
      make_stat_plot(
        SO = SO,
        feature = feature,
        clust_name = clust_name,
        geom2 = geom2,
        breaks = breaks
      )
    })

    return(patchwork::wrap_plots(p1, patchwork::wrap_plots(p2, ncol = 2), nrow = 1, widths = c(0.4, 0.6)))
    #return(list(pheno_clustering = p3_1, meta_boxplots = p3_2))
  })



  # p3_x <- lapply(qc_cols[which(!grepl("nCount_RNA_log|nFeature_RNA_log|pct_mt_log", qc_cols))], function(qcf) {
  #   p3_x <- feature_plot_stat(SO,
  #                             features = qcf,
  #                             meta_col = clust_name,
  #                             geom2 = geom2,
  #                             jitterwidth = 0.9) +
  #     ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey95"),
  #                    axis.title.y = ggplot2::element_blank())
  #   if (qcf != rev(qc_cols[which(!grepl("nCount_RNA_log|nFeature_RNA_log|pct_mt_log", qc_cols))])[1]) {
  #     p3_x <-
  #       p3_x +
  #       ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())
  #   }
  #   return(p3_x)
  # })

  qc_p3 <- qc_params_meta_cols(SO, meta_cols = "orig.ident", qc_cols = gsub("_log", "", qc_cols))


  message("save plots with: (replace out and, im_path)")
  show_command(ggplot2::ggsave(filename = "pheno_qc.png", plot = out[["pheno"]], device = png, path = im_path, width = 10, height = 7))
  show_command(for (i in names(out[["meta"]])) {
    ggplot2::ggsave(filename = paste0(i, "_qc.png"), plot = out[["meta"]][[i]], device = png, path = im_path, width = 12, height = 9)
  })
  show_command(ggplot2::ggsave(filename = "meta2_qc.png", plot = out[["meta2"]], device = png, path = im_path, width = 5, height = 8))

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = paste0("pheno_qc.", save_as),
      plot = qc_p1,
      device = save_as,
      path = save_path,
      width = 10,
      height = 7
    )
    for (i in names(qc_p2)) {
      ggplot2::ggsave(
        filename = paste0(i, "_qc.", save_as),
        plot = qc_p2[[i]],
        device = save_as,
        path = save_path,
        width = 12,
        height = 9
      )
    }
    ggplot2::ggsave(
      filename = paste0("meta2_qc.", save_as),
      plot = qc_p3,
      device = save_as,
      path = save_path,
      width = 4+0.5*length(unique(SO@meta.data$orig.ident)),
      height = 4+length(qc_cols)
    )
  }

  return(list(pheno = qc_p1, meta = qc_p2, meta2 = qc_p3))

}


make_stat_plot <- function(SO, feature, clust_name, geom2, breaks) {

  plot <- feature_plot_stat(SO,
                            features = feature,
                            col.pal = SO@misc$clustering_colors,
                            meta_col = clust_name,
                            geom2 = geom2,
                            jitterwidth = 0.3) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey95"),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())

  if (grepl("_log$", feature)) {
    plot <- plot +
      ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ expm1(.), breaks = breaks[intersect(which(breaks > min(expm1(SO@meta.data[[feature]]))),
                                                                                                     which(breaks < max(expm1(SO@meta.data[[feature]]))))]))
  }
  return(plot)
}

show_command <- function(expr) {
  print(substitute(expr))
}
