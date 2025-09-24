#' Title
#'
#' @param SO
#' @param features
#' @param assay
#' @param geom1
#' @param geom2
#' @param jitterwidth
#' @param label.size
#' @param meta_col
#' @param plot.expr.freq to plot labels of expression frequencies
#' @param plot.non.expr do plot dots of non-expressers (0 UMIs)
#' @param cells
#' @param downsample
#' @param split.by
#' @param pt.size
#' @param feature.aliases
#' @param feature_cut
#' @param feature_cut_expr
#' @param feature_ex
#' @param font.family
#' @param theme
#' @param ...
#' @param expr.freq.hjust
#' @param legend.title
#' @param col.pal
#' @param col.pal.dir
#' @param expr.freq.decimals decimal precision of frequency labels
#' @param expr.freq.pct plot expression frequency in percent?
#'
#' @return
#' @export
#'
#' @examples
feature_plot_stat <- function(SO,
                              features,
                              meta_col,
                              assay = "RNA",
                              geom1 = c("jitter", "point"),
                              geom2 = c("boxplot", "violin", "none"),
                              plot_first = c("geom1", "geom2"),
                              jitterwidth = 0.2,
                              dodgewidth = 0.9,
                              label.size = 4,
                              plot.expr.freq = F,
                              expr.freq.decimals = 2,
                              expr.freq.pct = F,
                              expr.freq.hjust = 0.3,
                              plot.non.expr = T,
                              plot.strip = T,
                              cells = NULL,
                              downsample = 1,
                              legend.title = "SO.split",
                              split.by = NULL,
                              pt.size = 0.5,
                              col.pal = colrr::col_pal("custom"),
                              col.na = "grey50",
                              col_pal_args = list(missing_fct_to_na = T),
                              feature.aliases = NULL,
                              feature_cut = NULL,
                              feature_cut_expr = 0,
                              feature_ex = NULL,
                              theme = ggplot2::theme_bw(),
                              geom2_args = list(outlier.shape = NA),
                              facetting_args = list(scales = "free_y",
                                                    axes = "all",
                                                    axis.labels = "margins")) {

  if (!requireNamespace("colrr", quietly = T)) {
    devtools::install_github("Close-your-eyes/colrr")
  }

  # add option to plot facet labels in italics

  if (missing(SO)) {
    stop("Seurat object list or feature vector is missing.")
  }
  if (missing(meta_col)) {
    stop("meta_col reauired.")
  }
  if (!is.null(split.by)) {
    message("split.by requires testing.")
  }

  geom1 <- rlang::arg_match(geom1)
  geom2 <- rlang::arg_match(geom2)
  plot_first <- rlang::arg_match(plot_first)

  SO <- check.SO(SO = SO, assay = assay)
  assay <- Seurat::DefaultAssay(SO[[1]])
  features <- check.features(SO = SO, features = features, meta.data = T)
  if (length(meta_col) > 1) {
    stop("Please provide only one meta_col.")
  }

  meta_col <- check.features(SO = SO, features = meta_col, rownames = F)
  cells <- check.and.get.cells(SO = SO,
                               assay = assay,
                               cells = cells,
                               feature_cut = feature_cut,
                               feature_cut_expr = feature_cut_expr,
                               feature_ex = feature_ex,
                               downsample = downsample,
                               included_only = T)

  data <- get_data(SO,
                   feature = features,
                   reduction = NULL,
                   assay = assay,
                   layer = "data",
                   cells = cells,
                   meta_col = meta_col,
                   split_feature = split.by)
  data <- dplyr::bind_rows(data, .id = "feature_split")

  # split.by requires testing - include in pivoting etc and geom_text
  #data <- tidyr::pivot_longer(data, dplyr::all_of(features), names_to = "feature_split", values_to = "expr")

  if (plot.expr.freq) {
    stat <-
      data %>%
      dplyr::mutate(max.feat.expr = max(feature), .by = feature_split) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("feature_split", meta_col, "SO.split", "max.feat.expr")))) %>%
      dplyr::summarise(pct.expr = sum(feature > 0)/dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(pct.expr.adjust.pct = dplyr::case_when(pct.expr == 0 ~ "0 %",
                                                           pct.expr > 0 & pct.expr < 0.01 ~ "> 1 %",
                                                           pct.expr >= 0.01 ~ paste0(round(pct.expr*100, expr.freq.decimals), " %"))) %>%
      dplyr::mutate(pct.expr.adjust = dplyr::case_when(pct.expr == 0 ~ "0",
                                                       pct.expr > 0 & pct.expr < 0.01 ~ "> 0.01",
                                                       pct.expr >= 0.01 ~ as.character(round(pct.expr, expr.freq.decimals))))
    names(stat)[which(names(stat) == "SO.split")] <- legend.title
  }
  names(data)[which(names(data) == "SO.split")] <- legend.title


  if (!plot.non.expr) {
    data <- dplyr::filter(data, feature > 0)
  }


  geom2 <- switch(geom2,
                  "violin" = ggplot2::geom_violin,
                  "boxplot" = ggplot2::geom_boxplot,
                  "none" = ggplot2::geom_blank)
  color_aes <- ifelse(length(SO) > 1, legend.title, meta_col)
  geom1pos <- if (geom1 == "jitter" && jitterwidth > 0) ggplot2::position_jitterdodge(jitter.width = jitterwidth, dodge.width = dodgewidth) else ggplot2::position_dodge(width = dodgewidth)

  col.pal <- colrr::make_col_pal(col_vec = col.pal,
                                 fct_lvls = if (is.factor(data[[color_aes]])) levels(data[[color_aes]]) else unique(data[[color_aes]]),
                                 missing_fct_to_na = ifelse("missing_fct_to_na" %in% names(col_pal_args), col_pal_args[["missing_fct_to_na"]], T),
                                 col_pal_args = col_pal_args[-which(names(col_pal_args) %in% c("name", "missing_fct_to_na"))])

  plot <-
    ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(meta_col), y = feature, color = !!rlang::sym(color_aes))) +
    ggplot2::scale_color_manual(values = col.pal, na.value = col.na) +
    theme +
    Gmisc::fastDoCall(ggplot2::facet_wrap, args = c(list(facets = ggplot2::vars(feature_split)), facetting_args))

  if (plot_first == "geom1") {
    plot <- plot +
      ggplot2::geom_point(size = pt.size, position = geom1pos) +
      suppressWarnings(Gmisc::fastDoCall(geom2, args = c(list(position = ggplot2::position_dodge(width = dodgewidth)),
                                                         geom2_args)))
  } else if (plot_first == "geom2") {
    plot <- plot +
      suppressWarnings(Gmisc::fastDoCall(geom2, args = c(list(position = ggplot2::position_dodge(width = dodgewidth)),
                                                         geom2_args))) +
      ggplot2::geom_point(size = pt.size, position = geom1pos)
  }



  ##ggforce::geom_sina()


  if (plot.expr.freq) {
    plot <- plot +
      ggplot2::geom_text(data = stat,
                         ggplot2::aes(label = !!rlang::sym(ifelse(expr.freq.pct,
                                                                  "pct.expr.adjust.pct",
                                                                  "pct.expr.adjust")),
                                      y = max.feat.expr + expr.freq.hjust),
                         position = ggplot2::position_dodge(width = 0.75),
                         size = label.size, show.legend = F)
    #expand_limits
  }

  if (length(SO) == 1) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }
  if (!plot.strip) {
    plot <- plot + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.x = ggplot2::element_blank())
  }

  return(plot)
}

