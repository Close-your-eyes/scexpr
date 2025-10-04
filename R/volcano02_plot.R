#' Make volcano plot
#'
#' Calculate with volcano01_calc. Plot with volcano02_plot.
#'
#' @param volc01_df data frame from volcano01_calc
#' @param x column from volc01_df for x axis
#' @param y column from volc01_df for y axis
#' @param minus_log10_y make y = -log10(y) (when p-value kind of thing)
#' @param feature column from volc01_df with feature (gene) names
#' @param x_symm make x axis symmetric
#' @param y_pseudo_log make y axis pseudo log with scales::pseudo_log_trans
#' @param pseudo_log_sigma sigma for pseudo log
#' @param features_to_color which features to plot with color on top
#' @param features_color_by which column to color them by;
#' numeric --> continuous, else: discrete
#' @param errorbar_low column from volc01_df with low errorbar size; only when
#' MetaVolcanoR
#' @param errorbar_up see errorbar_low
#' @param errorbar_size to ggplot2::geom_errorbar
#' @param errorbar_width ggplot2::geom_errorbar
#' @param col_pal color palette to colrr::col_pal
#' @param col_pal_dir to direction in colrr::col_pal
#' @param pt_size point size
#' @param pt_alpha point alpha
#' @param pt_col point color
#' @param pt_col_inf point color for those with Inf log fc
#' @param x_expnd expand (or reduce) x-axis by this at top and bottom
#' @param p_tick value for extra tick on y-axis (~p-value)
#' @param feature_exclude exclude features from plotting, vector
#' @param min_pct filter features by min_pct
#' @param cut_p provide cut_p and cut_fc and count number of genes above these
#' limits
#' @param cut_fc see cut_p
#' @param label_plot plot feature labels?
#' @param label_features select features for labeling; or pass "significant", then
#' all features below label_p_signif are labelled
#' @param label_topn_metric select features to label by top values in x and/or y
#' @param label_topn how many of top hits in metric? (roughly)
#' @param label_pt_col extra color for labelled points; NULL for no extra color
#' @param label_neg_pos_sep label points >0 and <0 separately? useful for uniform
#' repelling to left and right
#' @param label_p_signif p value (y-axis) cutoff for labeling, passed
#' "significant" to label_features to enable this
#' @param geom_text_repel_args arguments to ggrepel::geom_text_repel
#' @param theme ggplot theme
#'
#' @returns ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#'   volc <- scexpr::volcano01_calc(
#'     so2,
#'     neg_cells = c(15,16),
#'     pos_cells = c(17,18),
#'     method = "wilcox")
#'   volcano02_plot(volc[["df"]])
#' }
volcano02_plot <- function(volc01_df,
                           x = "avg_log2FC",
                           y = "p_val_adj",
                           feature = "feature",
                           pt_col = "grey60",
                           pt_col_inf = "cornflowerblue",
                           pt_size = 1,
                           pt_alpha = 0.8,
                           minus_log10_y = T,
                           x_symm = T,
                           x_expnd = 0.5,
                           y_pseudo_log = F,
                           pseudo_log_sigma = 1,
                           features_to_color = NULL,
                           features_color_by = NULL,
                           errorbar_low = NULL, # absolute coordinate of lower errorbar
                           errorbar_up = NULL, # absolute coordinate of upper errorbar
                           errorbar_size = 0.2,
                           errorbar_width = 0.2,
                           col_pal = "RdBu",
                           col_pal_dir = -1,
                           theme = colrr::theme_material(base_size = 14,
                                                         white = T),
                           p_tick = NULL,
                           feature_exclude = NULL,
                           min_pct = 0,
                           cut_p = NA,
                           cut_fc = NA,

                           label_plot = T,
                           label_features = NULL,
                           label_topn_metric = c("x", "y"),
                           label_topn = 30,
                           label_pt_col = "tomato2",
                           label_neg_pos_sep = F,
                           label_p_signif = 0.001,

                           geom_text_repel_args = list(fontface = "italic",
                                                       max.overlaps = 50)) {

  if (!requireNamespace("colrr", quietly = T)) {
    devtools::install_github("Close-your-eyes/colrr")
  }

  vd <- as.data.frame(volc01_df)

  x <- match.arg(x, colnames(vd))
  y <- match.arg(y, colnames(vd)) #c("adj.p.val", "p.val")

  label_topn_metric <- rlang::arg_match(label_topn_metric, multiple = T)
  att <- attributes(vd)


  if (!is.null(feature_exclude)) {
    # grepl(paste(feature_exclude, collapse = "|"), vd[,feature])
    print(paste0("The following features are excluded from the volcano plot: ", paste(grep(paste(feature_exclude, collapse = "|"), vd[,feature], value = T), collapse = ",")))
    vd <- vd[which(!vd[,feature] %in% feature_exclude),]
  }

  if (!is.null(features_to_color)) {
    if (any(!features_to_color %in% vd[,feature])) {
      message("features_to_color: ", paste(features_to_color[which(!features_to_color %in% vd[,feature])], collapse = ", "), " not found.")
    }
    features_to_color <- features_to_color[which(features_to_color %in% vd[,feature])]
    if (length(features_to_color) == 0) {
      features_to_color <- NULL
    }
  }

  if (!is.null(features_color_by) && !features_color_by %in% names(vd)) {
    message("features_color_by not found as column.")
    features_color_by <- NULL
  }

  if (!is.null(errorbar_low) && !errorbar_low %in% names(vd)) {
    message("errorbar_low not found as column. Both errorbar limit will be ignored.")
    errorbar_low <- NULL
    errorbar_up <- NULL
  }

  if (!is.null(errorbar_up) && !errorbar_up %in% names(vd)) {
    message("errorbar_up not found as column. Both errorbar limit will be ignored.")
    errorbar_low <- NULL
    errorbar_up <- NULL
  }

  if (min_pct > 0) {
    vd <- rbind(vd[intersect(which(vd[[paste0("pct.", att$neg_name)]] >= min_pct), which(vd[,x] < 0)),],
                vd[intersect(which(vd[[,paste0("pct.", att$pos_name)]] >= min_pct), which(vd[,x] > 0)),])
  }

  if (minus_log10_y) {
    vp <- ggplot2::ggplot(vd, ggplot2::aes(x = !!rlang::sym(x), y = -log10(!!rlang::sym(y))))
  } else {
    vp <- ggplot2::ggplot(vd, ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y)))
  }
  vp <-
    vp +
    ggplot2::geom_point(color = pt_col, alpha = pt_alpha, size = pt_size) +
    ggplot2::geom_point(data = vd[which(vd$infinite.FC == 1),], color = pt_col_inf, size = pt_size) +
    theme

  if (!is.null(features_color_by)) {

    if (is.null(features_to_color)) {
      features_to_color <- unique(vd[[feature]])
    }
    if (is.numeric(vd[[features_color_by]])) {
      # continuous
      fct_lvls <- NULL
    } else {
      fct_lvls <- vd[which(vd[,feature] %in% features_to_color),features_color_by]
      vd[[features_color_by]] <- as.factor(vd[[features_color_by]])
    }

    col_pal <- make_col_pal(col_vec = col_pal,
                            fct_lvls = fct_lvls,
                            col_pal_args = list(direction = col_pal_dir))
    vp <- vp +
      ggplot2::geom_point(
        data = vd[which(vd[[feature]] %in% features_to_color),],
        ggplot2::aes(color = !!rlang::sym(features_color_by)),
        size = pt_size)

    if (is.numeric(vd[[features_color_by]])) {
      vp <- vp + ggplot2::scale_color_gradientn(colors = col_pal)
    } else {
      vp <- vp + ggplot2::scale_color_manual(values = col_pal)
    }

    if (!is.null(errorbar_up)) {
      # checking one of errorbar_up, errorbar_low is enough
      vp <- vp +
        ggplot2::geom_errorbar(data = vd[which(vd[,feature] %in% features_to_color),],
                               ggplot2::aes(color = !!rlang::sym(features_color_by),
                                            xmin = !!rlang::sym(errorbar_low),
                                            xmax = !!rlang::sym(errorbar_up)),
                               size = errorbar_size, width = errorbar_width)
      ## 95 % conf-interval in case of metavolcanoR
    }
  }

  if (!is.null(p_tick) && p_tick > 0) {
    gg.brk <- ggplot2::ggplot_build(vp)[["layout"]][["panel_params"]][[1]][["y"]][["breaks"]]
    ord <- order(c(gg.brk, -log10(p_tick)))
    brk <- c(gg.brk, -log10(p_tick))
    lab <- c(gg.brk, paste0("p = ", p_tick))
    lab <- gsub("^0$", "", lab)
    brk <- brk[ord]
    lab <- lab[ord]
  } else {
    brk <- ggplot2::waiver()
    lab <- ggplot2::waiver()
  }

  if (y_pseudo_log) {
    vp <- vp + ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma),
                                           breaks = brk, labels = lab, limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  } else {
    vp <- vp + ggplot2::scale_y_continuous(breaks = brk, labels = lab, limits = vp[["layout"]][["panel_params"]][[1]][["y"]][["limits"]])
  }

  if (minus_log10_y) {
    vp <- vp + ggplot2::labs(x = bquote(bold(.(att$neg_name)) ~ " <-- " ~ log[2] ~ "FC" ~ " --> " ~ bold(.(att$pos_name))), y = bquote(-log[10]~.(rlang::sym(y))))
  } else {
    vp <- vp + ggplot2::labs(x = bquote(bold(.(att$neg_name)) ~ " <-- " ~ log[2] ~ "FC" ~ " --> " ~ bold(.(att$pos_name))), y = y)
  }

  if (x_symm) {
    vp <- vp + ggplot2::xlim(-round(max(abs(vd[,x]))) - x_expnd, round(max(abs(vd[,x]))) + x_expnd)
  } else {
    vp <- vp + ggplot2::xlim(round(min(vd[,x])) - x_expnd, round(max(vd[,x])) + x_expnd)
  }


  if (!any(c(is.na(cut_p), is.na(cut_fc)))) {
    vp <- vp + ggplot2::geom_hline(yintercept = -log10(cut_p), linetype = "dashed") + ggplot2::geom_vline(xintercept = c(-cut_fc, cut_fc), linetype = "dashed")

    vd.cut <-
      vd |>
      dplyr::filter(abs(!!rlang::sym(x)) >= cut_fc) |>
      dplyr::filter(!!rlang::sym(y) <= cut_p)

    vd.cut.sum <-
      vd.cut |>
      dplyr::mutate(sign = ifelse(!!rlang::sym(x) > 0, "plus", "minus")) |>
      dplyr::group_by(sign) |>
      dplyr::summarise(n.genes = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(xpos = ifelse(sign == "plus", cut_fc + 0.75*(max(abs(vd.cut[,x]))-cut_fc), -(cut_fc + 0.75*(max(abs(vd.cut[,x]))-cut_fc)))) |>
      dplyr::mutate(ypos = -log10(cut_p)*2)

    vp <-
      vp +
      ggplot2::geom_point(data = vd.cut, color = "black", size = pt_size) +
      ggplot2::geom_text(data = vd.cut.sum, ggplot2::aes(x = xpos, y = ypos, label = n.genes), inherit.aes = F, size = 5/14*base_size)

  }


  #### label section

  if (!is.null(feature_exclude)) {
    vd <- vd[which(!grepl(paste0(feature_exclude, collapse = "|"), vd[,feature])),]
  }

  if (is.null(label_features)) {
    if (length(label_topn_metric) == 1 && label_topn_metric == "y") {
      f_lab <- dplyr::top_n(vd, -label_topn, !!rlang::sym(y))
    } else if (length(label_topn_metric) == 2) {
      f_lab.p.val <- dplyr::top_n(vd, -label_topn, !!rlang::sym(y))
      f_lab.logfc <- dplyr::bind_rows(
        dplyr::top_n(vd, label_topn/2, !!rlang::sym(x)),
        dplyr::top_n(vd, -(label_topn/2), !!rlang::sym(x))
      )
      f_lab <- dplyr::bind_rows(f_lab.logfc, f_lab.p.val) |> dplyr::distinct()
    } else if (label_topn_metric == 1 && label_topn_metric == "x") {
      f_lab <- dplyr::bind_rows(
        dplyr::top_n(vd, label_topn/2, !!rlang::sym(x)),
        dplyr::top_n(vd, -(label_topn/2), !!rlang::sym(x))
      )
    }

  } else {

    if (length(label_features) == 1 && label_features == "significant") {
      label_features <- vd[which(as.numeric(vd[,y]) < label_p_signif), feature]
    }
    f_lab <- dplyr::filter(vd, !!rlang::sym(feature) %in% label_features)

  }

  if (!is.null(label_pt_col)) {
    vp <- vp + ggplot2::geom_point(data = f_lab, color = label_pt_col)
  }

  if (label_plot) {
    if (label_neg_pos_sep) {
      # only handle nudge_x automatically
      if (!"nudge_x" %in% names(geom_text_repel_args)) {
        geom_text_repel_args[["nudge_x"]] <- 0
      }
      geom_text_repel_args[["nudge_x"]] <- abs(geom_text_repel_args[["nudge_x"]])
      vp <- vp +
        Gmisc::fastDoCall(what = ggrepel::geom_text_repel,
                          args = c(list(mapping = ggplot2::aes(label = !!rlang::sym(feature)), data = dplyr::filter(f_lab, !!rlang::sym(x) > 0)),
                                   geom_text_repel_args))
      geom_text_repel_args[["nudge_x"]] <- -1*geom_text_repel_args[["nudge_x"]]
      vp <- vp +
        Gmisc::fastDoCall(what = ggrepel::geom_text_repel,
                          args = c(list(mapping = ggplot2::aes(label = !!rlang::sym(feature)), data = dplyr::filter(f_lab, !!rlang::sym(x) < 0)),
                                   geom_text_repel_args))
    } else {
      vp <- vp +
        Gmisc::fastDoCall(what = ggrepel::geom_text_repel,
                          args = c(list(mapping = ggplot2::aes(label = !!rlang::sym(feature)), data = f_lab),
                                   geom_text_repel_args))
    }
  }

  return(vp)
}
