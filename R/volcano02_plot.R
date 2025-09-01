#' Title
#'
#' @param vd
#' @param x
#' @param y
#' @param minus_log10_y
#' @param squish_max_zero_p
#' @param feature
#' @param x_label
#' @param y_label
#' @param x_symm
#' @param y_pseudo_log
#' @param pseudo_log_sigma
#' @param features.to.color
#' @param features.color.by
#' @param errorbar.low.col
#' @param errorbar.up.col
#' @param errorbar.size
#' @param errorbar.width
#' @param col.pal
#' @param col.pal.dir
#' @param col.type
#' @param pt_size
#' @param pt_alpha
#' @param font.size
#' @param font.family
#' @param ngn
#' @param pgn
#' @param x_expnd
#' @param p_tick
#' @param feature_exclude
#' @param min_pct
#' @param cut_p
#' @param cut_fc
#'
#' @returns
#' @export
#'
#' @examples
volcano02_plot <- function(vd,
                           x = "log2.fc",
                           y = "adj.p.val",
                           feature = "feature",
                           pt_size = 1,
                           pt_alpha = 0.8,
                           minus_log10_y = T,
                           x_symm = T,
                           x_expnd = 0.5,
                           y_pseudo_log = F,
                           pseudo_log_sigma = 1,
                           features.to.color = NULL, #which features to plot with color on top
                           features.color.by = NULL, #which column to color them by; numeric --> continuous, else: discrete
                           errorbar.low.col = NULL, # absolute coordinate of lower errorbar
                           errorbar.up.col = NULL, # absolute coordinate of upper errorbar
                           errorbar.size = 0.2,
                           errorbar.width = 0.2,
                           col.pal = "RdBu",
                           col.pal.dir = 1,
                           col.type = c("c", "d"), # continuous or discrete
                           font.size = 14,
                           p_tick = NULL,
                           feature_exclude = NULL,
                           min_pct = 0,
                           cut_p = NA,
                           cut_fc = NA,

                           label_plot = T,
                           label_features = NULL,
                           label_topn_metric = c("p_val", "fc", "both"),
                           label_topn = 30,
                           label_dot_color = "tomato2", # set NA to have no color
                           label_neg_pos_sep = F,
                           label_p_signif = 0.001,

                           geom_text_repel_args = list(fontface = "italic",
                                                       max.overlaps = 50)) {



  x <- match.arg(x, colnames(vd))
  y <- match.arg(y, colnames(vd)) #c("adj.p.val", "p.val")
  col.type <- match.arg(col.type, c("c", "d"))
  label_topn_metric <- rlang::arg_match(label_topn_metric)
  att <- attributes(vd)


  if (!is.null(feature_exclude)) {
    # grepl(paste(feature_exclude, collapse = "|"), vd[,feature])
    print(paste0("The following features are excluded from the volcano plot: ", paste(grep(paste(feature_exclude, collapse = "|"), vd[,feature], value = T), collapse = ",")))
    vd <- vd[which(!vd[,feature] %in% feature_exclude),]
  }

  if (!is.null(features.to.color)) {
    if (any(!features.to.color %in% vd[,feature])) {
      message("features.to.color: ", paste(features.to.color[which(!features.to.color %in% vd[,feature])], collapse = ", "), " not found.")
    }
    features.to.color <- features.to.color[which(features.to.color %in% vd[,feature])]
    if (length(features.to.color) == 0) {
      features.to.color <- NULL
    }
  }

  if (!is.null(features.color.by) && !features.color.by %in% names(vd)) {
    message("features.color.by not found as column.")
    features.color.by <- NULL
  }

  if (!is.null(errorbar.low.col) && !errorbar.low.col %in% names(vd)) {
    message("errorbar.low.col not found as column. Both errorbar limit will be ignored.")
    errorbar.low.col <- NULL
    errorbar.up.col <- NULL
  }

  if (!is.null(errorbar.up.col) && !errorbar.up.col %in% names(vd)) {
    message("errorbar.up.col not found as column. Both errorbar limit will be ignored.")
    errorbar.low.col <- NULL
    errorbar.up.col <- NULL
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
    ggplot2::geom_point(color = "#999999", alpha = pt_alpha, size = pt_size) +
    ggplot2::geom_point(data = vd[which(vd$infinite.FC == 1),], color = "cornflowerblue", size = pt_size) +
    ggplot2::theme_bw(base_size = font.size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


  if (!is.null(features.color.by) && !is.null(features.to.color)) {
    # select color scale
    if (col.type == "c") {
      if (length(col.pal) == 1 && !col.pal %in% grDevices::colors()) {
        col.pal <- colrr::col_pal(name = col.pal, direction = col.pal.dir)
      }
    } else if (col.type == "d") {
      if (length(col.pal) == 1 && !col.pal %in% grDevices::colors()) {
        col.pal <- colrr::col_pal(name = col.pal, direction = col.pal.dir, n = nlevels(as.factor(vd[which(vd[,feature] %in% features.to.color),features.color.by])))
      }
    }
    if (col.type == "c") {
      vp <-
        vp +
        ggplot2::geom_point(data = vd[which(vd[,feature] %in% features.to.color),], ggplot2::aes(color = !!rlang::sym(features.color.by)), size = pt_size) +
        ggplot2::scale_color_gradientn(colors = col.pal)
    }
    if (col.type == "d") {
      vd[,features.color.by] <- as.factor(vd[,features.color.by])
      vp <-
        vp +
        ggplot2::geom_point(data = vd[which(vd[,feature] %in% features.to.color),], ggplot2::aes(color = !!rlang::sym(features.color.by)), size = pt_size) +
        ggplot2::scale_color_manual(values = col.pal)
    }

    if (!is.null(errorbar.up.col)) {
      # checking one of errorbar.up.col, errorbar.low.col is enough
      vp <-
        vp +
        ggplot2::geom_errorbar(data = vd[which(vd[,feature] %in% features.to.color),], ggplot2::aes(color = !!rlang::sym(features.color.by),
                                                                                                    xmin = !!rlang::sym(errorbar.low.col),
                                                                                                    xmax = !!rlang::sym(errorbar.up.col)),
                               size = errorbar.size, width = errorbar.width)
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
      ggplot2::geom_text(data = vd.cut.sum, ggplot2::aes(x = xpos, y = ypos, label = n.genes), inherit.aes = F, size = 5/14*font.size)

  }


  #### label section

  if (!is.null(feature_exclude)) {
    vd <- vd[which(!grepl(paste0(feature_exclude, collapse = "|"), vd[,feature])),]
  }

  if (is.null(label_features)) {
    if (label_topn_metric == "p_val") {
      f_lab <- dplyr::top_n(vd, -label_topn, !!rlang::sym(y))
    } else if (label_topn_metric == "both") {
      f_lab.p.val <- vdplyr::top_n(vd, -label_topn, !!rlang::sym(y))
      f_lab.logfc <- dplyr::bind_rows(dplyr::top_n(vd, label_topn/2, !!rlang::sym(x)), dplyr::top_n(vd, -(label_topn/2), !!rlang::sym(x)))
      f_lab <- dplyr::bind_rows(f_lab.logfc, f_lab.p.val) |> dplyr::distinct()
    } else if (label_topn_metric == "fc") {
      f_lab <- dplyr::bind_rows(dplyr::top_n(vd, label_topn/2, !!rlang::sym(x)), dplyr::top_n(vd, -(label_topn/2), !!rlang::sym(x)))
    }

  } else {

    if (length(label_features) == 1 && label_features == "significant") {
      label_features <- vd[which(as.numeric(vd[,y]) < label_p_signif), feature]
    }
    f_lab <- dplyr::filter(vd, !!rlang::sym(feature) %in% label_features)

  }


  vp <- vp + ggplot2::geom_point(data = f_lab, color = label_dot_color)
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
