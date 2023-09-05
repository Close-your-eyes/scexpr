#' Title
#'
#' @param SO Seurat object or data frame
#' @param meta.col
#' @param label_size can one number or a vector for each label
#' @param label_inside_radius
#' @param legend.position
#' @param col_pal
#' @param order_pieces NULL to not change order, T for decreasing order of F for increasing order (with respect to size of pieces)
#' @param avoid_label_overlap strategy to elegantly avoid label overlap
#' @param border_color color of borders of pie pieces
#' @param outside_radius when avoid_label_overlap is outside; the radius where to place label;
#' 1 roughly equals the outer radius of the circle
#' @param label_outside
#' @param label_inside
#' @param label_outside_radius
#' @param label_angle
#' @param label_rel_percent
#' @param print_pct_sign
#' @param label_decimals
#' @param legend_title
#' @param theme
#' @param theme_args
#'
#' @return
#' @export
#'
#' @examples
freq_pie_chart <- function(SO,
                           meta.col,
                           pie_inside_radius = 0.3,
                           label_outside = c("none", "abs", "rel"),
                           label_inside = c("rel", "abs", "none"),
                           label_rel_cutoff = 0,
                           label_size = 5,
                           label_inside_radius = 0.75,
                           label_outside_radius = 1.1,
                           label_angle_inside = "circle", # circle or numeric
                           label_angle_outside = "circle", # circle or numeric
                           legend.position = "right",
                           col_pal = scexpr::col_pal("custom"),
                           border_color = "white",
                           order_pieces = T,
                           avoid_label_overlap = c("not", "alternating_shift", "outside"),
                           outside_radius = 1.1,
                           label_rel_percent = F,
                           print_pct_sign = T,
                           label_decimals = 2,
                           legend_title = NULL,
                           theme = ggplot2::theme_bw(),
                           theme_args = list(legend.background = ggplot2::element_blank(),
                                             legend.key.size = ggplot2::unit(0.3, "cm"),
                                             legend.key = ggplot2::element_blank(), ## keeps background of legend symbols transparent
                                             panel.grid = ggplot2::element_blank(),
                                             axis.title = ggplot2::element_blank(),
                                             axis.text = ggplot2::element_blank(),
                                             axis.ticks = ggplot2::element_blank())
                           #avoid_overlap_min_frac = 0.05
) {

  ## geom_textpath for labels?

  if (!requireNamespace("ggforce", quietly = T)) {
    utils::install.packages("ggforce")
  }
  if (!requireNamespace("farver", quietly = T)) {
    utils::install.packages("farver")
  }
  if (methods::is(SO, "Seurat")) {
    SO <- SO@meta.data
  }

  if (length(legend.position) > 2) {stop("legend.position should have length 1 being top, bottom, left, right; or length 2 indicating the corner where legend is to place.")}

  ## ... make that a list and check for default elements (see below)
  ## them call theme with Gmisc::fastDoCall

  avoid_label_overlap <- match.arg(avoid_label_overlap, c("not", "alternating_shift", "outside"))
  label_outside <- match.arg(label_outside, c("none", "abs", "rel"))
  label_inside <- match.arg(label_inside, c("rel", "abs", "none"))

  if (label_angle_inside != "circle") {
    if (!is.numeric(label_angle_inside)) {
      stop("label_angle_inside has to be 'circle' or numeric.")
    }
  }
  if (label_angle_outside != "circle") {
    if (!is.numeric(label_angle_outside)) {
      stop("label_angle_outside has to be 'circle' or numeric.")
    }
  }

  # https://stackoverflow.com/questions/16184188/ggplot-facet-piechart-placing-text-in-the-middle-of-pie-chart-slices (ggforce)
  tab <- table(SO[,meta.col,drop=T], exclude = c())
  tab <- stats::setNames(as.numeric(tab), names(tab))
  tab <- data.frame(abs = unname(tab), rel = as.numeric(tab/sum(tab)), cluster = factor(names(tab), levels = names(tab)))
  if (!is.null(order_pieces)) {
    tab <- tab[order(tab[,"rel"], decreasing = order_pieces), ]
  }
  tab$start_angle <- c(0,cumsum(tab[,"rel"]))[-(length(tab[,"rel"]) + 1)]*pi*2
  tab$end_angle <- c(cumsum(tab[,"rel"]))*pi*2
  tab$mid_angle <-  0.5*(tab[,"start_angle"] + tab[,"end_angle"])

  # text angle equal to angle of circle but readable
  if (label_angle_inside == "circle") {
    # relevant if order_pieces = T or order_pieces = F ??
    tab$text_angle_inside <- ifelse(tab$mid_angle > pi, 270 - tab$mid_angle*180/pi, 90 - tab$mid_angle*180/pi)
  } else {
    tab$text_angle_inside <- label_angle_inside
  }
  if (label_angle_outside == "circle") {
    # relevant if order_pieces = T or order_pieces = F ??
    tab$text_angle_outside <- ifelse(tab$mid_angle > pi, 270 - tab$mid_angle*180/pi, 90 - tab$mid_angle*180/pi)
  } else {
    tab$text_angle_outside <- label_angle_outside
  }

  if (length(label_outside_radius) > nrow(tab) || nrow(tab) %% length(label_outside_radius) != 0) {
    message("label_outside_radius has incompatible length. Will be set to 1.1.")
    label_outside_radius <- 1.1
  }
  if (length(label_inside_radius) > nrow(tab) || nrow(tab) %% length(label_inside_radius) != 0) {
    message("label_inside_radius has incompatible length. Will be set to 0.75.")
    label_inside_radius <- 0.75
  }
  # optional: adjust position of text labels
  tab$label_inside_radius <- label_inside_radius
  tab$label_outside_radius <- label_outside_radius
  tab[,"frac_lag"] <- as.numeric(lag(tab[,"rel"], default = 0)) # drop attributes; for rle below
  tab[,"frac_lag_diff"] <- tab[,"rel"] - tab[,"frac_lag"]

  ## not used yet
  tab[,"frac_lag_diff_series"] <- cumsum(abs(tab[,"frac_lag_diff",drop=T]) <= 0.05)

  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  rel_series <- rle(abs(tab[,"frac_lag_diff"]) <= 0.05)
  if (avoid_label_overlap == "alternating_shift") {
    tab$label_inside_radius[unlist(seq2(lag(cumsum(rel_series$lengths)+1)[-1], cumsum(rel_series$lengths)[-1]))] <-
      ifelse(unlist(seq2(lag(cumsum(rel_series$lengths)+1)[-1], cumsum(rel_series$lengths)[-1])) %% 2 == 0,
             tab$label_inside_radius[unlist(seq2(lag(cumsum(rel_series$lengths)+1)[-1], cumsum(rel_series$lengths)[-1]))] - 0.1,
             tab$label_inside_radius[unlist(seq2(lag(cumsum(rel_series$lengths)+1)[-1], cumsum(rel_series$lengths)[-1]))] + 0.1)
  } else if (avoid_label_overlap == "outside") {
    tab$label_inside_radius[unlist(seq2(lag(cumsum(rel_series$lengths)+1)[-1], cumsum(rel_series$lengths)[-1]))] <- outside_radius
  }

  if (length(col_pal) != nlevels(tab[,"cluster"])) {

    if (is.null(names(col_pal))) {
      if (length(col_pal) < length(unique(tab[,"cluster"]))) {
        col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
        warning("Number of colors provided not sufficient for number of factor levels. Falling back to scales::hue_pal().")
      } else {
        col_pal <- col_pal[1:length(unique(tab[,"cluster"]))]
      }
      names(col_pal) <- tab$cluster_cols
    } else {
      if (length(col_pal) > length(unique(tab[,"cluster"])) && all(names(col_pal) %in% unique(tab[,"cluster"]))) {
        col_pal <- col_pal[unique(as.character(tab[,"cluster"]))]
        tab$cluster_cols <- col_pal[tab$cluster]
      } else {
        warning("Number of colors provided not matching the number of factor levels in meta.col. Falling back to scales::hue_pal().")
        col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
        names(col_pal) <- tab$cluster_cols
      }
    }

  } else {
    if (!is.null(names(col_pal)) && !all(names(col_pal) %in% unique(tab[,"cluster"]))) {
      warning("Not all names of col_pal found in factor levels of meta.col. Falling back to scales::hue_pal().")
      col_pal <- scales::hue_pal()(length(unique(tab[,"cluster"])))
      names(col_pal) <- tab$cluster_cols
    }
  }
  tab$cluster_cols <- col_pal[as.character(tab[,"cluster"])]

  plot <-
    ggplot2::ggplot(tab, ggplot2::aes(x0 = 0, y0 = 0, r0 = pie_inside_radius, r = 1, start = start_angle, end = end_angle, fill = cluster)) +
    ggforce::geom_arc_bar(colour = border_color) +
    theme +
    labs(fill = legend_title) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::coord_fixed(ratio = 1)

  if (length(legend.position) == 1) {
    plot <-
      plot +
      Gmisc::fastDoCall(ggplot2::theme, args = c(theme_args, list(legend.position = legend.position)))
  } else {
    plot <-
      plot +
      Gmisc::fastDoCall(ggplot2::theme, args = c(theme_args, list(legend.justification = c(legend.position[1], legend.position[2]),
                                                                  legend.position = c(legend.position[1], legend.position[2]))))
  }

  if (is.na(plot[["theme"]][["panel.background"]][["fill"]])) {
    plot <- plot + theme(panel.background = element_rect(fill = "white"))
  }

  # in case label are outside of the circle, adjust their color (black or white) to the panel.background
  tab$text_color_inside <- ifelse(farver::decode_colour(tab$cluster_cols, to = "hcl")[, "l"] > 50, "black", "white")
  tab[which(tab$label_inside_radius >= 1), "text_color_inside"] <- ifelse(farver::decode_colour(plot[["theme"]][["panel.background"]][["fill"]], to = "hcl")[,"l"] > 50, "black", "white")

  tab$text_color_outside <- ifelse(farver::decode_colour(tab$cluster_cols, to = "hcl")[, "l"] > 50, "black", "white")
  tab[which(tab$label_outside_radius >= 1), "text_color_outside"] <- ifelse(farver::decode_colour(plot[["theme"]][["panel.background"]][["fill"]], to = "hcl")[,"l"] > 50, "black", "white")

  if (label_inside == "rel") {
    if (label_rel_percent) {
      ## problem with decimals and > 1 % may arise
      temp_labels <- round(tab[,"rel"]*100, label_decimals)
      temp_labels[which(temp_labels == 0 & tab[,"rel"] > 0)] <- "< 1"
      tab$label_inside_text <- temp_labels
      if (any(tab$rel < label_rel_cutoff)) {
        tab$label_inside_text[which(tab$rel < label_rel_cutoff)] <- ""
      }
      if (print_pct_sign) {
        for (i in 1:length(tab$label_inside_text)) {
          if (tab$label_inside_text[i] != "") {
            if (tab$rel[i] < 0.01 & tab$rel[i] > 0) {
              tab$label_inside_text[i] <- "< 1 %"
            } else if (tab$rel[i] > 0.01 & tab$rel[i] < 0.99) {
              tab$label_inside_text[i] <- paste0(tab$label_inside_text[i] , " %")
            } else if (tab$rel[i] > 0.99 & tab$rel[i] < 1) {
              tab$label_inside_text[i] <- "> 99 %"
            } else {
              tab$label_inside_text[i] <- paste0(tab$label_inside_text[i], " %")
            }
          }
        }
      } else {
        # or use sprinf fun
        tab$label_inside_text <- format(tab$label_inside_text, nsmall = label_decimals)
      }
    } else {
      tab$label_inside_text <- format(round(tab[,"rel"], label_decimals), nsmall = label_decimals)
    }
  } else if (label_inside == "abs") {
    tab$label_inside_text <- tab[,"abs"]
  }

  if (label_outside == "rel") {
    if (label_rel_percent) {
      ## problem with decimals and > 1 % may arise
      #tab$label_outside_text <- format(round(tab[,"rel",drop=T]*100, label_decimals), nsmall = label_decimals)
      temp_labels <- round(tab[,"rel"]*100, label_decimals)
      temp_labels[which(temp_labels == 0 & tab[,"rel"] > 0)] <- "< 1"
      tab$label_outside_text <- temp_labels
      if (any(tab$rel < label_rel_cutoff)) {
        tab$label_outside_text[which(tab$rel < label_rel_cutoff)] <- ""
      }
      if (print_pct_sign) {
        for (i in 1:length(tab$label_outside_text)) {
          if (tab$label_outside_text[i] != "") {
            if (tab$rel[i] < 0.01 & tab$rel[i] > 0) {
              tab$label_outside_text[i] <- "< 1 %"
            } else if (tab$rel[i] > 0.01 & tab$rel[i] < 0.99) {
              tab$label_outside_text[i] <- paste0(tab$label_outside_text[i] , " %")
            } else if (tab$rel[i] > 0.99 & tab$rel[i] < 1) {
              tab$label_outside_text[i] <- "> 99 %"
            } else {
              tab$label_outside_text[i] <- paste0(tab$label_outside_text[i], " %")
            }
          }
        }
      } else {
        tab$label_outside_text <- format(tab$label_outside_text, nsmall = label_decimals)
      }
    } else {
      tab$label_outside_text <- format(round(tab[,"rel"], label_decimals), nsmall = label_decimals)
    }
  } else if (label_outside == "abs") {
    tab$label_outside_text <- tab[,"abs"]
  }

  if (label_inside != "none") {
    plot <-
      plot +
      ggplot2::geom_text(data = tab, ggplot2::aes(color = I(text_color_inside),
                                                  x = label_inside_radius*sin(mid_angle),
                                                  y = label_inside_radius*cos(mid_angle),
                                                  angle = text_angle_inside,
                                                  label = label_inside_text),
                         size = label_size)
  }

  if (label_outside != "none") {
    plot <-
      plot +
      ggplot2::geom_text(data = tab, ggplot2::aes(color = I(text_color_outside),
                                                  x = label_outside_radius*sin(mid_angle),
                                                  y = label_outside_radius*cos(mid_angle),
                                                  angle = text_angle_outside,
                                                  label = label_outside_text),
                         size = label_size,
                         hjust = 1)
  }

  return(list(plot = plot, data = tab))
}
