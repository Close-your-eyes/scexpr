#' Plot composition of one category by another
#'
#' @param SO Seurat object or data frame, meta.data of SO is used
#' @param x_cat column to use as x axis
#' @param fill_cat column to use fill category for bar segments
#' @param color color of bars
#' @param plot_rel_labels plot relative labels of bar segment sizes
#' @param min_label_freq minimum frequency per x_cat to plot labels
#' @param label_rel_pct plot relative label as percent (TRUE) or as fraction (FALSE);
#' this also affects the axis if y = "rel"
#' @param label_rel_decimals number of decimal places for percent values
#' @param col_pal color palette for bar filling
#' @param plot_abs_labels plot absolute labels of bar segment sizes
#' @param summarize_all_x add an extra column that summarizes (averages) all others
#' @param label_rel_nudge list of x and y values to nudge relative labels;
#' list may have names of x_cat but names only apply if there is max. 1 label per x_cat,
#' e.g. when label_only_largest = T
#' @param label_abs_nudge see label_rel_nudge
#' @param flip flip x and y axis
#' @param y use absolute (n) or relative (%) y axis
#' @param plot_total_rel_labels plot total labels on top of bar showing how much the
#' x_cat occupies from total
#' @param plot_total_abs_labels see plot_total_rel_labels
#' @param label_size_rel vector of size for relative labels
#' @param label_size_abs see label_size_rel
#' @param label_total_rel_nudge x and y nudging of total label on top
#' @param label_total_abs_nudge see label_total_rel_nudge
#' @param label_position overwrite calculated label positions, which by default are centered in bar segments
#' @param label_only_largest label only the largest bar segment per x_cat,
#' e.g if there are two groups only per x_cat
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
composition_barplot <- function(SO,
                                x_cat, #x_var
                                fill_cat, #fill_var
                                y = c("rel", "abs"),
                                col_pal = colrr::col_pal("custom"),
                                color = "black",
                                plot_rel_labels = F,
                                label_only_largest = F,
                                plot_abs_labels = F,
                                plot_total_rel_labels = F,
                                plot_total_abs_labels = F,
                                #plot_sigma = F,
                                #sigma_just = ifelse(flip, 2, -1),
                                min_label_freq = 0.05,
                                label_rel_pct = T,
                                label_rel_decimals = 2,
                                summarize_all_x = F,
                                label_rel_nudge = list(c(0,0)),
                                label_abs_nudge = list(c(0,0)),
                                label_size_rel = 4,
                                label_size_abs = 4,
                                label_total_rel_nudge = c(0,0),
                                label_total_abs_nudge = c(0,0),
                                label_position = NULL,
                                flip = F) {


  if (methods::is(SO, "Seurat")) {
    SO <- SO@meta.data
  }

  if (!x_cat %in% names(SO)) {
    stop("x_cat not found in SO.")
  }
  if (!fill_cat %in% names(SO)) {
    stop("fill_cat not found in SO.")
  }

  if (label_rel_pct) {
    fctr <- 100
  } else {
    fctr <- 1
  }

  if (missing(label_rel_decimals)) {
    if (label_rel_pct) {
      label_rel_decimals <- 0
    } else {
      label_rel_decimals <- 2
    }
  }

  y = match.arg(y, c("rel", "abs"))

  table <-
    SO %>%
    dplyr::count(!!rlang::sym(x_cat), !!rlang::sym(fill_cat)) %>%
    dplyr::left_join(dplyr::count(SO, !!rlang::sym(x_cat), name = "x_total"), by = x_cat) %>%
    dplyr::left_join(dplyr::count(SO, !!rlang::sym(fill_cat), name = "fill_total"), by = fill_cat) %>%
    dplyr::mutate(rel_x = n/x_total) %>%
    dplyr::mutate(rel_fill = n/fill_total) # rel_fill is not used below

  if (summarize_all_x) {
    table <- dplyr::bind_rows(table,
                              SO %>%
                                dplyr::count(!!rlang::sym(fill_cat)) %>%
                                dplyr::mutate(x_total = nrow(SO)) %>%
                                dplyr::mutate(fill_total = nrow(SO)) %>%
                                dplyr::mutate(rel_x = n/x_total) %>%
                                dplyr::mutate(rel_fill = n/fill_total) %>%
                                dplyr::mutate(!!x_cat := "all"))
  }

  table0 <- NULL
  if (plot_total_rel_labels || plot_total_abs_labels) {
    table0 <-
      SO %>%
      dplyr::count(!!rlang::sym(x_cat)) %>%
      dplyr::mutate(pct = n/sum(n))

    if (summarize_all_x) {
      temp <- data.frame(x_cat = "all", n = nrow(SO), pct = 1)
      names(temp)[1] <- x_cat
      table0 <- dplyr::bind_rows(table0, temp)
    }
    if (y == "abs") {
      table0$total_labels_ypos <- max(table0$n)*1.05
    }
    if (y == "rel") {
      if (label_rel_pct) {
        table0$total_labels_ypos <- 108
      } else {
        table0$total_labels_ypos <- 1.08
      }
    }

    table0$pct2 <- paste0(brathering::round2(table0$pct*100, label_rel_decimals), " %")
    table0$pct_round <- brathering::round2(table0$pct, label_rel_decimals)
    table0$label_color = "black" # could be made relative to background color
  }

  # check this somewhen; in some cases label_color is not defined and then an error occurs below
  if (!is.null(names(col_pal)) && any(!unique(table[,fill_cat,drop=T]) %in% names(col_pal))) {
    col_pal <- unname(col_pal)
  } else {
    # new
    if (is.null(names(col_pal))) {
      names(col_pal) <- sort(unique(table[,fill_cat,drop=T]))
    }
    table$bar_segment_cols <- col_pal[as.character(table[,fill_cat,drop=T])]
    table$label_color <- ifelse(farver::decode_colour(table$bar_segment_cols, to = "hcl")[, "l"] > 50, "black", "white")
  }

  ## use cumsum approach to define label position in middle of bars
  ## this allow to remove label by min_label_freq and still have the correct coordinate
  ## coordinates became wrong with position = ggplot2::position_stack(vjust = 0.5) when some values were removed
  table <-
    table %>%
    dplyr::group_by(!!rlang::sym(x_cat)) %>%
    dplyr::arrange(dplyr::desc(!!rlang::sym(fill_cat)), .by_group = T) %>%
    dplyr::mutate(rel_x_cumsum = cumsum(rel_x)) %>%
    dplyr::mutate(rel_x_cumsum_lag = dplyr::lag(rel_x_cumsum, default = 0)) %>%
    dplyr::mutate(label_ypos = rel_x_cumsum_lag + (rel_x_cumsum-rel_x_cumsum_lag)/2) %>%
    dplyr::mutate(n_cumsum = cumsum(n)) %>%
    dplyr::mutate(n_cumsum_lag = dplyr::lag(n_cumsum, default = 0)) %>%
    dplyr::mutate(n_label_ypos = n_cumsum_lag + (n_cumsum-n_cumsum_lag)/2) %>%
    dplyr::mutate(rel_x_fctr = rel_x*fctr) %>%
    dplyr::mutate(rel_x_fctr_round = brathering::round2(rel_x_fctr, label_rel_decimals)) %>%
    dplyr::mutate(rel_x_fctr_pct = paste0(brathering::round2(rel_x_fctr, label_rel_decimals), " %")) %>%
    dplyr::mutate(label_ypos_fctr = label_ypos*fctr) %>%
    tibble::as_tibble()

  if (!is.null(levels(SO[,x_cat]))) {
    if (summarize_all_x) {
      table[,x_cat] <- factor(table[,x_cat,drop=T], levels = c("all", levels(SO[,x_cat])))
    } else {
      table[,x_cat] <- factor(table[,x_cat,drop=T], levels = levels(SO[,x_cat]))
    }

  }

  if (y == "rel") {
    plot <-
      ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(x_cat), y = rel_x*fctr, fill = !!rlang::sym(fill_cat))) +
      ggplot2::geom_col(color = color) +
      ggplot2::scale_fill_manual(values = col_pal) +
      ggplot2::labs(y = ifelse(label_rel_pct, "frequency [%]", "frequency")) +
      ggplot2::scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)*fctr)
  } else if (y == "abs") {
    plot <-
      ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(x_cat), y = n, fill = !!rlang::sym(fill_cat))) +
      ggplot2::geom_col(color = color) +
      ggplot2::scale_fill_manual(values = col_pal)
  }


  if (plot_rel_labels || plot_abs_labels) {
    table_temp <- dplyr::filter(table, rel_x >= min_label_freq)
    if (nrow(table_temp) < nrow(table)) {
      message(nrow(table) - nrow(table_temp), " text labels removed due to min_label_freq.")
    }
    if (label_only_largest) {
      table_temp <-
        table_temp %>%
        dplyr::group_by(!!rlang::sym(x_cat)) %>%
        dplyr::slice_max(rel_x_fctr) %>%
        dplyr::ungroup()
    }
    table_temp <-
      table_temp %>%
      dplyr::arrange(!!rlang::sym(fill_cat), rel_x_cumsum)

    if (!missing(label_rel_nudge) || !missing(label_abs_nudge) || !is.null(label_position)) {
      # print to show order of entries in the table which relevant for nudging labels
      message("order for label nudging or positioning:")
      print(table_temp[,c(1:5, which(names(table_temp) == "n_label_ypos"), which(names(table_temp) == "label_ypos_fctr"))], n = Inf)
    }

    if (y == "rel") {
      y_plot <- "label_ypos_fctr"
    } else if (y == "abs") {
      y_plot <- "n_label_ypos"
    }

    if (!is.null(label_position)) {
      if (!is.null(names(label_position))) {
        if (!anyDuplicated(names(label_position)) && !anyDuplicated(table_temp[[x_cat]])) {
          ## use names of label_position
          # order
          label_position <- label_position[as.character(table_temp[[x_cat]])]
          # handle missing values
          label_position <- purrr::map2_dbl(label_position, table_temp[[y_plot]], ~ if (is.na(.x)) .y else .x)
        } else {
          message("names of label_position cannot be used due to duplicates. will stick to the order provided.")
        }
      }
      if (length(label_position) != nrow(table_temp)) {
        message("length of label_position does not match number of labels to position.")
      } else {
        table_temp[[y_plot]] <- label_position
      }
    }

    if (plot_rel_labels) {
      plot <- plot_label_fun(nudge = label_rel_nudge,
                             size = label_size_rel,
                             label_var = ifelse(label_rel_pct, "rel_x_fctr_pct", "rel_x_fctr_round"),
                             name = "label_rel_nudge",
                             table_temp = table_temp,
                             x_cat = x_cat,
                             y_plot = y_plot,
                             label_color = label_color,
                             plot = plot)
    }
    if (plot_abs_labels) {
      plot <- plot_label_fun(nudge = label_abs_nudge,
                             size = label_size_abs,
                             label_var = "n",
                             name = "label_abs_nudge",
                             table_temp = table_temp,
                             x_cat = x_cat,
                             y_plot = y_plot,
                             label_color = label_color,
                             plot = plot)
    }
  }

  if (plot_total_rel_labels) {
    plot <-
      plot +
      ggplot2::geom_text(data = table0,
                         ggplot2::aes(color = I(label_color), label = if(label_rel_pct) {pct2} else {pct_round},
                                      x = !!rlang::sym(x_cat), y = total_labels_ypos),
                         nudge_x = label_total_rel_nudge[1], nudge_y = label_total_rel_nudge[2],
                         size = 4, inherit.aes = F)
  }
  if (plot_total_abs_labels) {
    plot <-
      plot +
      ggplot2::geom_text(data = table0,
                         ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(x_cat), y = total_labels_ypos),
                         nudge_x = label_total_abs_nudge[1], nudge_y = label_total_abs_nudge[2],
                         size = 4, inherit.aes = F)
  }

  if (flip) {
    plot <-
      plot +
      ggplot2::coord_flip()
    # if (plot_sigma) {
    #   plot <- plot + ggplot2::annotate("text", label = "\u03A3", x = Inf, y = total_labels_ypos*fctr, vjust = sigma_just, size = label_size)
    # }
  } #else {
    # if (plot_sigma) {
    #   plot <- plot + ggplot2::annotate("text", label = "\u03A3", x = -Inf, y = total_labels_ypos*fctr, hjust = sigma_just, size = label_size)
    # }
  #}

  return(list(data = table, plot = plot, data_total = table0))
}

plot_label_fun <- function(nudge, size, label_var, name, table_temp, x_cat, y_plot, label_color, plot) {

  if (!is.null(names(nudge))) {
    if (!anyDuplicated(names(nudge)) && !anyDuplicated(table_temp[[x_cat]])) {
      ## use names of nudge
      # order
      nudge <- nudge[as.character(table_temp[[x_cat]])]
      # handle missing values
      nudge <- purrr::map(nudge, ~ if (is.null(.x)) c(0,0) else .x)
    } else {
      message("names of ", name, " cannot be used due to duplicates. will stick to the order provided.")
      # brathering::recycle only when names of nudge not used
      nudge <- brathering::recycle(nudge, table_temp[,1,drop=T])
    }
  } else {
    nudge <- brathering::recycle(nudge, table_temp[,1,drop=T])
  }
  size <- brathering::recycle(size, nudge)
  for (j in seq_along(nudge)) {
    plot <-
      plot +
      ggplot2::geom_text(data = table_temp[j,],
                         ggplot2::aes(color = I(label_color), label = !!rlang::sym(label_var), x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
                         nudge_x = nudge[[j]][1], nudge_y = nudge[[j]][2],
                         size = size[j])
  }
  return(plot)
}


#   if (!is.null(names(label_rel_nudge))) {
#     if (!anyDuplicated(names(label_rel_nudge)) && !anyDuplicated(table_temp[[x_cat]])) {
#       ## use names of label_rel_nudge
#       # order
#       label_rel_nudge <- label_rel_nudge[as.character(table_temp[[x_cat]])]
#       # handle missing values
#       label_rel_nudge <- purrr::map(label_rel_nudge, ~ if (is.null(.x)) c(0,0) else .x)
#     } else {
#       message("names of label_rel_nudge cannot be used due to duplicates. will stick to the order provided.")
#       # brathering::recycle only when names of not used
#       label_rel_nudge <- brathering::recycle(label_rel_nudge, table_temp[,1,drop=T])
#     }
#   } else {
#     label_rel_nudge <- brathering::recycle(label_rel_nudge, table_temp[,1,drop=T])
#   }
#   label_size_rel <- brathering::recycle(label_size_rel, label_rel_nudge)
#   for (j in seq_along(label_rel_nudge)) {
#     plot <-
#       plot +
#       ggplot2::geom_text(data = table_temp[j,],
#                          ggplot2::aes(color = I(label_color), label = rel_x_fctr_pct, x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
#                          nudge_x = label_rel_nudge[[j]][1], nudge_y = label_rel_nudge[[j]][2],
#                          size = label_size_rel[j])
#   }
# }

#   if (plot_abs_labels) {
#     label_abs_nudge <- brathering::recycle(label_abs_nudge, table_temp[,1,drop=T])
#     label_size_abs <- brathering::recycle(label_size_abs, label_abs_nudge)
#     if (!is.null(names(label_abs_nudge))) {
#       if (!anyDuplicated(names(label_abs_nudge)) && !anyDuplicated(table_temp[[x_cat]])) {
#         label_abs_nudge <- label_abs_nudge[as.character(table_temp[[x_cat]])]
#       } else {
#         message("names of label_abs_nudge cannot be used due to duplicates. will stick to the order provided.")
#       }
#     }
#     for (j in seq_along(label_abs_nudge)) {
#       plot <-
#         plot +
#         ggplot2::geom_text(data = table_temp[j,],
#                            ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
#                            nudge_x = label_abs_nudge[[j]][1], nudge_y = label_abs_nudge[[j]][2],
#                            size = label_size_abs[j])
#     }
#   }
# }
