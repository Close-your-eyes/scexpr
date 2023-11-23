#' Plot composition of one category by another
#'
#' @param SO Seurat object or data frame
#' @param x_cat
#' @param fill_cat
#' @param border_color
#' @param plot_rel_labels
#' @param min_label_freq
#' @param label_rel_pct
#' @param label_rel_pct_decimals
#' @param col_pal
#' @param plot_abs_labels
#' @param summarize_all_x
#' @param label_rel_nudge
#' @param label_abs_nudge
#' @param flip
#' @param label_size
#' @param plot_total_labels
#' @param total_labels_ypos
#' @param plot_sigma
#' @param y
#' @param label_only_largest
#' @param sigma_just
#'
#' @return
#' @export
#'
#' @examples
composition_barplot <- function(SO,
                                x_cat,
                                fill_cat,
                                y = c("rel", "abs"),
                                col_pal = scexpr::col_pal("custom"),
                                border_color = "black",
                                plot_rel_labels = F,
                                label_only_largest = F,
                                plot_abs_labels = F,
                                plot_total_labels = F,
                                total_labels_ypos = 1.05,
                                plot_sigma = F,
                                sigma_just = ifelse(flip, 2, -1),
                                min_label_freq = 0.05,
                                label_rel_pct = T,
                                label_rel_pct_decimals = 2,
                                summarize_all_x = F,
                                label_rel_nudge = c(0,0),
                                label_abs_nudge = c(0,0),
                                label_size = 4, # vector of label_size only considered if label_rel_nudge and label_abs_nudge are lists (or one of them)
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

  y = match.arg(y, c("rel", "abs"))
  data_total <- NULL

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

  if (plot_total_labels) {
    table0 <-
      SO %>%
      dplyr::count(!!rlang::sym(x_cat)) %>%
      dplyr::mutate(pct = n/sum(n))

    if (summarize_all_x) {
      temp <- data.frame(x_cat = "all", n = nrow(SO), pct = 1)
      names(temp)[1] <- x_cat
      table0 <- dplyr::bind_rows(table0, temp)
      table0$label_ypos <- total_labels_ypos
      table0$label_color = "black" # could be made relative to background color
    }
    data_total <- table0
  }



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
    dplyr::arrange(desc(!!rlang::sym(fill_cat)), .by_group = T) %>%
    dplyr::mutate(rel_x_cumsum = cumsum(rel_x)) %>%
    dplyr::mutate(rel_x_cumsum_lag = dplyr::lag(rel_x_cumsum, default = 0)) %>%
    dplyr::mutate(label_ypos = rel_x_cumsum_lag + (rel_x_cumsum-rel_x_cumsum_lag)/2) %>%
    dplyr::mutate(n_cumsum = cumsum(n)) %>%
    dplyr::mutate(n_cumsum_lag = dplyr::lag(n_cumsum, default = 0)) %>%
    dplyr::mutate(n_label_ypos = n_cumsum_lag + (n_cumsum-n_cumsum_lag)/2) %>%
    dplyr::mutate(rel_x_fctr = rel_x*fctr) %>%
    dplyr::mutate(rel_x_fctr_pct = paste0(smart.round2(rel_x_fctr, label_rel_pct_decimals), " %")) %>%
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
      ggplot2::geom_col(color = border_color) +
      ggplot2::scale_fill_manual(values = col_pal) +
      ggplot2::labs(y = ifelse(label_rel_pct, "frequency [%]", "frequency")) +
      scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1)*fctr)
  } else if (y == "abs") {
    plot <-
      ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(x_cat), y = n, fill = !!rlang::sym(fill_cat))) +
      ggplot2::geom_col(color = border_color) +
      ggplot2::scale_fill_manual(values = col_pal)
  }

  if (plot_rel_labels || plot_abs_labels) {
    table_temp <- table %>% dplyr::filter(rel_x >= min_label_freq)

    if (y == "rel") {
      y_plot <- "label_ypos_fctr"
    } else if (y == "abs") {
      y_plot <- "n_label_ypos"
    }

    if (label_only_largest) {
      table_temp <-
        table_temp %>%
        dplyr::group_by(!!rlang::sym(x_cat)) %>%
        dplyr::slice_max(rel_x_fctr) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(clusters, rel_x_cumsum)
    } else {
      table_temp <-
        table %>%
        dplyr::arrange(clusters, rel_x_cumsum)
    }

    if (plot_rel_labels) {
      ## individual nudging:
      # add nudge_x and nudge_y to table
      # order the table first by x_cat and cumsum
      # this allows to provide a list of vectors for nudge_x and nudge_y
      if (methods::is(label_rel_nudge, "list")) {
        if (length(label_rel_nudge) != nrow(table_temp)) {
          message("Length of label_rel_nudge and number of labels are not equal.")
          if (label_size != 1 && length(label_size) != length(label_abs_nudge)) {
            message("Length of label_size and number of labels are not equal. Will use the first entry of label.size for all.")
          }
        } else {
          if (length(label_size) != length(label_rel_nudge)) {
            label_size <- rep(label_size[1], length(label_rel_nudge))
          }
          for (j in seq_along(label_rel_nudge)) {
            plot <-
              plot +
              ggplot2::geom_text(data = table_temp[j,],
                                 ggplot2::aes(color = I(label_color), label = rel_x_fctr_pct, x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
                                 nudge_x = label_rel_nudge[[j]][1], nudge_y = label_rel_nudge[[j]][2],
                                 size = label_size[j])
          }
        }
      } else {
        plot <-
          plot +
          ggplot2::geom_text(data = table_temp,
                             ggplot2::aes(color = I(label_color), label = rel_x_fctr_pct, x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
                             nudge_x = label_rel_nudge[1], nudge_y = label_rel_nudge[2],
                             size = label_size)
      }
    }

    if (plot_abs_labels) {
      if (methods::is(label_abs_nudge, "list")) {
        if (length(label_abs_nudge) != nrow(table_temp)) {
          message("Length of label_rel_nudge and number of labels are not equal.")
          if (label_size != 1 && length(label_size) != length(label_abs_nudge)) {
            message("Length of label_size and number of labels are not equal. Will use the first entry of label.size for all.")
          }
        } else {
          if (length(label_size) != length(label_abs_nudge)) {
            label_size <- rep(label_size[1], length(label_abs_nudge))
          }
          for (j in seq_along(label_abs_nudge)) {
            plot <-
              plot +
              ggplot2::geom_text(data = table_temp[j,],
                                 ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
                                 nudge_x = label_abs_nudge[[j]][1], nudge_y = label_abs_nudge[[j]][2],
                                 size = label_size[j])
          }
        }
      } else {
        plot <-
          plot +
          ggplot2::geom_text(data = table_temp,
                             ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(x_cat), y = !!rlang::sym(y_plot)),
                             nudge_x = label_abs_nudge[1], nudge_y = label_abs_nudge[2],
                             size = label_size)
      }
    }
  }

  if (plot_total_labels) {
    # TODO: adjust to plotting n (abs values, not rel)
    if (plot_rel_labels) {
      plot <-
        plot +
        ggplot2::geom_text(data = table0,
                           ggplot2::aes(color = I(label_color), label = paste0(smart.round2(pct*fctr, label_rel_pct_decimals), " %"), x = !!rlang::sym(x_cat), y = label_ypos*fctr),
                           nudge_x = label_rel_nudge[1], nudge_y = label_rel_nudge[2],
                           size = label_size, inherit.aes = F)
    }
    if (plot_abs_labels) {
      plot <-
        plot +
        ggplot2::geom_text(data = table0,
                           ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(x_cat), y = label_ypos*fctr),
                           nudge_x = label_abs_nudge[1], nudge_y = label_abs_nudge[2],
                           size = label_size, inherit.aes = F)
    }
  }

  if (flip) {
    plot <-
      plot +
      ggplot2::coord_flip()
    if (plot_sigma) {
      plot <- plot + annotate("text", label = "\u03A3", x = Inf, y = total_labels_ypos*fctr, vjust = sigma_just, size = label_size)
    }
  } else {
    if (plot_sigma) {
      plot <- plot + annotate("text", label = "\u03A3", x = -Inf, y = total_labels_ypos*fctr, hjust = sigma_just, size = label_size)
    }
  }

  return(list(data = table, plot = plot, data_total = data_total))
}


#https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart.round1 <- function(x) {
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}

smart.round2 <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}


