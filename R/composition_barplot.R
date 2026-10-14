#' Plot composition of one category by another
#'
#' Count observations for each combination of two metadata columns and display
#' their composition as stacked bars. Each row of the input represents one
#' observation. Bars can show within-category frequencies or absolute counts.
#'
#' @param obj A Seurat object or a data frame. For a Seurat object, the
#'   `meta.data` slot supplies the observations.
#' @param xvar Character scalar naming the column used for the x-axis categories.
#'   Numeric columns are not supported; convert them to character or factor.
#' @param fillvar Character scalar naming the column used for bar segments.
#' @param yvar Character scalar: `"rel"` for frequencies within each x-axis category
#'   or `"abs"` for observation counts.
#' @param col_pal Palette name or color vector passed to
#'   [colrr::make_col_pal()]. The default,
#'   `"..auto.."`, uses `obj@misc$metacolors[[fillvar]]` when available and
#'   otherwise uses the `"custom"` palette.
#' @param geom_col_args Named list of arguments passed to [ggplot2::geom_col()].
#'   Label positions assume the default stacked bars.
#' @param plot_rel_labels Logical. Display within-category frequencies inside
#'   bar segments.
#' @param label_only_largest Logical. Restrict segment labels to the largest
#'   eligible segment in each x-axis category. Ties retain all largest segments.
#' @param plot_abs_labels Logical. Display observation counts inside bar segments.
#' @param plot_total_rel_labels Logical. Display each x-axis category's share
#'   of all observations above its bar.
#' @param plot_total_abs_labels Logical. Display the total observation count
#'   for each x-axis category above its bar.
#' @param min_label_freq Numeric scalar. Minimum within-category fraction for a
#'   segment to receive a label, regardless of `yvar`. Defaults to `0.05`.
#' @param label_rel_pct Logical. Display relative labels as percentages when
#'   `TRUE` and as fractions otherwise. Also controls the scale of the relative
#'   y-axis when `yvar = "rel"`.
#' @param label_rel_decimals Number of decimal places for relative labels.
#'   When omitted, uses zero for percentages and two for fractions.
#' @param label_rel_rm_zero Logical. Remove the leading zero from relative
#'   segment labels displayed as fractions.
#' @param summarize_all_x Logical. Add a bar named `"___all___"` pooling all
#'   observations. Its frequencies reflect pooled counts, rather than an
#'   unweighted average of the individual bars.
#' @param label_rel_nudge List of numeric pairs giving x and y offsets for
#'   relative segment labels. Recycled to the number of eligible labels.
#'   Names can match x-axis categories when both the names and the eligible
#'   x-axis categories are unique; missing named entries use zero offsets.
#' @param label_abs_nudge List of x and y offset pairs for absolute segment
#'   labels, with the same recycling and naming rules as `label_rel_nudge`.
#' @param label_size_rel Numeric vector of relative segment label sizes,
#'   recycled to the number of labels and passed to [ggplot2::geom_text()].
#' @param label_size_abs Numeric vector of absolute segment label sizes,
#'   recycled in the same way as `label_size_rel`.
#' @param label_total_rel_nudge Numeric pair giving x and y offsets for all
#'   relative total labels.
#' @param label_total_abs_nudge Numeric pair giving x and y offsets for all
#'   absolute total labels.
#' @param label_position Optional numeric vector replacing the vertical segment
#'   label positions, in the plotted y-axis units. Its length must match the
#'   number of eligible labels. Unique names can match x-axis categories when
#'   each category has at most one eligible label; missing named positions use
#'   the calculated centers. Applies to both relative and absolute labels.
#' @param flip Logical. Flip the axes using [ggplot2::coord_flip()].
#' @param bw_cutoff Numeric cutoff passed to [brathering::bw_txt()] to select
#'   black or white segment label text from the segment background color.
#' @param theme A ggplot2 theme object. Defaults to
#'   `colrr::theme_material(white = TRUE)`.
#' @param theme_args Named list of additional modifications passed to
#'   [ggplot2::theme()] after applying `theme`. By default, x-axis text is
#'   rotated by 45 degrees and right-aligned.
#'
#' @details
#' Segment labels are centered using cumulative segment sizes before filtering
#' by `min_label_freq`. If both label types are enabled, they share the same
#' centers unless different nudges are supplied. Explicit nudges or positions
#' print the eligible label rows in the order used for positioning.
#'
#' Total labels use a common height above the bars. The pooled bar, when
#' requested, has a relative total of 100 percent, or one when fractions are used.
#'
#' Missing values in `fillvar` currently prevent segment labels from rendering.
#'
#' @return A named list containing:
#' \describe{
#'   \item{plot}{A ggplot object.}
#'   \item{data}{An ungrouped tibble of segment counts (`n`), category totals
#'     (`x_total` and `fill_total`), within-category frequencies (`rel_x` and
#'     `rel_fill`), and columns used for colors and label positioning.}
#'   \item{data_total}{A data frame of x-axis category counts, relative totals,
#'     and total-label positions when either total-label option is enabled;
#'     otherwise `NULL`.}
#' }
#' @export
#' @importFrom rlang :=
#'
#' @examples
#' metadata <- data.frame(
#'   sample = c("sample_1", "sample_1", "sample_2", "sample_2"),
#'   cell_type = c("T cell", "B cell", "T cell", "T cell")
#' )
#' result <- composition_barplot(
#'   metadata, xvar = "sample", fillvar = "cell_type",
#'   col_pal = "custom", plot_rel_labels = TRUE
#' )
#' result$plot
#' result$data
composition_barplot <- function(obj,
                                xvar,
                                fillvar,
                                yvar = c("rel", "abs"),
                                col_pal = "..auto..",
                                geom_col_args = list(color = "black"),
                                plot_rel_labels = F,
                                label_only_largest = F,
                                plot_abs_labels = F,
                                plot_total_rel_labels = F,
                                plot_total_abs_labels = F,
                                min_label_freq = 0.05,
                                label_rel_pct = T,
                                label_rel_decimals = 2,
                                label_rel_rm_zero = T,
                                summarize_all_x = F,
                                label_rel_nudge = list(c(0,0)),
                                label_abs_nudge = list(c(0,0)),
                                label_size_rel = 4,
                                label_size_abs = 4,
                                label_total_rel_nudge = c(0,0),
                                label_total_abs_nudge = c(0,0),
                                label_position = NULL,
                                flip = F,
                                bw_cutoff = 50,
                                theme = colrr::theme_material(white = T),
                                theme_args = list(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))) {

  scexpr:::.ensure_package("colrr")
  scexpr:::.ensure_package("brathering")

  if (methods::is(obj, "Seurat") && col_pal[1] == "..auto.." && "metacolors" %in% names(obj@misc) && fillvar %in% names(obj@misc$metacolors)) {
    col_pal <- obj@misc$metacolors[[fillvar]]
  } else if (col_pal[1] == "..auto..") {
    col_pal <- "custom"
  }

  if (methods::is(obj, "Seurat")) {
    obj <- obj@meta.data
  }

  if (!xvar %in% names(obj)) {
    stop("xvar not found in obj.")
  }
  if (!fillvar %in% names(obj)) {
    stop("fillvar not found in obj.")
  }
  if (is.numeric(obj[[xvar]])) {
    stop("xvar must not be numeric.")
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

  yvar <- rlang::arg_match(yvar)

  table <- obj |>
    dplyr::count(!!rlang::sym(xvar), !!rlang::sym(fillvar)) |>
    dplyr::left_join(dplyr::count(obj, !!rlang::sym(xvar), name = "x_total"), by = xvar) |>
    dplyr::left_join(dplyr::count(obj, !!rlang::sym(fillvar), name = "fill_total"), by = fillvar) |>
    dplyr::mutate(rel_x = n/x_total) |>
    dplyr::mutate(rel_fill = n/fill_total) # rel_fill is not used below

  # cl <- compare_labels(obj[[xvar]], obj[[fillvar]])
  #   browser()

  if (summarize_all_x) {
    table <- dplyr::bind_rows(table,
                              obj |>
                                dplyr::count(!!rlang::sym(fillvar)) |>
                                dplyr::mutate(x_total = nrow(obj)) |>
                                dplyr::mutate(fill_total = nrow(obj)) |>
                                dplyr::mutate(rel_x = n/x_total) |>
                                dplyr::mutate(rel_fill = n/fill_total) |>
                                dplyr::mutate(!!xvar := "___all___")) # show that it is 'special' and make it appear first
  }

  table0 <- NULL
  if (plot_total_rel_labels || plot_total_abs_labels) {
    table0 <- obj |>
      dplyr::count(!!rlang::sym(xvar)) |>
      dplyr::mutate(pct = n/sum(n))

    if (summarize_all_x) {
      temp <- data.frame(xvar = "___all___", n = nrow(obj), pct = 1)
      names(temp)[1] <- xvar
      table0 <- dplyr::bind_rows(table0, temp)
    }
    if (yvar == "abs") {
      table0$total_labels_ypos <- max(table0$n)*1.05
    }
    if (yvar == "rel") {
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

  # always works?
  # non existing factor levels ... problem
  col_pal <- colrr::make_col_pal(col_vec = col_pal,
                                 fct_lvls = sort(unique(table[,fillvar,drop=T])))

  # check this somewhen; in some cases label_color is not defined and then an error occurs below
  if (!is.null(names(col_pal)) && any(!unique(table[,fillvar,drop=T]) %in% names(col_pal))) {
    col_pal <- unname(col_pal)
  } else {
    # new
    if (is.null(names(col_pal))) {
      names(col_pal) <- sort(unique(table[,fillvar,drop=T]))
    }
    table$bar_segment_cols <- col_pal[as.character(table[,fillvar,drop=T])]
    table$label_color <- brathering::bw_txt(table$bar_segment_cols, cutoff = bw_cutoff)
  }

  ## use cumsum approach to define label position in middle of bars
  ## this allow to remove label by min_label_freq and still have the correct coordinate
  ## coordinates became wrong with position = ggplot2::position_stack(vjust = 0.5) when some values were removed


  table <- table |>
    # make sure NA come first in order, as they will be plotted first - needed for cumsum and correct label position
    dplyr::mutate(sort_helper = ifelse(is.na(!!rlang::sym(fillvar)), 0, 1)) |>
    dplyr::group_by(!!rlang::sym(xvar)) |>
    dplyr::arrange(sort_helper, dplyr::desc(!!rlang::sym(fillvar)), .by_group = T) |>
    dplyr::mutate(rel_x_cumsum = cumsum(rel_x)) |>
    dplyr::mutate(rel_x_cumsum_lag = dplyr::lag(rel_x_cumsum, default = 0)) |>
    dplyr::mutate(label_ypos = rel_x_cumsum_lag + (rel_x_cumsum-rel_x_cumsum_lag)/2) |>
    dplyr::mutate(n_cumsum = cumsum(n)) |>
    dplyr::mutate(n_cumsum_lag = dplyr::lag(n_cumsum, default = 0)) |>
    dplyr::mutate(n_label_ypos = n_cumsum_lag + (n_cumsum-n_cumsum_lag)/2) |>
    dplyr::mutate(rel_x_fctr = rel_x*fctr) |>
    dplyr::mutate(rel_x_fctr_round = brathering::round2(rel_x_fctr, label_rel_decimals)) |>
    dplyr::mutate(rel_x_fctr_round = ifelse(
      label_rel_rm_zero,
      gsub("^0", "", as.character(rel_x_fctr_round)),
      rel_x_fctr_round
    )) |>
    dplyr::mutate(rel_x_fctr_pct = paste0(brathering::round2(rel_x_fctr, label_rel_decimals), " %")) |>
    dplyr::mutate(label_ypos_fctr = label_ypos*fctr) |>
    tibble::as_tibble() |>
    dplyr::ungroup()

  if (!is.null(levels(obj[,xvar]))) {
    if (summarize_all_x) {
      table[,xvar] <- factor(table[,xvar,drop=T], levels = c("___all___", levels(obj[,xvar])))
    } else {
      table[,xvar] <- factor(table[,xvar,drop=T], levels = levels(obj[,xvar]))
    }
  }

  if (yvar == "rel") {
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(xvar), y = rel_x*fctr, fill = !!rlang::sym(fillvar))) +
      ggplot2::labs(y = ifelse(label_rel_pct, "frequency [%]", "frequency"))
    breaks <- c(0,0.25,0.5,0.75,1)*fctr
  } else if (yvar == "abs") {
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = !!rlang::sym(xvar), y = n, fill = !!rlang::sym(fillvar)))
    breaks <- ggplot2::waiver()
  }
  plot <- plot +
    do.call(ggplot2::geom_col, args = geom_col_args) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::scale_y_continuous(breaks = breaks, expand = ggplot2::expansion(mult = 0.01))

  if (plot_rel_labels || plot_abs_labels) {

    table_temp <- dplyr::filter(table, rel_x >= min_label_freq)
    if (nrow(table_temp) < nrow(table)) {
      message(nrow(table) - nrow(table_temp), " text labels removed due to min_label_freq.")
    }
    if (label_only_largest) {
      table_temp <- dplyr::slice_max(table_temp, order_by = rel_x_fctr, by = !!rlang::sym(xvar))
    }
    table_temp <- dplyr::arrange(table_temp, !!rlang::sym(fillvar), rel_x_cumsum)

    if (!missing(label_rel_nudge) || !missing(label_abs_nudge) || !is.null(label_position)) {
      # print to show order of entries in the table which relevant for nudging labels
      message("order for label nudging or positioning:")
      print(table_temp[,c(1:5, which(names(table_temp) == "n_label_ypos"), which(names(table_temp) == "label_ypos_fctr"))], n = Inf)
    }

    if (yvar == "rel") {
      y_plot <- "label_ypos_fctr"
    } else if (yvar == "abs") {
      y_plot <- "n_label_ypos"
    }

    if (!is.null(label_position)) {
      if (!is.null(names(label_position))) {
        if (!anyDuplicated(names(label_position)) && !anyDuplicated(table_temp[[xvar]])) {
          ## use names of label_position
          # order
          label_position <- label_position[as.character(table_temp[[xvar]])]
          # handle missing values
          label_position <- purrr::map2_dbl(label_position, table_temp[[y_plot]], ~ if (is.na(.x)) {.y} else {.x})
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
                             xvar = xvar,
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
                             xvar = xvar,
                             y_plot = y_plot,
                             label_color = label_color,
                             plot = plot)
    }
  }

  if (plot_total_rel_labels) {
    plot <- plot +
      ggplot2::geom_text(data = table0,
                         ggplot2::aes(color = I(label_color), label = if(label_rel_pct) {pct2} else {pct_round},
                                      x = !!rlang::sym(xvar), y = total_labels_ypos),
                         nudge_x = label_total_rel_nudge[1], nudge_y = label_total_rel_nudge[2],
                         size = 4, inherit.aes = F)
  }
  if (plot_total_abs_labels) {
    plot <- plot +
      ggplot2::geom_text(data = table0,
                         ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(xvar), y = total_labels_ypos),
                         nudge_x = label_total_abs_nudge[1], nudge_y = label_total_abs_nudge[2],
                         size = 4, inherit.aes = F)
  }

  if (flip) {
    plot <- plot + ggplot2::coord_flip()
  }

  plot <- plot +
    theme +
    do.call(ggplot2::theme, args = theme_args)

  return(list(plot = plot, data = table, data_total = table0))
}

plot_label_fun <- function(nudge,
                           size,
                           label_var,
                           name,
                           table_temp,
                           xvar,
                           y_plot,
                           label_color,
                           plot) {

  if (!is.null(names(nudge))) {
    if (!anyDuplicated(names(nudge)) && !anyDuplicated(table_temp[[xvar]])) {
      ## use names of nudge
      # order
      nudge <- nudge[as.character(table_temp[[xvar]])]
      # handle missing values
      nudge <- purrr::map(nudge, ~ if (is.null(.x)) {c(0,0)} else {.x})
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
    plot <- plot +
      ggplot2::geom_text(data = table_temp[j,],
                         ggplot2::aes(
                           color = I(label_color),
                           label = !!rlang::sym(label_var),
                           x = !!rlang::sym(xvar),
                           y = !!rlang::sym(y_plot)
                         ),
                         nudge_x = nudge[[j]][1], nudge_y = nudge[[j]][2],
                         size = size[j])
  }
  return(plot)
}


#   if (!is.null(names(label_rel_nudge))) {
#     if (!anyDuplicated(names(label_rel_nudge)) && !anyDuplicated(table_temp[[xvar]])) {
#       ## use names of label_rel_nudge
#       # order
#       label_rel_nudge <- label_rel_nudge[as.character(table_temp[[xvar]])]
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
#                          ggplot2::aes(color = I(label_color), label = rel_x_fctr_pct, x = !!rlang::sym(xvar), y = !!rlang::sym(y_plot)),
#                          nudge_x = label_rel_nudge[[j]][1], nudge_y = label_rel_nudge[[j]][2],
#                          size = label_size_rel[j])
#   }
# }

#   if (plot_abs_labels) {
#     label_abs_nudge <- brathering::recycle(label_abs_nudge, table_temp[,1,drop=T])
#     label_size_abs <- brathering::recycle(label_size_abs, label_abs_nudge)
#     if (!is.null(names(label_abs_nudge))) {
#       if (!anyDuplicated(names(label_abs_nudge)) && !anyDuplicated(table_temp[[xvar]])) {
#         label_abs_nudge <- label_abs_nudge[as.character(table_temp[[xvar]])]
#       } else {
#         message("names of label_abs_nudge cannot be used due to duplicates. will stick to the order provided.")
#       }
#     }
#     for (j in seq_along(label_abs_nudge)) {
#       plot <-
#         plot +
#         ggplot2::geom_text(data = table_temp[j,],
#                            ggplot2::aes(color = I(label_color), label = n, x = !!rlang::sym(xvar), y = !!rlang::sym(y_plot)),
#                            nudge_x = label_abs_nudge[[j]][1], nudge_y = label_abs_nudge[[j]][2],
#                            size = label_size_abs[j])
#     }
#   }
# }
