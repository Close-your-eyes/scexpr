feature_plot_gene <- function(plot,
                              freqs,
                              pt_size = 0.3,
                              pt_size_fct = 1,
                              col_expr = "tomato2",
                              col_non_expr = "grey85",
                              col_binary = F,
                              col_split = "grey30",
                              freq_plot = T,
                              plot_all_across_split = F,
                              freq_pos = c(-Inf, Inf),
                              freq_size = 4) {

  shapeby <- tryCatch(rlang::sym(attr(plot[["data"]], "shape_feature")), error = function(e) NULL)

  if (plot_all_across_split && "split_feature" %in% names(attributes(plot[["data"]]))) {
    plot <- plot + ggplot2::geom_point(data = ~dplyr::select(., -split_feature |> dplyr::filter(cells == 1),
                                       mapping = ggplot2::aes(shape = !!shapeby),
                                       size = pt_size,
                                       color = col_split))
  }

  if (col_binary) {
    if (any(dplyr::filter(plot[["data"]], cells == 1)[["feature"]] < 0)) {
      message("col_binary (+/-) may not be meaningful as there are negative expression values for ", attr(plot[["data"]], "feature"), ".")
    }
    plot[["data"]][["valuebin"]] <- ifelse(plot[["data"]][["feature"]] > 0, "+", "-")
    plot <- plot + ggplot2::geom_point(
      data = ~dplyr::filter(., cells == 1),
      mapping = ggplot2::aes(shape = !!shapeby, color = valuebin),
      size = pt_size*pt_size_fct) +
      ggplot2::scale_color_manual(values = c(col_non_expr, col_expr))
  } else {
    # non-expressers
    plot <- plot + ggplot2::geom_point(
      data = ~dplyr::filter(., cells == 1 & feature == 0),
      mapping = ggplot2::aes(shape = !!shapeby),
      size = pt_size,
      color = col_non_expr)
    # expressers
    plot <- plot + ggplot2::geom_point(
      data = ~dplyr::filter(., cells == 1 & feature > 0),
      mapping = ggplot2::aes(color = feature, shape = !!shapeby),
      size = pt_size*pt_size_fct)
  }


  if (freq_plot) {
    if (nlevels(plot[["data"]][["SO.split"]]) > 1 && !"split_feature" %in% names(attributes(plot[["data"]]))) {
      freqs2 <- freqs[["freq.expr.by.SO"]]
    } else if (nlevels(plot[["data"]][["SO.split"]]) > 1 && "split_feature" %in% names(attributes(plot[["data"]]))) {
      freqs2 <- freqs[["freq.expr.by.split.SO"]]
    } else if ("split_feature" %in% names(attributes(plot[["data"]])) && nlevels(plot[["data"]][["SO.split"]]) == 1) {
      freqs2 <- freqs[["freq.expr.by.split"]]
    } else {
      freqs2 <- freqs[["freq.expr"]]
    }

    plot <- plot + ggrepel::geom_text_repel(
      data = freqs2,
      size = freq_size,
      ggplot2::aes(
        label = freq2,
        x = freq_pos[1],
        y = freq_pos[2]
      )
    )
  }

  return(plot)
}
