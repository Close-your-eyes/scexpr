feature_plot_meta <- function(plot,
                              order_discr_explicit = NULL,
                              pt_size = 0.3,
                              col_split = "grey30",
                              plot_all_across_split = F,
                              shape_by = NULL) {

  shapeby <- tryCatch(rlang::sym(attr(plot[["data"]], "shape_feature")), error = function(e) NULL)

  if (plot_all_across_split && "split_feature" %in% names(attributes(plot[["data"]]))) {
    plot <- plot + ggplot2::geom_point(data = ~dplyr::select(., -split_feature |> dplyr::filter(cells == 1),
                                                             mapping = ggplot2::aes(shape = !!shapeby),
                                                             size = pt_size,
                                                             color = col_split))
  }

  if (is.null(order_discr_explicit)) {
    # order_discr_explicit T or F: ordering has been done in .get.data
    # use rownames(data) %in% ... to preserve random order

    plot <- plot +
      ggplot2::geom_point(data = ~dplyr::filter(., cells == 1),
                          mapping = ggplot2::aes(color = feature, shape = !!shapeby),
                          size = pt_size)
  } else {
    # order is explicit
    if (length(unique(order_discr_explicit)) != length(order_discr_explicit)) {
      order_discr_explicit <- unique(order_discr_explicit)
      message("order_discr_explicit is made unique: ", paste(order_discr_explicit, collapse = ", "))
    }

    if (any(!gsub("\\$$", "", (gsub("^\\^", "", order_discr_explicit))) %in% levels(as.factor(plot[["data"]][["value"]])))) {
      order_discr_explicit <- order_discr_explicit[which(gsub("\\$$", "", (gsub("^\\^", "", order_discr_explicit))) %in% levels(as.factor(plot[["data"]][["value"]])))]
      message("order_discr_explicit reduced to existing levels: ", paste(gsub("\\$$", "", (gsub("^\\^", "", order_discr_explicit))), collapse = ", "))
    }

    for (i in order_discr_explicit[which(grepl("^\\^", order_discr_explicit))]) {
      i <- gsub("^\\^", "", i)
      plot <- plot + ggplot2::geom_point(
        data = ~dplyr::filter(., cells == 1 & value == i),
        mapping = ggplot2::aes(color = feature, shape = !!shapeby),
        size = pt_size
      )
    }

    for (i in order_discr_explicit[intersect(which(!grepl("^\\^", order_discr_explicit)), which(!grepl("\\$$", order_discr_explicit)))]) {
      plot <- plot + ggplot2::geom_point(
        data = ~dplyr::filter(., cells == 1 & value == i),
        mapping = ggplot2::aes(color = feature, shape = !!shapeby),
        size = pt_size
      )
    }
    plot <- plot + ggplot2::geom_point(data = ~dplyr::filter(., cells == 1 & value %in% gsub("\\$$", "", (gsub("^\\^", "", order_discr_explicit)))),
                                       mapping = ggplot2::aes(color = feature, shape = !!shapeby),
                                       size = pt_size)

    for (i in order_discr_explicit[which(grepl("\\$$", order_discr_explicit))]) {
      i <- gsub("\\$$", "", i)
      plot <- plot + ggplot2::geom_point(
        data = ~dplyr::filter(., cells == 1 & value == i),
        ggplot2::aes(color = feature, shape = !!shapeby),
        size = pt_size
      )
    }
  }
  return(plot)
}
