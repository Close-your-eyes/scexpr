check.SO <- function(SO,
                     assay,
                     meta.col = NULL,
                     cells = NULL,
                     length = NULL) {
  if (!is.list(SO)) {
    SO <- list(SO)
  }

  if (!is.null(length) && length(SO) != length) {
    stop("Please provide exactly ", length, " SO.")
  }
  if (!any(unlist(lapply(SO, class)) == "Seurat")) {
    stop("All SO have to be Seurat objects (class == Seurat).")
  }
  if (!is.null(meta.col) && length(check.features(SO, meta.col, rownames = F)) == 0) {
    stop("meta.col not found in all objects.")
  }

  assay <- match.arg(assay, Reduce(intersect, lapply(SO, function(x) names(x@assays))))

  SO <- lapply(SO, function(x) {
    Seurat::DefaultAssay(x) <- assay
    return(x)
  })

  if (is.null(names(SO)) && length(SO) > 1) {
    message("List of SO has no names. Naming them numerically in order as provided.")
    names(SO) <- as.character(seq_along(SO))
  } else if (is.null(names(SO))) {
    names(SO) <- as.character(seq_along(SO))
  }

  ## check if data has been scaled (compare to counts)
  equal_slots <- unlist(lapply(SO, function(x) {
    if (utils::compareVersion(as.character(x@version), "4.9.9") == 1) {
      GetAssayData_args1 <- list(object = x, layer = "data", assay = assay)
      GetAssayData_args2 <- list(object = x, layer = "counts", assay = assay)
    } else {
      GetAssayData_args1 <- list(object = x, slot = "data", assay = assay)
      GetAssayData_args2 <- list(object = x, slot = "counts", assay = assay)
    }
    identical(Gmisc::fastDoCall(what = Seurat::GetAssayData, args = GetAssayData_args1),
              Gmisc::fastDoCall(what = Seurat::GetAssayData, args = GetAssayData_args2))

  }))
  if (any(equal_slots) && length(SO) > 1) {
    message("Data slot in at least one SO is equal to the counts slot. You may want to normalize.")
  } else if (any(equal_slots)) {
    message("Data slot is equal to the counts slot. You may want to normalize.")
  }

  if (!is.null(cells)) {
    allcells <- unique(unlist(purrr::map(SO, Seurat::Cells)))
    notfound <- cells[which(!cells %in% allcells)]
    if (length(notfound) > 0) {
      message(length(notfound), " of ", length(cells), " cells not found.")
    }
    cells <- cells[which(cells %in% allcells)]
    if (length(cells) == 0) {
      stop("None of cells found.")
    }
    assign("cells", cells, env = parent.frame())
  }

  if (!is.null(length) && length == 1) {
    return(SO[[1]])
  } else {
    return(SO)
  }
}

check.reduction <- function(SO,
                            reduction,
                            dims = c(1,2)) {

  if (!is.list(SO)) {
    SO <- list(SO)
  }

  red_list <- lapply(SO, function(x) names(x@reductions))
  if (any(lengths(red_list) == 0)) {
    stop("At least one SO has no reduction.")
  }

  common_red <- Reduce(intersect, red_list)
  if (length(common_red) == 0) {
    stop("No common reduction found in SOs.")
  }

  red <- grep(reduction, common_red, ignore.case = T, value = T)
  if (length(red) == 0) {
    message("reduction not found in SOs., Changing to an arbitrary common one in SOs.")
    red <- sample(common_red, 1)
  } else if (length(red) > 1) {
    red <- common_red[which.min(utils::adist(reduction, common_red, ignore.case = T))]
  }

  key <- SO[[1]]@reductions[[red]]@key
  red <- stats::setNames(gsub("_$", "", key), nm = red)
  return(red)
}

check.and.get.cells <- function(SO,
                                assay = c("RNA", "SCT"),
                                cells = NULL,
                                feature_cut = NULL,
                                feature_cut_expr = 0,
                                feature_ex  = NULL,
                                downsample = 1,
                                included_only = F) {

  if (!is.list(SO)) {
    SO <- list(SO)
  }

  feature_cut <- check.features(SO = SO, features = feature_cut, meta.data = F)
  feature_ex <- check.features(SO = SO, features = feature_ex, meta.data = F)

  if (downsample == 0) {
    stop("downsample cannot be 0.")
  }
  if (!is.null(feature_cut) && length(check.features(SO, feature_cut)) == 0) {
    stop("feature_cut not found in every SO.")
  }
  if (!is.null(feature_ex) && length(check.features(SO, feature_ex)) == 0) {
    stop("feature_ex not found in every SO.")
  }
  if (length(feature_cut) != length(feature_cut_expr) && length(feature_cut_expr) != 1) {
    stop("feature_cut and feature_cut_expr need to have the same length.")
  }
  if (length(feature_cut) != length(feature_cut_expr) && length(feature_cut_expr) == 1) {
    feature_cut_expr <- rep(feature_cut_expr, length(feature_cut))
  }

  assay <- match.arg(assay, c("RNA", "SCT"))

  if (any(duplicated(cells))) {
    message("Duplicates found in cells. Made unique now.")
    cells <- unique(cells)
  }

  # check if cell names are unique across SOs
  all.cells <- purrr::map(SO, Seurat::Cells)
  if (any(duplicated(unlist(all.cells)))) {
    message("Duplicated cell names found across SO.")
    # Get all pairwise combinations of indices
    combs <- combn(seq_along(all.cells), 2, simplify = F)
    # Check for intersection
    intersections <- lapply(combs, function(pair) {
      common <- intersect(all.cells[[pair[1]]], all.cells[[pair[2]]])
      list(pair = pair, intersect = common, has_intersect = length(common) > 0)
    })
    # Filter only those with intersection
    intersections_with_common <- Filter(function(x) x$has_intersect, intersections)
    message("Pairs with intersecting cell names:")
    print(sapply(intersections_with_common, "[", "pair"))
  }

  all.cells <- unique(unlist(all.cells))
  if (is.null(cells)) {
    cells <- all.cells
  } else {
    cells_l <- length(cells)
    cells <- intersect(cells, all.cells)
    if (length(cells) == 0) {
      stop("None of cells found in SO.")
    }
    if (length(cells) < cells_l) {
      message("Not all cells found in SO(s). Reduced to those which exist.")
    }
  }

  # create a vector of cells which identifies how to plot them; 0 indicates exclusion
  all.cells <- stats::setNames(rep(1, length(all.cells)), all.cells)
  ## cells excluded by arbitrary selection
  all.cells[!names(all.cells) %in% cells] <- 0

  ## cells excluded due to cutoff feature
  if (!is.null(feature_cut)) {
    compare_fun <- function(x,y) {x > y}
    exclude.cells <- names(which(Matrix::colSums(sweep(do.call(cbind, lapply(SO, function(x) {
      if (utils::compareVersion(as.character(x@version), "4.9.9") == 1) {
        GetAssayData_args <- list(object = x, layer = "data", assay = assay)
      } else {
        GetAssayData_args <- list(object = x, slot = "data", assay = assay)
      }
      Gmisc::fastDoCall(what = Seurat::GetAssayData, args = GetAssayData_args)[feature_cut,,drop=F]
    })), 1, feature_cut_expr, compare_fun)) < length(feature_cut)))
    all.cells[which(names(all.cells) %in% exclude.cells)] <- 0
  }
  if (length(all.cells) == 0) {
    stop("No cells left after filtering for feature_cut.")
  }
  ## cells excluded due to expression of an unwanted gene
  if (!is.null(feature_ex)) {
    exclude.cells <- unlist(lapply(SO, function(x) names(which(Matrix::colSums(do.call(cbind, lapply(SO, function(x) {
      if (utils::compareVersion(as.character(x@version), "4.9.9") == 1) {
        GetAssayData_args <- list(object = x,
                                  layer = "data",
                                  assay = assay)
      } else {
        GetAssayData_args <- list(object = x,
                                  slot = "data",
                                  assay = assay)
      }
      Gmisc::fastDoCall(
        what = Seurat::GetAssayData,
        args = GetAssayData_args
      )[feature_ex,,drop=F]

    }))) > 0))))
    all.cells[which(names(all.cells) %in% exclude.cells)] <- 0
  }

  if (length(all.cells) == 0) {
    stop("No cells left after filtering for feature_ex")
  }

  # downsample - completely remove cells to speed up plotting which may alter the statistics though
  if (downsample < 1) {
    downsample <- downsample*length(all.cells)
  } else if (downsample > 1) {
    downsample <- min(downsample, length(all.cells))
  }
  if (downsample == 0) {
    stop("downsample has become 0.")
  }
  if (downsample != 1) {
    all.cells <- all.cells[sample(seq_along(all.cells), size = downsample)]
  }

  if (included_only) {
    return(names(all.cells[which(all.cells == 1)]))
  } else {
    return(all.cells)
  }

}


check.aliases <- function(feature, feature.aliases, data) {
  if (!is.null(feature.aliases) && feature %in% names(feature.aliases)) {
    new_name <- as.character(feature.aliases[which(names(feature.aliases) == feature)])
    message(feature, " changed to ", new_name)
    attr(data, "feature") <- new_name
    #names(data)[which(names(data) == "feat")] <- new_name
  }
  return(data)
}

get.freqs2 <- function(data) {
  library(zeallot)
  freq.expr.by.split.SO <- dplyr::group_by(data, split_feature, SO.split)
  freq.expr.by.split <- dplyr::group_by(data, split_feature)
  freq.expr.by.SO <- dplyr::group_by(data, SO.split)
  freq.expr <- data

  d1 <- attr(data, "dim1")
  d2 <- attr(data, "dim2")

  c(freq.expr.by.split.SO,
    freq.expr.by.split,
    freq.expr.by.SO,
    freq.expr) %<-% purrr::map(list(freq.expr.by.split.SO,
                                    freq.expr.by.split,
                                    freq.expr.by.SO,
                                    freq.expr),
                               ~.x |>
                                 dplyr::filter(cells == 1) |>
                                 dplyr::summarise(freq = sum(feature>0)/dplyr::n()*100,
                                                  !!paste0(d1, "_min") := min(!!rlang::sym(d1)),
                                                  !!paste0(d1, "_max") := max(!!rlang::sym(d1)),
                                                  !!paste0(d2, "_min") := min(!!rlang::sym(d2)),
                                                  !!paste0(d2, "_max") := max(!!rlang::sym(d2)),
                                                  .groups = "drop") |>
                                 dplyr::mutate(freq2 = dplyr::case_when(freq>0 & freq<1 ~"<1 %",
                                                                        freq>1 & freq<99 ~paste0(round(freq, 0), " %"),
                                                                        freq>99 & freq<100 ~">99 %",
                                                                        .default = paste0(freq, " %"))))

  return(list(freq.expr.by.split.SO = freq.expr.by.split.SO,
              freq.expr.by.split = freq.expr.by.split,
              freq.expr.by.SO = freq.expr.by.SO,
              freq.expr = freq.expr))


}

get.freqs <- function(data,
                      cells,
                      reduction,
                      dims,
                      split_by,
                      SO.split) {






  freqs <- do.call(rbind, lapply(split(unique(data[,c("split_by", "SO.split")]), seq(nrow(unique(data[,c("split_by", "SO.split")])))), function(r) {

    b <- r[,"split_by"]
    c <- r[,"SO.split"]

    ref.rows <- Reduce(intersect, list(which(data[["id"]] %in% names(cells[which(cells == 1)])), which(data$split_by == b), which(data$SO.split == c)))
    freq.expr.by.split_by.SO <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    ref.rows <- Reduce(intersect, list(which(data[["id"]] %in% names(cells[which(cells == 1)])), which(data$split_by == b)))
    freq.expr.by.split <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    ref.rows <- Reduce(intersect, list(which(data[["id"]] %in% names(cells[which(cells == 1)])), which(data$SO.split == c)))
    freq.expr.by.SO <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    ref.rows <- which(data[["id"]] %in% names(cells[which(cells == 1)]))
    freq.expr <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    freqs <- c(freq.expr.by.split_by.SO, freq.expr.by.split, freq.expr.by.SO, freq.expr)

    freqs <- unlist(lapply(freqs, function(d) {
      if (is.na(d)) {
        d <- paste0(0, " %")
      } else {
        if (d < 1 & d > 0) {
          d <- "< 1 %"
        } else if (d > 1 & d < 99) {
          d <- paste0(round(d, 0), " %")
        } else if (d > 99 & d < 100) {
          d <- "> 99 %"
        } else {
          d <- paste0(d, " %")
        }
      }
    }))

    # get ranges of x and y axis to place annotation in individual facets
    ref.rows <- which(data$SO.split == c)
    xrng <- range(data[ref.rows,paste0(reduction, "_", dims[1])])
    yrng <- range(data[ref.rows,paste0(reduction, "_", dims[2])])

    return(data.frame(
      split_by = b,
      SO.split = c,
      freq.expr.by.split_by.SO = freqs[1],
      freq.expr.by.split = freqs[2],
      freq.expr.by.SO = freqs[3],
      freq.expr = freqs[4],
      xmin = xrng[1],
      xmax = xrng[2],
      ymin = yrng[1],
      ymax = yrng[2]
    ))
  }))

  return(freqs)
}


get_title <- function(feature_ex = NULL,
                      feature_cut = NULL,
                      feature_cut_expr = 0,
                      freq = NULL,
                      feature_italic = T,
                      feature = "FCMR",
                      name_anno = "{feature} ({freq}) in {feature_cut_ex}",
                      markdown = F) {
  if (!is.null(feature_ex)) {
    feature_ex <- paste(feature_ex, collapse = "&")
  }
  if (!is.null(feature_cut)) {
    feature_cut <- paste(feature_cut, collapse = "&")
  }
  if (length(feature_cut_expr) > 1) {
    feature_cut_expr <- paste(feature_cut_expr, collapse = "&")
  }

  if (markdown) {
    if (feature_italic) {feature <- paste0("*", feature, "*")}
    if (!is.null(feature_cut)) {feature_cut <- paste0("*", feature_cut, "*>", feature_cut_expr)} #else {feature_cut <- ""}
    if (!is.null(feature_ex)) {feature_ex <- paste0("*", feature_ex, "*=0")} #else {feature_ex <- ""}
    feature_cut_ex <- paste(c(feature_cut, feature_ex), collapse = " & ")

    title <- glue::glue(name_anno, .null = "")
    # replacing default left overs when feature_ex and feature_cut are NULL
    title <- gsub(" in ", "", title)
    title <- gsub("\\(\\)", "", title)
    title <- gsub(" {2,}", " ", title)
    title <- trimws(title)

    #title <- trimws(paste0(feature, feature_cut_ex_freq))
  } else {
    args <- list(x = feature, f = levels(as.factor(freqs$freq2)), y = feature_cut, z = feature_cut_expr, q = feature_ex)
    title <-
      if (!is.null(feature_cut) && !is.null(feature_ex)) {
        if (plot.freq.title && !is.null(freqs)) {
          if (feature_italic) {
            substitute(paste(italic(x), " in ", italic(y), ">", z, " & ", italic(q), " = 0", " (", f, ")", sep = ""), args)
          } else {
            substitute(paste(x, " in ", italic(y), ">", z, " & ", italic(q), " = 0", " (", f, ")", sep = ""), args)
          }
        } else {
          if (feature_italic) {
            substitute(paste(italic(x), " in ", italic(y), ">", z, " & ", italic(q), " = 0", sep = ""), args)
          } else {
            substitute(paste(x, " in ", italic(y), ">", z, " & ", italic(q), " = 0", sep = ""), args)
          }
        }
      } else if (!is.null(feature_cut)) {
        if (plot.freq.title && !is.null(freqs)) {
          if (feature_italic) {
            substitute(paste(italic(x), " in ", italic(y), ">", z, " (", f, ")", sep = ""), args)
          } else {
            substitute(paste(x, " in ", italic(y), ">", z, " (", f, ")", sep = ""), args)
          }
        } else {
          if (feature_italic) {
            substitute(paste(italic(x), " in ", italic(y), ">", z, sep = ""), args)
          } else {
            substitute(paste(x, " in ", italic(y), ">", z, sep = ""), args)
          }
        }
      } else if (!is.null(feature_ex)) {
        if (plot.freq.title && !is.null(freqs)) {
          if (feature_italic) {
            substitute(paste(italic(x), " in ", italic(q), " = 0", " (", f, ")", sep = ""), args)
          } else {
            substitute(paste(x, " in ", italic(q), " = 0", " (", f, ")", sep = ""), args)
          }
        } else {
          if (feature_italic) {
            substitute(paste(italic(x), " in ", italic(q), " = 0", sep = ""), args)
          } else {
            substitute(paste(x, " in ", italic(q), " = 0", sep = ""), args)
          }
        }
      } else {
        if (plot.freq.title && !is.null(freqs)) {
          if (feature_italic) {
            substitute(paste(italic(x), " (", f, ")", sep = ""), args)
          } else {
            substitute(paste(x, " (", f, ")", sep = ""), args)
          }
        } else {
          if (feature_italic) {
            substitute(paste(italic(x), sep = ""), args)
          } else {
            substitute(paste(x, sep = ""), args)
          }
        }
      }
    # use plot <- plot + ggplot2::labs(title = substitute(paste(x, sep = ""), list(x = title)))
  }

  return(title)
}


get_legend_text <- function(data) {

  if (is.numeric(data[["feature"]])) {
    if (all(data[which(data[["cells"]] == 1),"feature"] == 0)) {
      scale.max <- 0
      scale.min <- 0
      scale.mid <- 0
    } else {

      scale.max <- max(data[intersect(which(data[["cells"]] == 1), which(is.finite(data[["feature"]]))), "feature"], na.rm = T)
      scale.min <- min(data[Reduce(intersect, list(which(data[["feature"]] != 0), which(data[["cells"]] == 1), which(is.finite(data[["feature"]])))), "feature"], na.rm = T) # != 0 for module scores
      scale.mid <- scale.min + ((scale.max - scale.min) / 2)

      scale.max <- as.numeric(format(brathering::floor2(scale.max, 0.1), nsmall = 1))
      scale.min <- as.numeric(format(brathering::ceiling2(scale.min, 0.1), nsmall = 1))
      scale.mid <- as.numeric(format(round(scale.mid, 1), nsmall = 1))
    }
    return(list(scale.min, scale.mid, scale.max))
  }
  return(NULL)
}

get_col_pal <- function(data,
                        col_pal_c_args,
                        col_pal_d_args) {

  if (is.numeric(data[["feature"]])) {
    if (length(col_pal_c_args[["name"]]) == 1 && !col_pal_c_args[["name"]] %in% grDevices::colors()) {
      col.pal <- do.call(what = colrr::col_pal, args = col_pal_c_args)
    } else {
      col.pal <- col_pal_c_args[["name"]]
    }
  } else {
    if (length(col_pal_d_args[["name"]]) == 1 && !col_pal_d_args[["name"]] %in% grDevices::colors()) {
      col.pal <- do.call(what = colrr::col_pal, args = c(list(n = nlevels(as.factor(data[["feature"]]))), col_pal_d_args))
    } else {
      col.pal <- col_pal_d_args[["name"]]
    }
  }
  return(col.pal)
}


add_color_scale <- function(plot,
                            col.pal,
                            col_legend_args,
                            col_steps = "auto",
                            legendbreaks = "auto",
                            legendlabels = "auto",
                            col_steps_nice = T,
                            col_na = "grey50",
                            col_binary = F) {


  if (col_binary) {
    if (col_legend_args[["title"]] == "..auto..") {
      col_legend_args[["title"]] <- "UMI"
    }
    guide_fun <- ggplot2::guide_legend
  } else {

    if (is.numeric(plot[["data"]][["feature"]])) {
      if (col_legend_args[["title"]] == "..auto..") {
        if (attr(plot[["data"]], "feature_type") == "gene") {
          if (!is.null(attr(plot[["data"]], "layer"))) {
            if (attr(plot[["data"]], "layer") == "data") {
              col_legend_args[["title"]] <- "log (UMI)"
            } else if (attr(plot[["data"]], "layer") == "counts") {
              col_legend_args[["title"]] <- "UMI"
            } else {
              col_legend_args[["title"]] <- attr(plot[["data"]], "layer")
            }
          }
        } else if (attr(plot[["data"]], "feature_type") == "meta") {
          col_legend_args[["title"]] <- NA # omit by default as names shows up through names_anno
        }

      }
      c(scale.min, scale.mid, scale.max) %<-% get_legend_text(data = plot[["data"]])
      qmin <- attr(plot[["data"]], "qmin")
      qmax <- attr(plot[["data"]], "qmax")

      if (length(unique(c(scale.min, scale.mid, scale.max))) > 1) {

        # change this
        if (qmin > 0) {min.lab <- paste0(scale.min, " (q", round(qmin*100, 0), ")")} else {min.lab <- scale.min}
        if (qmax < 1) {max.lab <- paste0(scale.max, " (q", round(qmax*100, 0), ")")} else {max.lab <- scale.max}

        if (legendbreaks == "minmidmax") {
          legendbreaks <- c(scale.min, scale.mid, scale.max)
        }
        if (legendlabels == "auto" && (qmin > 0 || qmax < 1)) {
          legendlabels <- c(min.lab, scale.mid, max.lab)
        } else if (legendlabels != "auto" && (qmin > 0 || qmax < 1)) {
          message("qmax and/or qmin nor shown in legend.")
        }
        scale_color <- fcexpr:::colorscale_heuristic(colorscale_values = plot[["data"]][["feature"]],
                                                     values_zscored = F,
                                                     colorsteps = col_steps,
                                                     legendbreaks = legendbreaks,
                                                     legendlabels = legendlabels,
                                                     fill = col.pal,
                                                     colorsteps_nice = col_steps_nice,
                                                     type = "color",
                                                     col_na = col_na,
                                                     qmin = attr(plot[["data"]], "qmin"),
                                                     qmax = attr(plot[["data"]], "qmax"),
                                                     scale.min = scale.min)

        if (grepl("coloursteps", scale_color[["guide"]])) {
          guide_fun <- ggplot2::guide_colorsteps
        } else {
          guide_fun <- ggplot2::guide_colorbar
        }

        plot <- plot + scale_color
      } else {
        # if no expressers are found: breaks and labels would be of different lengths
        plot <- plot + ggplot2::scale_color_gradientn(colors = col.pal,
                                                      limits = c(scale.min, scale.max),
                                                      breaks = c(scale.min, scale.mid, scale.max),
                                                      na.value = col_na)
        guide_fun <- ggplot2::guide_colorbar
      }

    } else {

      if (col_legend_args[["title"]] == "..auto..") {
        #col_legend_args[["title"]] <- attr(plot[["data"]], "feature")
        col_legend_args[["title"]] <- NA # or ""; omit by default as it shows up in annotation or title, may take too much space as legend title
      }
      if ("barheight" %in% names(col_legend_args)) {
        col_legend_args[["barheight"]] <- 1
      }
      plot <- plot + ggplot2::scale_color_manual(values = col.pal,
                                                 na.value = col_na)
      guide_fun <- ggplot2::guide_legend
    }
  }


  plot <- plot + ggplot2::guides(color = Gmisc::fastDoCall(guide_fun, args = col_legend_args))

  return(plot)
}

add_facet <- function(plot,
                      facet_grid_row_var = NULL,
                      facet_scales = c("fixed", "free", "free_x", "free_y"),
                      nrow_inner = NULL,
                      ncol_inner = NULL) {
  if (!is.null(ncol_inner) && !is.null(nrow_inner)) {
    message("ncol_inner set NULL as ncol_inner and nrow_inner are provided.")
    ncol_inner <- NULL
  }
  facet_scales <- rlang::arg_match(facet_scales)
  if (is.null(facet_grid_row_var)) {
    wrap_by <- function(...) {ggplot2::facet_wrap(
      ggplot2::vars(...),
      labeller = ggplot2::label_wrap_gen(multi_line = F),
      scales = facet_scales,
      nrow = nrow_inner,
      ncol = ncol_inner
    )}
  } else {
    plot[["data"]][["facet_grid_row"]] <- facet_grid_row_var
    wrap_by <- function(...) {ggplot2::facet_grid(
      cols = ggplot2::vars(...),
      rows = ggplot2::vars(facet_grid_row),
      labeller = ggplot2::label_wrap_gen(multi_line = F),
      scales = facet_scales
    )}
  }

  if ("split_feature" %in% names(attributes(plot[["data"]]))) {
    split_by <- "split_feature"
  } else {
    split_by <- NULL
  }

  if (nlevels(plot[["data"]][["SO.split"]]) == 1 && !is.null(split_by)) {
    plot <- plot + wrap_by(!!rlang::sym(split_by))
  } else if (nlevels(plot[["data"]][["SO.split"]]) > 1 && is.null(split_by)) {
    plot <- plot + wrap_by(SO.split)
  } else if (nlevels(plot[["data"]][["SO.split"]]) > 1 && !is.null(split_by)) {
    plot <- plot + wrap_by(SO.split, !!rlang::sym(split_by))
  } else {
    # no facetting
  }
  return(plot)
}

add_axes_expansion <- function(plot,
                               axes_lim_set = list(),
                               axes_lim_expand = list()) {

  if (length(axes_lim_set) > 0 || length(axes_lim_expand) > 0) {

    if (length(axes_lim_expand) > 0) {
      if (length(axes_lim_set) > 0) {
        message("axes_lim_expand provided, axes_lim_set ignored.")
      }
      if (is.null(names(axes_lim_expand))) {
        message("axes_lim_expand needs names, e.g. list(x = c(2,2), y = c(3,3))")
        axes_lim_expand <- NULL
      }
      axes_lims <- brathering::gg_lims(plot)
      if (any(!names(axes_lim_expand) %in% names(axes_lims))) {
        message("names of axes_lim_expand must be: ", paste(names(axes_lims), collapse = ","))
        axes_lim_expand <- axes_lim_expand[names(axes_lims)]
      }
      if (!is.null(axes_lim_expand)) {
        for (i in names(axes_lims)) {
          axes_lims[[i]] <- axes_lims[[i]] + axes_lim_expand[[i]]
        }
      }
      axes_lim_set <- axes_lims
    }
    plot <- plot + Gmisc::fastDoCall(ggplot2::expand_limits, args = axes_lim_set)
  }
  return(plot)
}


add_labels <- function(plot = plot,
                       reduction = "UMAP",
                       dims = c(1,2),
                       label_filter_cells = T,
                       label_center_fun = c("median", "mean"),
                       label_nudge = c(0,0),
                       label_repel = F,
                       label_multi_try = F,
                       label_multi_max = 3,
                       label_args = list(
                         label.colour = NA,
                         fill = "white",
                         size = 4,
                         color = "black",
                         label.padding = ggplot2::unit(rep(0.1,4), "lines")),
                       finalize_plotting = F) {

  data <- plot[["data"]] #|> tidyr::pivot_wider(names_from = feat, values_from = value)
  label_feature <- ifelse("label_feature" %in% names(attributes(plot[["data"]])), "label_feature", "feature")

  if (is.numeric(data[[label_feature]])) {
    message("Labels not plotted as ", label_feature, " is numeric.")
    return(plot)
  }

  label_center_fun <- rlang::arg_match(label_center_fun)
  label_center_fun <- match.fun(label_center_fun)

  dimcol1 <- attr(plot[["data"]], "dim1")
  dimcol2 <- attr(plot[["data"]], "dim2")

  # plot labels only on included cells
  if (label_filter_cells) {
    data <- dplyr::filter(data, cells == 1)
  }

  label_df_multi <- data.frame()
  if (label_multi_try) {
    # low p: multimodality
    dip_p <-
      data |>
      dplyr::group_by(!!rlang::sym(label_feature), SO.split) |>
      dplyr::slice_sample(n = 5e4) |>
      dplyr::summarise(
        xp = diptest::dip.test(.data[[dimcol1]])$p.value,
        yp = diptest::dip.test(.data[[dimcol2]])$p.value,
        .groups = "drop") |>
      dplyr::filter(xp < 0.01 | yp < 0.01)


    plotranges <- purrr::map(brathering::gg_lims(plot), diff)
    names(plotranges) <- c(dimcol1, dimcol2)

    ## calculate separate avg coordinates for multi clusters
    # collapse them if too few cells or if too close
    dtach <- !"mclust" %in% .packages()
    null <- capture.output(library(mclust))
    label_df_multi <- purrr::pmap_dfr(.l = asplit(dip_p, 2), .f = function(label_feature, SO.split, xp, yp) {
      set.seed(as.numeric(xp))
      datasub <-
        data |>
        dplyr::filter(SO.split %in% !!SO.split & label_feature %in% !!label_feature) |>
        dplyr::slice_sample(n = 5000)
      xp <- as.numeric(xp)
      yp <- as.numeric(yp)

      if (xp < 0.01 && yp >= 0.01) {
        dimcolavgs <- get_dim_avg_multi(
          data = datasub,
          dim1 = dimcol1,
          dim2 = dimcol2,
          plotranges = plotranges,
          label_center_fun = label_center_fun
        )
        dimcol1_avg <- dimcolavgs[[1]]
        dimcol2_avg <- dimcolavgs[[2]]
      }
      if (yp < 0.01 && xp >= 0.01) {
        dimcolavgs <- get_dim_avg_multi(
          data = datasub,
          dim1 = dimcol2,
          dim2 = dimcol1,
          plotranges = plotranges,
          label_center_fun = label_center_fun
        )
        dimcol1_avg <- dimcolavgs[[2]]
        dimcol2_avg <- dimcolavgs[[1]]
      }
      if (xp < 0.01 && yp < 0.01) {
        mcl <- mclust::Mclust(
          dimcol12_data <- datasub[,c(dimcol1, dimcol2)],
          G = 2:label_multi_max,
          modelNames = "VVV",
          verbose = F
        )
        #dimcol12_data$clust <- mcl[["classification"]]
        #brathering::plot2(dimcol12_data, color = "clust")
        dimcol12_data_split <- split(dimcol12_data, mcl[["classification"]])
        fractions <- purrr::map_int(dimcol12_data_split,nrow)/nrow(dimcol12_data)
        dimcol12_data_split <- dimcol12_data_split[which(fractions>=0.2)]
        if (!length(dimcol12_data_split)) {
          dimcol12_data_split <- list(dimcol12_data)
        }
        dimcol1_avg <- purrr::map_dbl(dimcol12_data_split, ~label_center_fun(.x[[dimcol1]]))
        dimcol2_avg <- purrr::map_dbl(dimcol12_data_split, ~label_center_fun(.x[[dimcol2]]))
        new_avg <- collapse_close_points(x = dimcol1_avg, y = dimcol2_avg, threshold = mean(0.2*unlist(plotranges)))
        dimcol1_avg <- new_avg[["x"]]
        dimcol2_avg <- new_avg[["y"]]
      }
      return(stats::setNames(data.frame(
        dimcol1_avg = dimcol1_avg,
        dimcol2_avg = dimcol2_avg,
        label_feature = label_feature,
        SO.split = SO.split), c(dimcol1, dimcol2, "label", "SO.split")))
    })
    if (dtach) {
      detach("package:mclust", unload = T)
    }
  }

  # not split up clusters
  label_df <-
    data |>
    dplyr::group_by(!!rlang::sym(label_feature), SO.split) |>
    dplyr::summarise(!!dimcol1 := label_center_fun(!!rlang::sym(dimcol1)),
                     !!dimcol2 := label_center_fun(!!rlang::sym(dimcol2)),
                     .groups = "drop") |>
    dplyr::rename("label" = !!rlang::sym(label_feature)) |>
    as.data.frame()

  # join not split up and multi
  if (nrow(label_df_multi) > 0) {
    label_df <-
      dplyr::anti_join(label_df, label_df_multi, by = c("SO.split", "label")) |>
      dplyr::bind_rows(label_df_multi)
  }
  rownames(label_df) <- NULL

  label_df <-
    label_df |>
    # make for optional get_repel_coords, aligned with optional contour_labels
    dplyr::mutate("feature" = label) |>
    # reorder for optional get_repel_coords
    dplyr::select(!!rlang::sym(dimcol1), !!rlang::sym(dimcol2), label, feature, SO.split)

  if (!finalize_plotting) {
    assign("label_label", label_df, env = parent.frame())
  } else {
    if (label_repel) {
      label_df <- purrr::map_dfr(split(label_df, label_df$SO.split), brathering::get_repel_coords) #brathering::
    }
    label_nudge <- check_nudge_list(nudge_list = label_nudge,
                                    label_df = label_df)

    # nudge labels by name
    for (i in names(label_nudge)) {
      label_df[which(label_df$label == i),dimcol1] <- label_df[which(label_df$label == i),dimcol1] + label_nudge[[i]][1]
      label_df[which(label_df$label == i),dimcol2] <- label_df[which(label_df$label == i),dimcol2] + label_nudge[[i]][2]
    }

    plot <- plot +
      Gmisc::fastDoCall(what = ggtext::geom_richtext,
                        args = c(list(
                          data = label_df,
                          mapping = ggplot2::aes(label = label)), label_args))
  }

  return(plot)
}

get_dim_avg_multi <- function(data,
                              dim1,
                              dim2,
                              plotranges,
                              min_split_frac = 0.2,
                              min_range_frac = 0.2,
                              label_center_fun) {

  null <- capture.output(library(mclust))
  mcl <- mclust::Mclust(data_dim1 <- data[[dim1]], verbose = F)
  cluster_split <- split(data_dim1, mcl[["classification"]])
  # filter low fraction splits
  fractions <- lengths(cluster_split)/length(data_dim1)
  cluster_split <- cluster_split[which(fractions>=min_split_frac)]
  if (!length(cluster_split)) {
    cluster_split <- list(data_dim1)
  }
  dim1_avg <- purrr::map_dbl(cluster_split, label_center_fun)
  # filter clusters in close proximity
  dim1_avg <- collapse_close_values(dim1_avg, min_range_frac*plotranges[[dim1]])
  dim2_avg <- rep(label_center_fun(data[[dim2]]), length(dim1_avg))
  return(list(dim1_avg, dim2_avg))
}

collapse_close_values <- function(vec, threshold) {
  if (length(vec) <= 1) return(vec)

  # Compute the distance matrix
  d <- stats::dist(vec)

  # Perform hierarchical clustering
  hc <- stats::hclust(d, method = "single")

  # Cut the tree into clusters with max distance < threshold
  clusters <- stats::cutree(hc, h = threshold)

  # Collapse each cluster to its mean
  collapsed <- tapply(vec, clusters, mean)

  # Return the collapsed vector
  as.numeric(collapsed)
}

collapse_close_points <- function(x, y, threshold, return_cluster = F) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  if (length(x) <= 1) return(data.frame(x = x, y = y))

  # Combine x and y into a matrix of points
  coords <- cbind(x, y)

  # Compute Euclidean distances
  d <- stats::dist(coords)

  # Perform hierarchical clustering
  hc <- stats::hclust(d, method = "complete")

  # Cut tree at distance = threshold to form clusters
  clusters <- stats::cutree(hc, h = threshold)

  if (return_cluster) {
    return(clusters)
  }
  # Collapse each cluster to mean x and mean y
  collapsed <- stats::aggregate(coords, by = list(cluster = clusters), FUN = mean)

  # Return result as data.frame
  result <- collapsed[, c("x", "y")]
  return(result)
}

add_contour <- function(plot,
                        label_center_fun = c("median", "mean"),
                        contour_rm_outlier = F,
                        contour_rm_lowfreq_subcluster = F,
                        contour_filter_cells = T,
                        contour_col_pal_args = list(name = "custom"),
                        contour_args = list(contour_var = "ndensity", breaks = 0.3, linewidth = 1),
                        contour_label_nudge = c(0,0),
                        contour_label_args = list(size = 4),
                        contour_same_across_split = T,
                        contour_expr_freq = F,
                        contour_ggnewscale = F,
                        contour_multi_try = F,
                        contour_multi_max = 3,
                        contour_fun = ggplot2::geom_density2d,
                        contour_path_label = NULL,
                        finalize_plotting_expr_freq_labels = F) {
  dtach <- !"mclust" %in% .packages()
  null <- capture.output(library(mclust))
  # use all cells across split_by to assign contours? (in other words: same contour in every facet based on all cells, irrespective of split_by)

  label_center_fun <- rlang::arg_match(label_center_fun)
  label_center_fun <- match.fun(label_center_fun)

  dimcol1 <- attr(plot[["data"]], "dim1")
  dimcol2 <- attr(plot[["data"]], "dim2")

  data <- NULL
  try(expr = {
    # when contours were calculated for another feature already
    # parent frame 8 should be feature plot 2
    data <- get("contour_data", envir = parent.frame(8))
  }, silent = T)

  if (is.null(data)) {
    data <- plot[["data"]]
    #contour_feature <- attr(plot[["data"]], "contour_feature")
    if (contour_filter_cells) {
      data <- dplyr::filter(data, cells == 1)
    }

    data <- dplyr::mutate(data, group_id = as.character(dplyr::cur_group_id()), .by = c(contour_feature, SO.split))

    # filter 2D outliers by contour_feature level
    # instead of filtering below
    if (contour_rm_outlier) {
      # if min 5 algos indicate an outlier it is removed
      data <- purrr::map_dfr(split(data, data$group_id),
                             ~.x[which(rowSums(brathering::outlier(.x[,c(dimcol1, dimcol2)], methods = "dbscan")) == 0),])
    }

    data2 <- data.frame()
    if (contour_multi_try) {
      # check for multimodal clusters: first with diptest
      # low p: multimodality
      dip_p <-
        data |>
        dplyr::group_by(contour_feature, SO.split) |>
        dplyr::slice_sample(n = 5e4) |>
        dplyr::summarise(
          xp = diptest::dip.test(.data[[dimcol1]])$p.value,
          yp = diptest::dip.test(.data[[dimcol2]])$p.value,
          .groups = "drop") |>
        dplyr::filter(xp < 0.01 | yp < 0.01)

      plotranges <- purrr::map(brathering::gg_lims(plot), diff)
      names(plotranges) <- c(dimcol1, dimcol2)

      ## generate a separate id for split clusters to make separate contours
      # collapse them if too few cells or if too close
      data2 <- purrr::pmap_dfr(.l = asplit(dip_p, 2), .f = function(contour_feature, SO.split, xp, yp) {
        set.seed(as.numeric(xp))
        datasub <-
          data |>
          dplyr::filter(contour_feature %in% !!contour_feature & SO.split %in% !!SO.split) |>
          dplyr::slice_sample(n = 5000)
        xp <- as.numeric(xp)
        yp <- as.numeric(yp)
        if (xp < 0.01 || yp < 0.01) {
          mcl <- mclust::Mclust(
            datasub[,c(dimcol1, dimcol2)],
            modelNames = "VVV",
            G = 2:contour_multi_max,
            verbose = F
          )
          datasub_split <- split(datasub, mcl[["classification"]])

          # filter low fraction splits
          # these data will be missing for contour calculation
          # maybe remove that if filtering outliers above works
          if (contour_rm_lowfreq_subcluster) {
            ncells <- purrr::map_dbl(datasub_split, nrow)
            fractions <- ncells/nrow(datasub)
            datasub_split <- datasub_split[which(fractions>=0.2)]
            if (!length(datasub_split)) {
              datasub_split <- list(datasub)
            }
          }

          dimcol1_avg <- purrr::map_dbl(datasub_split, ~label_center_fun(.x[[dimcol1]]))
          dimcol2_avg <- purrr::map_dbl(datasub_split, ~label_center_fun(.x[[dimcol2]]))
          clusters <- stack(collapse_close_points(
            x = dimcol1_avg,
            y = dimcol2_avg,
            threshold = mean(0.2*unlist(plotranges)),
            return_cluster = T
          ))
          datasub_split <- purrr::map(stats::setNames(unique(clusters$values), unique(clusters$values)), function(x) {
            dplyr::bind_rows(datasub_split[clusters[which(clusters$values == x), "ind", drop = T]])
          })
          datasub_split <-
            dplyr::bind_rows(datasub_split, .id = "split") |>
            dplyr::mutate(contour_feature_split = paste0(contour_feature, "_", split))
        }
        return(datasub_split)
      })
    }

    # join rows of unimodal and multimodal clusters
    data <- dplyr::bind_rows(data |>
                               dplyr::filter(!contour_feature %in% unique(data2[["contour_feature"]])) |>
                               dplyr::mutate(contour_feature_split = NA),
                             data2) |>
      dplyr::mutate(contour_feature_split = ifelse(
        is.na(contour_feature_split),
        as.character(contour_feature),
        contour_feature_split
      ))

    # assign for next feature
    assign("contour_data", data, envir = parent.frame(8))
  }

  if (dtach) {
    detach("package:mclust", unload = T)
  }

  ## keep this as a list? or spare it?
  if (contour_same_across_split) {
    data <- purrr::map_dfr(stats::setNames(unique(data[["split_feature"]]), unique(data[["split_feature"]])),
                           function(x) data |> dplyr::mutate(split_by = x))
  }

  ## handle names of contour_col_pal
  if (is.factor(data[["contour_feature"]])) {
    colpal_names <- levels(data[["contour_feature"]])
  } else {
    colpal_names <- sort(as.character(unique(data[["contour_feature"]])))
  }
  if (length(contour_col_pal_args[["name"]]) == 1 && !contour_col_pal_args[["name"]] %in% grDevices::colors()) {
    # colors from col_pal and assignment to levels by order
    contour_col_pal <- do.call(colrr::col_pal, args = c(list(n = nlevels(as.factor(data[["contour_feature"]]))), contour_col_pal_args))
    contour_col_pal <- stats::setNames(contour_col_pal, colpal_names)
  } else {
    # colors explicitly provided
    contour_col_pal <- contour_col_pal_args[["name"]]
    if (length(contour_col_pal) != length(unique(data[["contour_feature"]]))) {
      message("contour_col_pal must have same length as groups to draw contours for.")
      contour_col_pal <- do.call(colrr::col_pal, args = c(list(n = nlevels(as.factor(data[["contour_feature"]])), name = "custom"), contour_col_pal_args[-which(names(contour_col_pal_args) == "name")]))
      contour_col_pal <- stats::setNames(contour_col_pal, colpal_names)
    }
    if (is.null(names(contour_col_pal)) || length(intersect(names(contour_col_pal), colpal_names)) != length(names(contour_col_pal))) {
      message("names of contour_col_pal do not match elements of contour_feature")
      contour_col_pal <- stats::setNames(contour_col_pal, colpal_names)
    }
  }

  # more than one list to contour_args (different settings for different contour for factor levels of contour_feature)
  # only when contour_ggnewscale is FALSE

  if (list_depth(contour_args) == 1 && !contour_ggnewscale) {
    contour_args <- rep(list(contour_args), length(unique(data[["contour_feature"]])))

    for (i in 1:length(lengths(contour_args))) {
      contour_args[[i]] <- contour_args[[i]][which(!names(contour_args[[i]]) %in% c("data", "mapping"))]
      if (!"contour_var" %in% names(contour_args[[i]])) {
        contour_args[[i]] <- c(contour_var = "ndensity", contour_args[[i]])
      }
      if (!"breaks" %in% names(contour_args[[i]])) {
        contour_args[[i]] <- c(breaks = 0.3, contour_args[[i]])
      }
      if (!"linewidth" %in% names(contour_args[[i]])) {
        contour_args[[i]] <- c(linewidth = 1, contour_args[[i]])
      }
    }
    names(contour_args) <- colpal_names
  } else if (list_depth(contour_args) > 1 && length(lengths(contour_args)) != length(unique(data[["contour_feature"]])) && !contour_ggnewscale) {
    message("Length of list of contour_args lists does not match length of factor levels of contour_feature.")
    contour_args <- rep(list(contour_args[[1]]), length(unique(data[["contour_feature"]])))
    names(contour_args) <- colpal_names
  } else if (!contour_ggnewscale) {
    if (is.null(names(contour_args)) || length(intersect(names(contour_args), colpal_names)) != length(names(contour_args))) {
      message("list of contour_args should have matching names of factor levels of contour_feature.")
      names(contour_args) <- colpal_names
    }
  } else if (contour_ggnewscale && list_depth(contour_args) > 1) {
    message("If contour_ggnewscale is TRUE only one list of contour_args is needed.")
    contour_args <- contour_args[[1]]
  }

  # calc expression freq per cluster
  if ((contour_expr_freq || !identical(ggplot2::geom_density_2d, contour_fun)) &&
      is.numeric(data[["feature"]]) &&
      !is.null(contour_path_label)) {
    # filter for largest subcluster of multimodal clusters
    # plot[["data"]] is unfiltered

  }

  if (contour_ggnewscale) {
    plot <-
      plot +
      ggnewscale::new_scale_color() +
      Gmisc::fastDoCall(ggplot2::geom_density2d, args = c(contour_args, list(data = data, mapping = ggplot2::aes(color = !!rlang::sym(contour_feature))))) +
      ggplot2::scale_color_manual(values = contour_col_pal)
  } else {

    for (i in unique(data[["contour_feature"]])) {
      for (j in unique(data[which(data[["contour_feature"]] == i), "contour_feature_split",drop=T])) {
        if (!identical(ggplot2::geom_density_2d, contour_fun) && !is.null(contour_path_label)) {
          # geomtextpath::geom_labeldensity2d does not allow for arbitrary label, hence: calculate paths by hand with
          # brathering::contour_est and use ...
          # ... geomtextpath::geom_textpath or geomtextpath::geom_labelpath
          # expr_by_cluster: in case this is used as path variable
          expr_by_cluster <-
            plot[["data"]] |>
            dplyr::group_by(contour_feature, SO.split) |>
            dplyr::summarise(pct = sum(feature > 0)/dplyr::n()*100, .groups = "drop") |>
            dplyr::mutate(pct = ifelse(pct < 1, ifelse(pct == 0, "0 %", "< 1 %"), paste0(round(pct,0), " %")))

          data2 <- data[which(data[["contour_feature_split"]] == j),,drop = F]
          contdata <- brathering::contour_est(as.data.frame(data2[,c(dimcol1, dimcol2)]))
          lvls <- unique(contdata[["level_norm"]])
          contdata <-
            dplyr::filter(contdata, level_norm == lvls[which.min(abs(contour_args[[i]][["breaks"]] - lvls))]) |>
            dplyr::bind_cols(data2[1,names(which(apply(data2, 2, function(x) length(unique(x))) == 1)), drop = F]) |>
            dplyr::left_join(expr_by_cluster, by = c("contour_feature", "SO.split"))
          suppressWarnings(expr = {
            plot <-
              plot +
              Gmisc::fastDoCall(contour_fun, args = c(
                contour_args[[i]],
                list(data = contdata,
                     mapping = ggplot2::aes(label = !!rlang::sym(contour_path_label)),
                     color = contour_col_pal[i])))
          })

        } else {
          #ggplot2::geom_density2d
          #contour_fun <- match.fun(contour_fun)
          plot <-
            plot +
            Gmisc::fastDoCall(contour_fun, args = c(
              contour_args[[i]],
              list(data = data[which(data[["contour_feature_split"]] == j),,drop = F],
                   color = contour_col_pal[i])))
        }

        # plot <-
        #   plot +
        #   Gmisc::fastDoCall(geomtextpath::geom_labeldensity2d, args = c(
        #     contour_args[[i]],
        #     list(data = data[which(data[["contour_feature_split"]] == j),,drop = F],
        #          color = contour_col_pal[i])))
      }
    }
  }

  if (contour_expr_freq && is.numeric(data[["feature"]])) {
    # filter for largest subcluster of multimodal clusters
    # plot[["data"]] is unfiltered
    # just recalculate for simplicity even though may have been done above
    expr_by_cluster <-
      plot[["data"]] |>
      dplyr::group_by(contour_feature, SO.split) |>
      # label = pct
      dplyr::summarise(label = sum(feature > 0)/dplyr::n()*100, .groups = "drop") |>
      dplyr::mutate(label = ifelse(label < 1, ifelse(label == 0, "0 %", "< 1 %"), paste0(round(label,0), " %")))

    # data may have been subject to filtering
    largest_subcluster <-
      data |>
      dplyr::group_by(contour_feature, SO.split, contour_feature_split) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::slice_max(order_by = n, n = 1, with_ties = F, by = c(contour_feature, SO.split))
    label_df <-
      data |>
      dplyr::semi_join(largest_subcluster, by = c("contour_feature_split", "SO.split")) |>
      dplyr::group_by(contour_feature, SO.split) |>
      dplyr::summarise(!!dimcol1 := label_center_fun(!!rlang::sym(dimcol1)),
                       !!dimcol2 := label_center_fun(!!rlang::sym(dimcol2)),
                       .groups = "drop") |>
      dplyr::left_join(expr_by_cluster, by = c("contour_feature", "SO.split")) |>
      dplyr::rename("feature" = contour_feature)

    if (!finalize_plotting_expr_freq_labels) {
      assign("label_contour", label_df, env = parent.frame())
    } else {
      contour_label_nudge <- check_nudge_list(nudge_list = contour_label_nudge,
                                              label_df = label_df)

      # nudge labels by name
      for (i in names(contour_label_nudge)) {
        label_df[which(label_df$label == i),dimcol1] <- label_df[which(label_df$label == i),dimcol1] + contour_label_nudge[[i]][1]
        label_df[which(label_df$label == i),dimcol2] <- label_df[which(label_df$label == i),dimcol2] + contour_label_nudge[[i]][2]
      }
      if (contour_ggnewscale) {
        plot <- plot +
          Gmisc::fastDoCall(ggtext::geom_richtext,
                            args = c(contour_label_args,
                                     list(data = label_df,
                                          show.legend = F,
                                          mapping = ggplot2::aes(
                                            x = !!rlang::sym(dimcol1),
                                            y = !!rlang::sym(dimcol2),
                                            label = label,
                                            color = !!rlang::sym(names(group_labels)[1])))))
      } else {
        plot <- plot +
          Gmisc::fastDoCall(ggtext::geom_richtext,
                            args = c(contour_label_args,
                                     list(data = label_df,
                                          mapping = ggplot2::aes(
                                            x = !!rlang::sym(dimcol1),
                                            y = !!rlang::sym(dimcol2),
                                            label = label))))
      }
    }
  }
  return(plot)
}

co_add_feature_and_contour_labels <- function(plot,
                                              label_label,
                                              label_contour,
                                              label_nudge,
                                              label_repel,
                                              label_args,
                                              contour_label_nudge,
                                              contour_label_args) {


  dimcol1 <- attr(plot[["data"]], "dim1")
  dimcol2 <- attr(plot[["data"]], "dim2")
  ### add nudging after repel
  if (!is.null(label_label) && !is.null(label_contour)) {
    # repel expr freq by contour and feature label
    label_df <- dplyr::bind_rows(label = label_label, contour = label_contour, .id = "group")
    label_df <- label_df[,c(2:ncol(label_df), 1)]
    label_df <- purrr::map_dfr(split(label_df, label_df$SO.split), brathering::get_repel_coords) #brathering::
    if ("split_feature" %in% names(plot[["data"]])) {
      label_df <- dplyr::left_join(
        label_df,
        dplyr::distinct(plot[["data"]], label_feature, split_feature),
        by = c("label" = "label_feature"))
    }
    label_nudge <- check_nudge_list(nudge_list = label_nudge,
                                    label_df = label_df)
    contour_label_nudge <- check_nudge_list(nudge_list = contour_label_nudge,
                                            label_df = label_df)
    for (i in names(label_nudge)) {
      label_df[which(label_df$label == i),dimcol1] <- label_df[which(label_df$label == i),dimcol1] + label_nudge[[i]][1]
      label_df[which(label_df$label == i),dimcol2] <- label_df[which(label_df$label == i),dimcol2] + label_nudge[[i]][2]
    }
    for (i in names(contour_label_nudge)) {
      label_df[which(label_df$label == i),dimcol1] <- label_df[which(label_df$label == i),dimcol1] + contour_label_nudge[[i]][1]
      label_df[which(label_df$label == i),dimcol2] <- label_df[which(label_df$label == i),dimcol2] + contour_label_nudge[[i]][2]
    }
    label_df <- split(label_df, label_df$group)
    plot <- plot +
      Gmisc::fastDoCall(what = ggtext::geom_richtext,
                        args = c(list(
                          data = label_df[["label"]],
                          mapping = ggplot2::aes(label = label)
                        ), label_args))
    plot <- plot +
      Gmisc::fastDoCall(what = ggtext::geom_richtext,
                        args = c(list(
                          data = label_df[["contour"]],
                          mapping = ggplot2::aes(label = label)
                        ), label_args))

  } else if (!is.null(label_label)) {

    label_df <- label_label
    if (label_repel) {
      label_df <- purrr::map_dfr(split(label_df, label_df$SO.split), brathering::get_repel_coords) #brathering::
    }

    if ("split_feature" %in% names(plot[["data"]])) {
      label_df <- dplyr::left_join(
        label_df,
        dplyr::distinct(plot[["data"]], label_feature, split_feature),
        by = c("label" = "label_feature"))
    }

    label_nudge <- check_nudge_list(nudge_list = label_nudge,
                                    label_df = label_df)
    # nudge labels by name
    for (i in names(label_nudge)) {
      label_df[which(label_df$label == i),dimcol1] <- label_df[which(label_df$label == i),dimcol1] + label_nudge[[i]][1]
      label_df[which(label_df$label == i),dimcol2] <- label_df[which(label_df$label == i),dimcol2] + label_nudge[[i]][2]
    }

    plot <- plot +
      Gmisc::fastDoCall(what = ggtext::geom_richtext,
                        args = c(list(
                          data = label_df,
                          mapping = ggplot2::aes(label = label)), label_args))

  } else if (!is.null(label_contour)) {

    label_df <- label_contour
    contour_label_nudge <- check_nudge_list(nudge_list = contour_label_nudge,
                                            label_df = label_df)
    for (i in names(contour_label_nudge)) {
      label_df[which(label_df$label == i),dimcol1] <- label_df[which(label_df$label == i),dimcol1] + contour_label_nudge[[i]][1]
      label_df[which(label_df$label == i),dimcol2] <- label_df[which(label_df$label == i),dimcol2] + contour_label_nudge[[i]][2]
    }

    plot <-
      plot +
      Gmisc::fastDoCall(what = ggtext::geom_richtext,
                        args = c(list(
                          data = label_contour,
                          mapping = ggplot2::aes(label = label)
                        ), contour_label_args))
  }
  return(plot)
}


check_nudge_list <- function(nudge_list,
                             label_df) {
  name <- deparse(substitute(nudge_list))
  if (!methods::is(nudge_list, "list") || (is.null(names(nudge_list)) && length(nudge_list) > 0)) {
    #nudge_list <- list(nudge_list)
    message(name, " has to be a named list with names being those of labels.")
    nudge_list <- replicate(length(label_df[["feature"]]), c(0,0), simplify = F)
    names(nudge_list) <- label_df[["feature"]]
  }
  if (any(lengths(nudge_list) != 2)) {
    message(name, " should contain vectors of length 2 only.")
    nudge_list <- nudge_list[which(lengths(nudge_list) == 2)]
  }
  if (length(nudge_list) > 0) {
    len_bef <- length(nudge_list)
    nudge_list <- nudge_list[which(names(nudge_list) %in% label_df[["feature"]])]
    if (length(nudge_list) < len_bef) {
      message(name, ": not all names found in label column of data.")
    }
  }
  return(nudge_list)
}

list_depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,list_depth,thisdepth=thisdepth+1))))
  }
}


check.features <- function(SO,
                           features,
                           rownames = T,
                           meta.data = T,
                           meta.data.numeric = F) {

  features <- unique(features)
  features <- features[which(!is.na(features))]
  features <- trimws(features)
  features <- features[which(features != "")]

  if (!length(features)) {return(NULL)}
  if (is.null(features)) {return(NULL)}
  if (!rownames && !meta.data) {return((NULL))}
  if (!is.list(SO)) {
    SO <- list(SO)
  }

  features.out <- feat.check(SO = SO,
                             features = features,
                             rownames = rownames,
                             meta.data = meta.data,
                             ignore.case = T,
                             meta.data.numeric = meta.data.numeric)

  if (length(features.out) == 0) {
    stop("Non of the provided features has not been found in every SO. No features left to plot.")
  }

  ### TO DO WITH meta.data.numeric
  hits <- hit.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = T)
  if (any(hits > 1)) {
    message(paste(names(hits)[which(hits > 1)], collapse = ","), " found more than once in at least one SO when ignoring case. So, case is being considered.")
    features.out <- feat.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = F)
    if (length(features.out) == 0) {stop("Non of the provided features has not been found in every SO. No features left to plot.")}
    hits <- hit.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = F)
    if (any(hits > 1)) {
      stop(paste0(paste(names(hits)[which(hits > 1)], collapse = ","), " still found more than once in at least one SO when not ignoring case. Please fix this (check rownames and meta.data)."))
    } else {
      ignore.case <- F
    }
  } else {
    ignore.case <- T
  }

  if (length(features.out) < length(features)) {message("Features not found in every SO: ", paste(setdiff(features, features.out), collapse = ","))}
  if (length(features.out) > length(features)) {stop("More features after check.features than before?!")}
  return(features.out)
}

hit.check <- function(SO, features, rownames, meta.data, ignore.case) {
  if (ignore.case) {
    sapply(features, function(x) {
      max(sapply(SO, function(y) {
        if (rownames & !meta.data) {
          length(which(tolower(x) == tolower(rownames(y))))
        } else if (!rownames & meta.data) {
          length(which(tolower(x) == tolower(names(y@meta.data))))
        } else if (rownames & meta.data) {
          length(which(tolower(x) == c(tolower(rownames(y)), tolower(names(y@meta.data)))))
        }
      }))
    })
  } else {
    sapply(features, function(x) {
      max(sapply(SO, function(y) {
        if (rownames & !meta.data) {
          length(which(x == rownames(y)))
        } else if (!rownames & meta.data) {
          length(which(x == names(y@meta.data)))
        } else if (rownames & meta.data) {
          length(which(x == c(rownames(y), names(y@meta.data))))
        }
      }))
    })
  }

}


feat.check <- function(SO, features, rownames, meta.data, ignore.case, meta.data.numeric) {
  unlist(lapply(features, function(x) {
    Reduce(intersect, lapply(SO, function(y) {
      if (rownames && !meta.data && !meta.data.numeric) {
        search <- rownames(y)
      } else if (rownames && !meta.data && meta.data.numeric) {
        search <- c(rownames(y), unique(c(names(which(sapply(y@meta.data, is.numeric))), names(which(sapply(y@meta.data, is.logical))))))
      } else if (!rownames && meta.data) {
        search <- names(y@meta.data)
      } else if (rownames && meta.data) {
        search <- c(names(y@meta.data), rownames(y))
      }
      grep(paste0("^",x,"$"), search, value = T, ignore.case = ignore.case)
    }))
  }))
}
