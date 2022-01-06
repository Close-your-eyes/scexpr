#' Title
#'
#' @param SO
#' @param features
#' @param assay
#' @param dims
#' @param cells
#' @param downsample
#' @param make.cellnames.unique
#' @param pt.size
#' @param pt.size.expr.factor
#' @param order
#' @param min.q.cutoff
#' @param max.q.cutoff
#' @param reduction
#' @param split.by
#' @param shape.by
#' @param combine
#' @param ncol.combine
#' @param nrow.combine
#' @param feature.aliases
#' @param title
#' @param cutoff.feature
#' @param cutoff.expression
#' @param exclusion.feature
#' @param binary.expr
#' @param col.expresser
#' @param legend.position
#' @param legend.title.text.size
#' @param legend.barheight
#' @param legend.barwidth
#' @param legend.text.size
#' @param legend.nrow
#' @param legend.ncol
#' @param legend.key.size
#' @param shape.legend.size
#' @param hide.shape.legend
#' @param col.legend.size
#' @param col.pal.c
#' @param col.pal.d
#' @param col.excluded.cells
#' @param col.non.expresser
#' @param col.pal.rev
#' @param theme
#' @param plot.axis.labels
#' @param plot.panel.grid
#' @param title.font.size
#' @param plot.freq.of.expr.title
#' @param plot.freq.of.expr.annotation
#' @param plot.legend.title
#' @param plot.title
#' @param plot.order
#' @param annotation.position
#' @param annotation.font.size
#' @param plot.cutoff.feature.inset
#' @param inset.position
#' @param inset.size
#' @param strip.font.size
#' @param plot.strip.selection
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
feature_plot <- function(SO,
                         features,
                         assay = c("RNA", "SCT"),
                         dims = c(1,2),
                         cells = NULL,
                         downsample = 1,
                         make.cellnames.unique = F,
                         pt.size = 1,
                         pt.size.expr.factor = 1,
                         order = T,
                         min.q.cutoff = 0,
                         max.q.cutoff = 1,
                         reduction = "tSNE",
                         split.by = NULL,
                         shape.by = NULL,
                         combine = T,
                         ncol.combine = NULL,
                         nrow.combine = NULL,
                         feature.aliases = NULL,
                         title = NULL,

                         cutoff.feature = NULL,
                         cutoff.expression = 0,
                         exclusion.feature = NULL,

                         binary.expr = F,
                         col.expresser = "#E69F00",

                         legend.position = "right",
                         legend.title.text.size = 14,
                         legend.barheight = 3,
                         legend.barwidth = 0.5,
                         legend.text.size = 10,
                         legend.nrow = NULL,
                         legend.ncol = NULL,
                         legend.key.size = 0.3,
                         shape.legend.size = 3,
                         hide.shape.legend = F,
                         col.legend.size = 3,

                         col.pal.c = "spectral",
                         col.pal.d = "custom",
                         col.excluded.cells = "#ededed",
                         col.non.expresser = "#cfcfcf",
                         col.pal.rev = F,
                         theme = NULL, #theme_bw(), #### 2021 01 05 change that
                         plot.axis.labels = F,
                         plot.panel.grid = F,
                         title.font.size = 14,
                         plot.freq.of.expr.title = T,
                         plot.freq.of.expr.annotation = T,
                         plot.legend.title = F,
                         plot.title = T,
                         plot.order = NULL,
                         annotation.position = c(-Inf, Inf),
                         annotation.font.size = 4,

                         plot.cutoff.feature.inset = F,
                         inset.position = c(0.75,0.1),
                         inset.size = c(0.15,0.15),

                         strip.font.size = 14,
                         plot.strip.selection = NA,
                         ...) {

  # tidy eval syntax: https://rlang.r-lib.org/reference/nse-force.html https://ggplot2.tidyverse.org/reference/aes.html#quasiquotation
  # numeric but discrete data columns from meta.data - how to tell that it is not continuous


  if (missing(features) | missing(SO)) {stop("Seurat object list or feature vector is missing.")}
  if (!class(features) %in% c("factor", "character")) {stop("Please provide a character vector for features.")}
  if (!is.null(ncol.combine) && !is.null(nrow.combine)) {stop("Please only select one, ncol.combine or nrow.combine.")}
  if (length(dims) != 2 | class(dims) != "numeric") {stop("dims has to be a numeric vector of length 2, e.g. c(1,2).")}
  assay <- match.arg(assay, c("RNA", "SCT"))
  if (!is.list(SO)) {SO <- list(SO)}
  if (is.null(names(SO))) {names(SO) <- as.character(seq_along(SO))}
  if (length(features) == 1) {combine <- F}


  # check if reduction is in all objects
  if (!all(unlist(lapply(SO, function(y) {tolower(reduction) %in% tolower(names(y@reductions))})))) {
    reduction <- base::Reduce(base::intersect, lapply(SO, function (x) {names(x@reductions)}))
    if (length(reduction) == 0) {stop("No common reduction found in SOs.")}
    if (length(reduction) > 1 && "tsne" %in% reduction) {
      reduction <- "tsne"
    } else if (length(reduction) > 1 && "umap" %in% reduction) {
      reduction <- "umap"
    } else {
      reduction <- sample(reduction, 1)
    }
  }

  SO <- lapply(SO, function(x) {
    if (!assay %in% Assays(x)) {stop("Assay not found in every SO.")}
    DefaultAssay(x) <- assay
    return(x)
  })

  if (!is.null(split.by) && length(.check.features(SO, split.by, rownames = F)) == 0) {stop("split.by not found in all objects.")}
  if (!is.null(shape.by) && length(.check.features(SO, shape.by, rownames = F)) == 0) {stop("shape.by not found in all objects.")}

  # removing features which could not be found
  features <- .check.features(SO, unique(features))
  if (!is.null(cutoff.feature) && length(.check.features(SO, cutoff.feature)) == 0) {stop("cutoff.feature not found in every SO.")}
  if (!is.null(exclusion.feature) && length(.check.features(SO, exclusion.feature)) == 0) {stop("exclusion.feature not found in every SO.")}

  if (length(cutoff.feature) != length(cutoff.expression) && length(cutoff.expression) != 1) {stop("Unequal lengths of cutoff.feature and cutoff.expression.")}
  if (length(cutoff.expression) == 1) {cutoff.expression <- rep(cutoff.expression, length(cutoff.feature))}

  # check if cell names are unique across SOs
  all.cells <- unlist(lapply(SO, function(x) {Seurat::Cells(x)}))
  if (length(SO) > 1 && !all(!duplicated(all.cells)) && !make.cellnames.unique) {stop("Cell names are not unique across SOs. Please fix that manually or pass make.cellnames.unique = T when calling feature.plot. Cells are then renamed with the prefix SO_i_. Where i the index of the SO in the list. Consider this when selecting cells for plotting.")}
  if (length(SO) > 1 && !all(!duplicated(all.cells)) && make.cellnames.unique) {
    names.temp <- names(SO)
    SO <- lapply(seq_along(SO), function(x) {
      Seurat::RenameCells(SO[[x]], add.cell.id = paste0("SO_", x))
      # RenameCells adds another "_" by itself
    })
    names(SO) <- names.temp
    all.cells <- unlist(lapply(SO, function(y) {Seurat::Cells(y)}))
  }
  if (is.null(cells)) {cells <- all.cells}

  # create a vector of cells which identifies how to plot them; 0 indicates exclusion
  cells.plot <- setNames(rep(1, length(all.cells)), all.cells)
  ## cells excluded by arbitrary selection
  cells.plot[!names(cells.plot) %in% cells] <- 0
  ## cells excluded by
  if (!is.null(cutoff.feature)) {
    for (y in seq_along(SO)) {
      for (z in seq_along(cutoff.feature)) {
        # test again with multiple SOs, I think okay though
        # the one-liner is only possible with '<=' (think of cutoff.expression == 0);
        #cells.plot[which(names(cells.plot) %in% names(which(Seurat::GetAssayData(SO[[y]], slot = "data", assay = assay)[cutoff.feature[z],] <= cutoff.expression[z])))] <- 0
        cells.plot[intersect(which(!names(cells.plot) %in% names(which(Seurat::GetAssayData(SO[[y]], slot = "data", assay = assay)[cutoff.feature[z],] > cutoff.expression[z]))),
                             which(names(cells.plot) %in% Seurat::Cells(SO[[y]])))] <- 0
      }
    }
  }


  ## cells excluded due to expression of an unwanted gene
  if (!is.null(exclusion.feature)) {
    for (y in seq_along(SO)) {
      for (z in exclusion.feature) {
        cells.plot[which(names(cells.plot) %in% names(which(Seurat::GetAssayData(SO[[y]], slot = "data", assay = assay)[z,] > 0)))] <- 0
      }
    }
  }

  # downsample
  if (downsample < 1) {
    downsample <- downsample*length(cells.plot)
  } else if (downsample == 1) {
    downsample <- length(cells.plot)
  }
  ### check!
  cells.downsample <- cells.plot[sample(seq_along(cells.plot), size = downsample)]

  # get data by features and objects and plot
  plots <- lapply(features, function(x) {
    data <- do.call(rbind, lapply(unname(SO), function(y) {

      if (x %in% rownames(Seurat::GetAssayData(y, slot = "data", assay = assay)) && !x %in% names(y@meta.data)) {
        data <- data.frame(t(as.matrix(Seurat::GetAssayData(y, slot = "data", assay = assay)[x,,drop = F])), check.names = F)
      } else if (!x %in% rownames(Seurat::GetAssayData(y, slot = "data", assay = assay)) && x %in% names(y@meta.data)) {
        data <- data.frame(y@meta.data[,x,drop=F], stringsAsFactors = F, check.names = F)
      } else {
        stop("Feature found in meta data and features.")
      }
      data <- cbind(data, Seurat::Embeddings(y, reduction = names(y@reductions)[grepl(reduction, names(y@reductions), ignore.case = T)]))
      data <- data[which(rownames(data) %in% names(cells.downsample)),]

      if (is.null(split.by)) {data[,"split.by"] <- "1"} else {data[,"split.by"] <- y@meta.data[,split.by]}
      if (!is.null(shape.by)) {data[,shape.by] <- as.factor(as.character(y@meta.data[,shape.by]))}

      return(data)
    }))

    names <- sapply(SO, function(y) {length(which(colnames(Seurat::GetAssayData(y, slot = "data", assay = assay)) %in% names(cells.downsample)))})
    data[,"SO.split"] <- rep(names(names), names)
    if (length(SO) > 1) {
      SO.split <- "SO.split"
    } else {
      SO.split <- NULL
    }

    shape.by <- tryCatch(sym(shape.by), error = function (e) NULL) # neccessary to make sym(shape.by) here, for !!shape.by to work; not possible within aes()

    if (is.numeric(data[,1])) {
      if (all(data[,1] == 0)) {print(paste0("No expresser found for ", x, "."))}
      data$binary.expr <- ifelse(data[,1] > 0, "+", "-")
      if (binary.expr & any(data[,1] < 0)) {
        print("binary.expr (+/-) may not be meaningful as there are negative expression values.")
      }
    }

    # correct reduction.name with respect to case
    reduction.name <- gsub("_[0-9]", "", names(data)[grepl(reduction, names(data), ignore.case = T)])[1]

    # use squishing to dampen extreme values - this will produce actually wrong limits on the legend
    if (is.numeric(data[,1])) {
      if (all(!data[,1] < 0, na.rm = T)) {
        # expression is always greater than 0 and non-expresser are excluded
        data[,1][which(data[,1] > 0)] <- scales::squish(data[,1][which(data[,1] > 0)], range = c(quantile(data[,1][which(data[,1] > 0)], min.q.cutoff, na.rm = T), quantile(data[,1][which(data[,1] > 0)], max.q.cutoff, na.rm = T)))
      } else {
        # e.g. for module scores below 0
        data[,1] <- scales::squish(data[,1], range = c(quantile(data[,1], min.q.cutoff, na.rm = T), quantile(data[,1], max.q.cutoff, na.rm = T)))
      }
    }

    # generate legend labels
    if (is.numeric(data[,1])) {
      if (all(data[,1] == 0)) {
        scale.max <- 0
        scale.min <- 0
        scale.mid <- 0
      } else {
        scale.max <- max(data[which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])),1])
        scale.min <- min(data[intersect(which(data[,1] != 0), which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)]))),1]) # != 0 for module scores
        scale.mid <- scale.min + ((scale.max - scale.min) / 2)

        scale.max <- as.numeric(format(ceiling_any(scale.max, 0.1), nsmall = 1))
        scale.min <- as.numeric(format(floor_any(scale.min, 0.1), nsmall = 1))
        scale.mid <- as.numeric(format(round(scale.mid, 1), nsmall = 1))
      }
    }

    # select color scale
    if (is.numeric(data[,1])) {
      if (length(col.pal.c) == 1) {
        scale.to.use <- col_pal(name = col.pal.c, reverse = col.pal.rev)
      } else {
        scale.to.use <- col.pal.c
      }
    } else {
      if (length(col.pal.d) == 1) {
        scale.to.use <- col_pal(name = col.pal.d, reverse = col.pal.rev, n = nlevels(as.factor(data[,1])))
      } else {
        scale.to.use <- col.pal.d
      }
    }

    # plot
    plot <- ggplot(data, aes(x = !!sym(paste0(reduction.name, "_", dims[1])), y = !!sym(paste0(reduction.name, "_", dims[2]))))
    # plot excluded cells
    plot <- plot + geom_point(aes(shape = !!shape.by), size = pt.size, colour = col.excluded.cells, data = data[which(rownames(data) %in% names(cells.plot[which(cells.plot == 0)])),])

    # different procedure for gene feature or meta.data feature
    if (all(unlist(lapply(unname(SO), function(y) {x %in% rownames(Seurat::GetAssayData(y, slot = "data", assay = assay))})))) {

      # calculate freq. of expresser by split.by variable (without dplyrs group.by)
      freqs <- do.call(rbind, lapply(split(unique(data[,c("split.by", "SO.split")]), seq(nrow(unique(data[,c("split.by", "SO.split")])))), function(r) {

        b <- r[,"split.by"]
        c <- r[,"SO.split"]

        ref.rows <- Reduce(intersect, list(which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])), which(data$split.by == b), which(data$SO.split == c)))
        freq.expr.by.split.by.SO <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

        ref.rows <- Reduce(intersect, list(which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])), which(data$split.by == b)))
        freq.expr.by.split <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

        ref.rows <- Reduce(intersect, list(which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])), which(data$SO.split == c)))
        freq.expr.by.SO <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

        ref.rows <- which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)]))
        freq.expr <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

        freqs <- c(freq.expr.by.split.by.SO, freq.expr.by.split, freq.expr.by.SO, freq.expr)

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
        xrng <- range(data[ref.rows,paste0(reduction.name, "_", dims[1])])
        yrng <- range(data[ref.rows,paste0(reduction.name, "_", dims[2])])

        data.frame(split.by = b, SO.split = c, freq.expr.by.split.by.SO = freqs[1], freq.expr.by.split = freqs[2], freq.expr.by.SO = freqs[3], freq.expr = freqs[4], xmin = xrng[1], xmax = xrng[2], ymin = yrng[1], ymax = yrng[2])
      }))

      aliases.list <- .check.aliases(x, feature.aliases, data)
      x <- aliases.list[[1]]
      data <- aliases.list[[2]]

      if (order) {
        data <- data[order(data[,1]),]
      }

      # plot expresser and non-expresser
      if (binary.expr) {
        plot <- plot + geom_point(data = data[which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])),], aes(shape = !!shape.by, colour = binary.expr), size = pt.size*pt.size.expr.factor) + scale_color_manual(values = c(col.non.expresser, col.expresser))
      } else {
        # non-expresser
        plot <- plot + geom_point(data = data[intersect(which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])), which(data[,1] == 0)),], aes(shape = !!shape.by), size = pt.size, colour = col.non.expresser)
        # expresser
        plot <- plot + geom_point(data = data[intersect(which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])), which(data[,1] > 0)),], aes(colour = !!sym(x), shape = !!shape.by), size = pt.size*pt.size.expr.factor)
      }

      if (length(SO) > 1) {
        if (is.null(split.by)) {label = "freq.expr.by.SO"} else {label = "freq.expr.by.split.by.SO"}
      } else if (!is.null(split.by) && length(SO) == 1) {
        label = "freq.expr.by.split"
      } else {
        label = "freq.expr"
      }
      if (plot.freq.of.expr.annotation) {
        plot <- plot + ggrepel::geom_text_repel(data = freqs, size = annotation.font.size, aes(label = !!sym(label), x = xmin + abs(xmin - xmax) * annotation.position[1], y = ymin + abs(ymin - ymax) * annotation.position[2]))
      }

      plot.colourbar <- !binary.expr
      make.italic <- T

    } else if (all(unlist(lapply(unname(SO), function(y) {x %in% names(y@meta.data)})))) {


      if (is.null(legend.ncol) && !is.null(legend.nrow)) {
        legend.ncol <- ceiling(length(unique(data[,1]))/legend.nrow)
      } else if (is.null(legend.nrow) && !is.null(legend.ncol)) {
        legend.nrow <- ceiling(length(unique(data[,1]))/legend.ncol)
      } else if (is.null(legend.nrow) && is.null(legend.ncol)) {
        legend.ncol <- 1
      } else {
        # error may pop up below
      }

      freqs <- NULL

      # feature alias
      aliases.list <- .check.aliases(x, feature.aliases, data)
      x <- aliases.list[[1]]
      data <- aliases.list[[2]]

      if (is.null(plot.order)) {
        # plot non-excluded cells, sample to make it random
        plot <- plot + geom_point(data = data[which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])),][sample(nrow(data[which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])),])), ], aes(colour = !!sym(x), shape = !!shape.by), size = pt.size)
      } else {
        names(scale.to.use) <- plot.order
        if (length(pt.size) == 1) {
          pt.size <- rep(pt.size, length(plot.order))
        }
        for (i in seq_along(plot.order)) {
          plot <- plot + geom_point(data = data[intersect(which(rownames(data) %in% names(cells.plot[which(cells.plot == 1)])), which(data[,x] == plot.order[i])),], aes(colour = !!sym(x), shape = !!shape.by), size = pt.size[i])
        }
      }

      plot.colourbar <- is.numeric(data[,1])
      make.italic <- F
    }

    if (!binary.expr) {
      if (is.numeric(data[,1])) {
        if (min.q.cutoff > 0) {min.lab <- paste0(scale.min, " (q", round(min.q.cutoff*100, 0), ")")} else {min.lab <- scale.min}
        if (max.q.cutoff < 1) {max.lab <- paste0(scale.max, " (q", round(max.q.cutoff*100, 0), ")")} else {max.lab <- scale.max}
        if (length(unique(c(scale.min, scale.mid, scale.max))) > 1) {
          plot <- plot + scale_colour_gradientn(colours = scale.to.use, limits = c(scale.min, scale.max), breaks = c(scale.min, scale.mid, scale.max), labels = c(min.lab,  scale.mid, max.lab))
        } else {
          # if no expresser are found: breaks and labels would be of different lengths
          plot <- plot + scale_colour_gradientn(colours = scale.to.use, limits = c(scale.min, scale.max), breaks = c(scale.min, scale.mid, scale.max))
        }
      } else {
        plot <- plot + scale_color_manual(values = scale.to.use)
      }
    }

    # modify different element of the plot
    plot <- plot + theme
    if (!plot.panel.grid) {plot <- plot + theme(panel.grid = element_blank())}
    if (!plot.axis.labels) {plot <- plot + theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())}

    # prepare title
    if (plot.title) {
      if (is.null(title)) {
        if (!is.null(exclusion.feature)) {
          exclusion.feature <- paste(exclusion.feature, collapse = "&")
        }
        if (!is.null(cutoff.feature)) {
          cutoff.feature <- paste(cutoff.feature, collapse = "&")
        }
        if (length(cutoff.expression) > 1) {
          cutoff.expression <- paste(cutoff.expression, collapse = "&")
        }

        args <- list(x = x, f = levels(as.factor(freqs$freq.expr)), y = cutoff.feature, z = cutoff.expression, q = exclusion.feature)
        title <-
          if (!is.null(cutoff.feature) && !is.null(exclusion.feature)) {
            if (plot.freq.of.expr.title && !is.null(freqs)) {
              if (make.italic) {
                substitute(paste(italic(x), " in ", italic(y), ">", z, " & ", italic(q), " = 0", " (", f, ")", sep = ""), args)
              } else {
                substitute(paste(x, " in ", italic(y), ">", z, " & ", italic(q), " = 0", " (", f, ")", sep = ""), args)
              }
            } else {
              if (make.italic) {
                substitute(paste(italic(x), " in ", italic(y), ">", z, " & ", italic(q), " = 0", sep = ""), args)
              } else {
                substitute(paste(x, " in ", italic(y), ">", z, " & ", italic(q), " = 0", sep = ""), args)
              }
            }
          } else if (!is.null(cutoff.feature)) {
            if (plot.freq.of.expr.title && !is.null(freqs)) {
              if (make.italic) {
                substitute(paste(italic(x), " in ", italic(y), ">", z, " (", f, ")", sep = ""), args)
              } else {
                substitute(paste(x, " in ", italic(y), ">", z, " (", f, ")", sep = ""), args)
              }
            } else {
              if (make.italic) {
                substitute(paste(italic(x), " in ", italic(y), ">", z, sep = ""), args)
              } else {
                substitute(paste(x, " in ", italic(y), ">", z, sep = ""), args)
              }
            }
          } else if (!is.null(exclusion.feature)) {
            if (plot.freq.of.expr.title && !is.null(freqs)) {
              if (make.italic) {
                substitute(paste(italic(x), " in ", italic(q), " = 0", " (", f, ")", sep = ""), args)
              } else {
                substitute(paste(x, " in ", italic(q), " = 0", " (", f, ")", sep = ""), args)
              }
            } else {
              if (make.italic) {
                substitute(paste(italic(x), " in ", italic(q), " = 0", sep = ""), args)
              } else {
                substitute(paste(x, " in ", italic(q), " = 0", sep = ""), args)
              }
            }
          } else {
            if (plot.freq.of.expr.title && !is.null(freqs)) {
              if (make.italic) {
                substitute(paste(italic(x), " (", f, ")", sep = ""), args)
              } else {
                substitute(paste(x, " (", f, ")", sep = ""), args)
              }
            } else {
              if (make.italic) {
                substitute(paste(italic(x), sep = ""), args)
              } else {
                substitute(paste(x, sep = ""), args)
              }
            }
          }
      }
    }


    if (make.italic) {
      legend.title <- substitute(paste(italic(x)), list(x = x))
    } else {
      legend.title <-substitute(paste(x), list(x = x))
    }

    # legend options
    if (length(legend.position) == 1) {
      plot <- plot + theme(legend.position = legend.position)
      if (legend.position %in% c("top", "bottom")) {
        temp <- legend.barheight
        legend.barheight <- legend.barwidth
        legend.barwidth <- temp
      }
    } else {
      plot <- plot + theme(legend.justification = c(legend.position[1], legend.position[2]), legend.position = c(legend.position[1], legend.position[2]))
    }

    plot <- plot + guides(shape = if (hide.shape.legend) {"none"} else {guide_legend(override.aes = list(size = shape.legend.size))},
                          colour = if (plot.colourbar) {
                            guide_colourbar(barwidth = legend.barwidth,
                                            barheight = legend.barheight,
                                            label.theme = element_text(size = legend.text.size, family = theme_get()[["text"]][["family"]]),
                                            title.theme = element_text(size = legend.title.text.size, family = theme_get()[["text"]][["family"]]),
                                            title = switch(plot.legend.title, legend.title, NULL))
                          } else {
                            guide_legend(override.aes = list(size = col.legend.size),
                                         nrow = legend.nrow,
                                         ncol = legend.ncol,
                                         label.theme = element_text(size = legend.text.size, family = theme_get()[["text"]][["family"]]),
                                         title.theme = element_text(size = legend.title.text.size, family = theme_get()[["text"]][["family"]]),
                                         title = switch(plot.legend.title, legend.title, NULL))
                          })
    plot <-
      plot +
      ggtitle(substitute(paste(x, sep = ""), list(x = title))) +
      theme(plot.title = element_text(size = title.font.size),
            strip.text.x = element_text(size = strip.font.size),
            legend.background = element_blank(),
            legend.key.size = unit(legend.key.size, "cm"),
            ...)
    '
            legend.key.width = unit(legend.barwidth, "cm"),
            legend.key.height = unit(legend.barheight, "cm"),
            legend.title = element_text(size = legend.title.text.size),
            legend.text = element_text(size = legend.text.size),
            ...)'

    # define facets and plot freq.of.expr annotation
    wrap_by <- function(...) {facet_wrap(vars(...), labeller = label_wrap_gen(multi_line = F), scales = "free")}
    if (is.null(SO.split) && !is.null(split.by)) {
      plot <- plot + wrap_by(split.by)
    } else if (!is.null(SO.split) && is.null(split.by)) {
      plot <- plot + wrap_by(SO.split)
    } else if (!is.null(SO.split) && !is.null(split.by)) {
      plot <- plot + wrap_by(SO.split, split.by)
    }

    # plot cutoff feature inset plot (https://stackoverflow.com/questions/5219671/it-is-possible-to-create-inset-graphs)
    if (!is.null(cutoff.feature) && plot.cutoff.feature.inset) {
      if (length(SO) > 1 || !is.null(split.by)) {
        print("The inset plot does not appear on every facet.")
      }
      if (length(legend.position) > 1) {
        print("You may set the legend.position to left, right, bottom or top to indicate it is valid for every facet.")
      }

      inset.data <- do.call(rbind, lapply(names(SO), function(y) {
        data.frame(cbind(as.matrix(t(GetAssayData(SO[[y]], slot = "data", assay = assay)[cutoff.feature,,drop = F]))))
      }))

      if (cutoff.expression == 0) {
        cutoff.expression.plot <- floor_any(min(inset.data[,cutoff.feature][which(inset.data[,cutoff.feature] > 0)]), 0.1)
      } else {
        cutoff.expression.plot <- cutoff.expression
      }
      inset <- ggplot(inset.data, aes(!!sym(cutoff.feature))) +
        geom_density(adjust = 1) +
        geom_vline(xintercept = cutoff.expression.plot, color = "red") +
        theme_bw() +
        theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.x = element_text(size = 9, face = "italic"), plot.margin = unit(c(0,0,0,0), "cm"), plot.background = element_rect(fill = "transparent")) +
        facet_wrap(vars(!!cutoff.feature))
      plot <- cowplot::ggdraw() + cowplot::draw_plot(plot) + cowplot::draw_plot(inset, x = inset.position[1], y = inset.position[2], width = inset.size[1], height = inset.size[2])
    }
    return (plot)
  })

  if (!all(is.na(plot.strip.selection))) {
    for (i in 1:length(plots)) {
      if (!i %in% plot.strip.selection) {
        plots[[i]] <- plots[[i]] + theme(strip.text.x = element_blank(), strip.background = element_blank())
      }
    }
  }

  if (combine) {plots <- cowplot::plot_grid(plotlist = plots, ncol = ncol.combine, nrow = nrow.combine, align = "hv", axis = "tblr")}
  if (length(plots) == 1 && !combine) {
    plots <- plots[[1]]
  }
  return(plots)
}

.check.aliases <- function(x, feature.aliases, data) {
  if (!is.null(feature.aliases) && x %in% names(feature.aliases)) {
    print(paste0(x, " changed to ", as.character(feature.aliases[which(names(feature.aliases) == x)])))
    names(data)[1] <- as.character(feature.aliases[which(names(feature.aliases) == x)])
    x <- names(data)[1]
  }
  return(list(x, data))
}


.check.features <- function(SO, features, rownames = T, meta.data = T) {
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
  feat.check <- function(SO, features, rownames, meta.data, ignore.case) {
    unlist(lapply(features, function(x) {
      Reduce(intersect, lapply(SO, function(y) {
        if (rownames & !meta.data) {
          search <- rownames(y)
        } else if (!rownames & meta.data) {
          search <- names(y@meta.data)
        } else if (rownames & meta.data) {
          search <- c(names(y@meta.data), rownames(y))
        }
        grep(paste0("^",x,"$"), search, value = T, ignore.case = ignore.case)
      }))
    }))
  }

  if (is.null(features)) {return(NULL)}
  if (!is.vector(features)) {stop("Features must be a vector of strings.")}
  if (!rownames && !meta.data) {stop("Dont be stupid, rownnames = F and meta.data = F ?!")}
  if (!is.list(SO)) {SO <- list(SO)}
  features <- unique(features)

  features.out <- feat.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = T)
  if (length(features.out) == 0) {stop("No features has not been found in every SO. No features left to plot.")}
  hits <- hit.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = T)
  if (any(hits > 1)) {
    print(paste0(paste(names(hits)[which(hits > 1)], collapse = ","), " found more than once in at least one SO when ignoring case. So, case is being considered."))
    features.out <- feat.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = F)
    if (length(features.out) == 0) {stop("No features has not been found in every SO. No features left to plot.")}
    hits <- hit.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = F)
    if (any(hits > 1)) {
      stop(paste0(paste(names(hits)[which(hits > 1)], collapse = ","), " still found more than once in at least one SO when not ignoring case. Please fix this (check rownames and meta.data)."))
    } else {
      ignore.case <- F
    }
  } else {
    ignore.case <- T
  }

  if (length(features.out) < length(features)) {print(paste0("Features not found in every SO: ", paste(setdiff(features, features.out), collapse = ",")))}
  if (length(features.out) > length(features)) {stop("More features after check.features than before?!")}
  return(features.out)
}
