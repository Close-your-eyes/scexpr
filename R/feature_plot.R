#' Title
#'
#' @param SO
#' @param features
#' @param assay
#' @param dims
#' @param cells
#' @param downsample
#' @param make.cells.unique
#' @param pt.size
#' @param pt.size.expr.factor
#' @param order for meta.col: remains T if var is continous but becomes F if var is integer (~probably discrete)
#' @param min.q.cutoff
#' @param max.q.cutoff
#' @param reduction
#' @param split.by
#' @param shape.by
#' @param combine
#' @param ncol.combine
#' @param nrow.combine
#' @param nrow.inner
#' @param ncol.inner
#' @param feature.aliases
#' @param title
#' @param title.font.size
#' @param cutoff.feature
#' @param cutoff.expression
#' @param exclusion.feature
#' @param binary
#' @param col.expresser
#' @param legend.position
#' @param legend.title.text.size
#' @param legend.text.size
#' @param legend.barheight
#' @param legend.barwidth
#' @param legend.nrow
#' @param legend.ncol
#' @param legend.key.size
#' @param legend.shape.size
#' @param legend.col.size
#' @param hide.shape.legend
#' @param font.family
#' @param col.pal.c
#' @param col.pal.d
#' @param col.excluded.cells
#' @param col.non.expresser
#' @param col.pal.rev
#' @param theme
#' @param plot.axis.labels
#' @param plot.panel.grid
#' @param plot.freq.title
#' @param plot.freq
#' @param plot.legend.title
#' @param plot.title
#' @param order.discrete
#' @param freq.position
#' @param freq.font.size
#' @param plot.cutoff.feature.inset
#' @param inset.position
#' @param inset.size
#' @param strip.font.size
#' @param strip.selection
#' @param ...
#' @param plot.labels
#' @param label.size
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
feature_plot <- function(SO,
                         features,
                         assay = c("RNA", "SCT"),
                         dims = c(1,2),
                         cells = NULL,
                         downsample = 1,
                         make.cells.unique = F,
                         pt.size = 1,
                         pt.size.expr.factor = 1,
                         order = T,
                         min.q.cutoff = 0,
                         max.q.cutoff = 1,
                         reduction = "tsne",
                         split.by = NULL,
                         shape.by = NULL,
                         combine = T,
                         ncol.combine = NULL,
                         nrow.combine = NULL,
                         nrow.inner = NULL,
                         ncol.inner = NULL,
                         feature.aliases = NULL,
                         binary = F,

                         title = NULL,
                         title.font.size = 14,
                         cutoff.feature = NULL,
                         cutoff.expression = 0,
                         exclusion.feature = NULL,

                         legend.position = "right",
                         legend.title.text.size = 14,
                         legend.text.size = 10,
                         legend.barheight = 3,
                         legend.barwidth = 0.5,
                         legend.nrow = NULL,
                         legend.ncol = NULL,
                         legend.key.size = 0.3,
                         legend.shape.size = 3,
                         legend.col.size = 3,
                         hide.shape.legend = F,

                         font.family = "sans",
                         col.pal.c = "spectral",
                         col.pal.d = "custom",
                         col.excluded.cells = "grey95",
                         col.non.expresser = "grey85",
                         col.expresser = "tomato2",
                         col.pal.rev = F,

                         theme = ggplot2::theme_bw(),
                         plot.axis.labels = F,
                         plot.panel.grid = F,

                         plot.freq.title = T,
                         plot.freq = T,
                         plot.legend.title = F,
                         plot.title = T,

                         order.discrete = NULL, # "^" as first element, "$" as last to indicate plotting first or last of respective factor levels, first plotted ones are buried ## or order.discrete = T to plot most in order of abundance
                         freq.position = c(-Inf, Inf),
                         freq.font.size = 4,

                         plot.cutoff.feature.inset = F,
                         inset.position = c(0.75,0.1),
                         inset.size = c(0.15,0.15),

                         strip.font.size = 14,
                         strip.selection = NA,

                         plot.labels = NULL,
                         label.size = 12,
                         ...) {

  # speed = F,
  # tidy eval syntax: https://rlang.r-lib.org/reference/nse-force.html https://ggplot2.tidyverse.org/reference/aes.html#quasiquotation
  # numeric but discrete data columns from meta.data - how to tell that it is not continuous

  if (missing(SO)) {stop("Seurat object list or feature vector is missing.")}
  if (length(features) == 1) {combine <- F}
  if (combine && (!is.null(ncol.combine) && !is.null(nrow.combine))) {stop("Please only select one, ncol.combine or nrow.combine. Leave the other NULL.")}
  if (!is.null(ncol.inner) && !is.null(nrow.inner)) {stop("Please only select one, ncol.inner or nrow.inner. Leave the other NULL.")}
  if (!is.null(legend.nrow) && !is.null(legend.ncol)) {stop("Please only select one, legend.nrow or legend.ncol. Leave the other NULL.")}
  if (length(dims) != 2 || class(dims) != "numeric") {stop("dims has to be a numeric vector of length 2, e.g. c(1,2).")}
  if (!is.null(plot.labels)) {
    plot.labels <- match.arg(plot.labels, c("text", "label"))
  }

  assay <- match.arg(assay, c("RNA", "SCT"))
  if (max.q.cutoff > 1) {
    message("max.q.cutoff and min.q.cutoff are divided by 100. Please provide values between 0 and 1.")
    max.q.cutoff <- max.q.cutoff/100
    min.q.cutoff <- min.q.cutoff/100
  }

  SO <- .check.SO(SO = SO, assay = assay, split.by = split.by,shape.by = shape.by)
  reduction <- .check.reduction(SO = SO, reduction = reduction, dims = dims)
  features <- .check.features(SO = SO, features = unique(features))
  cells <- .check.and.get.cells(SO = SO,
                                assay = assay,
                                cells = cells,
                                make.cells.unique = make.cells.unique,
                                cutoff.feature = cutoff.feature,
                                cutoff.expression = cutoff.expression,
                                exclusion.feature = exclusion.feature,
                                downsample = downsample,
                                make.cells.unique.warning = 1)

  if (length(SO) > 1) {
    SO.split <- "SO.split"
  } else {
    SO.split <- NULL
  }

  plots <- lapply(features, function(x) {

    data <- .get.data(SO = SO, assay = assay, cells = names(cells), split.by = split.by, shape.by = shape.by, reduction = names(reduction), feature = x, min.q.cutoff = min.q.cutoff, max.q.cutoff = max.q.cutoff, order = order)
    # necessary to make sym(shape.by) here, for !!shape.by to work; not possible within ggplot2::aes()
    # make it after .get.data
    shape.by <- tryCatch(rlang::sym(shape.by), error = function (e) NULL)

    # generate legend labels
    # shorten the cell selection
    if (is.numeric(data[,1])) {
      if (all(data[which(rownames(data) %in% names(which(cells == 1))),1] == 0)) {
        scale.max <- 0
        scale.min <- 0
        scale.mid <- 0
      } else {
        scale.max <- max(data[which(rownames(data) %in% names(which(cells == 1))),1])
        scale.min <- min(data[intersect(which(data[,1] != 0), which(rownames(data) %in% names(which(cells == 1)))),1]) # != 0 for module scores
        scale.mid <- scale.min + ((scale.max - scale.min) / 2)

        scale.max <- as.numeric(format(ceiling_any(scale.max, 0.1), nsmall = 1))
        scale.min <- as.numeric(format(floor_any(scale.min, 0.1), nsmall = 1))
        scale.mid <- as.numeric(format(round(scale.mid, 1), nsmall = 1))
      }
    }

    # select color scale
    if (is.numeric(data[,1])) {
      if (length(col.pal.c) == 1 && !col.pal.c %in% grDevices::colors()) {
        col.pal <- col_pal(name = col.pal.c, reverse = col.pal.rev)
      } else {
        col.pal <- col.pal.c
      }
    } else {
      if (length(col.pal.d) == 1 && !col.pal.d %in% grDevices::colors()) {
        col.pal <- col_pal(name = col.pal.d, reverse = col.pal.rev, n = nlevels(as.factor(data[,1])))
      } else {
        col.pal <- col.pal.d
      }
    }

   ' if (speed) {
      if (!requireNamespace("scattermore", quietly = T)) {
        devtools::install_github("exaexa/scattermore")
      }
      if (nlevels(as.factor(data$SO.split)) > 1 && nlevels(as.factor(data$split.by)) > 1) {
        subtitle <- "both"
      } else if (nlevels(as.factor(data$SO.split)) > 1) {
        subtitle <- "SO"
      } else if (nlevels(as.factor(data$split.by)) > 1) {
        subtitle <- "split"
      } else {
        subtitle <- "none"
      }
      plot <- lapply(levels(as.factor(data$SO.split)), function(y) {
        lapply(levels(as.factor(data$split.by)), function(z) {
          data <- data[intersect(which(data$SO.split == y), which(data$split.by == z)),]
          if (nrow(data) > 0) {
            plot <- ggplot2::ggplot() +  #no real data to plot
              suppressWarnings(scattermore::geom_scattermost(data[names(which(cells == 0)),c(paste0(reduction, "_", dims[1]), paste0(reduction, "_", dims[2]))],    #manually plot and color the points
                                                             col = col.excluded.cells,
                                                             pointsize = pt.size*2)) +
              suppressWarnings(scattermore::geom_scattermost(data[intersect(which(cells == 1), which(data[,1] == 0)),c(paste0(reduction, "_", dims[1]), paste0(reduction, "_", dims[2]))],    #manually plot and color the points
                                                             col = col.non.expresser,
                                                             pointsize = pt.size*2)) +
              scattermore::geom_scattermost(data[intersect(which(cells == 1), which(data[,1] > 0)),c(paste0(reduction, "_", dims[1]), paste0(reduction, "_", dims[2]))],    #manually plot and color the points
                                            col = col_pal("spectral")[1+99*(data[intersect(which(cells == 1), which(data[,1] > 0)),1]-min(data[intersect(which(cells == 1), which(data[,1] > 0)),1]))/(max(data[intersect(which(cells == 1), which(data[,1] > 0)),1])-min(data[intersect(which(cells == 1), which(data[,1] > 0)),1]))],
                                            pointsize = pt.size*2) +
              ggplot2::labs(subtitle = switch(subtitle, "both" = paste0(y, "  ", z),"SO" = y,"split" = z,"none" = ""), title = x)

            #geom_point(data=data.frame(x=double(0)), aes(x,x,color=x)) +
            #scale_color_gradientn(limits = c(min(data[intersect(which(cells == 1), which(data[,1] > 0)),1]), max(data[intersect(which(cells == 1), which(data[,1] > 0)),1])), colors = col_pal("spectral"), name=x) +
            plot <- plot + theme
            if (!plot.panel.grid) {plot <- plot + ggplot2::theme(panel.grid = ggplot2::element_blank())}
            if (!plot.axis.labels) {plot <- plot + ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank())}
            return(plot)
          } else {
            return(NULL)
          }
        })
      })
      plot <- purrr::flatten(plot)
      plot <- plot[!sapply(plot, is.null)]
      plot <- cowplot::plot_grid(plotlist = plot, ncol = ncol.combine, nrow = nrow.combine, align = "hv", axis = "tblr")
      return(plot)
    }'

    # excluded cells
    plot <-
      ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(paste0(reduction, "_", dims[1])), y = !!rlang::sym(paste0(reduction, "_", dims[2])))) +
      ggplot2::geom_point(data = data[names(which(cells == 0)),], ggplot2::aes(shape = !!shape.by), size = pt.size, color = col.excluded.cells)

    # different procedure for gene feature or meta.data feature
    if (x %in% rownames(SO[[1]])) {

      freqs <- .get.freqs(data = data, cells = cells, reduction = reduction, dims = dims, split.by = split.by, SO.split = SO.split)
      aliases.list <- .check.aliases(x, feature.aliases, data)
      x <- aliases.list[[1]]
      data <- aliases.list[[2]]

      # plot expressers and non-expressers
      if (binary) {
        if (any(data[,1] < 0)) {
          message("binary (+/-) may not be meaningful as there are negative expression values for ", feature, ".")
        }
        data[,1] <- ifelse(data[,1] > 0, "+", "-")
        plot <- plot + ggplot2::geom_point(data = data[which(cells == 1),], ggplot2::aes(shape = !!shape.by, color = !!rlang::sym(x)), size = pt.size*pt.size.expr.factor) + ggplot2::scale_color_manual(values = c(col.non.expresser, col.expresser))
      } else {
        # non-expressers
        plot <- plot + ggplot2::geom_point(data = data[intersect(which(cells == 1), which(data[,1] == 0)),], ggplot2::aes(shape = !!shape.by), size = pt.size, color = col.non.expresser)
        # expressers
        plot <- plot + ggplot2::geom_point(data = data[intersect(which(cells == 1), which(data[,1] > 0)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size*pt.size.expr.factor)
      }

      if (length(SO) > 1) {
        if (is.null(split.by)) {label = "freq.expr.by.SO"} else {label = "freq.expr.by.split.by.SO"}
      } else if (!is.null(split.by) && length(SO) == 1) {
        label = "freq.expr.by.split"
      } else {
        label = "freq.expr"
      }
      if (plot.freq) {
        plot <- plot + ggrepel::geom_text_repel(data = freqs, family = font.family, size = freq.font.size, ggplot2::aes(label = !!rlang::sym(label), x = xmin + abs(xmin - xmax) * freq.position[1], y = ymin + abs(ymin - ymax) * freq.position[2]))
      }

      plot.colorbar <- !binary
      make.italic <- T

    } else {

      freqs <- NULL
      aliases.list <- .check.aliases(x, feature.aliases, data)
      x <- aliases.list[[1]]
      data <- aliases.list[[2]]

      # shuffle to make plotting order random
      # only if order == TRUE and data is integer/discrete
      # or if is not numeric
      if ((!order && as.integer(data[,1]) == data[,1]) || !is.numeric(data[,1])) {
        data <- data[sample(x = 1:nrow(data), size = nrow(data), replace = F),]
      }

      if (is.null(order.discrete)) {
        plot <- plot + ggplot2::geom_point(data = data[names(which(cells == 1)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)
      } else if (is.logical(order.discrete) && order.discrete) {
        temp_name <- make.names(as.character(Sys.time()))
        data <-
          data %>%
          tibble::rownames_to_column(temp_name) %>%
          dplyr::group_by(!!rlang::sym(x)) %>%
          dplyr::mutate(count = dplyr::n()) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(-count) %>%
          tibble::column_to_rownames(temp_name)
        # use rownames %in% ... to preserve order
        plot <- plot + ggplot2::geom_point(data = data[rownames(data) %in% names(which(cells == 1)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)
      } else {
        if (length(unique(order.discrete)) != length(order.discrete)) {
          order.discrete <- unique(order.discrete)
          message("order.discrete is made unique: ", paste(order.discrete, collapse = ", "))
        }

        if (any(!gsub("\\$$", "", (gsub("^\\^", "", order.discrete))) %in% levels(as.factor(data[,1])))) {
          order.discrete <- order.discrete[which(gsub("\\$$", "", (gsub("^\\^", "", order.discrete))) %in% levels(as.factor(data[,1])))]
          message("order.discrete reduced to existing levels: ", paste(gsub("\\$$", "", (gsub("^\\^", "", order.discrete))), collapse = ", "))
        }

        # optionally: pt.size as names of order.discrete
        for (i in order.discrete[which(grepl("^\\^", order.discrete))]) {
          i <- gsub("^\\^", "", i)
          plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(data[,x] == i)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size) #[i]
        }

        for (i in order.discrete[intersect(which(!grepl("^\\^", order.discrete)), which(!grepl("\\$$", order.discrete)))]) {
          plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(data[,x] == i)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size) #[i]
        }
        plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(!data[,x] %in% gsub("\\$$", "", (gsub("^\\^", "", order.discrete))))),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)

        for (i in order.discrete[which(grepl("\\$$", order.discrete))]) {
          i <- gsub("\\$$", "", i)
          plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(data[,x] == i)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size) #[i]
        }

      }

      if (!is.null(plot.labels)) {
        if (is.numeric(data[,1])) {
          message("Labels not plotted as ", x, " is numeric.")
        } else {
          label_df <- do.call(rbind, lapply(unique(data[,1]), function(z) data.frame(label = z, avg1 = mean(data[which(data[,1] == z), paste0(reduction, "_", dims[1])]), avg2 = mean(data[which(data[,1] == z), paste0(reduction, "_", dims[2])]))))
          names(label_df)[c(2,3)] <- c(paste0(reduction, "_", dims[1]), paste0(reduction, "_", dims[2]))
          if (plot.labels == "text") {
            plot <- plot + ggplot2::geom_text(data = label_df, ggplot2::aes(label = label), size = label.size, family = font.family,) #...
          }
          if (plot.labels == "label") {
            plot <- plot + ggplot2::geom_label(data = label_df, ggplot2::aes(label = label), size = label.size, family = font.family) #...
          }
        }
      }

      plot.colorbar <- is.numeric(data[,1])
      make.italic <- F
    }

    if (!binary) {
      if (is.numeric(data[,1])) {
        if (min.q.cutoff > 0) {min.lab <- paste0(scale.min, " (q", round(min.q.cutoff*100, 0), ")")} else {min.lab <- scale.min}
        if (max.q.cutoff < 1) {max.lab <- paste0(scale.max, " (q", round(max.q.cutoff*100, 0), ")")} else {max.lab <- scale.max}
        if (length(unique(c(scale.min, scale.mid, scale.max))) > 1) {
          plot <- plot + ggplot2::scale_color_gradientn(colors = col.pal, limits = c(scale.min, scale.max), breaks = c(scale.min, scale.mid, scale.max), labels = c(min.lab,  scale.mid, max.lab))
        } else {
          # if no expressers are found: breaks and labels would be of different lengths
          plot <- plot + ggplot2::scale_color_gradientn(colors = col.pal, limits = c(scale.min, scale.max), breaks = c(scale.min, scale.mid, scale.max))
        }
      } else {
        plot <- plot + ggplot2::scale_color_manual(values = col.pal)
      }
    }

    if (plot.title && is.null(title)) {
      title <- .get.title(exclusion.feature = exclusion.feature, cutoff.feature = cutoff.feature, cutoff.expression = cutoff.expression, freqs = freqs, make.italic = make.italic, plot.freq.title = plot.freq.title, feature = x)
    }

    ##### theme and legend #####

    # modify different element of the plot
    plot <- plot + theme
    if (!plot.panel.grid) {plot <- plot + ggplot2::theme(panel.grid = ggplot2::element_blank())}
    if (!plot.axis.labels) {plot <- plot + ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank())}


    # legend options
    if (make.italic) {
      legend.title <- substitute(paste(italic(x)), list(x = x))
    } else {
      legend.title <- substitute(paste(x), list(x = x))
    }
    if (length(legend.position) == 1) {
      plot <- plot + ggplot2::theme(legend.position = legend.position)
      if (legend.position %in% c("top", "bottom")) {
        temp <- legend.barheight
        legend.barheight <- legend.barwidth
        legend.barwidth <- temp
      }
    } else {
      if (length(SO) > 1 || !is.null(split.by)) {
        message("You may set the legend.position to left, right, bottom or top to indicate it is valid for every facet.")
      }
      plot <- plot + ggplot2::theme(legend.justification = c(legend.position[1], legend.position[2]), legend.position = c(legend.position[1], legend.position[2]))
    }

    # theme_get()[["text"]][["family"]]
    plot <-
      plot +
      ggplot2::guides(
        shape = if (hide.shape.legend) {
          "none"
        } else {
          ggplot2::guide_legend(override.aes = list(size = legend.shape.size),
                                nrow = legend.nrow,
                                ncol = legend.ncol,
                                label.theme = ggplot2::element_text(size = legend.text.size, family = font.family),
                                title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family),
                                title = switch(plot.legend.title, legend.title, NULL))
        },
        color = if (plot.colorbar) {
          ggplot2::guide_colorbar(barwidth = legend.barwidth,
                                  barheight = legend.barheight,
                                  label.theme = ggplot2::element_text(size = legend.text.size, family = font.family),
                                  title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family),
                                  title = switch(plot.legend.title, legend.title, NULL))
        } else {
          ggplot2::guide_legend(override.aes = list(size = legend.col.size),
                                nrow = legend.nrow,
                                ncol = legend.ncol,
                                label.theme = ggplot2::element_text(size = legend.text.size, family = font.family),
                                title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family),
                                title = switch(plot.legend.title, legend.title, NULL))
        })

    plot <-
      plot +
      ggplot2::ggtitle(substitute(paste(x, sep = ""), list(x = title))) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title.font.size, family = font.family),
                     strip.text.x = ggplot2::element_text(size = strip.font.size, family = font.family),
                     legend.background = ggplot2::element_blank(),
                     legend.key.size = ggplot2::unit(legend.key.size, "cm"),
                     legend.key = ggplot2::element_blank(), ## keeps background of legend symbols transparent
                     ...)


    # define facets and plot freq.of.expr annotation

    wrap_by <- function(...) {ggplot2::facet_wrap(ggplot2::vars(...), labeller = ggplot2::label_wrap_gen(multi_line = F), scales = "free", nrow = nrow.inner, ncol = ncol.inner)}
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
        message("The inset plot does not appear on every facet.")
      }
      inset.data <- do.call(rbind, lapply(names(SO), function(y) {
        data.frame(cbind(as.matrix(t(Seurat::GetAssayData(SO[[y]], slot = "data", assay = assay)[cutoff.feature,,drop = F]))))
      }))
      if (cutoff.expression == 0) {
        cutoff.expression.plot <- floor_any(min(inset.data[,cutoff.feature][which(inset.data[,cutoff.feature] > 0)]), 0.1)
      } else {
        cutoff.expression.plot <- cutoff.expression
      }
      inset <- ggplot2::ggplot(inset.data, ggplot2::aes(!!rlang::sym(cutoff.feature))) +
        ggplot2::geom_density(adjust = 1) +
        ggplot2::geom_vline(xintercept = cutoff.expression.plot, color = "red") +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), strip.text.x = ggplot2::element_text(size = 9, face = "italic", family = font.family), plot.margin = ggplot2::unit(c(0,0,0,0), "cm"), plot.background = ggplot2::element_rect(fill = "transparent")) +
        ggplot2::facet_wrap(ggplot2::vars(!!cutoff.feature))
      plot <- cowplot::ggdraw() + cowplot::draw_plot(plot) + cowplot::draw_plot(inset, x = inset.position[1], y = inset.position[2], width = inset.size[1], height = inset.size[2])
    }
    return (plot)
  })


  if (!all(is.na(strip.selection))) {
    for (i in 1:length(plots)) {
      if (!i %in% strip.selection) {
        plots[[i]] <- plots[[i]] + ggplot2::theme(strip.text.x = ggplot2::element_blank(), strip.background = ggplot2::element_blank())
      }
    }
  }

  if (combine) {
    plots <- cowplot::plot_grid(plotlist = plots, ncol = ncol.combine, nrow = nrow.combine, align = "hv", axis = "tblr")
  }
  if (length(plots) == 1 && !combine) {
    plots <- plots[[1]]
  }
  return(plots)
}


.check.SO <- function(SO,
                      assay = c("RNA", "SCT"),
                      split.by = NULL,
                      shape.by = NULL,
                      length = NULL) {
  if (!is.list(SO)) {
    SO <- list(SO)
  }

  if (!is.null(length)) {
    if (length(SO) > length) {
      stop("Please provide exactly ", length, " SO.")
    }
  }

  assay <- match.arg(assay, c("RNA", "SCT"))

  if (is.null(names(SO)) && length(SO) > 1) {
    message("List of SO has no names. Naming them numerically in order as provided.")
    names(SO) <- as.character(seq_along(SO))
  }

  if (is.null(names(SO))) {
    names(SO) <- as.character(seq_along(SO))
  }

  if (!any(unlist(lapply(SO, class)) == "Seurat")) {
    stop("All SO have to Seurat objects (class == Seurat).")
  }

  if (!is.null(split.by) && length(.check.features(SO, split.by, rownames = F)) == 0) {
    stop("split.by not found in all objects.")
  }
  if (!is.null(shape.by) && length(.check.features(SO, shape.by, rownames = F)) == 0) {
    stop("shape.by not found in all objects.")
  }

  SO <- lapply(SO, function(x) {
    if (!assay %in% Seurat::Assays(x)) {
      stop("Assay not found in every SO.")
    }
    Seurat::DefaultAssay(x) <- assay
    return(x)
  })

  ## check if data has been scaled (compare to counts)
  check <- unlist(lapply(SO, function(x) identical(Seurat::GetAssayData(x, assay = assay, slot = "data"), Seurat::GetAssayData(x, assay = assay, slot = "counts"))))
  if (any(check)) {
    warning("data slot in at least one SO does not seem to contain normalized data since it is equal to the counts slot. You may want to normalize.")
  }

  if (!is.null(length) && length == 1) {
    return(SO[[1]])
  } else {
    return(SO)
  }
}

.check.reduction <- function(SO,
                             reduction,
                             dims = c(1,2)) {

  if (!is.list(SO)) {
    SO <- list(SO)
  }
  common_red <- Reduce(intersect, lapply(SO, function(x) names(x@reductions)))
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

.check.and.get.cells <- function (SO,
                                  assay = c("RNA", "SCT"),
                                  make.cells.unique = F,
                                  cells = NULL,
                                  cutoff.feature = NULL,
                                  cutoff.expression = NULL,
                                  exclusion.feature  = NULL,
                                  downsample = 1,
                                  make.cells.unique.warning = 1,
                                  return.included.cells.only = F) {

  if (!is.list(SO)) {
    SO <- list(SO)
  }

  if (!is.null(cutoff.feature) && length(.check.features(SO, cutoff.feature)) == 0) {
    stop("cutoff.feature not found in every SO.")
  }
  if (!is.null(exclusion.feature) && length(.check.features(SO, exclusion.feature)) == 0) {
    stop("exclusion.feature not found in every SO.")
  }

  if (length(cutoff.feature) != length(cutoff.expression) && length(cutoff.expression) != 1) {
    stop("cutoff.feature and cutoff.expression need to have the same length.")
  }
  if (length(cutoff.feature) != length(cutoff.expression) && length(cutoff.expression) == 1) {
    cutoff.expression <- rep(cutoff.expression, length(cutoff.feature))
  }

  assay <- match.arg(assay, c("RNA", "SCT"))

  if (downsample == 0) {
    stop("downsample cannot be 0.")
  }

  if (any(duplicated(cells))) {
    message("Duplicates found in cells.")
    cells <- unique(cells)
  }

  # check if cell names are unique across SOs
  all.cells <- unlist(lapply(SO, function(x) Seurat::Cells(x)))

  if (length(SO) > 1 && any(duplicated(all.cells)) && is.null(cells)) {
    # when cells is not provided it does not matter
    make.cells.unique <- T
  }

  if (length(SO) > 1 && any(duplicated(all.cells)) && !make.cells.unique) {
    if (make.cells.unique.warning == 1) {
      stop("Cell names are not unique across SOs. Please fix that manually with Seurat::RenameCells or pass make.cells.unique = T when calling this function. Cells are then renamed with the prefix SO_i_. Where i the index of the SO in the list. Consider this renaming of cells when selecting cells for plotting.")
    } else {
      stop("Cell names are not unique across SOs. Please fix that manually with Seurat::RenameCells.")
    }
  }

  if (length(SO) > 1 && any(duplicated(all.cells)) && make.cells.unique) {
    names.temp <- names(SO)
    SO <- lapply(seq_along(SO), function(x) Seurat::RenameCells(SO[[x]], add.cell.id = paste0("SO_", x)))
    names(SO) <- names.temp
    all.cells <- unlist(lapply(SO, function(x) Seurat::Cells(x)))
    #print("Cells have been prefixed with 'SO_i_' .")
    assign("SO", SO, envir = parent.frame()) # assigns in parent environment (https://stackoverflow.com/questions/10904124/global-and-local-variables-in-r?rq=1)
  }

  if (is.null(cells)) {
    cells <- all.cells
  } else {
    if (length(intersect(cells, all.cells)) == 0) {
      if (make.cells.unique) {
        stop("Non of cells found in SO. Cells have been prefixed with 'SO_i_' as cell names were not unique across SOs. Please consider when selecting a subset of cells for plotting. (e.g. manually prefix them with 'SO_i_' where i is the index of the SO.")
      } else {
        stop("Non of cells found in SO.")
      }
    }
    if (length(intersect(cells, all.cells)) < length(cells)) {
      message("Not all cells found in SO(s). Reduced to those which exist.")
    }
    cells <- intersect(cells, all.cells)
    if (length(cells) == 0) {
      stop("None of cells found.")
    }
  }

  # create a vector of cells which identifies how to plot them; 0 indicates exclusion
  cells.plot <- stats::setNames(rep(1, length(all.cells)), all.cells)
  ## cells excluded by arbitrary selection
  cells.plot[!names(cells.plot) %in% cells] <- 0
  ## cells excluded by
  if (!is.null(cutoff.feature)) {
    compare_fun <- function(x,y) {x > y}
    exclude.cells <- names(which(Matrix::colSums(sweep(Seurat::GetAssayData(x, slot = "data", assay = assay)[cutoff.feature,,drop=F], 1, cutoff.expression, compare_fun)) < length(cutoff.feature)))
    cells.plot[which(names(cells.plot) %in% exclude.cells)] <- 0
  }
  if (length(cells.plot) == 0) {
    stop("No cells left after filtering for cutoff.feature.")
  }
  ## cells excluded due to expression of an unwanted gene
  if (!is.null(exclusion.feature)) {
    exclude.cells <- unlist(lapply(SO, function(x) names(which(Matrix::colSums(Seurat::GetAssayData(x, slot = "data", assay = assay)[exclusion.feature,,drop=F]) > 0))))
    cells.plot[which(names(cells.plot) %in% exclude.cells)] <- 0
  }

  if (length(cells.plot) == 0) {
    stop("No cells left after filtering for exclusion.feature")
  }
  # downsample - completely remove cells to speed up plotting which may alter the statistics though
  if (downsample < 1) {
    downsample <- downsample*length(cells.plot)
  } else if (downsample > 1) {
    downsample <- min(downsample, length(cells.plot))
  }
  if (downsample == 0) {
    stop("downsample has become 0.")
  }
  if (downsample != 1) {
    cells.plot <- cells.plot[sample(seq_along(cells.plot), size = downsample)]
  }

  if (return.included.cells.only) {
    return(names(cells.plot[which(cells.plot == 1)]))
  } else if (!return.included.cells.only) {
    return(cells.plot)
  }

}


.check.aliases <- function(x, feature.aliases, data) {
  if (!is.null(feature.aliases) && x %in% names(feature.aliases)) {
    message(x, " changed to ", as.character(feature.aliases[which(names(feature.aliases) == x)]))
    names(data)[1] <- as.character(feature.aliases[which(names(feature.aliases) == x)])
    x <- names(data)[1]
  }
  return(list(x, data))
}


.check.features <- function(SO,
                            features,
                            rownames = T,
                            meta.data = T) {

  if (is.null(features)) {return(NULL)}
  if (!rownames && !meta.data) {return((NULL))}
  if (!is.list(SO)) {
    SO <- list(SO)
  }
  features <- unique(features)

  features.out <- .feat.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = T)
  if (length(features.out) == 0) {stop("Non of the provided features has not been found in every SO. No features left to plot.")}
  hits <- .hit.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = T)
  if (any(hits > 1)) {
    message(paste(names(hits)[which(hits > 1)], collapse = ","), " found more than once in at least one SO when ignoring case. So, case is being considered.")
    features.out <- .feat.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = F)
    if (length(features.out) == 0) {stop("Non of the provided features has not been found in every SO. No features left to plot.")}
    hits <- .hit.check(SO = SO, features = features, rownames = rownames, meta.data = meta.data, ignore.case = F)
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

.hit.check <- function(SO, features, rownames, meta.data, ignore.case) {
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
.feat.check <- function(SO, features, rownames, meta.data, ignore.case) {
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

.get.data <- function(SO,
                      feature,
                      assay = c("RNA", "SCT"),
                      slot = "data",
                      cells = NULL,
                      split.by = NULL,
                      shape.by = NULL,
                      reduction = "tsne",
                      meta.col = NULL,
                      min.q.cutoff = 0,
                      max.q.cutoff = 1,
                      order = T) {

  assay <- match.arg(assay, c("RNA", "SCT"))
  if (max.q.cutoff > 1) {
    message("max.q.cutoff and min.q.cutoff are divided by 100. Please provide values between 0 and 1.")
    max.q.cutoff <- max.q.cutoff/100
    min.q.cutoff <- min.q.cutoff/100
  }

  data <- do.call(rbind, lapply(names(SO), function(y) {

    if (feature %in% rownames(SO[[y]]) && !feature %in% names(SO[[y]]@meta.data)) {

      data <- data.frame(t(as.matrix(Seurat::GetAssayData(SO[[y]], slot = slot, assay = assay)[feature,,drop = F])), check.names = F)

    } else if (!feature %in% rownames(SO[[y]]) && feature %in% names(SO[[y]]@meta.data)) {

      data <- data.frame(SO[[y]]@meta.data[,feature,drop=F], stringsAsFactors = F, check.names = F)

    } else {
      stop("Feature found in meta.data and rownames of SO.")
    }
    if (!is.null(reduction)) {
      data <- cbind(data, Seurat::Embeddings(SO[[y]], reduction = names(SO[[y]]@reductions)[which.min(utils::adist(reduction, names(SO[[y]]@reductions), ignore.case = T))]))
      #data <- cbind(data, Seurat::Embeddings(SO[[y]], reduction = names(SO[[y]]@reductions)[grepl(reduction, names(SO[[y]]@reductions), ignore.case = T)]))
    }
    if (!is.null(meta.col)) {
      data <- cbind(data, SO[[y]]@meta.data[,meta.col,drop=F])
    }

    if (is.null(split.by)) {
      data[,"split.by"] <- "1"
    } else {
      data[,"split.by"] <- SO[[y]]@meta.data[,split.by]
    }
    if (!is.null(shape.by)) {
      data[,shape.by] <- as.factor(as.character(SO[[y]]@meta.data[,as.character(shape.by)]))
    }

    data[,"SO.split"] <- y
    return(data)
  }))

  if (!is.null(cells)) {
    data <- data[cells,]
  }

  if (is.numeric(data[,1]) && all(data[,1] == 0)) {
    message("No expressers found for ", feature, ".")
  }

  # use squishing to dampen extreme values - this will produce actually wrong limits on the legend
  if (is.numeric(data[,1]) && (min.q.cutoff > 0 || max.q.cutoff < 1)) {
    if (all(data[,1] > 0)) {
      # expression is always greater than 0 and non-expresser are excluded
      data[,1][which(data[,1] > 0)] <- scales::squish(data[,1][which(data[,1] > 0)], range = c(stats::quantile(data[,1][which(data[,1] > 0)], min.q.cutoff), stats::quantile(data[,1][which(data[,1] > 0)], max.q.cutoff)))
    } else {
      # e.g. for module scores below 0
      data[,1] <- scales::squish(data[,1], range = c(stats::quantile(data[,1], min.q.cutoff), stats::quantile(data[,1], max.q.cutoff)))
    }
  }


  if (order) {
    # abs: in case negative values are contained in meta.col, any extreme away from 0 will be plotted on top
    if (is.numeric(data[,1])) {
      data <- data[order(abs(data[,1])),]
    } else {
      data <- data[order(data[,1]),]
    }

  }

  return(data)
}

.get.freqs <- function(data,
                       cells,
                       reduction,
                       dims,
                       split.by,
                       SO.split) {

  freqs <- do.call(rbind, lapply(split(unique(data[,c("split.by", "SO.split")]), seq(nrow(unique(data[,c("split.by", "SO.split")])))), function(r) {

    b <- r[,"split.by"]
    c <- r[,"SO.split"]

    ref.rows <- Reduce(intersect, list(which(rownames(data) %in% names(cells[which(cells == 1)])), which(data$split.by == b), which(data$SO.split == c)))
    freq.expr.by.split.by.SO <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    ref.rows <- Reduce(intersect, list(which(rownames(data) %in% names(cells[which(cells == 1)])), which(data$split.by == b)))
    freq.expr.by.split <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    ref.rows <- Reduce(intersect, list(which(rownames(data) %in% names(cells[which(cells == 1)])), which(data$SO.split == c)))
    freq.expr.by.SO <- sum(data[ref.rows,1] > 0) / length(data[ref.rows,1])*100

    ref.rows <- which(rownames(data) %in% names(cells[which(cells == 1)]))
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
    xrng <- range(data[ref.rows,paste0(reduction, "_", dims[1])])
    yrng <- range(data[ref.rows,paste0(reduction, "_", dims[2])])

    return(data.frame(split.by = b, SO.split = c, freq.expr.by.split.by.SO = freqs[1], freq.expr.by.split = freqs[2], freq.expr.by.SO = freqs[3], freq.expr = freqs[4], xmin = xrng[1], xmax = xrng[2], ymin = yrng[1], ymax = yrng[2]))
  }))

  return(freqs)
}

.get.title <- function(exclusion.feature,
                       cutoff.feature,
                       cutoff.expression,
                       freqs,
                       make.italic,
                       plot.freq.title,
                       feature) {
  if (!is.null(exclusion.feature)) {
    exclusion.feature <- paste(exclusion.feature, collapse = "&")
  }
  if (!is.null(cutoff.feature)) {
    cutoff.feature <- paste(cutoff.feature, collapse = "&")
  }
  if (length(cutoff.expression) > 1) {
    cutoff.expression <- paste(cutoff.expression, collapse = "&")
  }
  # italic bold==
  args <- list(x = feature, f = levels(as.factor(freqs$freq.expr)), y = cutoff.feature, z = cutoff.expression, q = exclusion.feature)
  title <-
    if (!is.null(cutoff.feature) && !is.null(exclusion.feature)) {
      if (plot.freq.title && !is.null(freqs)) {
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
      if (plot.freq.title && !is.null(freqs)) {
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
      if (plot.freq.title && !is.null(freqs)) {
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
      if (plot.freq.title && !is.null(freqs)) {
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
  return(title)
}

.get.contour.level <- function(kk, x,y,prob) {
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  stats::approx(c1, sz, xout = 1 - prob)$y
}

.get.contours <- function(x, cells = NULL, levels = 0.9) {

  if (!requireNamespace("reshape2", quietly = T)) {
    utils::install.packages("reshape2")
  }

  ## x only with x any y columm (dim reduction dims)
  # handle levels attribute differently

  ## handle special cells arguments before passing here in a separate function??
  # cells = all (outlining all cells)
  # cells = clusters (outlining each cluster individually)
  # cells = quick selection of clusters
  # cells = list of cell vectors with attr like color, linetype size, level

  # return a list of contours with levels as attr to pass to geom_contour
  # also pass on all attr from cells vector for plotting

  # optionally subset here
  if (!is.null(cells)) {
    x <- x[cells,]
  }

  kk <- MASS::kde2d(x[,1], x[,2])
  dimnames(kk$z) <- list(kk$x, kk$y)
  kk <- reshape2::melt(kk$z)
  names(kk) <- c(names(x)[1], names(x)[2], "value")
  attributes(kk) <- attributes(x)
  attr(kk, "levels") <- sapply(levels, function(z) .get.contour.level(kk = kk, x = x[,1], y = x[,2], prob = z))
  return(kk)
}

.prep.contour.cells <- function(SO, reduction, x) {
  #multiple SO - how to?

  x<-"all"
  attr(x, "levels") <- 0.9
  attr(x, "linetype") <- "dashed"



  if (length(x) == 1 && x == "all") {
    y <- Seurat::Embeddings(SO, reduction = names(reduction))
    attributes(y) <- c(attributes(y), attributes(x))
  }
  data <- cbind(Seurat::FetchData(SO, "SCT_snn_res.0.4"), Seurat::Embeddings(SO, reduction = "tsne"))

}

ceiling_any = function(x, accuracy, f = ceiling){f(x/ accuracy) * accuracy}

floor_any = function(x, accuracy, f = floor){f(x/ accuracy) * accuracy}
