#' Plot features of single cell transcriptomes on a dimension reduction map
#'
#' @param SO one or more Seurat object(s)
#' @param features vector of features to plot (genes or column names in meta data)
#' @param assay which assay to get expression data from
#' @param dims which dimensions of the selected dimension reduction to plot
#' @param cells a vector of cell names to include (not selected ones are plotted with color col.excluded.cells)
#' @param downsample downsample the number of cells (intended to speed up test plottings)
#' @param pt.size dot size per cells
#' @param pt.size.expr.factor factor of increased dot size for expressing cells
#' @param order for meta.col: remains T if var is continuous but becomes F if var is integer (~probably discrete)
#' @param order.abs do use absolute values for ordering (any value away from zero (+/) is treated equally)
#' @param shuffle do shuffle if order if FALSE; allows to define a definite order for plotting if set to to F
#' @param order.rev reverse the ordering to have lowest values on top (or zeros if order.abs = T)
#' @param min.q.cutoff decimal number (> 1, < 0) of lower quantile limit where to cut the color scale, intended to squish extremes
#' @param max.q.cutoff decimal number (> 1, < 0) of upper quantile limit where to cut the color scale, intended to squish extremes
#' @param reduction which reduction to use for plotting
#' @param split.by column in meta data to use to split plots
#' @param shape.by column in meta data to shape dots (cells) b<
#' @param combine combine multiple features to one plot (TRUE) or return a list with one entry per feature (FALSE)
#' @param ncol.combine number of columns in combined graphic (feature combined)
#' @param nrow.combine number of rows in combined graphic (feature combined)
#' @param nrow.inner number of plots per row within one feature (originating from multiple SO and/or split.by)
#' @param ncol.inner number of plots per column within one feature (originating from multiple SO and/or split.by)
#' @param feature.aliases vector aliases for features; e.g. c("FAIM3" = "FCMR", "MS4A1" = "CD20")
#' @param title force title for the plot instead of feature name which is the default
#' @param title.font.size font size of the title
#' @param cutoff.feature select a feature for cells to consider in plotting
#' @param cutoff.expression select the cutoff feature expression levels for filtering
#' @param exclusion.feature
#' @param binary
#' @param col.expresser
#' @param legend.position # values between 0 and -1
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
#' @param split.by.scales
#' @param na.rm
#' @param inf.rm
#' @param bury_NA
#' @param trajectory.slot name of Misc slot that contains information on how to plot a trajectory
#' @param trajectory.color
#' @param trajectory.size
#' @param trajectory.linetype
#' @param contour_feature
#' @param col.pal.contour
#' @param contour_args
#' @param plot.expr.freq.by.contour.group
#' @param use_ggnewscale_for_contour_colors
#' @param expand_limits
#' @param label.feature provide a column in meta.data of SO which should be used for labels
#' (e.g. a short form that does not disturb the plot so much)
#' @param legend.label.position
#' @param label.center.fun
#' @param label.nudge
#' @param na.value
#' @param contour.label.nudge
#' @param color.scale.labels
#' @param n.colorsteps number of steps (numeric) to divide color scale into; if null then ordinary continuous fill scale is chosen;
#' if of length 1 this is passed as n.breaks to scale_fill_stepsn, if length > 1 then passed as breaks to scale_fill_stepsn
#' @param nice.breaks passed to scale_color_stepsn if length(n.colorsteps) == 1
#' @param show.limits passed to scale_color_stepsn if length(n.colorsteps) > 1; show min and max limit on legend
#' @param dotplot
#' @param legend.decimals passed to scale_color_stepsn if length(n.colorsteps) > 1; number of decimals to round legend labels to
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
                         pt.size = 1,
                         pt.size.expr.factor = 1,
                         order = T,
                         order.abs = T,
                         shuffle = T,
                         order.rev = F,
                         min.q.cutoff = 0,
                         max.q.cutoff = 1,
                         reduction = "umap",
                         split.by = NULL,
                         split.by.scales = "fixed",
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
                         legend.label.position = "right",
                         hide.shape.legend = F,

                         font.family = "sans",
                         col.pal.c = "spectral",
                         col.pal.d = "custom",
                         col.excluded.cells = "grey95",
                         col.non.expresser = "grey85",
                         col.expresser = "tomato2",
                         col.pal.rev = F,
                         n.colorsteps = NULL,
                         nice.breaks = F,
                         show.limits = T,
                         legend.decimals = 1,

                         theme = ggplot2::theme_bw(),
                         plot.axis.labels = F,
                         plot.panel.grid = F,

                         plot.freq.title = NULL, # if length(SO) == 1 --> defaults to T, else F
                         plot.freq = NULL,  # if length(SO) > 1 --> defaults to T, else F
                         plot.legend.title = F,
                         plot.title = T,

                         order.discrete = T, # "^" as first element, "$" as last to indicate plotting first or last of respective factor levels, first plotted ones are buried ## or order.discrete = T to plot most in order of abundance
                         freq.position = c(-Inf, Inf),
                         freq.font.size = 4,

                         plot.cutoff.feature.inset = F,
                         inset.position = c(0.75,0.1),
                         inset.size = c(0.15,0.15),

                         strip.font.size = 14,
                         strip.selection = NA,

                         plot.labels = NULL,
                         label.feature = NULL,
                         label.size = 12,
                         label.center.fun = c("median", "mean"),
                         label.nudge = c(0,0), # nudging in x and y direction
                         na.rm = F,
                         inf.rm = F,
                         bury_NA = T,
                         na.value = "grey50",

                         trajectory.slot = NULL,
                         trajectory.color = "grey30",
                         trajectory.size = 0.75,
                         trajectory.linetype = "solid",

                         contour_feature = NULL,
                         col.pal.contour = "custom",
                         contour_args = list(contour_var = "ndensity", breaks = 0.3, linewidth = 1), # arguments to geom_density_2d
                         contour.label.nudge = c(0,0),
                         plot.expr.freq.by.contour.group = F,
                         use_ggnewscale_for_contour_colors = T, ## ggnewscale breaks the legend of dot colors; setting to F will avoid that but also does not allow to have a legend for contour lines
                         expand_limits = list(), # arguments to ggplot2::expand_limits
                         color.scale.labels = NULL,
                         ...) {


  ## add axis arrows, shortened:
  #p_blood <- p_blood + guides(x = ggh4x::guide_axis_truncated(trunc_lower = ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1], trunc_upper = ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1] + abs(min(c(ggplot_build(p)$layout$panel_params[[1]]$x.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1])) - max(c(ggplot_build(p)$layout$panel_params[[1]]$x.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1])))/4 ),
  #                          y = ggh4x::guide_axis_truncated(trunc_lower = ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1], trunc_upper = ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1] + abs(min(c(ggplot_build(p)$layout$panel_params[[1]]$y.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1])) - max(c(ggplot_build(p)$layout$panel_params[[1]]$y.range[2], ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1])))/4 ))


  ## label position calculation is not facetted!!

  # tidy eval syntax: https://rlang.r-lib.org/reference/nse-force.html https://ggplot2.tidyverse.org/reference/aes.html#quasiquotation
  # numeric but discrete data columns from meta.data - how to tell that it is not continuous

  # no axis extension: # https://stackoverflow.com/questions/48255449/how-can-i-set-exactly-the-limits-for-axes-in-ggplot2-r-plots

  if (missing(SO)) {stop("Seurat object list or feature vector is missing.")}
  if (length(features) == 1) {combine <- F}
  if (!is.null(label.feature) && length(label.feature) != 1) {stop("label.feature must have length 1.")}
  if (combine && (!is.null(ncol.combine) && !is.null(nrow.combine))) {stop("Please only select one, ncol.combine or nrow.combine. Leave the other NULL.")}
  if (!is.null(ncol.inner) && !is.null(nrow.inner)) {stop("Please only select one, ncol.inner or nrow.inner. Leave the other NULL.")}
  if (!is.null(legend.nrow) && !is.null(legend.ncol)) {stop("Please only select one, legend.nrow or legend.ncol. Leave the other NULL.")}
  if (length(dims) != 2 || !methods::is(dims, "numeric")) {stop("dims has to be a numeric vector of length 2, e.g. c(1,2).")}
  if (!is.null(plot.labels)) {plot.labels <- match.arg(plot.labels, c("text", "label"))}
  if ((length(order.discrete) %in% c(0,1)) && (is.null(order.discrete) || is.na(order.discrete))) {stop("order.discrete should be logical or a vector of factor levels in order.")}
  if (length(contour_feature) > 1) {stop("Only provide one contour_feature.")}
  if (length(legend.position) > 2) {stop("legend.position should have length 1 being top, bottom, left, right; or length 2 indicating the corner where legend is to place.")}
  if (list_depth(contour_args) > 1 && use_ggnewscale_for_contour_colors) {
    stop("Different setting for contour_args cannot be passed when use_ggnewscale_for_contour_colors=T. Set it to FALSE or only pass one setting for contour_args.")
  }

  label.center.fun <- match.arg(label.center.fun, c("median", "mean"))
  assay <- match.arg(assay, c("RNA", "SCT"))
  if (max.q.cutoff > 1) {
    #message("max.q.cutoff and min.q.cutoff are divided by 100. Please provide values between 0 and 1.")
    max.q.cutoff <- max.q.cutoff/100
    min.q.cutoff <- min.q.cutoff/100
  }

  if (is.null(plot.freq.title)) {
    plot.freq.title <- length(SO) == 1
  } else {
    if (!is.logical(plot.freq.title)) {
      stop("plot.freq.title has to be NULL, TRUE or FALSE.")
    }
  }

  if (is.null(plot.freq)) {
    plot.freq <- length(SO) > 1
  } else {
    if (!is.logical(plot.freq)) {
      stop("plot.freq has to be NULL, TRUE or FALSE.")
    }
  }

  SO <- .check.SO(SO = SO, assay = assay, split.by = split.by, shape.by = shape.by)
  reduction <- .check.reduction(SO = SO, reduction = reduction, dims = dims)
  features <- .check.features(SO = SO, features = unique(features))
  label.feature <- .check.features(SO = SO, features = label.feature, rownames = F)
  contour_feature <- .check.features(SO = SO, unique(contour_feature), rownames = F)
  cells <- .check.and.get.cells(SO = SO,
                                assay = assay,
                                cells = cells,
                                #make.cells.unique = make.cells.unique,
                                cutoff.feature = cutoff.feature,
                                cutoff.expression = cutoff.expression,
                                exclusion.feature = exclusion.feature,
                                downsample = downsample)
  #make.cells.unique.warning = 1)

  if (length(SO) > 1) {
    SO.split <- "SO.split"
  } else {
    SO.split <- NULL
  }

  plots <- lapply(features, function(x) {

    data <- .get.data(SO = SO,
                      assay = assay,
                      cells = names(cells),
                      split.by = split.by,
                      shape.by = shape.by,
                      reduction = names(reduction),
                      feature = x,
                      label.feature = label.feature,
                      min.q.cutoff = min.q.cutoff,
                      max.q.cutoff = max.q.cutoff,
                      order = order,
                      order.rev = order.rev,
                      order.abs = order.abs,
                      shuffle = shuffle,
                      order.discrete = order.discrete,
                      bury_NA = bury_NA,
                      na.rm = na.rm,
                      inf.rm = inf.rm,
                      trajectory.slot = trajectory.slot)

    if (!is.data.frame(data)) {
      ## when !is.null(trajectory.slot) and the slot has been found
      data_traj <- data[["data_traj"]]
      data <- data[["data"]]
      plot_traj <- T
    } else {
      plot_traj <- F
    }

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
        scale.max <- max(data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(is.finite(data[,1]))), 1], na.rm = T)
        scale.min <- min(data[Reduce(intersect, list(which(data[,1] != 0), which(rownames(data) %in% names(which(cells == 1))), which(is.finite(data[,1])))),1], na.rm = T) # != 0 for module scores
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
        plot <- plot + ggplot2::geom_point(data = data[intersect(rownames(data[which(data[,1] == 0),]), names(which(cells == 1))),], ggplot2::aes(shape = !!shape.by), size = pt.size, color = col.non.expresser)
        # expressers
        # leave intersect like that to maintain order of data
        plot <- plot + ggplot2::geom_point(data = data[intersect(rownames(data[which(data[,1] > 0),]), names(which(cells == 1))),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size*pt.size.expr.factor)
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

      if (is.logical(order.discrete)) {
        # order.discrete T or F: ordering has been done in .get.data
        # use rownames(data) %in% ... to preserve random order

        plot <- plot + ggplot2::geom_point(data = data[rownames(data) %in% names(which(cells == 1)),],
                                           ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by),
                                           size = pt.size)
      } else {
        if (length(unique(order.discrete)) != length(order.discrete)) {
          order.discrete <- unique(order.discrete)
          message("order.discrete is made unique: ", paste(order.discrete, collapse = ", "))
        }

        if (any(!gsub("\\$$", "", (gsub("^\\^", "", order.discrete))) %in% levels(as.factor(data[,1])))) {
          order.discrete <- order.discrete[which(gsub("\\$$", "", (gsub("^\\^", "", order.discrete))) %in% levels(as.factor(data[,1])))]
          message("order.discrete reduced to existing levels: ", paste(gsub("\\$$", "", (gsub("^\\^", "", order.discrete))), collapse = ", "))
        }

        for (i in order.discrete[which(grepl("^\\^", order.discrete))]) {
          i <- gsub("^\\^", "", i)
          plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(data[,x] == i)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)
        }

        for (i in order.discrete[intersect(which(!grepl("^\\^", order.discrete)), which(!grepl("\\$$", order.discrete)))]) {
          plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(data[,x] == i)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)
        }
        plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(!data[,x] %in% gsub("\\$$", "", (gsub("^\\^", "", order.discrete))))),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)

        for (i in order.discrete[which(grepl("\\$$", order.discrete))]) {
          i <- gsub("\\$$", "", i)
          plot <- plot + ggplot2::geom_point(data = data[intersect(which(rownames(data) %in% names(which(cells == 1))), which(data[,x] == i)),], ggplot2::aes(color = !!rlang::sym(x), shape = !!shape.by), size = pt.size)
        }
      }

      # put this below contour plotting so that labels are on top (see below)
      '      if (!is.null(plot.labels)) {
        if (is.numeric(data[,1])) {
          message("Labels not plotted as ", x, " is numeric.")
        } else {
          # potentially: use mclust, densityMclust(), (mixed gaussian model) to find multimodal clusters; x,y separately
          label_df <- do.call(rbind, lapply(unique(data[,1]), function(z) data.frame(label = z,
                                                                                     avg1 = mean(data[which(data[,1] == z), paste0(reduction, "_", dims[1])]),
                                                                                     avg2 = mean(data[which(data[,1] == z), paste0(reduction, "_", dims[2])]))))

          label_column <- if (!is.null(label.feature)) {
            temp <- unique(data[,c(1,which(colnames(data) == "label.feature")[1])])
            label_df[,1] <- as.character(label_df[,1])
            label_df[,1] <- stats::setNames(temp[,2,drop=T], temp[,1,drop=T])[label_df[,1]]
          }

          names(label_df)[c(2,3)] <- c(paste0(reduction, "_", dims[1]), paste0(reduction, "_", dims[2]))
          if (plot.labels == "text") {
            plot <- plot + ggplot2::geom_text(data = label_df, ggplot2::aes(label = label), size = label.size, family = font.family) #...
          }
          if (plot.labels == "label") {
            plot <- plot + ggplot2::geom_label(data = label_df, ggplot2::aes(label = label), size = label.size, family = font.family) #...
          }
        }
      }'

      plot.colorbar <- is.numeric(data[,1])
      make.italic <- F
    }

    if (!binary) {
      if (is.numeric(data[,1])) {
        if (min.q.cutoff > 0) {min.lab <- paste0(scale.min, " (q", round(min.q.cutoff*100, 0), ")")} else {min.lab <- scale.min}
        if (max.q.cutoff < 1) {max.lab <- paste0(scale.max, " (q", round(max.q.cutoff*100, 0), ")")} else {max.lab <- scale.max}

        if (length(unique(c(scale.min, scale.mid, scale.max))) > 1) {

          if (!is.null(n.colorsteps)) {
            if (length(n.colorsteps) == 1) {
              plot <-
                plot +
                ggplot2::scale_color_stepsn(colors = col.pal,
                                            n.breaks = n.colorsteps,
                                            nice.breaks = nice.breaks,
                                            na.value = na.value)
            } else {
              plot <-
                plot +
                ggplot2::scale_color_stepsn(colors = col.pal,
                                            breaks = n.colorsteps,
                                            show.limits = T,
                                            na.value = na.value,
                                            labels = function(x) round(x,legend.decimals))
            }
          } else {
            plot <- plot + ggplot2::scale_color_gradientn(colors = col.pal,
                                                          limits = c(scale.min, scale.max),
                                                          breaks = c(scale.min, scale.mid, scale.max),
                                                          labels = if (is.null(color.scale.labels)) {c(min.lab, scale.mid, max.lab)} else {color.scale.labels},
                                                          na.value = na.value)
          }
        } else {
          # if no expressers are found: breaks and labels would be of different lengths
          plot <- plot + ggplot2::scale_color_gradientn(colors = col.pal,
                                                        limits = c(scale.min, scale.max),
                                                        breaks = c(scale.min, scale.mid, scale.max),
                                                        na.value = na.value)
        }

      } else {
        plot <- plot + ggplot2::scale_color_manual(values = col.pal,
                                                   na.value = na.value)
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
      if (any(legend.position < 0) | any(legend.position > 1)) {
        message("legend.position should have values between 0 and 1 indicating left/right and bottom/top corners")
      }
      plot <- plot + ggplot2::theme(legend.justification = c(legend.position[1], legend.position[2]), legend.position = c(legend.position[1], legend.position[2]))
    }


    if (plot_traj) {
      plot <- plot + ggplot2::geom_segment(data = data_traj, ggplot2::aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                                           linewidth = trajectory.size,
                                           color = trajectory.color,
                                           linetype = trajectory.linetype,
                                           na.rm = TRUE)
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
                                label.theme = ggtext::element_markdown(size = legend.text.size, family = font.family),
                                title.theme = ggtext::element_markdown(size = legend.title.text.size, family = font.family),
                                label.position = legend.label.position,
                                title = switch(plot.legend.title, legend.title, NULL))
        },
        color = if (plot.colorbar) {
          ggplot2::guide_colorbar(barwidth = legend.barwidth,
                                  barheight = legend.barheight,
                                  label.theme = ggtext::element_markdown(size = legend.text.size, family = font.family),
                                  title.theme = ggtext::element_markdown(size = legend.title.text.size, family = font.family),
                                  label.position = legend.label.position,
                                  title = switch(plot.legend.title, legend.title, NULL))
        } else {
          ggplot2::guide_legend(override.aes = list(size = legend.col.size),
                                nrow = legend.nrow,
                                ncol = legend.ncol,
                                label.theme = ggtext::element_markdown(size = legend.text.size, family = font.family),
                                title.theme = ggtext::element_markdown(size = legend.title.text.size, family = font.family),
                                label.position = legend.label.position,
                                title = switch(plot.legend.title, legend.title, NULL))
        })

    #ggplot2::element_text
    #ggtext::element_markdown ## this allows to pass italics in the legend text with e.g. *CCR7* or <i>CCR7</i>

    plot <-
      plot +
      ggplot2::ggtitle(substitute(paste(x, sep = ""), list(x = title))) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title.font.size, family = font.family),
                     strip.text.x = ggplot2::element_text(size = strip.font.size, family = font.family),
                     legend.background = ggplot2::element_blank(),
                     legend.key.size = ggplot2::unit(legend.key.size, "cm"),
                     legend.key = ggplot2::element_blank(), ## keeps background of legend symbols transparent
                     ...)


    # plot contours, optionally
    if (!is.null(contour_feature)) {
      ## cells currently not considered
      contour_data <- .get.data(SO, feature = contour_feature, reduction = names(reduction))
      if (length(col.pal.contour) == 1 && !col.pal.contour %in% grDevices::colors()) {
        col.pal.contour <- col_pal(name = col.pal.contour, reverse = col.pal.rev, n = nlevels(as.factor(contour_data[,1])))
      }

      # order of col.pal.contour in case it has names??

      # more than one list to contour_args? (different settings for different contour for factor levels of contour_feature)
      # only when use_ggnewscale_for_contour_colors is FALSE
      if (list_depth(contour_args) == 1 && !use_ggnewscale_for_contour_colors) {
        contour_args <- rep(list(contour_args), length(unique(contour_data[,contour_feature])))

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
        names(contour_args) <- unique(contour_data[,contour_feature])

      } else if (length(lengths(contour_args)) != length(unique(contour_data[,contour_feature])) && !use_ggnewscale_for_contour_colors) {
        stop("Length of list of contour_args lists does not match length of factor levels of contour_feature.")
      } else {
        if (is.null(names(contour_args))) {
          stop("list of contour_args has to have names of factor levels of contour_feature.")
        }
      }

      if (use_ggnewscale_for_contour_colors) {
        plot <-
          plot +
          ggnewscale::new_scale_color() +
          Gmisc::fastDoCall(ggplot2::geom_density2d, args = c(contour_args, list(data = contour_data, mapping = ggplot2::aes(color = !!rlang::sym(contour_feature))))) +
          ggplot2::scale_color_manual(values = col.pal.contour)
      } else {
        for (i in 1:length(unique(contour_data[,contour_feature]))) {
          plot <-
            plot +
            Gmisc::fastDoCall(ggplot2::geom_density2d, args = c(contour_args[[unique(contour_data[,contour_feature])[i]]], list(data = contour_data[which(contour_data[,contour_feature] == unique(contour_data[,contour_feature])[i]),], color = col.pal.contour[i])))
        }
      }



      if (plot.expr.freq.by.contour.group) {
        group_labels <-
          dplyr::left_join(data %>% dplyr::mutate(ID = rownames(.)), contour_data %>% dplyr::mutate(ID = rownames(.)) %>% dplyr::select(dplyr::all_of(c("ID", contour_feature))), by = "ID") %>%
          dplyr::group_by(!!rlang::sym(contour_feature)) %>%
          dplyr::summarise(dr1_avg = mean(!!rlang::sym(paste0(reduction, "_", dims[1]))),
                           dr2_avg = mean(!!rlang::sym(paste0(reduction, "_", dims[2]))),
                           pct = (sum(!!rlang::sym(colnames(.)[1]) > 0)/dplyr::n())*100) %>%
          dplyr::mutate(pct = ifelse(pct < 1, ifelse(pct == 0, "0 %", "< 1 %"), paste0(round(pct,0), " %")))


        group_labels[,"dr1_avg"] <- group_labels[,"dr1_avg"] + contour.label.nudge[1]
        group_labels[,"dr2_avg"] <- group_labels[,"dr2_avg"] + contour.label.nudge[2]


        if (use_ggnewscale_for_contour_colors) {
          plot <-
            plot +
            ggplot2::geom_label(data = group_labels, ggplot2::aes(x = dr1_avg, y = dr2_avg, label = pct, color = !!rlang::sym(names(group_labels)[1])),
                                show.legend = F)
        } else {
          plot <-
            plot +
            ggplot2::geom_label(data = group_labels, ggplot2::aes(x = dr1_avg, y = dr2_avg, label = pct))
        }

      }
    }


    if (!is.null(plot.labels)) {
      if (is.numeric(data[,1])) {
        message("Labels not plotted as ", x, " is numeric.")
      } else {
        label.center.fun <- match.fun(label.center.fun)
        # potentially: use mclust, densityMclust(), (mixed gaussian model) to find multimodal clusters; x,y separately
        '        label_df <- do.call(rbind, lapply(unique(data[,1]), function(z) data.frame(label = z,
                                                                                   avg1 = label.center.fun(data[which(data[,1] == z), paste0(reduction, "_", dims[1])]),
                                                                                   avg2 = label.center.fun(data[which(data[,1] == z), paste0(reduction, "_", dims[2])]))))'

        dimcol1 <- paste0(reduction, "_", dims[1])
        dimcol2 <- paste0(reduction, "_", dims[2])
        label_df <-
          data %>%
          dplyr::group_by(!!sym(x), SO.split) %>%
          dplyr::summarise(!!dimcol1 := label.center.fun(!!sym(dimcol1)) + label.nudge[1],
                           !!dimcol2 := label.center.fun(!!sym(dimcol2)) + label.nudge[2],
                           .groups = "drop") %>%
          dplyr::rename("label" = !!sym(x)) %>%
          as.data.frame()

        if (!is.null(label.feature)) {
          temp <- unique(data[,c(1,which(colnames(data) == "label.feature")[1])])
          label_df[,1] <- as.character(label_df[,1])
          label_df[,1] <- stats::setNames(temp[,2,drop=T], temp[,1,drop=T])[label_df[,1]]
        }

        if (plot.labels == "text") {
          plot <- plot + ggplot2::geom_text(data = label_df, ggplot2::aes(label = label), size = label.size, family = font.family) #...
        }
        if (plot.labels == "label") {
          plot <- plot + ggplot2::geom_label(data = label_df, ggplot2::aes(label = label), size = label.size, family = font.family) #...
        }
      }
    }

    if (length(expand_limits) > 0) {
      plot <-
        plot +
        Gmisc::fastDoCall(ggplot2::expand_limits, args = expand_limits)
    }

    # define facets and plot freq.of.expr annotation
    wrap_by <- function(...) {ggplot2::facet_wrap(ggplot2::vars(...), labeller = ggplot2::label_wrap_gen(multi_line = F), scales = split.by.scales, nrow = nrow.inner, ncol = ncol.inner)}
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
    #plots <- cowplot::plot_grid(plotlist = plots, ncol = ncol.combine, nrow = nrow.combine, align = "hv", axis = "tblr")
    plots <- patchwork::wrap_plots(plots, ncol = ncol.combine, nrow = nrow.combine)
  }
  if (length(plots) == 1 && !combine) {
    plots <- plots[[1]]
  }
  return(plots)
}


.check.SO <- function(SO,
                      assay,
                      split.by = NULL,
                      shape.by = NULL,
                      meta.col = NULL,
                      length = NULL) {
  if (!is.list(SO)) {
    SO <- list(SO)
  }

  if (!is.null(length)) {
    if (length(SO) != length) {
      stop("Please provide exactly ", length, " SO.")
    }
  }

  assay <- match.arg(assay, Reduce(intersect, lapply(SO, function(x) names(x@assays))))

  SO <- lapply(SO, function(x) {
    Seurat::DefaultAssay(x) <- assay
    return(x)
  })

  if (is.null(names(SO)) && length(SO) > 1) {
    message("List of SO has no names. Naming them numerically in order as provided.")
    names(SO) <- as.character(seq_along(SO))
  }

  if (is.null(names(SO))) {
    names(SO) <- as.character(seq_along(SO))
  }

  if (!any(unlist(lapply(SO, class)) == "Seurat")) {
    stop("All SO have to be Seurat objects (class == Seurat).")
  }
  if (!is.null(split.by) && length(.check.features(SO, split.by, rownames = F)) == 0) {
    stop("split.by not found in all objects.")
  }
  if (!is.null(shape.by) && length(.check.features(SO, shape.by, rownames = F)) == 0) {
    stop("shape.by not found in all objects.")
  }
  if (!is.null(meta.col) && length(.check.features(SO, meta.col, rownames = F)) == 0) {
    stop("meta.col not found in all objects.")
  }

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

.check.and.get.cells <- function(SO,
                                 assay = c("RNA", "SCT"),
                                 #make.cells.unique = F,
                                 cells = NULL,
                                 cutoff.feature = NULL,
                                 cutoff.expression = NULL,
                                 exclusion.feature  = NULL,
                                 downsample = 1,
                                 #make.cells.unique.warning = 1,
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
    message("Duplicates found in cells. Cells is made unique though.")
    cells <- unique(cells)
  }

  # check if cell names are unique across SOs
  all.cells <- unlist(lapply(SO, function(x) Seurat::Cells(x)))

  # check if any cells intersects with not unique cell names in SOs
  if (length(SO) > 1 && any(duplicated(all.cells)) &&
      !is.null(cells) && any(cells %in% all.cells[which(duplicated(all.cells))])) {
    stop("Selected cells can not be identified unambiguously as they intersect/overlap with duplicate cell names (barcodes) across SOs. Please fix with Seurat::RenameCells and make cell names unique.")
  }

  '  if (length(SO) > 1 && any(duplicated(all.cells)) && !make.cells.unique) {
    if (make.cells.unique.warning == 1) {
      stop("Cell names are not unique across SOs. Please fix that manually with Seurat::RenameCells or pass make.cells.unique = T when calling this function. Cells are then renamed with the prefix SO_i_, where i is the index (starting with 1) of the SO in the list. Consider this way of renaming cells when selecting cells for plotting or other calculations.")
    } else {
      stop("Cell names are not unique across SOs. Please fix that manually with Seurat::RenameCells.")
    }
  }'

  if (length(SO) > 1 && any(duplicated(all.cells))) {
    names.temp <- names(SO)
    SO <- lapply(seq_along(SO), function(x) Seurat::RenameCells(SO[[x]], add.cell.id = paste0("SO_", x)))
    names(SO) <- names.temp
    all.cells <- unlist(lapply(SO, function(x) Seurat::Cells(x)))
    #print("Cells have been prefixed with 'SO_i_' .")
    assign("SO", SO, envir = parent.frame()) # assigns in parent environment (https://stackoverflow.com/questions/10904124/global-and-local-variables-in-r?rq=1)
    if (!is.null(cells)) {
      # change names of selected cells as well (above it was confirmed that they do not belong to duplicate names, so this is safe)
      cells <- grep(pattern = paste(paste0("SO_[[:digit:]]{1,}_", cells, "$"), collapse = "|"), x = all.cells, value = T)
    }
  }

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
  ## cells excluded by

  if (!is.null(cutoff.feature)) {
    compare_fun <- function(x,y) {x > y}
    exclude.cells <- names(which(Matrix::colSums(sweep(do.call(cbind, lapply(SO, function(x) Seurat::GetAssayData(x, slot = "data", assay = assay)[cutoff.feature,,drop=F])), 1, cutoff.expression, compare_fun)) < length(cutoff.feature)))
    all.cells[which(names(all.cells) %in% exclude.cells)] <- 0
  }
  if (length(all.cells) == 0) {
    stop("No cells left after filtering for cutoff.feature.")
  }
  ## cells excluded due to expression of an unwanted gene
  if (!is.null(exclusion.feature)) {
    exclude.cells <- unlist(lapply(SO, function(x) names(which(Matrix::colSums(do.call(cbind, lapply(SO, function(x) Seurat::GetAssayData(x, slot = "data", assay = assay)[exclusion.feature,,drop=F]))) > 0))))
    all.cells[which(names(all.cells) %in% exclude.cells)] <- 0
  }

  if (length(all.cells) == 0) {
    stop("No cells left after filtering for exclusion.feature")
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

  if (return.included.cells.only) {
    return(names(all.cells[which(all.cells == 1)]))
  } else if (!return.included.cells.only) {
    return(all.cells)
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
                            meta.data = T,
                            meta.data.numeric = F) {

  if (is.null(features)) {return(NULL)}
  if (!rownames && !meta.data) {return((NULL))}
  if (!is.list(SO)) {
    SO <- list(SO)
  }
  features <- unique(features)

  features.out <- .feat.check(SO = SO,
                              features = features,
                              rownames = rownames,
                              meta.data = meta.data,
                              ignore.case = T,
                              meta.data.numeric = meta.data.numeric)

  if (length(features.out) == 0) {
    stop("Non of the provided features has not been found in every SO. No features left to plot.")
  }

  ### TO DO WITH meta.data.numeric
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


.feat.check <- function(SO, features, rownames, meta.data, ignore.case, meta.data.numeric) {
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

.get.data <- function(SO,
                      feature,
                      label.feature = NULL,
                      assay = c("RNA", "SCT"),
                      slot = "data",
                      cells = NULL,
                      split.by = NULL,
                      shape.by = NULL,
                      meta.col = NULL,
                      reduction = "umap",
                      min.q.cutoff = 0,
                      max.q.cutoff = 1,
                      order = T,
                      order.rev = F,
                      order.abs = T,
                      shuffle = F,
                      order.discrete = T,
                      bury_NA = T,
                      na.rm = F,
                      inf.rm = F,
                      trajectory.slot = NULL) {


  SO <- .check.SO(SO = SO, assay = assay, split.by = split.by, shape.by = shape.by, meta.col = meta.col)
  assay <- match.arg(assay, names(SO[[1]]@assays))

  if (length(feature) > 1) {
    # if function is called to get many features
    # ordering becomes irrelevant
    # also no values should be removed then (NA or Inf)
    na.rm <- F
    inf.rm <- F
    order <- F
    order.discrete <- F
    shuffle <- F
  }

  if (max.q.cutoff > 1) {
    #message("max.q.cutoff and min.q.cutoff are divided by 100. Please provide values between 0 and 1.")
    max.q.cutoff <- max.q.cutoff/100
    min.q.cutoff <- min.q.cutoff/100
  }

  all_gene_features <- Reduce(intersect, lapply(SO, rownames))
  all_meta_features <- Reduce(intersect, lapply(SO, function(x) names(x@meta.data)))

  if (any(!feature %in% all_gene_features & !feature %in% all_meta_features)) {
    message("Feature(s) ", paste(feature[which(!feature %in% all_gene_features & !feature %in% all_meta_features)], collapse = ","), " not found.")
    feature <- feature[which(feature %in% all_gene_features & feature %in% all_meta_features)]
    if (length(feature) == 0) {
      stop("No feature left.")
    }
  }

  if (any(feature %in% all_gene_features & feature %in% all_meta_features)) {
    message("Feature(s) ", paste(feature[which(feature %in% all_gene_features & feature %in% all_meta_features)], collapse = ","), " found in expression and meta data. This is not handled at the moment and filtered. Consider renaming the column of meta data.")
    feature <- feature[which(!(feature %in% all_gene_features & feature %in% all_meta_features))]
    if (length(feature) == 0) {
      stop("No feature left.")
    }
  }

  gene_features <- feature[which(feature %in% all_gene_features)]
  meta_features <- feature[which(feature %in% all_meta_features)]

  data <- do.call(rbind, lapply(names(SO), function(y) {

    data <- cbind(data.frame(t(as.matrix(Seurat::GetAssayData(SO[[y]], slot = slot, assay = assay)[gene_features,,drop = F])), check.names = F),
                  data.frame(SO[[y]]@meta.data[,meta_features,drop=F], stringsAsFactors = F, check.names = F))

    if (!is.null(reduction)) {
      reduction <- unique(unlist(lapply(reduction, function(z) {
        names(SO[[y]]@reductions)[which.min(utils::adist(z, names(SO[[y]]@reductions), ignore.case = T))]
      })))
      for (i in reduction) {
        data <- cbind(data, Seurat::Embeddings(SO[[y]], reduction = i))
      }
    }


    ## check if these exist
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

  if (!is.null(label.feature)) {
    data <- cbind(data, do.call(rbind, lapply(names(SO), function(y) {
      data.frame(label.feature = SO[[y]]@meta.data[,label.feature,drop=T], stringsAsFactors = F, check.names = F)
    })))
    if (nrow(unique(data[,which(colnames(data) %in% c(feature, "label.feature"))])) !=
        nrow(unique(data[,which(colnames(data) %in% c(feature)), drop=F]))) {
      stop("label.feature entries do not exactly one feature entry each. This must be the case though.")
    }
  }

  # ensure that facet ordering is according to the order of SO objects provided
  data$SO.split <- factor(data$SO.split, levels = names(SO))

  if (!is.null(cells)) {
    data <- data[cells,]
  }

  if (length(feature == 1) && is.numeric(data[,1]) && all(data[,1] == 0)) {
    message("No expressers found for ", feature, ".")
  }

  if (length(feature == 1) && anyNA(data[,1])) {
    message(feature, ": NA found in data.")
    if (na.rm) {
      data <- data[which(!is.na(data[,1])),]
    }
  }
  if (length(feature == 1) && any(is.infinite(data[,1]))) {
    message(feature, ": Inf found in data.")
    if (inf.rm) {
      data <- data[which(!is.infinite(data[,1])),]
    }
  }


  # use squishing to dampen extreme values - this will produce actually wrong limits on the legend
  if (min.q.cutoff > 0 || max.q.cutoff < 1) {
    for (i in intersect(names(which(sapply(data, is.numeric))), c(gene_features, meta_features))) {
      if (all(data[,i] >= 0)) { # > 0 or >= 0 ?!
        # expression is always greater than 0 and non-expresser are excluded
        data[,i][which(data[,i] > 0)] <- scales::squish(data[,i][which(data[,i] > 0)], range = c(stats::quantile(data[,i][which(data[,i] > 0)], min.q.cutoff), stats::quantile(data[,i][which(data[,i] > 0)], max.q.cutoff)))
      } else {
        # e.g. for module scores below 0
        data[,i] <- scales::squish(data[,i], range = c(stats::quantile(data[,i], min.q.cutoff), stats::quantile(data[,i], max.q.cutoff)))
      }
    }
  }

  ## params:
  # order = T --> ordering by values (lowest are negative values, then 0, then positive ones)
  # order.abs T --> absolute values far away from zero are plotted on top
  # shuffle = T --> only considered when !order; will shuffle the data data.frame; if !order and !shuffle a custom order can be provided from outside
  # order.rev = T --> reverse the order so that lowest values (or zeros if order.abs = T) are on top

  if (order && is.numeric(data[,1])) {
    # per default order will put NAs on top
    # abs: in case negative values are contained in meta.col, any extreme away from 0 will be plotted on top
    if (order.abs) {
      data <- data[order(abs(data[,1]), decreasing = order.rev, na.last = !bury_NA),]
    } else {
      data <- data[order(data[,1], decreasing = order.rev, na.last = !bury_NA),]
    }
  } else if (shuffle) {
    data <- data[sample(x = 1:nrow(data), size = nrow(data), replace = F),]
  }

  ## this is only for meta features; combine with if else from above??
  if (!is.numeric(data[,1]) && is.logical(order.discrete) && order.discrete) {
    # NA is not considered by split()
    # replace it by a character value just for splitting, undo this afterwards
    #data[,1] <- factor(data[,1], exclude = c())

    if (anyNA(data[,1])) {
      na_replace <- "NA"
      while(na_replace %in% unique(data[,1])) {
        na_replace <- paste(c(na_replace, na_replace), collapse = "_")
      }
      level_order <- NULL
      if (is.factor(data[,1])) {
        level_order <- levels(data[,1])
        data[,1] <- as.character(data[,1])
      }
      data[which(is.na(data[,1])),1] <- na_replace
      data <- dplyr::bind_rows(split(data, data[,1])[names(sort(table(data[,1]), decreasing = T))])
      data[which(data[,1] == na_replace),1] <- NA
      if (bury_NA) {
        data <- rbind(data[which(is.na(data[,1])),], data[which(!is.na(data[,1])),])
      }
      if (!is.null(level_order)) {
        data[,1] <- factor(data[,1], levels = level_order)
      }
    } else {
      data <- dplyr::bind_rows(split(data, data[,1])[names(sort(table(data[,1]), decreasing = T))])
    }
  } else if (!is.numeric(data[,1]) && is.logical(order.discrete) && shuffle) {
    data <- data[sample(x = 1:nrow(data), size = nrow(data), replace = F),]
  }

  if (!is.null(trajectory.slot)) {
    data_traj <- do.call(rbind, lapply(names(SO), function(y) {
      # rbind will throw error if column names to not match
      data <- NULL
      if (trajectory.slot %in% names(Seurat::Misc(SO[[y]]))) {
        data <- Seurat::Misc(SO[[y]], trajectory.slot)[["df"]]
        data[,"SO.split"] <- y
      } else {
        message("Trajectory slot not found in Seurat::Misc.")
      }
      return(data)
    }))
    return(list(data = data, data_traj = data_traj))
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

'.prep.contour.cells <- function(SO, reduction, x) {
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
'
ceiling_any = function(x, accuracy, f = ceiling){f(x/ accuracy) * accuracy}

floor_any = function(x, accuracy, f = floor){f(x/ accuracy) * accuracy}

list_depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,list_depth,thisdepth=thisdepth+1))))
  }
}
