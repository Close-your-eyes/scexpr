#' Title
#'
#' @param SO Seurat object
#' @param meta.col which column from meta.data of SO to use as x-axis; if NULL current Idents(SO) are used
#' @param levels.calc which levels in meta.col to include in calculation;
#' all levels if NULL; the order provided defines order on x-axis; level selection
#' will affect scaling calculations
#' @param levels.plot which levels in meta.col to include in platting;
#' if NULL this equals to levels.calc; must be subset of levels.calc;
#' the order provided defines order on x-axis; this will not affect scaling
#' calculation but only select levels for plotting;
#' defining ordering will not work perfectly when more than 1 SO is provided and factor levels
#' in meta.cols are not unique
#' @param assay which assay to obtain expression values from
#' @param features optionally choose which features to plot (supervised)
#' @param normalization how to scale expression values; may be "scale" to
#' use base::scale and transform average expression value to a standardized
#' normal distribution or may be a numeric vector of length 2 for scaling
#' each feature from min (first value) to max (second value); in the latter case
#' c(-1,1) is most meaningful
#' @param topn.features if no features are selected, this will select how many
#' features to plot per level in meta.col; selection is done based on the metric
#' selected in topn.metric; respective features with greatest difference between
#' meta.col levels are selected (best DE features so to say)
#' @param topn.metric which differential expression metric to apply for feature
#' selection
#' @param min.pct filter features across levels in meta.col for a minimal fraction
#' of expressing cells; 0.1 stands for min 10 % expressing cells
#' @param max.padj filter features across levels in meta.col for a significance
#' level of differential expression (one level vs. all others)
#' @param title which title to add to plot; if NULL no title is plotted
#' @param title.font.size font size of title
#' @param font.family which font type (family) to use for plotting,
#' e.g. mono or sans
#' @param y.font.size font size of features names on y-axis
#' @param color color of stroke (border) around tiles or dots; "NA" means no stroke is plotted;
#' other choices may be black, white or any other color code
#' @param plot.feature.breaks
#' @param plot.sec.axis
#' @param legend.position
#' @param legend.direction
#' @param legend.barheight
#' @param legend.barwidth
#' @param legend.text.size
#' @param legend.title.text.size
#' @param legend.title.fill
#' @param legend.title.size
#' @param legend.title.position
#' @param legend.labels
#' @param feature.labels
#' @param feature.labels.nudge_x
#' @param feature.labels.axis.width
#' @param ...
#' @param order_features when features are provided should they be ordered automatically to generate a pretty heatmap;
#' not relevant when is.null(features)
#' @param plot_hlines_between_groups plot horizontal lines to separate groups of genes that are
#' differentially expressed by a group on the x-axis; if hlines is NULL the breaks are
#' determined automatically
#' @param break.ties if FALSE, ties are not broken; if TRUE a second and third metric
#' beside topn.metric from c("padj", "logFC", "auc") are used
#' @param hlines provide manual breaks of where to plot horizontal lines;
#' even though the y-axis is discrete, a numeric vector works here:
#' hlines = c(1.5, 3.5) will plot lines between the first and second, and, third and fourth
#' entry on the y-axis, respectively
#' @param hlines_args arguments to geom_hline like linewidth or linetype
#' @param feature_selection_strategy different strategies to auto-select features; rule-of-thumb:
#' for large overview heatmap (100s of features per group) choose 2, for small detailed heatmap (10 per group) choose 1;
#' when a selection of features is provided by the 'feature' argument and order.features is TURE, feature_selection_strategy will be set to 2
#' by default if no choice is made as this gives a better order of features on the y-axis; in general choices 1 and 2 may yield different
#' features on the y-axis; choice 1 is slightly preferred
#' @param fill
#' @param flip.axes
#' @param legend.title.hjust
#' @param n.colorsteps number of steps (numeric) to divide color scale into; if null then ordinary continuous fill scale is chosen;
#' if of length 1 this is passed as n.breaks to scale_fill_stepsn, if length > 1 then passed as breaks to scale_fill_stepsn
#' @param nice.breaks passed to scale_fill_stepsn if length(n.colorsteps) == 1
#' @param show.limits passed to scale_fill_stepsn if length(n.colorsteps) > 1; show min and max limit on legend
#' @param dotplot
#' @param legend.decimals passed to scale_fill_stepsn if length(n.colorsteps) > 1; number of decimals to round legend labels to
#'
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
heatmap_pseudobulk <- function(SO,
                               meta.col = NULL,
                               levels.calc = NULL,
                               levels.plot = NULL,
                               assay = c("RNA", "SCT"),
                               features = NULL,
                               feature_selection_strategy = c(2,1),
                               order_features = F,
                               normalization = "scale", # scale
                               topn.features = 10,
                               break.ties = c(T, F),
                               topn.metric = c("padj", "logFC", "auc"),
                               min.pct = 0.1,
                               max.padj = 0.05,
                               title = NULL,
                               title.font.size = 14,
                               font.family = "sans",

                               y.font.size = 10,
                               color = "NA",
                               fill = scexpr::col_pal(name = "RdBu", nbrew = 9, reverse = T),
                               n.colorsteps = NULL,
                               nice.breaks = F,
                               show.limits = T,
                               legend.decimals = 1,
                               dotplot = F,
                               plot.feature.breaks = T,
                               plot.sec.axis = F,

                               flip.axes = F,

                               legend.position = "right",
                               legend.direction = "vertical",
                               legend.barheight = 8,
                               legend.barwidth = 1,
                               legend.text.size = 10,
                               legend.title.text.size = 10,
                               legend.title.fill = NULL,
                               legend.title.size = NULL,
                               legend.title.position = "top",
                               legend.title.hjust = 0.5,
                               legend.labels = c("min", "", "max"),
                               feature.labels = NULL,
                               feature.labels.nudge_x = -0.1,
                               feature.labels.axis.width = 0.2,
                               plot_hlines_between_groups = F,
                               hlines = NULL,
                               hlines_args = list(),
                               ...) {

  # ... arguments to ggrepel, like nudge_y and scexpr::convert_gene_identifier

  if (!requireNamespace("devtools", quietly = T)) {
    utils::install.packages("devtools")
  }
  if (!requireNamespace("presto", quietly = T)) {
    devtools::install_github("immunogenomics/presto")
  }

  # set to 2 when features are provided
  # this was found to reliably yield a beautiful order
  if (!is.null(features) && length(feature_selection_strategy) == 2) {
    feature_selection_strategy <- 2
  }

  topn.metric <- match.arg(topn.metric, c("padj", "logFC", "auc"))
  feature_selection_strategy <- match.arg(as.character(feature_selection_strategy), c("2","1"))

  if (!is.logical(break.ties)) {
    stop("break.ties must be logical.")
  }
  break.ties <- match.arg(as.character(break.ties), c("TRUE", "FALSE"))
  break.ties <- as.logical(break.ties)

  if (topn.features < 1) {
    message("topn.features has to be minimum 1.")
    topn.features <- 1
  }


  if (plot_hlines_between_groups && !is.null(features) && !order_features && is.null(hlines)) {
    message("horizontal lines between group (hlines) can only be inferred automatically if order_features is set to TRUE.
            Will set plot_hlines_between_groups to FALSE now.
            Alternatively provide them with the hlines argument manually")
    plot_hlines_between_groups <- F
  }

  if (!plot_hlines_between_groups && !is.null(hlines)) {
    message("hlines provided but plot_hlines_between_groups is FALSE.")
  }


  SO <- .check.SO(SO = SO, assay = assay) #length = 1
  assay <- match.arg(assay, Reduce(intersect, lapply(SO, function(x) names(x@assays))))

  if (methods::is(normalization, "numeric")) {
    if (length(normalization) != 2) {
      stop("Please set normalization to 'scale' or a numeric vector of length 2 (e.g. c(-1,1)).")
    }
    if (length(legend.labels) > 3) {
      stop("legend.labels can have max length of 3 only")
    }
    if (length(legend.labels) == 2) {
      legend.labels <- c(legend.labels[1], "", legend.labels[2])
    }
  } else {
    normalization <- match.arg(normalization, choices = c("scale"))
  }


  ## checking meta.col, levels.plot and levels.calc
  if (is.null(meta.col)) {
    SO <- lapply(SO, function(x) {
      x@meta.data$idents <- Seurat::Idents(x)
      return(x)
    })
    meta.col <- rep("idents", length(SO))
  } else {
    if (length(meta.col) != length(SO)) {
      stop("meta.col has to have same length as SO. (One column name in meta.data for each Seurat object.)")
    }
    if (any(unlist(mapply(x = SO, y = meta.col, function(x,y) {
      !y %in% names(x@meta.data)
    }, SIMPLIFY = F)))) {
      stop("One of meta.col not found in respective SO.")
    }
    if (any(unlist(mapply(x = SO, y = meta.col, function(x,y) {
      is.numeric(x@meta.data[,y,drop=T])
    }, SIMPLIFY = F)))) {
      stop("One of meta.col is numeric. Please make it factor or character.")
    }
  }


  # this will only apply if length(SO)==1 and levels.plot and levels.calc are vectors
  if (!is.null(levels.plot) && !is.list(levels.plot)) {
    levels.plot <- list(levels.plot)
  }
  if (!is.null(levels.calc) && !is.list(levels.calc)) {
    levels.calc <- list(levels.calc)
  }


  if (!is.null(levels.plot) && !is.null(levels.calc)) {
    levels.plot <- mapply(x = levels.plot, y = levels.calc, function(x,y) {
      if (any(!x %in% y)) {
        print("levels.plot has more levels than levels.calc. levels.plot is reduced to levels.calc.")
        x <- x[which(x %in% y)]
      }
      return(x)
    }, SIMPLIFY = F)
  }

  if (!is.null(levels.calc)) {
    if (length(levels.calc) != length(SO)) {
      stop("levels.calc has to have same length as SO.")
    }
  } else {
    levels.calc <- rep(NA, length(SO))
  }
  if (!is.null(levels.plot)) {
    if (length(levels.plot) != length(SO)) {
      stop("levels.plot has to have same length as SO.")
    }
  } else {
    levels.plot <- rep(NA, length(SO))
  }

  # this will replace NULL with all available levels
  levels.calc <- purrr::pmap(list(x = SO, y = levels.calc, z = meta.col), function(x,y,z) {
    .check.levels(SO = x, meta.col = z, levels = y, append_by_missing = F)
  })
  levels.plot <- purrr::pmap(list(x = SO, y = levels.plot, z = meta.col), function(x,y,z) {
    .check.levels(SO = x, meta.col = z, levels = y, append_by_missing = F)
  })
  #levels.calc <- .check.levels(SO = SO, meta.col = meta.col, levels = levels.calc, append_by_missing = F)
  #levels.plot <- .check.levels(SO = SO, meta.col = meta.col, levels = levels.plot, append_by_missing = F)

  # subset levels for calculation
  SO <- purrr::pmap(list(x = SO, y = levels.calc, z = meta.col), function(x,y,z) {
    if (length(y) < length(unique(x@meta.data[,z,drop=T]))) {
      x <- subset(x, cells = rownames(x@meta.data[,z,drop=F][which(x@meta.data[,z] %in% y),,drop=F]))
    }
    return(x)
  })

  #SO <- subset(SO, cells = rownames(SO@meta.data[,meta.col,drop=F][which(SO@meta.data[,meta.col] %in% levels.calc),,drop=F]))

  # this problem was discovered by chance (":" are replaced by Seurat::AverageExpression by "_"); maybe other symbols will also cause problems.
  for (i in length(SO)) {
    if (any(grepl(":", unique(SO[[i]]@meta.data[,meta.col[[i]]])))) {
      message("Found colon (:) in factor levels of meta.col. Will replace those with underscores (_).")
      SO[[i]]@meta.data[,meta.col[[i]]] <- gsub(":", "_", SO[[i]]@meta.data[,meta.col[[i]]])
      levels.plot[[i]] <- gsub(":", "_", levels.plot[[i]])
      levels.calc[[i]] <- gsub(":", "_", levels.calc[[i]])
    }
  }

  # make cluster names unique if there is intersection
  ## this will prohibit to define cluster order
  if (length(SO) > 1) {
    if (length(Reduce(intersect, levels.plot)) > 1 || length(Reduce(intersect, levels.calc)) > 1) {
      SO <- purrr::pmap(list(x = SO, y = meta.col, z = names(SO)), function(x,y,z) {
        x@meta.data[,y] <- paste0(z, "___", x@meta.data[,y,drop=T])
        return(x)
      })
    }
    if (length(Reduce(intersect, levels.plot)) > 1) {
      levels.plot <- mapply(x = names(SO), y = levels.plot, function(x,y) paste0(x, "___", y), SIMPLIFY = F)
    }
    if (length(Reduce(intersect, levels.calc)) > 1) {
      levels.calc <- mapply(x = names(SO), y = levels.calc, function(x,y) paste0(x, "___", y), SIMPLIFY = F)
    }
  }

  ###

  legend.direction <- match.arg(legend.direction, c("horizontal", "vertical"))
  if (legend.direction == "horizontal") {
    temp <- legend.barwidth
    legend.barwidth <- legend.barheight
    legend.barheight <- temp
  }

  ## presto gives deviating results with respect to avgExpr and logFC (maybe due to approximation which makes calculation faster)
  ## nevertheless, in order to select marker genes by logFC or other statistics fast performance of presto is useful
  ## filter all features with zero expression across the data set

  #expr_mat <- do.call(cbind, lapply(SO, function(x) {Seurat::GetAssayData(x, slot = "data", assay = assay)}))
  #unlist(purrr::pmap(list(x = SO, y = meta.col, z = names(SO)), function(x,y,z) {paste0(z, "___", x@meta.data[,y,drop=T])})) ~ SO@meta.data[,meta.col,drop=T]
  wil_auc <-
    presto::wilcoxauc(X = do.call(cbind, lapply(SO, function(x) {Seurat::GetAssayData(x, slot = "data", assay = assay)})),
                      y = unlist(purrr::pmap(list(x = SO, y = meta.col), function(x,y) {x@meta.data[,y,drop=T]}))) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(sum(pct_in) > 0) %>% # filter non-expressed features
    dplyr::ungroup()

  # prep features variable here to have one common pipeline below
  if (!is.null(features)) {
    features <- .check.features(SO = SO, features = unique(features), meta.data = F)
    if (any(!features %in% wil_auc$feature)) {
      message("No expressers found for: ", paste(features[which(!features %in% wil_auc$feature)], collapse = ","), ". Will not be plotted.")
    }
  } else {
    order_features <- T
    features <- wil_auc$feature
  }

  if (order_features) {

    # select multiple columns for ordering, in order to break ties; especially if padj is used as primary element to order
    if (break.ties) {
      all_metrics <- c("padj", "logFC", "auc")
      second_metric <- all_metrics[which(!all_metrics %in% topn.metric)][1]
      third_metric <- all_metrics[which(!all_metrics %in% topn.metric)][2]
    }


    if (feature_selection_strategy == "2") {
      # this is better for when features are provided; yields better order on y-axis
      features2 <-
        wil_auc %>%
        dplyr::filter(feature %in% features) %>%
        dplyr::filter(pct_in >= min.pct) %>%
        dplyr::filter(padj <= max.padj) %>%
        dplyr::mutate(padj = -padj) %>% # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc"))
        dplyr::group_by(feature) %>%
        dplyr::slice_max(order_by = avgExpr, n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(group = factor(group, levels = unlist(levels.plot))) %>%
        dplyr::group_by(group)
      if (break.ties) {
        features3 <-
          features2 %>%
          dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = topn.features)
      } else {
        features3 <-
          features2 %>%
          dplyr::slice_max(order_by = !!rlang::sym(topn.metric), n = topn.features, with_ties = T)
      }
      features3 <-
        features3 %>%
        dplyr::ungroup() %>%
        dplyr::arrange(group, avgExpr) #!!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)

      ## and maybe add additional ones, but how?
      features_check <-
        features3 %>%
        dplyr::group_by(feature) %>%
        dplyr::slice_max(order_by = avgExpr, n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(group) %>%
        dplyr::count()

    } else if (feature_selection_strategy == "1") {
      ## alternative feature selection; this procedure avoid that a diff gene is "stolen" by another group (due to max avgExpr) where it does not lie within topn.features
      ## this would cause that a feature actually relevant (diff expressed by 2 group does not appear on the heatmap)
      ## it may happen that a few more features per group are plotted
      ## for very large heatmaps (overview with 100 of features) is may look ugly; then feature_selection_strategy 2 may be better suited
      features2 <-
        wil_auc %>%
        dplyr::filter(feature %in% features) %>%
        dplyr::filter(pct_in >= min.pct) %>%
        dplyr::filter(padj <= max.padj) %>%
        dplyr::mutate(padj = -padj) %>% # make negative so that slice_max works equally for all three metrics (c("padj", "logFC", "auc"))
        dplyr::mutate(group = factor(group, levels = unlist(levels.plot))) %>%
        dplyr::group_by(group)

      if (break.ties) {
        features3 <-
          features2 %>%
          dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = topn.features) %>%
          dplyr::ungroup()

        ## check for number of features per group
        ## and maybe add additional ones, but how?
        features_check <-
          features3 %>%
          dplyr::group_by(feature, .drop = F) %>%
          dplyr::slice_max(order_by = avgExpr, n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(group) %>%
          dplyr::count()

        # try to make sure topn.features is reach for every group
        max_round <- min(topn.features, 20)
        n <- 1
        while(any(features_check$n < topn.features) && n <= max_round) {
          for (i in features_check[which(features_check$n<topn.features),"group",drop=T]) {
            nplus <- topn.features - features_check[which(features_check$group == i),"n",drop=T]
            # in case factor level i is missing in features_check
            if (length(nplus) == 0) {
              nplus <- topn.features
            }
            features_plus <-
              dplyr::anti_join(features2, features3, by = dplyr::join_by(feature, group, avgExpr, logFC, statistic, auc, pval, padj, pct_in, pct_out)) %>%
              dplyr::filter(group == i) %>%
              dplyr::slice_max(order_by = tibble::tibble(!!rlang::sym(topn.metric), !!rlang::sym(second_metric), !!rlang::sym(third_metric)), n = nplus)

            features3 <-
              rbind(features3, features_plus) %>%
              dplyr::distinct()

            features_check <-
              features3 %>%
              dplyr::group_by(feature) %>%
              dplyr::slice_max(order_by = avgExpr, n = 1) %>%
              dplyr::ungroup() %>%
              dplyr::group_by(group) %>%
              dplyr::count()
          }
          n<-n+1
        }
        features3 <-
          features3 %>%
          dplyr::arrange(group, avgExpr) # avgExpr instead of metric?! # -!!rlang::sym(topn.metric), -!!rlang::sym(second_metric), -!!rlang::sym(third_metric)

      } else {

        features3 <-
          features2 %>%
          dplyr::slice_max(order_by = !!rlang::sym(topn.metric), n = topn.features, with_ties = T) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(group, avgExpr) # avgExpr instead of metric?! #-!!rlang::sym(topn.metric)

        features_check <-
          features3 %>%
          dplyr::group_by(feature) %>%
          dplyr::slice_max(order_by = avgExpr, n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(group) %>%
          dplyr::count()
      }
    }

    print(as.data.frame(features_check))

    ## continue
    if (is.null(hlines)) {
      hlines <- cumsum(rle(as.character(features3$group))[["lengths"]]) + 0.5
      hlines <- hlines[-length(hlines)]
    }

    features <-
      features3 %>%
      dplyr::pull(feature)
  }

  raw_tab <- do.call(cbind, mapply(x = SO, y = meta.col, function(x,y) {
    Seurat::AverageExpression(x, assays = assay, group.by = y, slot = "data", verbose = F)[[1]][features,,drop=F]
  }, SIMPLIFY = F))
  #raw_tab <- Seurat::AverageExpression(SO, assays = assay, group.by = meta.col, slot = "data", verbose = F)[[1]][features,,drop=F]

  if (methods::is(normalization, "character") && normalization == "scale") {
    tab <- row_scale(raw_tab, add_attr = F)
  } else if (methods::is(normalization, "numeric")) {
    tab <- scale_min_max(raw_tab, min = min(normalization), max = max(normalization), margin = 1)
  }

  htp <-
    as.data.frame(tab) %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::filter(Feature %in% features) %>%
    tidyr::pivot_longer(cols = -Feature, names_to = "cluster", values_to = "norm_avgexpr") %>%
    dplyr::filter(cluster %in% unlist(levels.plot)) %>%
    dplyr::left_join(wil_auc[,c(which(names(wil_auc) %in% c("feature", "group", "pct_in")))], by = c("Feature" = "feature", "cluster" = "group")) %>%
    dplyr::mutate(cluster = factor(cluster, levels = unlist(levels.plot))) %>%
    dplyr::mutate(Feature = factor(Feature, levels = unique(features)))

  scale.max <- as.numeric(format(floor_any(max(htp$norm_avgexpr), 0.1), nsmall = 1))
  scale.min <- as.numeric(format(ceiling_any(min(htp$norm_avgexpr), 0.1), nsmall = 1))
  #scale.mid <- as.numeric(format(round(scale.min + ((scale.max - scale.min) / 2), 1), nsmall = 1))
  scale.mid <- 0

  if (methods::is(normalization, "numeric")) {
    labels <- c("min", "int", "max")
  } else {
    labels <- c(scale.min, scale.mid, scale.max)
  }

  if (flip.axes) {
    heatmap.plot <-
      ggplot2::ggplot(htp, ggplot2::aes(y = cluster, x = Feature, fill = norm_avgexpr)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(family = font.family),
                     axis.text.x = ggplot2::element_text(size = y.font.size, face = "italic", family = font.family))

  } else {
    heatmap.plot <-
      ggplot2::ggplot(htp, ggplot2::aes(x = cluster, y = Feature, fill = norm_avgexpr)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(family = font.family),
                     axis.text.y = ggplot2::element_text(size = y.font.size, face = "italic", family = font.family))
  }

  heatmap.plot <-
    heatmap.plot +
    ggplot2::ggtitle(title) +
    ggplot2::theme(title = ggplot2::element_text(size = title.font.size, family = font.family),
                   axis.title = ggplot2::element_blank(),
                   legend.position = legend.position, legend.direction = legend.direction)

  if (dotplot) {
    heatmap.plot <- heatmap.plot + ggplot2::geom_point(ggplot2::aes(size = pct_in), shape = 21, color = color)
  } else {

    heatmap.plot <- heatmap.plot + ggplot2::geom_tile(color = color)
  }

  if (!is.null(n.colorsteps)) {
    if (length(n.colorsteps) == 1) {
      heatmap.plot <-
        heatmap.plot +
        ggplot2::scale_fill_stepsn(colors = fill,
                                   n.breaks = n.colorsteps,
                                   nice.breaks = nice.breaks)
    } else {
      heatmap.plot <-
        heatmap.plot +
        ggplot2::scale_fill_stepsn(colors = fill,
                                   breaks = n.colorsteps,
                                   show.limits = show.limits,
                                   labels = function(x) round(x, legend.decimals))
    }
    guide_fun <- ggplot2::guide_colorsteps
  } else {
    heatmap.plot <-
      heatmap.plot +
      ggplot2::scale_fill_gradientn(values = scales::rescale(c(scale.min, scale.mid, scale.max)),
                                    colors = fill,
                                    breaks = c(scale.min, scale.mid, scale.max), labels = labels)
    guide_fun <- ggplot2::guide_colourbar
  }

  heatmap.plot <-
    heatmap.plot +
    ggplot2::guides(fill = guide_fun(label.theme = ggplot2::element_text(size = legend.text.size, family = font.family),
                                     title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family),
                                     title.position = legend.title.position,
                                     title = legend.title.fill,
                                     title.hjust = legend.title.hjust,
                                     barwidth = legend.barwidth,
                                     barheight = legend.barheight,
                                     order = 1),
                    size = ggplot2::guide_legend(label.theme = ggplot2::element_text(size = legend.text.size, family = font.family),
                                                 title.theme = ggplot2::element_text(size = legend.title.text.size, family = font.family),
                                                 title.position = legend.title.position,
                                                 title = legend.title.size,
                                                 title.hjust = legend.title.hjust,
                                                 label.position = "bottom",
                                                 order = 2,
                                                 override.aes = list(color = "black")))


  if (plot.feature.breaks & !is.null(feature.labels)) {
    axis.df <- data.frame(y = 1:length(levels(htp$Feature)), Feature = levels(htp$Feature))
    axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = Feature)) +
      ggrepel::geom_text_repel(fontface = "italic", family = font.family, data = axis.df[which(axis.df$Feature %in% feature.labels),], ggplot2::aes(label = Feature), nudge_x = feature.labels.nudge_x, direction = "y", ...) +
      ggplot2::scale_x_continuous(limits = c(-0.1, 0), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0, length(levels(htp$Feature)) + 0.5), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::theme_void()
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
    heatmap.plot <- cowplot::plot_grid(axis, heatmap.plot, align = "h", axis = "tb", nrow = 1, rel_widths = c(feature.labels.axis.width,1)) ## check how to replace with patchwork
  }

  if (plot.sec.axis) {
    out <- convert_gene_identifier(idents = levels(htp$Feature), ident_in = "SYMBOL", ident_out = "GENENAME", ...)
    #out <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = levels(htp$Feature), keytype = "SYMBOL", column = i, multiVals = "first")
    axis.df <- data.frame(y = 1:length(levels(htp$Feature)), Feature = levels(htp$Feature)) %>% dplyr::left_join(out, by = c("Feature" = "SYMBOL"))
    axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = GENENAME)) +
      ggplot2::geom_text(hjust = 0, family = font.family) +
      ggplot2::scale_x_continuous(limits = c(0, 2), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0.4, length(levels(htp$Feature)) + 0.4), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::theme(panel.background = ggplot2::element_blank(), plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"), axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    heatmap.plot <- cowplot::plot_grid(heatmap.plot + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt")), axis, align = "h", axis = "tb", nrow = 1, rel_widths = c(0.6,0.4))
  }

  if (plot_hlines_between_groups) {
    if (flip.axes) {
      heatmap.plot <-
        heatmap.plot +
        Gmisc::fastDoCall(geom_vline, args = c(list(xintercept = hlines), hlines_args))
    } else {
      heatmap.plot <-
        heatmap.plot +
        Gmisc::fastDoCall(geom_hline, args = c(list(yintercept = hlines), hlines_args))
    }

  } else {
    hlines = NULL
  }


  if (missing(features) & is.null(feature.labels) & nlevels(htp$Feature) > 200 & plot.feature.breaks) {
    print("Large number of rows without selecting feature.labels. Consider setting plot.feature.breaks to FALSE or selecting feature.labels to plot.")
  }

  if (!plot.feature.breaks) {
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }

  return(list(plot = heatmap.plot, data = htp, complete_data = wil_auc, hlines = hlines))
}

.check.levels <- function(SO, meta.col, levels = NULL, append_by_missing = F) {
  if (any(is.null(levels)) || any(is.na(levels))) {
    levels <- as.character(unique(SO@meta.data[,meta.col,drop=T]))
    if (suppressWarnings(!any(is.na(as.numeric(levels))))) {
      levels <- as.character(sort(as.numeric(unique(levels))))
    }
  } else {
    if (any(!levels %in% unique(SO@meta.data[,meta.col,drop=T]))) {
      print(paste0("levels not found in meta.col: ", paste(levels[which(!levels %in% unique(SO@meta.data[,meta.col,drop=T]))], collapse = ", ")))
    }
    levels <- as.character(unique(levels[which(levels %in% unique(SO@meta.data[,meta.col,drop=T]))]))
  }
  if (append_by_missing) {
    levels <- c(levels, as.character(unique(SO@meta.data[,meta.col,drop=T])[which(!as.character(unique(SO@meta.data[,meta.col,drop=T])) %in% levels)]))
  }
  return(levels)
}
