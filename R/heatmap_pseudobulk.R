#' Plot a heatmap of average gene transcription of transcriptomes in different groups (clusters)
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
#' c(-1,1) is most meaningful; set to NULL to not have any normalization computed
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
#' @param color color of stroke (border) around tiles or dots; "NA" means no stroke is plotted; NA has
#' to be put in quotation mark ("NA"), such that geom_point accepts it.
#' other choices may be black, white or any other color code; when "auto" a grey70 is used
#' when dotplot = F and a the number of features is below 100.
#' @param plot.feature.breaks
#' @param plot.sec.axis
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
#' @param fill color scale to use for filling tiles or dots.
#' if auto, then roughly scexpr::col_pal(name = "RColorBrewer::RdBu", n = 11, direction = -1) is used.
#' if !is.null(n.colorsteps) this scale is adjusted to make sure the split between blue and red happens at zero.
#' when a color scale is provided manually a meaningful split of a diverging scale at zero has to be forced manually.
#' @param flip.axes flip x and y axes
#' @param n.colorsteps number of steps (numeric) to divide color scale into; if null then ordinary continuous fill scale is chosen;
#' if of length 1 this is passed as n.breaks to scale_fill_stepsn, if length > 1 then passed as breaks to scale_fill_stepsn
#' @param nice.breaks passed to scale_fill_stepsn if length(n.colorsteps) == 1
#' @param show.limits passed to scale_fill_stepsn if length(n.colorsteps) > 1; show min and max limit on legend
#' @param dotplot do not plot tiles but dots, the size of which then indicates the percentage of transcribing cells
#' @param legend.decimals passed to scale_fill_stepsn if length(n.colorsteps) > 1; number of decimals to round legend labels to
#' @param border_linewidth linewidth (geom_tile) or stroke (geom_point); defines the size of borders around tiles or points
#' @param convert_gene_identifier_args list of named arguments to scexpr::convert_gene_identifier; only needed
#' if plot.sec.axis is TRUE, most likely you may have to set species to Hs or Mm as needed, in case it cannot be guessed.
#' @param theme ggplot theme, has to be passed with brackets
#' @param theme_args arguments to theme
#' @param legend_fill_args arguments passed to ggplot2::guide_colorsteps or ggplot2::guide_colorbar, depend upon the color scale
#' @param legend_size_args arguments passed to ggplot2::guide_legend to modify the size legend, only applies if dotplot = T
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
                               assay = "RNA",
                               features = NULL,
                               feature_selection_strategy = c(2,1),
                               order_features = F,
                               normalization = "scale", # scale, c(-1,1), NULL
                               topn.features = 5,
                               break.ties = c(T, F),
                               topn.metric = c("padj", "logFC", "auc"),
                               min.pct = 0.1,
                               max.padj = 0.05,

                               color = "auto",
                               border_linewidth = 0.2,
                               fill = "auto",
                               n.colorsteps = "auto",
                               nice.breaks = F,
                               show.limits = T,
                               legend.decimals = 1,
                               dotplot = F,
                               plot.feature.breaks = T,
                               plot.sec.axis = F,
                               flip.axes = F,
                               legend.labels = c("min", "", "max"),
                               feature.labels = NULL,
                               feature.labels.nudge_x = -0.1,
                               feature.labels.axis.width = 0.2,
                               plot_hlines_between_groups = F,
                               hlines = NULL,
                               hlines_args = list(),
                               convert_gene_identifier_args = list(ident_in = "SYMBOL", ident_out = "GENENAME"),
                               theme = ggplot2::theme_classic(),
                               theme_args = list(axis.title = ggplot2::element_blank(),
                                                 axis.text.x = ggplot2::element_text(),
                                                 axis.text.y = ggplot2::element_text(size = 10, face = "italic"),
                                                 legend.position = "right",
                                                 legend.direction = "vertical"),
                               legend_fill_args = list(label.theme = ggplot2::element_text(size = 10),
                                                       title.theme = ggplot2::element_text(size = 10),
                                                       title.position = "top",
                                                       title = "..auto..",
                                                       title.hjust = 0.5,
                                                       barwidth = 1,
                                                       barheight = 8,
                                                       order = 1),
                               legend_size_args = list(label.theme = ggplot2::element_text(size = 10),
                                                       title.theme = ggplot2::element_text(size = 10),
                                                       title.position = "top",
                                                       title = "transcription\nfrequency [%]",
                                                       title.hjust = 0.5,
                                                       label.position = "bottom",
                                                       order = 2,
                                                       ncol = NULL,
                                                       nrow = NULL,
                                                       override.aes = list(color = "black")),
                               ...) {

  # ... arguments to ggrepel, like nudge_y

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

  if (methods::is(normalization, "character") && normalization == "scale" && "title" %in% names(legend_fill_args) && legend_fill_args[["title"]] == "..auto..") {
    legend_fill_args[["title"]] <- "transcription\nlevel [z-score]"
  } else if ("legend.title" %in% names(legend_fill_args) && legend_fill_args[["title"]] == "..auto..") {
    legend_fill_args[["title"]] <- "transcription\nlevel"
  }

  #legend.direction <- match.arg(legend.direction, c("horizontal", "vertical"))
  if ("legend.direction" %in% names(theme_args)) {
    if (theme_args$legend.direction == "horizontal") {
      temp <- legend.barwidth
      legend.barwidth <- legend.barheight
      legend.barheight <- temp
    }
  }


  SO <- .check.SO(SO = SO, assay = assay) #length = 1
  assay <- Seurat::DefaultAssay(SO[[1]])

  if (methods::is(normalization, "numeric")) {
    if (length(normalization) != 2) {
      stop("Please set normalization to 'scale' or a numeric vector of length 2 (e.g. c(-1,1)).")
    }
    if (length(legend.labels) > 3) {
      stop("max length of legend.labels is 3. e.g. c('min','int','max') or c('min','','max')")
    }
    if (length(legend.labels) == 2) {
      legend.labels <- c(legend.labels[1], "", legend.labels[2])
    }
  } else if (!is.null(normalization)) {
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
    if (any(unlist(mapply(x = SO, y = meta.col, function(x,y) !y %in% names(x@meta.data), SIMPLIFY = F)))) {
      stop("One of meta.col not found in respective SO.")
    }
    if (any(unlist(mapply(x = SO, y = meta.col, function(x,y) is.numeric(x@meta.data[,y,drop=T]), SIMPLIFY = F)))) {
      stop("One of meta.col across SO is numeric. Please make it factor or character.")
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

  # subset levels for calculation
  SO <- purrr::pmap(list(x = SO, y = levels.calc, z = meta.col), function(x,y,z) {
    if (length(y) < length(unique(x@meta.data[,z,drop=T]))) {
      x <- subset(x, cells = rownames(x@meta.data[,z,drop=F][which(x@meta.data[,z] %in% y),,drop=F]))
    }
    return(x)
  })

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

  # filter non-expressed features with rowSums
  # filter group wise for min_pct with tapply; umi_mat in apply has to be coverted to dense - too expensive
'  group_vec <- unlist(purrr::pmap(list(x = SO, y = meta.col), function(x,y) x@meta.data[,y,drop=T]), use.names = F)
  umi_mat <- do.call(cbind, purrr::map(SO, Seurat::GetAssayData, slot = "data", assay = assay))
  temp_tapply_fun <- function(vec, group_vec) {tapply(vec, group_vec, function(x) sum(x>0)/length(x))}
  out <- purrr::map(split(1:ncol(umi_mat), ceiling(seq_along(1:ncol(umi_mat))/1000)), function(inds) {
    tt <- apply(umi_mat[inds,], MARGIN = 1, FUN = temp_tapply_fun, group_vec = group_vec)
    presto_feat <- names(which(apply(tt, MARGIN = 2, function(x) any(x > min.pct))))
    return(presto_feat)
  }, .progress = T)'


  presto_feat <- unique(unlist(purrr::map(SO, function(x) names(which(Matrix::rowSums(Seurat::GetAssayData(x, slot = "data", assay = assay)) > 0)))))
  if (!is.null(features)) {
    features <- .check.features(SO = SO, features = unique(features), meta.data = F)
    if (any(!features %in% presto_feat)) {
      message("No expressers found for: ", paste(features[which(!features %in% presto_feat)], collapse = ","), ". Will not be plotted.")
    }
    presto_feat <- intersect(presto_feat, features)
  }

  # by default, presto gives deviating results with respect to avgExpr and logFC
  # use expm1 to match results Seurats procedure, see example below (comparison of presto and Seurat)
  wil_auc_raw <- presto::wilcoxauc(X = Gmisc::fastDoCall(cbind, lapply(SO, function(x) expm1(Seurat::GetAssayData(x, slot = "data", assay = assay)[presto_feat,]))),
                                   y = unlist(purrr::pmap(list(x = SO, y = meta.col), function(x,y) x@meta.data[,y,drop=T]), use.names = F))

  # prep features variable here to have one common pipeline below
  if (is.null(features)) {
    order_features <- T
    features <- wil_auc_raw$feature
  }

  ## when no features are provided (features = NULL), then order_features is always TRUE.
  ## when !is.null(features), order_features may be TRUE or FALSE

  if (order_features) {
    features3 <- feature_order_fun(wil_auc_raw = wil_auc_raw,
                                   topn.features = topn.features,
                                   break.ties = break.ties,
                                   topn.metric = topn.metric,
                                   feature_selection_strategy = feature_selection_strategy,
                                   features = features,
                                   min.pct = min.pct,
                                   max.padj = max.padj,
                                   levels.plot = levels.plot)
    if (is.null(hlines)) {
      hlines <- cumsum(rle(as.character(features3$group))[["lengths"]]) + 0.5
      hlines <- hlines[-length(hlines)]
    }
    features <- features3$feature
  }

  ## comparison of average expression by presto and Seurat:
  #tt <- scexpr:::.get.data(SO1, feature = c("GPX3", "cluster_names")) %>% dplyr::filter(cluster_names == "0")
  #out1 <- mean(tt$GPX3) # this is how presto does it
  #out2 <- mean(expm1(tt$GPX3)) # this is how Seurat does it

  # old procedure, but now we use presto with expm1
  #wil_auc_mat <- do.call(cbind, mapply(x = SO, y = meta.col, function(x,y) Seurat::AverageExpression(x, assays = assay, group.by = y, slot = "data", verbose = F)[[1]][features,,drop=F], SIMPLIFY = F))

  ## subset wil_auc by features
  wil_auc <- dplyr::filter(wil_auc_raw, feature %in% features)

  if (!is.null(normalization)) {
    'wil_auc_mat <-
      wil_auc %>%
      dplyr::select(feature, group, avgExpr) %>%
      tidyr::pivot_wider(names_from = group, values_from = avgExpr) %>%
      tibble::column_to_rownames("feature") %>%
      as.matrix()'

    if (methods::is(normalization, "character") && normalization == "scale") {
      wil_auc <-
        wil_auc %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(avgExpr = as.vector(scale(avgExpr)))
    } else if (methods::is(normalization, "numeric")) {
      wil_auc <-
        wil_auc %>%
        dplyr::group_by(feature) %>%
        dplyr::mutate(avgExpr = scale_min_max(avgExpr, min = min(normalization), max = max(normalization)))
    }
    wil_auc <- dplyr::ungroup(wil_auc)
  }

  wil_auc <-
    wil_auc %>%
    dplyr::filter(group %in% unlist(levels.plot)) %>%
    dplyr::mutate(group = factor(group, levels = unlist(levels.plot))) %>%
    dplyr::mutate(feature = factor(feature, levels = unique(features)))

  scale.max <- as.numeric(format(floor_any(max(wil_auc$avgExpr), 0.1), nsmall = 1))
  scale.min <- as.numeric(format(ceiling_any(min(wil_auc$avgExpr), 0.1), nsmall = 1))
  #scale.mid <- as.numeric(format(round(scale.min + ((scale.max - scale.min) / 2), 1), nsmall = 1))
  scale.mid <- 0

  if (methods::is(normalization, "numeric")) {
    labels <- c("min", "int", "max")
  } else {
    labels <- c(scale.min, scale.mid, scale.max)
  }

  if (color == "auto") {
    if (dotplot || nlevels(wil_auc$feature) > 100) {
      color <- "NA"
    } else {
      color <- "grey70"
    }
  }

  ## calculate n.colorsteps
  ## somehow this all assumes that z-scores are provided
  n_neg_bin <- abs(floor(min(wil_auc$avgExpr)))
  n_pos_bin1 <- floor(max(wil_auc$avgExpr))
  n_pos_bin2 <- ceiling(max(wil_auc$avgExpr))
  min_bin <- min(n_neg_bin, n_pos_bin1)
  if (n.colorsteps == "auto") {
    if (n_pos_bin2 == 1 && n_neg_bin == 1) {
      # fine steps for heatmap with low range
      n.colorsteps <- seq(-1, 1, 0.5)
    } else {
      # this may be subject to optimization
      #n.colorsteps <- seq(-min_bin,min_bin,1)
      n.colorsteps <- seq(-n_neg_bin,n_pos_bin1,1)
    }
  }

  if (fill == "auto") {
    fill <- scexpr::col_pal(name = "RColorBrewer::RdBu", n = 11, direction = -1)[-c(1,11)]
    # increase color scale length here with *10
    fill <- prismatic::color(grDevices::colorRampPalette(fill)(length(n.colorsteps)*10))

    # somehow, for scale_fill_stepsn a color scale with one blue color only
    # and 4 red yield an ugly color scale. we have to interpolate the color scale by factor then
    # but then pick n colors * 5 depend upond if bins of n.colorsteps are negative or positive

    # condition on n.colorsteps needed here?
    #if (length(n.colorsteps) > 1 || n.colorsteps == "auto") {
    if (n_neg_bin != n_pos_bin1) {
      if (n_neg_bin > n_pos_bin1) {
        # more negative bins then positive: remove colors from top of the scale
        rm_top_bins <- abs(c(n_pos_bin1-n_neg_bin))
        rm_top_bins <- (length(fill)-(rm_top_bins-1)):length(fill)
        fill <- fill[-rm_top_bins]
      } else if (n_neg_bin < n_pos_bin1) {
        # more positive bins then negative: remove colors from bottom of the scale
        # here multiply by 5 and accordingly to color scale bins, pick a respective number from bottom and top
        top_bins <- (length(fill)-(n_pos_bin2*5-1)):length(fill)
        fill <- fill[c(1:(n_neg_bin*5), top_bins)]
      }
    }
  }

  if (flip.axes) {
    if ("axis.text.x" %in% names(theme_args) && "axis.text.y" %in% names(theme_args)) {
      names(theme_args)[which(names(theme_args) == "axis.text.x")] <- "axis.text.x2"
      names(theme_args)[which(names(theme_args) == "axis.text.y")] <- "axis.text.x"
      names(theme_args)[which(names(theme_args) == "axis.text.x2")] <- "axis.text.y"
    } else if ("axis.text.x" %in% names(theme_args)) {
      names(theme_args)[which(names(theme_args) == "axis.text.x")] <- "axis.text.y"
    } else if ("axis.text.y" %in% names(theme_args)) {
      names(theme_args)[which(names(theme_args) == "axis.text.y")] <- "axis.text.x"
    }
  }

  if (flip.axes) {
    heatmap.plot <- ggplot2::ggplot(wil_auc, ggplot2::aes(y = group, x = feature, fill = avgExpr))
  } else {
    heatmap.plot <- ggplot2::ggplot(wil_auc, ggplot2::aes(x = group, y = feature, fill = avgExpr))
  }

  heatmap.plot <-
    heatmap.plot +
    theme +
    Gmisc::fastDoCall(ggplot2::theme, args = theme_args)

  if (dotplot) {
    heatmap.plot <- heatmap.plot + ggplot2::geom_point(ggplot2::aes(size = pct_in), shape = 21, color = color, stroke = border_linewidth)
  } else {
    heatmap.plot <- heatmap.plot + ggplot2::geom_tile(color = color, linewidth = border_linewidth)
  }

  if (!is.null(n.colorsteps)) {
    if (length(n.colorsteps) == 1 && n.colorsteps != "auto") {
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
    ggplot2::guides(fill = Gmisc::fastDoCall(guide_fun, args = legend_fill_args),
                    size = Gmisc::fastDoCall(ggplot2::guide_legend, args = legend_size_args))

  if (plot.feature.breaks & !is.null(feature.labels)) {
    axis.df <- data.frame(y = 1:length(levels(wil_auc$feature)), feature = levels(wil_auc$feature))
    axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = feature)) +
      ggrepel::geom_text_repel(fontface = "italic", data = axis.df[which(axis.df$feature %in% feature.labels),], ggplot2::aes(label = feature), nudge_x = feature.labels.nudge_x, direction = "y", ...) +
      ggplot2::scale_x_continuous(limits = c(-0.1, 0), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0, length(levels(wil_auc$feature)) + 0.5), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::theme_void()
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
    heatmap.plot <- cowplot::plot_grid(axis, heatmap.plot, align = "h", axis = "tb", nrow = 1, rel_widths = c(feature.labels.axis.width,1)) ## check how to replace with patchwork
  }

  if (plot.sec.axis) {
    out <- Gmisc::fastDoCall(scexpr::convert_gene_identifier, args = c(list(idents = levels(wil_auc$feature)), convert_gene_identifier_args))
    axis.df <- data.frame(y = 1:length(levels(wil_auc$feature)), feature = levels(wil_auc$feature)) %>% dplyr::left_join(out, by = c("feature" = "SYMBOL"))
    axis <- ggplot2::ggplot(axis.df, ggplot2::aes(x = 0, y = y, label = GENENAME)) +
      ggplot2::geom_text(hjust = 0) +
      ggplot2::scale_x_continuous(limits = c(0, 2), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
      ggplot2::scale_y_continuous(limits = c(0.4, length(levels(wil_auc$feature)) + 0.4), expand = c(0, 0), breaks = NULL, labels = NULL, name = NULL) +
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


  if (missing(features) && is.null(feature.labels) && nlevels(wil_auc$feature) > 200 && plot.feature.breaks) {
    print("Large number of rows without selecting feature.labels. Consider setting plot.feature.breaks to FALSE or selecting feature.labels to plot.")
  }

  if (!plot.feature.breaks) {
    heatmap.plot <- heatmap.plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }

  return(list(plot = heatmap.plot, data = wil_auc, complete_data = wil_auc_raw, hlines = hlines))
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

feature_order_fun <- function(wil_auc_raw,
                              topn.features,
                              break.ties,
                              topn.metric,
                              feature_selection_strategy,
                              features,
                              min.pct,
                              max.padj,
                              levels.plot) {
  # select multiple columns for ordering, in order to break ties; especially if padj is used as primary element to order
  if (break.ties) {
    all_metrics <- c("padj", "logFC", "auc")
    second_metric <- all_metrics[which(!all_metrics %in% topn.metric)][1]
    third_metric <- all_metrics[which(!all_metrics %in% topn.metric)][2]
  }

  if (feature_selection_strategy == "2") {
    # this is better for when features are provided; yields better order on y-axis
    features2 <-
      wil_auc_raw %>%
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
      wil_auc_raw %>%
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
  return(features3)
}




