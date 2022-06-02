#' Find groups (clusters) of cells (single cell transcriptomes) with highest gene expression correlation in two Seurat objects
#'
#' Sometimes one wants to avoid combining data sets into a single seurat object as this often requires to employ an integration procedure or a
#' batch correction and respective uncertainties. At some point then one may be interested in which clusters correspond to
#' each other in separate seurat object. This function quickly returns a graphic and data frames to judge best matching clusters.
#' Seurat::AverageExpression is used to derive average gene expressions per cluster. Selected gene (features) are then used
#' for correlation calculation.
#'
#' @param SO named list of exactly 2 Seurat objects
#' @param assay which assay to obtain expression values from
#' @param meta.cols character vector of length 2, indicating the column names of meta.data in SO[[1]] and SO[[2]] to use for comparison,
#' e.g. the clustering columns in both SOs
#' @param levels list of length 2 indicating the factor levels in meta.cols to include; this also defines axis orders;
#' can be incomplete e.g. if order matters only for some levels, if complement.levels = T the missing levels are added randomly
#' @param complement.levels add all missing levels in random order
#' @param features features to use for correlation calculation; will always be reduced to intersecting features between SOs;
#' if all, all intersecting features in assay are considered; if pca, intersecting rownames of feature.loadings in pca-reduction;
#' or a vector of features
#' @param avg.expr avg.expr from a previous return list of cluster_correlation_matrix; e.g. provide that to  only change axis
#' orders or similar but avoid repeated calculation
#' @param method which correlation metric to calculate, passed to stats::cor
#' @param corr.in.percent display correlation in percent (TRUE, T) or as a fraction (FALSE, F)
#' @param round.corr decimal places to round correlation values to
#' @param split.by split correlation calculation by a common column with common levels in SOs;
#' individual correlation coefficients are averaged after fishers z-transformation (atanh) and then
#' transformed back (tanh)
#' @param min.cells minimum number of cells per group in meta.cols (and if provided split.by)
#' to calculate and return a correlation coefficient
#'
#' @return list with (i) ggplot object of correlation matrix plot, (ii) the data frame to that plot and (iii) calculated average expressions from Seurat::AverageExpression
#' @export
#'
#' @examples
cluster_correlation_matrix <- function(SO,
                                       assay = c("RNA", "SCT"),
                                       meta.cols,
                                       levels = NULL, # NA possible
                                       complement.levels = F,
                                       features = c("all", "pca"),
                                       split.by = NULL,
                                       avg.expr,
                                       method = c("pearson", "kendall", "spearman"),
                                       corr.in.percent = FALSE,
                                       min.cells = 20,
                                       round.corr = 2) {

  if (!requireNamespace("reshape2", quietly = T)) {
    utils::install.packages("reshape2")
  }

  if (length(features) == 1 || length(features) == 2 && length(intersect(features, c("all", "pca"))) == 2) {
    features <- match.arg(features, choices = c("all", "pca"))
  }

  if (missing(meta.cols)) {
    stop("Please provide meta.cols for the first and second SO.")
  } else {
    if (length(meta.cols) != 2) {
      stop("meta.cols has to be of length 2.")
    }
  }
  if (!missing(avg.expr)) {
    if (!is.list(avg.expr)) {
      stop("avg.expr has to be a list of matrices.")
    }
    if (length(avg.expr) != 2) {
      stop("avg.expr has to be of length 2.")
    }
  }

  if (!is.null(levels)) {
    if (!is.list(levels)) {
      stop("levels has to be a list.")
    }
    if (length(levels) != 2) {
      stop("levels has to be of length 2.")
    }
    levels <- lapply(seq_along(SO), function(x) .check.levels(SO[[x]], meta.cols[x], levels = levels[[x]], append_by_missing = complement.levels))

  } else {
    levels <- lapply(seq_along(SO), function(x) .check.levels(SO[[x]], meta.cols[x], levels = levels, append_by_missing = F))
  }

  SO <- .check.SO(SO, assay = assay, length = 2)
  ## names for SO are very important below. missing names will cause errors.
  assay <- match.arg(assay, c("RNA", "SCT"))
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  meta.cols[1] <- .check.features(SO[[1]], features = meta.cols[1], rownames = F)
  meta.cols[2] <- .check.features(SO[[2]], features = meta.cols[2], rownames = F)

  if (length(features) == 1 && features == "all") {
    features <- Reduce(intersect, lapply(SO, function(x) rownames(x)))
  } else if (length(features) == 1 && features == "pca") {
    features <- Reduce(intersect, lapply(SO, function(x) rownames(x@reductions[["pca"]]@feature.loadings)))
  } else {
    features <- .check.features(SO, features, meta.data = F)
  }

  ### add option to split comparison by another covariate and then use the average correlation coeff:
  # https://www.researchgate.net/post/average_of_Pearson_correlation_coefficient_values
  # https://en.wikipedia.org/wiki/Fisher_transformation
  # https://stats.stackexchange.com/questions/8019/averaging-correlation-values

  #split.by<-"Pat"
  if (!is.null(split.by)) {
    split.by <- .check.features(SO, features = split.by, rownames = F)
    intersect_levels <- intersect(SO[[1]]@meta.data[,split.by], SO[[2]]@meta.data[,split.by])
    if (length(intersect_levels) == 0) {
      stop("No intersecting levels in ", split.by, " column found in SOs.")
    } else {
      message("Intersecting levels between SOs for meta column ", split.by, ": ", paste(intersect_levels, collapse = ", "), ".")
    }
    # intersect_levels for correct order
    SO[[1]] <- Seurat::SplitObject(SO[[1]], split.by = split.by)[intersect_levels]
    SO[[2]] <- Seurat::SplitObject(SO[[2]], split.by = split.by)[intersect_levels]
  } else {
    SO[[1]] <- list(SO[[1]])
    SO[[2]] <- list(SO[[2]])
  }


  ## create a table of number of cells per clusters
  fff <- function(so, meta.col) {
    z <- utils::stack(table(so@meta.data[,meta.col]))
    names(z) <- c("n_cells", "meta_col_level")
    z[,2] <- as.character(z[,2])
    return(z)
  }
  ## use purrr to loop through SOs and check how many observations (cells) per cluster are found
  n_cell_df <- purrr::map2(.x = SO, .y = meta.cols, .f = ~ purrr::map2(.x = .x, .y = .y, .f = ~ fff(.x,.y)))
  n_cell_df <- purrr::map(n_cell_df, dplyr::bind_rows, .id = split.by)
  # filter for relevant levels only
  n_cell_df <- purrr::map2(.x = n_cell_df, .y = levels, .f = ~ .x[which(.x$meta_col_level %in% .y),])
  n_cell_df <- dplyr::bind_rows(n_cell_df, .id = "SO")

  if (nrow(n_cell_df[which(n_cell_df$n_cells < min.cells),]) > 0) {
    message("These elements are excluded from correlation analysis due to low n_cells below ", min.cells, ".")
    print(n_cell_df[which(n_cell_df$n_cells < min.cells),])
  }

'    SO_blood@meta.data <-
    SO_blood@meta.data %>%
    tidyr::separate(orig.ident, into = c("Pat", "rep", "body_fluid"), remove = F)
  SO_urine@meta.data <-
    SO_urine@meta.data %>%
    tidyr::separate(orig.ident, into = c("Pat", "rep", "body_fluid"), remove = F)
SO <- list(SO_blood, SO_urine)

SO_blood_split <- Seurat::SplitObject(SO_blood, split.by = "orig.ident")
SO_urine_split <- Seurat::SplitObject(SO_urine, split.by = "orig.ident")
'

  ## works
  #avg.expr2 <- purrr::map2(.x = SO, .y = meta.cols, .f = ~ purrr::map2(.x = .x, .y = .y, .f = ~ Seurat::AverageExpression(.x, assays = assay, group.by = .y, slot = "data", verbose = F)[[assay]]))
  #identical(Seurat::AverageExpression(SO[[1]][[1]], assays = assay, group.by = meta.cols[[1]], slot = "data", verbose = F)[[assay]], avg.expr2[[1]][[1]])

  # for pearson correlation, also a linear correlation would be feasable
  # conduct linear model in case of pearson
  # https://sebastianraschka.com/faq/docs/pearson-r-vs-linear-regr.html

  if (missing(avg.expr)) {
    avg.expr <- purrr::map2(.x = SO, .y = meta.cols, .f = ~ purrr::map2(.x = .x, .y = .y, .f = ~ Seurat::AverageExpression(.x, assays = assay, features = features, group.by = .y, slot = "data", verbose = F)[[assay]]))
    # avg.expr <- lapply(seq_along(SO), function(x) Seurat::AverageExpression(SO[[x]], assays = assay, features = features, group.by = meta.cols[[x]], slot = "data", verbose = F)[[assay]])
  } else {
    avg.expr <- purrr::map(.x = avg.expr, .f = ~ purrr::map(.x = .x, .f = ~ .x[features,,drop=F]))
    #avg.expr <- lapply(avg.expr, function(x) x[features,,drop=F])
  }

  # filter for levels
  avg.expr <- purrr::map2(.x = avg.expr, .y = levels, function(x,y) {purrr::map(.x = x, .f = ~ .x[,which(colnames(.x) %in% y),drop=F])})

  if (!is.null(split.by)) {
    n_cell_df_nest <-
      n_cell_df %>%
      dplyr::filter(n_cells >= min.cells) %>%
      dplyr::group_by(SO, !!rlang::sym(split.by)) %>%
      dplyr::summarise(meta_col_levels = list(meta_col_level), .groups = "drop") %>%
      as.data.frame()

    ## filter clusters with n cells below threshold
    for (i in 1:nrow(n_cell_df_nest)) {
      avg.expr[[n_cell_df_nest[i,"SO"]]][[n_cell_df_nest[i,split.by]]] <-
        avg.expr[[n_cell_df_nest[i,"SO"]]][[n_cell_df_nest[i,split.by]]][,n_cell_df_nest[[i,"meta_col_levels"]]]
    }

    # calculate correlations
    cm <- purrr::map2(.x = avg.expr[[1]], .y = avg.expr[[2]], function(x,y) stats::cor(x = x, y = y, method = method))

    combs <- as.data.frame(tidyr::expand_grid(x = levels[[1]], y = levels[[2]]))
    cm.melt <- dplyr::bind_rows(lapply(split(combs, 1:nrow(combs)), function(comb) {
      corrs <- unlist(sapply(cm, function(z) {
        tryCatch(z[comb[["x"]], comb[["y"]]], error = function (e) NULL)
      }))
      if (!is.null(corrs)) {
        data.frame(Var1 = comb[["x"]], Var2 = comb[["y"]], value = tanh(mean(atanh(corrs))))
      } else {
        NULL
      }
    }))

  } else {
    avg.expr <- purrr::flatten(avg.expr)

    n_cell_df_nest <-
      n_cell_df %>%
      dplyr::filter(n_cells >= min.cells) %>%
      dplyr::group_by(SO) %>%
      dplyr::summarise(meta_col_levels = list(meta_col_level), .groups = "drop") %>%
      as.data.frame()

    ## filter clusters with n cells below threshold
    for (i in 1:nrow(n_cell_df_nest)) {
      avg.expr[[n_cell_df_nest[i,"SO"]]] <-
        avg.expr[[n_cell_df_nest[i,"SO"]]][,n_cell_df_nest[[i,"meta_col_levels"]]]
    }

    # calculate correlations
    cm <- stats::cor(x = avg.expr[[1]], y = avg.expr[[2]], method = method)
    cm.melt <- reshape2::melt(cm)
  }

  cm.melt$Var1 <- factor(as.character(cm.melt$Var1), levels = levels[[1]])
  cm.melt$Var2 <- factor(as.character(cm.melt$Var2), levels = levels[[2]])

  pp <- c(col_pal("RdBu", n = 11, reverse = T)[2:5], rep(col_pal("RdBu", n = 11, reverse = T)[6], 6), col_pal("RdBu", n = 11, reverse = T)[7:10])

  cm.plot <- ggplot2::ggplot(cm.melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(colour = "black") +
    #ggplot2::scale_fill_gradient2(high = "#BC3F2B", low = "#4C6CA6", mid = "#edf9ff", midpoint = median(cm.melt$value)) +
    ggplot2::scale_fill_gradientn(colors = pp) +
    ggplot2::labs(x = paste0(names(SO)[1], " ", meta.cols[1]), y = paste0(names(SO)[2], " ", meta.cols[2])) +
    ggplot2::theme_classic()

  if (corr.in.percent) {
    cm.plot <- cm.plot + ggplot2::geom_text(size = cm.plot[["theme"]][["text"]][["size"]] *(5/14), ggplot2::aes(label = paste0(round(value*100, 0), " %")))
  } else {
    cm.plot <- cm.plot + ggplot2::geom_text(size = cm.plot[["theme"]][["text"]][["size"]] *(5/14), ggplot2::aes(label = format(round(value, round.corr), nsmall = round.corr)))
  }

  # dot plot of genes
  '  avg.expr <- lapply(avg.expr, function(x) reshape2::melt(as.matrix(x)))
  names(avg.expr[[1]]) <- c("Gene.symbol", "cluster.x", "avg.expr.x")
  names(avg.expr[[2]]) <- c("Gene.symbol", "cluster.y", "avg.expr.y")

  dot.plot.matrix.df <- left_join(avg.expr[[1]], avg.expr[[2]], by = "Gene.symbol")

  dot.plot.matrix <- ggplot(dot.plot.matrix.df, aes(x = scales::squish(avg.expr.x, c(0.001,100000)), y = scales::squish(avg.expr.y, c(0.001,100000)))) +
    geom_point(size = 0.2) +
    theme_bw() +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_log10() +
    scale_y_log10() +
    facet_grid(rows = vars(cluster.x), cols = vars(cluster.y), scales = "free")'

  return(list(cm.plot = cm.plot, cm.melt = cm.melt, avg.expr = avg.expr))
}
