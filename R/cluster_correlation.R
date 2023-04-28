#' Find groups (clusters) of cells (single cell transcriptomes) with highest gene expression correlation in two Seurat objects
#'
#' Sometimes one wants to avoid combining data sets into a single seurat object as this often requires to employ an integration procedure or a
#' batch correction and respective uncertainties. At some point then one may be interested in which clusters correspond to
#' each other in separate seurat object. This function quickly returns a graphic and data frames to judge best matching clusters.
#' Seurat::AverageExpression is used to derive average gene expressions per cluster. Selected gene (features) are then used
#' for correlation calculation.
#'
#' @param SO named list of exactly 2 Seurat objects; you may also pass the same object twice to have an intra-comparison
#' @param meta.cols character vector of length 2, indicating the column names of meta.data in SO[[1]] and SO[[2]] to use for comparison,
#' e.g. the clustering columns in both SOs
#' @param features features to use for correlation calculation; will always be reduced to intersecting features between SOs;
#' if all, all intersecting features in assay are considered; if pca, intersecting rownames of feature.loadings in pca-reduction;
#' or a vector of features
#' @param assay which assay to obtain expression values from
#' @param levels list of length 2 indicating the factor levels in meta.cols to include; this also defines axis orders;
#' can be incomplete e.g. if order matters only for some levels, if complement.levels = T the missing levels are added randomly
#' @param complement.levels add all missing levels in random order
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
#' @param lower.triangle.only return the lower triangle of correlation matrix which will contain
#' unique values only
#' @param aspect.ratio aspect ratio of matrix plot
#'
#' @return list with (i) ggplot object of correlation matrix plot, (ii) the data frame to that plot and (iii) calculated average expressions from Seurat::AverageExpression
#' @export
#'
#' @examples
cluster_correlation <- function(SO,
                                meta.cols,
                                features = c("all", "pca"),
                                assay = c("RNA", "SCT"),
                                levels = NULL, # NA possible
                                complement.levels = F,
                                split.by = NULL,
                                avg.expr,
                                method = c("pearson","spearman", "kendall"),
                                corr.in.percent = FALSE,
                                min.cells = 20,
                                round.corr = 2,
                                lower.triangle.only = F,
                                aspect.ratio = 1,
                                mid.white.strech.length = 4) {

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
  if (!meta.cols[1] %in% names(SO[[1]]@meta.data) || !meta.cols[2] %in% names(SO[[2]]@meta.data)) {
    stop("One of meta.cols not found in SOs.")
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
  if (!is.null(split.by) && method != "pearson") {
    warning("Averaging correlation values may only be valid for pearson. See https://stats.stackexchange.com/questions/8019/averaging-correlation-values?noredirect=1&lq=1 .")
  }
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

  if (nrow(n_cell_df[which(n_cell_df[,"n_cells",drop=T] < min.cells),]) > 0) {
    message("These elements are excluded from correlation analysis due to low n_cells below ", min.cells, ".")
    print(n_cell_df[which(n_cell_df[,"n_cells",drop=T] < min.cells),])
  }

  ## works
  #avg.expr2 <- purrr::map2(.x = SO, .y = meta.cols, .f = ~ purrr::map2(.x = .x, .y = .y, .f = ~ Seurat::AverageExpression(.x, assays = assay, group.by = .y, slot = "data", verbose = F)[[assay]]))
  #identical(Seurat::AverageExpression(SO[[1]][[1]], assays = assay, group.by = meta.cols[[1]], slot = "data", verbose = F)[[assay]], avg.expr2[[1]][[1]])

  if (missing(avg.expr)) {
    avg.expr <- purrr::map2(.x = SO, .y = meta.cols, .f = ~ purrr::map2(.x = .x, .y = .y, .f = ~ Seurat::AverageExpression(.x, assays = assay, features = features, group.by = .y, slot = "data", verbose = F)[[assay]]))
  } else {
    avg.expr <- purrr::map(.x = avg.expr, .f = ~ purrr::map(.x = .x, .f = ~ .x[features,,drop=F]))
  }

  # filter for levels
  avg.expr <- purrr::map2(.x = avg.expr, .y = levels, function(x,y) {purrr::map(.x = x, .f = ~ .x[,which(colnames(.x) %in% y),drop=F])})


  if (!is.null(split.by)) {

    # currently not implemented for splitted corr analyses
    lms_or_ranks <- NULL

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
    if (method == "pearson") {
      # for pearson correlation, also a linear correlation is feasible
      # this also allows for analysis of residuals
      # https://sebastianraschka.com/faq/docs/pearson-r-vs-linear-regr.html

      # replace mathematical operators and special characters from level names to allow passing formula to lm
      if (any(c(!identical(make.names(colnames(avg.expr[[1]])), colnames(avg.expr[[1]])),
                !identical(make.names(colnames(avg.expr[[2]])), colnames(avg.expr[[2]]))))) {
        message("Mathematical operators of special characters found in factor levels. Those are replace by dots in linear model formulae.")
      }
      # run the make.names thing to have a common command below
      orig_colnames_1 <- stats::setNames(colnames(avg.expr[[1]]), nm = make.names(colnames(avg.expr[[1]])))
      orig_colnames_2 <- stats::setNames(colnames(avg.expr[[2]]), nm = make.names(colnames(avg.expr[[2]])))
      colnames(avg.expr[[1]]) <- make.names(colnames(avg.expr[[1]]))
      colnames(avg.expr[[2]]) <- make.names(colnames(avg.expr[[2]]))

      combs_raw <- expand.grid(colnames(avg.expr[[1]]), colnames(avg.expr[[2]]))
      combs_raw <- combs_raw[which(as.character(combs_raw[,1]) != as.character(combs_raw[,2])),]
      #orders <- apply(combs_raw, 1, function(x) order(c(x[1], x[2])), simplify = F)
      #combs <- split(combs_raw, 1:nrow(combs_raw))
      #combs_helper <- dplyr::bind_rows(mapply(x = combs, y = orders, function(x,y) x[y], SIMPLIFY = F))
      #names(combs) <- apply(combs_raw, 1, function(x) paste0(x[1], "___", x[2]))

      ## use reformulate to have the correct var names in lms
      lms_or_ranks <- purrr::map2(.x = as.character(combs_raw$Var1), .y = as.character(combs_raw$Var2), .f = function(x,y) {
        stats::lm(formula = stats::reformulate(x, y), data = data.frame(scale(avg.expr[[1]][,x,drop=F]), y = scale(avg.expr[[2]][,y,drop=F])))
      })

      names(lms_or_ranks) <- paste0(orig_colnames_2[sapply(purrr::map(purrr::map(lms_or_ranks, purrr::chuck, "model"), colnames), "[", 1)], "___",
                                    orig_colnames_1[sapply(purrr::map(purrr::map(lms_or_ranks, purrr::chuck, "model"), colnames), "[", 2)])

      # matrix was found to be correct, like when stats::cor produces it
      cm <- purrr::map(lms_or_ranks, purrr::chuck, "coefficients", 2)
      cm <- matrix(unlist(cm), nrow = ncol(avg.expr[[1]]))
      rownames(cm) <- orig_colnames_1[colnames(avg.expr[[1]])]
      colnames(cm) <- orig_colnames_2[colnames(avg.expr[[2]])]
    } else {
      cm <- stats::cor(x = avg.expr[[1]], y = avg.expr[[2]], method = method)
    }

    # order columns and rows to get the right elements filtered in case of lower.triangle.only
    cm <- cm[levels[[1]], levels[[2]]]

    if (lower.triangle.only) {
      if (ncol(cm) != nrow(cm)) {
        message("Correlation matrix is not quadratic. Returning the lower triangle may not yield intended results.")
      }
      cm[which(!lower.tri(cm))] <- NA
    }

    cm.melt <- reshape2::melt(cm)
    cm.melt <- cm.melt[which(!is.na(cm.melt$value)),]

    # switch Var1 and Var2 to make it consistent with pearson case above (lm) where the formula is y~x (y in front stemming from avg.expr[[2]])
    # names(cm.melt)[1:2] <- c("Var2", "Var1")
    # why? or transpose cm?

    if (method != "pearson") {
      lms_or_ranks <- lapply(1:nrow(cm.melt), function(x) {
        ranks <- data.frame(avg.expr[[1]][,as.character(cm.melt[x, "Var1"])], avg.expr[[2]][,as.character(cm.melt[x, "Var2"])])
        names(ranks) <- c(paste0(names(SO)[1], "___", cm.melt[x, "Var1"]), paste0(names(SO)[2], "___", cm.melt[x, "Var2"]))
        name1 <- paste0(names(SO)[1], "___", cm.melt[x, "Var1"], "_", "rank")
        name2 <- paste0(names(SO)[2], "___", cm.melt[x, "Var2"], "_", "rank")

        ranks <-
          ranks %>%
          dplyr::mutate({{name1}} := dplyr::dense_rank(!!rlang::sym(paste0(names(SO)[1], "___", as.character(cm.melt[x, "Var1"])))),
                        {{name2}} := dplyr::dense_rank(!!rlang::sym(paste0(names(SO)[2], "___", as.character(cm.melt[x, "Var2"]))))) %>%
          dplyr::mutate(rank_diff = !!rlang::sym(name1) - !!rlang::sym(name2)) %>%
          dplyr::mutate(abs_rank_diff = abs(rank_diff)) %>%
          dplyr::arrange(rank_diff)

        names(ranks)[3:4] <- paste0(c(cm.melt[x, "Var1"], cm.melt[x, "Var2"]), "_rank")
        return(ranks)
      })
      names(lms_or_ranks) <- paste0(cm.melt[1:nrow(cm.melt),"Var1"], "___", cm.melt[1:nrow(cm.melt),"Var2"])
    }
  }

  cm.melt$Var1 <- factor(as.character(cm.melt$Var1), levels = levels[[1]])
  cm.melt$Var2 <- factor(as.character(cm.melt$Var2), levels = levels[[2]])

  pp <- c(col_pal("RdBu", n = 11, reverse = T)[2:5],
          rep(col_pal("RdBu", n = 11, reverse = T)[6], mid.white.strech.length),
          col_pal("RdBu", n = 11, reverse = T)[7:10])

  cm.plot <- ggplot2::ggplot(cm.melt, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(colour = "black") +
    #ggplot2::scale_fill_gradient2(high = "#BC3F2B", low = "#4C6CA6", mid = "#edf9ff", midpoint = median(cm.melt$value)) +
    ggplot2::scale_fill_gradientn(colors = pp) +
    ggplot2::labs(x = paste0(names(SO)[1], " ", meta.cols[1]), y = paste0(names(SO)[2], " ", meta.cols[2])) +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = aspect.ratio)

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

  if (method == "pearson") {
    return(list(plot = cm.plot, data = cm.melt, avg.expr = avg.expr, lms = lms_or_ranks))
  } else {
    return(list(plot = cm.plot, data = cm.melt, avg.expr = avg.expr, ranks = lms_or_ranks))
  }

}
