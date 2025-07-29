#' Get data frame for plotting from Seurat object(s)
#'
#' @param SO one Seurat object or a list of multiple ones
#' @param feature vector of features to fetch (genes or column names in
#' meta data)
#' @param reduction reduction to fetch from reductions slot of SO
#' @param dims dimensions of the selected reduction to extract from SO
#' @param assay which assay to get expression data from
#' @param layer which layer/slot to get expression data from, e.g. counts or data
#' @param cells vector of cell names to select for regular plotting (=1),
#' deselected ones are plotted with color col_ex_cells (=0) in feature_plot
#' @param meta_col meta data columns to attach
#' @param qmin lower quantile of feature values where to cut the color scale
#' @param qmax upper quantile
#' @param order order rows from lowest to highest (negative at bottom)
#' @param order_rev reverse row order
#' @param order_abs absolute values for ordering: zeros first and anything
#' distant from zero on top
#' @param order_discr set TRUE to for order by abundance, only relevant for
#' discrete meta features
#' @param shuffle do shuffle rows (if order is FALSE) to randomize plotting
#' order; if !order and !shuffle the provided order is not altered
#' @param bury_NA put NA to the bottom which will bury them upon plotting
#' @param na_rm remove NA values? (only for meta feature)
#' @param inf_rm remove infinite values?
#' @param label_feature feature for plotting labels
#' @param contour_feature feature for plotting contours
#' @param split_feature feature for creating facets
#' @param shape_feature feature to shape plotted points by
#' @param trajectory_slot
#' @param downsample random downsample cells which will completely remove them
#' from returned data frame and hence influence statistics and plotting
#' @param feature_ex filter cells by an exclusion feature, cells with UMI>0
#' are excluded and get col_ex_cells as color, passed to scexpr::get_data
#' @param feature_cut filter cells by a one or more cutoff features,
#' cells with UMI>feature_cut_expr are excluded and get col_ex_cells
#' @param feature_cut_expr numeric vector of cutoff expressions
#'
#' @return a data frame
#' @export
#'
#' @examples
get_data <- function(SO,
                     feature,
                     reduction = "umap",
                     dims = c(1,2),
                     assay = "RNA",
                     layer = "data",
                     cells = NULL,
                     meta_col = NULL,
                     qmin = 0,
                     qmax = 1,
                     order = T,
                     order_rev = F,
                     order_abs = T,
                     order_discr = T,
                     shuffle = F,
                     bury_NA = T,
                     na_rm = F,
                     inf_rm = F,
                     label_feature = NULL,
                     contour_feature = NULL,
                     split_feature = NULL,
                     shape_feature = NULL,
                     trajectory_slot = NULL,
                     downsample = 1,
                     feature_cut = NULL,
                     feature_cut_expr = 0,
                     feature_ex = NULL) {

  if (missing(SO)) {stop("Seurat object missing.")}
  if (length(dims) != 2 || !methods::is(dims, "numeric")) {stop("dims has to be a numeric vector of length 2, e.g. c(1,2).")}

  # also check layer
  SO <- check.SO(
    SO = SO,
    assay = assay,
    meta.col = meta_col)
  assay <- Seurat::DefaultAssay(SO[[1]])
  cells <- check.and.get.cells(SO = SO,
                               assay = assay,
                               cells = cells,
                               feature_cut = feature_cut,
                               feature_cut_expr = feature_cut_expr,
                               feature_ex = feature_ex,
                               downsample = downsample)


  feature <- check.features(SO = SO, feature = unique(feature))
  label_feature <- check.features(SO = SO, features = label_feature, rownames = F)
  contour_feature <- check.features(SO = SO, features = contour_feature, rownames = F)
  split_feature <- check.features(SO = SO, features = split_feature, rownames = F)
  shape_feature <- check.features(SO = SO, features = shape_feature, rownames = F)
  reduction <- check.reduction(SO = SO, reduction = reduction, dims = dims)
  assay <- match.arg(assay, names(SO[[1]]@assays))


  if (qmax > 1) {
    #message("qmax and qmin are divided by 100. Please provide values between 0 and 1.")
    qmax <- qmax/100
    qmin <- qmin/100
  }

  all_gene_features <- Reduce(intersect, lapply(SO, rownames))
  all_meta_features <- Reduce(intersect, lapply(SO, function(x) names(x@meta.data)))
  missing_features <- feature[which(!feature %in% all_gene_features & !feature %in% all_meta_features)]
  double_features <- feature[which(feature %in% all_gene_features & feature %in% all_meta_features)]
  gene_features <- feature[which(feature %in% all_gene_features)]
  meta_features <- feature[which(feature %in% all_meta_features)]

  if (length(missing_features) > 0) {
    message("Feature(s) ", paste(missing_features, collapse = ", "), " not found.")
    feature <- setdiff(feature, missing_features)
  }
  if (length(double_features) > 0) {
    message("Feature(s) ", paste(double_features, collapse = ", "), " found in expression and meta data. This is not handled at the moment and filtered. Consider renaming the column of meta data.")
    feature <- setdiff(feature, double_features)
  }
  if (length(feature) == 0) {
    stop("No feature left.")
  }

  data <- purrr::map_dfr(SO, function(x) {

    if (utils::compareVersion(as.character(x@version), "4.9.9") == 1) {
      GetAssayData_args <- list(object = x,
                                layer = layer,
                                assay = assay)
    } else {
      GetAssayData_args <- list(object = x,
                                slot = layer,
                                assay = assay)
    }

    data <- cbind(data.frame(t(as.matrix(Gmisc::fastDoCall(
      what = Seurat::GetAssayData,
      args = GetAssayData_args
    )[gene_features,,drop = F])), check.names = F),
    data.frame(
      x@meta.data[,meta_features,drop = F],
      stringsAsFactors = F,
      check.names = F
    ))

    if (!is.null(reduction)) {
      reduction <- unique(unlist(lapply(reduction, function(z) {
        names(x@reductions)[which.min(utils::adist(z, names(x@reductions), ignore.case = T))]
      })))
      for (i in reduction) {
        data <- cbind(data, Seurat::Embeddings(x, reduction = i))
      }
    }


    # redundant to meta_features
    if (!is.null(meta_col)) {
      data <- cbind(data, x@meta.data[,meta_col,drop=F])
    }
    if (!is.null(label_feature)) {
      data <- cbind(data, data.frame("label_feature" = x@meta.data[[label_feature]], stringsAsFactors = F, check.names = F))
    }
    if (!is.null(contour_feature)) {
      data <- cbind(data, data.frame("contour_feature" = x@meta.data[[contour_feature]], stringsAsFactors = F, check.names = F))
    }
    if (!is.null(shape_feature)) {
      data <- cbind(data, data.frame("shape_feature" = x@meta.data[[shape_feature]], stringsAsFactors = F, check.names = F))
    }
    if (!is.null(split_feature)) {
      data <- cbind(data, data.frame("split_feature" = x@meta.data[[split_feature]], stringsAsFactors = F, check.names = F))
    } else {
      data[["split_feature"]] <- factor("1")
    }
    return(data)
  }, .id = "SO.split")

  # ensure that facet ordering is according to the order of SO objects provided
  data$SO.split <- factor(data$SO.split, levels = names(SO))
  data <- dplyr::relocate(data, SO.split, .after = dplyr::last_col())

  if (!is.null(cells)) {
    # this avoids introduction of NA for duplicate names
    data <- data[which(rownames(data) %in% names(cells)),]
    data[["cells"]] <- unname(cells[rownames(data)])
  } else {
    data[["cells"]] <- 1
  }
  if (nrow(data) == 0) {
    stop("None of cells found.")
  }

  # check for expressers
  for (i in gene_features) {
    if (all(data[,i,drop = T] == 0)) {
      message("No expressers found for ", i, ".")
    }
  }

  nacol <- apply(data[,meta_features,drop = F], 2, anyNA)
  if (any(nacol)) {
    message(paste(names(which(nacol)), collapse = ", "), ": NA found in data.")
  }

  infcol <- apply(data[,meta_features,drop = F], 2, function(x) any(is.infinite(x)))
  if (any(infcol)) {
    message(paste(names(which(infcol)), collapse = ", "), ": Inf found in data.")
  }

  data <- tibble::rownames_to_column(data, "id")
  # rm all feature but one
  # split data into list with one feature each
  data <- purrr::map(stats::setNames(lapply(seq_along(feature), function(i) feature[-i]), feature),
                     ~dplyr::select(data, -dplyr::all_of(.x)))
  #data <- purrr::map2(.x = data, .y = feature, ~tidyr::pivot_longer(.x, cols = dplyr::all_of(.y), names_to = "feat"))
  data <- purrr::map2(.x = data, .y = feature, ~dplyr::rename(.x, "feature" = !!rlang::sym(.y)))

  # does not work if features are of different type
  #data <- tidyr::pivot_longer(data, cols = dplyr::all_of(feature), names_to = "feat")
  #data <- split(data, data$feat)

  ### data now is a list for each feature

  if (na_rm) {
    data <- purrr::map(data, function(x) x[which(!is.na(x[["feature"]])),])
  }
  if (inf_rm) {
    data <- purrr::map(data, function(x) x[which(!is.infinite(x[["feature"]])),])
  }

  # use squishing to dampen extreme values - this will produce actually wrong limits on the legend
  if (qmin > 0 || qmax < 1) {
    data <- purrr::map(data, function(x) {
      if (is.numeric(x[["feature"]])) {
        if (all(x[["feature"]] >= 0)) { # > 0 or >= 0 ?!
          # expression is always greater than 0 and non-expresser are excluded
          inds <- which(x[["feature"]] > 0)
          x[["feature"]][inds] <- scales::squish(
            x[["feature"]][inds],
            range = c(stats::quantile(x[["feature"]][inds], qmin),
                      stats::quantile(x[["feature"]][inds], qmax))
          )
        } else {
          # e.g. for module scores below 0
          x[["feature"]] <- scales::squish(x[["feature"]], range = c(stats::quantile(x[["feature"]], qmin), stats::quantile(x[["feature"]], qmax)))
        }
      }
      return(x)
    })
  }

  ## params:
  # order = T --> ordering by values (lowest are negative values, then 0, then positive ones)
  # order_abs T --> absolute values far away from zero are plotted on top
  # shuffle = T --> only considered when !order; will shuffle the data data.frame; if !order and !shuffle a custom order can be provided from outside
  # order_rev = T --> reverse the order so that lowest values (or zeros if order_abs = T) are on top

  data <- purrr::map(data, function(x) {
    # gene feat or numeric meta_feat
    if (is.numeric(x[["feature"]])) {
      if (order) {
        # per default order will put NAs on top
        # abs: in case negative values are contained in meta_col, any extreme away from 0 will be plotted on top
        if (order_abs) {
          x <- x[order(abs(x[["feature"]]), decreasing = order_rev, na.last = !bury_NA),]
        } else {
          x <- x[order(x[["feature"]], decreasing = order_rev, na.last = !bury_NA),]
        }
      } else if (shuffle) {
        x <- x[sample(1:nrow(x)),]
      }
    } else if (!is.numeric(x[["feature"]]) && is.logical(order_discr) && order_discr) {
      ## this is only for meta features
      # put most frequent groups to back (plot first)
      # NA is not considered by split()
      # replace it by a character value just for splitting, undo this afterwards
      #x[["feature"]] <- factor(x[["feature"]], exclude = c())

      if (anyNA(x[["feature"]])) {
        na_replace <- "NA"
        while(na_replace %in% unique(x[["feature"]])) {
          na_replace <- paste(c(na_replace, na_replace), collapse = "_")
        }
        level_order <- NULL
        if (is.factor(x[["feature"]])) {
          level_order <- levels(x[["feature"]])
          x[["feature"]] <- as.character(x[["feature"]])
        }
        x[which(is.na(x[["feature"]])),1] <- na_replace
        x <- dplyr::bind_rows(split(x, x[["feature"]])[names(sort(table(x[["feature"]]), decreasing = T))])
        x[which(x[["feature"]] == na_replace),1] <- NA
        if (bury_NA) {
          x <- rbind(x[which(is.na(x[["feature"]])),], x[which(!is.na(x[["feature"]])),])
        }
        if (!is.null(level_order)) {
          x[["feature"]] <- factor(x[["feature"]], levels = level_order)
        }
      } else {
        x <- dplyr::bind_rows(split(x, x[["feature"]])[names(sort(table(x[["feature"]]), decreasing = T))])
      }
    } else if (!is.numeric(x[["feature"]]) && is.logical(order_discr) && shuffle) {
      x <- x[sample(1:nrow(x)),]
    }

    return(x)
  })

  for (i in seq_along(data)) {
    if (names(data)[i] %in% gene_features) {
      attr(data[[i]], "feature_type") <- "gene"
    } else {
      attr(data[[i]], "feature_type") <- "meta"
    }
    attr(data[[i]], "feature") <- feature[i]
    attr(data[[i]], "layer") <- layer
    attr(data[[i]], "qmin") <- qmin
    attr(data[[i]], "qmax") <- qmax
    attr(data[[i]], "dim1") <- paste0(reduction, "_", dims[1])
    attr(data[[i]], "dim2") <- paste0(reduction, "_", dims[2])
    if (!is.null(label_feature)) {
      attr(data[[i]], "label_feature") <- label_feature
    }
    if (!is.null(contour_feature)) {
      attr(data[[i]], "contour_feature") <- contour_feature
    }
    if (!is.null(split_feature)) {
      attr(data[[i]], "split_feature") <- split_feature
    }
    if (!is.null(shape_feature)) {
      attr(data[[i]], "shape_feature") <- shape_feature
    }
    attr(data[[i]], "feature_ex") <- feature_ex
    attr(data[[i]], "feature_cut_expr") <- feature_cut_expr
    attr(data[[i]], "feature_cut") <- feature_cut
  }

  # data <- dplyr::bind_rows(data)
  #
  # if (!is.null(trajectory_slot)) {
  #   data_traj <- purrr::map_dfr(SO, function(x) {
  #     # rbind will throw error if column names to not match
  #     data_traj <- NULL
  #     if (trajectory_slot %in% names(Seurat::Misc(x))) {
  #       data_traj <- Seurat::Misc(x, trajectory_slot)[["df"]]
  #     } else {
  #       message("Trajectory slot not found in Seurat::Misc.")
  #     }
  #     return(data_traj)
  #   }, .id = "SO.split")
  #   return(list(data = data, data_traj = data_traj))
  # }

  return(data)
}


