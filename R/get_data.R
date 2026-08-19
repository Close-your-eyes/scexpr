#' Prepare plotting data from one or more Seurat objects
#'
#' Extract expression values, cell-level metadata, and two-dimensional
#' embeddings from one or more Seurat objects. The result is split by requested
#' feature so that each element can be passed directly to plotting functions.
#'
#' @details
#' Gene features are read from `layer` in `assay`; metadata features are read
#' from each object's `meta.data` slot. A requested feature must be present in
#' every object. Features that are missing, or that occur both as a gene and as
#' a metadata column, are reported and omitted.
#'
#' For numeric features, `qmin` and `qmax` winsorize values for visualization.
#' For non-negative features, quantiles are calculated from positive values
#' only. This changes the returned plotting values and therefore should not be
#' used when the result will be used for downstream statistics.
#'
#' Row order determines plotting order: rows appearing later are typically
#' drawn on top. Ordering and shuffling are applied independently to each
#' requested feature.
#'
#' @param SO A Seurat object or a named list of Seurat objects. When multiple
#'   objects are supplied, each requested feature must be available in all of
#'   them.
#' @param feature A character vector of gene names or columns in `meta.data` to
#'   extract.
#' @param reduction A reduction name, or one reduction name per object. Names
#'   are matched case-insensitively; compatible reductions in multiple objects
#'   may be renamed internally to a common label. Set to `NULL` to omit
#'   embeddings.
#'   An object's `misc` slot can override the requested reduction through, in
#'   precedence order, `reduction_preferred`, `preferred_reduction`, and
#'   `reduction`. Write the preferred reduction into one of these slots.
#' @param dims A numeric vector of length two identifying the embedding
#'   dimensions used by downstream plotting code. Defaults to `c(1, 2)`.
#' @param assay Name of the assay from which gene expression is extracted.
#' @param layer Name of the assay layer (or slot), such as `"counts"` or
#'   `"data"`.
#' @param cells Optional cell selection. Cell names identify rows to retain;
#'   associated values are copied to the returned `cells` column, where
#'   plotting code conventionally uses `1` for selected and `0` for excluded
#'   cells.
#' @param meta_col Optional character vector of additional `meta.data` columns
#'   to append to every result.
#' @param qmin,qmax Lower and upper quantiles used to winsorize numeric feature
#'   values. Supply proportions between 0 and 1. If `qmax > 1`, both values are
#'   interpreted as percentages and divided by 100.
#' @param order Logical; if `TRUE`, order rows by feature value before
#'   returning them.
#' @param order_rev Logical; reverse the selected numeric ordering.
#' @param order_abs Logical; when ordering numeric features, order by absolute
#'   value so that values farthest from zero are drawn together. Ignored when
#'   `order = FALSE`.
#' @param order_discr Logical; for discrete features, group rows by decreasing
#'   category abundance. The most abundant groups are placed first and are
#'   therefore typically drawn behind smaller groups.
#' @param shuffle Logical; randomly shuffle rows when numeric ordering is
#'   disabled. For discrete features, shuffling is used only when
#'   `order_discr = FALSE`.
#' @param bury_NA Logical; place missing feature values first so that they are
#'   typically hidden beneath non-missing points.
#' @param na_rm Logical; remove rows whose requested feature is `NA`.
#' @param inf_rm Logical; remove rows whose requested feature is infinite.
#' @param label_feature Optional metadata column used for plot labels. It is
#'   returned as `label_feature`.
#' @param contour_feature Optional metadata column used for contours. It is
#'   returned as `contour_feature`.
#' @param split_feature Optional metadata column used for facets. It is
#'   returned as `split_feature`; when omitted, that column is a factor with a
#'   single level, `"1"`.
#' @param shape_feature Optional metadata column used to map point shapes. It
#'   is returned as `shape_feature`.
#' @param trajectory_slot Reserved for compatibility; currently unused.
#' @param downsample Cell downsampling specification passed to
#'   `check.and.get.cells()`. Removed cells are absent from the returned data
#'   and from any downstream calculations.
#' @param feature_cut Optional feature or features used to flag cells by an
#'   expression cutoff. Passed to `check.and.get.cells()`.
#' @param feature_cut_expr Numeric cutoff value or vector corresponding to
#'   `feature_cut`.
#' @param feature_ex Optional exclusion feature. Cells with non-zero expression
#'   are flagged by `check.and.get.cells()` for downstream plotting.
#' @param try_df Logical; if `FALSE` (the default), return one data frame per
#'   requested feature. If `TRUE`, attempt to join those frames into a single
#'   wide data frame.
#'
#' @return If `try_df = FALSE`, a named list with one data frame per retained
#'   feature. In each data frame, the requested values are stored in `feature`,
#'   cell identifiers in `id`, object membership in `SO.split`, and selection
#'   status in `cells`; requested embeddings and annotation columns are also
#'   included. Each element carries attributes describing the feature, layer,
#'   quantile limits, embedding dimensions, and optional plotting mappings.
#'   If `try_df = TRUE`, the feature-specific frames are joined into one data
#'   frame and each feature keeps its original column name.
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
#' \dontrun{
#' # Return one plotting data frame per feature.
#' data_list <- get_data(
#'   pbmc,
#'   feature = c("CD3E", "orig.ident"),
#'   reduction = "umap"
#' )
#'
#' # Combine the feature columns into a single wide data frame.
#' data_wide <- get_data(
#'   pbmc,
#'   feature = c("CD3E", "MS4A1"),
#'   qmin = 0.01,
#'   qmax = 0.99,
#'   try_df = TRUE
#' )
#' }
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
                     feature_ex = NULL,
                     try_df = F) {



  if (missing(SO)) {stop("Seurat object missing.")}
  if (length(dims) != 2 || !methods::is(dims, "numeric")) {stop("dims has to be a numeric vector of length 2, e.g. c(1,2).")}

  # also check layer
  SO <- scexpr:::check.SO(
    SO = SO,
    assay = assay,
    meta.col = meta_col)
  assay <- Seurat::DefaultAssay(SO[[1]])
  cells <- scexpr:::check.and.get.cells(SO = SO,
                                        assay = assay,
                                        cells = cells,
                                        feature_cut = feature_cut,
                                        feature_cut_expr = feature_cut_expr,
                                        feature_ex = feature_ex,
                                        downsample = downsample)

  feature <- scexpr:::check.features(SO = SO, features = feature)
  label_feature <- scexpr:::check.features(SO = SO, features = label_feature, rownames = F)
  contour_feature <- scexpr:::check.features(SO = SO, features = contour_feature, rownames = F)
  split_feature <- scexpr:::check.features(SO = SO, features = split_feature, rownames = F)
  shape_feature <- scexpr:::check.features(SO = SO, features = shape_feature, rownames = F)

  c(reduction, SO) %<-% check.reduction(SO = SO, reduction = reduction)
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

  if (length(missing_features)) {
    message("Feature(s) ", paste(missing_features, collapse = ", "), " not found.")
    feature <- setdiff(feature, missing_features)
  }
  if (length(double_features)) {
    message("Feature(s) ", paste(double_features, collapse = ", "), " found in expression and meta data. This is not handled at the moment and filtered. Consider renaming the column of meta data.")
    feature <- setdiff(feature, double_features)
  }
  if (!length(feature)) {
    stop("No feature left.")
  }

  data <- purrr::map_dfr(SO, function(x) {
    data <- cbind(
      get_layer(obj = x,
                assay = assay,
                layer = layer,
                features = gene_features,
                transpose = T,
                as = "df"),
      data.frame(x@meta.data[,meta_features,drop = F],
                 stringsAsFactors = F,
                 check.names = F)
    )

    if (!is.null(reduction)) {
      reduction <- unique(unlist(lapply(names(reduction), function(z) {
        names(x@reductions)[which.min(stringdist::stringdist(tolower(z), tolower(names(x@reductions))))]
      })))

      for (i in reduction) {
        # valid umap colnames can be umap_1 but also UMAP_1
        # is case of multi-SO this could cause incompatibility
        # so always to lower
        red <- Seurat::Embeddings(x, reduction = i)
        colnames(red) <- tolower(colnames(red))
        data <- cbind(data, red)
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
    # do it here to enable duplicate names, e.g. from multiple objects
    # duplicate rownames would be appended by ...1 and prohibit checking for cells below
    data <- tibble::rownames_to_column(data, "id")
    return(data)
  }, .id = "SO.split")


  # ensure that facet ordering is according to the order of SO objects provided
  data$SO.split <- factor(data$SO.split, levels = names(SO))
  data <- dplyr::relocate(data, SO.split, .after = dplyr::last_col())

  if (!is.null(cells)) {
    # this avoids introduction of NA for duplicate names
    data <- data[which(data$id %in% names(cells)),]
    data[["cells"]] <- unname(cells[data$id])
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

  #data <- tibble::rownames_to_column(data, "id")

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

  classes <- purrr::map(data, freeze_classes, cols = "feature")
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

        na_replace <- "_NA_"
        while(na_replace %in% unique(x[["feature"]])) {
          na_replace <- paste(c(na_replace, na_replace), collapse = "_")
        }
        level_order <- NULL
        if (is.factor(x[["feature"]])) {
          level_order <- levels(x[["feature"]])
          x[["feature"]] <- as.character(x[["feature"]])
        }
        x[which(is.na(x[["feature"]])),"feature"] <- na_replace
        x <- dplyr::bind_rows(split(x, x[["feature"]])[names(sort(table(x[["feature"]]), decreasing = T))])
        x[which(x[["feature"]] == na_replace),"feature"] <- NA
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

  data <- purrr::map2(data, classes, ~thaw_cols(.x, .y))

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
    if (!is.null(reduction)) {
      reduction <- tolower(reduction)
      attr(data[[i]], "dim1") <- paste0(reduction, "_", dims[1])
      attr(data[[i]], "dim2") <- paste0(reduction, "_", dims[2])
    }
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

  if (try_df) {
    for(i in names(data)) {
      names(data[[i]])[which(names(data[[i]]) == "feature")] <- i
    }
    data <- purrr::reduce(data, dplyr::left_join)
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


## to brathering
freeze_classes <- function(df, cols) {

  # keep only columns that exist, ignore others silently
  cols <- intersect(cols, names(df))

  # store classes of selected columns
  classes <- sapply(df[cols], class)

  return(classes = classes)
}


thaw_cols <- function(df, classes) {
  restore_class <- function(x, cls) {
    if (cls == "numeric")  return(as.numeric(x))
    if (cls == "integer")  return(as.integer(x))
    if (cls == "logical")  return(as.logical(x))
    if (cls == "factor")   return(as.factor(x))
    if (cls == "Date")     return(as.Date(x))
    if (cls == "POSIXct")  return(as.POSIXct(x))
    x
  }

  # Only restore columns that still exist
  cols_to_restore <- intersect(names(classes), names(df))

  df[cols_to_restore] <- Map(
    restore_class,
    df[cols_to_restore],
    classes[cols_to_restore]
  )

  return(df)
}

get_gene_features <- function(obj) {
  genes <- unique(unlist(purrr::map(SeuratObject::Assays(obj),
                                    ~rownames(get_layer(obj, assay = .x, layer = "counts")))))
  return(genes)
}

check.reduction <- function(SO,
                            reduction = NULL) {

  if (is.null(reduction)) return(list(reduction, SO))
  SO <- if (is.list(SO)) SO else list(SO)
  red_list <- purrr::map(SO, ~names(.x@reductions))
  if (any(!lengths(red_list))) stop("At least one SO has no reduction.")
  reduction <- brathering::recycle(reduction, SO)

  reduction <- purrr::map2_chr(SO, reduction, function(x,y) {
    if ("reduction_preferred" %in% names(x@misc)) {
      if (tolower(x@misc$reduction_preferred[1]) != tolower(y)) {
        message("reduction set to x@misc$reduction_preferred[1]: ", x@misc$reduction_preferred[1])
      }
      y <- x@misc$reduction_preferred[1]
    }
    if ("preferred_reduction" %in% names(x@misc)) {
      if (tolower(x@misc$preferred_reduction[1]) != tolower(y)) {
        message("reduction set to x@misc$preferred_reduction[1]: ", x@misc$preferred_reduction[1])
      }
      y <- x@misc$preferred_reduction[1]
    }
    if ("reduction" %in% names(x@misc)) {
      if (tolower(x@misc$reduction[1]) != tolower(y)) {
        message("reduction set to x@misc$reduction[1]: ", x@misc$reduction[1])
      }
      y <- x@misc$reduction[1]
    }
    return(y)
  })

  if (length(unique(reduction)) == 1) {
    reduction <- unique(reduction)
  }

  if (length(reduction) == 1) {
    common_red <- Reduce(intersect, red_list)
    if (!length(common_red)) {
      # only when multiple SO

      message("No common reduction found in SOs. Will pick best matches and rename to a common label.")
      red_match <- purrr::map(red_list, ~grep(reduction, .x, ignore.case = T, value = T))
      if (any(x <- !lengths(red_match))) {
        for (i in which(x)) {
          red_match[[i]] <- red_list[[i]][which.min(utils::adist(reduction, red_list[[i]], ignore.case = T))]
          # red_match[[i]] <- red_list[[i]][length(red_list[[i]])]
        }
      }
      if (any(x <- lengths(red_match) > 1)) {
        for (i in which(x)) {
          red_match[[i]] <- red_match[[i]][which.min(utils::adist(reduction, red_match[[i]], ignore.case = T))]
          # red_match[[i]] <- red_match[[i]][1]
        }
      }
      message("using: ", paste(unlist(red_match), collapse = ", "))
      # only keep required reduction to avoid name conflict
      SO <- purrr::map2(SO, red_match, ~keep_reduction(.x, .y))
      SO <- purrr::map2(SO, red_match, ~rename_reduction(.x, .y, reduction))
      common_red <- reduction
    }
    red <- find_reduction(query = reduction, available = common_red)

  } else if (length(reduction) > 1) {

    if (length(reduction) != length(SO)) stop("length(reduction) != length(SO)")
    reduction <- purrr::map2_chr(SO, reduction, ~find_reduction(query = .y, available = names(.x@reductions)))
    SO <- purrr::map2(SO, reduction, ~keep_reduction(.x, .y))

    if (all(grepl("umap", reduction))) {
      red <- "umap"
    } else if (all(grepl("tsne", reduction))) {
      red <- "tsne"
    } else {
      red <- "reduction"
    }
    SO <- purrr::map2(SO, reduction, ~rename_reduction(.x, .y, red))
  }

  key <- SO[[1]]@reductions[[red]]@key
  red <- stats::setNames(sub("_$", "", key), nm = red)
  return(list(red, SO))
}

find_reduction <- function(query, available) {
  matches <- grep(query, available, ignore.case = TRUE, value = TRUE)

  if (!length(matches)) {
    matches <- grep("umap", available, ignore.case = TRUE, value = TRUE)
  }
  if (!length(matches)) {
    matches <- grep("tsne", available, ignore.case = TRUE, value = TRUE)
  }
  if (!length(matches)) {
    message("Reduction not found; using closest match.")
    matches <- available[which.min(utils::adist(query, available, ignore.case = TRUE)[1, ])]
  }
  if (length(matches) > 1) {
    matches <- matches[which.min(utils::adist(query, matches, ignore.case = TRUE)[1, ])]
  }

  return(matches)
}

keep_reduction <- function(x, reduction) {
  x@reductions <- x@reductions[reduction]
  return(x)
}
