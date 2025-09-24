#' Get data frame for plotting from Seurat object(s)
#'
#' From a wide data frame with any kind of expression values,
#' prepare a list of data frames ready to pass to feature_plot_data.
#'
#' @param feature columns that are features to plot as color scale
#' @param reduction columns that are dimension reduction
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
#'
#' @return a list of data frames
#' @export
#'
#' @examples
prep_data_df_wide <- function(data,
                              feature = NULL,
                              reduction = c("umap_1", "umap_2"),
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
                              shape_feature = NULL) {


  if (qmax > 1) {
    #message("qmax and qmin are divided by 100. Please provide values between 0 and 1.")
    qmax <- qmax/100
    qmin <- qmin/100
  }

  data <- as.data.frame(data)

  reduction <- grep(paste(reduction, collapse = "|"), names(data), value = T, ignore.case = T)
  if (any(!reduction %in% names(data))) {
    stop("At least one of reduction not found.")
  }
  if (length(reduction) > 2) {
    stop("length of reduction > 2: ", paste(reduction, collapse = ", "))
  }

  if (is.null(feature)) {
    feature <- setdiff(names(data), reduction)
  }

  data <- tibble::rownames_to_column(data, "id")
  data[["SO.split"]] <- 1
  data[["SO.split"]] <- factor(data[["SO.split"]], levels = 1)
  data[["cells"]] <- 1 # option to provide from outside?

  if (!is.null(label_feature)) {
    if (!label_feature %in% names(data)) {
      message("label_feature not found.")
    } else {
      data[["label_feature"]] <- as.character(data[[label_feature]])
    }
  }
  if (!is.null(contour_feature)) {
    if (!contour_feature %in% names(data)) {
      message("contour_feature not found.")
    } else {
      data[["contour_feature"]] <- as.character(data[[contour_feature]])
    }
  }
  if (!is.null(shape_feature)) {
    if (!shape_feature %in% names(data)) {
      message("shape_feature not found.")
    } else {
      data[["shape_feature"]] <- as.character(data[[shape_feature]])
    }
  }
  if (!is.null(split_feature)) {
    if (!split_feature %in% names(data)) {
      message("split_feature not found.")
    } else {
      data[["split_feature"]] <- as.character(data[[split_feature]])
    }
  } else {
    data[["split_feature"]] <- factor("1")
  }


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
        x[["feature"]] <- scales::squish(x[["feature"]], range = c(stats::quantile(x[["feature"]], qmin), stats::quantile(x[["feature"]], qmax)))
        
        ## dont do this here as non scrnaseq data are provided
        # if (all(x[["feature"]] >= 0)) { # > 0 or >= 0 ?!
        #   # expression is always greater than 0 and non-expresser are excluded
        #   inds <- which(x[["feature"]] > 0)
        #   x[["feature"]][inds] <- scales::squish(
        #     x[["feature"]][inds],
        #     range = c(stats::quantile(x[["feature"]][inds], qmin),
        #               stats::quantile(x[["feature"]][inds], qmax))
        #   )
        # } else {
        #   # e.g. for module scores below 0
        #   x[["feature"]] <- scales::squish(x[["feature"]], range = c(stats::quantile(x[["feature"]], qmin), stats::quantile(x[["feature"]], qmax)))
        # }
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

    if (is.numeric(data[[i]][["feature"]])) {
      attr(data[[i]], "feature_type") <- "gene"
    } else {
      attr(data[[i]], "feature_type") <- "meta"
    }

    attr(data[[i]], "feature") <- feature[i]
    attr(data[[i]], "layer") <- "data"
    attr(data[[i]], "qmin") <- qmin
    attr(data[[i]], "qmax") <- qmax
    attr(data[[i]], "dim1") <- reduction[1]
    attr(data[[i]], "dim2") <- reduction[2]
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
    attr(data[[i]], "feature_ex") <- NULL
    attr(data[[i]], "feature_cut_expr") <- 0
    attr(data[[i]], "feature_cut") <- NULL
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


