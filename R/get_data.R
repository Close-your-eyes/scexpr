get_data <- function(SO,
                     feature,
                     reduction = "umap",
                     assay = c("RNA", "SCT"),
                     slot = "data",
                     cells = NULL,
                     split_by = NULL,
                     shape_by = NULL,
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
                     trajectory_slot = NULL) {


  SO <- .check.SO(
    SO = SO,
    assay = assay,
    split.by = split_by,
    shape.by = shape_by,
    meta.col = meta_col,
    cells = cells
  )
  assay <- match.arg(assay, names(SO[[1]]@assays))

  # if (length(feature) > 1) {
  #   # if function is called to get many features
  #   # ordering becomes irrelevant
  #   # also no values should be removed then (NA or Inf)
  #   na_rm <- F
  #   inf_rm <- F
  #   order <- F
  #   order_discr <- F
  #   shuffle <- F
  # }

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

    data <- cbind(data.frame(t(as.matrix(Seurat::GetAssayData(
      object = x,
      slot = slot,
      assay = assay
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

    if (!is.null(meta_col)) {
      data <- cbind(data, x@meta.data[,meta_col,drop=F])
    }
    if (is.null(split_by)) {
      data[,"split_by"] <- "1"
    } else {
      data[,"split_by"] <- x@meta.data[,as.character(split_by)]
    }
    if (!is.null(shape_by)) {
      data[,shape_by] <- as.factor(as.character(x@meta.data[,as.character(shape_by)]))
    }

    if (!is.null(label_feature)) {
      data <- cbind(data, data.frame(label_feature = x@meta.data[,label_feature,drop=T], stringsAsFactors = F, check.names = F))
      if (nrow(unique(data[,which(colnames(data) %in% c(feature, "label_feature"))])) !=
          nrow(unique(data[,which(colnames(data) %in% c(feature)), drop=F]))) {
        stop("label_feature entries do not exactly one feature entry each. This must be the case though.")
      }
    }

    return(data)
  }, .id = "SO.split")

  # ensure that facet ordering is according to the order of SO objects provided
  data$SO.split <- factor(data$SO.split, levels = names(SO))
  data <- dplyr::relocate(data, SO.split, .after = dplyr::last_col())

  if (!is.null(cells)) {
    data <- data[cells,]
  }

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
  data <- tidyr::pivot_longer(data, cols = dplyr::all_of(feature), names_to = "feat")
  data <- split(data, data$feat)

  if (na_rm) {
    data <- purrr::map(data, function(x) x[which(!is.na(x[,"feat"])),])
  }
  if (inf_rm) {
    data <- purrr::map(data, function(x) x[which(!is.infinite(x[,"feat"])),])
  }

  # use squishing to dampen extreme values - this will produce actually wrong limits on the legend
  if (qmin > 0 || qmax < 1) {
    data <- purrr::map(data, function(x) {
      if (is.numeric(x[,"feat"])) {
        if (all(x[,"feat"] >= 0)) { # > 0 or >= 0 ?!
          # expression is always greater than 0 and non-expresser are excluded
          x[,"feat"][which(x[,"feat"] > 0)] <- scales::squish(
            x[,"feat"][which(x[,"feat"] > 0)],
            range = c(stats::quantile(x[,"feat"][which(x[,"feat"] > 0)], qmin),
                      stats::quantile(x[,"feat"][which(x[,"feat"] > 0)], qmax))
          )
        } else {
          # e.g. for module scores below 0
          x[,"feat"] <- scales::squish(x[,"feat"], range = c(stats::quantile(x[,"feat"], qmin), stats::quantile(x[,"feat"], qmax)))
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
    if (is.numeric(x[,"feat"])) {
      if (order) {
        # per default order will put NAs on top
        # abs: in case negative values are contained in meta_col, any extreme away from 0 will be plotted on top
        if (order_abs) {
          x <- x[order(abs(x[,"feat"]), decreasing = order_rev, na.last = !bury_NA),]
        } else {
          x <- x[order(x[,"feat"], decreasing = order_rev, na.last = !bury_NA),]
        }
      } else if (shuffle) {
        x <- x[sample(1:nrow(x)),]
      }
    } else if (!is.numeric(x[,"feat"]) && is.logical(order_discr) && order_discr) {
      ## this is only for meta features
      # put most frequent groups to back (plot first)
      # NA is not considered by split()
      # replace it by a character value just for splitting, undo this afterwards
      #x[,"feat"] <- factor(x[,"feat"], exclude = c())

      if (anyNA(x[,"feat"])) {
        na_replace <- "NA"
        while(na_replace %in% unique(x[,"feat"])) {
          na_replace <- paste(c(na_replace, na_replace), collapse = "_")
        }
        level_order <- NULL
        if (is.factor(x[,"feat"])) {
          level_order <- levels(x[,"feat"])
          x[,"feat"] <- as.character(x[,"feat"])
        }
        x[which(is.na(x[,"feat"])),1] <- na_replace
        x <- dplyr::bind_rows(split(x, x[,"feat"])[names(sort(table(x[,"feat"]), decreasing = T))])
        x[which(x[,"feat"] == na_replace),1] <- NA
        if (bury_NA) {
          x <- rbind(x[which(is.na(x[,"feat"])),], x[which(!is.na(x[,"feat"])),])
        }
        if (!is.null(level_order)) {
          x[,"feat"] <- factor(x[,"feat"], levels = level_order)
        }
      } else {
        x <- dplyr::bind_rows(split(x, x[,"feat"])[names(sort(table(x[,"feat"]), decreasing = T))])
      }
    } else if (!is.numeric(x[,"feat"]) && is.logical(order_discr) && shuffle) {
      x <- x[sample(1:nrow(x)),]
    }

    return(x)
  })

  data <- dplyr::bind_rows(data)

  if (!is.null(trajectory_slot)) {
    data_traj <- purrr::map_dfr(SO, function(x) {
      # rbind will throw error if column names to not match
      data_traj <- NULL
      if (trajectory_slot %in% names(Seurat::Misc(x))) {
        data_traj <- Seurat::Misc(x, trajectory_slot)[["df"]]
      } else {
        message("Trajectory slot not found in Seurat::Misc.")
      }
      return(data_traj)
    }, .id = "SO.split")
    return(list(data = data, data_traj = data_traj))
  }

  return(data)
}


