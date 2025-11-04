#' Check relationship (correlation) of clusters from one or two Seurat objects
#'
#' Check https://github.com/skinnider/dismay and https://github.com/tpq/propr
#' propr could not be made for 2 matrices yet.
#'
#' @param objs one Seurat or a named list of 2 Seurats
#' @param meta_cols one or two columns from meta.data
#' @param split split correlation estimation by another variable like
#' donor: column from meta data; must exist in both objs with shared levels
#' when objs is a list of two
#' @param features which features (genes) to use for correlation: all intersecting,
#' intersecting pca features, or a custom vector
#' @param assay assay to pull expression matrix from
#' @param layer layer to pull expression matrix from
#' @param method correlation metric: from psych pkg: pearson, spearman, kendall;
#' from dismay
#' @param heatmap_long_df_args arguments to fcexpr::heatmap_long_df for plotting
#' results
#' @param heatmap_ordering_args arguments to fcexpr::heatmap_ordering
#' @param lower_tri for a symmetric matrix when only one obj is provided:
#' plot the lower triangle only? (remove redundancy)
#' @param min_cells min cell number in both groups to calc correlation
#' @param avg_expression_args args to scexpr::avg_expression
#'
#' @returns
#' @export
#'
#' @examples
cluster_correlation2 <- function(objs,
                                 meta_cols,
                                 split = NULL,
                                 features = c("all", "pca"), # pairwise DEG?
                                 assay = "RNA",
                                 layer = "data",
                                 method = c("pearson","spearman", "kendall", "kendall_zi"),
                                 heatmap_long_df_args = list(
                                   fill = colrr::col_pal("spectral", direction = -1),
                                   colorsteps = "..auto.."
                                 ),
                                 heatmap_ordering_args = list(
                                   feature_order = "hclust",
                                   group_order = "hclust"),
                                 lower_tri = F,
                                 min_cells = 1,
                                 avg_expression_args = list(fun = Matrix::rowMeans,
                                                            fun2 = base::identity)) {

  # so1 <- readRDS("/Volumes/CMS_SSD_2TB/R_scRNAseq/2025_EnghardKlocke_DominoTx/data/SO_processed/SO_DTX_tissue_SCT_harmony_1_3000_25_251020_162905.rds")
  # so2 <- readRDS("/Volumes/CMS_SSD_2TB/R_scRNAseq/2025_Muto_GSE151302_PMID33850129/data/SO_processed/SO_healthy_scnuc_Muto_SCT_harmony_1_3000_25_251001_195449.rds")
  #
  # objs <- list(so1, so2)
  # meta_cols <- c("cluster", "celltype_muto")
  # split = "orig.ident"
  # features = c("all", "pca")
  # assay = "RNA"
  # layer = "data"
  # method <- "pearson"
  #
  # objs <- so1
  # meta_cols <- "cluster"

  if (!requireNamespace("colrr", quietly = T)) {
    devtools::install_github("Close-your-eyes/colrr")
  }
  if (!requireNamespace("fcexpr", quietly = T)) {
    devtools::install_github("Close-your-eyes/fcexpr")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  c(objs,
    meta_cols,
    split,
    assay,
    features,
    split_intersect,
    method) %<-% checks(
      objs = objs,
      meta_cols = meta_cols,
      split = split,
      assay = assay,
      features = features,
      method)

  avg_expr <- purrr::map2(
    objs,
    meta_cols,
    ~Gmisc::fastDoCall(what = avg_expression,
                       args = c(list(obj = .x,
                                     group = .y,
                                     split = split,
                                     assay = assay,
                                     layer = layer,
                                     features = features),
                                avg_expression_args)

    )
  )

  if (!is.null(split)) {
    # cell numbers
    # one more column when !is.null(split)

    ncellcol <- paste0("n_cells_", names(objs))
    cells_split <- purrr::map(c(1,2), function(x) {
      count_cells(obj = objs[[x]],
                  groups = stats::setNames(purrr::compact(c(meta_cols[x], split)),
                                           c(names(objs)[x], split)),
                  colname = paste0("n_cells_", names(objs)[x]))
    })

    # iterate over split_intersect, non-intersecting split-level are ignored
    if (method %in% c("pearson","spearman", "kendall")) {
      corrobj <- purrr::map(split_intersect,
                            ~psych::corr.test(x = avg_expr[[1]][[.x]],
                                              y = avg_expr[[2]][[.x]],
                                              method = method))
    }

    if (method == "kendall_zi") {
      corrobj <- list()
      corrobj[["r"]] <- purrr::map(split_intersect,
                                   ~brathering::kendall_zi_cross(x = avg_expr[[1]][[.x]],
                                                                 y = avg_expr[[2]][[.x]],
                                                                 mc.cores = 8))
    }
    corr_mats <- purrr::map(corrobj, `[[`, "r")
    # average corr coeffs
    # make dfs for joining. this guarantees that respective corr coeffs are grouped together, missing levels in either split level become irrelevant
    newcol <- paste0("r_", method)
    # calc average
    # wide version (does not permit to add cell counts easily)
    # corr_dfs <- purrr::map(names(corr_mats), ~brathering::mat_to_df_long(x = corr_mats[[.x]],
    #                                                                      rownames_to = names(objs)[1],
    #                                                                      colnames_to = names(objs)[2],
    #                                                                      values_to = .x))
    # corr_df <- purrr::reduce(corr_dfs, dplyr::left_join, by = names(objs)[c(1,2)])
    # corr_df[[newcol]] <- purrr::map_dbl(asplit(corr_df[,-c(1,2)], 1), brathering::combine_corrcoeff)
    # long version
    corr_df <- purrr::map_dfr(corr_mats,
                              ~brathering::mat_to_df_long(x = .x,
                                                          rownames_to = names(objs)[1],
                                                          colnames_to = names(objs)[2],
                                                          values_to = method),
                              .id = split) |>
      # add min_cells here
      dplyr::left_join(cells_split[[1]], by = c(names(objs)[1], split)) |>
      dplyr::left_join(cells_split[[2]], by = c(names(objs)[2], split)) |>
      dplyr::mutate(!!method := ifelse(!!rlang::sym(ncellcol[1]) < min_cells, NA, !!rlang::sym(method))) |>
      dplyr::mutate(!!method := ifelse(!!rlang::sym(ncellcol[2]) < min_cells, NA, !!rlang::sym(method))) |>
      dplyr::mutate(!!newcol := brathering::combine_corrcoeff(!!rlang::sym(method)),
                    .by = c(!!rlang::sym(names(objs)[1]), !!rlang::sym(names(objs)[2]))) |>
      tidyr::pivot_wider(names_from = !!rlang::sym(split), values_from = c(!!rlang::sym(method), !!rlang::sym(ncellcol1), !!rlang::sym(ncellcol2)))

    # back to matrix
    corr_mat <- brathering::df_long_to_mat(
      df = corr_df,
      to_rows = names(objs)[1],
      to_cols = names(objs)[2],
      values = newcol
    )

  } else {
    # no split
    corrobj <- psych::corr.test(
      avg_expr[[1]],
      y = avg_expr[[2]],
      method = method
    )
    corr_mat <- corrobj[["r"]]
  }

  # cell count w/o split
  ncellcol <- paste0("n_cells_", names(objs))
  cells <- purrr::map(c(1,2), function(x) {
    count_cells(obj = objs[[x]],
                groups = stats::setNames(purrr::compact(c(meta_cols[x])),
                                         c(names(objs)[x])),
                colname = paste0("n_cells_", names(objs)[x]))
  })

  rdf <- brathering::mat_to_df_long(corr_mat,
                                    rownames_to = names(objs)[1],
                                    colnames_to = names(objs)[2],
                                    values_to = method) |>
    dplyr::left_join(cells[[1]], by = names(objs)[1]) |>
    dplyr::left_join(cells[[2]], by = names(objs)[2])

  # ordering before NA may be introduced
  rdf <- Gmisc::fastDoCall(what = fcexpr::heatmap_ordering,
                           args = c(list(df = rdf,
                                         groups = names(rdf)[1],
                                         features = names(rdf)[2],
                                         values = method),
                                    heatmap_ordering_args))

  if (lower_tri) {
    # if (ncol(corr_mat) != nrow(corr_mat)) {
    #   message("Correlation matrix is not quadratic. Returning the lower triangle may not yield intended results.")
    # }
    rdfmat <- brathering::df_long_to_mat(rdf, to_rows = names(objs)[1], to_cols = names(objs)[2], values = method)
    rdfmat[which(!lower.tri(rdfmat))] <- NA
    rdf <- brathering::mat_to_df_long(rdfmat,
                                      rownames_to = names(objs)[1],
                                      colnames_to = names(objs)[2],
                                      values_to = method) |>
      dplyr::left_join(cells[[1]], by = names(objs)[1]) |>
      dplyr::left_join(cells[[2]], by = names(objs)[2])
    #dplyr::mutate(!!rlang::sym(names(objs)[1]) := factor(!!rlang::sym(names(objs)[1]), levels = levels(rdf[[names(objs)[1]]]))) |>
    #dplyr::mutate(!!rlang::sym(names(objs)[2]) := factor(!!rlang::sym(names(objs)[2]), levels = levels(rdf[[names(objs)[2]]])))
  }

  # set corr value to NA when min_cells is not met by both groups
  # option for min_cells in both vs. in at least one
  rdf_plot <- rdf |>
    dplyr::mutate(!!rlang::sym(method) := ifelse(!!rlang::sym(ncellcol[1]) < min_cells, NA, !!rlang::sym(method))) |>
    dplyr::mutate(!!rlang::sym(method) := ifelse(!!rlang::sym(ncellcol[2]) < min_cells, NA, !!rlang::sym(method)))

  plot <- Gmisc::fastDoCall(what = fcexpr::heatmap_long_df,
                            args = c(list(df = rdf_plot,
                                          groups = names(rdf)[1],
                                          features = names(rdf)[2],
                                          values = method,
                                          heatmap_ordering_args = list(feature_order = "none",
                                                                       group_order = "none")),
                                     heatmap_long_df_args))

  ret <- list(corrobj = corrobj, plot = plot, corr_df_plot = rdf_plot, corr_df = rdf)
  if (!is.null(split)) {
    ret <- c(list(corr_df_split = corr_df), ret)
  }
  return(ret)
}

checks <- function(objs, meta_cols, split, assay, features, method) {

  if (methods::is(objs, "Seurat")) {
    objs <- list(objs, objs)
  }
  objs <- check.SO(objs, assay = assay)
  if (length(meta_cols) == 1) {
    meta_cols <- c(meta_cols, meta_cols)
  }
  if (!meta_cols[1] %in% names(objs[[1]]@meta.data) || !meta_cols[2] %in% names(objs[[2]]@meta.data)) {
    stop("One of meta_cols not found in respective objs.")
  }

  if (!is.list(objs)) {
    stop("objs must be list of seurat objects.")
  }
  if (length(objs) != 2) {
    stop("objs must be of length two.")
  }

  # if (!is.null(split) && method != "pearson") {
  #   message("Averaging correlation values may only be valid for pearson. See https://stats.stackexchange.com/questions/8019/averaging-correlation-values?noredirect=1&lq=1 .")
  # }


  split_intersect <- NULL
  if (!is.null(split)) {
    split <- check.features(objs, features = split, rownames = F)
    split_intersect <- intersect(objs[[1]]@meta.data[[split]], objs[[2]]@meta.data[[split]])
    if (length(split_intersect) == 0) {
      stop("No intersecting levels in ", split, " column found in SOs.")
    }

    if (length(split_intersect) != length(unique(objs[[1]]@meta.data[[split]])) || length(split_intersect) != length(unique(objs[[2]]@meta.data[[split]]))) {
      message("not all levels of split in obj1 and obj2 overlap or intersect.")
      message("Intersecting levels: ", paste(split_intersect, collapse = ", "), ".")
    }
    names(split_intersect) <- split_intersect
  }
  assay <- match.arg(assay, names(objs[[1]]@assays))
  method = match.arg(method, c("pearson","spearman", "kendall", "kendall_zi"))

  if (length(features) %in% c(1,2) && length(intersect(features, c("all", "pca"))) == 2) {
    features <- match.arg(features, choices = c("all", "pca"))
  }

  if (length(features) == 1 && features == "all") {
    features <- Reduce(intersect, purrr::map(objs, rownames))
  } else if (length(features) == 1 && features == "pca") {
    features2 <- names(check.reduction(objs, reduction = "pca"))
    if (features != features2) {
      message("using pca: ", features2)
      features <- features2
    }
    features <- Reduce(intersect,  purrr::map(objs, ~rownames(.x@reductions[[features]]@feature.loadings)))
  } else {
    # features <- check.features(objs, features, meta.data = F)
  }
  return(list(objs, meta_cols, split, assay, features, split_intersect, method))
}

count_cells <- function(obj, groups, colname = "n") {

  out <- dplyr::summarise(obj@meta.data, n = dplyr::n(), .by = unname(groups))
  if (!is.null(names(groups))) {
    names(out)[1:(ncol(out)-1)] <- names(groups)
  }
  names(out)[ncol(out)] <- colname
  return(out)
}
