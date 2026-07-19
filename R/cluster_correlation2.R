#' Correlate average expression profiles between Seurat clusters or groups
#'
#' Computes correlations between groups of cells defined by metadata columns in
#' one Seurat object or between two Seurat objects. For each group, average
#' expression is calculated first, then correlations are computed between group
#' average-expression profiles. Optionally, correlations can be computed within
#' shared levels of a splitting variable and then combined across split levels.
#'
#' @param objs A Seurat object, or a named list of exactly two Seurat objects.
#'   If a single Seurat object is supplied, it is compared against itself.
#' @param meta_cols Character vector of length one or two. Metadata column(s)
#'   used to define the groups/clusters to compare. If length one, the same
#'   column is used for both objects.
#' @param split Optional character scalar. Metadata column used to compute
#'   correlations separately within each shared split level, such as donor,
#'   batch, sample, or condition. The column must exist in both objects and must
#'   have at least one shared level.
#' @param features Character scalar or character vector. If `"pca"`, uses the
#'   intersecting features from the PCA loadings of both objects. If `"all"`,
#'   uses all intersecting features between the two objects. If a character
#'   vector of feature names is supplied, those features are used directly.
#' @param assay Character scalar. Assay from which expression values are read.
#'   The assay must exist in both Seurat objects.
#' @param layer Character scalar. Assay layer used for expression values,
#'   usually `"data"`, `"counts"`, or `"scale.data"`.
#' @param method Character scalar. Correlation method. One of `"spearman"`,
#'   `"pearson"`, `"kendall"`, or `"kendall_zi"`. The first three are passed to
#'   `psych::corr.test()`. `"kendall_zi"` uses
#'   `brathering::kendall_zi_cross()`. kendall_zi not tested well and only
#'   used when split is provided.
#' @param heatmap_long_df_args Named list of arguments passed to
#'   `fcexpr::heatmap_long_df()` when generating the heatmap.
#' @param heatmap_ordering_args Named list of arguments passed to
#'   `fcexpr::heatmap_ordering()` before plotting.
#' @param min_cells Integer. Minimum number of cells required in both compared
#'   groups. Correlations involving groups with fewer cells are set to `NA` in
#'   the plotting data.
#' @param avg_expression_args Named list of additional arguments passed to
#'   `scexpr::avg_expression()`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{corrobj}{Raw correlation object or list of correlation objects.}
#'   \item{plot}{Heatmap plot object.}
#'   \item{corr_df_plot}{Long-format correlation data used for plotting.}
#'   \item{corr_df}{Long-format correlation data before plotting filters.}
#'   \item{corr_df_split}{Only when `split` is supplied; split-level correlation data.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clustcorr <- cluster_correlation2(
#'   objs = so,
#'   meta_cols = "SCT_harmony_snn_res.0.1",
#'   method = "spearman"
#' )
#'
#' clustcorr2 <- cluster_correlation2(
#'   objs = list(reference = so1, query = so2),
#'   meta_cols = c("cluster", "celltype"),
#'   split = "donor",
#'   features = "pca",
#'   min_cells = 10
#' )
#' }
cluster_correlation2 <- function(objs,
                                 meta_cols,
                                 split = NULL,
                                 features = c("pca", "all"), # pairwise DEG?
                                 assay = "RNA",
                                 layer = "data",
                                 method = c("spearman", "pearson", "kendall", "kendall_zi"),
                                 heatmap_long_df_args = list(
                                   fill = colrr::col_pal("spectral", direction = -1),
                                   values_zscored = F,
                                   colorsteps = 10,
                                   lower_tri = F,
                                   theme_args = list()),
                                 heatmap_ordering_args = list(
                                   feature_order = "hclust",
                                   group_order = "hclust"),
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
    pak::pak("Close-your-eyes/colrr")
  }
  if (!requireNamespace("fcexpr", quietly = T)) {
    pak::pak("Close-your-eyes/fcexpr")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    pak::pak("Close-your-eyes/brathering")
  }


  if (is.list(objs) && is.null(names(objs))) {
    nme <- as.character(deparse(substitute(objs)))
    names(objs) <- strsplit(gsub("list\\(|\\)", "", nme), ", ")[[1]]
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
      method = method)

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
    } else if (method == "kendall_zi") {
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
      tidyr::pivot_wider(names_from = !!rlang::sym(split), values_from = c(!!rlang::sym(method), !!rlang::sym(ncellcol[1]), !!rlang::sym(ncellcol[2])))

    # back to matrix
    corr_mat <- brathering::df_long_to_mat(
      df = corr_df,
      to_rows = names(objs)[1],
      to_cols = names(objs)[2],
      values = newcol
    )

  } else {
    # no split
    if (method == "kendall_zi") {
      stop("kendall_zi not implemented in non-split.")
    }
    corrobj <- psych::corr.test(
      x = avg_expr[[1]][[1]],
      y = avg_expr[[2]][[1]],
      method = method)
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

  dendroplot <- NULL
  hc <- NULL
  if (isSymmetric(corr_mat)) {
    dist_mat <- stats::as.dist(1 - corr_mat)
    hc <- stats::hclust(dist_mat)
    hcdata <- ggdendro::dendro_data(hc, type = "rectangle")
    dendroplot <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = ggdendro::segment(hcdata),
                            ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::geom_text(data = ggdendro::label(hcdata),
                         ggplot2::aes(x = x, y = y, label = label, hjust = 0),
                         size = 3) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_reverse(expand = c(0.2, 0)) +
      ggplot2::theme_void()
  }


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
  # set corr value to NA when min_cells is not met by both groups
  # option for min_cells in both vs. in at least one
  rdf_plot <- rdf |>
    dplyr::mutate(!!rlang::sym(method) := ifelse(!!rlang::sym(ncellcol[1]) < min_cells, NA, !!rlang::sym(method))) |>
    dplyr::mutate(!!rlang::sym(method) := ifelse(!!rlang::sym(ncellcol[2]) < min_cells, NA, !!rlang::sym(method)))

  plot <- Gmisc::fastDoCall(what = fcexpr::heatmap_long_df, #fcexpr::
                            args = c(list(df = rdf_plot,
                                          groups = names(rdf)[1],
                                          features = names(rdf)[2],
                                          values = method,
                                          heatmap_ordering_args = list(feature_order = "none",
                                                                       group_order = "none")),
                                     heatmap_long_df_args))

  ret <- list(
    corrobj = corrobj,
    plot = plot,
    corr_df_plot = rdf_plot,
    corr_df = rdf,
    dendroplot = dendroplot,
    hclust = hc
  )

  if (!is.null(split)) {
    ret <- c(list(corr_df_split = corr_df), ret)
  }
  return(ret)
}

checks <- function(objs,
                   meta_cols,
                   split,
                   assay,
                   features,
                   method) {

  if (methods::is(objs, "Seurat")) {
    objs <- list(objs, objs)
  }
  objs <- scexpr:::check.SO(objs, assay = assay, length = 2)
  assay <- rlang::arg_match(assay, names(objs[[1]]@assays))
  method <- rlang::arg_match(method, c("pearson","spearman", "kendall", "kendall_zi"))

  if (length(meta_cols) == 1) {
    meta_cols <- c(meta_cols, meta_cols)
  }
  if (!meta_cols[1] %in% names(objs[[1]]@meta.data) || !meta_cols[2] %in% names(objs[[2]]@meta.data)) {
    stop("One of meta_cols not found in respective objs.")
  }

  split_intersect <- NULL
  if (!is.null(split)) {
    split <- scexpr:::check.features(objs, features = split, rownames = F)
    split_intersect <- intersect(objs[[1]]@meta.data[[split]], objs[[2]]@meta.data[[split]])
    if (length(split_intersect) == 0) {
      stop("No intersecting levels in ", split, " column found in SOs.")
    }

    if (length(split_intersect) != length(unique(objs[[1]]@meta.data[[split]])) || length(split_intersect) != length(unique(objs[[2]]@meta.data[[split]]))) {
      message("not all levels of split in obj1 and obj2 overlap or intersect.")
      message("Intersecting levels: ", paste(split_intersect, collapse = ", "), ".")
    }
    split_intersect <- as.character(split_intersect)
    names(split_intersect) <- split_intersect
  }

  if (length(intersect(features, c("all", "pca"))) == 2) { # length(features) %in% c(1,2) &&
    features <- rlang::arg_match(features, c("all", "pca"))
  }

  if (length(features) == 1 && features == "all") {
    features <- Reduce(intersect, purrr::map(objs, scexpr:::get_gene_features))
  } else if (length(features) == 1 && features == "pca") {
    features2 <- names(scexpr:::check.reduction(objs, reduction = "pca")[[1]])
    if (length(features2)>1) {
      message("more than 1 pca found: ", paste(features2, collapse = ","), ". using first.")
      features2 <- features2[1]
    }
    if (features != features2) {
      message("using pca: ", features2)
      features <- features2
    }
    features <- Reduce(intersect,  purrr::map(objs, ~rownames(.x@reductions[[features]]@feature.loadings)))
  } else {
    # feature vector provided
    features <- intersect(Reduce(intersect, purrr::map(objs, scexpr:::get_gene_features)), features)
  }

  if (!length(features)) {
    stop("no features left.")
  } else {
    message("using ", length(features), " features.")
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
