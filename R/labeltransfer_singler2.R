#' Transfer cluster labels, i.e. celltype classification between objects
#'
#' How SingleR computes. Marker selection: For each label, SingleR identifies genes
#' that distinguish it from the others (de.method argument). Scoring: For a given
#' query cell, it computes correlation between that cell’s expression profile and
#' the reference expression profile restricted to those marker genes. Best
#' label: The label with the highest score is assigned. Fine-tuning: It optionally
#' refines the set of genes for the top labels and recomputes scores → final label.
#' How to interpret. Higher score = higher similarity between query cell and the
#' reference profile. Relative differences matter more than absolute values (e.g.
#' 0.6 vs 0.3 is decisive; 0.61 vs 0.59 is ambiguous). SingleR also reports:
#' score.max → the best score per cell. delta → the difference between top and
#' runner-up (margin of confidence). pruned.labels → removes labels when the margin
#' is too small (set to NA = “ambiguous”). Not strict probabilities. Scores are not
#' calibrated probabilities (they don’t sum to 1). They’re correlations, so you can
#' think of them as “confidence-like.” If you want probability-like values, you can
#' normalize across labels per cell (e.g. softmax), but this is an approximation.
#'
#' fine.tune may cause that the assigned label is not the top score.
#' for single cell data use de.method = "wilcox", for bulk as ref "classic".
#'
#' @param test_obj Seurat object or gene x cell matrix (sparse/dense) to which
#' labels should be transferred
#' @param ref_obj Seurat object, SummarizedExperiment or gene x cell matrix
#' @param ref_labels when ref_obj is Seurat: meta.data column, when
#' SummarizedExperiment: name of element in ref_obj_colData_listData, when
#' matrix: label vector of length ncol(ref_obj)
#' @param test_clusters NULL to run SingleR on single cells; to run on clusters,
#' when test_obj is Seurat: meta.data column, when matrix: vector of cluster
#' idents of length ncol(test_obj)
#' @param name_prefix prefix to all returned columns
#' @param singler_args arguments to SingleR::SingleR
#' @param get_layer_args arguments to scexpr::get_layer (applies when test_obj
#' or ref_obj is Seurat)
#'
#' @returns list of data frame, matrix, named vector, score plot of SingleR
#' results
#' @export
#'
#' @examples
#' \dontrun{
#' out <- labeltransfer_singler(test_obj = seu_test,
#'                              ref_obj = seu_ref,
#'                              ref_labels = seu_ref@meta.data$celltypes,
#'                              test_clusters = seu_test@meta.data$SCT_snn_res.0.8,
#'                              name_prefix = "enghard")
#' }
labeltransfer_singler2 <- function(test_obj,
                                  ref_obj,
                                  ref_labels,
                                  test_clusters = NULL,
                                  name_prefix = NULL,
                                  singler_args = list(de.method = "wilcox",
                                                      fine.tune = T,
                                                      aggr.ref = F),
                                  get_layer_args = list(layer = "data",
                                                        assay = "RNA")) {

  if (!requireNamespace("SingleR", quietly = T)) {
    BiocManager::install("SingleR")
  }
  if (!requireNamespace("fcexpr", quietly = T)) {
    devtools::install_github("Close-your-eyes/fcexpr")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  if (missing(ref_labels) || !length(ref_labels)) {
    stop("ref_labels are missing.")
  }

  # name_prefix <- deparse(substitute(ref_obj))

  #test: test_clusters is null - on cell level



  c(test_obj,
    test_clusters,
    test_clusters_name) %<-% prep_test_vars(test_obj = test_obj,
                                            test_clusters = test_clusters,
                                            get_layer_args = get_layer_args)

  c(ref_obj,
    ref_labels,
    ref_labels_name) %<-% prep_ref_vars(ref_obj = ref_obj,
                                        ref_labels = ref_labels,
                                        get_layer_args = get_layer_args)

  if (ref_labels_name == test_clusters_name) {
    stop("ref_labels and test_clusters cannot have the same name.")
  }
  prefix <- gsub("^__", "", paste0(name_prefix, "__", ref_labels_name, "__"))

  # BPPARAM <- BiocParallel::SerialParam()
  # num.threads = BiocParallel::bpnworkers(BPPARAM)
  # if (is.null(test_clusters)) {
  #   BPPARAM <- BiocParallel::MulticoreParam()
  #   num.threads = BiocParallel::bpnworkers(BPPARAM)
  #   message("no test_clusters.
  #           you may over-cluster your cells and provide result as test_clusters.
  #           this would flatten noise in single cells.
  #           cells belonging to the same very fine cluster are correlated and may hence be of same type.
  #           running SingleR on single cells takes a while. trying multicore though: ", num.threads, " threads.")
  # }

  BPPARAM <- BiocParallel::MulticoreParam()
  num.threads = BiocParallel::bpnworkers(BPPARAM)

  ## sc happens always
  labels_sc <- Gmisc::fastDoCall(SingleR::SingleR,
                                 args = c(list(
                                   test = test_obj,
                                   ref = ref_obj,
                                   labels = ref_labels,
                                   #clusters = test_clusters,
                                   BPPARAM = BPPARAM,
                                   num.threads = num.threads
                                 ), singler_args))
  inds_sc <- names(labels_sc@listData)[-which(names(labels_sc@listData) == "scores")]
  names(inds_sc) <- paste0(prefix, inds_sc)
  lablist_sc <- purrr::map(inds_sc, ~stats::setNames(labels_sc@listData[[.x]], labels_sc@rownames))
  labels_name <- grep("__labels$", names(lablist_sc), value = T)
  prunelabels_name <- grep("__pruned.labels$", names(lablist_sc), value = T)

  ## scores
  score_mat_sc <- labels_sc@listData[["scores"]]
  rownames(score_mat_sc) <- labels_sc@rownames
  score_df_sc <- prep_score_df(score_mat_sc, lablist_sc,
                               ref_labels_name,
                               labels_name, prunelabels_name,
                               rownames_to = "test_cells")
  names(score_df_sc)[-1] <- paste0(names(score_df_sc)[-1], "_sc")
  score_df_sc <- score_df_sc |>
  dplyr::left_join(
    score_df_sc |>
      dplyr::filter(is_label_sc) |>
      dplyr::select(!!rlang::sym("test_cells"), !!rlang::sym(paste0(ref_labels_name, "_sc"))) |>
      dplyr::rename(!!paste0(paste0(ref_labels_name, "_sc"), "_assigned_label") := !!rlang::sym(paste0(ref_labels_name, "_sc"))),
    by = dplyr::join_by(!!rlang::sym("test_cells")))

  if (!is.null(test_clusters)) {
    labels_clust <- Gmisc::fastDoCall(SingleR::SingleR,
                                      args = c(list(
                                        test = test_obj,
                                        ref = ref_obj,
                                        labels = ref_labels,
                                        clusters = test_clusters,
                                        BPPARAM = BPPARAM,
                                        num.threads = num.threads
                                      ), singler_args))

    ## clustered only if test_clusters provided
    inds_clust <- names(labels_clust@listData)[-which(names(labels_clust@listData) == "scores")]
    names(inds_clust) <- paste0(prefix, inds_clust)
    lablist_clust <- purrr::map(inds_clust, ~stats::setNames(labels_clust@listData[[.x]], labels_clust@rownames))


    score_mat_clust <- labels_clust@listData[["scores"]]
    rownames(score_mat_clust) <- labels_clust@rownames
    score_df_clust <- prep_score_df(score_mat_clust, lablist_clust,
                                 ref_labels_name,
                                 labels_name, prunelabels_name,
                                 rownames_to = test_clusters_name)
    names(score_df_clust)[-1] <- paste0(names(score_df_clust)[-1], "_clust")
  }








  ### fine.tune may cause that the assigned label is not the top score
  # color scores of assigned labels green in scores_plot; top scores in black
  # handle scores separately below to detect changes due to fine tuning



  if (any(score_df$is_max_score != score_df$is_label)) {
    message("maxore does not always match assigned labels. that means fine tuning affected label assignment.")
    print(dplyr::filter(score_df, is_max_score != is_label))
  } else {
    message("maxore always matches assigned labels.")
  }
  if (any(score_df$is_pruned_label != score_df$is_label)) {
    message("labels and pruned labels are not always the same. pruning removed labels.")
    print(dplyr::filter(score_df, is_pruned_label != is_label))
  } else {
    message("labels and pruned labels are equal.")
  }
  scores_plot_pseudobulk <- NULL
  if (!is.null(test_clusters)) {
    ## df for text labels on scores_plot
    score_df_textlabels <-
      score_df |>
      dplyr::filter(is_max_score | (is_label & !is_max_score)) |>
      tidyr::pivot_longer(cols = c(is_max_score, is_label)) |>
      dplyr::filter(value) |>
      dplyr::mutate(score_label = factor(name, levels = c("is_label", "is_max_score"))) |>
      dplyr::mutate(score2 = gsub("^0", "", format(round(score, 2), nsmall = 2)))

    scores_plot <- fcexpr::heatmap_long_df(
      df = score_df,
      groups = names(score_df)[2],
      features = names(score_df)[1],
      values = names(score_df)[3]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::geom_text(data = score_df_textlabels,
                         mapping = ggplot2::aes(label = score2, color = score_label)) +
      ggplot2::scale_color_manual(values = c("black", "hotpink"))

  } else {

    score_df_pseudobulk <- score_df |>
      dplyr::summarise(score = median(score), .by = c(!!rlang::sym(paste0(ref_labels_name, "_assigned_label")), !!rlang::sym(ref_labels_name)))
    score_df_textlabels <-
      score_df_pseudobulk |>
      dplyr::filter(!!rlang::sym(paste0(ref_labels_name, "_assigned_label")) == !!rlang::sym(ref_labels_name)) |>
      dplyr::mutate(score2 = gsub("^0", "", format(round(score, 2), nsmall = 2)))

    scores_plot_pseudobulk <- fcexpr::heatmap_long_df(
      df = score_df_pseudobulk,
      groups = names(score_df_pseudobulk)[1],
      features = names(score_df_pseudobulk)[2],
      values = names(score_df_pseudobulk)[3]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::geom_text(data = score_df_textlabels,
                         mapping = ggplot2::aes(label = score2))

    scores_plot <- ggplot2::ggplot(dplyr::filter(score_df, is_label),
                                   ggplot2::aes(x = !!rlang::sym(names(score_df)[2]), y = !!rlang::sym(names(score_df)[3]))) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.05, size = 0.3) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  labelore <- dplyr::filter(score_df, is_label)

  lablist[[gsub("^__", "", paste0(name_prefix, "__", paste0(ref_labels_name, "_ores")))]] <-
    stats::setNames(labelore$score, labelore[[test_clusters_name]])

  labdf <- NULL
  gttable <- NULL
  if (!is.null(test_clusters)) {

    labdf <- purrr::reduce(purrr::map(names(lablist), ~stats::setNames(utils::stack(lablist[[.x]]), c(.x, test_clusters_name))),
                           dplyr::left_join,
                           by = test_clusters_name) |>
      dplyr::select(2,1,5,4,3)
    rownames(labdf) <- labdf[[test_clusters_name]]

    labdfgt <- labdf
    names(labdfgt) <- gsub(prefix, "", names(labdfgt))
    labdfgt <- dplyr::mutate(labdfgt, scores = round(scores, 3), delta.next = round(delta.next, 3))

    gttable <- brathering::gt_tight(labdfgt) |>
      gt::tab_spanner(
        label = gsub("__$", "", prefix),
        columns = names(labdfgt)[-1],
        id = "prefix") |>
      gt::data_color(
        columns = scores,
        fn = scales::col_numeric(
          palette = c(colrr::col_pal("RdBu", direction = -1)),
          domain = range(labdfgt$scores))) |>
      gt::data_color(
        columns = delta.next,
        fn = scales::col_numeric(
          palette = c(colrr::col_pal("RdBu", direction = -1)),
          domain = range(labdfgt$delta.next)))
    # gt::tab_style(
    #   style = gt::cell_fill(color = "lightgreen"),
    #   locations = gt::cells_body(
    #     columns = c(),
    #     rows = c()
    #   )
    # )

    metadata <- stats::setNames(data.frame(test_clusters), test_clusters_name) |>
      dplyr::left_join(labdf, by = test_clusters_name) |>
      dplyr::select(-!!rlang::sym(test_clusters_name))
  } else {
    metadata <- as.data.frame(dplyr::bind_cols(lablist))
  }

  #rownames(metadata) <- test_colnames #labels@rownames
  rownames(metadata) <- labels@rownames

  return(list(
    labels_list = lablist,
    labels_df = labdf,
    metadata = metadata,
    scores = score_mat,
    scores_df = score_df,
    scores_plot = scores_plot,
    scores_plot_pseudobulk = scores_plot_pseudobulk,
    labels_gt = gttable,
    singler_obj = labels
  ))
}


prep_test_vars <- function(test_obj, test_clusters, get_layer_args) {
  test_clusters_name <- "test_cells" # default when no clause below applies
  if (methods::is(test_obj, "Seurat")) {
    if (!is.null(test_clusters)) {
      if (!test_clusters %in% names(test_obj@meta.data)) {
        stop(test_clusters, " not found in test_obj@meta.data.")
      }
      test_clusters_name <- test_clusters
      test_clusters <- test_obj@meta.data[[test_clusters]]
    }
    test_obj <- Gmisc::fastDoCall(get_layer, c(list(obj = test_obj), get_layer_args))
  } else if (methods::is(test_obj, "matrix") || methods::is(test_obj, "sparseMatrix")) {
    if (!is.null(test_clusters)) {
      test_clusters_name <- "test_clusters"
    }
  } else {
    stop("test_obj must be Seurat or (sparse) matrix")
  }
  if (!is.null(test_clusters) && !is.factor(test_clusters)) {
    test_clusters <- as.factor(test_clusters)
  }
  return(list(test_obj, test_clusters, test_clusters_name))
}

prep_ref_vars <- function(ref_obj, ref_labels, get_layer_args) {

  ref_labels_name <- "ref_labels"
  if (methods::is(ref_obj, "Seurat")) {
    if (!ref_labels %in% names(ref_obj@meta.data)) {
      message("labeltransfer_singler: ", ref_labels, " not found.")
      return(NULL)
    }
    ref_labels_name <- ref_labels
    ref_labels <- ref_obj@meta.data[[ref_labels]]
    ref_obj <-  Gmisc::fastDoCall(get_layer, c(list(obj = ref_obj), get_layer_args))
  } else if (methods::is(ref_obj, "SummarizedExperiment")) {
    if (!ref_labels %in% names(ref_obj@colData@listData)) {
      message("labeltransfer_singler: ", ref_labels, " not found.")
      return(NULL)
    }
    ref_labels_name <- ref_labels
    ref_labels <- ref_obj@colData@listData[[ref_labels]]
  } else if (!methods::is(ref_obj, "matrix") && !methods::is(ref_obj, "sparseMatrix")) {
    stop("ref_obj must be Seurat, SummarizedExperiment or (sparse) matrix")
  }

  return(list(ref_obj, ref_labels, ref_labels_name))
}

prep_score_df <- function(score_mat, lablist, ref_labels_name,
                          labels_name, prunelabels_name,
                          rownames_to) {
  score_df <- brathering::mat_to_df_long(
    score_mat,
    colnames_to = ref_labels_name,
    rownames_to = rownames_to,
    values_to = "score") |>
    dplyr::mutate(is_max_score = score == max(score), .by = !!rlang::sym(rownames_to)) |>
    dplyr::left_join(utils::stack(lablist[[labels_name]]) |>
                       dplyr::mutate(ind = as.character(ind)) |>
                       dplyr::rename(!!rlang::sym(ref_labels_name) := values,
                                     !!rlang::sym(rownames_to) := ind) |>
                       dplyr::mutate(is_label = T),
                     by = c(ref_labels_name, rownames_to)) |>
    dplyr::left_join(utils::stack(lablist[[prunelabels_name]]) |>
                       dplyr::mutate(ind = as.character(ind)) |>
                       dplyr::rename(!!rlang::sym(ref_labels_name) := values,
                                     !!rlang::sym(rownames_to) := ind) |>
                       dplyr::mutate(is_pruned_label = T),
                     by = c(ref_labels_name, rownames_to)) |>
    dplyr::mutate(is_label = ifelse(is.na(is_label), F, is_label),
                  is_pruned_label = ifelse(is.na(is_pruned_label), F, is_pruned_label))
}
