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
labeltransfer_singler <- function(test_obj,
                                  ref_obj,
                                  ref_labels,
                                  test_clusters = NULL,
                                  name_prefix = NULL,
                                  singler_args = list(de.method = "classic",
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

  test_clusters_name <- "test_cells" # default when no clause below applies
  test_colnames <- colnames(test_obj)
  if (is.null(test_colnames)) {
    test_colnames <- 1:ncol(test_obj)
  }
  if (methods::is(test_obj, "Seurat")) {
    if (!is.null(test_clusters)) {
      test_clusters_name <- test_clusters
      test_clusters <- test_obj@meta.data[[test_clusters]]
    }
    test_obj <- Gmisc::fastDoCall(get_layer, c(list(obj = test_obj), get_layer_args))
  } else if (!methods::is(ref_obj, "matrix") && !methods::is(ref_obj, "sparseMatrix")) {
    stop("test_obj must be Seurat or (sparse) matrix")
  } else {
    if (!is.null(test_clusters)) {
      test_clusters_name <- "test_clusters"
    }
  }


  if (!is.null(test_clusters) && !is.factor(test_clusters)) {
    test_clusters <- as.factor(test_clusters)
  }

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

  if (ref_labels_name == test_clusters_name) {
    stop("ref_labels and test_clusters cannot have the same name.")
  }

  if (is.null(test_clusters)) {
    BPPARAM <- BiocParallel::MulticoreParam()
    num.threads = BiocParallel::bpnworkers(BPPARAM)
    message("no test_clusters.
            you may over-cluster your cells and provide result as test_clusters.
            this would flatten noise in single cells.
            cells belonging to the same very fine cluster are correlated and may hence be of same type.
            running SingleR on single cells takes a while. trying multicore though: ", num.threads, " threads.")
  } else {
    BPPARAM <- BiocParallel::SerialParam()
    num.threads = BiocParallel::bpnworkers(BPPARAM)
  }

  labels <- Gmisc::fastDoCall(SingleR::SingleR,
                              args = c(list(
                                test = test_obj,
                                ref = ref_obj,
                                labels = ref_labels,
                                clusters = test_clusters,
                                BPPARAM = BPPARAM,
                                num.threads = num.threads
                              ), singler_args))

  lablist <- purrr::map(stats::setNames(c("scores", "labels", "delta.next", "pruned.labels"),
                                        c("scores", "labels", "delta.next", "pruned.labels")), function(x) {
                                          if (x == "scores") {
                                            y <- apply(labels@listData[[x]], 1, max)
                                          } else {
                                            y <- labels@listData[[x]]
                                          }
                                          return(stats::setNames(y, levels(test_clusters)))
                                        })

  # name_prefix + ref_labels as prefix
  newnames <- paste0(name_prefix, "__", paste0(ref_labels_name, "__", names(lablist)))
  newnames <- gsub("^__", "", newnames)
  names(lablist) <- newnames

  labdf <- NULL
  if (!is.null(test_clusters)) {
    labdf <- purrr::reduce(purrr::map(names(lablist), ~stats::setNames(stack(lablist[[.x]]), c(.x, test_clusters_name))),
                           dplyr::left_join,
                           by = test_clusters_name) |>
      dplyr::relocate(!!rlang::sym(test_clusters_name), 1)
    rownames(labdf) <- labdf[[test_clusters_name]]
    names(labdf)[-1] <- newnames

    metadata <- stats::setNames(data.frame(test_clusters), test_clusters_name) |>
      dplyr::left_join(labdf, by = test_clusters_name) |>
      dplyr::select(-!!rlang::sym(test_clusters_name))
  } else {
    metadata <- as.data.frame(dplyr::bind_cols(lablist))
  }
  rownames(metadata) <- test_colnames #labels@rownames


  ## scores
  score_mat <- labels@listData[["scores"]]
  rownames(score_mat) <- labels@rownames
  score_df <- brathering::mat_to_df_long(
    score_mat,
    colnames_to = ref_labels_name,
    rownames_to = test_clusters_name,
    values_to = "score"
  )
  scores_plot <- fcexpr::heatmap_long_df(
    df = score_df,
    groups = names(score_df)[2],
    features = names(score_df)[1],
    values = names(score_df)[3]) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::geom_text(data = dplyr::slice_max(score_df, order_by = score,
                                               n = 1,
                                               by = !!rlang::sym(names(score_df)[1])),
                       mapping = ggplot2::aes(label = round(score, 2)))

  return(list(
    labels_list = lablist,
    labels_df = labdf,
    metadata = metadata,
    scores = score_mat,
    scores_df = score_df,
    scores_plot = scores_plot
  ))
}
