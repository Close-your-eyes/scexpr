#' Transfer cluster labels, i.e. celltype classification based on gene sets or
#' modules
#'
#' In comparison to labeltransfer_singler a test object and modules (gene sets)
#' are starter material but not a ref object. Modules may be DE genes defining a
#' cell population in another dataset. Module score and test for enrichment with
#' fgsea are used to find matching populations (clusters) in test_obj.
#'
#' @param test_obj Seurat object
#' @param test_clusters meta.data column in test_obj
#' @param modules list of named gene sets / modules / pathways
#' @param AddModuleScore_UCell_args arguments to UCell::AddModuleScore_UCell
#'
#' @returns list
#' @export
#'
#' @examples
#' \dontrun{
#'   ## get modules/celltype signatrex from reference
#'   ref_modules <- scexpr::find_all_marker(ref, meta_col = "celltype_ref")
#'   ref_modules2 <- ref_modules |>
#'     dplyr::filter(avg_log2FC>1.5 & pct_in>30 & padj<0.00001) # filter for relevant markers
#'
#'   # optionally filter for variable features in test object
#'   # ref_modules2 <- ref_modules2 |> dplyr::filter(feature %in% scexpr::get_var_features(test))
#'   ref_modules_lst <- split(ref_modules2$feature, ref_modules2$group)
#'
#'   res <- labeltransfer_module_gsea(test,
#'                                    test_clusters = "SCT_pca_snn_res.2", # use an overclustering
#'                                    modules = ref_modules_lst,
#'                                    AddModuleScore_UCell_args = list(ncores = 8, name = "")) #leave name = ""
#' }
labeltransfer_module_gsea <- function(test_obj,
                                      test_clusters,
                                      modules,
                                      AddModuleScore_UCell_args = list(ncores = 1,
                                                                       name = "")) {
  if (!requireNamespace("UCell", quietly = T)) {
    BiocManager::install("UCell")
  }
  if (!requireNamespace("fgsea", quietly = T)) {
    devtools::install_github("alserglab/fgsea")
  }
  if (!requireNamespace("fcexpr", quietly = T)) {
    devtools::install_github("Close-your-eyes/fcexpr")
  }
  if (!requireNamespace("brathering", quietly = T)) {
    devtools::install_github("Close-your-eyes/brathering")
  }

  message("GSEA results may give a nice picture of similarities.")

  if (!methods::is(test_obj, "Seurat")) {
    stop("test_obj must be Seurat.")
  }
  if (missing(test_clusters)) {
    stop("test_clusters missing.")
  }
  if (missing(modules)) {
    stop("modules missing.")
  }
  if (!is.list(modules) || is.null(names(modules))) {
    stop("modules must named list.")
  }
  if (!test_clusters %in% names(test_obj@meta.data)) {
    stop("test_clusters not found in test_obj@meta.data.")
  }
  if (!"name" %in% names(AddModuleScore_UCell_args)) {
    AddModuleScore_UCell_args[["name"]] <- "_Ucell"
  }

  module_names <- paste0(names(modules), AddModuleScore_UCell_args[["name"]])

  # from seurat object returned: pull meta columns immediately
  score_df <- Gmisc::fastDoCall(UCell::AddModuleScore_UCell,
                                args = c(list(test_obj, features = modules),
                                         AddModuleScore_UCell_args))@meta.data[,c(test_clusters, module_names),drop = F]
  # avg scores by cluster
  score_df_avg <- dplyr::summarize(score_df,
                                   dplyr::across(
                                     dplyr::all_of(module_names),
                                     mean,
                                     na.rm = TRUE),
                                   .by = !!rlang::sym(test_clusters)) |>
    dplyr::arrange(!!rlang::sym(test_clusters))
  rownames(score_df_avg) <- score_df_avg[[test_clusters]]

  score_df_avg_long <- brathering::mat_to_df_long(
    score_df_avg[,-which(names(score_df_avg) == test_clusters)],
    colnames_to = "module",
    rownames_to = test_clusters,
    values_to = "score")

  scores_plot <- fcexpr::heatmap_long_df(score_df_avg_long,
                                         groups = "module",
                                         features = test_clusters,
                                         values = "score",
                                         values_zscored = F) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::geom_text(data = dplyr::slice_max(score_df_avg_long, order_by = score,
                                               n = 1,
                                               by = !!rlang::sym(test_clusters)),
                       mapping = ggplot2::aes(label = round(score, 2)))

  # get order from plot
  xyorder <- brathering::gg_get_axis_text(scores_plot)
  score_df_avg <- score_df_avg[xyorder$y,]

  score_best <- dplyr::slice_max(
    score_df_avg_long,
    score,
    n = 1,
    by = !!rlang::sym(test_clusters)) |>
    as.data.frame()
  rownames(score_best) <- score_best[[test_clusters]]
  score_best <- score_best[xyorder$y, ]

  conv <- stats::setNames(score_best$module, score_best[[test_clusters]])
  conv <- conv[xyorder$y]


  ## ---- do smth similar with gsea -----
  # s2n <- gsea_s2n_groupwise(test_obj, test_clusters)
  res <- gsea_groupwise(test_obj,
                        test_clusters,
                        fgseaMultilevel_args = list(pathways = modules))
  # res[["gsea"]] <- dplyr::mutate(res[["gsea"]], padj2 = -log10(padj))
  # res[["meta"]][[test_clusters]] <- score_df[[test_clusters]]
  # res[["meta"]] <- dplyr::relocate(res[["meta"]], !!rlang::sym(test_clusters), 1)
  resdf <- purrr::map_dfr(res, ~.x[["data"]]) |>
    dplyr::mutate(padj2 = -log(padj))

## plot ES as NES was NA sometimes
  scores_plot2 <- fcexpr::heatmap_long_df(resdf,
                                          groups = "pathway",
                                          features = test_clusters,
                                          values = "ES",
                                          dotsizes = "padj2",
                                          values_zscored = F,
                                          theme_args = list(panel.grid = ggplot2::element_blank(),
                                                            axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))) +
    ggplot2::geom_text(data = dplyr::slice_max(resdf,
                                               order_by = ES,
                                               n = 1,
                                               by = !!rlang::sym(test_clusters)),
                       mapping = ggplot2::aes(label = round(ES, 2)))
  xyorder2 <- brathering::gg_get_axis_text(scores_plot2)

  score_best2 <- dplyr::slice_max(
    resdf,
    NES,
    n = 1,
    by = !!rlang::sym(test_clusters)) |>
    as.data.frame()
  rownames(score_best2) <- score_best2[[test_clusters]]
  score_best2 <- score_best2[xyorder2$y, ]

  conv2 <- stats::setNames(score_best2$pathway, score_best2[[test_clusters]])
  conv2 <- conv2[xyorder2$y]

  score_best3 <- dplyr::left_join(
    dplyr::select(score_best, 1,2),
    dplyr::select(score_best2, 1,2),
    by = test_clusters)

  # gsea plot with order from module plot
  ### fails if ucell suffix is used?!
  suppressMessages(capture.output(scores_plot3 <- scores_plot2 +
                                    ggplot2::scale_x_discrete(limits = xyorder$x) +
                                    ggplot2::scale_y_discrete(limits = xyorder$y)))
  scores_plot3 <- patchwork::wrap_plots(scores_plot, scores_plot3)

  return(list(modulescore = list(score = score_df,
                                 score_avg_wide = score_df_avg,
                                 score_avg_long = score_df_avg_long,
                                 score_avg_long_plot = scores_plot,
                                 score_avg_best = score_best,
                                 conv = conv),
              gsea = list(#score = res[["meta"]],
                          score_avg_long = resdf,
                          score_avg_long_plot = scores_plot2,
                          score_avg_best = score_best2,
                          conv = conv2),
              module_vs_gsea = list(score_avg_best = score_best3,
                                    score_avg_long_plots = scores_plot3)))
}
