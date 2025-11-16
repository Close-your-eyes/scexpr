#' Assign clusters based on multiple module scores
#'
#' @param obj seurat object
#' @param modules list of modules/gene sets/pathways
#' @param method_score scoring method
#' @param method_cluster cluster method
#' @param AddModuleScore_UCell_args args to UCell::AddModuleScore_UCell when
#' method_score is module_score
#' @param fgseaMultilevel_args args to fgsea::fgseaMultilevel when method_score
#' is gsea
#' @param get_layer_args args to scexpr::getLayer in scexpr::gsea_groupwise
#' when method_score is gsea
#' @param n_times_modules_subcluster
#' @param gsea_overcluster_resolution
#'
#' @returns list
#' @export
#'
#' @examples
discretize_module_score_suggestive_clusters <- function(obj,
                                                        modules,
                                                        method_score = c("module_score", "gsea"),
                                                        method_cluster = c("louvain", "hclust"),
                                                        n_times_modules_subcluster = 1,
                                                        AddModuleScore_UCell_args = list(slot = "data"),
                                                        fgseaMultilevel_args = list(pathways = modules),
                                                        get_layer_args = list(),
                                                        gsea_overcluster_resolution = 2) {

  method_score <- rlang::arg_match(method_score)
  method_cluster <- rlang::arg_match(method_cluster)

  # method_score <- "module_score"
  # method_cluster <- "louvain"

  if (!is.list(modules)) {
    stop("modules must be a list of gene sets.")
  }

  if (length(modules) < 2) {
    stop("need at least 2 modules.")
  }

  if (!"name" %in% names(AddModuleScore_UCell_args)) {
    AddModuleScore_UCell_args[["name"]] <- "_UCell"
  }

  if (method_score == "module_score") {
    #meta_before <- names(obj@meta.data)
    obj <- Gmisc::fastDoCall(UCell::AddModuleScore_UCell,
                             args = c(list(obj = obj, features = modules),
                                      AddModuleScore_UCell_args))
    #meta_after <- names(obj@meta.data)
    #new_cols <- setdiff(meta_after, meta_before)
    new_cols <- paste0(names(modules), AddModuleScore_UCell_args[["name"]])
    score_df <- obj@meta.data[,new_cols]
  } else if (method_score == "gsea") {

    #meta_before <- names(obj@meta.data)
    # overcluster
    obj <- obj |>
      Seurat::FindNeighbors() |>
      Seurat::FindClusters(resolution = gsea_overcluster_resolution) # make argument
    #meta_after <- names(obj@meta.data)
    #new_cols <- setdiff(meta_after, meta_before) # n=1

    new_cols <- paste0(Seurat::DefaultAssay(obj), "_snn.res.", gsea_overcluster_resolution)
    new_cols <- gsub("\\.0$", "", new_cols)

    ## check output
    score_df <- gsea_groupwise(
      obj = obj,
      group = new_cols,
      fgseaMultilevel_args = fgseaMultilevel_args,
      get_layer_args = get_layer_args
    )
    score_df <- dplyr::select(score_df[["meta"]], dplyr::ends_with("_ES"))
    if (anyNA(score_df)) {
      stop("Some ES of GSEA are NA. Cannot continue with that")
    }
    obj <- Seurat::AddMetaData(obj, score_df)
    new_cols <- names(score_df)
  }

  score_mat <- scale(as.matrix(score_df))

  if (method_cluster == "louvain") {
    ncluster <- 0
    res <- 0.1
    # iterate louvain to match ncluster to len(modules)
    message("iterate louvain.")
    while(ncluster < n_times_modules_subcluster*length(modules)) {
      message(res)
      module_based_clusters <- fcexpr::get_louvain_cluster(exprs = score_mat,
                                                           FindClusters_args = list(resolution = res))
      ncluster <- length(unique(module_based_clusters[,1]))
      res <- res + 0.1
    }
    #res0 <- res
    incr <- 0.01
    while (ncluster > n_times_modules_subcluster*length(modules)) {
      res <- res - incr
      message(res)
      module_based_clusters <- fcexpr::get_louvain_cluster(exprs = score_mat,
                                                           FindClusters_args = list(resolution = res))
      ncluster <- length(unique(module_based_clusters[,1]))
      if (ncluster < n_times_modules_subcluster*length(modules)) {
        res <- res + incr
        incr <- incr/10
      }
    }

  } else if (method_cluster == "hclust") {
    module_based_clusters <- fcexpr::get_hclust_clusters(exprs = score_mat,
                                                         k = n_times_modules_subcluster*length(modules))
  }
  module_based_clusters <- as.data.frame(module_based_clusters)
  names(module_based_clusters) <- "cluster"
  module_based_clusters$cluster <- as.character(module_based_clusters$cluster)
  module_based_clusters$id <- rownames(score_mat)

  # assign cluster names
  dflong <- brathering::mat_to_df_long(x = score_mat,
                                       rownames_to = "id",
                                       colnames_to = "module",
                                       values_to = "score") |>
    dplyr::left_join(module_based_clusters, by = "id")


  assigned <- best_cluster_assignment_by_module_scores(df = dflong,
                                                       score = "score",
                                                       module = "module",
                                                       cluster = "cluster")
  ## use finer clusters?
  # assigned <- best_cluster_assignment_by_module_scores(df = dflong,
  #                                                      score = "score",
  #                                                      module = "module",
  #                                                      cluster = "cluster")

  assigned$df <- dplyr::mutate(assigned$df, is = assigned == module)

  p0 <- ggplot2::ggplot(assigned$df, ggplot2::aes(x = cluster, y = score, color = is)) +
    ggplot2::geom_jitter(width = 0.1) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(ggplot2::vars(module))

  obj@meta.data$id <- rownames(obj@meta.data)
  obj <- Seurat::AddMetaData(obj, assigned$df[,c("id", "assigned")] |> dplyr::distinct() |> tibble::column_to_rownames("id"))
  p1 <- feature_plot2(obj, features = c("assigned", new_cols))

  message("Use similar command as below to add assigned clusters to meta data. Make sure rownames of meta data exist as column 'id' and replace return with your return variable.")
  print("obj <- Seurat::AddMetaData(obj, return$assigned$df[,c('id', 'assigned')] |> dplyr::distinct() |> tibble::column_to_rownames('id'))")

  return(list(plots = list(p0 = p0, p1 = p1), assigned = assigned)) # obj = obj,
}


#' Assign modules/gene sets/pathways to anonymous clusters in optimal way
#'
#' Assign module scores that cover each cell with a reasonable maximum and
#' have re-calculated clustered assigned to them.
#' This function is called by discretize_module_score_suggestive_clusters().
#'
#' @param df long data frame
#' @param score module score column
#' @param module module name column
#' @param cluster cluster column which module is to be mapped to
#'
#' @returns
#'
#' @examples
best_cluster_assignment_by_module_scores0 <- function(df,
                                                      score,
                                                      module,
                                                      cluster) {

  #  no recycling as in best_cluster_assignment_by_module_scores

  # pt_scores <- so@meta.data[,grep("ucell", names(so@meta.data), ignore.case = T, value = T)]
  # pt_subcluster <- fcexpr::get_louvain_cluster(scale(as.matrix(pt_scores)), FindClusters_args = list(resolution = seq(0.1,2,0.1)))
  # pt_subcluster <- as.data.frame(pt_subcluster) |>
  #   dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  # rownames(pt_subcluster) <- rownames(pt_scores)
  # so <- Seurat::AddMetaData(so, pt_subcluster)

  # so_data <- scexpr::get_data(so, feature = c(grep("ucell", names(so@meta.data), ignore.case = T, value = T))) |>
  #   dplyr::bind_rows(.id = "pt_pathway") |>
  #   dplyr::rename("score" = feature) |>
  #   dplyr::left_join(scexpr::get_data(so, feature = names(pt_subcluster)[4])[[1]][,c("id", "feature")] |> dplyr::rename("cluster" = feature)) |>
  #   dplyr::mutate(pt_pathway2 = ifelse(pt_pathway %in% c("PT-S1_UCell", "PT-S2_UCell"), "PT-S12_UCell", pt_pathway)) |>
  #   dplyr::mutate(pt_pathway2 = ifelse(pt_pathway2 %in% c("PT-New3_UCell", "PT-New4_UCell"), "PT-New34_UCell", pt_pathway2))
  # best_cluster_assignment_by_module_scores(so_data,
  #                                          score = "score",
  #                                          module = "pt_pathway2",
  #                                          cluster = "cluster")
  # ggplot(so_data, aes(x = cluster, y = score)) +
  #   geom_jitter(width = 0.1) +
  #   geom_boxplot() +
  #   facet_wrap(vars(pt_pathway2))


  mean_scores <- dplyr::summarise(
    df,
    mean_score = median(!!rlang::sym(score), na.rm = TRUE),
    .by = c(!!rlang::sym(cluster), !!rlang::sym(module))
  )

  # 2. Make wide matrix
  score_matrix <- tidyr::pivot_wider(
    mean_scores,
    names_from = !!rlang::sym(module),
    values_from = mean_score,
    values_fill = 0
  )

  # Keep row (cluster) names
  clusters <- score_matrix[[cluster]]
  score_matrix <- as.matrix(score_matrix[ , -1])
  rownames(score_matrix) <- clusters

  # 3. Convert to cost matrix (since solve_LSAP minimizes)
  cost_matrix <- max(score_matrix) - score_matrix

  # 4. Solve assignment (Hungarian algorithm)
  assignment <- clue::solve_LSAP(cost_matrix)  # returns column indices

  # 5. Build mapping
  mapping <- data.frame(rownames(score_matrix), colnames(score_matrix)[assignment])
  names(mapping) <- c("cluster", "assigned")

  # 6. Apply mapping
  df <- dplyr::left_join(df, mapping, by = "cluster")

  return(list(df = df, mapping = mapping))
}



#' Assign modules/gene sets/pathways to anonymous clusters in optimal way
#'
#' Assign module scores that cover each cell with a reasonable maximum and
#' have re-calculated clustered assigned to them.
#'
#' @param df long data frame
#' @param score module score column
#' @param module module name column
#' @param cluster cluster column which module is to be mapped to
#'
#' @returns
#' @export
#'
#' @examples
best_cluster_assignment_by_module_scores <- function(df,
                                                     score,
                                                     module,
                                                     cluster) {


  mean_scores <- dplyr::summarise(
    df,
    mean_score = median(!!rlang::sym(score), na.rm = TRUE),
    .by = c(!!rlang::sym(cluster), !!rlang::sym(module))
  )

  # 2. Make wide matrix
  score_matrix <- tidyr::pivot_wider(
    mean_scores,
    names_from = !!rlang::sym(module),
    values_from = mean_score,
    values_fill = 0
  )

  # Keep row (cluster) names
  clusters <- score_matrix[[cluster]]
  score_matrix <- as.matrix(score_matrix[ , -1])
  rownames(score_matrix) <- clusters

  # initialization
  n_clusters <- nrow(score_matrix)
  n_modules <- ncol(score_matrix)
  remaining_clusters <- rownames(score_matrix)
  all_mappings <- list()
  round <- 1

  # iterative assignment
  while (length(remaining_clusters) > 0) {
    n_current <- min(length(remaining_clusters), n_modules)

    # subset the score matrix for remaining clusters
    current_scores <- score_matrix[remaining_clusters, , drop = FALSE]
    cost_matrix <- max(current_scores) - current_scores

    # run Hungarian algorithm
    assignment <- clue::solve_LSAP(cost_matrix[1:n_current, ])

    # build mapping for this batch
    mapping <- data.frame(
      cluster = remaining_clusters[1:n_current],
      assigned = colnames(current_scores)[assignment],
      iteration = round
    )
    names(mapping)[1] <- cluster

    all_mappings[[round]] <- mapping

    # remove assigned clusters
    remaining_clusters <- setdiff(remaining_clusters, mapping[[cluster]])
    round <- round + 1
  }

  # combine all mappings
  mapping_df <- dplyr::bind_rows(all_mappings)

  # optional: give reused pathways suffixes (_2, _3, etc.)
  mapping_df <- dplyr::mutate(mapping_df,
                              assigned_pathway_unique = ifelse(dplyr::n() == 1,
                                                               assigned,
                                                               paste0(assigned, "_", dplyr::row_number())),
                              .by = assigned)

  # merge with original data
  df <- dplyr::left_join(df, mapping_df, by = cluster)

  return(list(df = df, mapping = mapping_df))
}
