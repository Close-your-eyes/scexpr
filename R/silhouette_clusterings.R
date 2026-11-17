#' Compare silhouette widths across saved clusterings
#'
#' Score the metadata columns recorded in `SeuratObject::Misc(obj,
#' "clusterings")` using the same cells and embedding for every comparison.
#' Euclidean distances are calculated once and reused. Exact silhouette
#' widths require pairwise distances to fit in memory.
#'
#' @param obj A Seurat object containing saved clusterings.
#' @param reduction Name of a shared dimensional reduction, such as a PCA or
#'   Harmony reduction. If `NULL`, use the reduction identified by the saved
#'   SNN clustering names. If none match, use `"pca"` when available. Supply
#'   a name explicitly if the saved clusterings identify several reductions.
#' @param dims Optional integer vector selecting the dimensions of `reduction`
#'   to use for every clustering. By default, use all dimensions.
#'
#' @return A list with `summary`, `by_cluster`, and `plot`. `summary` contains
#'   one row per saved clustering with the number of cells and clusters, the
#'   mean silhouette width across cells, the fraction of negative widths,
#'   and the scoring status. Missing columns, incomplete labels, and
#'   partitions with fewer than two clusters or a cluster per cell receive
#'   missing scores and an explanatory status. `by_cluster` contains the
#'   size, mean width, and fraction of negative widths for each scored
#'   cluster. `plot` is a ggplot ordered by overall mean silhouette width.
#'
#' @seealso [silhouette_plot()] for a cell-level plot of one clustering.
#' @export
#'
#' @examples
#' \dontrun{
#' scores <- silhouette_clusterings(SO)
#' scores$summary
#' scores$plot
#' }
silhouette_clusterings <- function(obj, reduction = NULL, dims = NULL) {
  clusterings <- unique(SeuratObject::Misc(obj, "clusterings"))
  if (!is.character(clusterings) || !length(clusterings) ||
      anyNA(clusterings) || any(!nzchar(clusterings))) {
    stop("Misc(obj, 'clusterings') must contain metadata column names.",
         call. = FALSE)
  }

  reductions <- SeuratObject::Reductions(obj)
  if (is.null(reduction)) {
    matches <- vapply(reductions, function(name) {
      any(startsWith(clusterings, paste0(name, "_snn")))
    }, logical(1))
    reduction <- reductions[matches]
    if (!length(reduction)) {
      reduction <- intersect("pca", reductions)
    }
  }
  if (!is.character(reduction) || length(reduction) != 1L ||
      !reduction %in% reductions) {
    stop("Set `reduction` to the name of one shared embedding in `obj`.",
         call. = FALSE)
  }

  embeddings <- SeuratObject::Embeddings(obj, reduction = reduction)
  if (!is.null(dims)) {
    if (!is.vector(dims, mode = "numeric") || !length(dims) || anyNA(dims) ||
        any(dims < 1 | dims > ncol(embeddings) | dims != trunc(dims)) ||
        anyDuplicated(dims)) {
      stop("`dims` must contain distinct column indices in `reduction`.",
           call. = FALSE)
    }
    embeddings <- embeddings[, dims, drop = FALSE]
  }
  if (any(!is.finite(embeddings))) {
    stop("The selected embedding contains non-finite values.", call. = FALSE)
  }

  metadata <- obj[[]][rownames(embeddings), , drop = FALSE]
  labels <- lapply(clusterings, function(name) metadata[[name]])
  n_clusters <- vapply(labels, function(x) {
    if (is.null(x)) NA_integer_ else length(unique(x[!is.na(x)]))
  }, integer(1))
  summary <- data.frame(
    clustering = clusterings,
    n_cells = nrow(embeddings),
    n_clusters = n_clusters,
    mean_silhouette = NA_real_,
    fraction_negative = NA_real_,
    status = dplyr::case_when(
      is.na(n_clusters) ~ "Missing metadata column",
      vapply(labels, anyNA, logical(1)) ~ "Missing cell labels",
      n_clusters < 2L ~ "Fewer than two clusters",
      n_clusters == nrow(embeddings) ~ "Every cell is a singleton",
      TRUE ~ "Scored"
    )
  )

  scored <- which(summary$status == "Scored")
  cluster_scores <- vector("list", length(scored))
  if (length(scored)) {
    scexpr:::.ensure_package("cluster")
    distances <- stats::dist(embeddings)
    for (i in seq_along(scored)) {
      row <- scored[[i]]
      groups <- factor(labels[[row]])
      silhouette <- cluster::silhouette(as.integer(groups), dist = distances)
      widths <- silhouette[, "sil_width"]
      summary$mean_silhouette[[row]] <- mean(widths)
      summary$fraction_negative[[row]] <- mean(widths < 0)

      widths_by_cluster <- split(widths, groups)
      cluster_scores[[i]] <- data.frame(
        clustering = clusterings[[row]],
        cluster = names(widths_by_cluster),
        n_cells = lengths(widths_by_cluster),
        mean_silhouette = vapply(widths_by_cluster, mean, numeric(1)),
        fraction_negative = vapply(widths_by_cluster, function(x) mean(x < 0),
                                   numeric(1)),
        row.names = NULL
      )
    }
  }
  by_cluster <- dplyr::bind_rows(
    data.frame(clustering = character(), cluster = character(), n_cells = integer(),
               mean_silhouette = double(), fraction_negative = double()),
    cluster_scores
  )

  summary <- summary[order(summary$mean_silhouette, decreasing = TRUE), , drop = FALSE]
  rownames(summary) <- NULL
  has_score <- summary$status == "Scored"
  plot <- ggplot2::ggplot(summary, ggplot2::aes(y = .data$clustering)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_point(data = by_cluster, ggplot2::aes(x = .data$mean_silhouette),
                        colour = "#5989A5", alpha = 0.45) +
    ggplot2::geom_point(data = summary[has_score, ],
                        ggplot2::aes(x = .data$mean_silhouette),
                        colour = "#13354A", shape = 18, size = 3) +
    ggplot2::geom_text(data = summary[!has_score, ], ggplot2::aes(label = .data$status),
                       x = 0.04, hjust = 0, colour = "grey40", size = 3) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::scale_y_discrete(limits = rev(summary$clustering)) +
    ggplot2::labs(x = "Mean silhouette width", y = NULL,
                  subtitle = paste("Shared embedding", reduction),
                  caption = "Dots show individual cluster means. Diamonds show the mean across cells.") +
    ggplot2::theme_minimal()

  list(summary = summary, by_cluster = by_cluster, plot = plot)
}
