#' Simulate mixed pseudobulk samples from a Seurat object
#'
#' Cells are sampled within groups and their raw counts are summed to create
#' pseudobulk samples. Mixture proportions can be supplied explicitly or drawn
#' from a Dirichlet distribution.
#'
#' @param object A Seurat object.
#' @param group_labels Either the name of a column in `object[[]]`, or a vector
#'   of group labels. A label vector may be named by cell barcode; otherwise it
#'   must be in the same order as the columns of the selected count matrix.
#' @param proportions Optional numeric vector or matrix. A vector describes one
#'   mixture and should be named by group. A matrix has pseudobulk samples in
#'   rows and groups in columns. Rows are normalized to sum to one. Missing
#'   groups are assigned proportion zero.
#' @param n_samples Number of mixtures to simulate when `proportions = NULL`.
#' @param cells_per_sample Number of cells in each pseudobulk sample. Supply one
#'   integer or a vector of length `n_samples`.
#' @param concentration Dirichlet concentration used when
#'   `proportions = NULL`. Supply one positive number, one number per group, or
#'   a named vector. Values below 1 produce more extreme mixtures; values above
#'   1 produce more balanced mixtures.
#' @param assay Seurat assay containing raw counts. The default is the object's
#'   default assay.
#' @param layer Layer (Seurat v5) or slot (Seurat v4) containing raw counts.
#' @param replace Whether cells may be sampled more than once within a
#'   pseudobulk sample. Separate pseudobulk samples are always drawn
#'   independently.
#' @param seed Optional random seed.
#' @param sample_prefix Prefix used for generated pseudobulk sample names.
#' @param return_cell_ids Whether to return the sampled cell barcodes. With
#'   replacement, a barcode can occur more than once.
#'
#' @return A list containing:
#'   * `counts`: genes-by-pseudobulk raw-count matrix;
#'   * `requested_proportions`: normalized target proportions;
#'   * `realized_proportions`: proportions after integer cell allocation;
#'   * `sampled_cell_counts`: allocated cells per group;
#'   * `sampled_cells`: sampled barcodes, or `NULL`;
#'   * `settings`: assay, layer, groups, and sampling settings.
#'
#' @examples
#' # Random, relatively heterogeneous mixtures:
#' # sim <- simulate_seurat_pseudobulk(
#' #   object = seu,
#' #   group_labels = "seurat_clusters",
#' #   n_samples = 20,
#' #   cells_per_sample = 1000,
#' #   concentration = 0.5,
#' #   seed = 1
#' # )
#' # bulk_counts <- sim$counts
#'
#' # User-defined mixtures; columns must match group labels:
#' # mixes <- rbind(
#' #   mix_1 = c(B = 0.70, T = 0.20, Mono = 0.10),
#' #   mix_2 = c(B = 0.10, T = 0.30, Mono = 0.60)
#' # )
#' # sim <- simulate_seurat_pseudobulk(
#' #   seu, "cell_type", proportions = mixes,
#' #   cells_per_sample = 2000, seed = 1
#' # )
simulate_seurat_pseudobulk <- function(
    object,
    group_labels,
    proportions = NULL,
    n_samples = 10L,
    cells_per_sample = 1000L,
    concentration = 1,
    assay = NULL,
    layer = "counts",
    replace = TRUE,
    seed = NULL,
    sample_prefix = "pseudobulk",
    return_cell_ids = FALSE) {

  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required.", call. = FALSE)
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  if (!is.character(layer) || length(layer) != 1L || !nzchar(layer)) {
    stop("`layer` must be one non-empty character string.", call. = FALSE)
  }
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("`replace` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_cell_ids) ||
      length(return_cell_ids) != 1L || is.na(return_cell_ids)) {
    stop("`return_cell_ids` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(sample_prefix) ||
      length(sample_prefix) != 1L || !nzchar(sample_prefix)) {
    stop("`sample_prefix` must be one non-empty character string.",
         call. = FALSE)
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }

  # Seurat v5 stores expression in layers; Seurat v4 uses assay slots.
  selected_layer <- layer
  counts <- tryCatch(
    {
      if (utils::packageVersion("SeuratObject") >= "5.0.0") {
        available_layers <- SeuratObject::Layers(object[[assay]])
        if (!layer %in% available_layers) {
          split_matches <- available_layers[
            startsWith(available_layers, paste0(layer, "."))
          ]
          if (length(split_matches) > 1L) {
            stop(
              "Multiple layers match '", layer, "': ",
              paste(split_matches, collapse = ", "),
              ". Join them with SeuratObject::JoinLayers() before simulation."
            )
          }
          if (length(split_matches) == 1L) {
            selected_layer <- split_matches
          }
        }
        SeuratObject::LayerData(
          object = object,
          assay = assay,
          layer = selected_layer
        )
      } else {
        SeuratObject::GetAssayData(object = object, assay = assay, slot = layer)
      }
    },
    error = function(e) {
      stop(
        "Could not read assay '", assay, "', layer/slot '", layer, "': ",
        conditionMessage(e),
        ". For a Seurat v5 assay split into multiple count layers, join the ",
        "layers first or specify one exact layer.",
        call. = FALSE
      )
    }
  )

  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop("The selected count matrix is empty.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("The count matrix must have gene and cell names.", call. = FALSE)
  }
  if (anyDuplicated(colnames(counts))) {
    stop("Cell names in the count matrix must be unique.", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("The selected layer contains negative values; use raw counts.",
         call. = FALSE)
  }

  metadata <- object[[]]
  if (is.character(group_labels) && length(group_labels) == 1L &&
      group_labels %in% colnames(metadata)) {
    labels <- metadata[[group_labels]]
    names(labels) <- rownames(metadata)
    labels <- labels[match(colnames(counts), names(labels))]
  } else {
    labels <- group_labels
    if (!is.null(names(labels))) {
      cell_match <- match(colnames(counts), names(labels))
      if (anyNA(cell_match)) {
        stop("Named `group_labels` do not cover every cell in the count matrix.",
             call. = FALSE)
      }
      labels <- labels[cell_match]
    } else if (length(labels) != ncol(counts)) {
      stop(
        "`group_labels` must be a metadata column name, a named vector ",
        "covering all cells, or a vector of length ncol(counts).",
        call. = FALSE
      )
    }
  }

  keep <- !is.na(labels) & nzchar(as.character(labels))
  if (!all(keep)) {
    warning(sum(!keep), " cells with missing/empty group labels were removed.",
            call. = FALSE)
    counts <- counts[, keep, drop = FALSE]
    labels <- labels[keep]
  }
  if (ncol(counts) == 0L) {
    stop("No cells remain after removing missing group labels.", call. = FALSE)
  }

  if (is.factor(labels)) {
    groups <- levels(droplevels(labels))
  } else {
    groups <- unique(as.character(labels))
  }
  labels <- as.character(labels)
  n_groups <- length(groups)
  if (n_groups < 2L) {
    stop("At least two non-empty groups are required to create mixtures.",
         call. = FALSE)
  }

  validate_positive_integer <- function(x, name) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
        any(x < 1) || any(x != floor(x))) {
      stop("`", name, "` must contain positive integers.", call. = FALSE)
    }
    as.integer(x)
  }

  if (is.null(proportions)) {
    n_samples <- validate_positive_integer(n_samples, "n_samples")
    if (length(n_samples) != 1L) {
      stop("`n_samples` must have length one.", call. = FALSE)
    }

    if (!is.numeric(concentration) || anyNA(concentration) ||
        any(!is.finite(concentration)) || any(concentration <= 0)) {
      stop("`concentration` must contain positive finite numbers.",
           call. = FALSE)
    }
    if (!is.null(names(concentration))) {
      if (!all(groups %in% names(concentration))) {
        stop("Named `concentration` must contain every group.", call. = FALSE)
      }
      concentration <- concentration[groups]
    } else if (length(concentration) == 1L) {
      concentration <- rep(concentration, n_groups)
    } else if (length(concentration) != n_groups) {
      stop("`concentration` must have length 1 or one value per group.",
           call. = FALSE)
    }

    if (!is.null(seed)) {
      set.seed(seed)
    }
    gamma_draws <- matrix(
      stats::rgamma(n_samples * n_groups,
                    shape = rep(concentration, each = n_samples)),
      nrow = n_samples,
      ncol = n_groups,
      dimnames = list(NULL, groups)
    )
    gamma_sums <- rowSums(gamma_draws)
    if (any(!is.finite(gamma_sums)) || any(gamma_sums <= 0)) {
      stop("Dirichlet sampling failed; use less extreme concentration values.",
           call. = FALSE)
    }
    requested <- gamma_draws / gamma_sums
  } else {
    if (is.numeric(proportions) && is.null(dim(proportions))) {
      supplied_names <- names(proportions)
      requested <- matrix(
        as.numeric(proportions),
        nrow = 1L,
        dimnames = list(NULL, supplied_names)
      )
    } else {
      if (is.data.frame(proportions) &&
          !all(vapply(proportions, is.numeric, logical(1)))) {
        stop("Every column of `proportions` must be numeric.", call. = FALSE)
      }
      requested <- as.matrix(proportions)
      if (!is.numeric(requested)) {
        stop("`proportions` must be numeric.", call. = FALSE)
      }
    }
    if (nrow(requested) < 1L || ncol(requested) < 1L) {
      stop("`proportions` cannot be empty.", call. = FALSE)
    }

    if (is.null(colnames(requested))) {
      if (ncol(requested) != n_groups) {
        stop(
          "Unnamed `proportions` must have one column per group; named ",
          "columns may contain a subset of groups.",
          call. = FALSE
        )
      }
      colnames(requested) <- groups
    } else {
      if (anyDuplicated(colnames(requested))) {
        stop("Column names of `proportions` must be unique.", call. = FALSE)
      }
      unknown_groups <- setdiff(colnames(requested), groups)
      if (length(unknown_groups) > 0L) {
        stop(
          "Unknown groups in `proportions`: ",
          paste(unknown_groups, collapse = ", "),
          call. = FALSE
        )
      }
      complete_requested <- matrix(
        0,
        nrow = nrow(requested),
        ncol = n_groups,
        dimnames = list(rownames(requested), groups)
      )
      complete_requested[, colnames(requested)] <- requested
      requested <- complete_requested
    }

    if (anyNA(requested) || any(!is.finite(requested)) ||
        any(requested < 0)) {
      stop("`proportions` must be finite, non-negative, and non-missing.",
           call. = FALSE)
    }
    requested_sums <- rowSums(requested)
    if (any(requested_sums <= 0)) {
      stop("Every row of `proportions` must have a positive sum.",
           call. = FALSE)
    }
    requested <- requested / requested_sums
    n_samples <- nrow(requested)

    if (!is.null(seed)) {
      set.seed(seed)
    }
  }

  sample_ids <- rownames(requested)
  if (is.null(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    sample_ids <- paste0(sample_prefix, "_", seq_len(n_samples))
  }
  rownames(requested) <- sample_ids

  cells_per_sample <- validate_positive_integer(
    cells_per_sample, "cells_per_sample"
  )
  if (length(cells_per_sample) == 1L) {
    cells_per_sample <- rep(cells_per_sample, n_samples)
  } else if (length(cells_per_sample) != n_samples) {
    stop(
      "`cells_per_sample` must have length one or one value per pseudobulk sample.",
      call. = FALSE
    )
  }

  allocated <- t(vapply(
    seq_len(n_samples),
    function(i) {
      as.integer(stats::rmultinom(
        n = 1L,
        size = cells_per_sample[i],
        prob = requested[i, ]
      ))
    },
    integer(n_groups)
  ))
  dimnames(allocated) <- list(sample_ids, groups)

  if (!replace) {
    available <- table(factor(labels, levels = groups))
    too_many <- allocated > matrix(
      as.integer(available),
      nrow = n_samples,
      ncol = n_groups,
      byrow = TRUE
    )
    if (any(too_many)) {
      first_bad <- which(too_many, arr.ind = TRUE)[1L, ]
      stop(
        "Cannot sample ", allocated[first_bad[1], first_bad[2]],
        " cells from group '", groups[first_bad[2]], "' for sample '",
        sample_ids[first_bad[1]], "' without replacement; only ",
        available[first_bad[2]], " cells are available.",
        call. = FALSE
      )
    }
  }

  pseudobulk_counts <- matrix(
    0,
    nrow = nrow(counts),
    ncol = n_samples,
    dimnames = list(rownames(counts), sample_ids)
  )
  sampled_cells <- if (return_cell_ids) {
    stats::setNames(vector("list", n_samples), sample_ids)
  } else {
    NULL
  }

  group_pools <- lapply(groups, function(g) which(labels == g))
  names(group_pools) <- groups

  for (i in seq_len(n_samples)) {
    selected <- integer(0)
    for (k in seq_len(n_groups)) {
      number_to_draw <- allocated[i, k]
      if (number_to_draw == 0L) {
        next
      }
      pool <- group_pools[[k]]
      selected <- c(
        selected,
        pool[sample.int(length(pool), number_to_draw, replace = replace)]
      )
    }

    pseudobulk_counts[, i] <- Matrix::rowSums(
      counts[, selected, drop = FALSE]
    )
    if (return_cell_ids) {
      sampled_cells[[i]] <- colnames(counts)[selected]
    }
  }

  realized <- allocated / rowSums(allocated)


  list(
    counts = pseudobulk_counts,
    requested_proportions = requested[,names(proportions),drop = F],
    realized_proportions = realized[,names(proportions),drop = F],
    sampled_cell_counts = allocated[,names(proportions),drop = F],
    sampled_cells = sampled_cells,
    settings = list(
      assay = assay,
      layer = selected_layer,
      groups = groups,
      cells_per_sample = cells_per_sample,
      replace = replace,
      seed = seed,
      concentration = if (is.null(proportions)) concentration else NULL
    )
  )
}
