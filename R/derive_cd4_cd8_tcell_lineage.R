#' Derive CD4 or CD8 T-cell lineage from expression and module scores
#'
#' Assigns CD4/CD8 lineage labels in two stages. First, expression of
#' \code{CD4}, \code{CD8A}, and \code{CD8B} is averaged within each group and
#' used to generate high-confidence seed labels. Second, markers derived from
#' the seed populations are converted to per-cell UCell scores and used to
#' classify cells with initially mixed or unknown lineage.
#'
#' @param obj A Seurat object containing normalized expression values and the
#'   grouping variable specified by \code{group}.
#'
#' @param group A single character string naming a column in
#'   \code{obj@meta.data}. Expression is averaged within this grouping variable.
#'   This will typically identify clonotypes, clones, or another collection of
#'   cells expected to share T-cell lineage.
#'
#' @param min_expression Non-negative numeric scalar giving the minimum average
#'   expression required for the initial lineage assignment. Groups for which
#'   both CD4 and combined CD8 expression are less than or equal to this value
#'   receive the seed label \code{"unknown"}. The appropriate value depends on
#'   the assay and expression scale. Default is \code{0}.
#'
#' @param cd4cd8_ratio Numeric scalar greater than \code{1} specifying the
#'   minimum CD4-to-CD8 or CD8-to-CD4 expression ratio required for an initial
#'   lineage assignment. For example, a value of \code{3} requires one lineage
#'   signal to be at least three times the other. Groups below this separation
#'   are initially labeled \code{"mix"}. Default is \code{3}.
#'
#' @param colname A single character string naming the metadata column in which
#'   the initial group-level lineage labels are stored. Final per-cell
#'   predictions are stored in a column named
#'   \code{paste0(colname, "_predict")}. Default is \code{"Tlin"}.
#'
#' @param ncores Positive integer giving the number of CPU cores supplied to
#'   \code{\link[UCell]{AddModuleScore_UCell}}. Default is \code{6}.
#'
#' @details
#' The initial group-level assignment uses
#' \deqn{\log_2\left(\frac{CD4 + \epsilon}{CD8 + \epsilon}\right),}
#' where \eqn{CD8 = (CD8A + CD8B)/2} and \eqn{\epsilon = 10^{-4}}.
#' Groups are assigned \code{"CD4"} or \code{"CD8"} when the corresponding
#' expression ratio exceeds \code{cd4cd8_ratio}. Groups with insufficient
#' expression are labeled \code{"unknown"}, and the remaining groups are
#' labeled \code{"mix"}.
#'
#' Differential markers are then identified between the confidently seeded CD4
#' and CD8 populations. After removing \code{CD4}, \code{CD8A}, and
#' \code{CD8B} from their corresponding marker sets, the remaining markers are
#' used to calculate per-cell UCell scores.
#'
#' Score thresholds and score-difference margins are estimated from the 20th
#' percentile of the corresponding seed population. A minimum score-difference
#' margin of \code{0.02} is enforced. The final per-cell labels have the
#' following interpretations:
#'
#' \describe{
#'   \item{\code{"CD4"}}{The CD4 score passes its threshold and exceeds the CD8
#'     score by the required margin.}
#'   \item{\code{"CD8"}}{The CD8 score passes its threshold and exceeds the CD4
#'     score by the required margin.}
#'   \item{\code{"mix"}}{Both scores pass their thresholds, but neither lineage
#'     has the required advantage.}
#'   \item{\code{"ambiguous"}}{Some lineage signal is present, but the evidence
#'     is insufficient for CD4, CD8, or mixed classification.}
#'   \item{\code{"unknown"}}{Neither score passes its threshold, or one or both
#'     scores are missing.}
#' }
#'
#' High-confidence initial CD4 and CD8 labels are preserved in the prediction
#' column. Only cells initially labeled \code{"mix"} or \code{"unknown"} are
#' reassigned from their UCell scores.
#'
#' At least 20 initially labeled CD4 cells and 20 initially labeled CD8 cells
#' are required to estimate the score thresholds. The function stops if either
#' reference population is too small.
#'
#' The complete per-cell classification table, including scores, thresholds,
#' support flags, and intermediate assignments, is stored in
#' \code{obj@misc$cd4_cd8_classification}.
#'
#' @return The input Seurat object with:
#' \itemize{
#'   \item group-level seed labels in the metadata column specified by
#'     \code{colname};
#'   \item final per-cell predictions in
#'     \code{paste0(colname, "_predict")};
#'   \item UCell score columns added by
#'     \code{\link[UCell]{AddModuleScore_UCell}}; and
#'   \item the complete classification table in
#'     \code{obj@misc$cd4_cd8_classification}.
#' }
#'
#' @note
#' Marker selection and score calibration are performed on the same object.
#' The resulting labels should therefore be treated as data-driven annotations,
#' not as independently validated predictions. Where possible, inspect known
#' lineage markers and validate the thresholds using independent annotations
#' or held-out samples.
#'
#'   | Label       | Interpretation                                      |
#'   | ----------- | --------------------------------------------------- |
#'   | `CD4`       | Sufficient CD4 score and clear CD4 advantage        |
#'   | `CD8`       | Sufficient CD8 score and clear CD8 advantage        |
#'   | `mix`       | Both signatures high, but neither clearly dominates |
#'   | `ambiguous` | Some signal exists, but evidence is insufficient    |
#'   | `unknown`   | Neither lineage signature is sufficiently high      |
#'
#' The function depends on the project-specific helper functions
#' \code{avg_expression()}, \code{join_meta_data()}, \code{find_all_marker()},
#' \code{subset2()}, and \code{get_data()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' t_cells <- derive_cd4_cd8_tcell_lineage(
#'   obj = t_cells,
#'   group = "clonotype_id",
#'   min_expression = 0,
#'   cd4cd8_ratio = 3,
#'   colname = "Tlin",
#'   ncores = 6
#' )
#'
#' table(t_cells$Tlin, useNA = "ifany")
#' table(t_cells$Tlin_predict, useNA = "ifany")
#'
#' classification <- t_cells@misc$cd4_cd8_classification
#' }
derive_cd4_cd8_tcell_lineage <- function(obj,
                                         group,
                                         min_expression = 0,
                                         cd4cd8_ratio = 3,
                                         colname = "Tlin",
                                         ncores = 6,
                                         skip_module_scores = T) {


  obj <- infer_cd4_cd8_lineage(
    obj = obj,
    group = group,
    min_expression = min_expression,
    cd4cd8_ratio = cd4cd8_ratio,
    colname = colname,
    run_missing_group_single = T)

  obj <- vote_cd4_cd8_lineage_snn(
    object = obj,
    label_col = colname)

  if (!skip_module_scores) {
    obj <- add_cd4_cd8_scores(obj = obj,
                              label_col = colname)
  }


  # tt <- get_data(obj, c("CD4score_UCell", "CD8score_UCell", colname, group), reduction = NULL, try_df = T) |>
  #   dplyr::select(id, dplyr::ends_with("UCell"), dplyr::all_of(c(group, colname))) |>
  #   dplyr::mutate(diff = CD4score_UCell-CD8score_UCell)
  #
  # # Reference labels produced from the initial CD4/CD8 expression rule
  # reference_cd4 <- tt[[colname]] == "CD4"
  # reference_cd8 <- tt[[colname]] == "CD8"
  #
  # if (sum(reference_cd4, na.rm = TRUE) < 20L ||
  #     sum(reference_cd8, na.rm = TRUE) < 20L) {
  #   stop("Too few confidently seeded CD4 or CD8 cells to calibrate thresholds.")
  # }
  #
  # # Minimum lineage-specific score:
  # # 20th percentile means approximately 80% of the corresponding
  # # confident reference population passes this threshold.
  # cd4_threshold <- unname(stats::quantile(
  #   tt$CD4score_UCell[reference_cd4],
  #   probs = 0.20,
  #   na.rm = TRUE
  # ))
  #
  # cd8_threshold <- unname(stats::quantile(
  #   tt$CD8score_UCell[reference_cd8],
  #   probs = 0.20,
  #   na.rm = TRUE
  # ))
  #
  # # Score advantage expected in confident reference cells.
  # # These are directional:
  # # CD4 references should have CD4score - CD8score > 0.
  # # CD8 references should have CD8score - CD4score > 0.
  # cd4_margin <- unname(stats::quantile(
  #   tt$CD4score_UCell[reference_cd4] -
  #     tt$CD8score_UCell[reference_cd4],
  #   probs = 0.20,
  #   na.rm = TRUE
  # ))
  #
  # cd8_margin <- unname(stats::quantile(
  #   tt$CD8score_UCell[reference_cd8] -
  #     tt$CD4score_UCell[reference_cd8],
  #   probs = 0.20,
  #   na.rm = TRUE
  # ))
  #
  # # A negative calibrated margin would allow the wrong score to win.
  # # Clamp it to zero, or preferably to a user-supplied positive floor.
  # min_margin <- 0.02
  #
  # cd4_margin <- max(cd4_margin, min_margin)
  # cd8_margin <- max(cd8_margin, min_margin)
  #
  # prediction_col <- paste0(colname, "_predict")
  #
  # #   | Label       | Interpretation                                      |
  # #   | ----------- | --------------------------------------------------- |
  # #   | `CD4`       | Sufficient CD4 score and clear CD4 advantage        |
  # #   | `CD8`       | Sufficient CD8 score and clear CD8 advantage        |
  # #   | `mix`       | Both signatures high, but neither clearly dominates |
  # #   | `ambiguous` | Some signal exists, but evidence is insufficient    |
  # #   | `unknown`   | Neither lineage signature is sufficiently high      |
  #
  #
  # tt <- tt |>
  #   dplyr::mutate(
  #     cd4_advantage = CD4score_UCell - CD8score_UCell,
  #     cd8_advantage = CD8score_UCell - CD4score_UCell,
  #
  #     cd4_supported =
  #       CD4score_UCell >= cd4_threshold &
  #       cd4_advantage >= cd4_margin,
  #
  #     cd8_supported =
  #       CD8score_UCell >= cd8_threshold &
  #       cd8_advantage >= cd8_margin,
  #
  #     both_high =
  #       CD4score_UCell >= cd4_threshold &
  #       CD8score_UCell >= cd8_threshold,
  #
  #     neither_high =
  #       CD4score_UCell < cd4_threshold &
  #       CD8score_UCell < cd8_threshold,
  #
  #     predicted_lineage = dplyr::case_when(
  #       is.na(CD4score_UCell) | is.na(CD8score_UCell) ~ "unknown",
  #       cd4_supported                                ~ "CD4",
  #       cd8_supported                                ~ "CD8",
  #       neither_high                                 ~ "unknown",
  #       both_high                                    ~ "mix",
  #       TRUE                                         ~ "ambiguous"
  #     ),
  #
  #     # Preserve the high-confidence seed calls and predict only
  #     # initially unknown/mixed cells.
  #     !!rlang::sym(prediction_col) := dplyr::if_else(
  #       .data[[colname]] %in% c("CD4", "CD8"),
  #       .data[[colname]],
  #       predicted_lineage
  #     )
  #   )
  #
  # # collapse clonotypes with different annotations to one lineage
  # # if scores are consistently higher for either lineage
  # tt2 <- tt |>
  #   dplyr::distinct(cl_name, Tlin_predict) |>
  #   dplyr::add_count(cl_name) |>
  #   dplyr::filter(n>1) |>
  #   tidyr::drop_na() |>
  #   dplyr::distinct(!!rlang::sym(group))
  #
  # for (i in tt2[[group]]) {
  #   if (mean(tt[which(tt[[group]] == i),"CD4score_UCell"] > tt[which(tt[[group]] == i),"CD8score_UCell"]) > 0.9) {
  #     tt[which(tt[[group]] == i),prediction_col] <- "CD4"
  #   } else if (mean(tt[which(tt[[group]] == i),"CD4score_UCell"] < tt[which(tt[[group]] == i),"CD8score_UCell"]) > 0.9) {
  #     tt[which(tt[[group]] == i),prediction_col] <- "CD8"
  #   } else {
  #     tt[which(tt[[group]] == i),prediction_col] <- "ambiguous"
  #   }
  # }
  #
  # obj@misc$cd4_cd8_classification <- tt
  # obj <- Seurat::AddMetaData(obj, tt |> tibble::column_to_rownames("id") |> dplyr::select(dplyr::all_of(prediction_col)))


  return(obj)
}



infer_cd4_cd8_lineage <- function(
    obj,
    group = NULL,
    min_expression = 0,
    cd4cd8_ratio = 3,
    colname = "cd4cd8lin",
    run_missing_group_single = T) {


  features <- c("CD4", "CD8A", "CD8B")
  pseudocount <- 1e-4

  classify_lineage <- function(df) {
    df |>
      dplyr::mutate(
        CD8 = (CD8A + CD8B) / 2,
        log_ratio = log2((CD4 + pseudocount) / (CD8 + pseudocount)),
        !!rlang::sym(colname) := dplyr::case_when(
          CD4 <= min_expression & CD8 <= min_expression ~ "unknown",
          log_ratio >= log2(cd4cd8_ratio)               ~ "CD4",
          log_ratio <= -log2(cd4cd8_ratio)              ~ "CD8",
          .default =                                      "mix"
        )
      )
  }

  if (!is.null(group)) {
    if (length(group) != 1L || !group %in% names(obj@meta.data)) {
      stop("`group` must be NULL or the name of a column in obj@meta.data.")
    }


    group_na <- is.na(obj@meta.data[[group]])
    message(round(mean(group_na)*100, 2), " % of cases in group column are NA.")

    group_nonna <- length(stats::na.omit(unique(obj@meta.data[[group]])))
    message(group_nonna, " unique non-NA cases in group column found.")

    # Classify using average expression within each group
    df <- avg_expression(obj, group = group, features = features)[[1]] |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column(group) |>
      classify_lineage()
    print(table(df[[colname]]))

    obj <- join_meta_data(obj, df[, c(group, colname)], by = group, verbose = F)

    # grouping column was NA: classify every cell independently
    if (run_missing_group_single) {
      message("running single mode with NA-rows in group column.")
      df <- get_data(obj, cells = cells2(obj)[which(group_na)], features, reduction = NULL, try_df = T) |>
        dplyr::filter(cells == 1) |>
        classify_lineage() |>
        dplyr::select(id, dplyr::all_of(colname))
      print(table(df[[colname]]))

      obj <- join_meta_data(obj, df = df, by = "id", cols = colname, verbose = F)
    } else {
      obj[[colname]][is.na(obj[[colname]][, 1])] <- "unknown"
    }

  } else {
    # No grouping column: classify every cell independently
    df <- get_data(obj, features, reduction = NULL, try_df = T)|>
      classify_lineage() |>
      dplyr::select(id, dplyr::all_of(colname)) |>
      tibble::column_to_rownames("id")
    print(table(df[[colname]]))

    obj <- Seurat::AddMetaData(obj, df)
  }



  return(obj)
}


#' Assign unknown T cells by weighted voting on a Seurat SNN graph
#'
#' Known labels are never changed. Only cells whose current label is in
#' `unknown` (or is NA) are eligible for assignment. SNN edge weights are used
#' as vote weights, and only known cells are allowed to vote.
#'
#' @return The Seurat object with four added metadata columns: the voted label,
#'   vote confidence, number of labeled neighbors, and fraction of neighboring
#'   SNN weight connected to labeled cells.
vote_cd4_cd8_lineage_snn <- function(
    object,
    label_col = "cd4cd8lin",
    graph_name = NULL,
    classes = c("CD4", "CD8"),
    unknown = "unknown",
    output_col = paste0(label_col, "_snn"),
    min_confidence = 0.8,
    min_labeled_neighbors = 3L,
    min_labeled_weight_fraction = 0,
    unresolved_label = "unknown"
) {
  stopifnot(
    length(label_col) == 1L,
    length(output_col) == 1L,
    length(classes) >= 2L,
    !anyDuplicated(classes),
    min_confidence >= 0,
    min_confidence <= 1,
    min_labeled_neighbors >= 0,
    min_labeled_weight_fraction >= 0,
    min_labeled_weight_fraction <= 1
  )

  meta <- object[[]]
  if (!label_col %in% colnames(meta)) {
    stop("Metadata column not found: ", label_col)
  }

  available_graphs <- SeuratObject::Graphs(object)
  if (is.null(graph_name)) {
    preferred <- paste0(SeuratObject::DefaultAssay(object), "_snn")
    snn_graphs <- grep("snn$", available_graphs, value = TRUE, ignore.case = TRUE)
    if (preferred %in% available_graphs) {
      graph_name <- preferred
    } else if (length(snn_graphs) == 1L) {
      graph_name <- snn_graphs
    } else {
      stop(
        "Could not choose one SNN graph. Supply `graph_name`; available graphs: ",
        paste(available_graphs, collapse = ", ")
      )
    }
  }
  if (!graph_name %in% available_graphs) {
    stop("Graph not found: ", graph_name)
  }

  labels <- as.character(meta[[label_col]])
  cells <- rownames(meta)
  known <- !is.na(labels) & labels %in% classes
  candidates <- is.na(labels) | labels %in% unknown

  if (!any(known)) stop("No known cells have labels in `classes`.")
  if (!any(candidates)) {
    warning("No unknown/NA candidate cells were found; labels are copied unchanged.")
  }

  W <- object[[graph_name]]
  if (is.null(rownames(W)) || !all(cells %in% rownames(W))) {
    stop("The SNN graph does not contain all cells in the Seurat object.")
  }
  W <- W[cells, cells, drop = FALSE]
  diag(W) <- 0
  W <- Matrix::drop0(W)

  # Sparse one-hot matrix: rows are cells, columns are the trusted classes.
  class_matrix <- Matrix::sparseMatrix(
    i = which(known),
    j = match(labels[known], classes),
    x = 1,
    dims = c(length(cells), length(classes)),
    dimnames = list(cells, classes)
  )

  # Each score is the sum of SNN weights from a cell to known cells of a class.
  scores <- as.matrix(W %*% class_matrix)
  labeled_weight <- rowSums(scores)
  total_weight <- Matrix::rowSums(W)
  labeled_neighbors <- as.numeric((W > 0) %*% as.numeric(known))

  best_score <- apply(scores, 1L, max)
  confidence <- ifelse(labeled_weight > 0, best_score / labeled_weight, NA_real_)
  labeled_weight_fraction <- ifelse(
    total_weight > 0,
    labeled_weight / total_weight,
    NA_real_
  )

  # Do not break exact ties arbitrarily.
  unique_winner <- rowSums(scores == best_score) == 1L & best_score > 0
  winner <- max.col(scores, ties.method = "first")
  voted_class <- classes[winner]

  assignable <- candidates &
    unique_winner &
    !is.na(confidence) & confidence >= min_confidence &
    labeled_neighbors >= min_labeled_neighbors &
    !is.na(labeled_weight_fraction) &
    labeled_weight_fraction >= min_labeled_weight_fraction

  final_label <- labels
  final_label[candidates] <- unresolved_label
  final_label[assignable] <- voted_class[assignable]

  result <- data.frame(
    final_label,
    confidence,
    labeled_neighbors,
    labeled_weight_fraction,
    row.names = cells,
    check.names = FALSE
  )
  colnames(result) <- c(
    output_col,
    paste0(output_col, "_confidence"),
    paste0(output_col, "_n_labeled_neighbors"),
    paste0(output_col, "_labeled_weight_fraction")
  )

  object <- SeuratObject::AddMetaData(object, metadata = result)

  print(composition_barplot(
    object,
    y = "abs",
    x_cat = label_col,
    fill_cat = output_col
  )$plot)

  return(object)
}



#' Add CD4 and CD8 Gene-Signature Scores
#'
#' Identifies marker genes distinguishing CD4 and CD8 cells and calculates
#' UCell module scores for each signature. Marker genes are restricted to
#' features detected in more than 30 percent of cells with an adjusted
#' p-value below `1e-5`. Up to 20 top markers per group are selected using
#' both `avg_log2FC` and `logFC`.
#'
#' The canonical marker genes `CD4`, `CD8A`, and `CD8B` are excluded from
#' their corresponding signatures. The resulting CD4 score, CD8 score, and
#' their difference are stored in the object's cell-level metadata as
#' `CD4score_UCell`, `CD8score_UCell`, and `diff_CD4_CD8_score`, respectively.
#'
#' @param obj A Seurat object containing normalized expression data and a
#'   cell-level metadata column with CD4/CD8 lineage labels.
#' @param label_col Character scalar giving the name of the metadata column
#'   containing the labels `"CD4"` and `"CD8"`. Defaults to `"cd4cd8lin"`.
#' @param ncores Positive integer specifying the number of CPU cores passed to
#'   [UCell::AddModuleScore_UCell()]. Defaults to `8`.
#'
#' @return The input Seurat object with `CD4score_UCell`, `CD8score_UCell`,
#'   and `diff_CD4_CD8_score` added to its cell-level metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- add_cd4_cd8_scores(
#'   object,
#'   label_col = "cd4cd8lin",
#'   ncores = 4
#' )
#'
#' head(object[[]][
#'   c("CD4score_UCell", "CD8score_UCell", "diff_CD4_CD8_score")
#' ])
#' }
add_cd4_cd8_scores <- function(obj,
                               label_col = "cd4cd8lin",
                               ncores = 8) {

  tmark <- find_all_marker(subset2(
    obj,
    subset = !!rlang::sym(label_col) %in% c("CD4", "CD8")
  ), meta_col = label_col) |>
    dplyr::filter(pct_in>30 & padj<1e-5)

  tmark <- dplyr::bind_rows(dplyr::slice_max(tmark, avg_log2FC, n = 20, by = group),
                            dplyr::slice_max(tmark, logFC, n = 20, by = group)) |>
    dplyr::distinct()
  tmark <- split(tmark$feature, tmark$group)

  message("running AddModuleScore_UCell.")
  obj <- UCell::AddModuleScore_UCell(obj = obj,
                                     features = list(CD4score = setdiff(tmark[["CD4"]], "CD4"),
                                                     CD8score = setdiff(tmark[["CD8"]], c("CD8A", "CD8B"))),
                                     ncores = ncores,
                                     force.gc = T)

  tt <- get_data(obj, c("CD4score_UCell", "CD8score_UCell"), reduction = NULL, try_df = T) |>
    dplyr::mutate(diff_CD4_CD8_score = CD4score_UCell-CD8score_UCell) |>
    dplyr::select(id, diff_CD4_CD8_score) |>
    tibble::column_to_rownames("id")

  obj <- Seurat::AddMetaData(obj, tt)

  return(obj)
}
