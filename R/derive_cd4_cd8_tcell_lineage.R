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
                                         ncores = 6) {

  # run w/o group? intention: group = clonotypes
  #
  #   df <- get_data(obj, c("CD4", "CD8A", "CD8B", group), reduction = NULL, try_df = T)
  #   df <- df |>
  #     dplyr::summarise(dplyr::across(c(CD4, CD8A, CD8B), sum), .by = !!rlang::sym(group)) |>
  #     tidyr::drop_na() |>
  #     dplyr::mutate(CD8 = (CD8A+CD8B)/2) |>
  #     dplyr::mutate(!!rlang::sym(colname) := dplyr::case_when(CD4 == 0 & CD8 == 0 ~ "unknown",
  #                                                             CD4/CD8>cd4cd8_ratio ~ "CD4",
  #                                                             CD8/CD4>cd4cd8_ratio ~ "CD8",
  #                                                             .default = "mix"))

  if (missing(group)) {
    stop("group (meta.data column) is needed.")
  }

  avgexpr <- avg_expression(obj, group = group, features = c("CD4", "CD8A", "CD8B"))[[1]] |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(group)
  pseudocount <- 1e-4
  df <- avgexpr |>
    dplyr::mutate(CD8 = (CD8A+CD8B)/2) |>
    dplyr::mutate(log_ratio = log2((CD4 + pseudocount) / (CD8 + pseudocount))) |>
    dplyr::mutate(!!rlang::sym(colname) := dplyr::case_when(
      CD4 <= min_expression & CD8 <= min_expression ~ "unknown",
      log_ratio >= log2(cd4cd8_ratio)           ~ "CD4",
      log_ratio <= -log2(cd4cd8_ratio)          ~ "CD8",
      TRUE                                        ~ "mix"))

  obj <- join_meta_data(obj, df[,c(group, colname)], by = group)
  obj@meta.data[[colname]][which(is.na(obj@meta.data[[colname]]))] <- "unknown"

  tmark <- find_all_marker(subset2(obj, subset = !!rlang::sym(colname) %in% c("CD4", "CD8")), meta_col = colname) |>
    dplyr::filter(pct_in>30 & padj<1e-5)
  tmark <- dplyr::bind_rows(dplyr::slice_max(tmark, avg_log2FC, n = 20, by = group),
                            dplyr::slice_max(tmark, logFC, n = 20, by = group)) |>
    dplyr::distinct()
  tmark <- split(tmark$feature, tmark$group)

  obj <- UCell::AddModuleScore_UCell(obj = obj,
                                     features = list(CD4score = setdiff(tmark[["CD4"]], "CD4"),
                                                     CD8score = setdiff(tmark[["CD8"]], c("CD8A", "CD8B"))),
                                     ncores = ncores,
                                     force.gc = T)

  tt <- get_data(obj, c("CD4score_UCell", "CD8score_UCell", colname, group), reduction = NULL, try_df = T) |>
    dplyr::select(id, dplyr::ends_with("UCell"), dplyr::all_of(c(group, colname))) |>
    dplyr::mutate(diff = CD4score_UCell-CD8score_UCell)

  # summ <- tt |> dplyr::summarise(
  #   CD4score_UCell = mean(CD4score_UCell),
  #   CD8score_UCell = mean(CD8score_UCell),
  #   diff = mean(diff),
  #   .by = !!rlang::sym(colname)
  # )
  #
  # avgcd8score <- summ |> dplyr::filter(!!rlang::sym(colname) == "CD8") |> dplyr::pull(CD8score_UCell)
  # avgcd4score <- summ |> dplyr::filter(!!rlang::sym(colname) == "CD4") |> dplyr::pull(CD4score_UCell)
  #
  # cd4thresh <- quantile(
  #   tt$CD4score_UCell[tt[[colname]] == "CD4"],
  #   0.20,
  #   na.rm = TRUE
  # )
  #
  # cd8thresh <- quantile(
  #   tt$CD8score_UCell[tt[[colname]] == "CD8"],
  #   0.20,
  #   na.rm = TRUE
  # )
  #
  # cd8lm <- stats::lm(diff~CD8score_UCell, dplyr::filter(tt, !!rlang::sym(colname) == "CD8"))
  # cd4lm <- stats::lm(diff~CD4score_UCell, dplyr::filter(tt, !!rlang::sym(colname) == "CD4"))
  #
  # tt$cd8avgdiff <- stats::predict(cd8lm, newdata = tt)
  # tt$cd4avgdiff <- stats::predict(cd4lm, newdata = tt)
  #
  # if ("unknown" %in% tt[[colname]]) {
  #   # colname2 <- paste0(colname, "_predict")
  #   # tt[[colname2]] <- tt[[colname]]
  #   tt <- split(tt, tt[[colname]])
  #   tt$unknown <- tt$unknown |>
  #     dplyr::mutate(!!rlang::sym(colname) := dplyr::case_when(CD8score_UCell<=cd8thresh & (CD4score_UCell>cd4thresh | diff>cd4avgdiff)~"CD4",
  #                                                             CD8score_UCell>cd8thresh & (CD4score_UCell<=cd4thresh | diff>cd8avgdiff)~"CD8",
  #                                                             CD8score_UCell<=cd8thresh & (CD4score_UCell>cd4thresh | diff>cd4avgdiff) & CD8score_UCell>cd8thresh & (CD4score_UCell<=cd4thresh | diff>cd8avgdiff)~"mix",
  #                                                             CD8score_UCell>cd8thresh & CD4score_UCell>cd4thresh~"mix",
  #                                                             CD4score_UCell<=0.2*cd4thresh & CD8score_UCell<=0.2*cd8thresh~"unknown",
  #                                                             .default = !!rlang::sym(colname)))
  #   tt <- dplyr::bind_rows(tt)
  # }
  # tt <- tt |> tibble::column_to_rownames("id")

  # Reference labels produced from the initial CD4/CD8 expression rule
  reference_cd4 <- tt[[colname]] == "CD4"
  reference_cd8 <- tt[[colname]] == "CD8"

  if (sum(reference_cd4, na.rm = TRUE) < 20L ||
      sum(reference_cd8, na.rm = TRUE) < 20L) {
    stop("Too few confidently seeded CD4 or CD8 cells to calibrate thresholds.")
  }

  # Minimum lineage-specific score:
  # 20th percentile means approximately 80% of the corresponding
  # confident reference population passes this threshold.
  cd4_threshold <- unname(stats::quantile(
    tt$CD4score_UCell[reference_cd4],
    probs = 0.20,
    na.rm = TRUE
  ))

  cd8_threshold <- unname(stats::quantile(
    tt$CD8score_UCell[reference_cd8],
    probs = 0.20,
    na.rm = TRUE
  ))

  # Score advantage expected in confident reference cells.
  # These are directional:
  # CD4 references should have CD4score - CD8score > 0.
  # CD8 references should have CD8score - CD4score > 0.
  cd4_margin <- unname(stats::quantile(
    tt$CD4score_UCell[reference_cd4] -
      tt$CD8score_UCell[reference_cd4],
    probs = 0.20,
    na.rm = TRUE
  ))

  cd8_margin <- unname(stats::quantile(
    tt$CD8score_UCell[reference_cd8] -
      tt$CD4score_UCell[reference_cd8],
    probs = 0.20,
    na.rm = TRUE
  ))

  # A negative calibrated margin would allow the wrong score to win.
  # Clamp it to zero, or preferably to a user-supplied positive floor.
  min_margin <- 0.02

  cd4_margin <- max(cd4_margin, min_margin)
  cd8_margin <- max(cd8_margin, min_margin)

  prediction_col <- paste0(colname, "_predict")

  #   | Label       | Interpretation                                      |
  #   | ----------- | --------------------------------------------------- |
  #   | `CD4`       | Sufficient CD4 score and clear CD4 advantage        |
  #   | `CD8`       | Sufficient CD8 score and clear CD8 advantage        |
  #   | `mix`       | Both signatures high, but neither clearly dominates |
  #   | `ambiguous` | Some signal exists, but evidence is insufficient    |
  #   | `unknown`   | Neither lineage signature is sufficiently high      |


  tt <- tt |>
    dplyr::mutate(
      cd4_advantage = CD4score_UCell - CD8score_UCell,
      cd8_advantage = CD8score_UCell - CD4score_UCell,

      cd4_supported =
        CD4score_UCell >= cd4_threshold &
        cd4_advantage >= cd4_margin,

      cd8_supported =
        CD8score_UCell >= cd8_threshold &
        cd8_advantage >= cd8_margin,

      both_high =
        CD4score_UCell >= cd4_threshold &
        CD8score_UCell >= cd8_threshold,

      neither_high =
        CD4score_UCell < cd4_threshold &
        CD8score_UCell < cd8_threshold,

      predicted_lineage = dplyr::case_when(
        is.na(CD4score_UCell) | is.na(CD8score_UCell) ~ "unknown",
        cd4_supported                                ~ "CD4",
        cd8_supported                                ~ "CD8",
        neither_high                                 ~ "unknown",
        both_high                                    ~ "mix",
        TRUE                                         ~ "ambiguous"
      ),

      # Preserve the high-confidence seed calls and predict only
      # initially unknown/mixed cells.
      !!rlang::sym(prediction_col) := dplyr::if_else(
        .data[[colname]] %in% c("CD4", "CD8"),
        .data[[colname]],
        predicted_lineage
      )
    )

  # collapse clonotypes with different annotations to one lineage
  # if scores are consistently higher for either lineage
  tt2 <- tt |>
    dplyr::distinct(cl_name, Tlin_predict) |>
    dplyr::add_count(cl_name) |>
    dplyr::filter(n>1) |>
    tidyr::drop_na() |>
    dplyr::distinct(!!rlang::sym(group))

  for (i in tt2[[group]]) {
    if (mean(tt[which(tt[[group]] == i),"CD4score_UCell"] > tt[which(tt[[group]] == i),"CD8score_UCell"]) > 0.9) {
      tt[which(tt[[group]] == i),prediction_col] <- "CD4"
    } else if (mean(tt[which(tt[[group]] == i),"CD4score_UCell"] < tt[which(tt[[group]] == i),"CD8score_UCell"]) > 0.9) {
      tt[which(tt[[group]] == i),prediction_col] <- "CD8"
    } else {
      tt[which(tt[[group]] == i),prediction_col] <- "ambiguous"
    }
  }

  obj@misc$cd4_cd8_classification <- tt
  obj <- Seurat::AddMetaData(obj, tt |> tibble::column_to_rownames("id") |> dplyr::select(dplyr::all_of(prediction_col)))
  return(obj)
}
