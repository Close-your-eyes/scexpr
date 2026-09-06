#' Join Metadata to a Seurat Object
#'
#' Adds columns from a data frame to the cell-level metadata of a Seurat
#' object using a left join.
#'
#' @param obj A Seurat object.
#' @param df A data frame containing metadata to join.
#' @param by A character vector specifying the join columns. Passed to
#'   [dplyr::left_join()]. If `NULL`, the join uses all variables shared by
#'   both data frames.
#' @param rm_nonjoin remove non-joined columns
#'
#' @return The Seurat object with the joined data stored in its
#'   `meta.data` slot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' additional_metadata <- data.frame(
#'   sample_id = c("sample1", "sample2"),
#'   treatment = c("control", "treated")
#' )
#'
#' object <- join_meta_data(
#'   object,
#'   additional_metadata,
#'   by = "sample_id"
#' )
#' }
join_meta_data_legacy <- function(obj,
                                  df,
                                  by = NULL,
                                  rm_nonjoin = T) {


  original_n <- nrow(obj@meta.data)

  # Preserve join columns while removing columns that will be replaced.
  join_cols <- if (is.null(by)) {
    intersect(names(obj@meta.data), names(df))
  } else if (is.character(by)) {
    intersect(names(obj@meta.data), unname(by))
  } else {
    stop("by must be NULL or character.")
  }

  if (!length(join_cols)) {
    stop("no matching join columns found.")
  }

  message("join columns: ", paste(join_cols, collapse = ","))

  columns_to_remove <- setdiff(
    intersect(names(obj@meta.data), names(df)),
    join_cols)

  nonjoin_cols <- setdiff(names(df), join_cols)

  if (length(columns_to_remove)) {
    if (rm_nonjoin) {
      message("columns removed before join: ", paste(columns_to_remove, collapse = ","))
      obj@meta.data <- dplyr::select(obj@meta.data, -dplyr::any_of(columns_to_remove))
    } else {
      message("previously present columns in obj@meta.data and df: ", paste(columns_to_remove, collapse = ","))
    }
  }

  meta <- dplyr::left_join(obj@meta.data, df, by = join_cols)

  if (nrow(meta) != original_n) {
    stop(
      "`df` contains duplicate join keys that expanded the metadata from ",
      original_n, " to ", nrow(meta), " rows.",
      call. = FALSE
    )
  }

  check_na_col <- nonjoin_cols[1]
  if (length(check_na_col)) {
    if (anyNA(meta[[check_na_col]])) {
      message(sum(is.na(meta[[check_na_col]])), " rows with NA after join in ", check_na_col, " column.")
    }
  }


  rownames(meta) <- rownames(obj@meta.data)
  obj@meta.data <- meta
  return(obj)

}


#' Coalesce Selected Metadata into a Seurat Object
#'
#' Adds selected columns from a data frame to the cell-level metadata of a
#' Seurat object using `coalesce_join()`. Existing non-NA metadata values are
#' prioritized over values from `df`.
#'
#' @param obj A Seurat object.
#' @param df A data frame containing metadata to join.
#' @param by A character vector specifying the join columns. Passed to
#'   `coalesce_join()`.
#' @param cols A character vector specifying columns from `df` to add or
#'   coalesce. If `NULL`, all non-join columns in `df` are used.
#' @param verbose print all messages?
#'
#' @return The Seurat object with updated data stored in its `meta.data` slot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object <- join_meta_data(
#'   object,
#'   additional_metadata,
#'   by = "sample_id",
#'   cols = c("treatment", "batch")
#' )
#' }
join_meta_data <- function(obj,
                           df,
                           by = NULL,
                           cols = NULL,
                           verbose = T) {

  .ensure_package("brathering")

  original_meta <- obj@meta.data
  original_n <- nrow(original_meta)

  # Determine the left- and right-hand join columns.
  if (is.null(by)) {
    join_cols_left <- intersect(names(original_meta), names(df))
    join_cols_right <- join_cols_left
    join_by <- join_cols_left
  } else if (is.character(by)) {
    if (is.null(names(by))) {
      join_cols_left <- by
      join_cols_right <- by
    } else {
      join_cols_left <- names(by)
      join_cols_right <- unname(by)
    }

    join_by <- by
  } else {
    stop("`by` must be NULL or a character vector.", call. = FALSE)
  }

  if (!length(join_cols_left)) {
    stop("No matching join columns found.", call. = FALSE)
  }

  missing_left <- setdiff(join_cols_left, names(original_meta))
  missing_right <- setdiff(join_cols_right, names(df))

  if (length(missing_left)) {
    stop(
      "Join columns missing from `obj@meta.data`: ",
      paste(missing_left, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_right)) {
    stop(
      "Join columns missing from `df`: ",
      paste(missing_right, collapse = ", "),
      call. = FALSE
    )
  }

  # Select columns to add or coalesce.
  if (is.null(cols)) {
    cols <- setdiff(names(df), join_cols_right)
  } else {
    cols <- unique(as.character(cols))

    missing_cols <- setdiff(cols, names(df))

    if (length(missing_cols)) {
      stop(
        "Selected columns missing from `df`: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    cols <- setdiff(cols, join_cols_right)
  }

  if (verbose) {
    message(
      "join columns: ",
      paste(join_cols_left, collapse = ", ")
    )

    if (length(cols)) {
      message(
        "columns added or coalesced: ",
        paste(cols, collapse = ", ")
      )
    }
  }

  # Subsetting df ensures the original coalesce_join() only operates on
  # the selected metadata columns.
  df_selected <- dplyr::select(
    df,
    dplyr::all_of(unique(c(join_cols_right, cols)))
  )

  meta <- brathering::coalesce_join(
    x = original_meta,
    y = df_selected,
    by = join_by
  )

  if (nrow(meta) != original_n) {
    stop(
      "`df` contains duplicate join keys that expanded the metadata from ",
      original_n, " to ", nrow(meta), " rows.",
      call. = FALSE
    )
  }

  for (column in cols) {
    if (column %in% names(meta) && anyNA(meta[[column]])) {
      message(
        sum(is.na(meta[[column]])),
        " rows with NA after join in ",
        column,
        " column."
      )
    }
  }

  rownames(meta) <- rownames(original_meta)
  obj@meta.data <- meta

  return(obj)
}
