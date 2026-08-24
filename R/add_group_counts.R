#' Add Distinct Sample Counts to Seurat Metadata
#'
#' Calculates the number of distinct identifiers within each level of a
#' metadata variable and adds a formatted label to the metadata of a Seurat
#' object. For example, a disease value of `"Control"` with five distinct
#' samples becomes `"Control (n = 5)"`.
#'
#' @param obj A Seurat object.
#' @param group_col A single character string naming the metadata column used
#'   to define groups.
#' @param id_col A single character string naming the metadata column containing
#'   identifiers to count. Defaults to `"orig.ident"`.
#' @param output_col A single character string naming the new formatted metadata
#'   column. Defaults to `paste0(group_col, "2")`.
#'
#' @return The input Seurat object with `output_col` and `n` added to its
#'   metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sokid5 <- add_group_counts(
#'   obj = sokid5,
#'   group_col = "disease",
#'   id_col = "orig.ident",
#'   output_col = "disease2"
#' )
#' }
add_group_counts <- function(obj,
                             group_col,
                             id_col = "orig.ident",
                             output_col = paste0(group_col, "2")) {

  if (length(group_col)>1) {
    message("only one group_col allowed.")
    return(obj)
  }

  df <- obj@meta.data |>
    dplyr::summarise(
      n = dplyr::n_distinct(!!rlang::sym(id_col)),
      .by = !!rlang::sym(group_col)) |>
    dplyr::mutate(!!rlang::sym(output_col) := paste0(!!rlang::sym(group_col), " (n = ", n, ")"))
  print(df)
  obj <- join_meta_data(obj,
                        df = df,
                        by = group_col)


  return(obj)
}
