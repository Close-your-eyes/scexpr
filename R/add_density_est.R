#' Add density estimates to cell metadata
#'
#' Computes density estimates from a dimensional reduction embedding and adds the
#' density values to the object's metadata.
#'
#' @param obj A Seurat object.
#' @param reduction Character. Name of the dimensional reduction to use.
#' @param name Character. Name of the metadata column to add. Defaults to
#'   `paste0(reduction, "_dens")`.
#' @param density_est_args List of arguments passed to
#'   `brathering::density_est()`. Defaults to `list(n = 200, h = 1)`.
#' @param type use raw or scaled density values?
#'
#' @returns The input object with an added metadata column containing density
#'   estimates.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- add_density_est(obj, reduction = "umap")
#' }
add_density_est <- function(obj,
                            reduction,
                            name = paste0(reduction, "_dens"),
                            density_est_args = list(n = 200, h = 1),
                            type = c("scaled", "raw")) {
  .ensure_package("brathering")

  if (missing(reduction)) {
    stop("reduction missing.")
  }
  if (!reduction %in% names(obj@reductions)) {
    stop("reduction not found in obj.")
  }
  type <- rlang::arg_match(type)
  density_est_args[["x"]] <- obj@reductions[[reduction]]@cell.embeddings
  obj@meta.data[[name]] <- do.call(brathering::density_est, density_est_args)[[type]]
  return(obj)
}
