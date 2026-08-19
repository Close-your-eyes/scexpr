#' Split a Seurat Object into Cell Chunks
#'
#' Splits the cells of a Seurat object into chunks and returns one subsetted
#' Seurat object per chunk. Optionally reduces the object with
#' [Seurat::DietSeurat()] before splitting.
#'
#' @param obj A Seurat object.
#' @param chunks Integer. Number of chunks to create. Defaults to `2`.
#' @param size Optional integer specifying the desired number of cells per
#'   chunk. Passed to [brathering::split_chunks()].
#' @param shuffle Logical. Whether to randomly shuffle cells before splitting.
#'   Defaults to `FALSE`.
#' @param diet Logical. Whether to apply [Seurat::DietSeurat()] before
#'   splitting. Defaults to `FALSE`.
#' @param diet_args Named list of additional arguments passed to
#'   [Seurat::DietSeurat()]. The input object is supplied automatically as
#'   `object`. Defaults to retaining the RNA counts layer and dropping `misc`.
#'
#' @return A list of Seurat objects, with one object for each cell chunk.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' chunks <- split_seurat(pbmc_small, chunks = 2)
#'
#' chunks <- split_seurat(
#'   pbmc_small,
#'   chunks = 3,
#'   shuffle = TRUE,
#'   diet = TRUE
#' )
#' }
split_seurat <- function(obj,
                         chunks = 2,
                         size = NULL,
                         shuffle = F,
                         diet = F,
                         diet_args = list(misc = F, layers = "counts", assays = "RNA")) {

  if (diet) {
    diet_args <- c(list(object = obj), diet_args)
    obj <- Gmisc::fastDoCall(Seurat::DietSeurat, diet_args)
  }
  cells <- brathering::split_chunks(
    x = scexpr::cells2(obj),
    chunks = chunks,
    size = size,
    shuffle = shuffle
  )
  return(purrr::map(cells, ~subset2(obj, cells = .x)))
}
