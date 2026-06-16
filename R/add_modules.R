#' Add UCell module scores to a Seurat object
#'
#' Computes UCell scores for a collection of gene modules and adds the
#' resulting scores to the metadata of a Seurat object.
#'
#' Gene modules are filtered by size before scoring, processed in chunks to
#' reduce memory usage, and optionally filtered based on the standard
#' deviation of their scores across cells.
#'
#' @param obj A \code{\link[Seurat:Seurat-class]{Seurat}} object.
#' @param modules A named list of character vectors containing gene symbols
#'   for each module. Module names are used as prefixes for the generated
#'   metadata columns.
#' @param minlen Minimum number of genes required for a module to be scored.
#'   Modules with fewer genes are excluded. Default is \code{10}.
#' @param maxlen Maximum number of genes allowed in a module. Modules with
#'   more genes are excluded. Default is \code{500}.
#' @param chunksize Number of modules to process simultaneously in each call
#'   to \code{\link[UCell:AddModuleScore_UCell]{AddModuleScore_UCell}}.
#'   Reducing this value may decrease memory usage. Default is \code{300}.
#' @param min_sd Minimum standard deviation required for a module score to be
#'   retained. Modules with score variability below this threshold are
#'   discarded. Default is \code{0.04}.
#' @param ncores Number of CPU cores to use for UCell score computation.
#'   Passed to \code{UCell::AddModuleScore_UCell()}. Default is \code{8}.
#' @param ... Additional arguments passed to
#'   \code{\link[UCell:AddModuleScore_UCell]{AddModuleScore_UCell}}.
#' @param downsample_explore downsaple size to find relevant sets in a first pass
#'
#' @return A \code{\link[Seurat:Seurat-class]{Seurat}} object with UCell module
#'   scores added to the metadata.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters modules according to \code{minlen} and \code{maxlen}.
#'   \item Splits the remaining modules into chunks of size
#'     \code{chunksize}.
#'   \item Computes UCell scores for each chunk.
#'   \item Removes module scores with a standard deviation below
#'     \code{min_sd}.
#'   \item Adds the retained scores to the Seurat object's metadata.
#' }
#'
#' Large numbers of modules can require substantial memory. Adjust
#' \code{chunksize} and \code{ncores} based on available computational
#' resources.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' modules <- list(
#'   Tcell = c("CD3D", "CD3E", "TRBC1", "IL7R"),
#'   Myeloid = c("LYZ", "S100A8", "S100A9", "FCN1")
#' )
#'
#' obj <- add_modules(
#'   obj = obj,
#'   modules = modules,
#'   minlen = 3,
#'   ncores = 4
#' )
#' }
add_modules <- function(obj,
                        modules,
                        minlen = 10,
                        maxlen = 500,
                        chunksize = 1000,
                        min_sd = 0.04,
                        ncores = 8,
                        downsample_explore = 0.05,
                        ...) {

  # sds <- apply(meta, 2, sd)
  # plot(sort(sds))
  #
  # maxs <- apply(meta, 2, max)
  # plot(sort(maxs))
  #
  # frac <- colSums(meta>0)/nrow(meta)
  # plot(sort(frac))
  #
  # means <- apply(meta, 2, mean)
  # plot(sort(means))
  # kw <- matrixTests::col_kruskalwallis(meta, so$MajorCluster2)
  # brathering::plot2(data.frame(means = means, frac = frac))

  # use kruskall for filtering? kw <- matrixTests::col_kruskalwallis(meta, so$MajorCluster2)

  stopifnot(inherits(obj, "Seurat"))
  stopifnot(is.list(modules))
  stopifnot(is.numeric(minlen), minlen > 0)
  if (is.null(names(modules))) {
    stop("'modules' must be a named list.")
  }

  modules <- modules[which(lengths(modules) < maxlen & lengths(modules) > minlen)]
  chunks <- split(modules, ceiling(seq_along(modules) / chunksize))

  ## first check which modules are informative (stddev > min_sd)
  if (downsample_explore < 1) {
    message(length(chunks), " chunks will be calculated.")
    obj2 <- downsample_seurat(obj, downsample = downsample_explore)
    message(length(Seurat::Cells(obj)), " cells downsampled to ", length(Seurat::Cells(obj2)), ".")
    modules2 <- purrr::map(chunks, function(x) {
      res <- UCell::AddModuleScore_UCell(
        obj2,
        features = x,
        BPPARAM = BiocParallel::MulticoreParam(workers=ncores, progressbar = T),
        chunk.size = 1000,
        name = "_xxUCellxx",
        ...)@meta.data |>
        dplyr::select(dplyr::ends_with("_xxUCellxx"))
      sds <- apply(res, 2, sd)
      return(names(res[,which(sds>min_sd)]))
    })
    modules <- modules[gsub("_xxUCellxx", "", purrr::list_c(modules2))]
  }

  ## then run relevant on full data
  chunks2 <- split(modules, ceiling(seq_along(modules) / chunksize))
  message("Calculating ", length(modules), " modules on full data.")
  mods <- purrr::map(chunks2, function(x) {
    res <- UCell::AddModuleScore_UCell(
      obj,
      features = x,
      BPPARAM = BiocParallel::MulticoreParam(workers=ncores, progressbar = T),
      chunk.size = 1000,
      name = "_UCell",
      ...)@meta.data |>
      dplyr::select(dplyr::ends_with("_UCell"))
    sds <- apply(res, 2, sd)
    return(res[,which(sds>min_sd)])
  })

  mods <- dplyr::bind_cols(mods)
  obj <- Seurat::AddMetaData(obj, mods)
  return(obj)
}
