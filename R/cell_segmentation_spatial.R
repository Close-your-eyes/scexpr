#' Compute Cell Hull Polygons from Spatial Mask
#'
#' Extracts convex hull polygons for each cell (grouped by `z`) from a spatial
#' mask. The mask is optionally filtered to a region defined by cell coordinates
#' or a custom hull.
#'
#' @param obj A spatial object containing image and coordinate data (e.g.,
#'   Seurat-like object with `@misc[["spatial"]]` and `@reductions$spatial`).
#' @param mask A mask array representing segmented regions. Defaults to
#'   `obj@misc[["spatial"]][["staining_image_mask"]]`.
#' @param cell_coords A matrix or data frame of cell coordinates (typically
#'   `obj@reductions$spatial@cell.embeddings`). Used for filtering and defining
#'   the region of interest.
#' @param hull_filter Optional data frame of coordinates (`x`, `y`) defining a
#'   polygon to restrict the mask. If `NULL`, a convex hull is computed from
#'   `cell_coords`.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Subsets the mask using \code{get_mask_subset()}
#'   \item Groups pixels by cell ID (`z`)
#'   \item Computes a convex hull per group using \code{grDevices::chull()}
#' }
#'
#' The resulting output can be used for plotting segmentation boundaries.
#'
#' @return A data frame with columns `x`, `y`, and `z`, representing polygon
#'   vertices for each cell hull.
#'
#' @export
#'
#' @examples
get_cell_hulls <- function(obj,
                           mask = obj@misc[["spatial"]][["staining_image_mask"]],
                           cell_coords = obj@reductions$spatial@cell.embeddings,
                           hull_filter = NULL) {

  mask <- get_mask_subset(obj = obj,
                          mask = mask,
                          cell_coords = cell_coords,
                          hull_filter = hull_filter)
  # brathering::plot2(mask, size = 1)
  # brathering::plot2(cell_coords, size = 1)
  # brathering::plot2(hulls, size = 1)

  hulls <- dplyr::slice(mask, grDevices::chull(x,y), .by = z)

  return(hulls)
}



#' Subset Spatial Mask to Region of Interest
#'
#' Converts a raw image or mask array into a data frame and optionally restricts
#' it to a region defined by cell coordinates or a custom polygon.
#'
#' @param obj seurat object
#' @param mask A raw image or mask array to convert and subset.
#' @param cell_coords Optional matrix or data frame of cell coordinates. If
#'   provided and `hull_filter` is `NULL`, a convex hull is computed from these
#'   coordinates to define the region of interest.
#' @param hull_filter Optional data frame with columns `x`, `y` defining a polygon
#'   used to filter the mask. Overrides `cell_coords` if provided.
#' @param subset_factor Numeric. Downsampling factor for the mask to improve
#'   performance. Default is `1` (no downsampling).
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Optionally computes a convex hull from `cell_coords`
#'   \item Converts the mask array into a data frame using
#'     \code{raw_img_array_to_df()}
#'   \item Filters pixels inside the hull (if provided)
#'   \item Optionally filters by cell IDs (`z`) matching `cell_coords`
#' }
#'
#' If `cell_coords` has more than two columns, only the first two are used.
#'
#' @return A data frame containing at least `x`, `y`, and optionally `z` and
#'   `color` columns, depending on the input mask.
#'
#' @export
#'
#' @examples
get_mask_subset <- function(obj,
                            mask,
                            cell_coords = NULL,
                            hull_filter = NULL,
                            subset_factor = 1) {

  # optional cell_coords to increase speed of raw_img_array_to_df
  if (!is.null(cell_coords) && is.null(hull_filter)) {
    if (ncol(cell_coords) > 2) {
      message("Using cols 1,2 of cell_coords.")
    }
    cell_coords <- cell_coords[,c(1,2)]
    cell_coords <- as.data.frame(cell_coords)
    if (!all(colnames(cell_coords) != c("x", "y"))) {
      message("cols 1,2 of cell_coords are interpreted as x,y.")
    }
    colnames(cell_coords) <- c("x","y")
    hull_filter <- dplyr::slice(cell_coords, grDevices::chull(x, y))
    # brathering::plot2(hull_filter)
  }

  # mask matrix to df

  mask <- raw_img_array_to_df(raw_array = mask,
                              subset_factor = subset_factor,
                              hull_inside = T,
                              hull_filter = hull_filter)
  # all(rownames(cell_coords) %in% Seurat::Cells(obj))
  # all(Seurat::Cells(obj) %in% rownames(cell_coords))
  # setdiff(mask$z, rownames(cell_coords))
  # brathering::plot2(mask2[,c(1,2)])

  if (!is.null(cell_coords) && !is.null(rownames(cell_coords)) && "z" %in% names(mask)) {
    # check for mistmatch levels?
    mask <- dplyr::filter(mask, z %in% rownames(cell_coords))
  }

  return(mask)
}


#' Extract and Prepare Spatial Image Data for Plotting
#'
#' Converts a spatial image array into a data frame suitable for plotting,
#' optionally restricting it to a region of interest and applying a color
#' transformation.
#'
#' @param obj seurat object
#' @param img_array A raw image array. Defaults to
#'   `obj@misc[["spatial"]][["staining_image"]]`.
#' @param cell_coords Optional matrix or data frame of cell coordinates used to
#'   define a region of interest.
#' @param hull_filter Optional polygon (`x`, `y`) used to restrict the image.
#' @param subset_factor Numeric. Downsampling factor for performance. Default is `1`.
#' @param col_conv_fun Optional function to transform color values (e.g.,
#'   grayscale conversion).
#'
#' @details
#' This function:
#' \itemize{
#'   \item Calls \code{get_mask_subset()} to subset and reshape the image
#'   \item Optionally transforms the `color` column using `col_conv_fun`
#' }
#'
#' The output is typically used with \code{ggplot2::geom_raster()}.
#'
#' @return A data frame with columns `x`, `y`, and `color`.
#'
#' @export
#'
#' @examples
get_img <- function(obj,
                    img_array = obj@misc[["spatial"]][["staining_image"]],
                    cell_coords = obj@reductions$spatial@cell.embeddings,
                    hull_filter = NULL,
                    subset_factor = 1,
                    col_conv_fun = NULL) {

  img <- get_mask_subset(obj = obj,
                         mask = img_array,
                         cell_coords = cell_coords,
                         hull_filter = hull_filter,
                         subset_factor = subset_factor)

  if (!is.null(col_conv_fun)) {
    img$color <- col_conv_fun(img$color)
  }

  #hull <- dplyr::slice(img, grDevices::chull(x,y))

  return(img)
}
