#' Convert Image Array to Data Frame for Plotting
#'
#' Converts a 2D or 3D image array into a data frame with `x`/`y` pixel
#' coordinates, suitable for visualization with \code{ggplot2}. Supports RGB
#' images (3D arrays) as well as single-channel masks (2D arrays).
#'
#' @param raw_array A numeric or integer array representing an image.
#'   Must be either:
#'   \itemize{
#'     \item A 3D array with RGB channels in the third dimension
#'     \item A 2D array representing a single-channel mask (e.g., segmentation IDs)
#'   }
#' @param subset_factor Numeric. Downsampling factor for rows and columns to
#'   reduce memory usage and improve performance. Default is `1` (no downsampling).
#' @param hull_filter Optional data frame with columns `x` and `y` defining a
#'   polygon. Used to subset the output to points inside (or outside) the polygon.
#' @param hull_inside Logical. If `TRUE` (default), keeps only points inside or on
#'   the polygon defined by `hull_filter`. If `FALSE`, keeps only points outside.
#' @param xrange Optional numeric vector of length 2 specifying the x-range to
#'   retain (min, max).
#' @param yrange Optional numeric vector of length 2 specifying the y-range to
#'   retain (min, max).
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Converts the input array to integer values
#'   \item Subsamples pixels using `subset_factor`
#'   \item Expands pixel indices into `x`, `y` coordinate pairs
#'   \item Extracts RGB channels (for 3D arrays) or assigns values to `z` (for 2D arrays)
#'   \item Optionally filters points using:
#'     \itemize{
#'       \item Bounding box filters (`xrange`, `yrange`)
#'       \item Polygon inclusion via \code{sp::point.in.polygon()}
#'     }
#'   \item For RGB images, combines channels into a hex color string using
#'     \code{grDevices::rgb()}
#' }
#'
#' Filtering is applied after subsampling to improve efficiency.
#'
#' @return A data frame with:
#' \itemize{
#'   \item `x`, `y`: pixel coordinates
#'   \item `color`: hex color values (for RGB input), OR
#'   \item `z`: pixel values (for 2D input, returned as character)
#' }
#'
#' @export
#'
#' @examples
raw_img_array_to_df <- function(raw_array,
                                subset_factor = 1,
                                hull_filter = NULL,
                                hull_inside = T,
                                xrange = NULL,
                                yrange = NULL) {

  # raw array either 3d with rgb in third dim or 2d

  # set filters like hull or range to save memory and speed up

  # hull def
  # hull <- dplyr::slice(data, grDevices::chull(x, y))

  # ranges as simple filters, less precise as hull

  ints <- as.integer(raw_array)
  dim(ints) <- dim(raw_array)

  #ints <- aperm(ints, c(2,1,3)) # swap dims
  #subset with loss of original ranges
  #ints <- ints[seq(500,4700,3),seq(2100,10700,3),]
  #ints <- ints[seq(1,dim(ints)[1], 4), seq(1,dim(ints)[2], 4), ]

  # Convert to raster
  # r <- as.raster(ints / 255)
  # plot(r)
  #
  # range(so@reductions[["spatial"]]@cell.embeddings[,1])
  # range(so@reductions$spatial@cell.embeddings[,2])

  x_idx <- seq(1, dim(ints)[1], by = subset_factor)
  y_idx <- seq(1, dim(ints)[2], by = subset_factor)

  if (length(dim(ints)) == 3) {
    df <- data.frame(
      expand.grid(x = x_idx, y = y_idx),
      r = as.vector(ints[x_idx, y_idx, 1]),
      g = as.vector(ints[x_idx, y_idx, 2]),
      b = as.vector(ints[x_idx, y_idx, 3]))
  } else if (length(dim(ints)) == 2) {
    df <- data.frame(
      expand.grid(x = x_idx, y = y_idx),
      z = as.character(as.vector(ints[x_idx, y_idx])))
  } else {
    stop("raw_array must be 3d with 3 layer in third dim (=rgb), or 2d.")
  }


  if (!is.null(xrange)) {
    df <- dplyr::filter(df, dplyr::between(x, min(xrange), max(xrange)))
  }
  if (!is.null(yrange)) {
    df <- dplyr::filter(df, dplyr::between(y, min(yrange), max(yrange)))
  }

  if (!is.null(hull_filter)) {
    inside <- sp::point.in.polygon(
      point.x = df$x,
      point.y = df$y,
      pol.x = hull_filter$x,
      pol.y = hull_filter$y)

    if (hull_inside) {
      df  <- df[inside > 0, ] # inside and edge
    } else {
      df <- df[inside == 0, ]  # outside
    }
  }
  if (length(dim(ints)) == 3) {
    df$color <- grDevices::rgb(df$r, df$g, df$b, maxColorValue=255)
    df <- df[,-which(names(df) %in% c("r", "g", "b"))]
  }

  # geom_polygon(data = hull, aes(x, y), fill = NA, color = "red", linewidth = 1) +
  return(df)
}
