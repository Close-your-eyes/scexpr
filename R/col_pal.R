#' Color palettes
#'
#' Wrapper function around paletteer and two additional color palettes. See paletteer::palettes_c_names and paletteer::palettes_d_names.
#'
#' @param name name of the palette
#' @param n number of colors to return; may not work for every palette
#' @param nbrew number of color from brewer palettes
#' @param direction reverse palette with -1
#'
#' @return a color palette as character vector
#' @export
#'
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' col_pal(name = "custom", n = 10)
#' }
col_pal <- function(name = NULL,
                    n = NULL,
                    direction = c(1,-1)) {

  if (!requireNamespace("paletteer", quietly = T)) {
    utils::install.packages("paletteer")
  }

  paletteers <- dplyr::bind_rows(paletteer::palettes_c_names %>% dplyr::mutate(type2 = "continuous"),
                                 paletteer::palettes_d_names %>% dplyr::mutate(type2 = "discrete")) %>%
    dplyr::mutate(command = paste0(package, "::", palette))

  if (is.null(name)) {
    message("Select one palette by palette or command. Additional ones are 'custom' and 'ggplot' or 'hue'.")
    return(paletteers)
  }

  direction <- match.arg(direction, choices = c(1,-1))

  if (name %in% c("ggplot", "ggplot2", "hue", "hue_pal", "huepal")) {
    if (is.null(n)) {
      n <- 100
    }
    pal_select <- prismatic::color(scales::hue_pal()(n))
    if (direction == -1) {
      pal_select <- rev(pal_select)
    }
  } else if (name == "custom") {
    pal_select <- c("grey65", "darkgoldenrod1", "cornflowerblue", "forestgreen", "tomato2", "mediumpurple1", "turquoise3", "lightgreen", "navy", "plum1",
                    "red4", "khaki1", "tan4", "cadetblue1", "olivedrab3", "darkorange2", "burlywood2", "violetred3", "aquamarine3",
                    "grey30", "lavender", "yellow", "grey10", "pink3", "turquoise4", "darkkhaki", "magenta", "blue", "green", "blueviolet", "red",
                    "darkolivegreen", "orchid1", "springgreen", "dodgerblue4", "deepskyblue", "palevioletred4", "gold4", "maroon1", "lightyellow", "greenyellow", "purple4")
    if (!is.null(n)) {
      pal_select <- pal_select[1:n]
    }
    if (direction == -1) {
      pal_select <- rev(pal_select)
    }
  } else {
    if (grepl("::", name)) {
      pal_select <- paletteers %>% dplyr::filter(tolower(command) == tolower(name))
    } else {
      pal_select <- paletteers %>% dplyr::filter(tolower(palette) == tolower(name))
    }

    if (nrow(pal_select) == 0) {
      stop("Palette not found.")
    } else if (nrow(pal_select) > 1) {
      stop("Name is ambiguous. Please specify by command.")
    }

    if (pal_select$type2 == "discrete") {
      type <- "discrete"
      if (is.null(n)) {
        n <- pal_select$length
      }
      if (n > pal_select$length) {
        message("n = ", n, " larger than number of discrete color in palette (", pal_select$length, "). Going to interpolate to provide ", n, " colors.")
        type <- "continuous"
      }
      pal_return <- paletteer::paletteer_d(pal_select$command, n = n, type = type, direction = direction)
    } else if (pal_select$type2 == "continuous") {
      if (is.null(n)) {
        n <- 100
      }
      pal_return <- paletteer::paletteer_c(pal_select$command, n = n, direction = direction)
    }
  }
  return(pal_return)
}
