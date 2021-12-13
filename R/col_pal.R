#' Color palettes
#'
#' @param name name of the palette
#'
#' @return a color palette
#' @export
#'
#' @examples
#' \dontrun{
#' }
col_pal <- function(name) {


  # https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
  # https://stackoverflow.com/questions/5392061/algorithm-to-check-similarity-of-colors

  name <- match.arg(name, c("custom"))
  if (name == "custom") {
    return(c("grey65", "darkgoldenrod1", "cornflowerblue", "forestgreen", "tomato2", "mediumpurple3", "turquoise3", "lightgreen", "navy", "red4",
             "plum1", "tan4", "khaki2", "cadetblue2", "olivedrab1", "orange", "purple4", "yellowgreen", "violetred3", "rosybrown3", "aquamarine3",
             "grey25", "maroon4"))
  }
}
