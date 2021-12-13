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
             "plum1", "tan4", "khaki1", "cadetblue2", "olivedrab1", "orange", "purple4", "yellowgreen", "violetred3", "rosybrown3", "aquamarine3",
             "grey35", "lavender", "bisque", "mintcream", "yellow", "lightyellow2", "grey10"))
  }
}


'
 library(scales)

  tt <- c("grey65", "darkgoldenrod1", "cornflowerblue", "forestgreen", "tomato2", "mediumpurple3", "turquoise3", "lightgreen", "navy", "red4",
          "plum1", "tan4", "khaki1", "cadetblue2", "olivedrab1", "orange", "purple4", "yellowgreen", "violetred3", "rosybrown3", "aquamarine3",
          "grey35", "lavender", "bisque", "mintcream", "yellow", "lightyellow2", "grey10")
show_col(tt)

  install.packages("schemr")
  library(schemr)
  install.packages("ColorNameR")
  library(ColorNameR)

  r_color <- colors()
  r_color <- r_color[which(!grepl("^grey|^gray", r_color))]
  out <- as.data.frame(schemr::rgb_to_lab(t(col2rgb(r_color))))
  rownames(out) <- r_color

  euclidean <- function(x, mat) {
    sqrt(sum((as.numeric(mat[as.character(x[[1]]),]) - as.numeric(mat[x[2],]))^2))
  }
  table <- data.frame(dist = unlist(pbapply::pblapply(combn(r_color[1:100], 2, simplify = F), mat = out, euclidean)))

  table$col1 <- sapply(combn(r_color[1:100], 2, simplify = F), "[", 1)
  table$col2 <- sapply(combn(r_color[1:100], 2, simplify = F), "[", 2)
  library(tidyverse)
  library(viridis)

  ggplot(table %>% dplyr::filter(col1 != col2), aes(x=col1,y=col2,fill=dist)) +
    geom_tile() +
    scale_fill_viridis() +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_blank())


  table <- data.frame(dist = unlist(pbapply::pblapply(split(expand.grid(tt, r_color), 1:nrow(expand.grid(tt, r_color))), mat = out, euclidean)))
'