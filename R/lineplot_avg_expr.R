#' Make a lineplot of average gene expression across groups
#'
#' @param obj seurat object
#' @param features gene features
#' @param assay which assay
#' @param group groups on x-axis
#' @param split another group as color scale
#'
#' @returns ggplot
#' @export
#'
#' @examples
#' \dontrun{
#' gg <- lineplot_avg_expr(obj = so, features = c("CD3", "NKG7"), split = "cluster")
#' }
lineplot_avg_expr <- function(obj,
                              features,
                              assay = "RNA",
                              group = "orig.ident",
                              split = NULL,
                              return_df = F) {

  obj <- scexpr:::check.SO(SO = obj,
                           assay = assay,
                           #length = 1,
                           meta.col = c(group, split))

  if (missing(features)) {
    stop("features missing.")
  }

  features <- scexpr:::check.features(SO = obj,
                                      features = features)

  avgexpr <- purrr::map_dfr(obj, ~scexpr::avg_expression(
    obj = .x,
    assay = assay,
    features = features,
    group = group,
    split = split,
    return_as = "df"))


  if (is.null(split)) {
    line_mapping <- NULL
    point_mapping <- NULL
  } else {
    line_mapping <- ggplot2::aes(group = !!rlang::sym(split), color = !!rlang::sym(split))
    point_mapping <- ggplot2::aes(color = !!rlang::sym(split))
  }


  fct_lvl <- unlist(purrr::map(obj, function(x) {
    if (is.factor(x@meta.data[[split]])) {
      levels(x@meta.data[[split]])
    } else {
      unique(x@meta.data[[split]])
    }
  }))
  if (anyDuplicated(fct_lvl)) {
    message(("same split factor level in different obj."))
  }
  fct_lvl <- unique(fct_lvl)
  avgexpr[[split]] <- factor(avgexpr[[split]], levels = fct_lvl)

  p <- ggplot2::ggplot(avgexpr, ggplot2::aes(x = group, y = expr)) +
    ggplot2::geom_line(mapping = line_mapping) +
    ggplot2::geom_point(mapping = point_mapping, size = 2) +
    #ggplot2::scale_color_manual(values = so@misc$cluster_color) +
    #ggplot2::scale_fill_manual(values = so@misc$cluster_color) +
    #colrr::theme_material() +
    #colrr::guides_default() +
    ggplot2::labs(y = "avg expr [UMI]") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   strip.text.x = ggplot2::element_text(margin = ggplot2::margin(2,0,2,0, unit = "pt"))) +
    ggplot2::facet_wrap(ggplot2::vars(feature), scales = "free_y", axes = "all",
                        axis.labels = "margins")

  return(list(data = avgexpr, plot = p))
}
