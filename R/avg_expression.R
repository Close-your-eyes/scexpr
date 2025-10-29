#' Average expression
#'
#' When layer is data and assay RNA or SCT, values are transformed to linear
#' space with expm1 before application of fun. This is the default behavior as
#' in Seurat::AverageExpression. So, averaging happens in linear space. To transform
#' average values back to log space, provide log1p to fun2.
#' For simple summation of counts as in Seurat::AggregateExpression, make
#' layer = "counts" and fun = Matrix::rowSums.
#' This function does not change levels of group as Seurats functions do.
#' What a pain!
#'
#' @param obj seurat object
#' @param group group col from meta.data, if NULL, no split (so overall average)
#' @param assay assay
#' @param features select feature subset
#' @param cells select cell subset
#' @param layer which layer
#' @param fun aggregation function
#' @param fun2 function applied after fun, e.g. log1p. base:identity is neutral
#' fun and does nothing
#' @param split column in meta.data to split by before aggregation
#'
#' @returns matrix with average expression values or list thereof when
#' !is.null(split)
#' @export
#'
#' @examples
avg_expression <- function(obj,
                           group = NULL,
                           assay = "RNA",
                           layer = "data",
                           features = NULL,
                           cells = NULL,
                           split = NULL,
                           fun = Matrix::rowMeans,
                           fun2 = base::identity) {

  if (!is.null(group)) {
    f <- obj@meta.data[[group]]
  } else {
    # return x as x below, no splitting
    f <- rep("1", nrow(obj@meta.data))
  }

  if (is.null(split)) {
    f <- list(f)
    obj <- list(obj)
  } else {
    f <- split(f, obj@meta.data[[split]])
    obj <- Seurat::SplitObject(obj, split.by = split)
  }


  transformer <- if (layer == "data" && assay %in% c("RNA", "SCT")) expm1 else identity

  out <- purrr::map2(obj, f, function(obj, f) {
    obj <- brathering::split_mat(x = get_layer(obj = obj,
                                               assay = assay,
                                               layer = layer,
                                               features = features,
                                               cells = cells),
                                 f = f,
                                 byrow = F)

    obj <- purrr::map(obj, ~fun2(fun(transformer(.x))))

    obj <- do.call(cbind, obj)
    return(obj)
  })

  if (length(out) == 1) {
    out <- out[[1]]
  }
  return(out)
}
