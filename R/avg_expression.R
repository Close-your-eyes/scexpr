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
#' @param return_as how or what to return
#' @param na_to_zero impute NA values to zero?
#'
#' @returns matrix with average expression values or list thereof when
#' !is.null(split)
#' @export
#'
#' @examples
#' \dontrun{
#' avgexpr <- avg_expression(so, "celltype")
#' }
avg_expression <- function(obj,
                           group = NULL,
                           assay = "RNA",
                           layer = "data",
                           features = NULL,
                           cells = NULL,
                           split = NULL,
                           fun = Matrix::rowMeans,
                           fun2 = base::identity,
                           return_as = c("list", "df"),
                           na_to_zero = F) {

  return_as <- rlang::arg_match(return_as)

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
    f <- f[names(obj)]
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
    #obj <- purrr::map(obj, ~transformer(.x))
    #tt <- as.matrix(obj[["PDC"]])
    obj <- do.call(cbind, obj)
    return(obj)
  })

  if (na_to_zero) {
    ## impute NA with zeros
    out <- purrr::map(out, function(x) {
      x[which(is.na(x))] <- 0
      return(x)
    })
  }

  if (return_as == "df") {
    out <- purrr::map(
      .x = out,
      .f = brathering::mat_to_df_long,
      rownames_to = "feature",
      colnames_to = "group",
      values_to = "expr"
    )
    out <- dplyr::bind_rows(out, .id = split)
  }

  return(out)
}
