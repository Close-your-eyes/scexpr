#' PF-log1p-PF Depth Normalization
#'
#' Normalizes a gene-by-cell count matrix using proportional fitting,
#' a \code{log1p} transformation, and within-cell centering. Use it instead of
#' ScaleData to generate scale.data slot with HVF only (--> dense matrix). Then
#' PF is only used for dim red.
#' - https://x.com/lpachter/status/2064795979568783813
#' - https://x.com/lpachter/status/2064796022572986638/photo/1
#'
#' @param x A nonnegative numeric or sparse matrix with genes in rows and
#'   cells in columns. Every cell must have a positive total count.
#'
#' @return A dense numeric matrix with the same dimensions as \code{x}.
#'   Columns correspond to cells and have means approximately equal to zero.
#'
#' @details
#' The procedure consists of three steps:
#'
#' \enumerate{
#'   \item Scale each cell to the mean cell depth.
#'   \item Apply \code{log1p} to the scaled counts.
#'   \item Subtract the mean transformed value of each cell.
#' }
#'
#' The final centering step converts the result to a dense matrix, which
#' may require substantial memory for large single-cell datasets.
#'
#' @references
#' Booeshaghi AS, Hallgrímsdóttir IB, Gálvez-Merchán Á, Pachter L.
#' Depth normalization for single-cell genomics count data.
#' \doi{10.1101/2022.05.06.490859}
#'
#' @export
#'
#' @examples
#' counts <- Matrix::Matrix(
#'   matrix(
#'     c(
#'       1, 0, 4,
#'       0, 2, 1,
#'       3, 0, 0,
#'       0, 1, 2
#'     ),
#'     nrow = 4,
#'     ncol = 3,
#'     byrow = TRUE
#'   ),
#'   sparse = TRUE
#' )
#'
#' normalized <- proportional_fit_pachterlab(counts)
#' normalized
#' colMeans(normalized)
proportional_fit_pachterlab_matrix <- function(x) {
  # Scale cells to mean depth
  log1ppf <- pf(mtx = x)

  # Apply log1p only to stored nonzero entries
  log1ppf@x <- log1p(log1ppf@x)

  # Center each cell across genes
  cell_mean <- Matrix::colMeans(log1ppf)

  pflog1ppf <- sweep(
    as.matrix(log1ppf),
    MARGIN = 2,
    STATS = cell_mean,
    FUN = "-"
  )

  return(pflog1ppf)
}

#' PF-log1p-PF Normalization with Covariate Regression
#'
#' Applies proportional-fitting normalization to selected features of a Seurat
#' assay, optionally regresses cell-level covariates from every feature, and
#' writes the result to an assay layer such as \code{scale.data}.
#'
#' @param object A Seurat object.
#' @param assay A single character string specifying the assay to process.
#'   Defaults to \code{"RNA"}.
#' @param layer A single character string specifying the input assay layer.
#'   Defaults to \code{"data"}.
#' @param output_layer A single character string specifying the layer in which
#'   to store the result. Defaults to \code{"scale.data"}.
#' @param features An optional character vector containing the features to
#'   process. When \code{NULL}, the variable features defined for
#'   \code{assay} are used. Feature order is preserved in the output.
#' @param vars.to.regress An optional character vector containing names of
#'   cell-level variables in \code{object[[]]} to regress from every selected
#'   feature. For example, \code{c("nCount_RNA", "percent.mt")}. Variables are
#'   included jointly in a linear model. Numeric variables are used directly;
#'   factors and character variables are expanded into indicator variables.
#'   When \code{NULL}, no regression is performed.
#' @param do.center Logical indicating whether to center each feature across
#'   cells after optional regression. Defaults to \code{TRUE}.
#' @param do.scale Logical indicating whether to scale each feature across
#'   cells. Centered features are divided by their sample standard deviations.
#'   If \code{do.center = FALSE}, features are divided by their root mean
#'   squares. Defaults to \code{TRUE}.
#' @param scale.max An optional positive numeric value specifying the maximum
#'   absolute scaled value. Scaled values are restricted to
#'   \code{c(-scale.max, scale.max)}. Set to \code{NULL} to disable clipping.
#'   This argument is ignored when \code{do.scale = FALSE}. Defaults to
#'   \code{10}.
#'
#' @return The modified Seurat object. The specified \code{output_layer}
#'   contains a dense numeric matrix with selected features in rows and cells
#'   in columns.
#'
#' @details
#' If \code{features = NULL}, variable features must already have been defined,
#' usually by running \code{Seurat::FindVariableFeatures()}. Supplying
#' \code{features} explicitly does not require variable features to be defined.
#'
#' The procedure consists of the following steps:
#'
#' \enumerate{
#'   \item Extract the requested feature-by-cell matrix from \code{layer}.
#'   \item Calculate each cell's depth using only the selected features.
#'   \item Scale each cell to the mean cell depth.
#'   \item Apply \code{log1p} to the scaled expression values.
#'   \item Center each cell across the selected features.
#'   \item Optionally regress the variables in \code{vars.to.regress} from
#'     every feature.
#'   \item Optionally center, scale, and clip each feature across cells.
#' }
#'
#' For feature \eqn{g} and cell \eqn{i}, covariate regression fits
#'
#' \deqn{
#' y_{gi} = \beta_{0g} +
#' \sum_{k=1}^{p} \beta_{kg} x_{ik} +
#' \epsilon_{gi},
#' }
#'
#' where \eqn{y_{gi}} is the transformed expression value and \eqn{x_{ik}} is
#' covariate \eqn{k} for cell \eqn{i}. The fitted covariate component is
#' removed, and the residual is retained:
#'
#' \deqn{
#' r_{gi} = y_{gi} -
#' \left(
#' \widehat{\beta}_{0g} +
#' \sum_{k=1}^{p}\widehat{\beta}_{kg}x_{ik}
#' \right).
#' }
#'
#' Regression is performed simultaneously for all requested variables using a
#' QR-based linear-model fit.
#'
#' When both \code{do.center} and \code{do.scale} are \code{TRUE}, residuals
#' are converted to feature-wise standardized values:
#'
#' \deqn{
#' z_{gi} = \frac{r_{gi} - \bar{r}_g}{s_g}.
#' }
#'
#' The output matrix is dense and may require substantial memory for large
#' datasets.
#'
#' @section Important considerations:
#' Variables in \code{vars.to.regress} should represent unwanted sources of
#' variation. If a regressed variable is correlated with the biological
#' condition of interest, regression can remove genuine biological signal.
#'
#' The original normalization procedure is intended for nonnegative count
#' data. Seurat's \code{"data"} layer is commonly already normalized and
#' log-transformed. To apply the procedure to raw counts, use
#' \code{layer = "counts"}.
#'
#' Cell depth is calculated using only the selected features. Consequently,
#' using variable features can produce different scaling factors than using
#' all genes.
#'
#' @references
#' Booeshaghi AS, Hallgrímsdóttir IB, Gálvez-Merchán Á, Pachter L.
#' Depth normalization for single-cell genomics count data.
#' \doi{10.1101/2022.05.06.490859}
#'
#' @seealso
#' \code{Seurat::FindVariableFeatures()},
#' \code{Seurat::ScaleData()},
#' \code{SeuratObject::VariableFeatures()},
#' \code{SeuratObject::LayerData()}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use variable features
#' object <- Seurat::FindVariableFeatures(
#'   object,
#'   assay = "RNA"
#' )
#'
#' object <- proportional_fit_pachterlab(object)
#'
#' # Use explicitly selected features
#' object <- proportional_fit_pachterlab(
#'   object,
#'   features = c("CD3D", "CD3E", "IL7R", "CCR7")
#' )
#'
#' # Use raw counts and regress technical covariates
#' object <- proportional_fit_pachterlab(
#'   object,
#'   assay = "RNA",
#'   layer = "counts",
#'   vars.to.regress = c("nCount_RNA", "percent.mt")
#' )
#'
#' # Retrieve the result
#' normalized <- SeuratObject::LayerData(
#'   object,
#'   assay = "RNA",
#'   layer = "scale.data"
#' )
#' }
proportional_fit_pachterlab <- function(
    object,
    assay = "RNA",
    layer = "data",
    output_layer = "scale.data",
    features = NULL,
    vars.to.regress = NULL,
    do.center = TRUE,
    do.scale = TRUE,
    scale.max = 10
) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }

  if (
    length(assay) != 1L ||
    !is.character(assay) ||
    is.na(assay) ||
    !nzchar(assay)
  ) {
    stop("`assay` must be a single nonempty string.", call. = FALSE)
  }

  if (!assay %in% names(object)) {
    stop(
      sprintf("Assay '%s' was not found in the Seurat object.", assay),
      call. = FALSE
    )
  }

  if (
    length(layer) != 1L ||
    !is.character(layer) ||
    is.na(layer) ||
    !nzchar(layer)
  ) {
    stop("`layer` must be a single nonempty string.", call. = FALSE)
  }

  if (
    length(output_layer) != 1L ||
    !is.character(output_layer) ||
    is.na(output_layer) ||
    !nzchar(output_layer)
  ) {
    stop(
      "`output_layer` must be a single nonempty string.",
      call. = FALSE
    )
  }

  if (
    length(do.center) != 1L ||
    !is.logical(do.center) ||
    is.na(do.center)
  ) {
    stop("`do.center` must be TRUE or FALSE.", call. = FALSE)
  }

  if (
    length(do.scale) != 1L ||
    !is.logical(do.scale) ||
    is.na(do.scale)
  ) {
    stop("`do.scale` must be TRUE or FALSE.", call. = FALSE)
  }

  if (
    !is.null(scale.max) &&
    (
      !is.numeric(scale.max) ||
      length(scale.max) != 1L ||
      !is.finite(scale.max) ||
      scale.max <= 0
    )
  ) {
    stop(
      "`scale.max` must be a positive finite number or NULL.",
      call. = FALSE
    )
  }

  available_layers <- SeuratObject::Layers(
    object = object,
    assay = assay
  )

  if (!layer %in% available_layers) {
    stop(
      sprintf(
        "Layer '%s' was not found in assay '%s'. Available layers: %s",
        layer,
        assay,
        paste(available_layers, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Use variable features by default
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(object[[assay]])

    if (!length(features)) {
      stop(
        sprintf(
          paste0(
            "No variable features are defined for assay '%s'. ",
            "Run FindVariableFeatures() or supply `features`."
          ),
          assay
        ),
        call. = FALSE
      )
    }
  } else {
    if (
      !is.character(features) ||
      !length(features) ||
      anyNA(features) ||
      any(!nzchar(features))
    ) {
      stop(
        "`features` must be a nonempty character vector without missing values.",
        call. = FALSE
      )
    }

    features <- unique(features)
  }

  # Extract requested features; Seurat matrices are genes × cells
  x <- SeuratObject::LayerData(
    object = object,
    assay = assay,
    layer = layer,
    features = features
  )

  missing_features <- setdiff(features, rownames(x))

  if (length(missing_features)) {
    stop(
      sprintf(
        "Features absent from layer '%s': %s",
        layer,
        paste(missing_features, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Preserve the requested feature ordering
  x <- x[features, , drop = FALSE]

  if (anyNA(x) || any(!is.finite(x))) {
    stop(
      "The input expression matrix contains missing or nonfinite values.",
      call. = FALSE
    )
  }

  if (any(x < 0)) {
    stop(
      "The input expression matrix must contain nonnegative values.",
      call. = FALSE
    )
  }

  # Scale cells to mean depth
  transformed <- pf(x)

  # Apply log1p; sparse matrices remain sparse
  transformed <- log1p(transformed)

  # Center every cell across the selected features
  transformed <- sweep(
    as.matrix(transformed),
    MARGIN = 2,
    STATS = Matrix::colMeans(transformed),
    FUN = "-"
  )

  # Optionally regress cell-level metadata from every feature
  if (!is.null(vars.to.regress)) {
    if (
      !is.character(vars.to.regress) ||
      !length(vars.to.regress) ||
      anyNA(vars.to.regress) ||
      any(!nzchar(vars.to.regress))
    ) {
      stop(
        paste0(
          "`vars.to.regress` must be NULL or a nonempty character ",
          "vector without missing values."
        ),
        call. = FALSE
      )
    }

    vars.to.regress <- unique(vars.to.regress)

    metadata <- object[[]]

    missing_variables <- setdiff(
      vars.to.regress,
      colnames(metadata)
    )

    if (length(missing_variables)) {
      stop(
        sprintf(
          "Variables absent from object metadata: %s",
          paste(missing_variables, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    missing_cells <- setdiff(
      colnames(transformed),
      rownames(metadata)
    )

    if (length(missing_cells)) {
      stop(
        sprintf(
          "%d expression-matrix cells were absent from object metadata.",
          length(missing_cells)
        ),
        call. = FALSE
      )
    }

    latent_data <- metadata[
      colnames(transformed),
      vars.to.regress,
      drop = FALSE
    ]

    if (anyNA(latent_data)) {
      stop(
        "`vars.to.regress` contains missing values.",
        call. = FALSE
      )
    }

    design <- tryCatch(
      stats::model.matrix(
        ~ .,
        data = latent_data
      ),
      error = function(error) {
        stop(
          sprintf(
            "Could not construct the regression model: %s",
            conditionMessage(error)
          ),
          call. = FALSE
        )
      }
    )

    design_qr <- qr(design)

    if (design_qr$rank >= nrow(design)) {
      stop(
        paste0(
          "The regression model has no residual degrees of freedom. ",
          "Reduce the number of covariates or factor levels."
        ),
        call. = FALSE
      )
    }

    # Each column of the response is one feature
    fit <- stats::lm.fit(
      x = design,
      y = t(transformed)
    )

    transformed <- t(fit$residuals)

    dimnames(transformed) <- list(
      features,
      colnames(x)
    )
  }

  # Center each feature across cells
  if (do.center) {
    transformed <- sweep(
      transformed,
      MARGIN = 1,
      STATS = rowMeans(transformed),
      FUN = "-"
    )
  }

  # Scale each feature across cells
  if (do.scale) {
    if (do.center) {
      denominator <- max(ncol(transformed) - 1L, 1L)

      feature_scale <- sqrt(
        rowSums(transformed^2) / denominator
      )
    } else {
      feature_scale <- sqrt(
        rowMeans(transformed^2)
      )
    }

    # Constant features remain zero after centering
    feature_scale[
      !is.finite(feature_scale) |
        feature_scale == 0
    ] <- 1

    transformed <- sweep(
      transformed,
      MARGIN = 1,
      STATS = feature_scale,
      FUN = "/"
    )

    if (!is.null(scale.max)) {
      transformed <- pmax(
        -scale.max,
        pmin(scale.max, transformed)
      )
    }
  }

  # Ensure the output retains feature and cell names
  dimnames(transformed) <- list(
    features,
    colnames(x)
  )

  # Write the dense result to the requested assay layer
  SeuratObject::LayerData(
    object = object,
    assay = assay,
    layer = output_layer
  ) <- transformed

  object
}


#' Scale Cells to the Mean Cell Depth
#'
#' Scales every column of a gene-by-cell matrix so that its total equals the
#' mean cell depth of the input matrix.
#'
#' @param mtx A nonnegative numeric or sparse matrix with features in rows and
#'   cells in columns.
#'
#' @return A matrix with the same dimensions and dimnames as \code{mtx}.
#'
#' @noRd
pf <- function(mtx) {
  depth <- Matrix::colSums(mtx)

  if (any(!is.finite(depth))) {
    stop("Cell depths must be finite.", call. = FALSE)
  }

  if (any(depth <= 0)) {
    stop(
      "Some cells have zero total expression among the selected features.",
      call. = FALSE
    )
  }

  mtx %*% Matrix::Diagonal(
    x = mean(depth) / depth
  )
}

