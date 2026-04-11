#' Train a Clustering Classifier Model from a Seurat Object
#'
#' Builds a multi-class classification model using features extracted from a
#' Seurat object. The function can use dimensional reductions, variable features,
#' or a user-specified feature set as input.
#'
#' @param obj A \code{Seurat} object containing expression data and metadata.
#' @param meta_col A character string specifying the metadata column used as the
#'   target variable (e.g., cluster labels). Default is \code{"orig.ident"}.
#' @param features A character vector specifying features to use. Can be:
#'   \itemize{
#'     \item The name of a dimensional reduction (e.g., \code{"pca"})
#'     \item \code{"var_feat"} to automatically use variable features
#'     \item A vector of gene/feature names
#'   }
#' @param model_args A named list of arguments passed to the model training
#'   function (e.g., \code{run_xgboost_multi_classifier}). Includes options for
#'   model type, cross-validation, and tuning.
#' @param get_layer_args A named list of arguments passed to \code{get_layer()},
#'   controlling which assay and layer are used to extract expression data.
#'
#' @returns A trained classification model object (typically from
#'   \code{run_xgboost_multi_classifier}).
#'
#' @details
#' The function constructs a feature matrix from either:
#' \itemize{
#'   \item A dimensional reduction stored in \code{obj@reductions}
#'   \item Variable features obtained via \code{scexpr::get_var_features()}
#'   \item A user-defined set of features extracted using \code{get_layer()}
#' }
#'
#' The resulting data frame is combined with the specified metadata column and
#' passed to a multi-class classifier training function.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## try parallel computing
#' num_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(12)
#' doParallel::registerDoParallel(cl)
#'
#' parallel::stopCluster(cl)
#' foreach::registerDoSEQ()
#' }
make_clustering_classifier_model <- function(obj,
                                             meta_col = "orig.ident",
                                             features = "var_feat",
                                             model_args = list(method_model = "glmnet",
                                                               train_search = "random",
                                                               method_train = "repeatedcv",
                                                               train_iter = 5,
                                                               train_repeat = 3,
                                                               tuneGrid = NULL,
                                                               tuneLength = 10),
                                             get_layer_args = list(assay = "RNA",
                                                                   layer = "data")) {

  meta_col <- scexpr:::check.features(obj, features = meta_col)

  if (length(features) == 1 && features %in% names(obj@reductions)) {
    df <- as.data.frame(obj@reductions[[features]]@cell.embeddings)
    message("using reduction: ", features)
  } else {
    if (length(features) == 1 && features == "var_feat") {
      features <- scexpr::get_var_features(obj)
    }
    df <- Gmisc::fastDoCall(what = get_layer,
                            args = c(list(obj = obj,
                                          features = features,
                                          as = "df",
                                          transpose = T),
                                     get_layer_args))
    nme1 <- names(df)
    names(df) <- gsub("-", "_", names(df))
    nme1 <- stats::setNames(names(df), nme1)
  }
  df <- cbind(obj@meta.data[,meta_col,drop = F], df)


  model <- Gmisc::fastDoCall(what = run_xgboost_multi_classifier,
                             args = c(list(df = df),
                                      model_args))
  return(c(model, list(feature_name_conv = nme1)))
}

