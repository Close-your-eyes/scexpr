#' Calculate a multi-classifier model with xgboost
#'
#' xgbTree, randomForest of generazed linear model (see examples)
#' try with https://www.tidymodels.org/ somewhen.
#'
#' @param df data frame with classes in first column (=character). otherwise
#' numeric columns only.
#' @param split train/test split fraction
#' @param preprocess preprocess columns of df, e.g. c("center", "scale"), see
#' ?caret::train
#' @param train_iter in case of cross-validation (method_train = cv) this is
#' n folds; passed to number of caret::trainControl
#' @param method_model which model from caret. function was made for xgbtree,
#' see all models: https://topepo.github.io/caret/available-models.html
#' @param method_train method arg of caret::trainControl
#' @param train_search search arg of caret::trainControl
#' @param tuneGrid tuneGrid arg of caret::train; set NULL and search=random to
#' have 30 combinations automatically picked; when search is random it is
#' the permitted range per parameter; when search is grid it gives all parameter
#' combinations (brute force all combinations with expand.grid or fine tune, e.g.
#' some pars depend on each other like eta and rounds)
#' @param tuneLength tuneLength arg of caret::train; when train search is random
#' this is the number of iterations per element of tuneGrid; tuneGrid then defines
#' the allowed parameter range; INTERACTION OF TUNEGRID AND TUNELENGTH NOT
#' CLEAR YET
#' @param seed random seed for reproducibility
#' @param ... more args to caret::train
#' @param train_repeat
#'
#' @returns list with train obj from caret and accessory info
#' @export
#'
#' @examples
#'\dontrun{
#' ## try parallel computing
#' num_cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(12)
#' doParallel::registerDoParallel(cl)
#'
#' parallel::stopCluster(cl)
#' foreach::registerDoSEQ()
#'
#' model <- run_xgboost_multi_classifier(tt,
#'                                       train_search = "random",
#'                                       method_model = "xgbTree",
#'                                       tuneGrid = NULL,
#'                                       tuneLength = 10)
#'
#' # random Forest
#' rf_model <- run_xgboost_multi_classifier(tt,
#'                                          train_search = "random",
#'                                          method_model = "rf",
#'                                          tuneGrid = NULL,
#'                                          tuneLength = 10,
#'                                          ntree = 200) # ntree rf specific
#'
#' rf_model <- run_xgboost_multi_classifier(tt,
#'                                          train_search = "grid",
#'                                          method_model = "rf",
#'                                          tuneGrid = data.frame(mtry = c(2, 4, 5)), # mtry bound by dimension in tt?!
#'                                          ntree = 200)
#'
#' # generalized linear models (with LASSO regularization)
#' # use preprocessing here by default: scaling
#' glm_model <- run_xgboost_multi_classifier(tt,
#'                                           train_search = "random",
#'                                           method_model = "glmnet",
#'                                           tuneGrid = NULL,
#'                                           tuneLength = 10,
#'                                           preprocess = c("center", "scale"))
#' # alpha = 1 → LASSO
#' # alpha = 0 → Ridge
#' # 0 < alpha < 1 → Elastic Net
#' # all regularized linear models
#' # Ordinary regression: minimize error only
#' # Regularized regression: minimize error + penalty on coefficients
#' # Ridge Regression (L2):
#' # Shrinks coefficients toward 0, Never sets them exactly to 0
#' # When to use
#' # Many correlated features, You want to keep all variables
#' # LASSO (L1):
#' # Can set some exactly to 0 → feature selection
#' # When to use
#' # You want a simpler model, Many irrelevant features
#' # Elastic Net (L1 + L2):
#' # mixture: Select features, but don’t be too aggressive
#'
#' # Correlated features:
#' # Ridge → spreads weight across them
#' # LASSO → picks one, drops others
#' # Elastic Net → keeps a group of them
#'
#' lasso_grid <- expand.grid(
#'   alpha = c(0,1),
#'   lambda = seq(0.001, 0.2, length = 30))
#' glm_model <- run_xgboost_multi_classifier(tt,
#'                                           train_search = "grid",
#'                                           method_model = "glmnet",
#'                                           tuneGrid = lasso_grid,
#'                                           preprocess = c("center", "scale"))
#'
#' train_iter:
#' # folds = data subsets
#' # k = 3 --> 3 subsetted groups f1, f2, f3
#' # then cross-validation:
#' # training on f1,f2 --> validation on f3
#' # training on f1,f3 --> validation on f2
#' # training on f2,f3 --> validation on f1
#' # Get a more reliable estimate of model performance
#' # Reduce overfitting
#' # tune grid:
#' # eta: learning rate, 0.01 → slow, safer, needs more trees, 0.3 → fast, risk of overfitting 0.1 → very common default
#' # nrounds: number of tree, depends on/interacts with eta: smaller eta needs more trees, 200 is too low for eta = 0.01
#' # max_depth: tree complexity, 3–5 → safer, less overfitting, 7–9 → more complex
#' # gamma, regularization, Minimum loss reduction to split. 0 → aggressive splitting, higher (1,5) → more conservative, useful default: 0.1
#' # colsample_bytree = Feature subsampling. 0.5 to 1 are useful
#' # subsample = row subsample. 0.5 to 1
#' # min_child_weight: Controls leaf size / regularization. 1 → flexible, larger → more conservative. from 1 to 10
#'
#' # a brute force grid with all combinations may not be meaningful due to dependencies of parameters, e.g. eta and rounds
#' # low eta → high nrounds
#' # deep trees → need more regularization
#' # small subsample → reduces overfitting
#'
#' # alternative: random parameter search with train_search = random, then tuneLength is used instead of tuneGrid
#' # if tuneGrid is not NULL then it defines the search space:
#' #' # for search = grid
#' # insane brute-force:
#' tuneGrid <- expand.grid(
#'   nrounds = c(100, 300, 500, 800, 1200, 2000),
#'   max_depth = c(2, 3, 5, 7, 9, 12),
#'   eta = c(0.005, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3),
#'   gamma = c(0, 0.1, 1, 5, 10),
#'   colsample_bytree = c(0.4, 0.6, 0.8, 1),
#'   min_child_weight = c(1, 3, 5, 10, 20),
#'   subsample = c(0.5, 0.7, 0.85, 1)
#' )
#' # sample from it
#' tuneGrid_small <- tuneGrid[sample(nrow(tuneGrid), 100), ]
#'
#' # less insane:
#' tuneGrid <- expand.grid(
#'   nrounds = c(100, 500, 1000),
#'   max_depth = c(3, 6, 9),
#'   eta = c(0.01, 0.05, 0.1, 0.3),
#'   gamma = c(0, 1, 5),
#'   colsample_bytree = c(0.6, 0.8, 1),
#'   min_child_weight = c(1, 5, 10),
#'   subsample = c(0.7, 1)
#' )
#' # even less:
#' tuneGrid = expand.grid(
#'   eta = c(0.01, 0.1),
#'   nrounds = c(100, 200),
#'   max_depth = c(6, 9),
#'   gamma = c(0, 1),
#'   colsample_bytree = c(0.6, 0.8),
#'   min_child_weight = c(1, 5),
#'   subsample = c(0.7, 1))
#' # define range for random search (unclear how tuneLength interacts)
#' tuneGrid = expand.grid(
#'   nrounds = seq(100, 1000, by = 100),
#'   max_depth = 3:10,
#'   eta = c(0.01, 0.05, 0.1, 0.3),
#'   gamma = c(0, 0.1, 1, 5),
#'   colsample_bytree = seq(0.5, 1, by = 0.1),
#'   min_child_weight = c(1, 3, 5, 10),
#'   subsample = seq(0.5, 1, by = 0.1))
#'
#' # run in parallel
#' # library(doParallel)
#' # cl <- makeCluster(4)
#' # registerDoParallel(cl)
#' # disable
#' # registerDoSEQ()
#'
#' ## method specific examples
#' # 1. GLMNET (LASSO / Elastic Net) wins
#' # Linear + sparse signal
#' # Only a few features matter
#' # Relationship is linear
#' # Many irrelevant features
#' n <- 800
#' p <- 40
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' colnames(X) <- paste0("x", 1:p)
#' # Sparse linear signal: only first 4 variables matter
#' eta <- 2 * X[,1] - 1.5 * X[,2] + 1.2 * X[,3] - 1 * X[,4]
#' prob <- 1 / (1 + exp(-eta))
#' y <- factor(ifelse(runif(n) < prob, "yes", "no"))
#' tt <- data.frame(y, X)
#' model <- run_xgboost_multi_classifier(tt,
#'                                       train_search = "random",
#'                                       method_model = "xgbTree",
#'                                       tuneGrid = NULL,
#'                                       tuneLength = 10,
#'                                       preprocess = c("center", "scale"))
#' rf_model <- run_xgboost_multi_classifier(tt,
#'                                          train_search = "random",
#'                                          method_model = "rf",
#'                                          tuneGrid = NULL,
#'                                          tuneLength = 10,
#'                                          preprocess = c("center", "scale"),
#'                                          ntree = 200) # ntree rf specific
#' glm_model <- run_xgboost_multi_classifier(tt,
#'                                           train_search = "random",
#'                                           method_model = "glmnet",
#'                                           tuneGrid = NULL,
#'                                           tuneLength = 10,
#'                                           preprocess = c("center", "scale"))
#' model$confmat_test$overall
#' model$confmat_test$table
#'
#' rf_model$confmat_test$overall
#' rf_model$confmat_test$table
#'
#' glm_model$confmat_test$overall
#' glm_model$confmat_test$table
#'
#'
#' # 2. Random Forest wins
#' # Nonlinear + interactions + noise
#' # Strong interactions
#' # Nonlinear relationships
#' # Moderate noise
#' n <- 800
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- runif(n)
#' x4 <- runif(n)
#' # Noise predictors
#' noise <- matrix(rnorm(n * 20), nrow = n, ncol = 20)
#' colnames(noise) <- paste0("z", 1:20)
#' # Nonlinear interaction / rule-based structure
#' y <- ifelse(
#'   (x1 > 0.7 & x2 < 0.4) | (x3 > 0.8) | (x1 * x4 > 0.5),
#'   "yes", "no"
#' )
#' # Add label noise
#' flip <- sample(1:n, size = round(0.1 * n))
#' y[flip] <- ifelse(y[flip] == "yes", "no", "yes")
#' y <- factor(y)
#' tt <- data.frame(y, x1, x2, x3, x4, noise)
#' model <- run_xgboost_multi_classifier(tt,
#'                                       train_search = "random",
#'                                       method_model = "xgbTree",
#'                                       tuneGrid = NULL,
#'                                       tuneLength = 10,
#'                                       preprocess = c("center", "scale"))
#' rf_model <- run_xgboost_multi_classifier(tt,
#'                                          train_search = "random",
#'                                          method_model = "rf",
#'                                          tuneGrid = NULL,
#'                                          tuneLength = 10,
#'                                          preprocess = c("center", "scale"),
#'                                          ntree = 200) # ntree rf specific
#' glm_model <- run_xgboost_multi_classifier(tt,
#'                                           train_search = "random",
#'                                           method_model = "glmnet",
#'                                           tuneGrid = NULL,
#'                                           tuneLength = 10,
#'                                           preprocess = c("center", "scale"))
#' model$confmat_test$overall
#' model$confmat_test$table
#'
#' rf_model$confmat_test$overall
#' rf_model$confmat_test$table
#'
#' glm_model$confmat_test$overall
#' glm_model$confmat_test$table
#'
#'
#' # 3. XGBoost wins
#' # Complex nonlinear structure + subtle patterns
#' # Smooth nonlinearities
#' # Additive + interaction effects
#' # Requires boosting precision
#'
#' n <- 1000
#' x1 <- runif(n, -3, 3)
#' x2 <- runif(n, -3, 3)
#' x3 <- runif(n, -3, 3)
#' x4 <- runif(n, -3, 3)
#' # Many small smooth effects (boosting-friendly)
#' eta <- 0.5 * sin(x1) +
#'   0.4 * cos(x2) +
#'   0.3 * x3 +
#'   0.2 * x4 +
#'   0.3 * sin(x1 * x2) +
#'   0.2 * x3 * x4
#' prob <- 1 / (1 + exp(-eta))
#' y <- factor(ifelse(runif(n) < prob, "yes", "no"), levels = c("yes", "no"))
#' tt <- data.frame(y, x1, x2, x3, x4)
#' model <- run_xgboost_multi_classifier(tt,
#'                                       train_search = "random",
#'                                       method_model = "xgbTree",
#'                                       tuneGrid = NULL,
#'                                       tuneLength = 10,
#'                                       preprocess = c("center", "scale"))
#' rf_model <- run_xgboost_multi_classifier(tt,
#'                                          train_search = "random",
#'                                          method_model = "rf",
#'                                          tuneGrid = NULL,
#'                                          tuneLength = 10,
#'                                          preprocess = c("center", "scale"),
#'                                          ntree = 200) # ntree rf specific
#' glm_model <- run_xgboost_multi_classifier(tt,
#'                                           train_search = "random",
#'                                           method_model = "glmnet",
#'                                           tuneGrid = NULL,
#'                                           tuneLength = 10,
#'                                           preprocess = c("center", "scale"))
#' model$confmat_test$overall
#' model$confmat_test$table
#'
#' rf_model$confmat_test$overall
#' rf_model$confmat_test$table
#'
#' glm_model$confmat_test$overall
#' glm_model$confmat_test$table
#'
#' final_model <- model[["trainobj"]][["finalModel"]]
#' lambda_opt <- model[["trainobj"]][["bestTune"]]$lambda
#' coefs <- coef(final_model, lambda_opt)
#' coefs <- purrr::map(coefs, ~as.data.frame(.x))
#' coefs <- purrr::map_dfc(names(coefs), ~setNames(coefs[[.x]], .x))
#'
#' }
run_xgboost_multi_classifier <- function(df,
                                         split = 0.7,
                                         preprocess = NULL,
                                         train_iter = 3,
                                         method_model = c("xgbTree", "rf", "glmnet"),
                                         method_train = c("cv", "repeatedcv", "boot", "boot632", "LOOCV", "none"),
                                         train_search = c("grid", "random"),
                                         train_repeat = 1,
                                         tuneGrid = expand.grid(
                                           eta = c(0.01, 0.1),
                                           nrounds = c(100, 200),
                                           max_depth = c(6, 9),
                                           gamma = c(0, 1),
                                           colsample_bytree = c(0.6, 0.8),
                                           min_child_weight = c(1, 5),
                                           subsample = c(0.7, 1)),
                                         tuneLength = 3,
                                         seed = 42,
                                         ...) {



  # tune grid:
  # eta: learning rate, 0.01 → slow, safer, needs more trees, 0.3 → fast, risk of overfitting 0.1 → very common default
  # nrounds: number of tree, depends on/interacts with eta: smaller eta needs more trees, 200 is too low for eta = 0.01
  # max_depth: tree complexity, 3–5 → safer, less overfitting, 7–9 → more complex
  # gamma, regularization, Minimum loss reduction to split. 0 → aggressive splitting, higher (1,5) → more conservative, useful default: 0.1
  # colsample_bytree = Feature subsampling. 0.5 to 1 are useful
  # subsample = row subsample. 0.5 to 1
  # min_child_weight: Controls leaf size / regularization. 1 → flexible, larger → more conservative. from 1 to 10

  # a brute force grid with all combinations may not be meaningful due to dependencies of parameters, e.g. eta and rounds
  # low eta → high nrounds
  # deep trees → need more regularization
  # small subsample → reduces overfitting

  # alternative: random parameter search with train_search = random, then tuneLength is used instead of tuneGrid
  # if tuneGrid is not NULL then it defines the search space:
  # grid <- expand.grid(
  #   nrounds = seq(100, 1000, by = 100),
  #   max_depth = 3:10,
  #   eta = c(0.01, 0.05, 0.1, 0.3),
  #   gamma = c(0, 0.1, 1, 5),
  #   colsample_bytree = seq(0.5, 1, by = 0.1),
  #   min_child_weight = c(1, 3, 5, 10),
  #   subsample = seq(0.5, 1, by = 0.1)
  # )
  # tuneLength: number of different values to try for each tuning parameter

  # manual grid:
  # tuneGrid = expand.grid(
  #   eta = c(0.01, 0.1, 0.3),
  #   nrounds = c(100, 200),
  #   max_depth = c(3, 6, 9),
  #   gamma = c(0, 1),
  #   colsample_bytree = c(0.6, 0.8),
  #   min_child_weight = c(1, 5),
  #   subsample = c(0.7, 1))

  # run in parallel:
  # library(doParallel)
  # cl <- makeCluster(4)
  # registerDoParallel(cl)
  # disable:
  # registerDoSEQ()

  if (compareVersion(as.character(packageVersion("xgboost")), "1.8") == 1) {
    message("https://stackoverflow.com/questions/79849114/new-version-of-xgboost-package-is-not-working-under-caret-environment")
    stop("install old version of xgboost 1.7.11.1 like so: install.packages('xgboost', repos = 'https://p3m.dev/cran/2025-12-01')")
  }

  if (!is.data.frame(df)) {
    stop("df must be data frame.")
  }

  if (!all(apply(df[,-1], 2, is.numeric))) {
    stop("all columns of df except first must be numeric.")
  }

  if (anyNA(df[,1])) {
    stop("NA in first column.")
  }

  if (length(unique((df[,1]))) == 1) {
    stop("Only one level in first column.")
  }

  if (any(make.names(df[,1]) != df[,1])) {
    message("first column of df must valid R names as by make.names. applying make.names.")
    df[,1] <- make.names(df[,1])
  }

  method_model <- rlang::arg_match(method_model)
  method_train <- rlang::arg_match(method_train)

  train_search <- rlang::arg_match(train_search)
  if (train_search == "grid") {
    message("tuneGrid: ", nrow(tuneGrid), " combinations.")
    tuneLength <- NULL
  }

  set.seed(seed)
  train_index <- caret::createDataPartition(y = df[,1,drop=T],
                                            p = split,
                                            list = FALSE)[,1]
  train_data <- df[train_index, ]
  test_data  <- df[-train_index, ]

  # folds = data subsets
  # k = 3 --> 3 subsetted groups f1, f2, f3
  # then cross-validation:
  # training on f1,f2 --> validation on f3
  # training on f1,f3 --> validation on f2
  # training on f2,f3 --> validation on f1
  # Get a more reliable estimate of model performance
  # Reduce overfitting

  trControl <- caret::trainControl(method = method_train,
                                   number = train_iter,
                                   search = train_search,
                                   repeats = train_repeat,
                                   verboseIter = TRUE,
                                   classProbs = TRUE)

  # if (method_model == "rf") {
  #   trControl <- caret::trainControl(method = method_train,
  #                                    number = train_iter)
  # }

  # seed: for very precise control use seeds argument in trainControl
  # this controls every single split during training
  # however, required seeds depends on method and resamples
  # needed seeds = (number of resamples) + 1
  # with method_train = cv and 5 folds: 5 + 1
  # maybe just set seed once before caret::tain for now

  #requireNamespaceQuietStop <- caret:::requireNamespaceQuietStop
  set.seed(seed)
  trainobj <- caret::train(
    as.formula(paste(names(df)[1], "~ .")),
    data = train_data,
    preProcess = preprocess,
    method = method_model,
    trControl = trControl,
    tuneLength = tuneLength,
    tuneGrid = tuneGrid,
    metric = "Accuracy",
    ...)


  # Predict classes
  # on test data
  pred_label_test <- stats::predict(trainobj, newdata = test_data)
  confmat_test <- caret::confusionMatrix(pred_label_test, as.factor(test_data[,1]))

  # on full data
  pred_label_full <- stats::predict(trainobj, newdata = df[,-1])
  confmat_full <- caret::confusionMatrix(pred_label_full, as.factor(df[,1]))
  label = df[,1,drop=T]


  feat_cor_mat <- stats::cor(as.matrix(df[,-1]))
  # caret::findCorrelation(feat_cor_mat)
  cor_feat <- brathering::mat_to_df_long(feat_cor_mat) |>
    dplyr::filter(cname != rname) |>
    dplyr::filter(value >= 0.9)


  feat_imp <- NULL
  if (methods::is(trainobj$finalModel, "xgb.Booster")) {
    feat_imp <- as.data.frame(xgboost::xgb.importance(model = trainobj$finalModel))
    message(nrow(feat_imp), " of ", length(trainobj[["coefnames"]]), " features used in model. Checkout feat_imp.")
  }



  return(list(
    confmat_test = confmat_test,
    trainobj = trainobj,
    feat_cor_mat = feat_cor_mat,
    cor_feat = cor_feat,
    feat_imp = feat_imp,
    confmat_full = confmat_full,
    pred_label_full = pred_label_full,
    label = label,
    mislabel = which(pred_label_full != label)
  ))
}



run_xgboost_multi_classifier_tidymodels <- function(
    df,
    split = 0.7,
    preprocess = NULL,
    train_iter = 3,
    train_method = c("cv", "boot", "none"),
    train_search = c("grid", "random"),
    tune_grid = NULL,
    tune_length = 10,
    seed = 42
) {

  base::set.seed(seed)


  # checks

  if (!base::is.data.frame(df)) {
    base::stop("df must be data frame.")
  }

  if (!base::all(base::apply(df[, -1], 2, base::is.numeric))) {
    base::stop("all columns except first must be numeric.")
  }

  df[[1]] <- base::as.factor(df[[1]])
  outcome_name <- base::names(df)[1]


  # split

  split_obj <- rsample::initial_split(
    df,
    prop = split,
    strata = rlang::sym(outcome_name)
  )

  train_data <- rsample::training(split_obj)
  test_data  <- rsample::testing(split_obj)


  # recipe

  rec <- recipes::recipe(
    stats::as.formula(base::paste(outcome_name, "~ .")),
    data = train_data
  )

  if (!base::is.null(preprocess) && base::is.function(preprocess)) {
    rec <- preprocess(rec)
  }


  # model

  model <- parsnip::boost_tree(
    trees = tune::tune(),
    tree_depth = tune::tune(),
    learn_rate = tune::tune(),
    loss_reduction = tune::tune(),
    sample_size = tune::tune(),
    mtry = tune::tune(),
    min_n = tune::tune()
  ) |>
    parsnip::set_engine("xgboost") |>
    parsnip::set_mode("classification")


  # workflow

  wf <- workflows::workflow() |>
    workflows::add_model(model) |>
    workflows::add_recipe(rec)


  # resampling

  train_method <- rlang::arg_match(train_method)

  resamples <- base::switch(
    train_method,
    cv = rsample::vfold_cv(
      train_data,
      v = train_iter,
      strata = rlang::sym(outcome_name)
    ),
    boot = rsample::bootstraps(train_data, times = train_iter),
    none = NULL
  )


  # tuning grid

  if (base::is.null(tune_grid)) {

    if (train_search == "grid") {
      grid <- dials::grid_regular(
        dials::parameters(wf),
        levels = tune_length
      )
    } else {
      grid <- dials::grid_random(
        dials::parameters(wf),
        size = tune_length
      )
    }

  } else {
    grid <- tune_grid
  }


  # fit / tune

  if (!base::is.null(resamples)) {

    tuned <- tune::tune_grid(
      wf,
      resamples = resamples,
      grid = grid,
      metrics = yardstick::metric_set(yardstick::accuracy),
      control = tune::control_grid(verbose = TRUE)
    )

    best_params <- tune::select_best(tuned, metric = "accuracy")

    final_wf <- workflows::finalize_workflow(wf, best_params)

    fit_obj <- parsnip::fit(final_wf, data = train_data)

  } else {

    fit_obj <- parsnip::fit(wf, data = train_data)
  }


  # predictions

  pred_test <- parsnip::predict(fit_obj, test_data) |>
    dplyr::bind_cols(test_data)

  pred_full <- parsnip::predict(fit_obj, df) |>
    dplyr::bind_cols(df)


  # confusion matrices

  confmat_test <- yardstick::conf_mat(
    pred_test,
    truth = rlang::sym(outcome_name),
    estimate = .pred_class
  )

  confmat_full <- yardstick::conf_mat(
    pred_full,
    truth = rlang::sym(outcome_name),
    estimate = .pred_class
  )


  # feature correlation

  feat_cor_mat <- stats::cor(base::as.matrix(df[, -1]))

  cor_feat <- base::as.data.frame(base::as.table(feat_cor_mat)) |>
    dplyr::filter(.data$Var1 != .data$Var2, .data$Freq >= 0.9)


  # feature importance

  feat_imp <- NULL

  fitted_engine <- workflows::extract_fit_engine(fit_obj)

  if (base::inherits(fitted_engine, "xgb.Booster")) {
    feat_imp <- xgboost::xgb.importance(model = fitted_engine)
  }


  # outputs

  return(base::list(
    confmat_test = confmat_test,
    model = fit_obj,
    confmat_full = confmat_full,
    pred_full = pred_full$.pred_class,
    label = df[[1]],
    mislabel = base::which(pred_full$.pred_class != df[[1]]),
    feat_cor_mat = feat_cor_mat,
    cor_feat = cor_feat,
    feat_imp = feat_imp
  ))
}

