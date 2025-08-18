#' Perform Regularized Regression Analysis
#'
#' @title Perform LASSO and Elastic Net Regression with Cross-Validation
#'
#' @description
#' This function performs regularized regression analysis using LASSO (Least Absolute
#' Shrinkage and Selection Operator) and/or Elastic Net regression methods. Both methods
#' are regularization techniques that prevent overfitting by adding penalty terms to the
#' loss function. LASSO uses L1 regularization for feature selection by shrinking
#' coefficients to zero, while Elastic Net combines L1 and L2 penalties to handle
#' multicollinearity and perform simultaneous feature selection. The function supports
#' both binary and multinomial classification tasks with comprehensive model evaluation
#' and result reporting.
#'
#' @param data A list object containing preprocessed data. Must be the output from the
#'   \code{perform_PreprocessingPeakData} function, containing the following elements:
#'   \itemize{
#'     \item{\code{FunctionOrigin}}: Character string indicating data source
#'     \item{\code{Metadata}}: Data frame with sample metadata including Group column
#'     \item{\code{data_scaledPCA_rsdFiltered_varFiltered}}: Matrix of preprocessed features
#'   }
#' @param method Character vector specifying regression method(s) to perform. Options:
#'   \itemize{
#'     \item{\code{"lasso"}}: LASSO regression only (L1 penalty, alpha = 1)
#'     \item{\code{"enet"}}: Elastic Net regression only (L1+L2 penalty, alpha = 0.5)
#'     \item{\code{c("lasso", "enet")}}: Both methods (recommended for comparison)
#'   }
#'   Default: \code{"enet"}
#' @param specify_response Character string specifying the response variable column name.
#'   If \code{NULL}, uses the Group column from metadata. Default: \code{NULL}
#' @param train_percent Numeric value between 1 and 99 specifying the percentage of data
#'   to use for training. Remaining data used for testing. Default: \code{80}
#' @param ref Character string specifying the reference level for the response variable.
#'   If \code{NULL}, uses the first factor level alphabetically. Default: \code{NULL}
#' @param lambda Character string specifying lambda selection criterion:
#'   \itemize{
#'     \item{\code{"1se"}}: Lambda within one standard error of minimum (conservative, fewer features)
#'     \item{\code{"min"}}: Lambda that minimizes cross-validation error (aggressive, more features)
#'   }
#'   Default: \code{"1se"}
#' @param remember Numeric value for reproducible results. Sets random seed using
#'   \code{set.seed(remember)}. If \code{NULL}, no seed is set. Default: \code{NULL}
#' @param verbose Logical indicating whether to print progress messages and results
#'   to console. Default: \code{TRUE}
#' @param cv_folds Integer specifying number of cross-validation folds for model
#'   selection. Must be between 3 and 20. Default: \code{10}
#' @param parallel Logical indicating whether to use parallel processing for
#'   cross-validation. Default: \code{FALSE}
#'
#' @return A list containing regression results with the following structure:
#' \describe{
#'   \item{\code{FunctionOrigin}}{Character string identifying the source function}
#'   \item{\code{ModelSummary}}{Data frame summarizing model performance metrics}
#'   \item{\code{DataSplit}}{List containing training/testing data split information}
#'   \item{\code{LASSO_Results}}{List of LASSO results (if method includes "lasso")}
#'   \item{\code{ElasticNet_Results}}{List of Elastic Net results (if method includes "enet")}
#'   \item{\code{ComparisonSummary}}{Data frame comparing methods (if both performed)}
#' }
#'
#' Each method-specific results list contains:
#' \itemize{
#'   \item{\code{Model}}: Fitted cv.glmnet object
#'   \item{\code{Predictions}}: Data frame with actual vs predicted values
#'   \item{\code{ConfusionMatrix}}: Complete confusion matrix object
#'   \item{\code{Performance}}: Data frame with accuracy, sensitivity, specificity, etc.
#'   \item{\code{Coefficients}}: Data frame with non-zero coefficients and odds ratios
#'   \item{\code{Lambda}}: Selected lambda value
#'   \item{\code{Alpha}}: Alpha parameter used
#'   \item{\code{ReferenceLevel}}: Reference level for classification
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(glmnet)
#' library(caret)
#' library(dplyr)
#'
#' # Perform both LASSO and Elastic Net regression
#' regression_results <- perform_Regression(
#'   data = preprocessed_data,
#'   method = c("lasso", "enet"),
#'   train_percent = 75,
#'   lambda = "1se",
#'   remember = 123,
#'   cv_folds = 10
#' )
#'
#' # View model comparison
#' print(regression_results$ModelSummary)
#' print(regression_results$ComparisonSummary)
#'
#' # Access LASSO results
#' lasso_coef <- regression_results$LASSO_Results$Coefficients
#' lasso_perf <- regression_results$LASSO_Results$Performance
#'
#' # Access Elastic Net results
#' enet_coef <- regression_results$ElasticNet_Results$Coefficients
#' enet_perf <- regression_results$ElasticNet_Results$Performance
#'
#' # View confusion matrices
#' print(regression_results$LASSO_Results$ConfusionMatrix$table)
#' print(regression_results$ElasticNet_Results$ConfusionMatrix$table)
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link[caret]{confusionMatrix}}
#'
#' @importFrom glmnet cv.glmnet coef.glmnet
#' @importFrom caret confusionMatrix
#' @importFrom dplyr n_distinct
#' @importFrom stats predict
#'
#' @export
perform_Regression <- function(
    data,
    method = "enet",
    specify_response = NULL,
    train_percent = 80,
    ref = NULL,
    lambda = "1se",
    remember = NULL,
    verbose = TRUE,
    cv_folds = 10,
    parallel = FALSE
) {

  # Validate inputs
  .validate_inputs(data, method, train_percent, lambda, remember, cv_folds)

  # Set seed for reproducibility
  if (!is.null(remember)) {
    set.seed(remember)
    if (verbose) cat("Random seed set to:", remember, "\n")
  }

  # Prepare data
  data_prep <- .prepare_data_regression(data, specify_response, ref, train_percent, verbose)

  # Initialize results structure
  results <- list(
    FunctionOrigin = "perform_Regression",
    ModelSummary = data.frame(),
    DataSplit = data_prep$split_info,
    ComparisonSummary = NULL
  )

  # Store performance metrics for comparison
  performance_comparison <- list()

  # Perform regression for each method
  for (method_type in method) {

    if (verbose) {
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat("PERFORMING", toupper(method_type), "REGRESSION\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
    }

    # Fit model and get results
    method_results <- .fit_regression_model(
      x_train = data_prep$x_train,
      y_train = data_prep$y_train,
      x_test = data_prep$x_test,
      y_test = data_prep$y_test,
      method_type = method_type,
      lambda = lambda,
      groups = data_prep$groups,
      ref = data_prep$ref_level,
      cv_folds = cv_folds,
      parallel = parallel,
      verbose = verbose
    )

    # Store results
    if (method_type == "lasso") {
      results$LASSO_Results <- method_results
    } else {
      results$ElasticNet_Results <- method_results
    }

    # Store for comparison
    performance_comparison[[method_type]] <- method_results$Performance
  }

  # Create model summary
  results$ModelSummary <- .create_model_summary(performance_comparison, method)

  # Create comparison summary if both methods were used
  if (length(method) > 1) {
    results$ComparisonSummary <- .create_comparison_summary(performance_comparison)
  }

  if (verbose) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("REGRESSION ANALYSIS COMPLETE\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    print(results$ModelSummary)
    if (!is.null(results$ComparisonSummary)) {
      cat("\nMETHOD COMPARISON:\n")
      print(results$ComparisonSummary)
    }
  }

  return(results)
}

# Helper function: Validate inputs
.validate_inputs <- function(data, method, train_percent, lambda, remember, cv_folds) {

  # Check data structure
  if (!is.list(data)) {
    stop("'data' must be a list object from perform_PreprocessingPeakData function.")
  }

  if (!"FunctionOrigin" %in% names(data) || data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("The supplied data did not come from 'perform_PreprocessingPeakData' function. Required metadata might be missing.")
  }

  required_elements <- c("Metadata", "data_scaledPCA_rsdFiltered_varFiltered")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    stop(paste("Missing required data elements:", paste(missing_elements, collapse = ", ")))
  }

  # Check method
  allowed_methods <- c("lasso", "enet")
  if (!all(method %in% allowed_methods)) {
    stop("'method' must be one or more of: ", paste(allowed_methods, collapse = ", "))
  }

  # Check train_percent
  if (!is.numeric(train_percent) || train_percent <= 0 || train_percent >= 100) {
    stop("'train_percent' must be a numeric value between 1 and 99.")
  }

  # Check lambda
  if (!lambda %in% c("min", "1se")) {
    stop("'lambda' must be either 'min' or '1se'.")
  }

  # Check remember
  if (!is.null(remember) && !is.numeric(remember)) {
    stop("'remember' parameter must be numeric or NULL.")
  }

  # Check cv_folds
  if (!is.numeric(cv_folds) || cv_folds < 3 || cv_folds > 20) {
    stop("'cv_folds' must be an integer between 3 and 20.")
  }
}

# Helper function: Prepare data for modeling - FIXED VERSION
.prepare_data_regression <- function(data, specify_response, ref, train_percent, verbose) {

  # Extract metadata and features
  metadata <- data$Metadata
  features <- data$data_scaledPCA_rsdFiltered_varFiltered

  # Convert features to matrix if not already
  if (!is.matrix(features)) {
    features <- as.matrix(features)
  }

  # Identify QC samples
  qc_indices <- metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices <- !qc_indices

  if (sum(non_qc_indices) == 0) {
    stop("No non-QC samples found in the data.")
  }

  # Prepare response variable
  if (!is.null(specify_response)) {
    if (!specify_response %in% colnames(metadata)) {
      stop("Specified response variable '", specify_response, "' not found in metadata.")
    }
    groups <- metadata[[specify_response]][non_qc_indices]
  } else {
    groups <- metadata$Group[non_qc_indices]
  }

  # Remove NA values
  na_indices <- is.na(groups)
  if (any(na_indices)) {
    warning(paste("Removing", sum(na_indices), "samples with NA response values."))
    non_qc_indices[which(non_qc_indices)[na_indices]] <- FALSE
    groups <- groups[!na_indices]
  }

  # Set reference level
  ref_level <- NULL
  if (!is.null(ref)) {
    if (!ref %in% unique(groups)) {
      stop("Reference level '", ref, "' not found in response variable.")
    }
    ref_level <- ref
    groups <- factor(groups, levels = c(ref, setdiff(unique(groups), ref)))
  } else {
    groups <- factor(groups)
    ref_level <- levels(groups)[1]
  }

  # Check for sufficient samples per group
  group_counts <- table(groups)
  if (any(group_counts < 3)) {
    warning("Some groups have fewer than 3 samples. This may affect model performance.")
  }

  # Extract non-QC feature data
  df <- features[non_qc_indices, , drop = FALSE]
  df <- df[!na_indices, , drop = FALSE]

  # Ensure df is a matrix
  if (!is.matrix(df)) {
    df <- as.matrix(df)
  }

  # Check for constant or near-constant features
  feature_var <- apply(df, 2, var, na.rm = TRUE)
  constant_features <- is.na(feature_var) | feature_var < 1e-10
  if (any(constant_features)) {
    warning(paste("Removing", sum(constant_features), "constant or near-constant features."))
    df <- df[, !constant_features, drop = FALSE]
  }

  # Check for missing values in features
  if (any(is.na(df))) {
    warning("Missing values detected in features. Consider imputation.")
    # Option 1: Remove samples with missing values
    complete_cases <- complete.cases(df)
    if (sum(complete_cases) < nrow(df)) {
      warning(paste("Removing", sum(!complete_cases), "samples with missing feature values."))
      df <- df[complete_cases, , drop = FALSE]
      groups <- groups[complete_cases]
    }
  }

  # Split data with stratification to ensure balanced groups
  n_samples <- nrow(df)
  train_size <- floor((train_percent / 100) * n_samples)

  if (train_size < length(levels(groups))) {
    stop("Training set too small. Reduce train_percent or increase sample size.")
  }

  # Stratified sampling to ensure all groups are represented
  train_indices <- c()
  for (level in levels(groups)) {
    level_indices <- which(groups == level)
    level_train_size <- max(1, floor(length(level_indices) * train_percent / 100))
    level_train_indices <- sample(level_indices, level_train_size)
    train_indices <- c(train_indices, level_train_indices)
  }

  # If we need more samples to reach train_size, add randomly from remaining
  remaining_indices <- setdiff(seq_len(n_samples), train_indices)
  if (length(train_indices) < train_size && length(remaining_indices) > 0) {
    additional_needed <- min(train_size - length(train_indices), length(remaining_indices))
    additional_indices <- sample(remaining_indices, additional_needed)
    train_indices <- c(train_indices, additional_indices)
  }

  x_train <- df[train_indices, , drop = FALSE]
  x_test <- df[-train_indices, , drop = FALSE]
  y_train <- groups[train_indices]
  y_test <- groups[-train_indices]

  # Final validation
  if (nrow(x_train) == 0 || nrow(x_test) == 0) {
    stop("Data split resulted in empty training or testing set.")
  }

  if (ncol(x_train) == 0) {
    stop("No features remaining after preprocessing.")
  }

  # Ensure factor levels are consistent
  y_train <- factor(y_train, levels = levels(groups))
  y_test <- factor(y_test, levels = levels(groups))

  if (verbose) {
    cat("Data preparation complete:\n")
    cat("  Total samples:", n_samples, "\n")
    cat("  Training samples:", nrow(x_train), "\n")
    cat("  Testing samples:", nrow(x_test), "\n")
    cat("  Features:", ncol(x_train), "\n")
    cat("  Groups:", paste(levels(groups), collapse = ", "), "\n")
    cat("  Reference level:", ref_level, "\n")
    cat("  Training group distribution:\n")
    print(table(y_train))
    cat("  Testing group distribution:\n")
    print(table(y_test))
  }

  return(list(
    x_train = x_train,
    x_test = x_test,
    y_train = y_train,
    y_test = y_test,
    groups = groups,
    ref_level = ref_level,
    split_info = list(
      total_samples = n_samples,
      train_samples = nrow(x_train),
      test_samples = nrow(x_test),
      features = ncol(x_train),
      groups = levels(groups),
      reference_level = ref_level,
      train_indices = train_indices
    )
  ))
}

# Helper function: Fit regression model - FIXED VERSION
.fit_regression_model <- function(x_train, y_train, x_test, y_test, method_type,
                                  lambda, groups, ref, cv_folds, parallel, verbose) {

  # Set alpha value
  alpha_val <- if (method_type == "lasso") 1 else 0.5

  # Determine family
  n_groups <- length(levels(groups))
  family_type <- if (n_groups == 2) "binomial" else "multinomial"

  if (verbose) {
    cat("  Alpha:", alpha_val, "\n")
    cat("  Family:", family_type, "\n")
    cat("  CV folds:", cv_folds, "\n")
    cat("  Training data dimensions:", paste(dim(x_train), collapse = " x "), "\n")
    cat("  Testing data dimensions:", paste(dim(x_test), collapse = " x "), "\n")
  }

  # Additional data validation before modeling
  if (any(is.na(x_train)) || any(is.na(y_train))) {
    stop("Training data contains missing values. Please clean data before modeling.")
  }

  if (any(is.na(x_test)) || any(is.na(y_test))) {
    stop("Testing data contains missing values. Please clean data before modeling.")
  }

  # Check for adequate sample size per fold
  min_group_size <- min(table(y_train))
  if (cv_folds > min_group_size) {
    cv_folds <- max(3, min_group_size - 1)
    warning(paste("Reducing cv_folds to", cv_folds, "due to small group sizes."))
  }

  # Fit model with enhanced error handling
  tryCatch({
    model_fit <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      alpha = alpha_val,
      family = family_type,
      type.measure = "class",
      nfolds = cv_folds,
      parallel = parallel,
      standardize = TRUE,  # Ensure standardization
      maxit = 100000       # Increase max iterations if needed
    )
  }, error = function(e) {
    stop("Model fitting failed: ", e$message,
         "\nDebugging info:",
         "\n  x_train dimensions: ", paste(dim(x_train), collapse = " x "),
         "\n  y_train length: ", length(y_train),
         "\n  y_train levels: ", paste(levels(y_train), collapse = ", "),
         "\n  Alpha: ", alpha_val,
         "\n  Family: ", family_type,
         "\n  CV folds: ", cv_folds)
  })

  # Select lambda
  selected_lambda <- if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min

  if (verbose) {
    cat("  Selected lambda:", selected_lambda, "\n")
  }

  # Make predictions with error handling
  predictions <- tryCatch({
    pred_result <- stats::predict(
      model_fit,
      s = selected_lambda,
      newx = x_test,
      type = "class"
    )
    as.vector(pred_result)
  }, error = function(e) {
    stop("Prediction failed: ", e$message)
  })

  # Ensure predictions have correct factor levels
  predictions <- factor(predictions, levels = levels(groups))

  # Create confusion matrix
  conf_matrix <- caret::confusionMatrix(
    predictions,
    factor(y_test, levels = levels(groups))
  )

  # Extract coefficients and calculate odds ratios
  coef_results <- .extract_coefficients(model_fit, selected_lambda, groups, family_type)

  # Create performance summary
  performance <- .create_performance_summary(conf_matrix, method_type, alpha_val, selected_lambda, ref)

  if (verbose) {
    cat("  Accuracy:", round(performance$Accuracy, 4), "\n")
    cat("  Non-zero coefficients:", nrow(coef_results), "\n")
  }

  return(list(
    Model = model_fit,
    Predictions = data.frame(
      Actual = as.character(y_test),
      Predicted = as.character(predictions),
      stringsAsFactors = FALSE
    ),
    ConfusionMatrix = conf_matrix,
    Performance = performance,
    Coefficients = coef_results,
    Lambda = selected_lambda,
    Alpha = alpha_val,
    ReferenceLevel = ref
  ))
}

# Helper function: Extract coefficients and odds ratios
.extract_coefficients <- function(model_fit, lambda, groups, family_type) {

  coef_model <- glmnet::coef.glmnet(model_fit, s = lambda)

  if (family_type == "binomial") {
    # Binary classification
    coef_df <- as.data.frame(as.matrix(coef_model))
    colnames(coef_df) <- "Coefficient"
    coef_df$Feature <- rownames(coef_df)
    coef_df <- coef_df[coef_df$Coefficient != 0, ]
    coef_df$OddsRatio <- exp(coef_df$Coefficient)
    coef_df <- coef_df[, c("Feature", "Coefficient", "OddsRatio")]

  } else {
    # Multinomial classification
    group_levels <- levels(groups)

    coef_list <- lapply(seq_along(coef_model), function(i) {
      coef_df <- as.data.frame(as.matrix(coef_model[[i]]))
      colnames(coef_df) <- "Coefficient"
      coef_df$Feature <- rownames(coef_df)
      coef_df <- coef_df[coef_df$Coefficient != 0, ]
      coef_df$OddsRatio <- exp(coef_df$Coefficient)
      coef_df$Group <- group_levels[i]
      coef_df[, c("Group", "Feature", "Coefficient", "OddsRatio")]
    })

    coef_df <- do.call(rbind, coef_list)
  }

  rownames(coef_df) <- NULL
  return(coef_df)
}

# Helper function: Create performance summary
.create_performance_summary <- function(conf_matrix, method, alpha, lambda, ref) {

  overall <- conf_matrix$overall

  performance <- data.frame(
    Method = method,
    Alpha = alpha,
    Lambda = lambda,
    Reference = ref,
    Accuracy = overall["Accuracy"],
    AccuracyLower = overall["AccuracyLower"],
    AccuracyUpper = overall["AccuracyUpper"],
    Kappa = overall["Kappa"],
    stringsAsFactors = FALSE
  )

  # Add class-specific metrics if available
  if (!is.null(conf_matrix$byClass)) {
    by_class <- conf_matrix$byClass
    if (is.matrix(by_class)) {
      # Multi-class case
      performance$MeanSensitivity <- mean(by_class[, "Sensitivity"], na.rm = TRUE)
      performance$MeanSpecificity <- mean(by_class[, "Specificity"], na.rm = TRUE)
      performance$MeanPrecision <- mean(by_class[, "Pos Pred Value"], na.rm = TRUE)
      performance$MeanF1 <- mean(by_class[, "F1"], na.rm = TRUE)
    } else {
      # Binary case
      performance$Sensitivity <- by_class["Sensitivity"]
      performance$Specificity <- by_class["Specificity"]
      performance$Precision <- by_class["Pos Pred Value"]
      performance$F1 <- by_class["F1"]
    }
  }

  rownames(performance) <- NULL
  return(performance)
}

# Helper function: Create model summary
.create_model_summary <- function(performance_list, methods) {

  summary_df <- do.call(rbind, performance_list)
  rownames(summary_df) <- NULL

  # Add ranking based on accuracy
  summary_df$AccuracyRank <- rank(-summary_df$Accuracy, ties.method = "min")

  return(summary_df)
}

# Helper function: Create comparison summary
.create_comparison_summary <- function(performance_list) {

  if (length(performance_list) < 2) {
    return(NULL)
  }

  methods <- names(performance_list)
  metrics <- c("Accuracy", "Kappa")

  comparison <- data.frame(
    Metric = metrics,
    stringsAsFactors = FALSE
  )

  for (method in methods) {
    perf <- performance_list[[method]]
    comparison[[paste0(method, "_Value")]] <- perf[metrics]
  }

  # Calculate differences
  if (length(methods) == 2) {
    method1 <- methods[1]
    method2 <- methods[2]
    comparison[[paste0("Difference_", method1, "_vs_", method2)]] <-
      comparison[[paste0(method1, "_Value")]] - comparison[[paste0(method2, "_Value")]]
  }

  return(comparison)
}
