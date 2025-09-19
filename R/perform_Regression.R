#' Perform Regularized Regression Analysis with Enhanced Features
#'
#' @title Perform LASSO and Elastic Net Regression with Cross-Validation for Categorical and Numeric Responses
#'
#' @description
#' This function performs regularized regression analysis using LASSO (Least Absolute
#' Shrinkage and Selection Operator) and/or Elastic Net regression methods with enhanced
#' capabilities. Both methods are regularization techniques that prevent overfitting by
#' adding penalty terms to the loss function. LASSO uses L1 regularization for feature
#' selection, while Elastic Net combines L1 and L2 penalties. The function supports:
#' \itemize{
#'   \item{Binary and multinomial classification (categorical responses)}
#'   \item{Continuous regression (numeric responses)}
#'   \item{Multiple alpha values for comprehensive Elastic Net tuning}
#'   \item{Automatic data source selection based on preprocessing parameters}
#'   \item{Robust error handling with detailed diagnostics}
#' }
#'
#' @param data A list object containing preprocessed data. Must be the output from the
#'   \code{perform_PreprocessingPeakData} function, containing the following elements:
#'   \itemize{
#'     \item{\code{FunctionOrigin}}: Character string indicating data source
#'     \item{\code{Metadata}}: Data frame with sample metadata including Group column
#'     \item{\code{data_scaledPCA_varFiltered}}: Matrix of preprocessed features
#'     \item{\code{data_scaledPCA_merged}}: Matrix of merged replicate features (optional)
#'     \item{\code{Parameters}}: List containing preprocessing parameters (optional)
#'   }
#' @param method Character vector specifying regression method(s) to perform. Options:
#'   \itemize{
#'     \item{\code{"lasso"}}: LASSO regression only (L1 penalty, alpha = 1)
#'     \item{\code{"enet"}}: Elastic Net regression only (L1+L2 penalty)
#'     \item{\code{c("lasso", "enet")}}: Both methods (recommended for comparison)
#'   }
#'   Default: \code{"enet"}
#' @param specify_response Character string specifying the response variable column name.
#'   Can be categorical (classification) or numeric (regression). Supported values:
#'   \itemize{
#'     \item{\code{NULL}}: Uses the Group column from metadata (default)
#'     \item{\code{"Group"}}: Uses the Group column explicitly
#'     \item{\code{"Group2"}}: Uses the Group2 column if available
#'     \item{\code{"Response"}}: Uses the Response column if available
#'     \item{Custom column name}: Any valid column name in metadata
#'   }
#'   Default: \code{NULL}
#' @param train_percent Numeric value between 1 and 99 specifying the percentage of data
#'   to use for training. Remaining data used for testing. Default: \code{80}
#' @param ref Character string specifying the reference level for categorical response
#'   variables. If \code{NULL}, uses the first factor level alphabetically. Ignored for
#'   numeric responses. Default: \code{NULL}
#' @param lambda Character string specifying lambda selection criterion:
#'   \itemize{
#'     \item{\code{"1se"}}: Lambda within one standard error of minimum (conservative, fewer features)
#'     \item{\code{"min"}}: Lambda that minimizes cross-validation error (aggressive, more features)
#'   }
#'   Default: \code{"1se"}
#' @param alpha Numeric value or vector specifying the elastic net mixing parameter(s):
#'   \itemize{
#'     \item{\code{Single value}}: Between 0 (ridge) and 1 (lasso). Default: \code{0.5}
#'     \item{\code{Vector}}: Multiple values for comprehensive tuning, e.g., \code{c(0, 0.25, 0.5, 0.75, 1)}
#'   }
#'   When multiple values provided, all will be evaluated and results stored separately.
#' @param remember Numeric value for reproducible results. Sets random seed using
#'   \code{set.seed(remember)}. If \code{NULL}, no seed is set. Default: \code{123}
#' @param verbose Logical indicating whether to print progress messages and results
#'   to console. Default: \code{TRUE}
#' @param cv_folds Integer specifying number of cross-validation folds for model
#'   selection. Must be between 3 and 20. Default: \code{10}
#' @param parallel Logical indicating whether to use parallel processing for
#'   cross-validation. Default: \code{FALSE}
#' @param standardize Logical indicating whether to standardize features before modeling.
#'   Default: \code{TRUE}
#' @param maxit Integer specifying maximum iterations for model convergence.
#'   Default: \code{100000}
#'
#' @return A list containing regression results with the following structure:
#' \describe{
#'   \item{\code{FunctionOrigin}}{Character string identifying the source function}
#'   \item{\code{ModelSummary}}{Data frame summarizing model performance metrics}
#'   \item{\code{DataSplit}}{List containing training/testing data split information}
#'   \item{\code{DataSource}}{Character string indicating which data matrix was used}
#'   \item{\code{ResponseType}}{Character string indicating "categorical" or "numeric"}
#'   \item{\code{LASSO_Results}}{List of LASSO results (if method includes "lasso")}
#'   \item{\code{ElasticNet_Results}}{List of Elastic Net results (if method includes "enet")}
#'   \item{\code{AlphaComparison}}{Data frame comparing multiple alpha values (if vector provided)}
#'   \item{\code{ComparisonSummary}}{Data frame comparing methods (if both performed)}
#'   \item{\code{ErrorLog}}{List of any warnings or errors encountered}
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Example 1: Categorical response with multiple alpha values
#' regression_results <- perform_Regression(
#'   data = preprocessed_data,
#'   method = c("lasso", "enet"),
#'   specify_response = "Group",
#'   alpha = c(0.1, 0.5, 0.9),
#'   train_percent = 75,
#'   lambda = "1se",
#'   remember = 123
#' )
#'
#' # Example 2: Numeric response
#' regression_results <- perform_Regression(
#'   data = preprocessed_data,
#'   method = "enet",
#'   specify_response = "Response",
#'   alpha = 0.5,
#'   train_percent = 80
#' )
#'
#' # View results
#' print(regression_results$ModelSummary)
#' print(regression_results$AlphaComparison)
#' }
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link[caret]{confusionMatrix}}
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent (2010), Journal of Statistical Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01. (for glmnet)
#' @references Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05. (for glmnet)
#' @references Kuhn, M. (2008), “Building predictive models in R using the caret package, ” Journal of Statistical Software, (doi:10.18637/jss.v028.i05). (for confusionMatrix )
#' @references Altman, D.G., Bland, J.M. (1994) “Diagnostic tests 1: sensitivity and specificity,” British Medical Journal, vol 308, 1552. (for confusionMatrix )
#' @references Altman, D.G., Bland, J.M. (1994) “Diagnostic tests 2: predictive values,” British Medical Journal, vol 309, 102. (for confusionMatrix )
#' @references Velez, D.R., et. al. (2008) “A balanced accuracy function for epistasis modeling in imbalanced datasets using multifactor dimensionality reduction.,” Genetic Epidemiology, vol 4, 306. (for confusionMatrix )
#' @references Chambers, J. M. and Hastie, T. J. (1992) Statistical Models in S. Wadsworth & Brooks/Cole. (for predict)
#'
#' @export
perform_Regression <- function(
    data,
    method = "enet",
    specify_response = NULL,
    train_percent = 80,
    ref = NULL,
    lambda = "1se",
    alpha = 0.5,
    remember = 123,
    verbose = TRUE,
    cv_folds = 10,
    parallel = FALSE,
    standardize = TRUE,
    maxit = 100000
) {

  # Initialize error logging
  error_log <- list()

  # Suppress browser() calls by setting debug options
  old_debug <- getOption("error")
  options(error = NULL)
  on.exit(options(error = old_debug), add = TRUE)

  tryCatch({

    # Validate inputs
    validation_result <- .validate_inputs_regression(data, method, train_percent, lambda, alpha, remember, cv_folds, specify_response)
    if (!validation_result$valid) {
      stop(validation_result$message)
    }

    # Set seed for reproducibility
    if (!is.null(remember)) {
      set.seed(remember)
      if (verbose) cat("Random seed set to:", remember, "\n")
    }

    # Determine data source and prepare data
    data_prep_result <- .prepare_data_regression(data, specify_response, ref, train_percent, verbose)
    if (!data_prep_result$success) {
      stop("Data preparation failed: ", data_prep_result$error)
    }
    data_prep <- data_prep_result$data

    # Initialize results structure
    results <- list(
      FunctionOrigin = "perform_Regression",
      ModelSummary = data.frame(),
      DataSplit = data_prep$split_info,
      DataSource = data_prep$data_source,
      ResponseType = data_prep$response_type,
      ComparisonSummary = NULL,
      AlphaComparison = NULL,
      ErrorLog = error_log
    )

    # Store performance metrics for comparison
    performance_comparison <- list()
    alpha_performance <- data.frame()

    # Handle multiple alpha values
    alpha_values <- if ("lasso" %in% method) c(1, alpha[alpha != 1]) else alpha
    alpha_values <- unique(alpha_values)

    # Perform regression for each method
    for (method_type in method) {

      if (verbose) {
        cat("\n", paste(rep("=", 60), collapse = ""), "\n")
        cat("PERFORMING", toupper(method_type), "REGRESSION\n")
        cat(paste(rep("=", 60), collapse = ""), "\n")
      }

      if (method_type == "lasso") {
        # LASSO with alpha = 1
        method_results <- .fit_regression_model_enhanced(
          x_train = data_prep$x_train,
          y_train = data_prep$y_train,
          x_test = data_prep$x_test,
          y_test = data_prep$y_test,
          method_type = method_type,
          alpha_val = 1,
          lambda = lambda,
          response_info = data_prep$response_info,
          cv_folds = cv_folds,
          parallel = parallel,
          standardize = standardize,
          maxit = maxit,
          verbose = verbose
        )

        results$LASSO_Results <- method_results
        performance_comparison[["lasso"]] <- method_results$Performance

        # Add to alpha comparison
        alpha_row <- data.frame(
          Method = "lasso",
          Alpha = 1,
          method_results$Performance[setdiff(names(method_results$Performance), c("Method", "Alpha"))],
          stringsAsFactors = FALSE
        )
        alpha_performance <- rbind(alpha_performance, alpha_row)

      } else if (method_type == "enet") {
        # Elastic Net with specified alpha value(s)
        enet_results_list <- list()

        for (alpha_val in alpha) {
          if (verbose && length(alpha) > 1) {
            cat("  Testing alpha =", alpha_val, "\n")
          }

          method_results <- .fit_regression_model_enhanced(
            x_train = data_prep$x_train,
            y_train = data_prep$y_train,
            x_test = data_prep$x_test,
            y_test = data_prep$y_test,
            method_type = method_type,
            alpha_val = alpha_val,
            lambda = lambda,
            response_info = data_prep$response_info,
            cv_folds = cv_folds,
            parallel = parallel,
            standardize = standardize,
            maxit = maxit,
            verbose = verbose && length(alpha) == 1
          )

          enet_results_list[[paste0("alpha_", alpha_val)]] <- method_results

          # Add to alpha comparison
          alpha_row <- data.frame(
            Method = "enet",
            Alpha = alpha_val,
            method_results$Performance[setdiff(names(method_results$Performance), c("Method", "Alpha"))],
            stringsAsFactors = FALSE
          )
          alpha_performance <- rbind(alpha_performance, alpha_row)
        }

        # Store best performing alpha for main results
        if (data_prep$response_type == "categorical") {
          best_alpha_idx <- which.max(sapply(enet_results_list, function(x) x$Performance$Accuracy))
        } else {
          best_alpha_idx <- which.min(sapply(enet_results_list, function(x) x$Performance$RMSE))
        }

        results$ElasticNet_Results <- enet_results_list[[best_alpha_idx]]
        results$ElasticNet_Results$AllAlphaResults <- enet_results_list
        performance_comparison[["enet"]] <- results$ElasticNet_Results$Performance

        # Create coefficient dimensions summary for all alphas
        if (method_type == "enet" && length(alpha) > 1) {
          coef_dims_list <- lapply(enet_results_list, function(x) {
            data.frame(
              Alpha = x$Alpha,
              n_coefficients = nrow(x$Coefficients),
              stringsAsFactors = FALSE
            )
          })

          results$ElasticNet_Results$CoefficientDimensions <- do.call(rbind, coef_dims_list)
          rownames(results$ElasticNet_Results$CoefficientDimensions) <- NULL
        }
      }
    }

    # Create alpha comparison summary
    if (nrow(alpha_performance) > 1) {
      results$AlphaComparison <- alpha_performance[order(alpha_performance$Method, alpha_performance$Alpha), ]
      rownames(results$AlphaComparison) <- NULL
    }

    # Create model summary
    results$ModelSummary <- .create_model_summary_enhanced(performance_comparison, method, data_prep$response_type)

    # Create comparison summary if both methods were used
    if (length(method) > 1) {
      results$ComparisonSummary <- .create_comparison_summary_enhanced(performance_comparison, data_prep$response_type)
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
      if (!is.null(results$AlphaComparison)) {
        cat("\nALPHA COMPARISON:\n")
        print(results$AlphaComparison)
      }
    }

    return(results)

  }, error = function(e) {
    error_msg <- paste("Error in perform_Regression:", e$message)
    if (verbose) {
      cat("\n", paste(rep("!", 60), collapse = ""), "\n")
      cat("ERROR OCCURRED:\n")
      cat(error_msg, "\n")
      cat(paste(rep("!", 60), collapse = ""), "\n")
    }

    # Return partial results with error information
    return(list(
      FunctionOrigin = "perform_Regression",
      Error = error_msg,
      ErrorDetails = as.character(e),
      ErrorLog = error_log,
      Success = FALSE
    ))
  })
}

# Enhanced input validation
.validate_inputs_regression <- function(data, method, train_percent, lambda, alpha, remember, cv_folds, specify_response) {

  # Check data structure
  if (!is.list(data)) {
    return(list(valid = FALSE, message = "'data' must be a list object from perform_PreprocessingPeakData function."))
  }

  if (!"FunctionOrigin" %in% names(data) || data$FunctionOrigin != "perform_PreprocessingPeakData") {
    return(list(valid = FALSE, message = "The supplied data did not come from 'perform_PreprocessingPeakData' function."))
  }

  required_elements <- c("Metadata")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    return(list(valid = FALSE, message = paste("Missing required data elements:", paste(missing_elements, collapse = ", "))))
  }

  # Check for at least one data matrix
  if (!("data_scaledPCA_varFiltered" %in% names(data)) && !("data_scaledPCA_merged" %in% names(data))) {
    return(list(valid = FALSE, message = "No valid data matrix found. Need either 'data_scaledPCA_varFiltered' or 'data_scaledPCA_merged'."))
  }

  # Check method
  allowed_methods <- c("lasso", "enet")
  if (!all(method %in% allowed_methods)) {
    return(list(valid = FALSE, message = paste("'method' must be one or more of:", paste(allowed_methods, collapse = ", "))))
  }

  # Check train_percent
  if (!is.numeric(train_percent) || train_percent <= 0 || train_percent >= 100) {
    return(list(valid = FALSE, message = "'train_percent' must be a numeric value between 1 and 99."))
  }

  # Check lambda
  if (!lambda %in% c("min", "1se")) {
    return(list(valid = FALSE, message = "'lambda' must be either 'min' or '1se'."))
  }

  # Check alpha
  if (!is.numeric(alpha) || any(alpha < 0) || any(alpha > 1)) {
    return(list(valid = FALSE, message = "'alpha' must be numeric value(s) between 0 and 1."))
  }

  # Check remember
  if (!is.null(remember) && !is.numeric(remember)) {
    return(list(valid = FALSE, message = "'remember' parameter must be numeric or NULL."))
  }

  # Check cv_folds
  if (!is.numeric(cv_folds) || cv_folds < 3 || cv_folds > 20) {
    return(list(valid = FALSE, message = "'cv_folds' must be an integer between 3 and 20."))
  }

  # Check specify_response if provided
  if (!is.null(specify_response)) {
    if (!is.character(specify_response) || length(specify_response) != 1) {
      return(list(valid = FALSE, message = "'specify_response' must be a single character string."))
    }
    if (!specify_response %in% colnames(data$Metadata)) {
      return(list(valid = FALSE, message = paste("Specified response variable '", specify_response, "' not found in metadata.")))
    }
  }

  return(list(valid = TRUE, message = "All inputs valid"))
}

# Enhanced data preparation function
.prepare_data_regression <- function(data, specify_response, ref, train_percent, verbose) {

  tryCatch({

    # Determine which data matrix to use
    data_source <- "data_scaledPCA_varFiltered"  # default

    if ("Parameters" %in% names(data) &&
        "auto_merge_replicates" %in% names(data$Parameters) &&
        data$Parameters$auto_merge_replicates == TRUE &&
        "data_scaledPCA_merged" %in% names(data)) {

      if (verbose) cat("Using merged replicate data: data_scaledPCA_merged\n")
      features <- data$data_scaledPCA_merged
      data_source <- "data_scaledPCA_merged"

      # Extract metadata
      metadata <- data$Metadata_merged

    } else {

      if (verbose) cat("Using filtered data: data_scaledPCA_varFiltered\n")
      features <- data$data_scaledPCA_varFiltered
      data_source <- "data_scaledPCA_varFiltered"

      # Extract metadata
      metadata <- data$Metadata
    }

    # Convert features to matrix if not already
    if (!is.matrix(features)) {
      features <- as.matrix(features)
    }

    # Identify QC samples
    qc_indices <- metadata$Group %in% c("SQC", "EQC", "QC")
    non_qc_indices <- !qc_indices

    if (sum(non_qc_indices) == 0) {
      return(list(success = FALSE, error = "No non-QC samples found in the data."))
    }

    # Prepare response variable
    if (!is.null(specify_response)) {
      if (!specify_response %in% colnames(metadata)) {
        return(list(success = FALSE, error = paste("Specified response variable '", specify_response, "' not found in metadata.")))
      }
      response_var <- metadata[[specify_response]][non_qc_indices]
    } else {
      response_var <- metadata$Group[non_qc_indices]
    }

    # Remove NA values
    na_indices <- is.na(response_var)
    if (any(na_indices)) {
      if (verbose) cat("Removing", sum(na_indices), "samples with NA response values.\n")
      non_qc_indices[which(non_qc_indices)[na_indices]] <- FALSE
      response_var <- response_var[!na_indices]
    }

    # Determine response type and prepare accordingly
    response_type <- "categorical"
    response_info <- list()

    if (is.numeric(response_var)) {
      response_type <- "numeric"
      response_info$type <- "numeric"
      response_info$min <- min(response_var, na.rm = TRUE)
      response_info$max <- max(response_var, na.rm = TRUE)
      response_info$mean <- mean(response_var, na.rm = TRUE)
      response_info$sd <- sd(response_var, na.rm = TRUE)

      if (verbose) {
        cat("Response variable is numeric:\n")
        cat("  Range:", round(response_info$min, 4), "to", round(response_info$max, 4), "\n")
        cat("  Mean (SD):", round(response_info$mean, 4), "(", round(response_info$sd, 4), ")\n")
      }

    } else {
      # Categorical response
      response_type <- "categorical"

      # Set reference level
      ref_level <- NULL
      if (!is.null(ref)) {
        if (!ref %in% unique(response_var)) {
          return(list(success = FALSE, error = paste("Reference level '", ref, "' not found in response variable.")))
        }
        ref_level <- ref
        response_var <- factor(response_var, levels = c(ref, setdiff(unique(response_var), ref)))
      } else {
        response_var <- factor(response_var)
        ref_level <- levels(response_var)[1]
      }

      # Check for sufficient samples per group
      group_counts <- table(response_var)
      if (any(group_counts < 3)) {
        if (verbose) cat("Warning: Some groups have fewer than 3 samples. This may affect model performance.\n")
      }

      response_info$type <- "categorical"
      response_info$levels <- levels(response_var)
      response_info$ref_level <- ref_level
      response_info$group_counts <- as.vector(group_counts)
      names(response_info$group_counts) <- names(group_counts)

      if (verbose) {
        cat("Response variable is categorical:\n")
        cat("  Groups:", paste(levels(response_var), collapse = ", "), "\n")
        cat("  Reference level:", ref_level, "\n")
        cat("  Group distribution:\n")
        print(group_counts)
      }
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
      if (verbose) cat("Removing", sum(constant_features), "constant or near-constant features.\n")
      df <- df[, !constant_features, drop = FALSE]
    }

    # Check for missing values in features
    if (any(is.na(df))) {
      if (verbose) cat("Warning: Missing values detected in features. Removing incomplete cases.\n")
      complete_cases <- complete.cases(df)
      if (sum(complete_cases) < nrow(df)) {
        if (verbose) cat("Removing", sum(!complete_cases), "samples with missing feature values.\n")
        df <- df[complete_cases, , drop = FALSE]
        response_var <- response_var[complete_cases]
      }
    }

    # Split data
    n_samples <- nrow(df)
    train_size <- floor((train_percent / 100) * n_samples)

    if (train_size < 3) {
      return(list(success = FALSE, error = "Training set too small. Need at least 3 samples."))
    }

    # Stratified sampling for categorical, random for numeric
    if (response_type == "categorical") {
      if (train_size < length(levels(response_var))) {
        return(list(success = FALSE, error = "Training set too small. Reduce train_percent or increase sample size."))
      }

      # Stratified sampling
      train_indices <- c()
      for (level in levels(response_var)) {
        level_indices <- which(response_var == level)
        level_train_size <- max(1, floor(length(level_indices) * train_percent / 100))
        level_train_indices <- sample(level_indices, level_train_size)
        train_indices <- c(train_indices, level_train_indices)
      }

      # Add more samples if needed
      remaining_indices <- setdiff(seq_len(n_samples), train_indices)
      if (length(train_indices) < train_size && length(remaining_indices) > 0) {
        additional_needed <- min(train_size - length(train_indices), length(remaining_indices))
        additional_indices <- sample(remaining_indices, additional_needed)
        train_indices <- c(train_indices, additional_indices)
      }

    } else {
      # Random sampling for numeric response
      train_indices <- sample(seq_len(n_samples), train_size)
    }

    x_train <- df[train_indices, , drop = FALSE]
    x_test <- df[-train_indices, , drop = FALSE]
    y_train <- response_var[train_indices]
    y_test <- response_var[-train_indices]

    # Final validation
    if (nrow(x_train) == 0 || nrow(x_test) == 0) {
      return(list(success = FALSE, error = "Data split resulted in empty training or testing set."))
    }

    if (ncol(x_train) == 0) {
      return(list(success = FALSE, error = "No features remaining after preprocessing."))
    }

    # Ensure factor levels are consistent for categorical
    if (response_type == "categorical") {
      y_train <- factor(y_train, levels = levels(response_var))
      y_test <- factor(y_test, levels = levels(response_var))
    }

    if (verbose) {
      cat("Data preparation complete:\n")
      cat("  Data source:", data_source, "\n")
      cat("  Response type:", response_type, "\n")
      cat("  Total samples:", n_samples, "\n")
      cat("  Training samples:", nrow(x_train), "\n")
      cat("  Testing samples:", nrow(x_test), "\n")
      cat("  Features:", ncol(x_train), "\n")

      if (response_type == "categorical") {
        cat("  Training group distribution:\n")
        print(table(y_train))
        cat("  Testing group distribution:\n")
        print(table(y_test))
      }
    }

    return(list(
      success = TRUE,
      data = list(
        x_train = x_train,
        x_test = x_test,
        y_train = y_train,
        y_test = y_test,
        response_info = response_info,
        response_type = response_type,
        data_source = data_source,
        split_info = list(
          total_samples = n_samples,
          train_samples = nrow(x_train),
          test_samples = nrow(x_test),
          features = ncol(x_train),
          response_type = response_type,
          train_indices = train_indices
        )
      )
    ))

  }, error = function(e) {
    return(list(success = FALSE, error = paste("Data preparation error:", e$message)))
  })
}

# Enhanced model fitting function
.fit_regression_model_enhanced <- function(x_train, y_train, x_test, y_test, method_type,
                                         alpha_val, lambda, response_info, cv_folds,
                                         parallel, standardize, maxit, verbose) {

  tryCatch({

    # Determine family and type.measure based on response type
    if (response_info$type == "categorical") {
      n_groups <- length(response_info$levels)
      family_type <- if (n_groups == 2) "binomial" else "multinomial"
      type_measure <- "class"
    } else {
      family_type <- "gaussian"
      type_measure <- "mse"
    }

    if (verbose) {
      cat("  Method:", method_type, "\n")
      cat("  Alpha:", alpha_val, "\n")
      cat("  Family:", family_type, "\n")
      cat("  Type measure:", type_measure, "\n")
      cat("  CV folds:", cv_folds, "\n")
    }

    # Additional data validation
    if (any(is.na(x_train)) || any(is.na(y_train))) {
      stop("Training data contains missing values.")
    }

    if (any(is.na(x_test)) || any(is.na(y_test))) {
      stop("Testing data contains missing values.")
    }

    # Check for adequate sample size per fold
    if (response_info$type == "categorical") {
      min_group_size <- min(table(y_train))
      if (cv_folds > min_group_size) {
        cv_folds <- max(3, min_group_size - 1)
        if (verbose) cat("  Reduced cv_folds to", cv_folds, "due to small group sizes.\n")
      }
    }

    # Fit model
    model_fit <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      alpha = alpha_val,
      family = family_type,
      type.measure = type_measure,
      nfolds = cv_folds,
      parallel = parallel,
      standardize = standardize,
      maxit = maxit
    )

    # Select lambda
    selected_lambda <- if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min

    if (verbose) {
      cat("  Selected lambda:", selected_lambda, "\n")
    }

    # Make predictions
    if (response_info$type == "categorical") {
      predictions <- predict(
        model_fit,
        s = selected_lambda,
        newx = x_test,
        type = "class"
      )
      predictions <- as.vector(predictions)
      predictions <- factor(predictions, levels = response_info$levels)
    } else {
      predictions <- predict(
        model_fit,
        s = selected_lambda,
        newx = x_test,
        type = "response"
      )
      predictions <- as.numeric(predictions)
    }

    # Create performance metrics
    if (response_info$type == "categorical") {
      # Confusion matrix for classification
      conf_matrix <- caret::confusionMatrix(
        predictions,
        factor(y_test, levels = response_info$levels)
      )

      performance <- .create_performance_summary_enhanced(
        conf_matrix = conf_matrix,
        method = method_type,
        alpha = alpha_val,
        lambda = selected_lambda,
        ref = response_info$ref_level,
        response_type = "categorical"
      )

      if (verbose) {
        cat("  Accuracy:", round(performance$Accuracy, 4), "\n")
      }

    } else {
      # Regression metrics
      conf_matrix <- NULL

      # Calculate regression metrics
      mse <- mean((y_test - predictions)^2)
      rmse <- sqrt(mse)
      mae <- mean(abs(y_test - predictions))

      # R-squared
      ss_res <- sum((y_test - predictions)^2)
      ss_tot <- sum((y_test - mean(y_test))^2)
      r_squared <- 1 - (ss_res / ss_tot)

      # Correlation
      correlation <- cor(y_test, predictions, use = "complete.obs")

      performance <- data.frame(
        Method = method_type,
        Alpha = alpha_val,
        Lambda = selected_lambda,
        MSE = mse,
        RMSE = rmse,
        MAE = mae,
        R_squared = r_squared,
        Correlation = correlation,
        stringsAsFactors = FALSE
      )

      if (verbose) {
        cat("  RMSE:", round(rmse, 4), "\n")
        cat("  R-squared:", round(r_squared, 4), "\n")
      }
    }

    # Extract coefficients
    coef_results <- .extract_coefficients_enhanced(
      model_fit,
      selected_lambda,
      response_info,
      family_type
    )

    # Create coefficient dimensions summary
    coef_dimensions <- data.frame(
      Alpha = alpha_val,
      n_coefficients = nrow(coef_results),
      stringsAsFactors = FALSE
    )

    if (verbose) {
      cat("  Non-zero coefficients:", nrow(coef_results), "\n")
    }

    return(list(
      Model = model_fit,
      Predictions = data.frame(
        Actual = if (response_info$type == "categorical") as.character(y_test) else y_test,
        Predicted = if (response_info$type == "categorical") as.character(predictions) else predictions,
        stringsAsFactors = FALSE
      ),
      ConfusionMatrix = conf_matrix,
      Performance = performance,
      Coefficients = coef_results,
      CoefficientDimensions = coef_dimensions,
      Lambda = selected_lambda,
      Alpha = alpha_val,
      ResponseType = response_info$type,
      ReferenceLevel = if (response_info$type == "categorical") response_info$ref_level else NULL
    ))

  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })
}

# Enhanced coefficient extraction
.extract_coefficients_enhanced <- function(model_fit, lambda, response_info, family_type) {

  coef_model <- glmnet::coef.glmnet(model_fit, s = lambda)

  if (family_type == "binomial") {
    # Binary classification
    coef_df <- as.data.frame(as.matrix(coef_model))
    colnames(coef_df) <- "Coefficient"
    coef_df$Feature <- rownames(coef_df)
    coef_df <- coef_df[coef_df$Coefficient != 0, ]
    coef_df$OddsRatio <- exp(coef_df$Coefficient)
    coef_df <- coef_df[, c("Feature", "Coefficient", "OddsRatio")]

  } else if (family_type == "multinomial") {
    # Multinomial classification
    group_levels <- response_info$levels

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

  } else {
    # Gaussian (regression)
    coef_df <- as.data.frame(as.matrix(coef_model))
    colnames(coef_df) <- "Coefficient"
    coef_df$Feature <- rownames(coef_df)
    coef_df <- coef_df[coef_df$Coefficient != 0, ]
    coef_df <- coef_df[, c("Feature", "Coefficient")]
  }

  rownames(coef_df) <- NULL
  return(coef_df)
}

# Enhanced performance summary creation
.create_performance_summary_enhanced <- function(conf_matrix = NULL, method, alpha, lambda, ref = NULL, response_type) {

  if (response_type == "categorical") {
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
  }

  rownames(performance) <- NULL
  return(performance)
}

# Enhanced model summary creation
.create_model_summary_enhanced <- function(performance_list, methods, response_type) {

  summary_df <- do.call(rbind, performance_list)
  rownames(summary_df) <- NULL

  # Add ranking based on primary metric
  if (response_type == "categorical") {
    summary_df$AccuracyRank <- rank(-summary_df$Accuracy, ties.method = "min")
  } else {
    summary_df$RMSERank <- rank(summary_df$RMSE, ties.method = "min")
  }

  return(summary_df)
}

# Enhanced comparison summary creation
.create_comparison_summary_enhanced <- function(performance_list, response_type) {

  if (length(performance_list) < 2) {
    return(NULL)
  }

  methods <- names(performance_list)

  if (response_type == "categorical") {
    metrics <- c("Accuracy", "Kappa")
  } else {
    metrics <- c("RMSE", "MAE", "R_squared", "Correlation")
  }

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

#' Extract Results from Regularized Regression Analysis
#'
#' @title Extract Specific Results from perform_Regression Output
#'
#' @description
#' This function extracts specific results from the output of \code{perform_Regression}.
#' It can extract all results for a specific alpha value or specific components
#' (like coefficients, performance metrics, etc.) for targeted analysis.
#' Supports both single method results and multi-alpha comparisons.
#'
#' @param results A list object containing results from \code{perform_Regression}.
#'   Must contain either \code{LASSO_Results}, \code{ElasticNet_Results}, or both.
#' @param method Character string specifying which method results to extract.
#'   Options: \code{"lasso"}, \code{"enet"}, or \code{"both"}. If \code{"both"},
#'   extracts from both methods if available. Default: \code{"enet"}
#' @param alpha Numeric value specifying which alpha value to extract results for.
#'   Only applicable when the method has multiple alpha results stored in
#'   \code{AllAlphaResults}. If \code{NULL}, extracts the best-performing alpha.
#'   Default: \code{NULL}
#' @param component Character string specifying which component to extract.
#'   Options:
#'   \itemize{
#'     \item{\code{NULL}}: Extract all components for the specified alpha
#'     \item{\code{"model"}}: Extract the fitted cv.glmnet model object
#'     \item{\code{"predictions"}}: Extract actual vs predicted values
#'     \item{\code{"confusion"}}: Extract confusion matrix (classification only)
#'     \item{\code{"performance"}}: Extract performance metrics
#'     \item{\code{"coefficients"}}: Extract non-zero coefficients and odds ratios
#'     \item{\code{"lambda"}}: Extract selected lambda value
#'     \item{\code{"summary"}}: Extract key summary information
#'   }
#'   Default: \code{NULL}
#' @param alpha_comparison Logical indicating whether to return the alpha comparison
#'   data frame instead of specific method results. Overrides other parameters
#'   except \code{component}. Default: \code{FALSE}
#' @param format Character string specifying output format for certain components.
#'   Options: \code{"list"}, \code{"dataframe"}, \code{"matrix"}. Default: \code{"list"}
#' @param verbose Logical indicating whether to print extraction details.
#'   Default: \code{TRUE}
#'
#' @return Depends on the parameters:
#' \describe{
#'   \item{All components (\code{component = NULL})}{List containing all result components}
#'   \item{Specific component}{The requested component (model, data frame, matrix, etc.)}
#'   \item{Alpha comparison (\code{alpha_comparison = TRUE})}{Data frame with alpha comparisons}
#'   \item{Both methods (\code{method = "both"})}{List with results from both methods}
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Assume results_enet_1se_model1 is your regression results
#'
#' # Extract all results for alpha = 0.9
#' alpha_09_results <- extract_RegressionResults(
#'   results = results_enet_1se_model1,
#'   method = "enet",
#'   alpha = 0.9
#' )
#'
#' # Extract only coefficients for alpha = 0.9
#' coefficients_09 <- extract_RegressionResults(
#'   results = results_enet_1se_model1,
#'   method = "enet",
#'   alpha = 0.9,
#'   component = "coefficients"
#' )
#'
#' # Extract performance metrics for best alpha
#' performance_best <- extract_RegressionResults(
#'   results = results_enet_1se_model1,
#'   method = "enet",
#'   component = "performance"
#' )
#'
#' # Extract alpha comparison table
#' alpha_comparison <- extract_RegressionResults(
#'   results = results_enet_1se_model1,
#'   alpha_comparison = TRUE
#' )
#'
#' # Extract predictions for alpha = 0.5 as data frame
#' predictions_df <- extract_RegressionResults(
#'   results = results_enet_1se_model1,
#'   method = "enet",
#'   alpha = 0.5,
#'   component = "predictions",
#'   format = "dataframe"
#' )
#' }
#'
#' @seealso \code{\link{perform_Regression}}
#'
#' @export
extract_RegressionResults <- function(
    results,
    method = "enet",
    alpha = NULL,
    component = NULL,
    alpha_comparison = FALSE,
    format = "list",
    verbose = TRUE
) {

  # Validate inputs
  validation_result <- .validate_extraction_inputs(results, method, component, format)
  if (!validation_result$valid) {
    stop(validation_result$message)
  }

  # Handle alpha comparison request
  if (alpha_comparison) {
    if (verbose) cat("Extracting alpha comparison data...\n")

    if ("AlphaComparison" %in% names(results) && !is.null(results$AlphaComparison)) {
      if (!is.null(component)) {
        if (component %in% colnames(results$AlphaComparison)) {
          return(results$AlphaComparison[[component]])
        } else {
          stop("Component '", component, "' not found in AlphaComparison.")
        }
      }
      return(results$AlphaComparison)
    } else {
      stop("No AlphaComparison found in results. This may be a single-alpha analysis.")
    }
  }

  # Handle both methods request
  if (method == "both") {
    if (verbose) cat("Extracting results from both methods...\n")

    both_results <- list()

    # Extract LASSO if available
    if ("LASSO_Results" %in% names(results)) {
      both_results$LASSO <- .extract_single_method(
        results$LASSO_Results, "lasso", alpha, component, format, verbose
      )
    }

    # Extract Elastic Net if available
    if ("ElasticNet_Results" %in% names(results)) {
      both_results$ElasticNet <- .extract_single_method(
        results$ElasticNet_Results, "enet", alpha, component, format, verbose
      )
    }

    if (length(both_results) == 0) {
      stop("No method results found in the provided results object.")
    }

    return(both_results)
  }

  # Handle single method request
  method_key <- if (method == "lasso") "LASSO_Results" else "ElasticNet_Results"

  if (!method_key %in% names(results)) {
    stop("Method '", method, "' results not found. Available methods: ",
         paste(gsub("_Results", "", names(results)[grepl("_Results$", names(results))]), collapse = ", "))
  }

  method_results <- results[[method_key]]

  return(.extract_single_method(method_results, method, alpha, component, format, verbose))
}

# Helper function to validate extraction inputs
.validate_extraction_inputs <- function(results, method, component, format) {

  # Check results structure
  if (!is.list(results)) {
    return(list(valid = FALSE, message = "'results' must be a list object from perform_Regression."))
  }

  if (!"FunctionOrigin" %in% names(results) || results$FunctionOrigin != "perform_Regression") {
    return(list(valid = FALSE, message = "Results object does not appear to be from perform_Regression function."))
  }

  # Check method
  if (!method %in% c("lasso", "enet", "both")) {
    return(list(valid = FALSE, message = "'method' must be one of: 'lasso', 'enet', 'both'"))
  }

  # Check component
  valid_components <- c("model", "predictions", "confusion", "performance",
                        "coefficients", "lambda", "summary")
  if (!is.null(component) && !component %in% valid_components) {
    return(list(valid = FALSE, message = paste("'component' must be one of:",
                                               paste(valid_components, collapse = ", "))))
  }

  # Check format
  if (!format %in% c("list", "dataframe", "matrix")) {
    return(list(valid = FALSE, message = "'format' must be one of: 'list', 'dataframe', 'matrix'"))
  }

  return(list(valid = TRUE, message = "All inputs valid"))
}

# Helper function to extract from a single method
.extract_single_method <- function(method_results, method_name, alpha, component, format, verbose) {

  # Determine if we need to extract from AllAlphaResults
  if (!is.null(alpha) && "AllAlphaResults" %in% names(method_results)) {

    if (verbose) cat("Extracting results for", method_name, "with alpha =", alpha, "...\n")

    # Find the matching alpha result
    alpha_key <- paste0("alpha_", alpha)

    if (!alpha_key %in% names(method_results$AllAlphaResults)) {
      available_alphas <- gsub("alpha_", "", names(method_results$AllAlphaResults))
      stop("Alpha value ", alpha, " not found. Available alphas: ",
           paste(available_alphas, collapse = ", "))
    }

    target_results <- method_results$AllAlphaResults[[alpha_key]]

  } else {
    # Use the main results (best performing alpha or single alpha)
    if (verbose && is.null(alpha)) {
      cat("Extracting best-performing", method_name, "results...\n")
    } else if (verbose) {
      cat("Extracting", method_name, "results (single alpha analysis)...\n")
    }

    target_results <- method_results
  }

  # Extract specific component if requested
  if (!is.null(component)) {
    return(.extract_component(target_results, component, format, verbose))
  }

  # Return all results
  return(target_results)
}

# Helper function to extract specific component
.extract_component <- function(results, component, format, verbose) {

  switch(component,
         "model" = {
           if (verbose) cat("Extracting model object...\n")
           if ("Model" %in% names(results)) {
             return(results$Model)
           } else {
             stop("Model not found in results.")
           }
         },

         "predictions" = {
           if (verbose) cat("Extracting predictions...\n")
           if ("Predictions" %in% names(results)) {
             pred_data <- results$Predictions
             return(.format_output(pred_data, format))
           } else {
             stop("Predictions not found in results.")
           }
         },

         "confusion" = {
           if (verbose) cat("Extracting confusion matrix...\n")
           if ("ConfusionMatrix" %in% names(results) && !is.null(results$ConfusionMatrix)) {
             if (format == "matrix") {
               return(results$ConfusionMatrix$table)
             } else if (format == "dataframe") {
               conf_table <- as.data.frame(results$ConfusionMatrix$table)
               return(conf_table)
             } else {
               return(results$ConfusionMatrix)
             }
           } else {
             stop("Confusion matrix not found in results (may be numeric response).")
           }
         },

         "performance" = {
           if (verbose) cat("Extracting performance metrics...\n")
           if ("Performance" %in% names(results)) {
             perf_data <- results$Performance
             return(.format_output(perf_data, format))
           } else {
             stop("Performance metrics not found in results.")
           }
         },

         "coefficients" = {
           if (verbose) cat("Extracting coefficients...\n")
           if ("Coefficients" %in% names(results)) {
             coef_data <- results$Coefficients
             return(.format_output(coef_data, format))
           } else {
             stop("Coefficients not found in results.")
           }
         },

         "lambda" = {
           if (verbose) cat("Extracting lambda value...\n")
           if ("Lambda" %in% names(results)) {
             return(results$Lambda)
           } else {
             stop("Lambda not found in results.")
           }
         },

         "summary" = {
           if (verbose) cat("Extracting summary information...\n")
           return(.create_extraction_summary(results))
         },

         stop("Unknown component: ", component)
  )
}

# Helper function to format output
.format_output <- function(data, format) {

  if (is.null(data)) return(NULL)

  switch(format,
         "list" = {
           if (is.data.frame(data)) {
             return(as.list(data))
           } else {
             return(data)
           }
         },

         "dataframe" = {
           if (is.data.frame(data)) {
             return(data)
           } else if (is.matrix(data)) {
             return(as.data.frame(data))
           } else if (is.list(data)) {
             return(data.frame(data, stringsAsFactors = FALSE))
           } else {
             return(data.frame(Value = data, stringsAsFactors = FALSE))
           }
         },

         "matrix" = {
           if (is.matrix(data)) {
             return(data)
           } else if (is.data.frame(data)) {
             return(as.matrix(data))
           } else {
             return(as.matrix(data))
           }
         },

         data
  )
}

# Helper function to create extraction summary
.create_extraction_summary <- function(results) {

  summary_list <- list()

  # Basic information
  if ("Alpha" %in% names(results)) summary_list$Alpha <- results$Alpha
  if ("Lambda" %in% names(results)) summary_list$Lambda <- results$Lambda
  if ("ResponseType" %in% names(results)) summary_list$ResponseType <- results$ResponseType
  if ("ReferenceLevel" %in% names(results)) summary_list$ReferenceLevel <- results$ReferenceLevel

  # Performance summary
  if ("Performance" %in% names(results)) {
    perf <- results$Performance
    if ("Accuracy" %in% names(perf)) summary_list$Accuracy <- perf$Accuracy
    if ("RMSE" %in% names(perf)) summary_list$RMSE <- perf$RMSE
    if ("R_squared" %in% names(perf)) summary_list$R_squared <- perf$R_squared
    if ("Kappa" %in% names(perf)) summary_list$Kappa <- perf$Kappa
  }

  # Model characteristics
  if ("Coefficients" %in% names(results)) {
    summary_list$NonZeroCoefficients <- nrow(results$Coefficients)
    if (nrow(results$Coefficients) > 1) {
      summary_list$TopFeatures <- head(results$Coefficients$Feature[
        results$Coefficients$Feature != "(Intercept)"], 5)
    }
  }

  # Prediction summary
  if ("Predictions" %in% names(results)) {
    pred <- results$Predictions
    if ("Actual" %in% names(pred) && "Predicted" %in% names(pred)) {
      if (is.numeric(pred$Actual)) {
        summary_list$PredictionRange <- range(pred$Predicted, na.rm = TRUE)
        summary_list$ActualRange <- range(pred$Actual, na.rm = TRUE)
      } else {
        summary_list$CorrectPredictions <- sum(pred$Actual == pred$Predicted)
        summary_list$TotalPredictions <- nrow(pred)
      }
    }
  }

  return(summary_list)
}

# Convenience function to list all available alphas
list_available_alphas <- function(results, method = "enet") {

  method_key <- if (method == "lasso") "LASSO_Results" else "ElasticNet_Results"

  if (!method_key %in% names(results)) {
    stop("Method '", method, "' results not found.")
  }

  method_results <- results[[method_key]]

  if ("AllAlphaResults" %in% names(method_results)) {
    alpha_keys <- names(method_results$AllAlphaResults)
    alphas <- as.numeric(gsub("alpha_", "", alpha_keys))
    return(sort(alphas))
  } else {
    if ("Alpha" %in% names(method_results)) {
      return(method_results$Alpha)
    } else {
      return(NULL)
    }
  }
}

# Convenience function to get the best performing alpha
get_best_alpha <- function(results, method = "enet", metric = "auto") {

  if ("AlphaComparison" %in% names(results) && !is.null(results$AlphaComparison)) {
    alpha_comp <- results$AlphaComparison
    method_data <- alpha_comp[alpha_comp$Method == method, ]

    if (nrow(method_data) == 0) {
      stop("No data found for method '", method, "'")
    }

    # Auto-select metric based on response type
    if (metric == "auto") {
      if ("Accuracy" %in% names(method_data)) {
        metric <- "Accuracy"
        best_idx <- which.max(method_data[[metric]])
      } else if ("RMSE" %in% names(method_data)) {
        metric <- "RMSE"
        best_idx <- which.min(method_data[[metric]])
      } else {
        stop("Cannot determine appropriate metric automatically")
      }
    } else {
      if (!metric %in% names(method_data)) {
        stop("Metric '", metric, "' not found in alpha comparison")
      }

      if (metric %in% c("RMSE", "MAE", "MSE")) {
        best_idx <- which.min(method_data[[metric]])
      } else {
        best_idx <- which.max(method_data[[metric]])
      }
    }

    best_alpha <- method_data$Alpha[best_idx]
    return(list(alpha = best_alpha, metric = metric, value = method_data[[metric]][best_idx]))
  } else {
    stop("No alpha comparison available. This may be a single-alpha analysis.")
  }
}
