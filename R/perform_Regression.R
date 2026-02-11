#' Perform Regularized Regression Analysis
#'
#' @title Perform Regularized Regression with Cross-Validation for Categorical and Numeric Responses
#'
#' @description
#' This function performs regularized regression analysis using Elastic Net regression,
#' which includes LASSO (alpha = 1) and Ridge (alpha = 0) as special cases.
#' Regularization techniques prevent overfitting by adding penalty terms to the
#' loss function. The function supports:
#' \itemize{
#'   \item{Binary and multinomial classification (categorical responses)}
#'   \item{Continuous regression (numeric responses)}
#'   \item{Multiple alpha values for comprehensive parameter tuning}
#'   \item{Automatic data source selection based on preprocessing parameters}
#'   \item{Robust error handling with detailed diagnostics}
#' }
#'
#' The function properly handles multinomial models by normalizing coefficients
#' to the reference group, making all comparisons interpretable as odds ratios
#' relative to the baseline. Note: P-values are not provided as they are
#' statistically invalid after variable selection in regularized regression.
#'
#' @param data A list object containing preprocessed data. Must be the output from the
#'   \code{perform_PreprocessingPeakData} function, containing the following elements:
#'   \itemize{
#'     \item{\code{FunctionOrigin}}: Character string indicating data source
#'     \item{\code{Metadata}}: Data frame with sample metadata including Group column
#'     \item{\code{data_scaledNONPLS_varFiltered}}: Matrix of preprocessed features
#'     \item{\code{data_scaledNONPLS_merged}}: Matrix of merged replicate features (optional)
#'     \item{\code{Parameters}}: List containing preprocessing parameters (optional)
#'   }
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
#'     \item{\code{0}}: Ridge regression (L2 penalty only)
#'     \item{\code{1}}: LASSO regression (L1 penalty only)
#'     \item{\code{0 < alpha < 1}}: Elastic Net (combination of L1 and L2)
#'     \item{\code{Single value}}: Tests one specific alpha value. Default: \code{0.5}
#'     \item{\code{Vector}}: Tests multiple alpha values, e.g., \code{c(0, 0.25, 0.5, 0.75, 1)}
#'   }
#'   When multiple values are provided, all will be evaluated and the best-performing
#'   model will be selected automatically based on cross-validation performance.
#' @param remember Numeric value for reproducible results. Sets random seed using
#'   \code{set.seed(remember)}. If \code{NULL}, no seed is set. Default: \code{123}
#' @param verbose Logical indicating whether to print progress messages and results
#'   to console. Default: \code{TRUE}
#' @param cv_folds Integer specifying number of cross-validation folds for model
#'   selection. Must be between 3 and 20. Default: \code{10}
#' @param parallel Logical indicating whether to use parallel processing for
#'   cross-validation. Requires a parallel backend to be registered (e.g., via
#'   \code{doParallel::registerDoParallel()}). Default: \code{FALSE}
#' @param standardize Logical indicating whether to standardize features before modeling.
#'   Default: \code{TRUE}
#' @param maxit Integer specifying maximum iterations for model convergence.
#'   Default: \code{100000}
#' @param type_measure Character string specifying the loss function to use for
#'   cross-validation. The default value depends on the response type:
#'   \itemize{
#'     \item \code{NULL}: Automatically selects based on response type (default)
#'     \item \strong{For categorical responses:}
#'       \itemize{
#'         \item \code{"class"}: Misclassification error (default for categorical)
#'         \item \code{"auc"}: Area under ROC curve (binomial only)
#'         \item \code{"deviance"}: Multinomial deviance
#'       }
#'     \item \strong{For numeric responses:}
#'       \itemize{
#'         \item \code{"mse"}: Mean squared error (default for numeric)
#'         \item \code{"mae"}: Mean absolute error
#'         \item \code{"deviance"}: Gaussian deviance
#'       }
#'   }
#'   Note: "auc" is only valid for binary classification. For multinomial
#'   classification with >2 classes, "auc" will automatically fall back to "class".
#'   Default: \code{NULL}
#' @param custom_contrasts Named list of custom group comparisons for multinomial
#'   models. Each element should be a list with 'group1' and 'group2' specifying
#'   the groups to compare. For example:
#'   \code{list(severe_vs_mild = list(group1 = "Severe", group2 = "Mild"))}.
#'   Only applies to multinomial classification. Default: \code{NULL}
#'
#' @return A list containing regression results with the following structure:
#' \describe{
#'   \item{\code{FunctionOrigin}}{Character string identifying the source function}
#'   \item{\code{ModelSummary}}{Data frame summarizing model performance metrics for all tested alpha values}
#'   \item{\code{DataSplit}}{List containing training/testing data split information}
#'   \item{\code{DataSource}}{Character string indicating which data matrix was used}
#'   \item{\code{ResponseType}}{Character string indicating "categorical" or "numeric"}
#'   \item{\code{SampleDistribution}}{Data frame showing the distribution of samples across
#'     training and testing sets. For categorical responses, shows counts for each group.
#'     Includes a "Total" row summing all samples.}
#'   \item{\code{BestModel}}{List containing results for the best-performing alpha value:
#'     \itemize{
#'       \item{\code{Model}}: The fitted cv.glmnet object
#'       \item{\code{Performance}}: Data frame with performance metrics
#'       \item{\code{Coefficients}}: Data frame of non-zero coefficients (including Intercepts) with Feature names,
#'         Coefficient values, and Odds Ratios. For multinomial models, includes a
#'         Comparison column indicating which groups are being compared (e.g., "Group1 vs Reference")
#'       \item{\code{N_Coefficients}}: Integer count of total non-zero coefficients (including Intercepts)
#'       \item{\code{ConfusionMatrix}}: Confusion matrix (categorical only)
#'       \item{\code{Predictions}}: Data frame with actual vs predicted values
#'       \item{\code{Lambda}}: Selected lambda value
#'       \item{\code{Alpha}}: Alpha value used
#'     }
#'   }
#'   \item{\code{AllAlphaResults}}{Named list containing results for all tested alpha values.
#'     Each element is named "alpha_X" where X is the alpha value (e.g., "alpha_0.5", "alpha_1").
#'     Structure is identical to BestModel.}
#'   \item{\code{AlphaComparison}}{Data frame comparing all tested alpha values, including
#'     performance metrics and the number of non-zero coefficients for each alpha}
#'   \item{\code{ErrorLog}}{List of any warnings or errors encountered}
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Example 1: Single alpha value (Elastic Net with alpha = 0.5)
#' results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Group",
#'   alpha = 0.5,
#'   train_percent = 80,
#'   lambda = "1se",
#'   remember = 123
#' )
#'
#' # Example 2: LASSO regression (alpha = 1)
#' lasso_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Group",
#'   alpha = 1,
#'   lambda = "1se"
#' )
#'
#' # Example 3: Ridge regression (alpha = 0)
#' ridge_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Group",
#'   alpha = 0,
#'   lambda = "1se"
#' )
#'
#' # Example 4: Multiple alpha values for comprehensive tuning
#' multi_alpha_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Group",
#'   alpha = c(0, 0.25, 0.5, 0.75, 1),  # Tests Ridge to LASSO
#'   train_percent = 75,
#'   lambda = "1se",
#'   remember = 123
#' )
#'
#' # Example 5: Numeric response variable
#' numeric_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Response",
#'   alpha = 0.5,
#'   train_percent = 80
#' )
#'
#' # Example 6: Custom contrasts for multinomial classification
#' contrast_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Group",
#'   alpha = 0.5,
#'   custom_contrasts = list(
#'     severe_vs_mild = list(group1 = "Severe", group2 = "Mild"),
#'     moderate_vs_mild = list(group1 = "Moderate", group2 = "Mild")
#'   )
#' )
#'
#' # View results
#' print(results$ModelSummary)           # Performance summary
#' print(results$AlphaComparison)        # Compare all alpha values
#' print(results$BestModel$Coefficients) # Best model coefficients
#'
#' # Extract specific alpha results
#' alpha_05 <- extract_RegressionResults(results, alpha = 0.5)
#' lasso <- extract_RegressionResults(results, alpha = 1)
#' }
#'
#' # Example 7: Binary classification optimizing AUC instead of accuracy
#' auc_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Group",  # Must have exactly 2 levels
#'   alpha = 0.5,
#'   type_measure = "auc"
#' )
#'
#' # Example 8: Numeric response using MAE instead of MSE
#' mae_results <- perform_Regression(
#'   data = preprocessed_data,
#'   specify_response = "Response",  # Numeric variable
#'   alpha = 0.5,
#'   type_measure = "mae"
#' )
#'
#' # View results
#' print(results$ModelSummary) # Performance summary
#'
#' @seealso \code{\link[glmnet]{cv.glmnet}}, \code{\link[caret]{confusionMatrix}},
#'   \code{\link{extract_RegressionResults}}
#'
#' @note
#' \strong{Alpha Parameter Guide:}
#' \itemize{
#'   \item{\strong{alpha = 0}: Ridge regression - keeps all features but shrinks coefficients}
#'   \item{\strong{alpha = 1}: LASSO regression - performs feature selection by setting some coefficients to zero}
#'   \item{\strong{0 < alpha < 1}: Elastic Net - balances feature selection and coefficient shrinkage}
#' }
#'
#' \strong{Why no p-values?} P-values and standard confidence intervals are not
#' provided because they are statistically invalid after variable selection in
#' regularized regression. The selection process induces bias that standard
#' inference methods cannot account for. Valid alternatives include:
#' \itemize{
#'   \item{Cross-validated prediction performance (provided in Performance metrics)}
#'   \item{Coefficient magnitudes and directions (provided in Coefficients)}
#'   \item{Replication in independent datasets}
#'   \item{Specialized post-selection inference methods (e.g., selective inference)}
#' }
#' For scientific reporting, focus on cross-validated performance metrics
#' (Accuracy, RMSE, etc.) rather than individual coefficient p-values.
#'
#' @references
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for
#' Generalized Linear Models via Coordinate Descent, Journal of Statistical
#' Software, Vol. 33(1), 1-22, doi:10.18637/jss.v033.i01.
#'
#' Zou, H. and Hastie, T. (2005) Regularization and variable selection via the
#' elastic net, Journal of the Royal Statistical Society: Series B, 67(2), 301-320.
#'
#' Tibshirani, R. (1996) Regression shrinkage and selection via the lasso,
#' Journal of the Royal Statistical Society: Series B, 58(1), 267-288.
#'
#' Kuhn, M. (2008) Building predictive models in R using the caret package,
#' Journal of Statistical Software, doi:10.18637/jss.v028.i05.
#'
#' Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011) Regularization
#' Paths for Cox's Proportional Hazards Model via Coordinate Descent,
#' Journal of Statistical Software, Vol. 39(5), 1-13, doi:10.18637/jss.v039.i05.
#'
#' @export
# perform_Regression <- function(
    #     data,
#     specify_response = NULL,
#     train_percent    = 80,
#     ref              = NULL,
#     lambda           = "1se",
#     alpha            = 0.5,
#     remember         = 123,
#     verbose          = TRUE,
#     cv_folds         = 10,
#     parallel         = FALSE,
#     standardize      = TRUE,
#     maxit            = 100000,
#     type_measure     = NULL,
#     custom_contrasts = NULL
# ) {
#
#   error_log <- list()
#   old_debug <- getOption("error")
#   options(error = NULL)
#   on.exit(options(error = old_debug), add = TRUE)
#
#   tryCatch({
#
#     # --- 1. Validation and Setup ---
#     validation_result <- .validate_inputs_regression(data, train_percent,
#                                                      lambda, alpha, remember,
#                                                      cv_folds, specify_response,
#                                                      parallel, type_measure)
#     if (!validation_result$valid) stop(validation_result$message)
#
#     if (!is.null(remember)) {
#       set.seed(remember)
#       if (verbose) cat("Random seed set to:", remember, "\n")
#     }
#
#     # Prepare Data
#     # data_prep_result <- .prepare_data_regression(data, specify_response, ref,
#     #                                              train_percent, verbose, cv_folds)
#     data_prep_result <- .prepare_data_regression(data, specify_response, ref,
#                                                  train_percent, verbose, cv_folds, remember)
#     if (!data_prep_result$success) {
#       stop("Data preparation failed: ", data_prep_result$error)
#     }
#     data_prep <- data_prep_result$data
#
#     # Use adjusted cv_folds if it was modified during data prep
#     if ("cv_folds" %in% names(data_prep)) {
#       cv_folds <- data_prep$cv_folds
#     }
#
#     # Create sample distribution summary
#     sample_distribution <- .create_sample_distribution(
#       y_train       = data_prep$y_train,
#       y_test        = data_prep$y_test,
#       response_type = data_prep$response_type,
#       train_percent = train_percent
#     )
#
#     # Initialize Results Container
#     results <- list(
#       FunctionOrigin     = "perform_Regression",
#       ModelSummary       = data.frame(),
#       DataSplit          = data_prep$split_info,
#       DataSource         = data_prep$data_source,
#       ResponseType       = data_prep$response_type,
#       SampleDistribution = sample_distribution,
#       AlphaComparison    = NULL,
#       ErrorLog           = error_log
#     )
#
#     performance_comparison <- list()
#     alpha_performance      <- data.frame()
#
#     # Handle multiple alphas
#     alpha_values <- unique(alpha)
#
#     # Initialize storage for all alpha results
#     results_list <- list()
#
#     # --- 2. Main Loop Over Alpha Values ---
#     for (alpha_val in alpha_values) {
#
#       if (verbose && length(alpha_values) > 1) {
#         cat("\n", paste(rep("=", 60), collapse = ""), "\n")
#         cat("TESTING ALPHA =", alpha_val, "\n")
#         cat(paste(rep("=", 60), collapse = ""), "\n")
#       }
#
#       method_results <- .fit_regression_model(
#         x_train          = data_prep$x_train,
#         y_train          = data_prep$y_train,
#         x_test           = data_prep$x_test,
#         y_test           = data_prep$y_test,
#         alpha_val        = alpha_val,
#         lambda           = lambda,
#         response_info    = data_prep$response_info,
#         cv_folds         = cv_folds,
#         parallel         = parallel,
#         standardize      = standardize,
#         maxit            = maxit,
#         type_measure     = type_measure,
#         verbose          = verbose,
#         custom_contrasts = custom_contrasts
#       )
#
#       # Store results
#       results_list[[paste0("alpha_", alpha_val)]] <- method_results
#
#       # Add to alpha comparison
#       alpha_row <- data.frame(
#         Alpha            = alpha_val,
#         method_results$Performance,
#         N_Coefficients   = method_results$N_Coefficients,
#         stringsAsFactors = FALSE
#       )
#       alpha_performance <- rbind(alpha_performance, alpha_row)
#     }
#
#     # Determine best performing alpha
#     if (data_prep$response_type == "categorical") {
#       best_alpha_idx <- which.min(sapply(results_list,
#                                          function(x) min(x$Model$cvm)
#       ))
#     } else {
#       best_alpha_idx <- which.min(sapply(results_list, function(x) x$Performance$RMSE))
#     }
#
#     # Store best result as main result
#     results$BestModel       <- results_list[[best_alpha_idx]]
#     results$AllAlphaResults <- results_list
#
#     # Create alpha comparison summary
#     if (nrow(alpha_performance) >= 1) {
#       results$AlphaComparison           <- alpha_performance[order(alpha_performance$Alpha), ]
#       rownames(results$AlphaComparison) <- NULL
#     }
#
#     # Create model summary
#     results$ModelSummary <- .create_model_summary(alpha_performance, data_prep$response_type)
#
#     if (verbose) {
#       cat("\n", paste(rep("=", 60), collapse = ""), "\n")
#       cat("REGRESSION ANALYSIS COMPLETE\n")
#       cat(paste(rep("=", 60), collapse = ""), "\n")
#       print(results$ModelSummary)
#       if (!is.null(results$AlphaComparison)) {
#         cat("\nALPHA COMPARISON:\n")
#         print(results$AlphaComparison)
#       }
#
#       if (!is.null(results$SampleDistribution)) {
#         cat("\nSAMPLE DISTRIBUTION:\n")
#         print(results$SampleDistribution, row.names = FALSE)
#       }
#     }
#
#     class(results) <- c("perform_Regression", "list")
#     return(results)
#
#   }, error = function(e) {
#     error_msg <- paste("Error in perform_Regression:", e$message)
#     if (verbose) {
#       cat("\n", paste(rep("!", 60), collapse = ""), "\n")
#       cat("ERROR OCCURRED:\n")
#       cat(error_msg, "\n")
#       cat(paste(rep("!", 60), collapse = ""), "\n")
#     }
#
#     return(list(
#       FunctionOrigin = "perform_Regression",
#       Error          = error_msg,
#       ErrorDetails   = as.character(e),
#       ErrorLog       = error_log,
#       Success        = FALSE
#     ))
#   })
# }
perform_Regression <- function(
    data,
    specify_response = NULL,
    train_percent    = 80,
    ref              = NULL,
    lambda           = "1se",
    alpha            = 0.5,
    remember         = 123,
    verbose          = TRUE,
    cv_folds         = 10,
    parallel         = FALSE,
    standardize      = TRUE,
    maxit            = 100000,
    type_measure     = NULL,
    custom_contrasts = NULL
) {

  error_log <- list()
  old_debug <- getOption("error")
  options(error = NULL)
  on.exit(options(error = old_debug), add = TRUE)

  tryCatch({

    # --- 1. Validation and Setup ---
    validation_result <- .validate_inputs_regression(data, train_percent,
                                                     lambda, alpha, remember,
                                                     cv_folds, specify_response,
                                                     parallel, type_measure)
    if (!validation_result$valid) stop(validation_result$message)

    if (!is.null(remember)) {
      set.seed(remember)
      if (verbose) cat("Random seed set to:", remember, "\n")
    }

    # Prepare Data
    data_prep_result <- .prepare_data_regression(data, specify_response, ref,
                                                 train_percent, verbose, cv_folds, remember)
    if (!data_prep_result$success) {
      stop("Data preparation failed: ", data_prep_result$error)
    }
    data_prep <- data_prep_result$data

    # Use adjusted cv_folds if it was modified during data prep
    if ("cv_folds" %in% names(data_prep)) {
      cv_folds <- data_prep$cv_folds
    }

    # Create sample distribution summary
    sample_distribution <- .create_sample_distribution(
      y_train       = data_prep$y_train,
      y_test        = data_prep$y_test,
      response_type = data_prep$response_type,
      train_percent = train_percent
    )

    # Initialize Results Container
    results <- list(
      FunctionOrigin     = "perform_Regression",
      ModelSummary       = data.frame(),
      DataSplit          = data_prep$split_info,
      DataSource         = data_prep$data_source,
      ResponseType       = data_prep$response_type,
      SampleDistribution = sample_distribution,
      AlphaComparison    = NULL,
      ErrorLog           = error_log
    )

    performance_comparison <- list()
    alpha_performance      <- data.frame()

    # Handle multiple alphas
    alpha_values <- unique(alpha)

    # Initialize storage for all alpha results
    results_list <- list()

    # --- FIX: Generate Fold IDs Once ---
    # This ensures that every alpha value uses the EXACT same cross-validation folds.
    # This guarantees consistency regardless of the number/order of alphas tested.
    if (!is.null(remember)) set.seed(remember)

    if (data_prep$response_type == "categorical") {
      # Stratified folds for categorical data (using caret to ensure balance)
      # list = FALSE returns a vector of integers (1 to k)
      fold_ids <- caret::createFolds(data_prep$y_train, k = cv_folds, list = FALSE)
    } else {
      # Random balanced folds for numeric data
      n_train  <- nrow(data_prep$x_train)
      fold_ids <- sample(rep(seq_len(cv_folds), length.out = n_train))
    }
    # -----------------------------------

    # --- 2. Main Loop Over Alpha Values ---
    for (alpha_val in alpha_values) {

      if (verbose && length(alpha_values) > 1) {
        cat("\n", paste(rep("=", 60), collapse = ""), "\n")
        cat("TESTING ALPHA =", alpha_val, "\n")
        cat(paste(rep("=", 60), collapse = ""), "\n")
      }

      method_results <- .fit_regression_model(
        x_train          = data_prep$x_train,
        y_train          = data_prep$y_train,
        x_test           = data_prep$x_test,
        y_test           = data_prep$y_test,
        alpha_val        = alpha_val,
        lambda           = lambda,
        response_info    = data_prep$response_info,
        cv_folds         = cv_folds,
        fold_ids         = fold_ids,  # <--- PASSING THE FIXED FOLD IDs
        parallel         = parallel,
        standardize      = standardize,
        maxit            = maxit,
        type_measure     = type_measure,
        verbose          = verbose,
        custom_contrasts = custom_contrasts
      )

      # Store results
      results_list[[paste0("alpha_", alpha_val)]] <- method_results

      # Add to alpha comparison
      alpha_row <- data.frame(
        Alpha            = alpha_val,
        method_results$Performance,
        N_Coefficients   = method_results$N_Coefficients,
        stringsAsFactors = FALSE
      )
      alpha_performance <- rbind(alpha_performance, alpha_row)
    }

    # Determine best performing alpha
    if (data_prep$response_type == "categorical") {
      best_alpha_idx <- which.min(sapply(results_list,
                                         function(x) min(x$Model$cvm)
      ))
    } else {
      best_alpha_idx <- which.min(sapply(results_list, function(x) x$Performance$RMSE))
    }

    # Store best result as main result
    results$BestModel       <- results_list[[best_alpha_idx]]
    results$AllAlphaResults <- results_list

    # Create alpha comparison summary
    if (nrow(alpha_performance) >= 1) {
      results$AlphaComparison           <- alpha_performance[order(alpha_performance$Alpha), ]
      rownames(results$AlphaComparison) <- NULL
    }

    # Create model summary
    results$ModelSummary <- .create_model_summary(alpha_performance, data_prep$response_type)

    if (verbose) {
      cat("\n", paste(rep("=", 60), collapse = ""), "\n")
      cat("REGRESSION ANALYSIS COMPLETE\n")
      cat(paste(rep("=", 60), collapse = ""), "\n")
      print(results$ModelSummary)
      # if (!is.null(results$AlphaComparison)) {
      #   cat("\nALPHA COMPARISON:\n")
      #   print(results$AlphaComparison)
      # }

      if (!is.null(results$SampleDistribution)) {
        cat("\nSAMPLE DISTRIBUTION:\n")
        print(results$SampleDistribution, row.names = FALSE)
      }
    }

    class(results) <- c("perform_Regression", "list")
    return(results)

  }, error = function(e) {
    error_msg <- paste("Error in perform_Regression:", e$message)
    if (verbose) {
      cat("\n", paste(rep("!", 60), collapse = ""), "\n")
      cat("ERROR OCCURRED:\n")
      cat(error_msg, "\n")
      cat(paste(rep("!", 60), collapse = ""), "\n")
    }

    return(list(
      FunctionOrigin = "perform_Regression",
      Error          = error_msg,
      ErrorDetails   = as.character(e),
      ErrorLog       = error_log,
      Success        = FALSE
    ))
  })
}

# Input validation
.validate_inputs_regression <- function(data, train_percent, lambda, alpha, remember, cv_folds, specify_response, parallel, type_measure = NULL) {
  # Check data structure
  if (!is.list(data)) {
    return(list(valid = FALSE, message = "'data' must be a list object from perform_PreprocessingPeakData function."))
  }

  if (!"FunctionOrigin" %in% names(data) || data$FunctionOrigin != "perform_PreprocessingPeakData") {
    return(list(valid = FALSE, message = "The supplied data did not come from 'perform_PreprocessingPeakData' function."))
  }

  required_elements <- c("Metadata")
  missing_elements  <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    return(list(valid = FALSE, message = paste("Missing required data elements:", paste(missing_elements, collapse = ", "))))
  }

  # Check for at least one data matrix
  if (!("data_scaledNONPLS_varFiltered" %in% names(data)) && !("data_scaledNONPLS_merged" %in% names(data))) {
    return(list(valid = FALSE, message = "No valid data matrix found. Need either 'data_scaledNONPLS_varFiltered' or 'data_scaledNONPLS_merged'."))
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

  # Check type_measure
  if (!is.null(type_measure)) {
    valid_measures <- c("class", "auc", "deviance", "mse", "mae")
    if (!type_measure %in% valid_measures) {
      return(list(valid = FALSE,
                  message = paste("'type_measure' must be one of:",
                                  paste(valid_measures, collapse = ", "),
                                  "or NULL for automatic selection.")))
    }
  }

  # Check Parallel Backend
  if (isTRUE(parallel)) {
    # Check if a parallel backend is registered (e.g., doParallel/doSNOW)
    if (!foreach::getDoParRegistered()) {
      warning("Parallel processing set to TRUE, but no parallel backend is registered. ",
              "Cross-validation will run sequentially. ",
              "Use doParallel::registerDoParallel() before running this function.")
    }
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

# Data preparation function
# .prepare_data_regression <- function(data, specify_response, ref, train_percent, verbose, cv_folds) {
#
#   tryCatch({
#     if ("Metadata_merged" %in% names(data)) {
#       metadata <- data$Metadata_merged
#     } else {
#       warning("Metadata_merged not found, using Metadata instead")
#       metadata <- data$Metadata
#     }
#
#     # Determine which data matrix to use
#     data_source <- "data_scaledNONPLS_varFiltered"  # default
#
#     if ("Parameters" %in% names(data) &&
#         "merge_replicates" %in% names(data$Parameters) &&
#         data$Parameters$merge_replicates == TRUE &&
#         "data_scaledNONPLS_merged" %in% names(data)) {
#
#       if (verbose) cat("Using merged replicate data: data_scaledNONPLS_merged\n")
#       features <- data$data_scaledNONPLS_merged
#       data_source <- "data_scaledNONPLS_merged"
#
#       # Extract metadata
#       metadata <- data$Metadata_merged
#
#     } else {
#
#       if (verbose) cat("Using filtered data: data_scaledNONPLS_varFiltered\n")
#       features <- data$data_scaledNONPLS_varFiltered
#       data_source <- "data_scaledNONPLS_varFiltered"
#
#       # Extract metadata
#       metadata <- data$Metadata
#     }
#
#     # Convert features to matrix if not already
#     if (!is.matrix(features)) {
#       features <- as.matrix(features)
#     }
#
#     # Identify QC samples
#     qc_indices <- metadata$Group %in% c("SQC", "EQC", "QC")
#     non_qc_indices <- !qc_indices
#
#     if (sum(non_qc_indices) == 0) {
#       return(list(success = FALSE, error = "No non-QC samples found in the data."))
#     }
#
#     # Prepare response variable
#     if (!is.null(specify_response)) {
#       if (!specify_response %in% colnames(metadata)) {
#         return(list(success = FALSE, error = paste("Specified response variable '", specify_response, "' not found in metadata.")))
#       }
#       response_var <- metadata[[specify_response]][non_qc_indices]
#     } else {
#       response_var <- metadata$Group[non_qc_indices]
#     }
#
#     # Remove NA values
#     na_indices <- is.na(response_var)
#     if (any(na_indices)) {
#       if (verbose) cat("Removing", sum(na_indices), "samples with NA response values.\n")
#       non_qc_indices[which(non_qc_indices)[na_indices]] <- FALSE
#       response_var <- response_var[!na_indices]
#     }
#
#     # Determine response type and prepare accordingly
#     response_type <- "categorical"
#     response_info <- list()
#
#     if (is.numeric(response_var)) {
#       response_type <- "numeric"
#       response_info$type <- "numeric"
#       response_info$min <- min(response_var, na.rm = TRUE)
#       response_info$max <- max(response_var, na.rm = TRUE)
#       response_info$mean <- mean(response_var, na.rm = TRUE)
#       response_info$sd <- sd(response_var, na.rm = TRUE)
#
#       if (verbose) {
#         cat("Response variable is numeric:\n")
#         cat("  Range:", round(response_info$min, 4), "to", round(response_info$max, 4), "\n")
#         cat("  Mean (SD):", round(response_info$mean, 4), "(", round(response_info$sd, 4), ")\n")
#       }
#
#     } else {
#       # Categorical response
#       response_type <- "categorical"
#
#       # Set reference level
#       ref_level <- NULL
#       if (!is.null(ref)) {
#         if (!ref %in% unique(response_var)) {
#           return(list(success = FALSE, error = paste("Reference level '", ref, "' not found in response variable.")))
#         }
#         ref_level <- ref
#         response_var <- factor(response_var, levels = c(ref, setdiff(unique(response_var), ref)))
#       } else {
#         response_var <- factor(response_var)
#         ref_level <- levels(response_var)[1]
#       }
#
#       # Check for sufficient samples per group
#       group_counts <- table(response_var)
#       if (any(group_counts < 3)) {
#         if (verbose) cat("Warning: Some groups have fewer than 3 samples. This may affect model performance.\n")
#       }
#
#       response_info$type <- "categorical"
#       response_info$levels <- levels(response_var)
#       response_info$ref_level <- ref_level
#       response_info$group_counts <- as.vector(group_counts)
#       names(response_info$group_counts) <- names(group_counts)
#
#       if (verbose) {
#         cat("Response variable is categorical:\n")
#         cat("  Groups:", paste(levels(response_var), collapse = ", "), "\n")
#         cat("  Reference level:", ref_level, "\n")
#         cat("  Group distribution:\n")
#         print(group_counts)
#       }
#     }
#
#     # Extract non-QC feature data
#     df <- features[non_qc_indices, , drop = FALSE]
#     df <- df[!na_indices, , drop = FALSE]
#
#     # Ensure df is a matrix
#     if (!is.matrix(df)) {
#       df <- as.matrix(df)
#     }
#
#     # Check for constant or near-constant features
#     feature_var <- apply(df, 2, var, na.rm = TRUE)
#     constant_features <- is.na(feature_var) | feature_var < 1e-10
#     if (any(constant_features)) {
#       if (verbose) cat("Removing", sum(constant_features), "constant or near-constant features.\n")
#       df <- df[, !constant_features, drop = FALSE]
#     }
#
#     # Check for missing values in features
#     if (any(is.na(df))) {
#       if (verbose) cat("Warning: Missing values detected in features. Removing incomplete cases.\n")
#       complete_cases <- complete.cases(df)
#       if (sum(complete_cases) < nrow(df)) {
#         if (verbose) cat("Removing", sum(!complete_cases), "samples with missing feature values.\n")
#         df <- df[complete_cases, , drop = FALSE]
#         response_var <- response_var[complete_cases]
#       }
#     }
#
#     # Split data
#     n_samples <- nrow(df)
#     train_size <- floor((train_percent / 100) * n_samples)
#
#     if (train_size < 3) {
#       return(list(success = FALSE, error = "Training set too small. Need at least 3 samples."))
#     }
#
#     # Split data using caret for robust stratification
#     if (response_type == "categorical") {
#       train_indices <- caret::createDataPartition(
#         y = response_var,
#         p = train_percent / 100,
#         list = FALSE,
#         times = 1
#       )
#       train_indices <- as.vector(train_indices)
#     } else {
#       # Simple random sampling for numeric
#       train_indices <- sample(seq_len(n_samples), train_size)
#     }
#
#     x_train <- df[train_indices, , drop = FALSE]
#     x_test <- df[-train_indices, , drop = FALSE]
#     y_train <- response_var[train_indices]
#     y_test <- response_var[-train_indices]
#
#     # Validate cv_folds for categorical responses
#     if (response_type == "categorical") {
#       min_group_size <- min(table(y_train))
#       if (cv_folds > min_group_size) {
#         cv_folds_original <- cv_folds
#         cv_folds <- max(3, min_group_size)  # Ensure at least 3 folds
#         if (verbose) {
#           cat("Warning: Reducing cv_folds from", cv_folds_original, "to", cv_folds,
#               "due to small group size (min =", min_group_size, "samples)\n")
#         }
#       }
#     }
#
#     # Final validation
#     if (nrow(x_train) == 0 || nrow(x_test) == 0) {
#       return(list(success = FALSE, error = "Data split resulted in empty training or testing set."))
#     }
#
#     if (ncol(x_train) == 0) {
#       return(list(success = FALSE, error = "No features remaining after preprocessing."))
#     }
#
#     # Ensure factor levels are consistent for categorical
#     if (response_type == "categorical") {
#       y_train <- factor(y_train, levels = levels(response_var))
#       y_test <- factor(y_test, levels = levels(response_var))
#     }
#
#     if (verbose) {
#       cat("Data preparation complete:\n")
#       cat("  Data source:", data_source, "\n")
#       cat("  Response type:", response_type, "\n")
#       cat("  Total samples:", n_samples, "\n")
#       cat("  Training samples:", nrow(x_train), "\n")
#       cat("  Testing samples:", nrow(x_test), "\n")
#       cat("  Features:", ncol(x_train), "\n")
#
#       if (response_type == "categorical") {
#         cat("  Training group distribution:\n")
#         print(table(y_train))
#         cat("  Testing group distribution:\n")
#         print(table(y_test))
#       }
#     }
#
#     return(list(
#       success = TRUE,
#       data = list(
#         x_train = x_train,
#         x_test = x_test,
#         y_train = y_train,
#         y_test = y_test,
#         response_info = response_info,
#         response_type = response_type,
#         data_source = data_source,
#         cv_folds = cv_folds,
#         split_info = list(
#           total_samples = n_samples,
#           train_samples = nrow(x_train),
#           test_samples = nrow(x_test),
#           features = ncol(x_train),
#           response_type = response_type,
#           train_indices = train_indices
#         )
#       )
#     ))
#
#   }, error = function(e) {
#     return(list(success = FALSE, error = paste("Data preparation error:", e$message)))
#   })
# }

.prepare_data_regression <- function(data, specify_response, ref, train_percent, verbose, cv_folds, remember = NULL) {

  tryCatch({
    if ("Metadata_merged" %in% names(data)) {
      metadata <- data$Metadata_merged
    } else {
      warning("Metadata_merged not found, using Metadata instead")
      metadata <- data$Metadata
    }

    # Determine which data matrix to use
    data_source <- "data_scaledNONPLS_varFiltered"  # default

    if ("Parameters" %in% names(data) &&
        "merge_replicates" %in% names(data$Parameters) &&
        data$Parameters$merge_replicates == TRUE &&
        "data_scaledNONPLS_merged" %in% names(data)) {

      if (verbose) cat("Using merged replicate data: data_scaledNONPLS_merged\n")
      features    <- data$data_scaledNONPLS_merged
      data_source <- "data_scaledNONPLS_merged"

      # Extract metadata
      metadata <- data$Metadata_merged

    } else {

      if (verbose) cat("Using filtered data: data_scaledNONPLS_varFiltered\n")
      features    <- data$data_scaledNONPLS_varFiltered
      data_source <- "data_scaledNONPLS_varFiltered"

      # Extract metadata
      metadata <- data$Metadata
    }

    # Convert features to matrix if not already
    if (!is.matrix(features)) {
      features <- as.matrix(features)
    }

    # Identify QC samples
    qc_indices     <- metadata$Group %in% c("SQC", "EQC", "QC")
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
      response_var                                      <- response_var[!na_indices]
    }

    # Determine response type and prepare accordingly
    response_type <- "categorical"
    response_info <- list()

    if (is.numeric(response_var)) {
      response_type      <- "numeric"
      response_info$type <- "numeric"
      response_info$min  <- min(response_var, na.rm = TRUE)
      response_info$max  <- max(response_var, na.rm = TRUE)
      response_info$mean <- mean(response_var, na.rm = TRUE)
      response_info$sd   <- sd(response_var, na.rm = TRUE)

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
        ref_level    <- levels(response_var)[1]
      }

      # Check for sufficient samples per group
      group_counts <- table(response_var)
      if (any(group_counts < 3)) {
        if (verbose) cat("Warning: Some groups have fewer than 3 samples. This may affect model performance.\n")
      }

      response_info$type                <- "categorical"
      response_info$levels              <- levels(response_var)
      response_info$ref_level           <- ref_level
      response_info$group_counts        <- as.vector(group_counts)
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
    feature_var       <- apply(df, 2, var, na.rm = TRUE)
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
        df           <- df[complete_cases, , drop = FALSE]
        response_var <- response_var[complete_cases]
      }
    }

    # Split data
    n_samples  <- nrow(df)
    train_size <- floor((train_percent / 100) * n_samples)

    if (train_size < 3) {
      return(list(success = FALSE, error = "Training set too small. Need at least 3 samples."))
    }

    # === FIX: Set seed right before splitting ===
    if (!is.null(remember)) {
      set.seed(remember)
    }

    # Split data using caret for robust stratification
    if (response_type == "categorical") {
      train_indices <- caret::createDataPartition(
        y     = response_var,
        p     = train_percent / 100,
        list  = FALSE,
        times = 1
      )
      train_indices <- as.vector(train_indices)
    } else {
      # Simple random sampling for numeric
      train_indices <- sample(seq_len(n_samples), train_size)
    }

    x_train <- df[train_indices, , drop = FALSE]
    x_test  <- df[-train_indices, , drop = FALSE]
    y_train <- response_var[train_indices]
    y_test  <- response_var[-train_indices]

    # Validate cv_folds for categorical responses
    if (response_type == "categorical") {
      min_group_size <- min(table(y_train))
      if (cv_folds > min_group_size) {
        cv_folds_original <- cv_folds
        cv_folds          <- max(3, min_group_size)  # Ensure at least 3 folds
        if (verbose) {
          cat("Warning: Reducing cv_folds from", cv_folds_original, "to", cv_folds,
              "due to small group size (min =", min_group_size, "samples)\n")
        }
      }
    }

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
      y_test  <- factor(y_test, levels = levels(response_var))
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
      data    = list(
        x_train       = x_train,
        x_test        = x_test,
        y_train       = y_train,
        y_test        = y_test,
        response_info = response_info,
        response_type = response_type,
        data_source   = data_source,
        cv_folds      = cv_folds,
        split_info = list(
          total_samples = n_samples,
          train_samples = nrow(x_train),
          test_samples  = nrow(x_test),
          features      = ncol(x_train),
          response_type = response_type,
          train_indices = train_indices
        )
      )
    ))

  }, error = function(e) {
    return(list(success = FALSE, error = paste("Data preparation error:", e$message)))
  })
}

# Model fitting function
# .fit_regression_model <- function(x_train, y_train, x_test, y_test,
#                                   alpha_val, lambda, response_info, cv_folds,
#                                   parallel, standardize, maxit, type_measure = NULL,
#                                   verbose, custom_contrasts = NULL) {
#
#   # # Determine family
#   # if (response_info$type == "categorical") {
#   #   family_type <- if (length(response_info$levels) == 2) "binomial" else "multinomial"
#   #   type_measure <- "class"
#   # } else {
#   #   family_type <- "gaussian"
#   #   type_measure <- "mse"
#   # }
#   # Determine family and type.measure
#   if (response_info$type == "categorical") {
#     family_type <- if (length(response_info$levels) == 2) "binomial" else "multinomial"
#
#     # Set default or validate type_measure for categorical
#     if (is.null(type_measure)) {
#       type_measure <- "class"  # default for categorical
#     } else {
#       # Validate type_measure for categorical responses
#       valid_categorical <- c("class", "auc", "deviance")
#       if (!type_measure %in% valid_categorical) {
#         warning("type_measure '", type_measure, "' is not valid for categorical responses. ",
#                 "Using 'class' instead. Valid options: ", paste(valid_categorical, collapse = ", "))
#         type_measure <- "class"
#       }
#
#       # Special check: auc only works for binomial
#       if (type_measure == "auc" && family_type == "multinomial") {
#         warning("type_measure 'auc' is only valid for binary classification. ",
#                 "Using 'class' for multinomial response.")
#         type_measure <- "class"
#       }
#     }
#
#   } else {
#     # Numeric response
#     family_type <- "gaussian"
#
#     # Set default or validate type_measure for numeric
#     if (is.null(type_measure)) {
#       type_measure <- "mse"  # default for numeric
#     } else {
#       # Validate type_measure for numeric responses
#       valid_numeric <- c("mse", "mae", "deviance")
#       if (!type_measure %in% valid_numeric) {
#         warning("type_measure '", type_measure, "' is not valid for numeric responses. ",
#                 "Using 'mse' instead. Valid options: ", paste(valid_numeric, collapse = ", "))
#         type_measure <- "mse"
#       }
#     }
#   }
#
#   if (verbose) cat("  Using type.measure:", type_measure, "\n")
#
#   # Fit glmnet
#   model_fit <- glmnet::cv.glmnet(
#     x            = x_train,
#     y            = y_train,
#     alpha        = alpha_val,
#     family       = family_type,
#     type.measure = type_measure,
#     nfolds       = cv_folds,
#     parallel     = parallel,
#     standardize  = standardize,
#     maxit        = maxit
#   )
#
#   # selected_lambda <- if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min
#   selected_lambda <- if (lambda == "1se") {
#     if (is.na(model_fit$lambda.1se)) model_fit$lambda.min else model_fit$lambda.1se
#   } else {
#     model_fit$lambda.min
#     message("lambda == 1se failed (lambda is NA, used min instead")
#   }
#
#   # --- Prediction & Performance Logic ---
#
#   if (response_info$type == "categorical") {
#     # Class predictions
#     preds_class <- predict(model_fit, s = selected_lambda, newx = x_test, type = "class")
#     preds_class <- factor(as.vector(preds_class), levels = response_info$levels)
#
#     # Store Predictions
#     prediction_df <- data.frame(
#       Actual           = y_test,
#       Predicted        = preds_class,
#       stringsAsFactors = FALSE
#     )
#
#     conf_matrix <- caret::confusionMatrix(preds_class, factor(y_test, levels = response_info$levels))
#
#     accuracy    <- conf_matrix$overall["Accuracy"]
#     kappa       <- conf_matrix$overall["Kappa"]
#
#     performance <- data.frame(Accuracy = accuracy, Kappa = kappa)
#
#   } else {
#     # Numeric predictions
#     preds <- predict(model_fit, s = selected_lambda, newx = x_test, type = "response")
#     preds <- as.vector(preds)
#
#     # Store Predictions
#     prediction_df <- data.frame(
#       Actual           = y_test,
#       Predicted        = preds,
#       stringsAsFactors = FALSE
#     )
#
#     # Metrics
#     mse  <- mean((y_test - preds)^2)
#     rmse <- sqrt(mse)
#     mae  <- mean(abs(y_test - preds))
#
#     # Calculate R-squared
#     sst       <- sum((y_test - mean(y_test))^2)
#     sse       <- sum((y_test - preds)^2)
#     r_squared <- 1 - (sse / sst)
#
#     performance <- data.frame(
#       MSE       = mse,
#       RMSE      = rmse,
#       MAE       = mae,
#       R_squared = r_squared
#     )
#
#     conf_matrix <- NULL
#   }
#
#   # Extract Coefficients & Count
#   coef_extract <- .extract_coefficients(model_fit, selected_lambda, response_info, family_type, custom_contrasts)
#   coef_results <- coef_extract$Coefficients
#   n_coefs      <- coef_extract$N_Coefficients
#
#   if (verbose) cat("  Non-zero coefficients (incl. intercept):", n_coefs, "\n")
#
#   return(list(
#     Model           = model_fit,
#     Performance     = performance,
#     ConfusionMatrix = conf_matrix,
#     Predictions     = prediction_df,
#     Coefficients    = coef_results,
#     N_Coefficients  = n_coefs,
#     Lambda          = selected_lambda,
#     Alpha           = alpha_val
#   ))
# }
.fit_regression_model <- function(x_train, y_train, x_test, y_test,
                                  alpha_val, lambda, response_info, cv_folds,
                                  fold_ids, # <--- NEW ARGUMENT
                                  parallel, standardize, maxit, type_measure = NULL,
                                  verbose, custom_contrasts = NULL) {

  # Determine family and type.measure
  if (response_info$type == "categorical") {
    family_type <- if (length(response_info$levels) == 2) "binomial" else "multinomial"

    # Set default or validate type_measure for categorical
    if (is.null(type_measure)) {
      type_measure <- "class"  # default for categorical
    } else {
      # Validate type_measure for categorical responses
      valid_categorical <- c("class", "auc", "deviance")
      if (!type_measure %in% valid_categorical) {
        warning("type_measure '", type_measure, "' is not valid for categorical responses. ",
                "Using 'class' instead. Valid options: ", paste(valid_categorical, collapse = ", "))
        type_measure <- "class"
      }

      # Special check: auc only works for binomial
      if (type_measure == "auc" && family_type == "multinomial") {
        warning("type_measure 'auc' is only valid for binary classification. ",
                "Using 'class' for multinomial response.")
        type_measure <- "class"
      }
    }

  } else {
    # Numeric response
    family_type <- "gaussian"

    # Set default or validate type_measure for numeric
    if (is.null(type_measure)) {
      type_measure <- "mse"  # default for numeric
    } else {
      # Validate type_measure for numeric responses
      valid_numeric <- c("mse", "mae", "deviance")
      if (!type_measure %in% valid_numeric) {
        warning("type_measure '", type_measure, "' is not valid for numeric responses. ",
                "Using 'mse' instead. Valid options: ", paste(valid_numeric, collapse = ", "))
        type_measure <- "mse"
      }
    }
  }

  if (verbose) cat("  Using type.measure:", type_measure, "\n")

  # Fit glmnet with FIXED FOLDS
  model_fit <- glmnet::cv.glmnet(
    x            = x_train,
    y            = y_train,
    alpha        = alpha_val,
    family       = family_type,
    type.measure = type_measure,
    foldid       = fold_ids, # <--- USING FIXED FOLDS HERE
    # nfolds     = cv_folds, # nfolds is ignored if foldid is provided
    parallel     = parallel,
    standardize  = standardize,
    maxit        = maxit
  )

  selected_lambda <- if (lambda == "1se") {
    if (is.na(model_fit$lambda.1se)) model_fit$lambda.min else model_fit$lambda.1se
  } else {
    model_fit$lambda.min
  }

  # If 1se failed and we fell back to min, warn logic could go here, but the check above handles the NA assignment.

  # --- Prediction & Performance Logic ---

  if (response_info$type == "categorical") {
    # Class predictions
    preds_class <- predict(model_fit, s = selected_lambda, newx = x_test, type = "class")
    preds_class <- factor(as.vector(preds_class), levels = response_info$levels)

    # Store Predictions
    prediction_df <- data.frame(
      Actual           = y_test,
      Predicted        = preds_class,
      stringsAsFactors = FALSE
    )

    conf_matrix <- caret::confusionMatrix(preds_class, factor(y_test, levels = response_info$levels))

    accuracy    <- conf_matrix$overall["Accuracy"]
    kappa       <- conf_matrix$overall["Kappa"]

    performance <- data.frame(Accuracy = accuracy, Kappa = kappa)

  } else {
    # Numeric predictions
    preds <- predict(model_fit, s = selected_lambda, newx = x_test, type = "response")
    preds <- as.vector(preds)

    # Store Predictions
    prediction_df <- data.frame(
      Actual           = y_test,
      Predicted        = preds,
      stringsAsFactors = FALSE
    )

    # Metrics
    mse  <- mean((y_test - preds)^2)
    rmse <- sqrt(mse)
    mae  <- mean(abs(y_test - preds))

    # Calculate R-squared
    sst       <- sum((y_test - mean(y_test))^2)
    sse       <- sum((y_test - preds)^2)
    r_squared <- 1 - (sse / sst)

    performance <- data.frame(
      MSE       = mse,
      RMSE      = rmse,
      MAE       = mae,
      R_squared = r_squared
    )

    conf_matrix <- NULL
  }

  # Extract Coefficients & Count
  coef_extract <- .extract_coefficients(model_fit, selected_lambda, response_info, family_type, custom_contrasts)
  coef_results <- coef_extract$Coefficients
  n_coefs      <- coef_extract$N_Coefficients

  if (verbose) cat("  Non-zero coefficients (incl. intercept):", n_coefs, "\n")

  return(list(
    Model           = model_fit,
    Performance     = performance,
    ConfusionMatrix = conf_matrix,
    Predictions     = prediction_df,
    Coefficients    = coef_results,
    N_Coefficients  = n_coefs,
    Lambda          = selected_lambda,
    Alpha           = alpha_val
  ))
}

# Coefficient extraction
.extract_coefficients <- function(model_fit, lambda, response_info,
                                  family_type, custom_contrasts = NULL) {

  coef_model <- glmnet::coef.glmnet(model_fit, s = lambda)

  if (family_type == "multinomial") {
    group_levels <- response_info$levels
    ref_level    <- response_info$ref_level

    # Convert sparse matrices to dense matrix
    coef_matrix           <- do.call(cbind, lapply(coef_model, as.matrix))
    colnames(coef_matrix) <- group_levels

    # CORRECT: Normalize to Reference Level
    ref_col_idx       <- which(colnames(coef_matrix) == ref_level)
    ref_betas         <- coef_matrix[, ref_col_idx]
    normalized_matrix <- coef_matrix - ref_betas

    results_list <- list()

    # Get Coefficients for Groups vs Reference
    target_groups <- setdiff(group_levels, ref_level)

    for (grp in target_groups) {
      betas <- normalized_matrix[, grp]
      # Include (Intercept) if non-zero: Removed the filter `& rownames(normalized_matrix) != "(Intercept)"`
      active_idx <- which(betas != 0)

      if (length(active_idx) > 0) {
        df <- data.frame(
          Comparison       = paste0(grp, " vs ", ref_level),
          Feature          = rownames(normalized_matrix)[active_idx],
          Coefficient      = betas[active_idx],
          OddsRatio        = exp(betas[active_idx]),
          stringsAsFactors = FALSE
        )
        results_list[[grp]] <- df
      }
    }

    # User-specified custom contrasts (optional)
    if (!is.null(custom_contrasts)) {
      for (contrast_name in names(custom_contrasts)) {
        contrast_spec <- custom_contrasts[[contrast_name]]

        # Expect format: list(group1 = "GroupA", group2 = "GroupB")
        if (all(c("group1", "group2") %in% names(contrast_spec))) {
          grp1 <- contrast_spec$group1
          grp2 <- contrast_spec$group2

          # Validate groups exist
          if (grp1 %in% group_levels && grp2 %in% group_levels) {
            contrast_betas <- normalized_matrix[, grp1] - normalized_matrix[, grp2]
            # Include (Intercept) if non-zero: Removed the filter `& rownames(normalized_matrix) != "(Intercept)"`
            active_idx <- which(contrast_betas != 0)

            if (length(active_idx) > 0) {
              df_contrast <- data.frame(
                Comparison       = paste0(grp1, " vs ", grp2),
                Feature          = rownames(normalized_matrix)[active_idx],
                Coefficient      = contrast_betas[active_idx],
                OddsRatio        = exp(contrast_betas[active_idx]),
                stringsAsFactors = FALSE
              )
              results_list[[contrast_name]] <- df_contrast
            }
          }
        }
      }
    }

    final_df <- do.call(rbind, results_list)
    rownames(final_df) <- NULL

    # Count of non-zero coefficients (total rows in the combined comparison df)
    n_coefs <- nrow(final_df)

    return(list(
      Coefficients   = final_df,
      N_Coefficients = n_coefs
    ))


  } else if (family_type == "binomial") {
    # Binary classification
    coef_df           <- as.data.frame(as.matrix(coef_model))
    colnames(coef_df) <- "Coefficient"
    coef_df$Feature   <- rownames(coef_df)
    coef_df           <- coef_df[coef_df$Coefficient != 0, ] # Intercept is retained if non-zero
    coef_df$OddsRatio <- exp(coef_df$Coefficient)
    coef_df           <- coef_df[, c("Feature", "Coefficient", "OddsRatio")]

    n_coefs <- nrow(coef_df)

    return(list(
      Coefficients   = coef_df,
      N_Coefficients = n_coefs
    ))

  } else {
    # Gaussian (regression)
    coef_df           <- as.data.frame(as.matrix(coef_model))
    colnames(coef_df) <- "Coefficient"
    coef_df$Feature   <- rownames(coef_df)
    coef_df           <- coef_df[coef_df$Coefficient != 0, ] # Intercept is retained if non-zero
    coef_df           <- coef_df[, c("Feature", "Coefficient")]

    n_coefs <- nrow(coef_df)

    return(list(
      Coefficients   = coef_df,
      N_Coefficients = n_coefs
    ))
  }
}

# Model summary creation
.create_model_summary <- function(alpha_performance, response_type) {
  # Simply return the alpha performance dataframe with rankings
  if (response_type == "categorical") {
    alpha_performance$AccuracyRank <- rank(-alpha_performance$Accuracy, ties.method = "min")
  } else {
    alpha_performance$RMSERank <- rank(alpha_performance$RMSE, ties.method = "min")
  }
  return(alpha_performance)
}

# Create sample distribution summary
.create_sample_distribution <- function(y_train, y_test, response_type, train_percent) {

  if (response_type == "categorical") {
    # Get counts for each group
    train_counts <- table(y_train)
    test_counts  <- table(y_test)

    # Get all unique levels
    all_levels <- union(names(train_counts), names(test_counts))

    # Create data frame
    dist_df <- data.frame(
      `Dependent Variable` = all_levels,
      stringsAsFactors     = FALSE,
      check.names          = FALSE
    )

    # Add training counts
    dist_df[[paste0("No. of Training Samples (", train_percent, "%)")]] <-
      sapply(all_levels, function(x) {
        if (x %in% names(train_counts)) train_counts[x] else 0
      })

    # Add testing counts
    dist_df[[paste0("No. of Testing Samples (", 100 - train_percent, "%)")]] <-
      sapply(all_levels, function(x) {
        if (x %in% names(test_counts)) test_counts[x] else 0
      })

    # Add total column
    dist_df$Total <- dist_df[[2]] + dist_df[[3]]

    # Add total row
    total_row <- data.frame(
      `Dependent Variable` = "Total",
      stringsAsFactors     = FALSE,
      check.names          = FALSE
    )
    total_row[[paste0("No. of Training Samples (", train_percent, "%)")]]      <- sum(dist_df[[2]])
    total_row[[paste0("No. of Testing Samples (", 100 - train_percent, "%)")]] <- sum(dist_df[[3]])
    total_row$Total                                                            <- sum(dist_df$Total)

    dist_df <- rbind(dist_df, total_row)

  } else {
    # For numeric response, just show totals
    dist_df <- data.frame(
      `Dependent Variable` = "Numeric Response",
      stringsAsFactors     = FALSE,
      check.names          = FALSE
    )

    dist_df[[paste0("No. of Training Samples (", train_percent, "%)")]]      <- length(y_train)
    dist_df[[paste0("No. of Testing Samples (", 100 - train_percent, "%)")]] <- length(y_test)
    dist_df$Total                                                            <- length(y_train) + length(y_test)
  }

  return(dist_df)
}

#' Extract Results from Regularized Regression Analysis
#'
#' @title Extract Specific Results from perform_Regression Output
#'
#' @description
#' This function extracts specific results from the output of \code{perform_Regression}.
#' It can extract all results for a specific alpha value or specific components
#' (like coefficients, performance metrics, etc.) for targeted analysis.
#' When multiple alpha values were tested, this function allows you to access
#' results for any specific alpha or retrieve the best-performing model.
#'
#' @param results A list object containing results from \code{perform_Regression}.
#'   Must contain \code{BestModel}, \code{AllAlphaResults}, and \code{AlphaComparison}.
#' @param alpha Numeric value specifying which alpha value to extract results for.
#'   Must match one of the alpha values used in the original analysis. Common values:
#'   \itemize{
#'     \item{\code{NULL}}: Extracts the best-performing alpha (default)
#'     \item{\code{0}}: Extracts Ridge regression results
#'     \item{\code{1}}: Extracts LASSO regression results
#'     \item{\code{0 < alpha < 1}}: Extracts Elastic Net results for that specific alpha
#'   }
#'   Default: \code{NULL}
#' @param component Character string specifying which component to extract.
#'   Options:
#'   \itemize{
#'     \item{\code{NULL}}: Extract all components for the specified alpha (default)
#'     \item{\code{"model"}}: Extract the fitted cv.glmnet model object
#'     \item{\code{"predictions"}}: Extract actual vs predicted values
#'     \item{\code{"confusion"}}: Extract confusion matrix (classification only)
#'     \item{\code{"performance"}}: Extract performance metrics (Accuracy, RMSE, etc.)
#'     \item{\code{"coefficients"}}: Extract non-zero coefficients and odds ratios
#'     \item{\code{"lambda"}}: Extract selected lambda value
#'     \item{\code{"n_coefficients"}}: Extract the count of non-zero coefficients
#'     \item{\code{"alpha_comparison"}}: Extract comparison table across all alpha values
#'     \item{\code{"summary"}}: Extract key summary information
#'   }
#'   Default: \code{NULL}
#' @param format Character string specifying output format for certain components.
#'   Options:
#'   \itemize{
#'     \item{\code{"list"}}: Return as list (default)
#'     \item{\code{"dataframe"}}: Convert to data frame where possible
#'     \item{\code{"matrix"}}: Convert to matrix where possible
#'   }
#'   Default: \code{"list"}
#' @param verbose Logical indicating whether to print extraction details to console.
#'   Default: \code{TRUE}
#'
#' @return Depends on the parameters specified:
#' \describe{
#'   \item{All components (\code{component = NULL})}{
#'     List containing all result components for the specified alpha:
#'     Model, Performance, ConfusionMatrix, Predictions, Coefficients,
#'     N_Coefficients, Lambda, Alpha
#'   }
#'   \item{Specific component}{
#'     The requested component in the specified format (model object, data frame,
#'     matrix, numeric value, etc.)
#'   }
#'   \item{Alpha comparison (\code{component = "alpha_comparison"})}{
#'     Data frame comparing performance metrics across all tested alpha values
#'   }
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Assume you ran regression with multiple alphas:
#' results <- perform_Regression(
#'   data = preprocessed_data,
#'   alpha = c(0, 0.25, 0.5, 0.75, 1)
#' )
#'
#' # Example 1: Extract best-performing model (automatic selection)
#' best_model <- extract_RegressionResults(
#'   results = results,
#'   alpha = NULL  # NULL = best model
#' )
#'
#' # Example 2: Extract LASSO results (alpha = 1)
#' lasso_results <- extract_RegressionResults(
#'   results = results,
#'   alpha = 1
#' )
#'
#' # Example 3: Extract Ridge results (alpha = 0)
#' ridge_results <- extract_RegressionResults(
#'   results = results,
#'   alpha = 0
#' )
#'
#' # Example 4: Extract only coefficients for alpha = 0.5
#' coefficients_05 <- extract_RegressionResults(
#'   results = results,
#'   alpha = 0.5,
#'   component = "coefficients"
#' )
#'
#' # Example 5: Extract performance metrics for best alpha
#' best_performance <- extract_RegressionResults(
#'   results = results,
#'   alpha = NULL,
#'   component = "performance"
#' )
#'
#' # Example 6: Extract alpha comparison table
#' alpha_comparison <- extract_RegressionResults(
#'   results = results,
#'   component = "alpha_comparison"
#' )
#'
#' # Example 7: Extract predictions as data frame
#' predictions_df <- extract_RegressionResults(
#'   results = results,
#'   alpha = 1,
#'   component = "predictions",
#'   format = "dataframe"
#' )
#'
#' # Example 8: Get number of selected features for LASSO
#' n_features <- extract_RegressionResults(
#'   results = results,
#'   alpha = 1,
#'   component = "n_coefficients"
#' )
#'
#' # Example 9: Compare all alpha values side by side
#' comparison <- extract_RegressionResults(
#'   results = results,
#'   component = "alpha_comparison"
#' )
#' print(comparison)
#' }
#'
#' @seealso \code{\link{perform_Regression}}, \code{\link{list_available_alphas}},
#'   \code{\link{get_best_alpha}}
#'
#' @note
#' \strong{Extracting Results for Different Alpha Values:}
#' \itemize{
#'   \item{Use \code{alpha = NULL} to get the best model automatically}
#'   \item{Use \code{alpha = 1} to get LASSO results}
#'   \item{Use \code{alpha = 0} to get Ridge results}
#'   \item{Use \code{0 < alpha < 1} to get Elastic Net results}
#'   \item{Use \code{component = "alpha_comparison"} to compare all tested alphas}
#' }
#'
#' \strong{Helper Functions:}
#' \itemize{
#'   \item{\code{list_available_alphas(results)}: Lists all alpha values tested}
#'   \item{\code{get_best_alpha(results)}: Returns the best alpha value and its performance}
#' }
#'
#' @export
extract_RegressionResults <- function(
    results,
    alpha     = NULL,
    component = NULL,
    format    = "list",
    verbose   = TRUE
) {

  # Validate inputs
  if (!is.list(results)) {
    stop("'results' must be a list object from perform_Regression.")
  }

  if (!"FunctionOrigin" %in% names(results) || results$FunctionOrigin != "perform_Regression") {
    stop("Results object does not appear to be from perform_Regression function.")
  }

  # Validate component
  valid_components <- c("model", "predictions", "confusion", "performance",
                        "coefficients", "lambda", "n_coefficients",
                        "alpha_comparison", "summary")
  if (!is.null(component) && !component %in% valid_components) {
    stop("'component' must be one of: ", paste(valid_components, collapse = ", "))
  }

  # Check format
  if (!format %in% c("list", "dataframe", "matrix")) {
    stop("'format' must be one of: 'list', 'dataframe', 'matrix'")
  }

  # Handle alpha comparison request
  if (!is.null(component) && component == "alpha_comparison") {
    if (verbose) cat("Extracting alpha comparison data...\n")

    if ("AlphaComparison" %in% names(results) && !is.null(results$AlphaComparison)) {
      return(results$AlphaComparison)
    } else {
      stop("No AlphaComparison found in results.")
    }
  }

  # Determine which alpha to extract
  if (is.null(alpha)) {
    # Extract best model
    if (verbose) cat("Extracting best-performing model...\n")

    if (!"BestModel" %in% names(results)) {
      stop("No BestModel found in results.")
    }

    target_results <- results$BestModel

  } else {
    # Extract specific alpha
    if (verbose) cat("Extracting results for alpha =", alpha, "...\n")

    if (!"AllAlphaResults" %in% names(results)) {
      stop("No AllAlphaResults found. This may be a single-alpha analysis.")
    }

    alpha_key <- paste0("alpha_", alpha)

    if (!alpha_key %in% names(results$AllAlphaResults)) {
      available_alphas <- gsub("alpha_", "", names(results$AllAlphaResults))
      stop("Alpha value ", alpha, " not found. Available alphas: ",
           paste(available_alphas, collapse = ", "))
    }

    target_results <- results$AllAlphaResults[[alpha_key]]
  }

  # Extract specific component if requested
  if (!is.null(component)) {
    return(.extract_component(target_results, component, format, verbose))
  }

  # Return all results
  return(target_results)
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

  # Check component
  valid_components <- c("model", "predictions", "confusion", "performance",
                        "coefficients", "lambda", "n_coefficients", "summary")
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

# Helper function to extract from results
.extract_single_result <- function(results, alpha, component, format, verbose) {

  # Handle alpha comparison request
  if (!is.null(component) && component == "alpha_comparison") {
    if ("AlphaComparison" %in% names(results)) {
      return(results$AlphaComparison)
    }
    stop("No AlphaComparison found in results.")
  }

  # Determine which alpha to extract
  if (is.null(alpha)) {
    # Extract best model
    if (verbose) cat("Extracting best-performing model...\n")
    target_results <- results$BestModel
  } else {
    # Extract specific alpha
    if (verbose) cat("Extracting results for alpha =", alpha, "...\n")

    alpha_key <- paste0("alpha_", alpha)

    if (!"AllAlphaResults" %in% names(results)) {
      stop("No AllAlphaResults found. This may be a single-alpha analysis.")
    }

    if (!alpha_key %in% names(results$AllAlphaResults)) {
      available_alphas <- gsub("alpha_", "", names(results$AllAlphaResults))
      stop("Alpha value ", alpha, " not found. Available alphas: ",
           paste(available_alphas, collapse = ", "))
    }

    target_results <- results$AllAlphaResults[[alpha_key]]
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

         "n_coefficients" = {
           if (verbose) cat("Extracting non-zero coefficient count...\n")
           if ("N_Coefficients" %in% names(results)) {
             return(results$N_Coefficients)
           } else {
             stop("N_Coefficients not found in results.")
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
  if ("Alpha" %in% names(results)) summary_list$Alpha                   <- results$Alpha
  if ("Lambda" %in% names(results)) summary_list$Lambda                 <- results$Lambda
  if ("ResponseType" %in% names(results)) summary_list$ResponseType     <- results$ResponseType
  if ("ReferenceLevel" %in% names(results)) summary_list$ReferenceLevel <- results$ReferenceLevel

  # Performance summary
  if ("Performance" %in% names(results)) {
    perf <- results$Performance
    if ("Accuracy" %in% names(perf)) summary_list$Accuracy   <- perf$Accuracy
    if ("RMSE" %in% names(perf)) summary_list$RMSE           <- perf$RMSE
    if ("R_squared" %in% names(perf)) summary_list$R_squared <- perf$R_squared
    if ("Kappa" %in% names(perf)) summary_list$Kappa         <- perf$Kappa
  }

  # Model characteristics
  if ("N_Coefficients" %in% names(results)) { # <--- Updated to use N_Coefficients
    summary_list$NonZeroCoefficients <- results$N_Coefficients
  }

  if ("Coefficients" %in% names(results)) {
    feature_coefs <- results$Coefficients[results$Coefficients$Feature != "(Intercept)", , drop=FALSE]
    if (nrow(feature_coefs) > 0) {
      # Use head to get top features based on absolute magnitude (a reasonable proxy)
      feature_coefs$AbsCoef    <- abs(feature_coefs$Coefficient)
      feature_coefs            <- feature_coefs[order(feature_coefs$AbsCoef, decreasing = TRUE), ]
      summary_list$TopFeatures <- head(feature_coefs$Feature, 5)
    }
  }

  # Prediction summary
  if ("Predictions" %in% names(results)) {
    pred <- results$Predictions
    if ("Actual" %in% names(pred) && "Predicted" %in% names(pred)) {
      if (is.numeric(pred$Actual)) {
        summary_list$PredictionRange <- range(pred$Predicted, na.rm = TRUE)
        summary_list$ActualRange     <- range(pred$Actual, na.rm = TRUE)
      } else {
        summary_list$CorrectPredictions <- sum(pred$Actual == pred$Predicted)
        summary_list$TotalPredictions   <- nrow(pred)
      }
    }
  }

  return(summary_list)
}

# Convenience function to list all available alphas
list_available_alphas <- function(results) {

  if (!"AllAlphaResults" %in% names(results)) {
    # Single alpha analysis
    if ("BestModel" %in% names(results) && "Alpha" %in% names(results$BestModel)) {
      return(results$BestModel$Alpha)
    }
    stop("No alpha information found in results.")
  }

  alpha_keys <- names(results$AllAlphaResults)
  alphas <- as.numeric(gsub("alpha_", "", alpha_keys))
  return(sort(alphas))
}

# Convenience function to get the best performing alpha
get_best_alpha <- function(results, metric = "auto") {

  if (!"AlphaComparison" %in% names(results) || is.null(results$AlphaComparison)) {
    stop("No AlphaComparison found. This may be a single-alpha analysis.")
  }

  alpha_comp <- results$AlphaComparison

  # Auto-select metric based on what's available
  if (metric == "auto") {
    if ("Accuracy" %in% names(alpha_comp)) {
      metric   <- "Accuracy"
      best_idx <- which.max(alpha_comp[[metric]])
    } else if ("RMSE" %in% names(alpha_comp)) {
      metric   <- "RMSE"
      best_idx <- which.min(alpha_comp[[metric]])
    } else {
      stop("Cannot determine appropriate metric automatically")
    }
  } else {
    if (!metric %in% names(alpha_comp)) {
      stop("Metric '", metric, "' not found in alpha comparison")
    }

    if (metric %in% c("RMSE", "MAE", "MSE")) {
      best_idx <- which.min(alpha_comp[[metric]])
    } else {
      best_idx <- which.max(alpha_comp[[metric]])
    }
  }

  best_alpha <- alpha_comp$Alpha[best_idx]
  return(list(alpha = best_alpha, metric = metric, value = alpha_comp[[metric]][best_idx]))
}

# S3 Methods
#' @export
print.perform_Regression <- function(x, ...) {
  cat("=== Regularized Regression (Lasso/Elastic Net) ===\n")
  cat("Response Type: ", x$ResponseType, "\n")
  cat("Data Source:   ", x$DataSource, "\n")
  cat("Training Split:", round(x$DataSplit$train_samples / x$DataSplit$total_samples * 100, 1), "%\n")

  if (!is.null(x$BestModel)) {
    cat("Best Alpha:    ", x$BestModel$Alpha, "\n")
    if ("N_Coefficients" %in% names(x$BestModel)) {
      cat("Non-Zero Coefs:", x$BestModel$N_Coefficients, "\n")
    }
  }
  invisible(x)
}

#' @export
summary.perform_Regression <- function(object, ...) {
  ans <- list(
    summary    = object$ModelSummary,
    comparison = object$ComparisonSummary,
    type       = object$ResponseType
  )
  class(ans) <- "summary.perform_Regression"
  return(ans)
}

#' @export
print.summary.perform_Regression <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("Regression Model Summary\n")
  cat("---------------------------------------\n")
  print(x$summary, row.names = FALSE)

  if(!is.null(x$comparison)) {
    cat("\n-- Method Comparison --\n")
    print(x$comparison, row.names = FALSE)
  }
  invisible(x)
}
