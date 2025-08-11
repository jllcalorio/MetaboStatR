#' Perform Least Absolute Shrinkage and Selection Operator (LASSO) and Elastic Net Regression
#'
#' @description
#' This function performs either or both the LASSO and/or Elastic Net Regression. Both are regularization
#' techniques used to prevent overfitting in linear regression models by adding penalty terms to the loss function.
#' LASSO uses L1 regularization to shrink the coefficients of less important features towards zero, potentially eliminating them entirely.
#' Elastic Net uses a combination of L1 and L2 regularization to address multicollinearity and perform feature selection simultaneously.
#'
#' @param data List. This list must be a result from the `perform_PreprocessingPeakData` function.
#' @param method String or vector.
#'   \itemize{
#'     \item "lasso": Analyses data using LASSO only.
#'     \item "enet": Analyses data using Elastic Net only.
#'     \item c("lasso", "enet"): Analyses data using both LASSO Elastic Net. Needs to be a vector since 'Ridge Regression' will be added.
#'     }
#'     Defaults to "enet".
#' @param specify_response String. A column
#' @param train_percent Numeric. The percent to be used in the training data set. The remaining will be used in testing. Selection will be done in a random manner. For consistent results, suggests to use the `remember` parameter.
#' @param ref String. The reference value. Defaults to the 1st level/category of the dependent variable if `NULL`.
#' @param lambda String. The lambda value.
#'   \itemize{
#'     \item "1se"
#'     \item "min"
#'     }
#'     Defaults to "1se" for fewer features, lower risk of over fitting.
#' @param remember Numeric. This value will be used in the `base::set.seed(remember)` function.
#'
#' @returns Data frame of regression results containing coefficients and odds ratios.
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform the function
#' myregression <- performRegression(
#'   data     = results_from_perform_PreprocessingPeakData_function,
#'   remember = 1234
#'   )
#'
#' # View the results
#' View(myregression[["LASSO_coefOddsRatio"]])
#' View(myregression[["ElasticNet_coefOddsRatio"]])
#' myregression[["LASSO_confusionMatrix"]]$overall[1]
#' myregression[["ElasticNet_confusionMatrix"]]$overall[1]
#' }
#'
performRegression <- function(
    data,
    method           = "enet",
    specify_response = NULL,
    train_percent    = 80,
    ref              = NULL,
    lambda           = "1se",
    remember         = NULL
) {

  # Parameter checks
  if (!base::is.numeric(train_percent) || train_percent <= 0 || train_percent > 100) {
    stop("train_percent must be a numeric value between 0 and 100.")
  }

  if (data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("The supplied data did not come from 'perform_PreprocessingPeakData' function. Required metadata might be missing. Set 'ignore = TRUE' if your data is correctly formatted.")
  }

  allowed_scales <- c("min", "1se")
  if (!(lambda %in% allowed_scales)) {
    errors <- base::c(errors, "lambda: Must be either min or 1se. min = lambda value that achieves the minimum cross-validated error, 1se = lambda value which the cross-validated error is within one standard error of the minimum.")
  }

  allowed_methods <- c("lasso", "enet")
  if (!base::is.null(method)) {
    if (!base::all(method %in% allowed_methods)) {
      stop("method is invalid. Must be either 'lasso', 'enet', or both.")
    }
  } else { # if NULL
    stop("method must not be NULL. Must be either 'lasso', 'enet', or both.")
  }

  if (!base::is.null(remember)) { # perform base::set.seed if remember = any number
    if (base::is.numeric(remember)) {
      base::set.seed(remember)
    } else {
      stop("'remember' parameter must be numeric.")
    }
  }

  qc_indices      <- data$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices  <- !qc_indices
  groups          <- data$Metadata$Group[non_qc_indices] # This is to be replaced in the latter part of the script

  # Check if arrangeLevels have the same actual values of groups
  if (!base::is.null(ref)) {
    if (base::length(base::setdiff(ref, base::unique(groups))) > 0) {
      stop("Invalid ref: check for typos or missing groups.")
    } else {
      ref    <- ref[1] # The ref[1] ensures that it only takes the 1st value in the vector, in the case the user inputs more than 1
      groups <- base::factor(groups, levels = c(ref,
                                                base::setdiff(base::unique(data$Metadata$Group[non_qc_indices]), ref)))
    }
  } else {
    groups <- base::factor(groups)
  }

  num_groups <- dplyr::n_distinct(groups)

  df         <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ] # Extract non-QC data

  train_size <- base::floor((train_percent / 100) * base::dim(df)[1])
  train_rows <- base::sample(1:dim(df)[1], train_size)

  x_train    <- df[ train_rows, ] %>% base::as.matrix()
  x_test     <- df[-train_rows, ] %>% base::as.matrix()
  y_train    <- groups[ train_rows]
  y_test     <- groups[-train_rows]

  # #####################################################
  # #####################################################
  # #####################################################
  #
  # List of results
  regression_results                  <- base::list() # Create empty list

  regression_results$FunctionOrigin   <- "performRegression"

  regression_results$data             <- df

  regression_results$train_percent    <- train_percent %>% base::as.data.frame() %>% base::`colnames<-`(NULL)
  regression_results$data_indep_train <- x_train       %>% base::as.data.frame() %>% base::`colnames<-`(NULL)
  regression_results$data_indep_test  <- x_test        %>% base::as.data.frame() %>% base::`colnames<-`(NULL)
  regression_results$data_dep_train   <- y_train       %>% base::as.data.frame() %>% base::`colnames<-`(NULL)
  regression_results$data_dep_test    <- y_test        %>% base::as.data.frame() %>% base::`colnames<-`(NULL)

  # Loop through the chosen methods
  for (method_type in method) {

    cat("########################################################\n")
    cat("RESULTS FOR", base::toupper(method_type), "\n")

    alpha_val   <- base::ifelse(method_type == "lasso",          1, 0.5)
    family_type <- base::ifelse(num_groups  ==       2, "binomial", "multinomial")

    model_fit <- glmnet::cv.glmnet(
      x            = x_train,
      y            = y_train,
      alpha        = alpha_val,
      family       = family_type,
      type.measure = "class"
    )

    predicted <- stats::predict(
      model_fit,
      # s = if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min,
      s    = base::ifelse(lambda == "1se", model_fit$lambda.1se, model_fit$lambda.min),
      newx = x_test,
      type = "class"
    )

    confMat <- caret::confusionMatrix(
      base::factor(y_test),
      base::factor(predicted, levels = base::levels(groups))
    )

    base::print(confMat)

    coef_model <- glmnet::coef.glmnet(model_fit, s = if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min)

    coef.OR_list <- if (num_groups == 2) {
      coef_df              <- base::as.data.frame(base::as.matrix(coef_model))
      colnames(coef_df)    <- "Coefficient"
      coef_df$Feature      <- base::rownames(coef_df)
      coef_df              <- coef_df[coef_df$Coefficient != 0, ]
      coef_df$`Odds Ratio` <- base::exp(coef_df$Coefficient)

      # list(coef_df[, c("Coefficient", "Odds Ratio")])
      coef.OR_list <- coef_df[, c("Coefficient", "Odds Ratio")]

    } else {
      group_levels <- base::levels(groups)

      # Combine all group-specific coefficient tables into one data frame
      coef_list <- base::lapply(base::seq_along(coef_model), function(i) {
        coef_df              <- base::as.data.frame(base::as.matrix(coef_model[[i]]))
        colnames(coef_df)    <- "Coefficient"
        coef_df$Feature      <- base::rownames(coef_df)
        coef_df              <- coef_df[coef_df$Coefficient != 0, ]
        coef_df$`Odds Ratio` <- base::exp(coef_df$Coefficient)
        coef_df$Group        <- group_levels[i]
        coef_df[, c("Group", "Feature", "Coefficient", "Odds Ratio")]
      })

      base::do.call(rbind, coef_list) %>% base::`rownames<-`(NULL)
    }

    base::print(coef.OR_list)
    base::cat("Reference: ", base::ifelse(base::is.null(ref), base::levels(groups)[1], ref), "\n")                              # If reference is NULL, the reference value used is the 1st level category

    # Save results
    if (method_type == "lasso") {
      regression_results$LASSO_fit                    <- model_fit
      regression_results$LASSO_train_vs_predict       <- base::data.frame(Training = y_test, Predicted = base::as.vector(predicted))
      regression_results$LASSO_confMatrix_All         <- confMat
      regression_results$LASSO_confMatrix_Table       <- confMat$table   %>% base::as.data.frame()                              # Confusion Matrix table
      regression_results$LASSO_OverallStatistics      <- confMat$overall %>% base::as.data.frame() %>% base::`colnames<-`(NULL) # Overall statistics, from accuracy, CI, p-value etc.
      regression_results$LASSO_SensSpecsPosNeg        <- confMat$byClass %>% base::t()             %>% base::as.data.frame()    # Sensitivity, Specificity, Neg/Pos Pred Value
      regression_results$LASSO_coefOddsRatio          <- coef.OR_list    %>% base::`rownames<-`(NULL)
    } else {
      regression_results$ElasticNet_fit               <- model_fit
      regression_results$ElasticNet_train_vs_predict  <- base::data.frame(Training = y_test, Predicted = base::as.vector(predicted))
      regression_results$ElasticNet_confMatrix_All    <- confMat
      regression_results$ElasticNet_confMatrix_Table  <- confMat$table   %>% base::as.data.frame()                              # Confusion Matrix table
      regression_results$ElasticNet_OverallStatistics <- confMat$overall %>% base::as.data.frame() %>% base::`colnames<-`(NULL) # Overall statistics, from accuracy, CI, p-value etc.
      regression_results$ElasticNet_SensSpecsPosNeg   <- confMat$byClass %>% base::t() %>% base::as.data.frame()                # Sensitivity, Specificity, Neg/Pos Pred Value
      regression_results$ElasticNet_coefOddsRatio     <- coef.OR_list    %>% base::`rownames<-`(NULL)
    }
  }

  return(regression_results)
}
