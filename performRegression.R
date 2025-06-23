# Function for LASSO and Elastic Net
performRegression <- function(
    data,
    method        = c("lasso", "enet"),  # Choose either "lasso", "enet", or both
    train_percent = 80,   # percent of data to be used in training
    ref           = NULL, # The reference value. Defaults to the 1st level/category of the dependent variable
    lambda        = "1se", # The lambda value. c("1se", "min") Defaults to "1se" for fewer features, lower risk of over fitting
    remember      = NULL  # value used in set.seed
) {

  # Parameter checks
  if (!is.numeric(train_percent) || train_percent <= 0 || train_percent > 100) {
    stop("train_percent must be a numeric value between 0 and 100.")
  }

  if (data$Metadata$FunctionOrigin != "performPreprocessingPeakData") {
    stop("The supplied data did not come from 'performPreprocessingPeakData' function. Required metadata might be missing. Set 'ignore = TRUE' if your data is correctly formatted.")
  }

  allowed_scales <- c("min", "1se")
  if (!(lambda %in% allowed_scales)) {
    errors <- c(errors, "lambda: Must be either min or 1se. min = lambda value that achieves the minimum cross-validated error, 1se = lambda value which the cross-validated error is within one standard error of the minimum.")
  }

  allowed_methods <- c("lasso", "enet")
  if (!is.null(method)) {
    if (!all(method %in% allowed_methods)) {
      stop("method is invalid. Must be either 'lasso', 'enet', or both.")
    }
  } else { # if NULL
    stop("method must not be NULL. Must be either 'lasso', 'enet', or both.")
  }


  if (!is.null(remember)) { # perform set.seed if remember = any number
    if (is.numeric(remember)) {
      set.seed(remember)
    } else {
      stop("'remember' parameter must be numeric.")
    }
  }

  non_qc_indices  <- data$Metadata$Groups != "QC"
  groups          <- data$Metadata$Groups[non_qc_indices] # This is to be replaced in the latter part of the script

  # Check if arrangeLevels have the same actual values of groups
  if (!is.null(ref)) {
    if (length(setdiff(ref, unique(groups))) > 0) {
      stop("Invalid ref: check for typos or missing groups.")
    } else {
      ref <- ref[1] # The ref[1] ensures that it only takes the 1st value in the vector, in the case the user inputs more than 1
      groups <- factor(groups, levels = c(ref,
                                          setdiff(unique(data$Metadata$Groups[non_qc_indices]), ref)))
    }
  } else {
    groups <- factor(groups)
  }

  num_groups      <- dplyr::n_distinct(groups)

  # if (dplyr::n_distinct(groups) != 2) {
  #   stop("There are more than 2 groups.")
  # } else if (dplyr::n_distinct(groups) == 1) {
  #   stop("There is only 1 group.")
  # }

  df              <- data$data_scaledOPLSDA[non_qc_indices, ] # Extract non-QC data



  train_size      <- floor((train_percent / 100) * dim(df)[1])
  train_rows      <- sample(1:dim(df)[1], train_size)

  x_train         <- df[ train_rows, ] %>% as.matrix()
  x_test          <- df[-train_rows, ] %>% as.matrix()
  y_train         <- groups[ train_rows]
  y_test          <- groups[-train_rows]

  # #####################################################
  # #####################################################
  # #####################################################
  #
  # List of results
  regression_results = list() # Create empty list

  regression_results$Class            <- "performRegression"

  regression_results$data             <- df

  regression_results$train_percent    <- train_percent %>% as.data.frame() %>% `colnames<-`(NULL)
  regression_results$data_indep_train <- x_train %>% as.data.frame() %>% `colnames<-`(NULL)
  regression_results$data_indep_test  <- x_test %>% as.data.frame() %>% `colnames<-`(NULL)
  regression_results$data_dep_train   <- y_train %>% as.data.frame() %>% `colnames<-`(NULL)
  regression_results$data_dep_test    <- y_test %>% as.data.frame() %>% `colnames<-`(NULL)
  #
  #
  # #####################################################
  # #####################################################
  # #####################################################
  #
  # # LASSO Regression (alpha = 1)
  # lasso_fit <- glmnet::cv.glmnet(
  #   x            = x_train,
  #   y            = y_train,
  #   alpha        = 1, # must be 1 for LASSO
  #   family       = ifelse(num_groups == 2, "binomial", "multinomial"),
  #   type.measure = "class"
  # )
  #
  # lasso_predicted <- predict(
  #   lasso_fit,
  #   # lambda.1se for fewer metabolites
  #   # lambda.1se for more metabolites
  #   s    = lasso_fit$lambda.1se,
  #   newx = x_test, # Matrix of values which predictions are to be made
  #   type = "class" # "class" for binomial/multinomial logistic regression
  # ) # predictions return the class labels for classification
  #
  # # mean(y_test == lasso_predicted)
  #
  #
  #
  # confMat_lasso   <- caret::confusionMatrix(factor(y_test), factor(lasso_predicted, levels = c(ref,
  #                                                                                              setdiff(unique(data$Metadata$Groups[non_qc_indices]), ref))))
  # print("RESULTS FOR LASSO")
  # # print(paste0("LASSO Accuracy: ", round(confMat_lasso$overall[[1]], 3),
  # #              " (95% CI: ", round(confMat_lasso$overall[[3]], 3), "-", round(confMat_lasso$overall[[4]], 3), ")"))
  #
  #
  # print(confMat_lasso)
  #
  # # coef_lasso      <- as.matrix(coef(lasso_fit, s = lasso_fit$lambda.1se)) %>%
  # #   .[. != 0, , drop = FALSE]
  #
  # # Get coefficients
  # coef_lasso      <- glmnet::coef.glmnet(lasso_fit, s = lasso_fit$lambda.1se)
  # # Loop through the list (one per class/group)
  # coef.OR_list <- if (num_groups == 2) {
  #   coef_df <- as.data.frame(as.matrix(coef_lasso))
  #   colnames(coef_df) <- "Coefficient"
  #   coef_df$Feature <- rownames(coef_df)
  #   coef_df <- coef_df[coef_df$Coefficient != 0, ]
  #   coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
  #   list(coef_df[, c("Coefficient", "Odds Ratio")])
  # } else {
  #   lapply(coef_lasso, function(mat) {
  #     coef_df <- as.data.frame(as.matrix(mat))
  #     colnames(coef_df) <- "Coefficient"
  #     coef_df$Feature <- rownames(coef_df)
  #     coef_df <- coef_df[coef_df$Coefficient != 0, ]
  #     coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
  #     coef_df[, c("Coefficient", "Odds Ratio")]
  #   })
  # }
  #
  # print(coef.OR_list)
  # print(paste0("Reference: ", ifelse(is.null(ref), groups[1], ref)))
  #
  #
  #
  # regression_results$LASSO_fit              <- lasso_fit
  # regression_results$LASSO_train_vs_predict <- data.frame(Training = y_test, Predicted = as.vector(lasso_predicted))
  # regression_results$LASSO_confusionMatrix  <- confMat_lasso
  # regression_results$LASSO_coefOddsRatio    <- coef.OR_list
  #
  # #####################################################
  # #####################################################
  # #####################################################
  #
  # # Elastic Net Regression (alpha = .5)
  # enet_fit <- glmnet::cv.glmnet(
  #   x            = x_train,
  #   y            = y_train,
  #   alpha        = .5, # must be between 0 and 1 for Elastic Net
  #   family       = ifelse(num_groups == 2, "binomial", "multinomial"),
  #   type.measure = "class"
  # )
  #
  # enet_predicted <- predict(
  #   enet_fit,
  #   # lambda.1se for fewer metabolites
  #   # lambda.1se for more metabolites
  #   s    = enet_fit$lambda.1se,
  #   newx = x_test, # Matrix of values which predictions are to be made
  #   type = "class" # class for binomial/multinomial logistic regression
  # ) # predictions return the class labels for classification
  #
  #
  # # mean(y_test == enet_predicted)
  #
  # confMat_enet   <- caret::confusionMatrix(factor(y_test), factor(enet_predicted, levels = c(ref,
  #                                                                                            setdiff(unique(data$Metadata$Groups[non_qc_indices]), ref))))
  # print("########################################################")
  # print("########################################################")
  # print("########################################################")
  #
  # print("RESULTS FOR ELASTIC NET")
  # # print(paste0("Elastic Net Accuracy: ", round(confMat_enet$overall[[1]], 3),
  # #              " (95% CI: ", round(confMat_enet$overall[[3]], 3), "-", round(confMat_enet$overall[[4]], 3), ")"))
  #
  #
  # print(confMat_enet)
  #
  # # Get coefficients
  # coef_enet      <- glmnet::coef.glmnet(enet_fit, s = enet_fit$lambda.1se)
  # # Loop through the list (one per class/group)
  # coef.OR_list <- if (num_groups == 2) {
  #   coef_df <- as.data.frame(as.matrix(coef_lasso))
  #   colnames(coef_df) <- "Coefficient"
  #   coef_df$Feature <- rownames(coef_df)
  #   coef_df <- coef_df[coef_df$Coefficient != 0, ]
  #   coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
  #   list(coef_df[, c("Coefficient", "Odds Ratio")])
  # } else {
  #   lapply(coef_lasso, function(mat) {
  #     coef_df <- as.data.frame(as.matrix(mat))
  #     colnames(coef_df) <- "Coefficient"
  #     coef_df$Feature <- rownames(coef_df)
  #     coef_df <- coef_df[coef_df$Coefficient != 0, ]
  #     coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
  #     coef_df[, c("Coefficient", "Odds Ratio")]
  #   })
  # }
  #
  # print(coef.OR_list)
  # print(paste0("Reference: ", ifelse(is.null(ref), groups[1], ref)))
  #
  # regression_results$ElasticNet_fit              <- enet_fit
  # regression_results$ElasticNet_train_vs_predict <- data.frame(Training = y_test, Predicted = as.vector(enet_predicted))
  # regression_results$ElasticNet_confusionMatrix  <- confMat_enet
  # regression_results$ElasticNet_coefOddsRatio    <- coef.OR_list
  #
  # #####################################################
  # #####################################################
  # #####################################################
  #
  # return(regression_results)

  # Loop through the chosen methods
  for (method_type in method) {

    cat("########################################################\n")
    cat("RESULTS FOR", toupper(method_type), "\n")

    alpha_val   <- ifelse(method_type == "lasso", 1, 0.5)
    family_type <- ifelse(num_groups  == 2, "binomial", "multinomial")

    model_fit <- glmnet::cv.glmnet(
      x            = x_train,
      y            = y_train,
      alpha        = alpha_val,
      family       = family_type,
      type.measure = "class"
    )

    predicted <- predict(
      model_fit,
      # s = if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min,
      s    = ifelse(lambda == "1se", model_fit$lambda.1se, model_fit$lambda.min),
      newx = x_test,
      type = "class"
    )

    confMat <- caret::confusionMatrix(
      factor(y_test),
      factor(predicted, levels = levels(groups))
    )

    print(confMat)

    coef_model <- glmnet::coef.glmnet(model_fit, s = if (lambda == "1se") model_fit$lambda.1se else model_fit$lambda.min)

    coef.OR_list <- if (num_groups == 2) {
      coef_df              <- as.data.frame(as.matrix(coef_model))
      colnames(coef_df)    <- "Coefficient"
      coef_df$Feature      <- rownames(coef_df)
      coef_df              <- coef_df[coef_df$Coefficient != 0, ]
      coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)

      # list(coef_df[, c("Coefficient", "Odds Ratio")])
      coef.OR_list <- coef_df[, c("Coefficient", "Odds Ratio")]

    } else {
      # lapply(coef_model, function(mat) {
      #   coef_df              <- as.data.frame(as.matrix(mat))
      #   colnames(coef_df)    <- "Coefficient"
      #   coef_df$Feature      <- rownames(coef_df)
      #   coef_df              <- coef_df[coef_df$Coefficient != 0, ]
      #   coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
      #   coef_df[, c("Coefficient", "Odds Ratio")]
      # })

      # group_levels <- levels(groups)
      #
      # for (i in seq_along(coef_model)) {
      #   coef_df              <- as.data.frame(as.matrix(coef_model[[i]]))
      #   colnames(coef_df)    <- "Coefficient"
      #   coef_df$Feature      <- rownames(coef_df)
      #   coef_df              <- coef_df[coef_df$Coefficient != 0, ]
      #   coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
      #
      #   # Assign to separate data frames named by group level
      #   assign(paste0("coef_OR_", group_levels[i]), coef_df[, c("Coefficient", "Odds Ratio")], envir = .GlobalEnv)
      # }

      # group_levels <- levels(groups)
      #
      # # Create a named list of data frames
      # coef_list <- lapply(seq_along(coef_model), function(i) {
      #   coef_df              <- as.data.frame(as.matrix(coef_model[[i]]))
      #   colnames(coef_df)    <- "Coefficient"
      #   coef_df$Feature      <- rownames(coef_df)
      #   coef_df              <- coef_df[coef_df$Coefficient != 0, ]
      #   coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
      #   coef_df[, c("Coefficient", "Odds Ratio")]
      # })

      # names(coef_list) <- group_levels
      # coef_list

      group_levels <- levels(groups)

      # Combine all group-specific coefficient tables into one data frame
      coef_list <- lapply(seq_along(coef_model), function(i) {
        coef_df              <- as.data.frame(as.matrix(coef_model[[i]]))
        colnames(coef_df)    <- "Coefficient"
        coef_df$Feature      <- rownames(coef_df)
        coef_df              <- coef_df[coef_df$Coefficient != 0, ]
        coef_df$`Odds Ratio` <- exp(coef_df$Coefficient)
        coef_df$Group        <- group_levels[i]
        coef_df[, c("Group", "Feature", "Coefficient", "Odds Ratio")]
      })

      do.call(rbind, coef_list) %>% `rownames<-`(NULL)



    }

    print(coef.OR_list)
    cat("Reference: ", ifelse(is.null(ref), levels(groups)[1], ref), "\n") # If reference is NULL, the reference value used is the 1st level category

    # Save results
    if (method_type == "lasso") {
      regression_results$LASSO_fit                    <- model_fit
      regression_results$LASSO_train_vs_predict       <- data.frame(Training = y_test, Predicted = as.vector(predicted))
      regression_results$LASSO_confMatrix_All         <- confMat
      regression_results$LASSO_confMatrix_Table       <- confMat$table %>% as.data.frame() # Confusion Matrix table
      regression_results$LASSO_OverallStatistics      <- confMat$overall %>% as.data.frame() %>% `colnames<-`(NULL) # Overall statistics, from accuracy, CI, p-value etc.
      regression_results$LASSO_SensSpecsPosNeg        <- confMat$byClass %>% t() %>% as.data.frame() # Sensitivity, Specificity, Neg/Pos Pred Value
      regression_results$LASSO_coefOddsRatio          <- coef.OR_list %>% `rownames<-`(NULL)
    } else {
      regression_results$ElasticNet_fit               <- model_fit
      regression_results$ElasticNet_train_vs_predict  <- data.frame(Training = y_test, Predicted = as.vector(predicted))
      regression_results$ElasticNet_confMatrix_All    <- confMat
      regression_results$ElasticNet_confMatrix_Table  <- confMat$table %>% as.data.frame() # Confusion Matrix table
      regression_results$ElasticNet_OverallStatistics <- confMat$overall %>% as.data.frame() %>% `colnames<-`(NULL) # Overall statistics, from accuracy, CI, p-value etc.
      regression_results$ElasticNet_SensSpecsPosNeg   <- confMat$byClass %>% t() %>% as.data.frame() # Sensitivity, Specificity, Neg/Pos Pred Value
      regression_results$ElasticNet_coefOddsRatio     <- coef.OR_list %>% `rownames<-`(NULL)
    }
  }

  return(regression_results)

  # regression_results[[paste0(method_type, "_ConfMatrix")]] <- confMat$table # Confusion Matrix table
  # regression_results[[paste0(method_type, "_OverallStatistics")]] <- confMat$overall # Overall statistics, from accuracy, CI, p-value etc.
  # regression_results[[paste0(method_type, "_SensSpecsPosNeg")]] <- confMat$byClass # Sensitivity, Specificity, Neg/Pos Pred Value


}

# Usage # 1258
# myregression <- performRegression(mydata, remember = 1258); View(myregression[["LASSO_coefOddsRatio"]]); View(myregression[["ElasticNet_coefOddsRatio"]]); myregression[["LASSO_confusionMatrix"]]$overall[1]; myregression[["ElasticNet_confusionMatrix"]]$overall[1]
