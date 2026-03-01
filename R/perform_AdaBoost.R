#' Perform Adaptive Boosting (AdaBoost) Classification on Preprocessed Metabolomics Data
#'
#' Uses the AdaBoost.M1 algorithm from the \code{adabag} package. The number
#' of boosting iterations (\code{mfinal}) can be tuned automatically via
#' cross-validation or set manually. Variable importance is derived from the
#' weighted aggregation of base learner importances across all boosting rounds.
#'
#' @param prepped_data Output from perform_PreprocessingPeakData function
#' @param outcome_var Character string specifying the outcome variable column
#'        name in metadata (default: "Group")
#' @param features Character vector of feature/metabolite names to include.
#'        If NULL, uses all available features (default: NULL)
#' @param data_type Character string specifying which data matrix to use:
#'        "NONPLS" or "PLS" (default: "NONPLS")
#' @param use_merged Logical indicating whether to use merged replicates if
#'        available (default: TRUE)
#' @param tune_mfinal Logical indicating whether to automatically tune the
#'        number of boosting iterations via cross-validation (default: TRUE)
#' @param mfinal Integer or integer vector of iteration values to search when
#'        \code{tune_mfinal = TRUE}, or a single value when
#'        \code{tune_mfinal = FALSE} (default: \code{c(50, 100, 150, 200)})
#' @param cross Integer specifying the number of cross-validation folds used
#'        during tuning (default: 10)
#' @param maxdepth Integer specifying the maximum depth of the base decision
#'        tree learners (default: 1, i.e. decision stumps)
#' @param coeflearn Character string specifying the coefficient update rule:
#'        \code{"Breiman"} (default), \code{"Freund"}, or \code{"Zhu"}
#' @param seed Integer for reproducibility (default: 42)
#'
#' @return Object of class "perform_AdaBoost" (list) containing:
#'   \itemize{
#'     \item \code{model} – Final adaboost object
#'     \item \code{accuracy} – Cross-validated accuracy at optimal mfinal
#'     \item \code{confusion_matrix} – Confusion matrix with class error rates
#'     \item \code{probabilities} – Matrix of class probabilities per sample
#'     \item \code{variable_importance} – Named numeric vector of variable
#'           importance scores aggregated across boosting rounds
#'     \item \code{error_trajectory} – Data frame of training error per
#'           boosting iteration
#'     \item \code{tuning_results} – Data frame of mfinal vs CV accuracy
#'           (NULL if \code{tune_mfinal = FALSE})
#'     \item \code{optimal_mfinal} – Optimal or specified mfinal value
#'     \item \code{metadata} – Metadata used (QC samples removed)
#'     \item \code{data_matrix} – Data matrix used
#'     \item \code{outcome_var} – Name of outcome variable
#'     \item \code{features_used} – Character vector of features in model
#'     \item \code{parameters} – List of function call parameters
#'   }
#'
#' @examples
#' ada_results <- perform_AdaBoost(df_prepped, outcome_var = "Group")
#'
#' # Manual mfinal, no tuning
#' ada_results <- perform_AdaBoost(df_prepped, outcome_var = "Group",
#'                                  tune_mfinal = FALSE, mfinal = 100)
#'
#' # Deeper base learners with custom iteration grid
#' ada_results <- perform_AdaBoost(df_prepped,
#'                                  mfinal   = c(50, 100, 200, 300, 500),
#'                                  maxdepth = 2)
#' @export
perform_AdaBoost <- function(prepped_data,
                             outcome_var  = "Group",
                             features     = NULL,
                             data_type    = "NONPLS",
                             use_merged   = TRUE,
                             tune_mfinal  = TRUE,
                             mfinal       = c(50, 100, 150, 200),
                             cross        = 10,
                             maxdepth     = 1,
                             coeflearn    = "Breiman",
                             seed         = 42) {

  # ---- package checks ----
  if (!requireNamespace("adabag", quietly = TRUE))
    stop("Package 'adabag' is required but not installed.")
  if (!requireNamespace("rpart", quietly = TRUE))
    stop("Package 'rpart' is required but not installed.")

  # ---- input validation ----
  if (!inherits(prepped_data, "perform_PreprocessingPeakData"))
    stop("Input must be output from perform_PreprocessingPeakData function")

  if (!data_type %in% c("NONPLS", "PLS"))
    stop("data_type must be either 'NONPLS' or 'PLS'")

  coeflearn <- match.arg(coeflearn, c("Breiman", "Freund", "Zhu"))

  if (any(mfinal < 1) || any(mfinal != floor(mfinal)))
    stop("All mfinal values must be positive integers.")

  # ---- select data and metadata ----
  if (use_merged && !is.null(prepped_data$data_scaledNONPLS_merged)) {
    data_matrix <- if (data_type == "NONPLS") prepped_data$data_scaledNONPLS_merged
    else                        prepped_data$data_scaledPLS_merged
    metadata    <- prepped_data$Metadata_merged
    message("Using merged replicate data")
  } else {
    if (data_type == "NONPLS") {
      data_matrix <- prepped_data$data_scaledNONPLS_varFiltered %||%
        prepped_data$data_scaledNONPLS
    } else {
      data_matrix <- prepped_data$data_scaledPLS_varFiltered %||%
        prepped_data$data_scaledPLS
    }
    metadata <- prepped_data$Metadata
    message("Using non-merged data")
  }

  # ---- validate outcome variable ----
  if (!outcome_var %in% colnames(metadata))
    stop(paste0("Outcome variable '", outcome_var, "' not found in metadata. ",
                "Available columns: ", paste(colnames(metadata), collapse = ", ")))

  # ---- remove QC samples ----
  if ("Group_" %in% colnames(metadata)) {
    qc_idx <- which(metadata$Group_ == "QC")
    if (length(qc_idx) > 0) {
      message(paste0("Removing ", length(qc_idx), " QC samples from analysis"))
      data_matrix <- data_matrix[-qc_idx, , drop = FALSE]
      metadata    <- metadata[-qc_idx,    , drop = FALSE]
    }
  }

  # ---- feature selection ----
  if (!is.null(features)) {
    available   <- colnames(data_matrix)
    missing_ft  <- setdiff(features, available)
    if (length(missing_ft) > 0)
      warning(paste0("Features not found and excluded: ",
                     paste(missing_ft, collapse = ", ")))
    features_to_use <- intersect(features, available)
    if (length(features_to_use) == 0)
      stop("None of the specified features were found in the data matrix.")
    message(paste0("Using ", length(features_to_use), " out of ",
                   length(features), " requested features"))
    data_matrix <- data_matrix[, features_to_use, drop = FALSE]
  } else {
    features_to_use <- colnames(data_matrix)
  }

  # ---- handle missing values in outcome ----
  outcome <- as.factor(metadata[[outcome_var]])
  if (any(is.na(outcome))) {
    warning("Missing values found in outcome variable. Removing those samples.")
    valid_idx   <- !is.na(outcome)
    data_matrix <- data_matrix[valid_idx, , drop = FALSE]
    metadata    <- metadata[valid_idx,    , drop = FALSE]
    outcome     <- outcome[valid_idx]
  }

  # Sanitise column names — adabag uses formula interface; make.names()
  # is applied explicitly so train and predict-time names always match.
  colnames(data_matrix) <- make.names(colnames(data_matrix), unique = TRUE)
  features_to_use       <- make.names(features_to_use,       unique = TRUE)
  model_data            <- data.frame(outcome = outcome, data_matrix)

  # rpart control for base learners
  rpart_ctrl <- rpart::rpart.control(maxdepth = maxdepth)

  set.seed(seed)

  # ---- tune mfinal via k-fold CV ----
  tuning_results       <- NULL
  optimal_mfinal       <- NULL
  cv_trajectory        <- NULL   # CV error at every iteration (for trajectory plot)

  if (tune_mfinal) {
    message(paste0("Tuning mfinal via ", cross, "-fold cross-validation..."))

    n        <- nrow(model_data)
    fold_ids <- sample(rep(seq_len(cross), length.out = n))
    max_mf   <- max(mfinal)

    # Run each fold once up to max(mfinal), collecting per-iteration predictions
    # via errorevol() so CV error at every iteration is available cheaply.
    message("Computing per-iteration CV error trajectory (this may take a moment)...")
    fold_iter_errors <- lapply(seq_len(cross), function(f) {
      train <- model_data[fold_ids != f, , drop = FALSE]
      test  <- model_data[fold_ids == f, , drop = FALSE]
      fit   <- adabag::boosting(outcome ~ .,
                                data      = train,
                                mfinal    = max_mf,
                                coeflearn = coeflearn,
                                control   = rpart_ctrl)
      # errorevol returns error at each iteration on test set
      adabag::errorevol(fit, newdata = test)$error
    })

    # Average CV error across folds at each iteration
    cv_iter_errors <- rowMeans(do.call(cbind, fold_iter_errors))
    cv_trajectory  <- data.frame(
      Iteration = seq_len(max_mf),
      CV_Error  = cv_iter_errors
    )

    # Summarise CV accuracy at each candidate mfinal
    cv_acc <- vapply(mfinal, function(mf) 1 - cv_iter_errors[mf], numeric(1))

    tuning_results <- data.frame(
      mfinal      = mfinal,
      cv_accuracy = cv_acc * 100,
      cv_error    = (1 - cv_acc) * 100
    )

    optimal_mfinal <- mfinal[which.max(cv_acc)]
    message(paste0("Optimal mfinal: ", optimal_mfinal,
                   "  (CV accuracy: ", round(max(cv_acc) * 100, 2), "%)"))

  } else {
    optimal_mfinal <- mfinal[1]
    message(paste0("Using specified mfinal: ", optimal_mfinal))
  }

  # ---- build final model ----
  message("Building final AdaBoost model...")
  final_model <- adabag::boosting(outcome ~ .,
                                  data      = model_data,
                                  mfinal    = optimal_mfinal,
                                  coeflearn = coeflearn,
                                  control   = rpart_ctrl)

  # ---- predictions and probabilities ----
  pred_result <- predict(final_model, newdata = model_data)
  predicted   <- as.factor(pred_result$class)
  probs       <- pred_result$prob
  colnames(probs) <- levels(outcome)

  # ---- confusion matrix ----
  conf_table  <- table(Predicted = predicted, Actual = outcome)
  class_err   <- 1 - diag(conf_table) / colSums(conf_table)
  conf_matrix <- cbind(conf_table, class.error = round(class_err, 4))

  cv_accuracy <- if (!is.null(tuning_results)) {
    max(tuning_results$cv_accuracy)
  } else {
    mean(predicted == outcome) * 100
  }

  # ---- error trajectory across boosting iterations ----
  # adabag stores cumulative training error per iteration in $error
  error_trajectory <- data.frame(
    Iteration   = seq_along(final_model$error),
    Train_Error = final_model$error
  )
  # Merge CV trajectory if available (aligns on Iteration)
  if (!is.null(cv_trajectory)) {
    error_trajectory <- merge(error_trajectory, cv_trajectory, by = "Iteration")
  }

  # ---- variable importance ----
  # adabag provides $importance: a named numeric vector
  vi_raw <- final_model$importance
  vi_sorted <- sort(vi_raw[vi_raw > 0], decreasing = TRUE)

  # ---- output ----
  results <- list(
    model                = final_model,
    accuracy             = cv_accuracy,
    confusion_matrix     = conf_matrix,
    probabilities        = probs,
    variable_importance  = vi_sorted,
    error_trajectory     = error_trajectory,
    tuning_results       = tuning_results,
    cv_trajectory        = cv_trajectory,
    optimal_mfinal       = optimal_mfinal,
    metadata             = metadata,
    data_matrix          = data_matrix,
    outcome_var          = outcome_var,
    features_used        = features_to_use,
    parameters           = list(
      data_type  = data_type,
      use_merged = use_merged,
      maxdepth   = maxdepth,
      coeflearn  = coeflearn,
      cross      = cross,
      seed       = seed
    )
  )

  class(results) <- c("perform_AdaBoost", "list")

  message("AdaBoost analysis complete!")
  message(paste0("Optimal mfinal: ", optimal_mfinal,
                 "  |  Accuracy: ", round(cv_accuracy, 2), "%",
                 if (!is.null(tuning_results)) " (cross-validated)" else " (training)"))

  return(results)
}


# ---- internal null-coalescing helper (skip if already defined) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}


#' Create Plots from AdaBoost Results
#'
#' @param ada_results Output from perform_AdaBoost function
#' @param plot_type Character string specifying which plot to create:
#'   \itemize{
#'     \item \code{"trajectory"} – Error trajectory across boosting iterations
#'     \item \code{"varimp"} – Variable importance bar chart
#'     \item \code{"roc"} – ROC curve with AUC (binary outcomes only)
#'     \item \code{"confusion"} – Confusion matrix heatmap with class error rates
#'     \item \code{"tuning"} – mfinal tuning profile (requires
#'           \code{tune_mfinal = TRUE} in \code{perform_AdaBoost})
#'   }
#'   Default: \code{"trajectory"}
#' @param trajectory_type Character string controlling what is shown in the
#'        \code{"trajectory"} plot:
#'   \itemize{
#'     \item \code{"train"} – Training error only
#'     \item \code{"cv"} – CV error trajectory only (requires
#'           \code{tune_mfinal = TRUE})
#'     \item \code{"both"} – Training and CV error overlaid, with a vertical
#'           line marking the CV-optimal mfinal (requires
#'           \code{tune_mfinal = TRUE})
#'     \item \code{"train_mark"} – Training error with a vertical line marking
#'           the CV-optimal mfinal (requires \code{tune_mfinal = TRUE})
#'   }
#'   Default: \code{"both"}
#' @param top_n Integer specifying how many top features to show in the
#'        variable importance plot (default: 20)
#' @param bar_color Character string specifying bar fill color for the variable
#'        importance plot (default: \code{"#4575b4"})
#' @param line_size Numeric specifying line width for trajectory and tuning
#'        plots (default: 1)
#'
#' @return ggplot2 object
#'
#' @examples
#' ada_results <- perform_AdaBoost(df_prepped, outcome_var = "Group")
#'
#' plot_AdaBoost(ada_results, plot_type = "trajectory")
#' plot_AdaBoost(ada_results, plot_type = "trajectory", trajectory_type = "train")
#' plot_AdaBoost(ada_results, plot_type = "trajectory", trajectory_type = "cv")
#' plot_AdaBoost(ada_results, plot_type = "trajectory", trajectory_type = "train_mark")
#' plot_AdaBoost(ada_results, plot_type = "varimp", top_n = 15)
#' plot_AdaBoost(ada_results, plot_type = "roc")
#' plot_AdaBoost(ada_results, plot_type = "confusion")
#' plot_AdaBoost(ada_results, plot_type = "tuning")
#' @export
plot_AdaBoost <- function(ada_results,
                          plot_type       = "trajectory",
                          trajectory_type = "both",
                          top_n           = 20,
                          bar_color       = "#4575b4",
                          line_size       = 1) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.")

  if (!inherits(ada_results, "perform_AdaBoost"))
    stop("Input must be output from perform_AdaBoost function")

  plot_type       <- match.arg(plot_type,
                               c("trajectory", "varimp", "roc", "confusion", "tuning"))
  trajectory_type <- match.arg(trajectory_type,
                               c("both", "train", "cv", "train_mark"))

  switch(plot_type,
         trajectory = .plot_ada_trajectory(ada_results, trajectory_type, line_size),
         varimp     = .plot_ada_varimp(ada_results, top_n, bar_color),
         roc        = .plot_ada_roc(ada_results),
         confusion  = .plot_ada_confusion(ada_results),
         tuning     = .plot_ada_tuning(ada_results, line_size)
  )
}


# ---- internal: error trajectory across boosting iterations ----
.plot_ada_trajectory <- function(ada_results, trajectory_type, line_size) {

  td         <- ada_results$error_trajectory
  has_cv     <- "CV_Error" %in% colnames(td)
  opt_mfinal <- ada_results$optimal_mfinal

  # Guard: types that need CV data
  if (trajectory_type %in% c("cv", "both", "train_mark") && !has_cv)
    stop("trajectory_type = '", trajectory_type, "' requires CV error data. ",
         "Re-run perform_AdaBoost with tune_mfinal = TRUE.")

  base_subtitle <- paste0("mfinal = ", opt_mfinal,
                          "  |  depth = ", ada_results$parameters$maxdepth,
                          "  |  coeflearn = ", ada_results$parameters$coeflearn)

  # -- "train": training error only --
  if (trajectory_type == "train" || nrow(td) == 0) {
    # Guard: if Train_Error is empty (e.g. adabag version difference), show CV
    if (nrow(td) == 0 || all(is.na(td$Train_Error))) {
      message("Train_Error not available; falling back to CV error trajectory.")
      td_plot <- data.frame(Iteration = td$Iteration,
                            Error     = td$CV_Error %||% numeric(0))
      y_label <- "CV Error"
      title_  <- "AdaBoost — CV Error Trajectory (train error unavailable)"
    } else {
      td_plot <- td
      y_label <- "Training Error"
      title_  <- "AdaBoost — Training Error Trajectory"
    }
    p <- ggplot2::ggplot(td_plot,
                         ggplot2::aes(x = Iteration, y = Error)) +
      ggplot2::geom_line(color = "#d73027", linewidth = line_size) +
      ggplot2::labs(y = y_label, title = title_)

    # -- "cv": CV error only --
  } else if (trajectory_type == "cv") {
    p <- ggplot2::ggplot(td, ggplot2::aes(x = Iteration, y = CV_Error)) +
      ggplot2::geom_line(color = "#4575b4", linewidth = line_size) +
      ggplot2::geom_vline(xintercept = opt_mfinal,
                          linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::annotate("text", x = opt_mfinal, y = max(td$CV_Error),
                        label = paste0("Optimal = ", opt_mfinal),
                        hjust = -0.1, color = "red", size = 3.5) +
      ggplot2::labs(y = "CV Error",
                    title = "AdaBoost — CV Error Trajectory")

    # -- "train_mark": training error + vertical line at CV-optimal mfinal --
  } else if (trajectory_type == "train_mark") {
    p <- ggplot2::ggplot(td, ggplot2::aes(x = Iteration, y = Train_Error)) +
      ggplot2::geom_line(color = "#d73027", linewidth = line_size) +
      ggplot2::geom_vline(xintercept = opt_mfinal,
                          linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::annotate("text", x = opt_mfinal, y = max(td$Train_Error),
                        label = paste0("CV-optimal = ", opt_mfinal),
                        hjust = -0.1, color = "red", size = 3.5) +
      ggplot2::labs(y = "Training Error",
                    title = "AdaBoost — Training Error Trajectory")

    # -- "both": train + CV overlaid, vertical line at CV-optimal mfinal --
  } else {
    # Trim CV trajectory to match the length of Train_Error (final model may
    # have fewer iterations than max(mfinal) used during tuning).
    # Also guard against empty Train_Error vector.
    n_train <- nrow(td)

    if (n_train == 0 || all(is.na(td$Train_Error))) {
      message("Train_Error not available; showing CV error only.")
      p <- ggplot2::ggplot(td, ggplot2::aes(x = Iteration, y = CV_Error)) +
        ggplot2::geom_line(color = "#4575b4", linewidth = line_size) +
        ggplot2::labs(y = "CV Error",
                      title = "AdaBoost — CV Error Trajectory") +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Boosting Iteration", subtitle = base_subtitle) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5,
                                                color = "grey40", size = 9))
      return(p)
    }

    cv_vals <- if ("CV_Error" %in% colnames(td)) {
      td$CV_Error[seq_len(n_train)]
    } else {
      rep(NA_real_, n_train)
    }
    # Reshape to long for clean legend
    td_long <- rbind(
      data.frame(Iteration = td$Iteration,
                 Error     = td$Train_Error,
                 Type      = "Training"),
      data.frame(Iteration = td$Iteration,
                 Error     = cv_vals,
                 Type      = "CV (avg)")
    )
    td_long <- td_long[!is.na(td_long$Error), ]
    p <- ggplot2::ggplot(td_long,
                         ggplot2::aes(x = Iteration, y = Error,
                                      color = Type, linewidth = Type)) +
      ggplot2::geom_line() +
      ggplot2::scale_color_manual(values = c("Training" = "#d73027",
                                             "CV (avg)" = "#4575b4")) +
      ggplot2::scale_linewidth_manual(values = c("Training" = line_size,
                                                 "CV (avg)" = line_size),
                                      guide = "none") +
      ggplot2::geom_vline(xintercept = opt_mfinal,
                          linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::annotate("text", x = opt_mfinal, y = max(td_long$Error),
                        label = paste0("CV-optimal = ", opt_mfinal),
                        hjust = -0.1, color = "red", size = 3.5) +
      ggplot2::labs(y = "Error", color = NULL,
                    title = "AdaBoost — Training vs CV Error Trajectory")
  }

  p +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Boosting Iteration", subtitle = base_subtitle) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      legend.position = "top"
    )
}


# ---- internal: variable importance bar chart ----
.plot_ada_varimp <- function(ada_results, top_n, bar_color) {

  vi <- ada_results$variable_importance

  if (length(vi) == 0)
    stop("No variable importance scores available.")

  top_n <- min(top_n, length(vi))
  vi_df <- data.frame(
    Feature    = factor(names(vi)[seq_len(top_n)],
                        levels = rev(names(vi)[seq_len(top_n)])),
    Importance = as.numeric(vi[seq_len(top_n)])
  )

  ggplot2::ggplot(vi_df, ggplot2::aes(x = Feature, y = Importance)) +
    ggplot2::geom_bar(stat = "identity", fill = bar_color) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = NULL,
      y        = "Importance (weighted across boosting rounds)",
      title    = paste0("AdaBoost Variable Importance — Top ", top_n, " Features"),
      subtitle = paste0("mfinal = ", ada_results$optimal_mfinal,
                        "  |  coeflearn = ", ada_results$parameters$coeflearn)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text.y   = ggplot2::element_text(size = 9)
    )
}


# ---- internal: ROC curve (binary only) ----
.plot_ada_roc <- function(ada_results) {

  outcome <- as.factor(ada_results$metadata[[ada_results$outcome_var]])
  lvls    <- levels(outcome)

  if (length(lvls) != 2)
    stop("ROC curve is only supported for binary outcomes. ",
         "Detected classes: ", paste(lvls, collapse = ", "))

  probs     <- ada_results$probabilities
  pos_class <- lvls[2]

  if (!pos_class %in% colnames(probs))
    stop("Positive class '", pos_class, "' not found in probability matrix.")

  pos_prob  <- probs[, pos_class]
  truth_bin <- as.integer(outcome == pos_class)

  thresholds <- sort(unique(pos_prob), decreasing = TRUE)
  roc_df <- do.call(rbind, lapply(thresholds, function(thr) {
    pred <- as.integer(pos_prob >= thr)
    tp   <- sum(pred == 1 & truth_bin == 1)
    fp   <- sum(pred == 1 & truth_bin == 0)
    fn   <- sum(pred == 0 & truth_bin == 1)
    tn   <- sum(pred == 0 & truth_bin == 0)
    data.frame(Threshold = thr,
               TPR = tp / (tp + fn),
               FPR = fp / (fp + tn))
  }))
  roc_df <- rbind(data.frame(Threshold = 1, TPR = 0, FPR = 0), roc_df)

  auc_val <- round(
    sum(diff(roc_df$FPR) *
          (roc_df$TPR[-1] + roc_df$TPR[-nrow(roc_df)]) / 2), 4
  )

  ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = "#d73027", linewidth = 1.2) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50") +
    ggplot2::annotate("text", x = 0.75, y = 0.1,
                      label = paste0("AUC = ", auc_val), size = 5) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = "False Positive Rate (1 - Specificity)",
      y        = "True Positive Rate (Sensitivity)",
      title    = paste0("ROC Curve — ", ada_results$outcome_var,
                        "  (", paste(lvls, collapse = " vs "), ")"),
      subtitle = paste0("mfinal = ", ada_results$optimal_mfinal,
                        "  |  coeflearn = ", ada_results$parameters$coeflearn)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9)
    )
}


# ---- internal: confusion matrix heatmap ----
.plot_ada_confusion <- function(ada_results) {

  conf      <- ada_results$confusion_matrix
  class_err <- conf[, "class.error"]
  conf_vals <- conf[, colnames(conf) != "class.error", drop = FALSE]

  conf_df <- as.data.frame(as.table(conf_vals))
  colnames(conf_df) <- c("Predicted", "Actual", "Count")

  err_df <- data.frame(
    Predicted = rownames(conf_vals),
    Actual    = colnames(conf_vals),
    label     = paste0("err: ", round(class_err * 100, 1), "%")
  )

  ggplot2::ggplot(conf_df, ggplot2::aes(x = Actual, y = Predicted, fill = Count)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Count),
                       size = 5, color = "white", fontface = "bold") +
    ggplot2::geom_text(data = err_df,
                       ggplot2::aes(x = Actual, y = Predicted, label = label),
                       inherit.aes = FALSE, size = 3,
                       vjust = 2.2, color = "white") +
    ggplot2::scale_fill_gradient(low = "#4575b4", high = "#d73027") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title    = paste0("AdaBoost Confusion Matrix  (mfinal = ",
                        ada_results$optimal_mfinal, ")"),
      subtitle = paste0("coeflearn = ", ada_results$parameters$coeflearn,
                        "  |  Base learner depth = ",
                        ada_results$parameters$maxdepth),
      fill     = "Count"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text     = ggplot2::element_text(size = 11),
      axis.title    = ggplot2::element_text(size = 12)
    )
}


# ---- internal: mfinal tuning profile ----
.plot_ada_tuning <- function(ada_results, line_size) {

  if (is.null(ada_results$tuning_results))
    stop("No tuning results found. Re-run perform_AdaBoost with tune_mfinal = TRUE.")

  td         <- ada_results$tuning_results
  opt_mfinal <- ada_results$optimal_mfinal

  ggplot2::ggplot(td, ggplot2::aes(x = mfinal, y = cv_accuracy)) +
    ggplot2::geom_line(color = "#4575b4", linewidth = line_size) +
    ggplot2::geom_point(color = "#4575b4", size = 3) +
    ggplot2::geom_vline(xintercept = opt_mfinal,
                        linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::annotate("text",
                      x     = opt_mfinal,
                      y     = min(td$cv_accuracy),
                      label = paste0("Optimal = ", opt_mfinal),
                      hjust = -0.1, color = "red", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = td$mfinal) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = "Number of Boosting Iterations (mfinal)",
      y        = paste0("CV Accuracy (%)  [", ada_results$parameters$cross, "-fold]"),
      title    = "AdaBoost Tuning — Boosting Iterations",
      subtitle = paste0("coeflearn = ", ada_results$parameters$coeflearn,
                        "  |  Base learner depth = ",
                        ada_results$parameters$maxdepth)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
