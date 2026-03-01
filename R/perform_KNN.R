#' Perform K-Nearest Neighbors Classification on Preprocessed Metabolomics Data
#'
#' Uses KNN via the \code{caret} package with built-in cross-validation for
#' K tuning. Distance metric is controlled via the \code{distance} parameter
#' using the \code{kknn} method inside caret, which supports any Minkowski
#' distance. K can be tuned automatically via caret's \code{trainControl} or
#' set manually.
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
#' @param tune_k Logical indicating whether to automatically tune K via
#'        cross-validation through caret (default: TRUE)
#' @param k Integer or integer vector of K values to search when
#'        \code{tune_k = TRUE}, or a single value to use directly when
#'        \code{tune_k = FALSE} (default: \code{seq(1, 21, by = 2)})
#' @param distance Numeric specifying the Minkowski distance parameter:
#'        1 = Manhattan, 2 = Euclidean, higher values approach Chebyshev
#'        (default: 2)
#' @param kernel Character string specifying the kernel for weighted voting
#'        passed to \code{kknn}. One of \code{"rectangular"} (standard
#'        unweighted KNN), \code{"triangular"}, \code{"epanechnikov"},
#'        \code{"biweight"}, \code{"triweight"}, \code{"cos"}, \code{"inv"},
#'        \code{"gaussian"}, or \code{"optimal"} (default: \code{"optimal"})
#' @param cv_method Character string specifying caret's resampling method:
#'        \code{"cv"} (k-fold), \code{"LOOCV"}, or \code{"repeatedcv"}
#'        (default: \code{"cv"})
#' @param cross Integer specifying number of CV folds (used when
#'        \code{cv_method = "cv"} or \code{"repeatedcv"}, default: 10)
#' @param repeats Integer specifying number of repeats (used only when
#'        \code{cv_method = "repeatedcv"}, default: 3)
#' @param seed Integer for reproducibility (default: 42)
#'
#' @return Object of class "perform_KNN" (list) containing:
#'   \itemize{
#'     \item \code{model} – Final caret train object
#'     \item \code{accuracy} – Cross-validated accuracy at optimal K
#'     \item \code{confusion_matrix} – Confusion matrix with class error rates
#'     \item \code{probabilities} – Matrix of class probabilities per sample
#'     \item \code{tuning_results} – Data frame of K vs CV accuracy/Kappa
#'           (NULL if \code{tune_k = FALSE})
#'     \item \code{optimal_k} – Optimal or specified K value
#'     \item \code{metadata} – Metadata used (QC samples removed)
#'     \item \code{data_matrix} – Data matrix used
#'     \item \code{outcome_var} – Name of outcome variable
#'     \item \code{features_used} – Character vector of features in model
#'     \item \code{parameters} – List of function call parameters
#'   }
#'
#' @examples
#' knn_results <- perform_KNN(df_prepped, outcome_var = "Group")
#'
#' # Manhattan distance, LOOCV, no tuning
#' knn_results <- perform_KNN(df_prepped, outcome_var = "Group",
#'                             tune_k    = FALSE,
#'                             k         = 7,
#'                             distance  = 1,
#'                             cv_method = "LOOCV")
#'
#' # Repeated CV with wider K range
#' knn_results <- perform_KNN(df_prepped,
#'                             k         = seq(1, 31, by = 2),
#'                             cv_method = "repeatedcv",
#'                             cross     = 10,
#'                             repeats   = 5)
#' @export
perform_KNN <- function(prepped_data,
                        outcome_var = "Group",
                        features    = NULL,
                        data_type   = "NONPLS",
                        use_merged  = TRUE,
                        tune_k      = TRUE,
                        k           = seq(1, 21, by = 2),
                        distance    = 2,
                        kernel      = "optimal",
                        cv_method   = "cv",
                        cross       = 10,
                        repeats     = 3,
                        seed        = 42) {

  # ---- package checks ----
  for (pkg in c("caret", "kknn")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(paste0("Package '", pkg, "' is required but not installed."))
  }

  # ---- input validation ----
  if (!inherits(prepped_data, "perform_PreprocessingPeakData"))
    stop("Input must be output from perform_PreprocessingPeakData function")

  if (!data_type %in% c("NONPLS", "PLS"))
    stop("data_type must be either 'NONPLS' or 'PLS'")

  valid_kernels <- c("rectangular", "triangular", "epanechnikov", "biweight",
                     "triweight", "cos", "inv", "gaussian", "optimal")
  kernel    <- match.arg(kernel, valid_kernels)
  cv_method <- match.arg(cv_method, c("cv", "LOOCV", "repeatedcv"))

  if (any(k < 1) || any(k != floor(k)))
    stop("All k values must be positive integers.")

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

  set.seed(seed)

  # ---- caret trainControl ----
  cv_label <- switch(cv_method,
                     cv         = paste0(cross, "-fold CV"),
                     LOOCV      = "LOOCV",
                     repeatedcv = paste0(cross, "-fold CV × ", repeats, " repeats")
  )

  ctrl_args <- list(
    method          = cv_method,
    classProbs      = TRUE,
    summaryFunction = caret::multiClassSummary,
    savePredictions = "final"
  )
  if (cv_method %in% c("cv", "repeatedcv")) ctrl_args$number  <- cross
  if (cv_method == "repeatedcv")            ctrl_args$repeats <- repeats

  train_ctrl <- do.call(caret::trainControl, ctrl_args)

  # ---- tune grid: kknn uses kmax, distance, kernel ----
  # caret's kknn method tunes over kmax (= K), distance, kernel
  if (tune_k) {
    message(paste0("Tuning K via caret (", cv_label, ")..."))
    tune_grid <- expand.grid(
      kmax     = k,
      distance = distance,   # fixed distance, user-specified
      kernel   = kernel      # fixed kernel,   user-specified
    )
  } else {
    tune_grid <- expand.grid(
      kmax     = k[1],
      distance = distance,
      kernel   = kernel
    )
    message(paste0("Using specified K = ", k[1],
                   ", distance = ", distance,
                   ", kernel = ", kernel))
  }

  # ---- train model ----
  message("Training KNN model via caret...")
  final_model <- caret::train(
    x         = data_matrix,
    y         = outcome,
    method    = "kknn",
    trControl = train_ctrl,
    tuneGrid  = tune_grid,
    metric    = "Accuracy"
  )

  optimal_k <- final_model$bestTune$kmax

  # ---- tuning results ----
  tuning_results <- NULL
  if (tune_k && nrow(final_model$results) > 1) {
    tuning_results <- final_model$results[, c("kmax", "Accuracy", "Kappa")]
    colnames(tuning_results) <- c("k", "cv_accuracy", "cv_kappa")
    tuning_results$cv_accuracy <- tuning_results$cv_accuracy * 100
  }

  cv_accuracy <- max(final_model$results$Accuracy) * 100

  # ---- confusion matrix ----
  predicted  <- stats::predict(final_model, data_matrix)
  conf_table <- table(Predicted = predicted, Actual = outcome)
  class_err  <- 1 - diag(conf_table) / colSums(conf_table)
  conf_matrix <- cbind(conf_table, class.error = round(class_err, 4))

  # ---- class probabilities ----
  probs <- stats::predict(final_model, data_matrix, type = "prob")

  # ---- output ----
  results <- list(
    model            = final_model,
    accuracy         = cv_accuracy,
    confusion_matrix = conf_matrix,
    probabilities    = as.matrix(probs),
    tuning_results   = tuning_results,
    optimal_k        = optimal_k,
    metadata         = metadata,
    data_matrix      = data_matrix,
    outcome_var      = outcome_var,
    features_used    = features_to_use,
    parameters       = list(
      data_type  = data_type,
      use_merged = use_merged,
      distance   = distance,
      kernel     = kernel,
      cv_method  = cv_method,
      cross      = cross,
      repeats    = repeats,
      seed       = seed
    )
  )

  class(results) <- c("perform_KNN", "list")

  message("KNN analysis complete!")
  message(paste0("Optimal K: ", optimal_k,
                 "  |  CV accuracy: ", round(cv_accuracy, 2), "%",
                 "  (", cv_label, ")"))

  return(results)
}


# ---- internal null-coalescing helper (skip if already defined) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}


#' Create Plots from KNN Results
#'
#' @param knn_results Output from perform_KNN function
#' @param plot_type Character string specifying which plot to create:
#'   \itemize{
#'     \item \code{"tuning"} – K tuning profile (CV accuracy vs K)
#'     \item \code{"confusion"} – Confusion matrix heatmap with class error rates
#'     \item \code{"roc"} – ROC curve with AUC (binary outcomes only)
#'     \item \code{"boundary"} – Decision boundary for exactly 2 user-specified features
#'   }
#'   Default: \code{"tuning"}
#' @param point_size Numeric specifying size of points for \code{"boundary"}
#'        plot (default: 3)
#' @param boundary_features Character vector of exactly 2 feature names to use
#'        for the decision boundary plot (required when
#'        \code{plot_type = "boundary"})
#' @param grid_resolution Integer specifying the number of grid points per axis
#'        for the decision boundary surface (default: 150)
#'
#' @return ggplot2 object
#'
#' @examples
#' knn_results <- perform_KNN(df_prepped, outcome_var = "Group")
#'
#' plot_KNN(knn_results, plot_type = "tuning")
#' plot_KNN(knn_results, plot_type = "confusion")
#' plot_KNN(knn_results, plot_type = "roc")
#' plot_KNN(knn_results, plot_type = "boundary",
#'          boundary_features = c("Met001", "Met010"))
#' @export
plot_KNN <- function(knn_results,
                     plot_type         = "tuning",
                     point_size        = 3,
                     boundary_features = NULL,
                     grid_resolution   = 150) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.")

  if (!inherits(knn_results, "perform_KNN"))
    stop("Input must be output from perform_KNN function")

  plot_type <- match.arg(plot_type, c("tuning", "confusion", "roc", "boundary"))

  switch(plot_type,
         tuning    = .plot_knn_tuning(knn_results),
         confusion = .plot_knn_confusion(knn_results),
         roc       = .plot_knn_roc(knn_results),
         boundary  = .plot_knn_boundary(knn_results, boundary_features,
                                        grid_resolution, point_size)
  )
}


# ---- internal: K tuning profile ----
.plot_knn_tuning <- function(knn_results) {

  if (is.null(knn_results$tuning_results))
    stop("No tuning results found. Re-run perform_KNN with tune_k = TRUE ",
         "and more than one K value.")

  td    <- knn_results$tuning_results
  opt_k <- knn_results$optimal_k

  dist_label <- switch(as.character(knn_results$parameters$distance),
                       "1" = "Manhattan",
                       "2" = "Euclidean",
                       paste0("Minkowski (p=", knn_results$parameters$distance, ")"))

  ggplot2::ggplot(td, ggplot2::aes(x = k, y = cv_accuracy)) +
    ggplot2::geom_line(color = "#4575b4", linewidth = 1) +
    ggplot2::geom_point(color = "#4575b4", size = 3) +
    ggplot2::geom_vline(xintercept = opt_k,
                        linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::annotate("text",
                      x     = opt_k,
                      y     = min(td$cv_accuracy),
                      label = paste0("Optimal K = ", opt_k),
                      hjust = -0.1, color = "red", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = td$k) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = "K (Number of Neighbors)",
      y        = "CV Accuracy (%)",
      title    = "KNN Tuning — K Parameter",
      subtitle = paste0("Distance: ", dist_label,
                        "  |  Kernel: ", knn_results$parameters$kernel,
                        "  |  Method: ", knn_results$parameters$cv_method)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1)
    )
}


# ---- internal: confusion matrix heatmap ----
.plot_knn_confusion <- function(knn_results) {

  conf      <- knn_results$confusion_matrix
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
      title = paste0("KNN Confusion Matrix  (K = ", knn_results$optimal_k, ")"),
      fill  = "Count"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text  = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12)
    )
}


# ---- internal: ROC curve (binary only) ----
.plot_knn_roc <- function(knn_results) {

  outcome <- as.factor(knn_results$metadata[[knn_results$outcome_var]])
  lvls    <- levels(outcome)

  if (length(lvls) != 2)
    stop("ROC curve is only supported for binary outcomes. ",
         "Detected classes: ", paste(lvls, collapse = ", "))

  probs     <- knn_results$probabilities
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

  dist_label <- switch(as.character(knn_results$parameters$distance),
                       "1" = "Manhattan",
                       "2" = "Euclidean",
                       paste0("Minkowski (p=", knn_results$parameters$distance, ")"))

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
      title    = paste0("ROC Curve — ", knn_results$outcome_var,
                        "  (", paste(lvls, collapse = " vs "), ")"),
      subtitle = paste0("K = ", knn_results$optimal_k,
                        "  |  Distance: ", dist_label,
                        "  |  Kernel: ", knn_results$parameters$kernel)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9)
    )
}


# ---- internal: decision boundary (2 user-specified features) ----
.plot_knn_boundary <- function(knn_results, boundary_features,
                               grid_resolution, point_size) {

  if (is.null(boundary_features) || length(boundary_features) == 0) {
    message("boundary_features not specified. Automatically selecting the ",
            "top 2 features by caret variable importance.")
    vi <- tryCatch(
      caret::varImp(knn_results$model)$importance,
      error = function(e) NULL
    )
    boundary_features <- if (!is.null(vi)) {
      vi$Overall <- rowMeans(vi)
      rownames(vi)[order(vi$Overall, decreasing = TRUE)][seq_len(2)]
    } else {
      head(knn_results$features_used, 2)
    }
    message(paste0("Using features: ", paste(boundary_features, collapse = ", ")))
  }

  if (length(boundary_features) != 2)
    stop("Exactly 2 features must be supplied to boundary_features. ",
         "You supplied ", length(boundary_features), ".")

  missing_ft <- setdiff(boundary_features, knn_results$features_used)
  if (length(missing_ft) > 0)
    stop("The following features were not used in the KNN model: ",
         paste(missing_ft, collapse = ", "), ". ",
         "Choose from: ", paste(head(knn_results$features_used, 10), collapse = ", "),
         if (length(knn_results$features_used) > 10) " ..." else "")

  f1      <- boundary_features[1]
  f2      <- boundary_features[2]
  outcome <- as.factor(knn_results$metadata[[knn_results$outcome_var]])

  sub_x <- data.frame(
    x1 = knn_results$data_matrix[[f1]],
    x2 = knn_results$data_matrix[[f2]]
  )

  # Refit caret KNN on 2 features using optimal K + same parameters
  ctrl_2d <- caret::trainControl(method = "none", classProbs = TRUE)
  model_2d <- caret::train(
    x        = sub_x,
    y        = outcome,
    method   = "kknn",
    trControl = ctrl_2d,
    tuneGrid  = expand.grid(
      kmax     = knn_results$optimal_k,
      distance = knn_results$parameters$distance,
      kernel   = knn_results$parameters$kernel
    )
  )

  # Prediction grid
  x1_seq <- seq(min(sub_x$x1), max(sub_x$x1), length.out = grid_resolution)
  x2_seq <- seq(min(sub_x$x2), max(sub_x$x2), length.out = grid_resolution)
  grid   <- expand.grid(x1 = x1_seq, x2 = x2_seq)
  grid$predicted <- stats::predict(model_2d, grid)

  dist_label <- switch(as.character(knn_results$parameters$distance),
                       "1" = "Manhattan",
                       "2" = "Euclidean",
                       paste0("Minkowski p=", knn_results$parameters$distance))

  sub_x$outcome <- outcome

  ggplot2::ggplot() +
    ggplot2::geom_tile(data = grid,
                       ggplot2::aes(x = x1, y = x2, fill = predicted),
                       alpha = 0.25) +
    ggplot2::geom_point(data = sub_x,
                        ggplot2::aes(x = x1, y = x2, colour = outcome),
                        size = point_size) +
    ggplot2::scale_fill_brewer(palette = "Pastel1", guide = "none") +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = f1,
      y        = f2,
      colour   = knn_results$outcome_var,
      title    = paste0("KNN Decision Boundary  (K = ", knn_results$optimal_k, ")"),
      subtitle = paste0("Distance: ", dist_label,
                        "  |  Kernel: ", knn_results$parameters$kernel,
                        "  |  Grid: ", grid_resolution, " \u00d7 ", grid_resolution)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9)
    )
}
