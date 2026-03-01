#' Perform Support Vector Machine Classification on Preprocessed Metabolomics Data
#'
#' @param prepped_data Output from perform_PreprocessingPeakData function
#' @param outcome_var Character string specifying the outcome variable column name in metadata (default: "Group")
#' @param features Character vector of feature/metabolite names to include in analysis.
#'        If NULL, uses all available features (default: NULL)
#' @param data_type Character string specifying which data matrix to use: "NONPLS" or "PLS" (default: "NONPLS")
#' @param use_merged Logical indicating whether to use merged replicates if available (default: TRUE)
#' @param kernel Character string specifying the SVM kernel: "linear", "radial", "polynomial",
#'        or "sigmoid" (default: "radial")
#' @param tune_params Logical indicating whether to automatically tune cost and gamma via grid
#'        search (default: TRUE)
#' @param cost Numeric or numeric vector of cost values to search when tune_params = TRUE,
#'        or a single value to use directly (default: c(0.01, 0.1, 1, 10, 100))
#' @param gamma Numeric or numeric vector of gamma values to search when tune_params = TRUE
#'        (only for radial, polynomial, sigmoid kernels), or a single value to use directly.
#'        If NULL, defaults to 1/ncol(data) (default: NULL)
#' @param degree Integer or integer vector of degree values to search (only for polynomial kernel)
#'        when tune_params = TRUE, or a single value to use directly (default: c(2, 3, 4))
#' @param cross Integer specifying number of cross-validation folds for tuning (default: 10)
#' @param probability Logical indicating whether to compute class probabilities (default: TRUE)
#' @param seed Integer for reproducibility (default: 42)
#'
#' @return Object of class "perform_SVM" (list) containing:
#'   - model: The final svm model object
#'   - accuracy: Overall cross-validated accuracy
#'   - confusion_matrix: Confusion matrix with class error rates
#'   - tuning_results: Output from tune() if tune_params = TRUE, else NULL
#'   - optimal_params: List of optimal cost, gamma, and/or degree used
#'   - metadata: Metadata used (QC samples removed)
#'   - data_matrix: Data matrix used
#'   - outcome_var: Name of outcome variable
#'   - features_used: Character vector of features used in the model
#'   - kernel: Kernel used
#'   - parameters: List of function call parameters
#'
#' @examples
#' svm_results <- perform_SVM(df_prepped, outcome_var = "Group")
#' svm_results <- perform_SVM(df_prepped, outcome_var = "Group",
#'                             kernel = "linear", tune_params = TRUE,
#'                             cost = c(0.1, 1, 10))
#' svm_results <- perform_SVM(df_prepped, outcome_var = "Group2",
#'                             features = c("Met001", "Met002", "Met010"))
#' @export
perform_SVM <- function(prepped_data,
                        outcome_var   = "Group",
                        features      = NULL,
                        data_type     = "NONPLS",
                        use_merged    = TRUE,
                        kernel        = "radial",
                        tune_params   = TRUE,
                        cost          = c(0.01, 0.1, 1, 10, 100),
                        gamma         = NULL,
                        degree        = c(2, 3, 4),
                        cross         = 10,
                        probability   = TRUE,
                        seed          = 42) {

  # ---- package checks ----
  if (!requireNamespace("e1071", quietly = TRUE))
    stop("Package 'e1071' is required but not installed.")

  # ---- input validation ----
  if (!inherits(prepped_data, "perform_PreprocessingPeakData"))
    stop("Input must be output from perform_PreprocessingPeakData function")

  if (!data_type %in% c("NONPLS", "PLS"))
    stop("data_type must be either 'NONPLS' or 'PLS'")

  kernel <- match.arg(kernel, c("linear", "radial", "polynomial", "sigmoid"))

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
    available  <- colnames(data_matrix)
    missing_ft <- setdiff(features, available)
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

  model_data <- data.frame(outcome = outcome, data_matrix)

  # ---- default gamma ----
  default_gamma <- 1 / ncol(data_matrix)
  gamma_search  <- if (is.null(gamma)) {
    # Build a range around the default
    round(default_gamma * c(0.01, 0.1, 1, 10, 100), 8)
  } else {
    gamma
  }

  set.seed(seed)

  # ---- hyperparameter tuning ----
  tuning_results  <- NULL
  optimal_params  <- list()

  if (tune_params) {
    message("Tuning hyperparameters via ", cross, "-fold cross-validation...")

    tune_ranges <- if (kernel == "linear") {
      list(cost = cost)
    } else if (kernel == "polynomial") {
      list(cost = cost, gamma = gamma_search, degree = degree)
    } else {
      list(cost = cost, gamma = gamma_search)
    }

    tuning_results <- e1071::tune(
      e1071::svm,
      outcome ~ .,
      data       = model_data,
      kernel     = kernel,
      ranges     = tune_ranges,
      tunecontrol = e1071::tune.control(cross = cross)
    )

    best        <- tuning_results$best.parameters
    optimal_params$cost  <- best$cost
    if (kernel != "linear") optimal_params$gamma  <- best$gamma
    if (kernel == "polynomial") optimal_params$degree <- best$degree

    message(paste0("Optimal parameters: ",
                   paste(names(optimal_params),
                         round(unlist(optimal_params), 6),
                         sep = " = ", collapse = ", ")))

  } else {
    # Use single supplied values (first element if a vector was passed)
    optimal_params$cost  <- cost[1]
    if (kernel != "linear") optimal_params$gamma  <- if (length(gamma_search) == 1) gamma_search else default_gamma
    if (kernel == "polynomial") optimal_params$degree <- degree[1]
    message(paste0("Using supplied parameters: ",
                   paste(names(optimal_params),
                         round(unlist(optimal_params), 6),
                         sep = " = ", collapse = ", ")))
  }

  # ---- build final model ----
  message("Building final SVM model...")

  svm_args <- list(
    formula     = outcome ~ .,
    data        = model_data,
    kernel      = kernel,
    cost        = optimal_params$cost,
    probability = probability,
    cross       = cross   # internal CV accuracy on final model
  )
  if (kernel != "linear")      svm_args$gamma  <- optimal_params$gamma
  if (kernel == "polynomial")  svm_args$degree <- optimal_params$degree

  final_model <- do.call(e1071::svm, svm_args)

  # ---- confusion matrix ----
  predicted  <- stats::predict(final_model, model_data)
  conf_table <- table(Predicted = predicted, Actual = outcome)
  class_err  <- 1 - diag(conf_table) / colSums(conf_table)
  conf_matrix <- cbind(conf_table, class.error = round(class_err, 4))

  # cross-validated accuracy stored inside model (from cross = cross above)
  cv_accuracy <- final_model$tot.accuracy

  # ---- output ----
  results <- list(
    model           = final_model,
    accuracy        = cv_accuracy,
    confusion_matrix = conf_matrix,
    tuning_results  = tuning_results,
    optimal_params  = optimal_params,
    metadata        = metadata,
    data_matrix     = data_matrix,
    outcome_var     = outcome_var,
    features_used   = features_to_use,
    kernel          = kernel,
    parameters      = list(
      data_type   = data_type,
      use_merged  = use_merged,
      cross       = cross,
      probability = probability,
      seed        = seed
    )
  )

  class(results) <- c("perform_SVM", "list")

  message("SVM analysis complete!")
  message(paste0("Cross-validated accuracy: ", round(cv_accuracy, 2), "%"))

  return(results)
}


# ---- internal null-coalescing helper (not exported) ----
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Create Plots from SVM Results
#'
#' @param svm_results Output from perform_SVM function
#' @param plot_type Character string specifying which plot to create:
#'   \itemize{
#'     \item \code{"mds"} – MDS plot using distances derived from the SVM decision values
#'     \item \code{"confusion"} – Confusion matrix heatmap
#'     \item \code{"roc"} – ROC curve (requires probability = TRUE and binary outcome)
#'     \item \code{"tuning"} – Hyperparameter tuning error surface / profile
#'     \item \code{"boundary"} – Decision boundary plot for exactly 2 user-specified features
#'   }
#'   Default: \code{"mds"}
#' @param color_by Character string specifying metadata column to color points by
#'        (for \code{"mds"} plot, default: uses outcome_var from SVM)
#' @param point_size Numeric specifying size of points (for \code{"mds"} and \code{"boundary"} plots, default: 3)
#' @param add_labels Logical indicating whether to add sample labels (for \code{"mds"} plot, default: FALSE)
#' @param label_size Numeric specifying size of labels (for \code{"mds"} plot, default: 3)
#' @param boundary_features Character vector of exactly 2 feature names to use for the decision
#'        boundary plot (required when \code{plot_type = "boundary"}). Both features must be
#'        present in \code{svm_results$features_used}
#' @param grid_resolution Integer specifying the number of grid points per axis for the decision
#'        boundary surface (for \code{"boundary"} plot, default: 200)
#'
#' @return ggplot2 object
#'
#' @examples
#' svm_results <- perform_SVM(df_prepped, outcome_var = "Group")
#'
#' plot_SVM(svm_results, plot_type = "mds")
#' plot_SVM(svm_results, plot_type = "confusion")
#' plot_SVM(svm_results, plot_type = "roc")
#' plot_SVM(svm_results, plot_type = "tuning")
#'
#' # Decision boundary for two specific metabolites
#' plot_SVM(svm_results, plot_type = "boundary",
#'          boundary_features = c("Met001", "Met010"))
#'
#' # MDS coloured by a different metadata variable
#' plot_SVM(svm_results, plot_type = "mds", color_by = "Group2")
#' @export
plot_SVM <- function(svm_results,
                     plot_type          = "mds",
                     color_by           = NULL,
                     point_size         = 3,
                     add_labels         = FALSE,
                     label_size         = 3,
                     boundary_features  = NULL,
                     grid_resolution    = 200) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.")

  if (!inherits(svm_results, "perform_SVM"))
    stop("Input must be output from perform_SVM function")

  plot_type <- match.arg(plot_type, c("mds", "confusion", "roc", "tuning", "boundary"))

  switch(plot_type,
         mds       = .plot_svm_mds(svm_results, color_by, point_size, add_labels, label_size),
         confusion = .plot_svm_confusion(svm_results),
         roc       = .plot_svm_roc(svm_results),
         tuning    = .plot_svm_tuning(svm_results),
         boundary  = .plot_svm_boundary(svm_results, boundary_features, grid_resolution, point_size)
  )
}


# ---- internal: MDS from SVM decision values ----
.plot_svm_mds <- function(svm_results, color_by, point_size, add_labels, label_size) {

  color_by <- color_by %||% svm_results$outcome_var

  if (!color_by %in% colnames(svm_results$metadata))
    stop(paste0("Color variable '", color_by, "' not found in metadata. ",
                "Available columns: ", paste(colnames(svm_results$metadata), collapse = ", ")))

  # Use decision values as a proxy for a distance/score space
  dv <- svm_results$model$decision.values
  if (is.null(dv)) {
    # Fallback: PCA on the raw data matrix if decision values unavailable
    message("Decision values unavailable; using PCA on data matrix for MDS-like projection.")
    pca_res  <- stats::prcomp(svm_results$data_matrix, center = FALSE, scale. = FALSE)
    mds_data <- data.frame(
      Sample   = rownames(svm_results$data_matrix),
      X        = pca_res$x[, 1],
      Y        = pca_res$x[, 2],
      ColorVar = svm_results$metadata[[color_by]]
    )
    xlab <- paste0("PC1 - ", round(summary(pca_res)$importance[2, 1] * 100, 1), "%")
    ylab <- paste0("PC2 - ", round(summary(pca_res)$importance[2, 2] * 100, 1), "%")
    title <- "PCA Plot of SVM Feature Space"
  } else {
    dist_mat <- as.dist(stats::dist(dv))
    mds_res  <- cmdscale(dist_mat, eig = TRUE, x.ret = TRUE)
    mds_var  <- round(mds_res$eig / sum(abs(mds_res$eig)) * 100, 1)
    pts      <- mds_res$points
    mds_data <- data.frame(
      Sample   = rownames(pts),
      X        = pts[, 1],
      Y        = pts[, 2],
      ColorVar = svm_results$metadata[[color_by]]
    )
    xlab  <- paste0("MDS1 - ", mds_var[1], "%")
    ylab  <- paste0("MDS2 - ", mds_var[2], "%")
    title <- "MDS Plot of SVM Decision Value Space"
  }

  p <- ggplot2::ggplot(mds_data, ggplot2::aes(x = X, y = Y, color = ColorVar))

  if (add_labels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), size = label_size)
  } else {
    p <- p + ggplot2::geom_point(size = point_size)
  }

  p + ggplot2::theme_bw() +
    ggplot2::labs(x = xlab, y = ylab, color = color_by, title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
}


# ---- internal: confusion matrix heatmap ----
.plot_svm_confusion <- function(svm_results) {

  conf <- svm_results$confusion_matrix

  # Drop the class.error column for the tile plot, keep it for annotation
  class_err <- conf[, "class.error"]
  conf_vals <- conf[, colnames(conf) != "class.error", drop = FALSE]

  conf_df <- as.data.frame(as.table(conf_vals))
  colnames(conf_df) <- c("Predicted", "Actual", "Count")

  # Add class error labels on the diagonal
  n_classes <- length(class_err)
  err_df    <- data.frame(
    Predicted  = rownames(conf_vals),
    Actual     = colnames(conf_vals),
    label      = paste0("err: ", round(class_err * 100, 1), "%")
  )

  ggplot2::ggplot(conf_df, ggplot2::aes(x = Actual, y = Predicted, fill = Count)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Count), size = 5, color = "white", fontface = "bold") +
    ggplot2::geom_text(data = err_df,
                       ggplot2::aes(x = Actual, y = Predicted, label = label),
                       inherit.aes = FALSE, size = 3, vjust = 2.2, color = "white") +
    ggplot2::scale_fill_gradient(low = "#4575b4", high = "#d73027") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("SVM Confusion Matrix (", svm_results$kernel, " kernel)"),
      fill  = "Count"
    ) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text   = ggplot2::element_text(size = 11),
      axis.title  = ggplot2::element_text(size = 12)
    )
}


# ---- internal: ROC curve (binary only) ----
.plot_svm_roc <- function(svm_results) {

  if (!svm_results$parameters$probability)
    stop("ROC curve requires probability = TRUE. Re-run perform_SVM with probability = TRUE.")

  outcome  <- as.factor(svm_results$metadata[[svm_results$outcome_var]])
  lvls     <- levels(outcome)

  if (length(lvls) != 2)
    stop("ROC curve is only supported for binary outcomes. ",
         "Detected classes: ", paste(lvls, collapse = ", "))

  probs <- attr(
    stats::predict(svm_results$model,
                   data.frame(svm_results$data_matrix),
                   probability = TRUE),
    "probabilities"
  )

  pos_class <- lvls[2]
  pos_prob  <- probs[, pos_class]
  truth_bin <- as.integer(outcome == pos_class)

  # Manually compute ROC points
  thresholds <- sort(unique(pos_prob), decreasing = TRUE)
  roc_df <- do.call(rbind, lapply(thresholds, function(thr) {
    pred <- as.integer(pos_prob >= thr)
    tp   <- sum(pred == 1 & truth_bin == 1)
    fp   <- sum(pred == 1 & truth_bin == 0)
    fn   <- sum(pred == 0 & truth_bin == 1)
    tn   <- sum(pred == 0 & truth_bin == 0)
    data.frame(
      Threshold = thr,
      TPR = tp / (tp + fn),    # sensitivity
      FPR = fp / (fp + tn)     # 1 - specificity
    )
  }))
  roc_df <- rbind(data.frame(Threshold = 1, TPR = 0, FPR = 0), roc_df)

  # Trapezoidal AUC
  auc_val <- round(
    sum(diff(roc_df$FPR) * (roc_df$TPR[-1] + roc_df$TPR[-nrow(roc_df)]) / 2),
    4
  )

  ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line(color = "#d73027", linewidth = 1.2) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50") +
    ggplot2::annotate("text", x = 0.75, y = 0.1,
                      label = paste0("AUC = ", auc_val), size = 5) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x     = "False Positive Rate (1 - Specificity)",
      y     = "True Positive Rate (Sensitivity)",
      title = paste0("ROC Curve — ", svm_results$outcome_var,
                     " (", paste(lvls, collapse = " vs "), ")")
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
}


# ---- internal: decision boundary (2 user-specified features) ----
.plot_svm_boundary <- function(svm_results, boundary_features, grid_resolution, point_size) {

  # ---- validate boundary_features ----
  if (is.null(boundary_features) || length(boundary_features) == 0) {
    message("boundary_features not specified. Automatically selecting the ",
            "top 2 features by SVM weight / importance.")
    boundary_features <- tryCatch({
      w   <- t(svm_results$model$coefs) %*% svm_results$model$SV
      imp <- abs(colMeans(w))
      # SV column names may be sanitised — map back to features_used
      feat_sanitised <- make.names(svm_results$features_used, unique = TRUE)
      top2_sanitised <- names(sort(imp, decreasing = TRUE))[seq_len(2)]
      # Match sanitised names back to original feature names
      svm_results$features_used[match(top2_sanitised, feat_sanitised)]
    }, error = function(e) {
      head(svm_results$features_used, 2)
    })
    # Remove any NAs from failed match
    boundary_features <- boundary_features[!is.na(boundary_features)]
    if (length(boundary_features) < 2)
      boundary_features <- head(svm_results$features_used, 2)
    message(paste0("Using features: ",
                   paste(boundary_features, collapse = ", ")))
  }

  if (length(boundary_features) != 2)
    stop("Exactly 2 features must be supplied to boundary_features. ",
         "You supplied ", length(boundary_features), ".")

  missing_ft <- setdiff(boundary_features, svm_results$features_used)
  if (length(missing_ft) > 0)
    stop("The following features were not used in the SVM model: ",
         paste(missing_ft, collapse = ", "), ". ",
         "Choose from: ", paste(head(svm_results$features_used, 10), collapse = ", "),
         if (length(svm_results$features_used) > 10) " ..." else "")

  f1 <- boundary_features[1]
  f2 <- boundary_features[2]

  # ---- refit a 2-feature SVM with the same optimal params ----
  outcome    <- as.factor(svm_results$metadata[[svm_results$outcome_var]])
  sub_data   <- data.frame(
    outcome  = outcome,
    x1       = svm_results$data_matrix[[f1]],
    x2       = svm_results$data_matrix[[f2]]
  )

  svm_args <- list(
    formula     = outcome ~ x1 + x2,
    data        = sub_data,
    kernel      = svm_results$kernel,
    cost        = svm_results$optimal_params$cost,
    probability = FALSE
  )
  if (svm_results$kernel != "linear")
    svm_args$gamma  <- svm_results$optimal_params$gamma
  if (svm_results$kernel == "polynomial")
    svm_args$degree <- svm_results$optimal_params$degree

  model_2d <- do.call(e1071::svm, svm_args)

  # ---- build prediction grid ----
  x1_seq <- seq(min(sub_data$x1), max(sub_data$x1), length.out = grid_resolution)
  x2_seq <- seq(min(sub_data$x2), max(sub_data$x2), length.out = grid_resolution)
  grid   <- expand.grid(x1 = x1_seq, x2 = x2_seq)

  grid$predicted <- stats::predict(model_2d, grid)

  # ---- identify support vectors ----
  sv_idx  <- svm_results$model$index   # indices from original full model
  # keep only those that are within the row range of sub_data
  sv_idx  <- sv_idx[sv_idx <= nrow(sub_data)]
  sv_data <- sub_data[sv_idx, , drop = FALSE]

  # ---- plot ----
  ggplot2::ggplot() +
    # decision region background
    ggplot2::geom_tile(data  = grid,
                       ggplot2::aes(x = x1, y = x2, fill = predicted),
                       alpha = 0.25) +
    # all samples
    ggplot2::geom_point(data  = sub_data,
                        ggplot2::aes(x = x1, y = x2, colour = outcome),
                        size  = point_size) +
    # support vectors highlighted with a larger open circle
    ggplot2::geom_point(data  = sv_data,
                        ggplot2::aes(x = x1, y = x2),
                        shape = 21, size = point_size + 2,
                        colour = "black", fill = NA, stroke = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Pastel1", guide = "none") +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x      = f1,
      y      = f2,
      colour = svm_results$outcome_var,
      title  = paste0("SVM Decision Boundary (", svm_results$kernel, " kernel)"),
      subtitle = paste0("Support vectors shown as open circles  |  ",
                        "Grid resolution: ", grid_resolution, " \u00d7 ", grid_resolution)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40", size = 9)
    )
}


# ---- internal: tuning error surface / profile ----
.plot_svm_tuning <- function(svm_results) {

  if (is.null(svm_results$tuning_results))
    stop("No tuning results found. Re-run perform_SVM with tune_params = TRUE.")

  tune_perf <- svm_results$tuning_results$performances
  kernel    <- svm_results$kernel

  # Linear kernel: 1D profile over cost
  if (kernel == "linear") {
    p <- ggplot2::ggplot(tune_perf,
                         ggplot2::aes(x = cost, y = error)) +
      ggplot2::geom_line(color = "#4575b4", linewidth = 1) +
      ggplot2::geom_point(color = "#4575b4", size = 3) +
      ggplot2::geom_vline(xintercept = svm_results$optimal_params$cost,
                          linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::scale_x_log10() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Cost (log scale)",
        y     = "Cross-validation Error",
        title = "SVM Tuning — Cost (linear kernel)"
      ) +
      ggplot2::annotate("text",
                        x = svm_results$optimal_params$cost,
                        y = max(tune_perf$error),
                        label = paste0("Optimal C = ", svm_results$optimal_params$cost),
                        hjust = -0.1, color = "red", size = 3.5) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  } else {
    # 2D heat-map over cost × gamma (and facet by degree for polynomial)
    base_plot <- ggplot2::ggplot(
      tune_perf,
      ggplot2::aes(x = factor(cost), y = factor(gamma), fill = error)
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "#4575b4", high = "#d73027",
                                   name = "CV Error") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x     = "Cost (C)",
        y     = "Gamma",
        title = paste0("SVM Tuning — Cost × Gamma (", kernel, " kernel)")
      ) +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    p <- if (kernel == "polynomial" && "degree" %in% colnames(tune_perf)) {
      base_plot + ggplot2::facet_wrap(~ paste0("degree = ", degree))
    } else {
      base_plot
    }
  }

  p
}
