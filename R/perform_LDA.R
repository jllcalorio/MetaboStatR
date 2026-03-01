#' Perform Linear Discriminant Analysis on Preprocessed Metabolomics Data
#'
#' Uses \code{MASS::lda()} via \code{caret::train()} for unified cross-validation.
#' Because metabolomics data is typically high-dimensional (p >> n), the function
#' automatically applies PCA preprocessing prior to LDA (controlled by
#' \code{pca_preprocess}). Prior probabilities can be selected automatically
#' by comparing uniform vs proportional via caret CV, or set manually.
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
#' @param prior Either \code{"auto"} (select between uniform and proportional
#'        via caret CV), \code{"uniform"} (equal priors for all classes),
#'        \code{"proportional"} (class frequencies in training data), or a
#'        named numeric vector of manual prior probabilities that must sum to 1
#'        and whose names match the class levels (default: \code{"auto"})
#' @param pca_preprocess Logical indicating whether to reduce dimensions with
#'        PCA before LDA via caret's \code{preProcess}. Strongly recommended
#'        when features exceed samples (default: TRUE)
#' @param pca_variance_threshold Numeric in (0, 1] specifying the cumulative
#'        proportion of variance to retain when \code{pca_preprocess = TRUE}.
#'        Passed to caret as \code{thresh} in \code{preProcess} (default: 0.95)
#' @param cv_method Character string specifying caret's resampling method:
#'        \code{"cv"} (k-fold), \code{"LOOCV"}, or \code{"repeatedcv"}
#'        (default: \code{"cv"})
#' @param cross Integer specifying number of CV folds (used when
#'        \code{cv_method = "cv"} or \code{"repeatedcv"}, default: 10)
#' @param repeats Integer specifying number of repeats (used only when
#'        \code{cv_method = "repeatedcv"}, default: 3)
#' @param tol Numeric tolerance passed to \code{MASS::lda()} to handle
#'        near-singular covariance matrices (default: 1e-4)
#' @param seed Integer for reproducibility (default: 42)
#'
#' @return Object of class "perform_LDA" (list) containing:
#'   \itemize{
#'     \item \code{model} – Final caret train object (method = "lda")
#'     \item \code{accuracy} – Cross-validated accuracy from caret
#'     \item \code{confusion_matrix} – Confusion matrix with class error rates
#'     \item \code{probabilities} – Matrix of posterior class probabilities
#'     \item \code{ld_scores} – Data frame of sample scores on each LD axis
#'     \item \code{ld_coefficients} – Data frame of standardised LD coefficients
#'     \item \code{variance_explained} – Named numeric vector of proportion of
#'           between-group variance explained by each LD axis
#'     \item \code{prior_used} – Named numeric vector of priors used
#'     \item \code{pca_model} – caret preProcess object (NULL if
#'           \code{pca_preprocess = FALSE})
#'     \item \code{tuning_results} – Data frame comparing CV accuracy of
#'           uniform vs proportional prior (NULL if prior is not "auto")
#'     \item \code{metadata} – Metadata used (QC samples removed)
#'     \item \code{data_matrix} – Data matrix used
#'     \item \code{outcome_var} – Name of outcome variable
#'     \item \code{features_used} – Character vector of features in model
#'     \item \code{parameters} – List of function call parameters
#'   }
#'
#' @examples
#' lda_results <- perform_LDA(df_prepped, outcome_var = "Group")
#'
#' # LOOCV, manual uniform prior
#' lda_results <- perform_LDA(df_prepped, prior = "uniform",
#'                             cv_method = "LOOCV")
#'
#' # Manual prior probabilities
#' lda_results <- perform_LDA(df_prepped,
#'                             prior = c(Start = 0.4, End = 0.6))
#'
#' # Repeated CV, retain 90% PCA variance
#' lda_results <- perform_LDA(df_prepped,
#'                             pca_variance_threshold = 0.90,
#'                             cv_method = "repeatedcv",
#'                             cross = 10, repeats = 5)
#' @export
perform_LDA <- function(prepped_data,
                        outcome_var            = "Group",
                        features               = NULL,
                        data_type              = "NONPLS",
                        use_merged             = TRUE,
                        prior                  = "auto",
                        pca_preprocess         = TRUE,
                        pca_variance_threshold = 0.95,
                        cv_method              = "cv",
                        cross                  = 10,
                        repeats                = 3,
                        tol                    = 1e-4,
                        seed                   = 42) {

  # ---- package checks ----
  for (pkg in c("caret", "MASS")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(paste0("Package '", pkg, "' is required but not installed."))
  }

  # ---- input validation ----
  if (!inherits(prepped_data, "perform_PreprocessingPeakData"))
    stop("Input must be output from perform_PreprocessingPeakData function")

  if (!data_type %in% c("NONPLS", "PLS"))
    stop("data_type must be either 'NONPLS' or 'PLS'")

  cv_method <- match.arg(cv_method, c("cv", "LOOCV", "repeatedcv"))

  prior_is_auto   <- identical(prior, "auto")
  prior_is_string <- is.character(prior) && length(prior) == 1
  prior_is_vector <- is.numeric(prior) && !is.null(names(prior))

  if (!prior_is_auto && prior_is_string &&
      !prior %in% c("uniform", "proportional"))
    stop("prior must be 'auto', 'uniform', 'proportional', or a named numeric vector.")

  if (prior_is_vector && abs(sum(prior) - 1) > 1e-6)
    stop("Manual prior probabilities must sum to 1.")

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

  lvls <- levels(outcome)
  set.seed(seed)

  # ---- PCA preprocessing via caret::preProcess ----
  pca_model  <- NULL
  input_data <- as.data.frame(data_matrix)

  if (pca_preprocess) {
    pca_model  <- caret::preProcess(input_data,
                                    method = "pca",
                                    thresh  = pca_variance_threshold)
    input_data <- stats::predict(pca_model, input_data)
    n_pcs      <- ncol(input_data)
    message(paste0("PCA preprocessing: retaining ", n_pcs,
                   " PCs (>= ", pca_variance_threshold * 100, "% variance)"))
  }

  # ---- prior helper ----
  .build_prior <- function(type, outcome_vec) {
    if (type == "uniform") {
      setNames(rep(1 / length(lvls), length(lvls)), lvls)
    } else {
      tbl <- table(outcome_vec)
      setNames(as.numeric(tbl) / sum(tbl), names(tbl))
    }
  }

  # ---- caret trainControl ----
  .make_ctrl <- function() {
    args <- list(
      method          = cv_method,
      classProbs      = TRUE,
      summaryFunction = caret::multiClassSummary,
      savePredictions = "final"
    )
    if (cv_method %in% c("cv", "repeatedcv")) args$number  <- cross
    if (cv_method == "repeatedcv")            args$repeats <- repeats
    do.call(caret::trainControl, args)
  }

  cv_label <- switch(cv_method,
                     cv         = paste0(cross, "-fold CV"),
                     LOOCV      = "LOOCV",
                     repeatedcv = paste0(cross, "-fold CV \u00d7 ", repeats, " repeats")
  )

  # ---- auto prior selection via caret CV ----
  tuning_results <- NULL
  prior_used     <- NULL

  if (prior_is_auto) {
    message(paste0("Auto-selecting prior via caret (", cv_label, ")..."))

    cv_acc_by_prior <- vapply(c("uniform", "proportional"), function(pt) {
      pr <- .build_prior(pt, outcome)
      # Pass prior via modelLookup-compatible model wrapper
      fit <- tryCatch(
        caret::train(
          x         = input_data,
          y         = outcome,
          method    = "lda",
          trControl = .make_ctrl(),
          metric    = "Accuracy",
          prior     = pr,
          tol       = tol
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      max(fit$results$Accuracy, na.rm = TRUE)
    }, numeric(1))

    tuning_results <- data.frame(
      prior       = c("uniform", "proportional"),
      cv_accuracy = cv_acc_by_prior * 100
    )

    best_prior_type <- c("uniform", "proportional")[which.max(cv_acc_by_prior)]
    if (length(best_prior_type) == 0 || is.na(best_prior_type)) {
      warning("Prior CV selection failed; defaulting to 'proportional'.")
      best_prior_type <- "proportional"
    }
    prior_used <- .build_prior(best_prior_type, outcome)
    message(paste0("Selected prior: '", best_prior_type,
                   "'  (CV accuracy: ",
                   round(max(cv_acc_by_prior, na.rm = TRUE) * 100, 2), "%)"))

  } else if (prior_is_string) {
    prior_used <- .build_prior(prior, outcome)
    message(paste0("Using '", prior, "' prior"))

  } else {
    missing_lvls <- setdiff(lvls, names(prior))
    if (length(missing_lvls) > 0)
      stop("Manual prior is missing entries for classes: ",
           paste(missing_lvls, collapse = ", "))
    prior_used <- prior[lvls]
    message("Using manual prior probabilities")
  }

  # ---- build final caret model with selected prior ----
  message(paste0("Training final LDA model via caret (", cv_label, ")..."))
  final_model <- caret::train(
    x         = input_data,
    y         = outcome,
    method    = "lda",
    trControl = .make_ctrl(),
    metric    = "Accuracy",
    prior     = prior_used,
    tol       = tol
  )

  cv_accuracy <- max(final_model$results$Accuracy, na.rm = TRUE) * 100

  # ---- extract underlying MASS::lda object for scores/loadings ----
  lda_fit <- final_model$finalModel

  # ---- predictions and probabilities on full data ----
  probs_df   <- stats::predict(final_model, input_data, type = "prob")
  predicted  <- stats::predict(final_model, input_data)
  probs      <- as.matrix(probs_df)

  # ---- confusion matrix ----
  conf_table  <- table(Predicted = predicted, Actual = outcome)
  class_err   <- 1 - diag(conf_table) / colSums(conf_table)
  conf_matrix <- cbind(conf_table, class.error = round(class_err, 4))

  # ---- LD scores ----
  ld_scores_mat <- as.matrix(input_data) %*% lda_fit$scaling
  colnames(ld_scores_mat) <- paste0("LD", seq_len(ncol(lda_fit$scaling)))
  ld_scores_df  <- as.data.frame(ld_scores_mat)
  ld_scores_df  <- cbind(
    Sample    = rownames(data_matrix),
    ld_scores_df,
    Predicted = predicted,
    Actual    = outcome
  )

  # ---- variance explained ----
  svd_vals      <- lda_fit$svd
  var_explained <- setNames(svd_vals^2 / sum(svd_vals^2),
                            paste0("LD", seq_along(svd_vals)))

  # ---- LD coefficients: back-project PCA → original feature space ----
  coef_raw <- lda_fit$scaling   # PCs × LDs (or features × LDs if no PCA)

  if (pca_preprocess) {
    rot       <- pca_model$rotation[, rownames(coef_raw), drop = FALSE]
    coef_orig <- rot %*% coef_raw
    rownames(coef_orig) <- features_to_use
  } else {
    coef_orig <- coef_raw
  }

  ld_sd     <- apply(ld_scores_mat, 2, stats::sd)
  coef_std  <- sweep(coef_orig, 2, ld_sd, "/")
  ld_coef_df <- as.data.frame(coef_std)
  ld_coef_df$Feature <- rownames(coef_std)

  # ---- output ----
  results <- list(
    model              = final_model,
    accuracy           = cv_accuracy,
    confusion_matrix   = conf_matrix,
    probabilities      = probs,
    ld_scores          = ld_scores_df,
    ld_coefficients    = ld_coef_df,
    variance_explained = var_explained,
    prior_used         = prior_used,
    pca_model          = pca_model,
    tuning_results     = tuning_results,
    metadata           = metadata,
    data_matrix        = data_matrix,
    outcome_var        = outcome_var,
    features_used      = features_to_use,
    parameters         = list(
      data_type              = data_type,
      use_merged             = use_merged,
      pca_preprocess         = pca_preprocess,
      pca_variance_threshold = pca_variance_threshold,
      cv_method              = cv_method,
      cross                  = cross,
      repeats                = repeats,
      tol                    = tol,
      seed                   = seed
    )
  )

  class(results) <- c("perform_LDA", "list")

  message("LDA analysis complete!")
  message(paste0("CV accuracy: ", round(cv_accuracy, 2), "%",
                 "  (", cv_label, ")",
                 "  |  LD axes: ", length(var_explained)))

  return(results)
}


# ---- internal null-coalescing helper (skip if already defined) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}


#' Create Plots from LDA Results
#'
#' @param lda_results Output from perform_LDA function
#' @param plot_type Character string specifying which plot to create:
#'   \itemize{
#'     \item \code{"scores"} – Sample scores in discriminant space (LD1 vs LD2
#'           scatter with ellipses, or LD1 density if only one LD axis exists)
#'     \item \code{"coefficients"} – Top feature loadings on a chosen LD axis
#'     \item \code{"roc"} – ROC curve with AUC (binary outcomes only)
#'     \item \code{"confusion"} – Confusion matrix heatmap
#'     \item \code{"prior"} – Bar chart comparing CV accuracy of uniform vs
#'           proportional prior (requires \code{prior = "auto"})
#'   }
#'   Default: \code{"scores"}
#' @param color_by Character string specifying metadata column to color points
#'        by in the scores plot (default: uses outcome_var from LDA)
#' @param ld_axis Integer specifying which LD axis to show in the
#'        \code{"coefficients"} plot (default: 1)
#' @param top_n Integer specifying how many top features (by absolute loading)
#'        to show in the coefficients plot (default: 20)
#' @param point_size Numeric specifying size of points in scores plot (default: 3)
#' @param add_labels Logical indicating whether to add sample labels in scores
#'        plot (default: FALSE)
#' @param label_size Numeric specifying size of labels (default: 3)
#' @param bar_color Character string specifying bar fill color for the
#'        coefficients plot (default: \code{"#4575b4"})
#'
#' @return ggplot2 object
#'
#' @examples
#' lda_results <- perform_LDA(df_prepped, outcome_var = "Group")
#'
#' plot_LDA(lda_results, plot_type = "scores")
#' plot_LDA(lda_results, plot_type = "scores", color_by = "Group2",
#'          add_labels = TRUE)
#' plot_LDA(lda_results, plot_type = "coefficients", ld_axis = 1, top_n = 15)
#' plot_LDA(lda_results, plot_type = "roc")
#' plot_LDA(lda_results, plot_type = "confusion")
#' plot_LDA(lda_results, plot_type = "prior")
#' @export
plot_LDA <- function(lda_results,
                     plot_type  = "scores",
                     color_by   = NULL,
                     ld_axis    = 1,
                     top_n      = 20,
                     point_size = 3,
                     add_labels = FALSE,
                     label_size = 3,
                     bar_color  = "#4575b4") {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.")

  if (!inherits(lda_results, "perform_LDA"))
    stop("Input must be output from perform_LDA function")

  plot_type <- match.arg(plot_type,
                         c("scores", "coefficients", "roc", "confusion", "prior"))

  switch(plot_type,
         scores       = .plot_lda_scores(lda_results, color_by, point_size,
                                         add_labels, label_size),
         coefficients = .plot_lda_coefficients(lda_results, ld_axis, top_n, bar_color),
         roc          = .plot_lda_roc(lda_results),
         confusion    = .plot_lda_confusion(lda_results),
         prior        = .plot_lda_prior(lda_results)
  )
}


# ---- internal: LD scores plot ----
.plot_lda_scores <- function(lda_results, color_by, point_size,
                             add_labels, label_size) {

  color_by  <- color_by %||% lda_results$outcome_var
  scores_df <- lda_results$ld_scores
  var_exp   <- round(lda_results$variance_explained * 100, 1)
  n_ld      <- length(var_exp)

  if (!color_by %in% colnames(lda_results$metadata))
    stop(paste0("color_by variable '", color_by, "' not found in metadata. ",
                "Available columns: ",
                paste(colnames(lda_results$metadata), collapse = ", ")))

  scores_df$ColorVar <- lda_results$metadata[[color_by]]

  cv_label <- switch(lda_results$parameters$cv_method,
                     cv         = paste0(lda_results$parameters$cross, "-fold CV"),
                     LOOCV      = "LOOCV",
                     repeatedcv = paste0(lda_results$parameters$cross, "-fold CV \u00d7 ",
                                         lda_results$parameters$repeats, " repeats")
  )

  base_sub <- paste0(cv_label,
                     "  |  PCA: ", lda_results$parameters$pca_preprocess,
                     "  |  CV accuracy: ",
                     round(lda_results$accuracy, 2), "%")

  # 1 LD axis → density plot
  if (n_ld == 1) {
    p <- ggplot2::ggplot(scores_df,
                         ggplot2::aes(x = LD1, fill = ColorVar,
                                      color = ColorVar)) +
      ggplot2::geom_density(alpha = 0.4) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x        = paste0("LD1 (", var_exp[1], "% between-group variance)"),
        y        = "Density",
        fill     = color_by, color = color_by,
        title    = "LDA — LD1 Score Distribution",
        subtitle = base_sub
      )
  } else {
    p <- ggplot2::ggplot(scores_df,
                         ggplot2::aes(x = LD1, y = LD2, color = ColorVar))
    if (add_labels) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), size = label_size)
    } else {
      p <- p + ggplot2::geom_point(size = point_size)
    }
    p <- p +
      ggplot2::stat_ellipse(ggplot2::aes(group = ColorVar),
                            linetype = "dashed", linewidth = 0.6) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x        = paste0("LD1 (", var_exp["LD1"], "% between-group variance)"),
        y        = paste0("LD2 (", var_exp["LD2"], "% between-group variance)"),
        color    = color_by,
        title    = "LDA — Sample Scores in Discriminant Space",
        subtitle = base_sub
      )
  }

  p + ggplot2::theme(
    plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9)
  )
}


# ---- internal: LD coefficient / loading plot ----
.plot_lda_coefficients <- function(lda_results, ld_axis, top_n, bar_color) {

  coef_df <- lda_results$ld_coefficients
  ld_col  <- paste0("LD", ld_axis)
  var_exp <- round(lda_results$variance_explained * 100, 1)

  if (!ld_col %in% colnames(coef_df))
    stop(paste0("LD axis '", ld_col, "' not found. Available: ",
                paste(names(lda_results$variance_explained), collapse = ", ")))

  top_n    <- min(top_n, nrow(coef_df))
  coef_ord <- coef_df[order(abs(coef_df[[ld_col]]), decreasing = TRUE), ]
  coef_top <- coef_ord[seq_len(top_n), ]

  coef_top$Feature   <- factor(coef_top$Feature,
                               levels = coef_top$Feature[
                                 order(coef_top[[ld_col]])])
  coef_top$Direction <- ifelse(coef_top[[ld_col]] >= 0, "Positive", "Negative")

  ggplot2::ggplot(coef_top,
                  ggplot2::aes(x = Feature, y = .data[[ld_col]],
                               fill = Direction)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("Positive" = "#d73027",
                                          "Negative" = "#4575b4")) +
    ggplot2::geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = NULL,
      y        = paste0("Standardised LD", ld_axis, " Coefficient"),
      fill     = "Direction",
      title    = paste0("LDA — Top ", top_n,
                        " Feature Loadings on LD", ld_axis),
      subtitle = paste0("LD", ld_axis, " explains ",
                        var_exp[paste0("LD", ld_axis)],
                        "% of between-group variance")
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text.y   = ggplot2::element_text(size = 9)
    )
}


# ---- internal: ROC curve (binary only) ----
.plot_lda_roc <- function(lda_results) {

  outcome <- as.factor(lda_results$metadata[[lda_results$outcome_var]])
  lvls    <- levels(outcome)

  if (length(lvls) != 2)
    stop("ROC curve is only supported for binary outcomes. ",
         "Detected classes: ", paste(lvls, collapse = ", "))

  probs     <- lda_results$probabilities
  pos_class <- lvls[2]
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

  cv_label <- switch(lda_results$parameters$cv_method,
                     cv         = paste0(lda_results$parameters$cross, "-fold CV"),
                     LOOCV      = "LOOCV",
                     repeatedcv = paste0(lda_results$parameters$cross, "-fold CV \u00d7 ",
                                         lda_results$parameters$repeats, " repeats")
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
      title    = paste0("ROC Curve — ", lda_results$outcome_var,
                        "  (", paste(lvls, collapse = " vs "), ")"),
      subtitle = paste0(cv_label,
                        "  |  CV accuracy: ",
                        round(lda_results$accuracy, 2), "%")
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9)
    )
}


# ---- internal: confusion matrix heatmap ----
.plot_lda_confusion <- function(lda_results) {

  conf      <- lda_results$confusion_matrix
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
      title    = "LDA Confusion Matrix",
      subtitle = paste0("CV accuracy: ", round(lda_results$accuracy, 2), "%"),
      fill     = "Count"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text     = ggplot2::element_text(size = 11),
      axis.title    = ggplot2::element_text(size = 12)
    )
}


# ---- internal: prior comparison bar chart ----
.plot_lda_prior <- function(lda_results) {

  if (is.null(lda_results$tuning_results))
    stop("No prior tuning results found. ",
         "Re-run perform_LDA with prior = 'auto'.")

  td      <- lda_results$tuning_results
  best_pr <- td$prior[which.max(td$cv_accuracy)]

  ggplot2::ggplot(td, ggplot2::aes(x = prior, y = cv_accuracy,
                                   fill = prior == best_pr)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(cv_accuracy, 2), "%")),
                       vjust = -0.4, size = 4.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("TRUE"  = "#d73027",
                                          "FALSE" = "#4575b4"),
                               guide  = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = "Prior Type",
      y        = "CV Accuracy (%)",
      title    = "LDA — Prior Selection via Cross-Validation",
      subtitle = paste0("Selected: '", best_pr, "'")
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9)
    )
}
