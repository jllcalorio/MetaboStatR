#' Perform Decision Tree Classification on Preprocessed Metabolomics Data
#'
#' Uses conditional inference trees (ctree) from the \code{party} package.
#' Tree complexity is controlled via the \code{mincriterion} parameter, which
#' sets the threshold for the test statistic that must be exceeded before a
#' split is made (analogous to pruning in rpart). A higher value produces a
#' smaller, more conservative tree.
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
#' @param tune_mincriterion Logical indicating whether to automatically tune
#'        \code{mincriterion} via cross-validation (default: TRUE)
#' @param mincriterion Numeric value in (0, 1) or numeric vector of values to
#'        search when \code{tune_mincriterion = TRUE}. When
#'        \code{tune_mincriterion = FALSE}, the first element is used directly.
#'        (default: \code{c(0.75, 0.85, 0.90, 0.95, 0.99)})
#' @param cross Integer specifying number of cross-validation folds used during
#'        tuning (default: 10)
#' @param minsplit Integer specifying the minimum number of samples that must
#'        exist in a node for a split to be attempted (default: 20)
#' @param minbucket Integer specifying the minimum number of samples in any
#'        terminal node (default: 7)
#' @param seed Integer for reproducibility (default: 42)
#'
#' @return Object of class "perform_DT" (list) containing:
#'   \itemize{
#'     \item \code{model} – Final ctree model object
#'     \item \code{accuracy} – Cross-validated accuracy from tuning (or
#'           training accuracy if \code{tune_mincriterion = FALSE})
#'     \item \code{confusion_matrix} – Confusion matrix with class error rates
#'     \item \code{variable_importance} – Named numeric vector of variable
#'           importance scores (based on split statistic aggregation)
#'     \item \code{tuning_results} – Data frame of mincriterion vs CV error
#'           (NULL if \code{tune_mincriterion = FALSE})
#'     \item \code{optimal_mincriterion} – Optimal or specified mincriterion value
#'     \item \code{metadata} – Metadata used (QC samples removed)
#'     \item \code{data_matrix} – Data matrix used
#'     \item \code{outcome_var} – Name of outcome variable
#'     \item \code{features_used} – Character vector of features in model
#'     \item \code{parameters} – List of function call parameters
#'   }
#'
#' @examples
#' dt_results <- perform_DT(df_prepped, outcome_var = "Group")
#'
#' # Manual mincriterion, no tuning
#' dt_results <- perform_DT(df_prepped, outcome_var = "Group",
#'                           tune_mincriterion = FALSE, mincriterion = 0.95)
#'
#' # Specific features with custom tuning grid
#' dt_results <- perform_DT(df_prepped,
#'                           features      = c("Met001", "Met010", "Met050"),
#'                           mincriterion  = c(0.80, 0.90, 0.95, 0.99))
#' @export
perform_DT <- function(prepped_data,
                       outcome_var        = "Group",
                       features           = NULL,
                       data_type          = "NONPLS",
                       use_merged         = TRUE,
                       tune_mincriterion  = TRUE,
                       mincriterion       = c(0.75, 0.85, 0.90, 0.95, 0.99),
                       cross              = 10,
                       minsplit           = 20,
                       minbucket          = 7,
                       seed               = 42) {

  # ---- package checks ----
  if (!requireNamespace("party", quietly = TRUE))
    stop("Package 'party' is required but not installed.")

  # ---- input validation ----
  if (!inherits(prepped_data, "perform_PreprocessingPeakData"))
    stop("Input must be output from perform_PreprocessingPeakData function")

  if (!data_type %in% c("NONPLS", "PLS"))
    stop("data_type must be either 'NONPLS' or 'PLS'")

  if (any(mincriterion <= 0) || any(mincriterion >= 1))
    stop("All mincriterion values must be in the open interval (0, 1)")

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

  # Sanitise column names — ctree uses formula interface which calls
  # make.names() internally; we do it explicitly so predict-time names match.
  colnames(data_matrix) <- make.names(colnames(data_matrix), unique = TRUE)
  features_to_use       <- make.names(features_to_use,       unique = TRUE)
  model_data            <- data.frame(outcome = outcome, data_matrix)

  set.seed(seed)

  # ---- ctree control builder ----
  .make_controls <- function(mc) {
    party::ctree_control(
      mincriterion = mc,
      minsplit     = minsplit,
      minbucket    = minbucket
    )
  }

  # ---- tune mincriterion via k-fold CV ----
  tuning_results    <- NULL
  optimal_mc        <- NULL
  cv_accuracy_final <- NULL

  if (tune_mincriterion) {
    message(paste0("Tuning mincriterion via ", cross, "-fold cross-validation..."))

    n         <- nrow(model_data)
    fold_ids  <- sample(rep(seq_len(cross), length.out = n))

    cv_errors <- vapply(mincriterion, function(mc) {
      fold_errs <- vapply(seq_len(cross), function(k) {
        train <- model_data[fold_ids != k, , drop = FALSE]
        test  <- model_data[fold_ids == k, , drop = FALSE]
        fit   <- party::ctree(outcome ~ ., data = train,
                              controls = .make_controls(mc))
        preds <- stats::predict(fit, newdata = test)
        mean(preds != test$outcome)
      }, numeric(1))
      mean(fold_errs)
    }, numeric(1))

    tuning_results <- data.frame(
      mincriterion = mincriterion,
      cv_error     = cv_errors,
      cv_accuracy  = 1 - cv_errors
    )

    optimal_mc        <- mincriterion[which.min(cv_errors)]
    cv_accuracy_final <- (1 - min(cv_errors)) * 100

    message(paste0("Optimal mincriterion: ", optimal_mc,
                   "  (CV accuracy: ", round(cv_accuracy_final, 2), "%)"))

  } else {
    optimal_mc <- mincriterion[1]
    message(paste0("Using specified mincriterion: ", optimal_mc))
  }

  # ---- build final model ----
  message("Building final ctree model...")
  final_model <- party::ctree(outcome ~ .,
                              data     = model_data,
                              controls = .make_controls(optimal_mc))

  # ---- confusion matrix ----
  predicted  <- stats::predict(final_model, model_data)
  conf_table <- table(Predicted = predicted, Actual = outcome)
  class_err  <- 1 - diag(conf_table) / colSums(conf_table)
  conf_matrix <- cbind(conf_table, class.error = round(class_err, 4))

  # training accuracy if no CV was done
  if (is.null(cv_accuracy_final))
    cv_accuracy_final <- mean(predicted == outcome) * 100

  # ---- variable importance via split statistic aggregation ----
  # Primary: walk inner nodes accumulating test statistics.
  # Fallback 1: count how many times each variable is used as a split variable.
  # Fallback 2: return a zero-named vector (tree has no splits at all).
  .extract_varimp <- function(tree_obj) {
    varimp_vec <- setNames(numeric(length(features_to_use)), features_to_use)

    node_ids <- tryCatch(unique(party::where(tree_obj)), error = function(e) integer(0))
    if (length(node_ids) == 0) return(varimp_vec)

    nodes <- tryCatch(party::nodes(tree_obj, node_ids), error = function(e) list())

    has_splits <- FALSE
    for (nd in nodes) {
      if (!party::is.leaf(nd)) {
        has_splits <- TRUE
        tryCatch({
          split_var <- nd@psplit@variableID
          criterion <- max(nd@criterion@statistic, na.rm = TRUE)
          var_name  <- features_to_use[split_var]
          if (!is.na(var_name) && var_name %in% names(varimp_vec))
            varimp_vec[var_name] <- varimp_vec[var_name] + criterion
        }, error = function(e) NULL)
      }
    }

    # Fallback: split-count importance when statistic slots are inaccessible
    if (!has_splits || all(varimp_vec == 0)) {
      message("Statistic-based importance unavailable; using split-count fallback.")
      for (nd in nodes) {
        if (!party::is.leaf(nd)) {
          tryCatch({
            split_var <- nd@psplit@variableID
            var_name  <- features_to_use[split_var]
            if (!is.na(var_name) && var_name %in% names(varimp_vec))
              varimp_vec[var_name] <- varimp_vec[var_name] + 1
          }, error = function(e) NULL)
        }
      }
    }

    sort(varimp_vec[varimp_vec > 0], decreasing = TRUE)
  }

  var_importance <- tryCatch(
    .extract_varimp(final_model),
    error = function(e) {
      warning("Variable importance extraction failed. Returning empty vector.")
      setNames(numeric(0), character(0))
    }
  )

  # ---- output ----
  results <- list(
    model                = final_model,
    accuracy             = cv_accuracy_final,
    confusion_matrix     = conf_matrix,
    variable_importance  = var_importance,
    tuning_results       = tuning_results,
    optimal_mincriterion = optimal_mc,
    metadata             = metadata,
    data_matrix          = data_matrix,
    outcome_var          = outcome_var,
    features_used        = features_to_use,
    parameters           = list(
      data_type  = data_type,
      use_merged = use_merged,
      minsplit   = minsplit,
      minbucket  = minbucket,
      cross      = cross,
      seed       = seed
    )
  )

  class(results) <- c("perform_DT", "list")

  message("Decision Tree analysis complete!")
  message(paste0("Accuracy: ", round(cv_accuracy_final, 2), "%"))

  return(results)
}


# ---- internal null-coalescing helper (skip if already defined) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}


#' Create Plots from Decision Tree Results
#'
#' @param dt_results Output from perform_DT function
#' @param plot_type Character string specifying which plot to create:
#'   \itemize{
#'     \item \code{"tree"} – Rendered tree diagram using \code{party}'s built-in
#'           plot method (printed directly; returns invisibly)
#'     \item \code{"varimp"} – Horizontal bar chart of variable importance scores
#'     \item \code{"roc"} – ROC curve with AUC (binary outcomes only)
#'     \item \code{"confusion"} – Confusion matrix heatmap with class error rates
#'     \item \code{"tuning"} – mincriterion tuning profile (requires
#'           \code{tune_mincriterion = TRUE} in \code{perform_DT})
#'   }
#'   Default: \code{"tree"}
#' @param top_n Integer specifying how many top features to show in the variable
#'        importance plot (default: 20)
#' @param bar_color Character string specifying bar fill color for the variable
#'        importance plot (default: "#4575b4")
#'
#' @return A ggplot2 object for all plot types except \code{"tree"}, which
#'         renders via \code{party}'s own plot method and returns invisibly.
#'
#' @examples
#' dt_results <- perform_DT(df_prepped, outcome_var = "Group")
#'
#' plot_DT(dt_results, plot_type = "tree")
#' plot_DT(dt_results, plot_type = "varimp")
#' plot_DT(dt_results, plot_type = "varimp", top_n = 10)
#' plot_DT(dt_results, plot_type = "roc")
#' plot_DT(dt_results, plot_type = "confusion")
#' plot_DT(dt_results, plot_type = "tuning")
#' @export
plot_DT <- function(dt_results,
                    plot_type  = "tree",
                    top_n      = 20,
                    bar_color  = "#4575b4") {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.")

  if (!inherits(dt_results, "perform_DT"))
    stop("Input must be output from perform_DT function")

  plot_type <- match.arg(plot_type, c("tree", "varimp", "roc", "confusion", "tuning"))

  switch(plot_type,
         tree      = .plot_dt_tree(dt_results),
         varimp    = .plot_dt_varimp(dt_results, top_n, bar_color),
         roc       = .plot_dt_roc(dt_results),
         confusion = .plot_dt_confusion(dt_results),
         tuning    = .plot_dt_tuning(dt_results)
  )
}


# ---- internal: tree diagram ----
.plot_dt_tree <- function(dt_results) {
  message("Rendering ctree diagram via party's plot method...")
  graphics::plot(dt_results$model,
                 main = paste0("Conditional Inference Tree — ",
                               dt_results$outcome_var,
                               "  (mincriterion = ",
                               dt_results$optimal_mincriterion, ")"))
  invisible(dt_results$model)
}


# ---- internal: variable importance bar chart ----
.plot_dt_varimp <- function(dt_results, top_n, bar_color) {

  vi <- dt_results$variable_importance

  if (length(vi) == 0)
    stop("No variable importance scores available — the ctree has no splits. ",
         "Try lowering mincriterion (e.g. mincriterion = 0.75) or reducing ",
         "minsplit / minbucket, or providing more / different features.")

  top_n  <- min(top_n, length(vi))
  vi_df  <- data.frame(
    Feature    = factor(names(vi)[seq_len(top_n)],
                        levels = rev(names(vi)[seq_len(top_n)])),
    Importance = vi[seq_len(top_n)]
  )

  ggplot2::ggplot(vi_df, ggplot2::aes(x = Feature, y = Importance)) +
    ggplot2::geom_bar(stat = "identity", fill = bar_color) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x     = NULL,
      y     = "Aggregated Split Statistic",
      title = paste0("Variable Importance — Top ", top_n, " Features")
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9)
    )
}


# ---- internal: ROC curve (binary only) ----
.plot_dt_roc <- function(dt_results) {

  outcome <- as.factor(dt_results$metadata[[dt_results$outcome_var]])
  lvls    <- levels(outcome)

  if (length(lvls) != 2)
    stop("ROC curve is only supported for binary outcomes. ",
         "Detected classes: ", paste(lvls, collapse = ", "))

  # ctree predict with type = "prob" returns a list of probability matrices
  prob_list <- stats::predict(dt_results$model,
                              newdata = data.frame(dt_results$data_matrix),
                              type    = "prob")
  # Each element is a named vector; extract positive class probability
  pos_class <- lvls[2]
  pos_prob  <- vapply(prob_list, function(x) x[pos_class], numeric(1))
  truth_bin <- as.integer(outcome == pos_class)

  # Compute ROC points
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
      title = paste0("ROC Curve — ", dt_results$outcome_var,
                     "  (", paste(lvls, collapse = " vs "), ")")
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
}


# ---- internal: confusion matrix heatmap ----
.plot_dt_confusion <- function(dt_results) {

  conf      <- dt_results$confusion_matrix
  class_err <- conf[, "class.error"]
  conf_vals <- conf[, colnames(conf) != "class.error", drop = FALSE]

  conf_df   <- as.data.frame(as.table(conf_vals))
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
      title = paste0("Decision Tree Confusion Matrix",
                     "  (mincriterion = ", dt_results$optimal_mincriterion, ")"),
      fill  = "Count"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text  = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12)
    )
}


# ---- internal: mincriterion tuning profile ----
.plot_dt_tuning <- function(dt_results) {

  if (is.null(dt_results$tuning_results))
    stop("No tuning results found. ",
         "Re-run perform_DT with tune_mincriterion = TRUE.")

  td      <- dt_results$tuning_results
  opt_mc  <- dt_results$optimal_mincriterion

  ggplot2::ggplot(td, ggplot2::aes(x = mincriterion, y = cv_error)) +
    ggplot2::geom_line(color = "#4575b4", linewidth = 1) +
    ggplot2::geom_point(color = "#4575b4", size = 3) +
    ggplot2::geom_vline(xintercept = opt_mc,
                        linetype = "dashed", color = "red", linewidth = 0.8) +
    ggplot2::annotate("text",
                      x     = opt_mc,
                      y     = max(td$cv_error),
                      label = paste0("Optimal = ", opt_mc),
                      hjust = -0.1, color = "red", size = 3.5) +
    ggplot2::scale_x_continuous(breaks = td$mincriterion) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x     = "mincriterion",
      y     = "Cross-validation Error",
      title = "Decision Tree Tuning — mincriterion"
    ) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
