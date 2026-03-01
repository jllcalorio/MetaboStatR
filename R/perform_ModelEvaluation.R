# ==============================================================================
#  perform_ModelEvaluation()  +  plot_ModelEvaluation()
#
#  Packages required:
#    CRAN : kernelshap, shapviz, ggplot2, patchwork, ggrepel
#
#  kernelshap  — model-agnostic Kernel SHAP; works with all 6 model classes
#  shapviz     — tidy SHAP visualisation layer on top of kernelshap output
# ==============================================================================


# ---- internal null-coalescing helper (skip if already defined) ----
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}


# ==============================================================================
#  Internal helpers
# ==============================================================================

# Safe column name sanitiser — stores original→sanitized mapping
.sanitise_names <- function(x) make.names(x, unique = TRUE)

# Rebuild a data.frame with sanitised column names ONLY for formula-based models
.safe_df <- function(newdata, sanitise = FALSE) {
  nd <- as.data.frame(newdata)
  if (sanitise) colnames(nd) <- .sanitise_names(colnames(nd))
  nd
}

# Unified predict wrapper — returns a numeric probability matrix (n × classes)
.predict_prob <- function(model_obj, newdata) {

  cls  <- class(model_obj)[1]
  lvls <- levels(as.factor(model_obj$metadata[[model_obj$outcome_var]]))

  # Check if the trained model has sanitized feature names
  .model_uses_sanitized_names <- function(model_obj) {
    tryCatch({
      if (cls == "perform_RandomForest") {
        # Check the terms or xlevels
        trained_vars <- names(model_obj$model$forest$xlevels)
        if (length(trained_vars) > 0) {
          # If first variable starts with X followed by digit, it's sanitized
          return(grepl("^X[0-9]", trained_vars[1]))
        }
      } else if (cls == "perform_SVM") {
        # SVM from e1071 with formula also gets sanitized names
        if (!is.null(model_obj$model$terms)) {
          trained_vars <- attr(model_obj$model$terms, "term.labels")
          if (length(trained_vars) > 0) {
            return(grepl("^X[0-9]", trained_vars[1]))
          }
        }
      } else if (cls %in% c("perform_DT", "perform_AdaBoost")) {
        # Formula-based models always sanitize
        return(TRUE)
      } else if (cls == "perform_KNN") {
        # Check caret model's xNames
        if (!is.null(model_obj$model$xNames)) {
          return(grepl("^X[0-9]", model_obj$model$xNames[1]))
        }
      } else if (cls == "perform_LDA") {
        # Check caret model's xNames
        if (!is.null(model_obj$model$xNames)) {
          return(grepl("^X[0-9]", model_obj$model$xNames[1]))
        }
      }
      return(FALSE)
    }, error = function(e) FALSE)
  }

  needs_sanitized <- .model_uses_sanitized_names(model_obj)

  if (cls == "perform_RandomForest") {
    nd    <- .safe_df(newdata, sanitise = needs_sanitized)
    probs <- stats::predict(model_obj$model, newdata = nd, type = "prob")

  } else if (cls == "perform_SVM") {
    nd    <- .safe_df(newdata, sanitise = needs_sanitized)
    probs <- attr(
      stats::predict(model_obj$model, newdata = nd, probability = TRUE),
      "probabilities"
    )

  } else if (cls == "perform_DT") {
    nd          <- .safe_df(newdata, sanitise = TRUE)
    nd$outcome  <- factor(NA, levels = lvls)
    prob_list   <- stats::predict(model_obj$model, newdata = nd, type = "prob")
    probs       <- do.call(rbind, prob_list)

  } else if (cls == "perform_KNN") {
    nd    <- .safe_df(newdata, sanitise = needs_sanitized)
    probs <- as.matrix(stats::predict(model_obj$model, newdata = nd,
                                      type = "prob"))

  } else if (cls == "perform_AdaBoost") {
    nd         <- .safe_df(newdata, sanitise = TRUE)
    nd$outcome <- factor(NA, levels = lvls)
    res        <- stats::predict(model_obj$model, newdata = nd)
    probs      <- res$prob
    colnames(probs) <- lvls

  } else if (cls == "perform_LDA") {
    nd <- .safe_df(newdata, sanitise = needs_sanitized)
    if (!is.null(model_obj$pca_model))
      nd <- stats::predict(model_obj$pca_model, nd)
    probs <- as.matrix(stats::predict(model_obj$model, newdata = nd,
                                      type = "prob"))

  } else {
    stop("Unsupported model class: ", cls)
  }

  as.matrix(probs)
}

# Compute all 9 metrics for one model (binary or multi-class macro-average)
.compute_metrics <- function(model_obj) {

  outcome   <- as.factor(model_obj$metadata[[model_obj$outcome_var]])
  lvls      <- levels(outcome)
  newdata   <- model_obj$data_matrix
  probs     <- .predict_prob(model_obj, newdata)

  # Ensure probs has proper column names
  if (is.null(colnames(probs))) {
    colnames(probs) <- lvls
  }

  predicted <- lvls[apply(probs, 1, which.max)]
  predicted <- factor(predicted, levels = lvls)

  is_binary <- length(lvls) == 2

  # ---- AUC (hand-rolled, works for binary + macro-average multiclass) ----
  .auc_binary <- function(truth_bin, pos_prob) {
    ord  <- order(pos_prob, decreasing = TRUE)
    tbin <- truth_bin[ord]
    n_p  <- sum(tbin);  n_n <- sum(1 - tbin)
    if (n_p == 0 || n_n == 0) return(NA_real_)
    tp_cum <- cumsum(tbin)
    fp_cum <- cumsum(1 - tbin)
    sum(diff(c(0, fp_cum / n_n)) * (c(0, tp_cum / n_p)[-length(tp_cum) - 1] +
                                      c(0, tp_cum / n_p)[-1]) / 2)
  }

  if (is_binary) {
    pos_class <- lvls[2]
    # Safely extract positive class probabilities
    if (pos_class %in% colnames(probs)) {
      pos_prob <- probs[, pos_class]
    } else {
      # Fallback: use second column
      pos_prob <- probs[, 2]
    }
    auc <- .auc_binary(as.integer(outcome == pos_class), pos_prob)
  } else {
    auc_vals <- vapply(lvls, function(lv) {
      if (lv %in% colnames(probs)) {
        .auc_binary(as.integer(outcome == lv), probs[, lv])
      } else {
        # Find column by index
        idx <- which(lvls == lv)
        .auc_binary(as.integer(outcome == lv), probs[, idx])
      }
    }, numeric(1))
    auc <- mean(auc_vals, na.rm = TRUE)   # macro-average
  }

  # ---- Per-class TP/FP/TN/FN then macro-average ----
  per_class <- lapply(lvls, function(lv) {
    tp <- sum(predicted == lv & outcome == lv)
    fp <- sum(predicted == lv & outcome != lv)
    tn <- sum(predicted != lv & outcome != lv)
    fn <- sum(predicted != lv & outcome == lv)
    list(tp = tp, fp = fp, tn = tn, fn = fn)
  })

  .macro <- function(fn_metric) {
    vals <- vapply(per_class, fn_metric, numeric(1))
    mean(vals, na.rm = TRUE)
  }

  accuracy  <- mean(predicted == outcome)
  precision <- .macro(function(x) if ((x$tp + x$fp) == 0) NA else x$tp / (x$tp + x$fp))
  recall    <- .macro(function(x) if ((x$tp + x$fn) == 0) NA else x$tp / (x$tp + x$fn))
  fpr       <- .macro(function(x) if ((x$fp + x$tn) == 0) NA else x$fp / (x$fp + x$tn))
  f1        <- .macro(function(x) {
    p <- if ((x$tp + x$fp) == 0) NA else x$tp / (x$tp + x$fp)
    r <- if ((x$tp + x$fn) == 0) NA else x$tp / (x$tp + x$fn)
    if (is.na(p) || is.na(r) || (p + r) == 0) NA else 2 * p * r / (p + r)
  })
  npv       <- .macro(function(x) if ((x$tn + x$fn) == 0) NA else x$tn / (x$tn + x$fn))
  tnr       <- .macro(function(x) if ((x$tn + x$fp) == 0) NA else x$tn / (x$tn + x$fp))
  fnr       <- .macro(function(x) if ((x$tp + x$fn) == 0) NA else x$fn / (x$tp + x$fn))

  list(
    AUC       = auc,
    Accuracy  = accuracy,
    Precision = precision,
    Recall    = recall,
    FPR       = fpr,
    F1        = f1,
    NPV       = npv,
    TNR       = tnr,
    FNR       = fnr,
    predicted = predicted,
    probs     = probs,
    outcome   = outcome,
    lvls      = lvls
  )
}


# Friendly model name extractor
.model_name <- function(model_obj) {
  switch(class(model_obj)[1],
         perform_RandomForest = "Random Forest",
         perform_SVM          = paste0("SVM (", model_obj$kernel, ")"),
         perform_DT           = "Decision Tree",
         perform_KNN          = paste0("KNN (K=", model_obj$optimal_k, ")"),
         perform_AdaBoost     = paste0("AdaBoost (", model_obj$optimal_mfinal, " iter)"),
         perform_LDA          = "LDA",
         class(model_obj)[1]
  )
}


# ==============================================================================
#  perform_ModelEvaluation()
# ==============================================================================

#' Compute Performance Metrics and SHAP Values Across Multiple Models
#'
#' Accepts two or more model result objects from \code{perform_RandomForest},
#' \code{perform_SVM}, \code{perform_DT}, \code{perform_KNN},
#' \code{perform_AdaBoost}, or \code{perform_LDA} and returns a unified
#' evaluation object containing all nine performance metrics and SHAP values
#' for each model.
#'
#' SHAP values are computed via \code{kernelshap} (model-agnostic Kernel SHAP)
#' and stored as \code{shapviz} objects for direct plotting.
#'
#' @param ... Two or more model result objects. Each must be an output from one
#'        of the six supported \code{perform_*()} functions. Objects must all
#'        share the same \code{outcome_var} and compatible feature sets.
#' @param model_names Optional character vector of custom display names for
#'        each model. Must be the same length as the number of models supplied.
#'        If NULL, names are derived automatically (default: NULL)
#' @param compute_shap Logical indicating whether to compute SHAP values via
#'        \code{kernelshap} + \code{shapviz} (default: TRUE). Set to FALSE to
#'        skip SHAP (faster)
#' @param shap_nsim Integer specifying the number of Monte Carlo samples used
#'        by \code{kernelshap} per observation (default: 100)
#' @param shap_bg_n Integer specifying the number of background samples to use
#'        for the SHAP reference distribution. Larger values are more accurate
#'        but slower (default: 50)
#' @param shap_type Character string controlling which SHAP plot style is
#'        pre-computed and stored. One of \code{"bar"}, \code{"beeswarm"}, or
#'        \code{"both"} (default: \code{"both"}). This does not affect what is
#'        stored — all \code{shapviz} objects support both styles — but is
#'        passed through to the output for reference by \code{plot_ModelEvaluation}
#' @param seed Integer for reproducibility (default: 42)
#'
#' @return Object of class "perform_ModelEvaluation" (list) containing:
#'   \itemize{
#'     \item \code{metrics_table} – Data frame of all 9 metrics × models
#'     \item \code{models} – Named list of the input model objects
#'     \item \code{model_names} – Character vector of display names
#'     \item \code{shap} – Named list of \code{shapviz} objects, one per model
#'           (NULL if \code{compute_shap = FALSE})
#'     \item \code{roc_data} – Named list of ROC data frames per model
#'     \item \code{outcome_var} – Shared outcome variable name
#'     \item \code{parameters} – List of function call parameters
#'   }
#'
#' @examples
#' rf_res  <- perform_RandomForest(df_prepped, outcome_var = "Group")
#' svm_res <- perform_SVM(df_prepped,          outcome_var = "Group")
#' knn_res <- perform_KNN(df_prepped,          outcome_var = "Group")
#'
#' eval_res <- perform_ModelEvaluation(rf_res, svm_res, knn_res)
#'
#' # Custom names, no SHAP
#' eval_res <- perform_ModelEvaluation(rf_res, svm_res,
#'                                      model_names   = c("RF", "SVM"),
#'                                      compute_shap  = FALSE)
#' @export
perform_ModelEvaluation <- function(...,
                                    model_names   = NULL,
                                    compute_shap  = TRUE,
                                    shap_nsim     = 100,
                                    shap_bg_n     = 50,
                                    shap_type     = "both",
                                    seed          = 42) {

  shap_type <- match.arg(shap_type, c("both", "bar", "beeswarm"))

  # ---- package checks ----
  for (pkg in c("ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(paste0("Package '", pkg, "' is required but not installed."))
  }
  if (compute_shap) {
    for (pkg in c("kernelshap", "shapviz")) {
      if (!requireNamespace(pkg, quietly = TRUE))
        stop(paste0("Package '", pkg, "' is required for SHAP. ",
                    "Install via install.packages('", pkg, "') or set ",
                    "compute_shap = FALSE."))
    }
  }

  # ---- collect models ----
  model_list <- list(...)
  n_models   <- length(model_list)

  if (n_models < 2)
    stop("At least 2 model result objects must be supplied.")

  supported <- c("perform_RandomForest", "perform_SVM", "perform_DT",
                 "perform_KNN", "perform_AdaBoost", "perform_LDA")

  for (i in seq_len(n_models)) {
    if (!any(class(model_list[[i]]) %in% supported))
      stop(paste0("Model ", i, " is not a supported perform_*() result object."))
  }

  # ---- validate shared outcome_var ----
  outcome_vars <- vapply(model_list, `[[`, character(1), "outcome_var")
  if (length(unique(outcome_vars)) > 1)
    stop("All models must share the same outcome_var. Found: ",
         paste(unique(outcome_vars), collapse = ", "))
  outcome_var <- outcome_vars[1]

  # ---- auto model names ----
  if (is.null(model_names)) {
    model_names <- vapply(model_list, .model_name, character(1))
    # Deduplicate if same model type run twice
    if (anyDuplicated(model_names)) {
      dups <- which(duplicated(model_names) | duplicated(model_names, fromLast = TRUE))
      model_names[dups] <- paste0(model_names[dups], " (", seq_along(dups), ")")
    }
  } else {
    if (length(model_names) != n_models)
      stop("model_names must have the same length as the number of models (",
           n_models, ").")
  }

  names(model_list) <- model_names

  set.seed(seed)

  # ---- compute metrics for each model ----
  message("Computing performance metrics...")
  metrics_list <- lapply(seq_len(n_models), function(i) {
    message(paste0("  [", i, "/", n_models, "] ", model_names[i]))
    .compute_metrics(model_list[[i]])
  })
  names(metrics_list) <- model_names

  # ---- assemble metrics table ----
  metric_cols <- c("AUC", "Accuracy", "Precision", "Recall",
                   "FPR", "F1", "NPV", "TNR", "FNR")

  metrics_table <- do.call(rbind, lapply(seq_len(n_models), function(i) {
    m  <- metrics_list[[i]]
    df <- as.data.frame(lapply(metric_cols, function(mc) round(m[[mc]], 4)))
    colnames(df) <- metric_cols
    cbind(Model = model_names[i], df)
  }))
  rownames(metrics_table) <- NULL

  # ---- ROC data per model ----
  roc_data <- lapply(seq_len(n_models), function(i) {
    m    <- metrics_list[[i]]
    lvls <- m$lvls

    if (length(lvls) != 2) {
      message(paste0("  Note: ROC for '", model_names[i],
                     "' uses macro-average (multiclass)"))
    }

    # For binary: single curve; for multiclass: macro-average
    if (length(lvls) == 2) {
      pos_class <- lvls[2]
      pos_prob  <- m$probs[, pos_class]
      truth_bin <- as.integer(m$outcome == pos_class)
      thrs <- sort(unique(pos_prob), decreasing = TRUE)
      rd   <- do.call(rbind, lapply(thrs, function(thr) {
        p  <- as.integer(pos_prob >= thr)
        tp <- sum(p == 1 & truth_bin == 1)
        fp <- sum(p == 1 & truth_bin == 0)
        fn <- sum(p == 0 & truth_bin == 1)
        tn <- sum(p == 0 & truth_bin == 0)
        data.frame(TPR = tp / (tp + fn), FPR = fp / (fp + tn))
      }))
      rd <- rbind(data.frame(TPR = 0, FPR = 0), rd)
    } else {
      # Macro-average ROC
      roc_per_class <- lapply(lvls, function(lv) {
        pos_prob  <- m$probs[, lv]
        truth_bin <- as.integer(m$outcome == lv)
        thrs <- sort(unique(pos_prob), decreasing = TRUE)
        do.call(rbind, lapply(thrs, function(thr) {
          p  <- as.integer(pos_prob >= thr)
          tp <- sum(p == 1 & truth_bin == 1)
          fp <- sum(p == 1 & truth_bin == 0)
          fn <- sum(p == 0 & truth_bin == 1)
          tn <- sum(p == 0 & truth_bin == 0)
          data.frame(TPR = tp / (tp + fn), FPR = fp / (fp + tn))
        }))
      })
      fpr_grid <- seq(0, 1, length.out = 200)
      tpr_mat  <- sapply(roc_per_class, function(rc) {
        stats::approx(rc$FPR, rc$TPR, xout = fpr_grid,
                      method = "linear", rule = 2)$y
      })
      rd <- data.frame(FPR = fpr_grid, TPR = rowMeans(tpr_mat))
      rd <- rbind(data.frame(FPR = 0, TPR = 0), rd)
    }

    auc_val <- round(
      sum(diff(rd$FPR) * (rd$TPR[-1] + rd$TPR[-nrow(rd)]) / 2), 4
    )
    rd$Model <- model_names[i]
    rd$AUC   <- auc_val
    rd
  })
  names(roc_data) <- model_names

  # ---- SHAP values via kernelshap + shapviz ----
  shap_list <- NULL

  if (compute_shap) {
    message("Computing SHAP values via kernelshap...")
    shap_list <- lapply(seq_len(n_models), function(i) {
      message(paste0("  [", i, "/", n_models, "] ", model_names[i]))
      mod <- model_list[[i]]
      X   <- as.data.frame(mod$data_matrix)

      # Background dataset (subsample for speed)
      bg_idx <- sample(nrow(X), min(shap_bg_n, nrow(X)), replace = FALSE)
      bg     <- X[bg_idx, , drop = FALSE]

      # Unified predict function returning numeric matrix (n × classes)
      pfun <- function(object, newdata) {
        p <- .predict_prob(object, as.data.frame(newdata))
        p   # kernelshap handles multi-output automatically
      }

      ks <- tryCatch(
        kernelshap::kernelshap(
          object   = mod,
          X        = X,
          bg_X     = bg,
          pred_fun = pfun,
          verbose  = FALSE#,
          # nsim     = shap_nsim
        ),
        error = function(e) {
          warning(paste0("SHAP failed for '", model_names[i], "': ", e$message))
          NULL
        }
      )

      if (is.null(ks)) return(NULL)
      shapviz::shapviz(ks)
    })
    names(shap_list) <- model_names
    message("SHAP computation complete.")
  }

  # ---- output ----
  results <- list(
    metrics_table = metrics_table,
    models        = model_list,
    model_names   = model_names,
    shap          = shap_list,
    roc_data      = roc_data,
    outcome_var   = outcome_var,
    parameters    = list(
      compute_shap = compute_shap,
      shap_nsim    = shap_nsim,
      shap_bg_n    = shap_bg_n,
      shap_type    = shap_type,
      seed         = seed
    )
  )

  class(results) <- c("perform_ModelEvaluation", "list")

  message("\nModel evaluation complete!")
  message(paste0("Models evaluated: ", paste(model_names, collapse = ", ")))
  print(metrics_table[, c("Model", "AUC", "Accuracy", "F1")])

  return(results)
}


# ==============================================================================
#  plot_ModelEvaluation()
# ==============================================================================

#' Create Plots from Model Evaluation Results
#'
#' @param eval_results Output from \code{perform_ModelEvaluation}
#' @param plot_type Character string specifying which plot to create:
#'   \itemize{
#'     \item \code{"metrics_table"} – Heatmap table of all 9 metrics × models
#'     \item \code{"metrics_bar"} – Bar chart comparing one metric across models
#'     \item \code{"radar"} – Radar / spider chart of all metrics per model
#'     \item \code{"roc"} – Overlaid ROC curves for all models with AUC labels
#'     \item \code{"shap"} – SHAP importance plot; style controlled by
#'           \code{shap_type} parameter
#'   }
#'   Default: \code{"metrics_table"}
#' @param metric Character string specifying which metric to use for
#'        \code{"metrics_bar"} (default: \code{"AUC"}). One of:
#'        \code{"AUC"}, \code{"Accuracy"}, \code{"Precision"},
#'        \code{"Recall"}, \code{"FPR"}, \code{"F1"}, \code{"NPV"},
#'        \code{"TNR"}, \code{"FNR"}
#' @param shap_type Character string controlling the SHAP visualization style
#'        when \code{plot_type = "shap"}:
#'   \itemize{
#'     \item \code{"bar"} – Mean |SHAP| horizontal bar chart
#'     \item \code{"beeswarm"} – SHAP value distribution beeswarm summary plot
#'     \item \code{"both"} – Returns a named list with both \code{"bar"} and
#'           \code{"beeswarm"} plots for each selected model
#'   }
#'   Default: inherits from \code{eval_results$parameters$shap_type},
#'   falls back to \code{"both"}
#' @param model_select Character string or NULL specifying which model to show
#'        in SHAP plots. If NULL, produces plots for all models and returns a
#'        named list (default: NULL)
#' @param top_n Integer specifying the number of top features to show in SHAP
#'        plots (default: 15)
#' @param shap_class For multi-class outcomes, the class name to display SHAP
#'        values for. If NULL, uses the first class level (default: NULL)
#' @param bar_palette Character vector of colours for bars / radar fills.
#'        If NULL, uses a built-in palette (default: NULL)
#'
#' @return A ggplot2 object, or a named list of ggplot2 objects when
#'         \code{model_select = NULL} and \code{plot_type} is a SHAP plot
#'
#' @examples
#' eval_res <- perform_ModelEvaluation(rf_res, svm_res, knn_res)
#'
#' plot_ModelEvaluation(eval_res, plot_type = "metrics_table")
#' plot_ModelEvaluation(eval_res, plot_type = "metrics_bar", metric = "F1")
#' plot_ModelEvaluation(eval_res, plot_type = "radar")
#' plot_ModelEvaluation(eval_res, plot_type = "roc")
#' plot_ModelEvaluation(eval_res, plot_type = "shap")
#' plot_ModelEvaluation(eval_res, plot_type = "shap", shap_type = "bar")
#' plot_ModelEvaluation(eval_res, plot_type = "shap", shap_type = "beeswarm")
#' plot_ModelEvaluation(eval_res, plot_type = "shap", shap_type = "both",
#'                      model_select = "Random Forest", top_n = 20)
#' # Multi-class: specify which class to show SHAP for
#' plot_ModelEvaluation(eval_res, plot_type = "shap",
#'                      shap_class = "End")
#' @export
plot_ModelEvaluation <- function(eval_results,
                                 plot_type    = "metrics_table",
                                 metric       = "AUC",
                                 shap_type    = NULL,
                                 model_select = NULL,
                                 top_n        = 15,
                                 shap_class   = NULL,
                                 bar_palette  = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required but not installed.")

  if (!inherits(eval_results, "perform_ModelEvaluation"))
    stop("Input must be output from perform_ModelEvaluation function")

  plot_type <- match.arg(plot_type,
                         c("metrics_table", "metrics_bar", "radar",
                           "roc", "shap"))

  valid_metrics <- c("AUC", "Accuracy", "Precision", "Recall",
                     "FPR", "F1", "NPV", "TNR", "FNR")
  metric <- match.arg(metric, valid_metrics)

  # Inherit shap_type from eval object if not explicitly set
  resolved_shap_type <- shap_type %||%
    eval_results$parameters$shap_type %||% "both"
  resolved_shap_type <- match.arg(resolved_shap_type,
                                  c("both", "bar", "beeswarm"))

  # Default palette
  n_mod <- length(eval_results$model_names)
  default_pal <- c("#4575b4", "#d73027", "#1a9850", "#f46d43",
                   "#74add1", "#a50026", "#313695")
  pal <- bar_palette %||% default_pal[seq_len(n_mod)]

  switch(plot_type,
         metrics_table = .plot_eval_table(eval_results),
         metrics_bar   = .plot_eval_bar(eval_results, metric, pal),
         radar         = .plot_eval_radar(eval_results, pal),
         roc           = .plot_eval_roc(eval_results, pal),
         shap          = .plot_eval_shap(eval_results, resolved_shap_type,
                                         model_select, top_n, shap_class)
  )
}


# ---- internal: metrics heatmap table ----
.plot_eval_table <- function(eval_results) {

  mt      <- eval_results$metrics_table
  met_cols <- c("AUC", "Accuracy", "Precision", "Recall",
                "FPR", "F1", "NPV", "TNR", "FNR")

  # Lower-is-better metrics (shown in reverse colour scale)
  lower_better <- c("FPR", "FNR")

  mt_long <- do.call(rbind, lapply(met_cols, function(mc) {
    data.frame(
      Model  = mt$Model,
      Metric = mc,
      Value  = mt[[mc]],
      LB     = mc %in% lower_better
    )
  }))

  # Normalise within each metric for fill, flip direction for lower-is-better
  mt_long$Fill <- ave(mt_long$Value, mt_long$Metric, FUN = function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
  })
  mt_long$Fill <- ifelse(mt_long$LB, 1 - mt_long$Fill, mt_long$Fill)
  mt_long$Label <- sprintf("%.3f", mt_long$Value)
  mt_long$Metric <- factor(mt_long$Metric, levels = met_cols)
  mt_long$Model  <- factor(mt_long$Model,
                           levels = rev(eval_results$model_names))

  ggplot2::ggplot(mt_long,
                  ggplot2::aes(x = Metric, y = Model, fill = Fill)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 3.5,
                       fontface = "bold", color = "white") +
    ggplot2::scale_fill_gradient(low = "#4575b4", high = "#d73027",
                                 guide = "none") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x     = NULL, y = NULL,
      title = "Model Performance Metrics Comparison"
    ) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10)
    )
}


# ---- internal: single-metric bar chart ----
.plot_eval_bar <- function(eval_results, metric, pal) {

  mt <- eval_results$metrics_table
  mt$Model <- factor(mt$Model, levels = eval_results$model_names)

  lower_better <- metric %in% c("FPR", "FNR")
  best_model   <- mt$Model[if (lower_better) which.min(mt[[metric]])
                           else             which.max(mt[[metric]])]

  ggplot2::ggplot(mt, ggplot2::aes(x = Model, y = .data[[metric]],
                                   fill = Model == best_model)) +
    ggplot2::geom_bar(stat = "identity", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data[[metric]])),
                       vjust = -0.4, size = 4, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                               guide  = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, max(mt[[metric]], na.rm = TRUE) * 1.12),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = NULL, y = metric,
      title    = paste0("Model Comparison — ", metric),
      subtitle = paste0("Best: ", best_model,
                        if (lower_better) " (lower is better)"
                        else              " (higher is better)")
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      axis.text.x   = ggplot2::element_text(angle = 25, hjust = 1)
    )
}


# ---- internal: radar / spider chart ----
.plot_eval_radar <- function(eval_results, pal) {

  mt      <- eval_results$metrics_table
  # For radar: flip FPR and FNR so higher = better for all axes
  radar_cols <- c("AUC", "Accuracy", "Precision", "Recall",
                  "F1", "NPV", "TNR",
                  "1-FPR", "1-FNR")

  mt[["1-FPR"]] <- 1 - mt[["FPR"]]
  mt[["1-FNR"]] <- 1 - mt[["FNR"]]

  n_axes  <- length(radar_cols)
  angles  <- seq(0, 2 * pi, length.out = n_axes + 1)[-(n_axes + 1)]

  radar_df <- do.call(rbind, lapply(seq_len(nrow(mt)), function(i) {
    vals <- as.numeric(mt[i, radar_cols])
    data.frame(
      Model  = mt$Model[i],
      Metric = radar_cols,
      Value  = vals,
      Angle  = angles,
      x      = vals * cos(angles),
      y      = vals * sin(angles)
    )
  }))

  # Close each polygon
  first_pts <- radar_df[radar_df$Metric == radar_cols[1], ]
  first_pts$Metric <- "__close__"
  radar_df <- rbind(radar_df, first_pts)
  radar_df$Model <- factor(radar_df$Model, levels = eval_results$model_names)

  # Axis labels positions
  ax_labels <- data.frame(
    Metric = radar_cols,
    x      = 1.12 * cos(angles),
    y      = 1.12 * sin(angles)
  )

  # Grid circles
  grid_df <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1), function(r) {
    th <- seq(0, 2 * pi, length.out = 200)
    data.frame(x = r * cos(th), y = r * sin(th), r = r)
  }))

  ggplot2::ggplot() +
    ggplot2::geom_path(data = grid_df,
                       ggplot2::aes(x = x, y = y, group = r),
                       color = "grey80", linewidth = 0.3) +
    ggplot2::geom_polygon(data = radar_df[radar_df$Metric != "__close__", ],
                          ggplot2::aes(x = x, y = y, group = Model,
                                       color = Model, fill = Model),
                          alpha = 0.12, linewidth = 1) +
    ggplot2::geom_text(data = ax_labels,
                       ggplot2::aes(x = x, y = y, label = Metric),
                       size = 3, fontface = "bold") +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_fill_manual(values  = pal) +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::labs(
      title    = "Model Comparison — Radar Chart",
      subtitle = "All metrics scaled to [0, 1]  |  FPR & FNR inverted (higher = better)",
      color    = "Model", fill = "Model"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      legend.position = "bottom"
    )
}


# ---- internal: overlaid ROC curves ----
.plot_eval_roc <- function(eval_results, pal) {

  roc_all <- do.call(rbind, eval_results$roc_data)

  # AUC label positions (near top-right of each curve)
  auc_labels <- do.call(rbind, lapply(eval_results$model_names, function(nm) {
    rd <- eval_results$roc_data[[nm]]
    data.frame(
      Model = nm,
      FPR   = 0.62,
      TPR   = 0.10 + 0.08 * (which(eval_results$model_names == nm) - 1),
      Label = paste0(nm, "  AUC=", unique(rd$AUC))
    )
  }))

  roc_all$Model <- factor(roc_all$Model, levels = eval_results$model_names)

  ggplot2::ggplot(roc_all, ggplot2::aes(x = FPR, y = TPR,
                                        color = Model)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50") +
    ggplot2::geom_text(data = auc_labels,
                       ggplot2::aes(x = FPR, y = TPR, label = Label,
                                    color = Model),
                       size = 3.2, hjust = 0, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x        = "False Positive Rate (1 - Specificity)",
      y        = "True Positive Rate (Sensitivity)",
      color    = "Model",
      title    = "ROC Curves — Model Comparison",
      subtitle = paste0("Outcome: ", eval_results$outcome_var)
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40", size = 9),
      legend.position = "none"
    )
}


# ---- internal: SHAP plots (bar, beeswarm, or both) ----
.plot_eval_shap <- function(eval_results, shap_type, model_select,
                            top_n, shap_class) {

  if (!requireNamespace("shapviz", quietly = TRUE))
    stop("Package 'shapviz' is required for SHAP plots. ",
         "Install via install.packages('shapviz').")

  if (is.null(eval_results$shap))
    stop("No SHAP values found. Re-run perform_ModelEvaluation with ",
         "compute_shap = TRUE.")

  # Filter models
  if (!is.null(model_select)) {
    if (!model_select %in% eval_results$model_names)
      stop(paste0("model_select '", model_select, "' not found. ",
                  "Available: ",
                  paste(eval_results$model_names, collapse = ", ")))
    shap_subset <- eval_results$shap[model_select]
  } else {
    shap_subset <- eval_results$shap
  }

  # Remove NULLs (models where SHAP failed)
  shap_subset <- Filter(Negate(is.null), shap_subset)
  if (length(shap_subset) == 0)
    stop("No valid SHAP objects available to plot.")

  # Build one or two plots per model
  plots <- lapply(names(shap_subset), function(nm) {
    sv <- shap_subset[[nm]]

    # For multi-output shapviz (multiclass), select class
    if (inherits(sv, "mshapviz")) {
      cls_names <- names(sv)
      chosen    <- shap_class %||% cls_names[1]
      if (!chosen %in% cls_names) {
        warning(paste0("shap_class '", chosen, "' not found in '", nm,
                       "'. Using '", cls_names[1], "'."))
        chosen <- cls_names[1]
      }
      sv <- sv[[chosen]]
    }

    .make_shap_plot <- function(kind) {
      shapviz::sv_importance(sv, kind = kind, max_display = top_n) +
        ggplot2::labs(
          title    = paste0("SHAP — ", nm),
          subtitle = if (kind == "bar")
            paste0("Top ", top_n, " features by mean |SHAP|")
          else
            paste0("Top ", top_n,
                   " features — SHAP value distribution")
        ) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5,
                                                color = "grey40", size = 9)
        )
    }

    if (shap_type == "both") {
      list(bar      = .make_shap_plot("bar"),
           beeswarm = .make_shap_plot("beeswarm"))
    } else {
      .make_shap_plot(shap_type)
    }
  })
  names(plots) <- names(shap_subset)

  # Simplify return value when only one model and not "both"
  if (length(plots) == 1 && shap_type != "both") return(plots[[1]])
  if (length(plots) == 1 && shap_type == "both") return(plots[[1]])
  plots
}
