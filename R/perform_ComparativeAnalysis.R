#' Perform Comparative Statistical Analysis on Preprocessed Metabolomics Data
#'
#' @description
#' Conducts comprehensive comparative statistical analysis on preprocessed metabolomics data
#' by iterating over each metabolite and applying the 'run_diff' function.
#'
#' @details
#' This function serves as a wrapper for the \code{run_diff} function, applying
#' it to each metabolite (column) in a metabolomics data matrix.
#'
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates input data structure and parameters.
#'   \item Removes quality control samples from analysis.
#'   \item For each metabolite, it combines the data with the group metadata.
#'   \item Calls \code{run_diff} to select and perform the appropriate statistical test.
#'   \item Aggregates the results from all metabolites into a single data frame.
#'   \item Applies the specified p-value adjustment method *within* each call to \code{run_diff}
#'         for post-hoc tests.
#'   \item Generates optional visualization plots using \code{ggstatsplot}.
#' }
#'
#' **NOTE:** This function requires the \code{run_diff} function to be
#' loaded into the R environment (e.g., by sourcing \code{run_diff.R}).
#'
#' @param data List. Output from \code{perform_PreprocessingPeakData} function containing:
#'   \itemize{
#'     \item \code{data_scaledNONPLS_varFiltered}: Numeric matrix of processed metabolite data
#'     \item \code{Metadata}: Data frame with sample metadata including 'Group' column
#'   }
#' @param adjust_p_method Character. Method for p-value adjustment, passed to
#'   \code{run_diff}. Default is "BH". Options include:
#'   \itemize{
#'     \item "holm": Holm (1979) - Controls family-wise error rate
#'     \item "hochberg": Hochberg (1988) - Less conservative than Bonferroni
#'     \item "hommel": Hommel (1988) - More powerful than Hochberg
#'     \item "bonferroni": Classical Bonferroni correction
#'     \item "BH": Benjamini & Hochberg (1995) - Controls false discovery rate
#'     \item "BY": Benjamini & Yekutieli (2001) - More conservative FDR control
#'     \item "fdr": Alias for "BH"
#'     \item "none": No adjustment
#'   }
#' @param sort_p Logical. If \code{TRUE} (default), results are sorted by the
#'   adjusted post-hoc p-values in ascending order.
#' @param paired Logical. If \code{TRUE}, performs paired sample tests (passed to
#'   \code{run_diff}). Default is \code{FALSE}.
#' @param plot_metabolites Character vector. Names of metabolites to visualize. If provided,
#'   generates statistical plots using \code{ggstatsplot}. Default is \code{NULL} (no plots).
#' @param alpha Numeric. Significance threshold passed to \code{run_diff} for
#'   assumption tests (\code{alpha_normality}, \code{alpha_variance}) and the main
#'   test (\code{test_alpha}). Default is 0.05.
#' @param min_group_size Integer. Minimum required sample size per group. Default is 3.
#' @param verbose Logical. If \code{TRUE}, prints detailed progress information. Default is \code{FALSE}.
#' @param num_cores Integer or "max" specifying the number of cores to be used.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{results}: Data frame with statistical test results for each metabolite.
#'       Includes \code{posthoc_test_used} and a combined \code{interpretation} column.
#'     \item \code{plots}: List of ggplot objects (if \code{plot_metabolites} specified)
#'     \item \code{summary}: Summary statistics of the analysis
#'     \item \code{metadata}: Analysis metadata and parameters used
#'   }
#'
#' @examples
#' \dontrun{
#' # Assuming 'run_diff.R' has been sourced
#' # source("run_diff.R")
#'
#' # Basic two-group comparison
#' results <- perform_ComparativeAnalysis(
#'   data = preprocessed_data,
#'   adjust_p_method = "BH"
#' )
#'
#' # Paired comparison with plots
#' results <- perform_ComparativeAnalysis(
#'   data = preprocessed_data,
#'   paired = TRUE,
#'   plot_metabolites = c("metabolite_1", "metabolite_2"),
#'   verbose = TRUE
#' )
#'
#' # Multi-group comparison with strict correction
#' results <- perform_ComparativeAnalysis(
#'   data = preprocessed_data,
#'   adjust_p_method = "bonferroni",
#'   sort_p = TRUE
#' )
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{perform_PreprocessingPeakData}}, \code{run_diff},
#'   \code{\link[ggstatsplot]{ggbetweenstats}}
#'
#' @importFrom BiocParallel SnowParam SerialParam
#' @importFrom parallelly availableCores
#'
#' @export
#'
perform_ComparativeAnalysis <- function(
    data,
    adjust_p_method  = "BH",
    sort_p           = TRUE,
    paired           = FALSE,
    plot_metabolites = NULL,
    alpha            = 0.05,
    min_group_size   = 3,
    verbose          = FALSE,
    num_cores        = 1
) {

  # ── Restore serial param on exit ──────────────────────────────────────────
  on.exit(BiocParallel::register(BiocParallel::SerialParam()), add = TRUE)

  # Only register a parallel backend when the user actually wants > 1 core,
  # otherwise leave the default SerialParam in place.
  if (num_cores != 1) {
    avail <- max(1L, parallelly::availableCores(omit = 2L))
    workers <- if (identical(num_cores, "max")) avail else min(as.integer(num_cores), avail)
    BiocParallel::register(BiocParallel::SnowParam(workers = workers))
  }

  # ── Dependency check ──────────────────────────────────────────────────────
  if (!exists("run_diff", mode = "function")) {
    stop("The 'run_diff' function is not loaded. Please source the 'run_diff.R' file.")
  }

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.list(data) ||
      !all(c("data_scaledNONPLS_varFiltered", "Metadata") %in% names(data))) {
    stop("Input 'data' must be a list with 'data_scaledNONPLS_varFiltered' and 'Metadata'.")
  }
  if (!"Group" %in% colnames(data$Metadata)) {
    stop("Metadata must contain a 'Group' column.")
  }
  if (nrow(data$data_scaledNONPLS_varFiltered) != nrow(data$Metadata)) {
    stop("Number of rows in data and metadata must match.")
  }

  # ── Remove QC samples ─────────────────────────────────────────────────────
  qc_patterns <- c("SQC", "EQC", "QC", "BLANK", "blank", "Blank")
  groups_col  <- data$Metadata$Group
  keep        <- !groups_col %in% qc_patterns
  df          <- data$data_scaledNONPLS_varFiltered[keep, , drop = FALSE]
  groups      <- factor(as.character(groups_col)[keep])

  # ── Validate group sizes ───────────────────────────────────────────────────
  group_counts <- table(groups)
  if (any(group_counts < min_group_size)) {
    stop("Groups with insufficient sample size: ",
         paste(names(group_counts)[group_counts < min_group_size], collapse = ", "))
  }
  if (length(group_counts) < 2) stop("At least 2 groups required.")

  num_groups <- length(levels(groups))
  if (verbose) {
    cat("Analyzing", ncol(df), "metabolites across",
        num_groups, "groups using run_diff\n")
  }

  # ── Per-metabolite analysis ───────────────────────────────────────────────
  analyze_metabolite <- function(met_name) {
    x_vec <- df[[met_name]]

    # Build the data frame run_diff expects
    analysis_df <- data.frame(outcome = x_vec, group = groups,
                              stringsAsFactors = FALSE)

    # ── run_diff call — note first arg is 'x', NOT 'data' ──────────────────
    ac_result <- tryCatch(
      run_diff(
        x                     = analysis_df,   # FIX: was `data =`
        outcome               = "outcome",
        group                 = "group",
        paired                = paired,
        test_type             = "auto",
        test_alpha            = alpha,
        alpha_normality       = alpha,
        alpha_variance        = alpha,
        p_adjust_method       = adjust_p_method,
        min_n_threshold       = min_group_size,
        calculate_effect_size = TRUE,
        perform_posthoc       = TRUE,
        verbose               = FALSE,
        num_cores             = 1L            # avoid nested parallelism
      ),
      error = function(e) {
        if (verbose) cat("  [ERROR]", met_name, "->", conditionMessage(e), "\n")
        NULL
      }
    )

    # Return a consistent failure record if run_diff errored
    if (is.null(ac_result)) {
      return(list(
        test_used        = "Test Failed",
        omnibus_p        = NA_real_,
        agg_type         = NA_character_,
        agg_stats        = NULL,
        effect_size      = NA_real_,
        effect_size_metric = NA_character_,
        posthoc          = NULL,
        posthoc_test_used = NA_character_,
        interpretation   = NA_character_
      ))
    }

    # ── Extract aggregates ─────────────────────────────────────────────────
    agg_type <- ifelse(isTRUE(ac_result$parametric), "mean", "median")
    ds       <- ac_result$data_summary   # always a data.frame with col 'group'

    agg_stats <- if (agg_type == "mean") {
      setNames(as.numeric(ds$mean),   as.character(ds$group))
    } else {
      setNames(as.numeric(ds$median), as.character(ds$group))
    }

    # ── Extract post-hoc ──────────────────────────────────────────────────
    posthoc_df        <- NULL
    posthoc_test_name <- NA_character_
    interpretation_tx <- NA_character_

    ph <- ac_result$posthoc_result
    if (!is.null(ph) && nrow(ph) > 0) {
      posthoc_df <- data.frame(
        comparison = paste(ph$group1, "vs", ph$group2),
        p_adj      = as.numeric(ph$p.adj),
        stringsAsFactors = FALSE
      )
      posthoc_test_name <- as.character(ph$posthoc_test[1])
      interpretation_tx <- paste(ph$interpretation, collapse = "; ")
    }

    # ── Extract omnibus p-value robustly ──────────────────────────────────
    omnibus_p <- tryCatch(
      as.numeric(ac_result$test_result$p.value),
      error = function(e) NA_real_
    )

    list(
      test_used          = ac_result$test_used,
      omnibus_p          = omnibus_p,
      agg_type           = agg_type,
      agg_stats          = agg_stats,
      effect_size        = if (!is.null(ac_result$effect_size)) ac_result$effect_size$estimate else NA_real_,
      effect_size_metric = if (!is.null(ac_result$effect_size)) ac_result$effect_size$metric   else NA_character_,
      posthoc            = posthoc_df,
      posthoc_test_used  = posthoc_test_name,
      interpretation     = interpretation_tx
    )
  }

  # Run over all metabolites
  met_names    <- colnames(df)
  results_list <- lapply(met_names, function(m) {
    if (verbose) cat(".")
    analyze_metabolite(m)
  })
  if (verbose) cat("\nAnalysis complete. Formatting results.\n")
  names(results_list) <- met_names

  # ── Build results data frame ──────────────────────────────────────────────
  unique_groups <- levels(groups)

  safe_val <- function(r, field, default = NA) {
    v <- r[[field]]
    if (is.null(v) || length(v) == 0) default else v
  }

  results_df <- data.frame(
    test_used          = sapply(results_list, safe_val, "test_used",          NA_character_),
    omnibus_p_value    = sapply(results_list, safe_val, "omnibus_p",          NA_real_),
    effect_size        = sapply(results_list, safe_val, "effect_size",        NA_real_),
    statistic_type     = sapply(results_list, safe_val, "agg_type",           NA_character_),
    effect_size_metric = sapply(results_list, safe_val, "effect_size_metric", NA_character_),
    posthoc_test_used  = sapply(results_list, safe_val, "posthoc_test_used",  NA_character_),
    interpretation     = sapply(results_list, safe_val, "interpretation",     NA_character_),
    row.names          = met_names,
    stringsAsFactors   = FALSE,
    check.names        = FALSE
  )

  # Aggregate columns per group
  for (grp_name in unique_groups) {
    col_name <- paste0("aggregate_", grp_name)
    results_df[[col_name]] <- sapply(results_list, function(r) {
      agg <- r$agg_stats
      if (is.null(agg) || is.na(agg[grp_name])) NA_real_ else as.numeric(agg[grp_name])
    }, USE.NAMES = FALSE)
  }

  # Post-hoc p-value columns (multi-group only)
  if (num_groups > 2) {
    all_comparisons <- unique(unlist(lapply(results_list, function(r) {
      if (!is.null(r$posthoc) && nrow(r$posthoc) > 0) r$posthoc$comparison
    })))

    if (length(all_comparisons) == 0) {
      all_comparisons <- apply(combn(unique_groups, 2), 2,
                               function(x) paste(x[1], "vs", x[2]))
    }

    for (comp in all_comparisons) {
      col_name <- paste0("posthoc_p_", comp)
      results_df[[col_name]] <- sapply(results_list, function(r) {
        ph <- r$posthoc
        if (!is.null(ph) && nrow(ph) > 0 && comp %in% ph$comparison) {
          as.numeric(ph$p_adj[match(comp, ph$comparison)])
        } else {
          NA_real_
        }
      }, USE.NAMES = FALSE)
    }
  }

  # ── Sort ──────────────────────────────────────────────────────────────────
  posthoc_cols <- grep("^posthoc_p_", colnames(results_df), value = TRUE)
  if (sort_p) {
    if (length(posthoc_cols) > 0) {
      min_ph_p <- suppressWarnings(
        apply(results_df[, posthoc_cols, drop = FALSE], 1, function(row) {
          if (all(is.na(row))) NA_real_ else min(row, na.rm = TRUE)
        })
      )
      results_df <- results_df[order(min_ph_p, results_df$omnibus_p_value, na.last = TRUE), ]
    } else {
      results_df <- results_df[order(results_df$omnibus_p_value, na.last = TRUE), ]
    }
  }

  # ── Summary ───────────────────────────────────────────────────────────────
  summary_stats <- list(
    total_metabolites  = nrow(results_df),
    valid_results      = sum(!is.na(results_df$omnibus_p_value)),
    significant_omnibus = sum(results_df$omnibus_p_value < 0.05, na.rm = TRUE),
    significant_posthoc = if (length(posthoc_cols) > 0) {
      sum(apply(results_df[, posthoc_cols, drop = FALSE], 1,
                function(row) any(row < 0.05, na.rm = TRUE)), na.rm = TRUE)
    } else 0L,
    test_distribution  = table(results_df$test_used),
    posthoc_distribution = table(results_df$posthoc_test_used),
    group_sizes        = table(groups)
  )

  # ── Optional plots ────────────────────────────────────────────────────────
  plots_list <- list()
  if (!is.null(plot_metabolites) && requireNamespace("ggstatsplot", quietly = TRUE)) {
    if (verbose) cat("Generating plots...\n")
    for (met in plot_metabolites) {
      if (!met %in% colnames(df)) {
        warning("Metabolite '", met, "' not found for plotting.")
        next
      }
      ok        <- complete.cases(df[[met]])
      plot_data <- data.frame(group = groups[ok], value = df[[met]][ok])
      t_type    <- if (grepl("t-test|ANOVA", results_df[met, "test_used"],
                             ignore.case = TRUE)) "parametric" else "nonparametric"
      plots_list[[met]] <- ggstatsplot::ggbetweenstats(
        plot_data, x = group, y = value, type = t_type,
        title             = paste("Comparative Analysis:", met),
        p.adjust.method   = adjust_p_method,
        pairwise.comparisons = TRUE
      )
    }
  }

  # ── Return ────────────────────────────────────────────────────────────────
  out <- list(
    results  = results_df,
    plots    = plots_list,
    summary  = summary_stats,
    metadata = list(
      function_origin = "perform_ComparativeAnalysis",
      timestamp       = Sys.time(),
      num_groups      = num_groups,
      num_metabolites = ncol(df),
      num_samples     = nrow(df),
      adjust_p_method = adjust_p_method,
      paired          = paired,
      alpha           = alpha
    )
  )
  class(out) <- c("perform_ComparativeAnalysis", "list")
  return(out)
}


# ── S3 Methods ────────────────────────────────────────────────────────────────

#' @export
print.perform_ComparativeAnalysis <- function(x, ...) {
  cat("=== Comparative Statistical Analysis ===\n")
  cat("Metabolites Tested: ", nrow(x$results), "\n")
  cat("Significant (p<0.05):", x$summary$significant_omnibus, "\n")
  cat("Adjustment Method:  ", x$metadata$adjust_p_method, "\n")
  cat("Groups Compared:    ", x$metadata$num_groups, "\n")
  invisible(x)
}

#' @export
summary.perform_ComparativeAnalysis <- function(object, ...) {
  top_feats <- head(
    object$results[order(object$results$omnibus_p_value),
                   c("omnibus_p_value", "test_used", "interpretation")],
    5
  )
  ans <- list(stats = object$summary, top_features = top_feats, meta = object$metadata)
  class(ans) <- "summary.perform_ComparativeAnalysis"
  return(ans)
}

#' @export
print.summary.perform_ComparativeAnalysis <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("Comparative Analysis Summary\n")
  cat("---------------------------------------\n")
  cat("Total Metabolites:", x$stats$total_metabolites, "\n")
  cat("Sig. (Omnibus):   ", x$stats$significant_omnibus, "\n")
  if (x$stats$significant_posthoc > 0)
    cat("Sig. (Post-hoc):  ", x$stats$significant_posthoc, "\n")
  cat("\n-- Test Distribution --\n")
  print(x$stats$test_distribution)
  cat("\n-- Top 5 Significant Features --\n")
  print(x$top_features)
  invisible(x)
}