#' Perform Comparative Statistical Analysis on Preprocessed Metabolomics Data
#'
#' @description
#' Conducts comprehensive comparative statistical analysis on preprocessed metabolomics data
#' by iterating over each metabolite and applying the 'auto_compare' function.
#'
#' @details
#' This function serves as a wrapper for the \code{auto_compare} function, applying
#' it to each metabolite (column) in a metabolomics data matrix.
#'
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates input data structure and parameters.
#'   \item Removes quality control samples from analysis.
#'   \item For each metabolite, it combines the data with the group metadata.
#'   \item Calls \code{auto_compare} to select and perform the appropriate statistical test.
#'   \item Aggregates the results from all metabolites into a single data frame.
#'   \item Applies the specified p-value adjustment method *within* each call to \code{auto_compare}
#'         for post-hoc tests.
#'   \item Generates optional visualization plots using \code{ggstatsplot}.
#' }
#'
#' **NOTE:** This function requires the \code{auto_compare} function to be
#' loaded into the R environment (e.g., by sourcing \code{auto_compare.R}).
#'
#' @param data List. Output from \code{perform_PreprocessingPeakData} function containing:
#'   \itemize{
#'     \item \code{data_scaledNONPLS_varFiltered}: Numeric matrix of processed metabolite data
#'     \item \code{Metadata}: Data frame with sample metadata including 'Group' column
#'   }
#' @param adjust_p_method Character. Method for p-value adjustment, passed to
#'   \code{auto_compare}. Default is "BH". Options include:
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
#'   \code{auto_compare}). Default is \code{FALSE}.
#' @param plot_metabolites Character vector. Names of metabolites to visualize. If provided,
#'   generates statistical plots using \code{ggstatsplot}. Default is \code{NULL} (no plots).
#' @param alpha Numeric. Significance threshold passed to \code{auto_compare} for
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
#' # Assuming 'auto_compare.R' has been sourced
#' # source("auto_compare.R")
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
#' @seealso \code{\link{perform_PreprocessingPeakData}}, \code{auto_compare},
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

  # Parallel setup
  BPPARAM_to_use <- BiocParallel::SnowParam(workers = parallelly::availableCores(omit = 2))
  BiocParallel::register(BPPARAM_to_use)

  # Ensure SerialParam is restored even if the function errors
  on.exit({
    BiocParallel::register(BiocParallel::SerialParam())
  }, add = TRUE)

  # Check if auto_compare function exists
  if (!exists("auto_compare", mode = "function")) {
    stop("The 'auto_compare' function is not loaded. Please source the 'auto_compare.R' file.")
  }

  # Validate inputs
  if (!is.list(data) || !all(c("data_scaledNONPLS_varFiltered", "Metadata") %in% names(data))) {
    stop("Input 'data' must be a list with 'data_scaledNONPLS_varFiltered' and 'Metadata'")
  }
  if (!"Group" %in% colnames(data$Metadata)) {
    stop("Metadata must contain a 'Group' column")
  }
  if (nrow(data$data_scaledNONPLS_varFiltered) != nrow(data$Metadata)) {
    stop("Number of rows in data and metadata must match")
  }

  # Remove QC samples
  qc_patterns <- c("SQC", "EQC", "QC", "BLANK", "blank", "Blank")
  groups_col <- data$Metadata$Group
  qc_indices <- groups_col %in% qc_patterns
  df <- data$data_scaledNONPLS_varFiltered[!qc_indices, , drop = FALSE]
  groups <- factor(as.character(groups_col)[!qc_indices])

  # Validate group sizes
  group_counts <- table(groups)
  if (any(group_counts < min_group_size)) {
    stop("Groups with insufficient sample size: ",
         paste(names(group_counts)[group_counts < min_group_size], collapse = ", "))
  }
  if (length(group_counts) < 2) stop("At least 2 groups required")

  num_groups <- length(levels(groups))
  if (verbose) cat("Analyzing", ncol(df), "metabolites across", num_groups, "groups using auto_compare\n")

  # Analysis function (now uses auto_compare)
  analyze_metabolite <- function(x, metabolite_name) {
    tryCatch({
      # 1. Create the data frame required by auto_compare
      analysis_df <- data.frame(
        outcome = x,
        group = groups
      )

      # 2. Call auto_compare
      ac_result <- auto_compare(
        data = analysis_df,
        outcome = "outcome",
        group = "group",
        paired = paired,
        test_type = "auto",
        test_alpha = alpha,
        alpha_normality = alpha,
        alpha_variance = alpha,
        p_adjust_method = adjust_p_method,
        min_n_threshold = min_group_size,
        calculate_effect_size = TRUE,
        perform_posthoc = TRUE,
        verbose = FALSE, # Keep this FALSE to avoid spam
        num_cores = num_cores
      )

      # 3. Extract results into the expected format

      # Determine aggregate type based on test
      agg_type <- ifelse(ac_result$parametric, "mean", "median")

      # Extract aggregate stats
      if(agg_type == "mean") {
        agg_stats <- setNames(ac_result$data_summary$mean, ac_result$data_summary$group)
      } else {
        agg_stats <- setNames(ac_result$data_summary$median, ac_result$data_summary$group)
      }

      # Extract post-hoc results
      posthoc_df <- NULL
      posthoc_test_name <- NA_character_
      interpretation_text <- NA_character_

      if (!is.null(ac_result$posthoc_result) && nrow(ac_result$posthoc_result) > 0) {
        ph_res <- ac_result$posthoc_result

        # Get post-hoc p-values
        posthoc_df <- data.frame(
          comparison = paste(ph_res$group1, "vs", ph_res$group2),
          p_adj = ph_res$p.adj
        )

        # Get post-hoc test name
        posthoc_test_name <- ph_res$posthoc_test[1]

        # Get combined interpretation
        interpretation_text <- paste(ph_res$interpretation, collapse = "; ")
      }

      return(list(
        test_used = ac_result$test_used,
        omnibus_p = ac_result$test_result$p.value,
        agg_type = agg_type,
        agg_stats = agg_stats,
        effect_size = ac_result$effect_size$estimate,
        effect_size_metric = ac_result$effect_size$metric,
        posthoc = posthoc_df,
        posthoc_test_used = posthoc_test_name,   # NEW
        interpretation = interpretation_text  # NEW
      ))

    }, error = function(e) {
      if (verbose) cat("Error in", metabolite_name, ":", e$message, "\n")
      return(list(
        test_used = "Test Failed",
        omnibus_p = NA,
        agg_type = NA,
        agg_stats = NULL,
        effect_size = NA,
        effect_size_metric = NA,
        posthoc = NULL,
        posthoc_test_used = NA, # NEW
        interpretation = NA    # NEW
      ))
    })
  }

  # Run analysis
  results_list <- lapply(colnames(df), function(met) {
    if (verbose) cat(".") # Print progress
    analyze_metabolite(df[[met]], met)
  })
  if (verbose) cat("\nAnalysis complete. Formatting results.\n")
  names(results_list) <- colnames(df)

  # Format results
  unique_groups <- levels(groups)
  n_metabolites <- ncol(df)

  # Initialize result dataframe with consistent length vectors
  results_df <- data.frame(
    test_used = sapply(results_list, function(r) r$test_used),
    omnibus_p_value = sapply(results_list, function(r) ifelse(is.null(r$omnibus_p), NA, r$omnibus_p)),
    effect_size = sapply(results_list, function(r) ifelse(is.null(r$effect_size), NA, r$effect_size)),
    statistic_type = sapply(results_list, function(r) ifelse(is.null(r$agg_type), NA, r$agg_type)),
    effect_size_metric = sapply(results_list, function(r) ifelse(is.null(r$effect_size_metric), NA, r$effect_size_metric)),
    posthoc_test_used = sapply(results_list, function(r) ifelse(is.null(r$posthoc_test_used), NA, r$posthoc_test_used)),
    interpretation = sapply(results_list, function(r) ifelse(is.null(r$interpretation), NA, r$interpretation)),
    row.names = colnames(df),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Add aggregate columns - ensure consistent length
  for (group_name in unique_groups) {
    col_name <- paste0("aggregate_", group_name)
    results_df[[col_name]] <- sapply(results_list, function(r) {
      if (is.null(r$agg_stats) || is.null(r$agg_stats[[group_name]]) || is.na(r$agg_stats[[group_name]])) {
        return(NA)
      } else {
        return(r$agg_stats[[group_name]])
      }
    }, USE.NAMES = FALSE)
  }

  # Add post-hoc columns only if multi-group and significant
  if (num_groups > 2) {
    # Get all possible comparisons (handle case where no posthoc was done)
    all_comparisons <- unique(unlist(lapply(results_list, function(r) {
      if (!is.null(r$posthoc) && nrow(r$posthoc) > 0) {
        r$posthoc$comparison
      } else {
        character(0)
      }
    })))

    # Generate all pairwise combinations if no posthoc was done
    if (length(all_comparisons) == 0) {
      all_comparisons <- apply(combn(unique_groups, 2), 2, function(x) paste(x[1], "vs", x[2]))
    }

    for (comp in all_comparisons) {
      col_name <- paste0("posthoc_p_", comp) # Sanitize name
      results_df[[col_name]] <- sapply(results_list, function(r) {
        if (!is.null(r$posthoc) && nrow(r$posthoc) > 0 && comp %in% r$posthoc$comparison) {
          idx <- match(comp, r$posthoc$comparison)
          return(r$posthoc$p_adj[idx])
        } else {
          return(NA)
        }
      }, USE.NAMES = FALSE)
    }
  }

  # # Sort results
  # posthoc_cols <- grep("^posthoc_p_", colnames(results_df), value = TRUE)
  # if (sort_p) {
  #   # Sort by omnibus p-value if no post-hoc, or by min post-hoc p-value if they exist
  #   if (length(posthoc_cols) > 0) {
  #     min_posthoc_p <- apply(results_df[, posthoc_cols, drop = FALSE], 1, min, na.rm = TRUE)
  #     # Replace Inf from min(c(NA,NA), na.rm=T) with NA
  #     min_posthoc_p[is.infinite(min_posthoc_p)] <- NA
  #     results_df <- results_df[order(min_posthoc_p, results_df$omnibus_p_value, na.last = TRUE), ]
  #   } else {
  #     results_df <- results_df[order(results_df$omnibus_p_value, na.last = TRUE), ]
  #   }
  # }
  # Sort results
  posthoc_cols <- grep("^posthoc_p_", colnames(results_df), value = TRUE)
  if (sort_p) {
    # Sort by omnibus p-value if no post-hoc, or by min post-hoc p-value if they exist
    if (length(posthoc_cols) > 0) {
      # Suppress warnings and handle all-NA rows
      min_posthoc_p <- suppressWarnings(
        apply(results_df[, posthoc_cols, drop = FALSE], 1, function(row) {
          if (all(is.na(row))) {
            return(NA)
          } else {
            return(min(row, na.rm = TRUE))
          }
        })
      )
      results_df <- results_df[order(min_posthoc_p, results_df$omnibus_p_value, na.last = TRUE), ]
    } else {
      results_df <- results_df[order(results_df$omnibus_p_value, na.last = TRUE), ]
    }
  }

  # Generate summary
  summary_stats <- list(
    total_metabolites = nrow(results_df),
    valid_results = sum(!is.na(results_df$omnibus_p_value)),
    significant_omnibus = sum(results_df$omnibus_p_value < 0.05, na.rm = TRUE),
    significant_posthoc = if (length(posthoc_cols) > 0) {
      sum(apply(results_df[, posthoc_cols, drop = FALSE], 1, function(row) any(row < 0.05, na.rm = TRUE)), na.rm = TRUE)
    } else 0,
    test_distribution = table(results_df$test_used),
    posthoc_distribution = table(results_df$posthoc_test_used),
    group_sizes = table(groups)
  )

  # Generate plots if requested
  plots_list <- list()
  if (!is.null(plot_metabolites) && requireNamespace("ggstatsplot", quietly = TRUE)) {
    if (verbose) cat("Generating plots...\n")
    for (met in plot_metabolites) {
      if (met %in% colnames(df)) {
        plot_data <- data.frame(
          group = groups[complete.cases(df[[met]])],
          value = df[[met]][complete.cases(df[[met]])]
        )
        # Use the test type decided by auto_compare
        test_type <- if (grepl("t-test|ANOVA", results_df[met, "test_used"])) "parametric" else "nonparametric"

        plots_list[[met]] <- ggstatsplot::ggbetweenstats(
          plot_data, x = group, y = value, type = test_type,
          title = paste("Comparative Analysis:", met),
          p.adjust.method = adjust_p_method,
          pairwise.comparisons = TRUE
        )
      } else {
        warning(paste("Metabolite", met, "not found for plotting."))
      }
    }
  }

  # return(list(
  #   results = results_df,
  #   plots = plots_list,
  #   summary = summary_stats,
  #   metadata = list(
  #     function_origin = "perform_ComparativeAnalysis",
  #     timestamp = Sys.time(),
  #     num_groups = num_groups,
  #     num_metabolites = ncol(df),
  #     num_samples = nrow(df),
  #     adjust_p_method = adjust_p_method,
  #     paired = paired,
  #     alpha = alpha
  #   )
  # ))

  out <- list(
    results           = results_df,
    plots             = plots_list,
    summary           = summary_stats,
    metadata          = list(
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


# S3 Methods
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
  # Get top 5 significant features
  top_feats <- head(object$results[order(object$results$omnibus_p_value),
                                   c("omnibus_p_value", "test_used", "interpretation")], 5)

  ans <- list(
    stats        = object$summary,
    top_features = top_feats,
    meta         = object$metadata
  )
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
  if(x$stats$significant_posthoc > 0) {
    cat("Sig. (Post-hoc):  ", x$stats$significant_posthoc, "\n")
  }

  cat("\n-- Test Distribution --\n")
  print(x$stats$test_distribution)

  cat("\n-- Top 5 Significant Features --\n")
  print(x$top_features)
  invisible(x)
}
