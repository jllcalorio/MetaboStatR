#' Perform Comparative Statistical Analysis on Preprocessed Metabolomics Data
#'
#' @description
#' Conducts comprehensive comparative statistical analysis on preprocessed metabolomics data.
#' The function automatically selects appropriate statistical tests based on data characteristics
#' including normality, variance homogeneity, and sample independence. Supports both two-group
#' and multi-group comparisons with parametric and non-parametric alternatives.
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates input data structure and parameters
#'   \item Removes quality control samples from analysis
#'   \item Tests statistical assumptions (normality, variance homogeneity)
#'   \item Selects appropriate statistical tests automatically
#'   \item Applies multiple comparison corrections
#'   \item Generates optional visualization plots
#' }
#'
#' For two-group comparisons, the function chooses between:
#' \itemize{
#'   \item Paired/Independent t-test (parametric)
#'   \item Welch's t-test (unequal variances)
#'   \item Mann-Whitney U test (non-parametric)
#'   \item Wilcoxon signed-rank test (paired non-parametric)
#' }
#'
#' For multi-group comparisons:
#' \itemize{
#'   \item One-way ANOVA (parametric)
#'   \item Repeated measures ANOVA (paired)
#'   \item Kruskal-Wallis test (non-parametric)
#' }
#'
#' @param data List. Output from \code{perform_PreprocessingPeakData} function containing:
#'   \itemize{
#'     \item \code{data_scaledPCA_rsdFiltered_varFiltered}: Numeric matrix of processed metabolite data
#'     \item \code{Metadata}: Data frame with sample metadata including 'Group' column
#'   }
#' @param adjust_p_method Character. Method for p-value adjustment. Default is "BH".
#'   Options include:
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
#' @param sort_p Logical. If \code{TRUE} (default), results are sorted by adjusted p-values
#'   in ascending order.
#' @param paired Logical. If \code{TRUE}, performs paired sample tests. Default is \code{FALSE}.
#'   Note: Requires equal group sizes for multi-group comparisons.
#' @param plot_metabolites Character vector. Names of metabolites to visualize. If provided,
#'   generates statistical plots using \code{ggstatsplot}. Default is \code{NULL} (no plots).
#' @param alpha Numeric. Significance threshold for assumption tests. Default is 0.05.
#' @param min_group_size Integer. Minimum required sample size per group. Default is 3.
#' @param verbose Logical. If \code{TRUE}, prints detailed progress information. Default is \code{FALSE}.
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{results}: Data frame with statistical test results for each metabolite
#'     \item \code{plots}: List of ggplot objects (if \code{plot_metabolites} specified)
#'     \item \code{summary}: Summary statistics of the analysis
#'     \item \code{assumptions}: Results of assumption tests
#'     \item \code{metadata}: Analysis metadata and parameters used
#'   }
#'
#' @examples
#' \dontrun{
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
#' @author Your Name
#' @seealso \code{\link{perform_PreprocessingPeakData}}, \code{\link[stats]{t.test}},
#'   \code{\link[stats]{aov}}, \code{\link[ggstatsplot]{ggbetweenstats}}
#' @export
#'
perform_ComparativeAnalysis2 <- function(
    data,
    adjust_p_method = "BH",
    sort_p = TRUE,
    paired = FALSE,
    plot_metabolites = NULL,
    alpha = 0.05,
    min_group_size = 3,
    verbose = FALSE
) {

  # Validate inputs
  .validate_inputs(data, adjust_p_method, alpha, min_group_size)

  # Initialize results structure
  results_list <- .initialize_results()

  if (verbose) cat("Starting comparative analysis...\n")

  # Extract and validate data
  analysis_data <- .extract_analysis_data(data, min_group_size, verbose)
  df <- analysis_data$df
  groups <- analysis_data$groups
  num_groups <- length(unique(groups))

  if (verbose) cat("Number of groups detected:", num_groups, "\n")

  # Validate paired analysis requirements
  if (paired) {
    .validate_paired_analysis(groups, num_groups)
  }

  # Perform statistical analysis
  if (num_groups == 2) {
    if (verbose) cat("Performing two-group comparison...\n")
    stat_results <- .perform_two_group_analysis(df, groups, paired, alpha, verbose)
  } else {
    if (verbose) cat("Performing multi-group comparison...\n")
    stat_results <- .perform_multi_group_analysis(df, groups, paired, alpha, verbose)
  }

  # Process results
  final_results <- .process_results(stat_results, adjust_p_method, sort_p, verbose)

  # Store main results
  results_list$results <- final_results$results
  results_list$assumptions <- final_results$assumptions
  results_list$summary <- .generate_summary(final_results$results, groups)

  # Generate plots if requested
  if (!is.null(plot_metabolites)) {
    if (verbose) cat("Generating plots for", length(plot_metabolites), "metabolites...\n")
    results_list$plots <- .generate_plots(
      plot_metabolites, df, groups, final_results$results,
      adjust_p_method, verbose
    )
  }

  # Store metadata
  results_list$metadata <- list(
    function_origin = "perform_ComparativeAnalysis",
    timestamp = Sys.time(),
    num_groups = num_groups,
    num_metabolites = ncol(df),
    num_samples = nrow(df),
    adjust_p_method = adjust_p_method,
    paired = paired,
    alpha = alpha,
    r_version = R.version.string
  )

  if (verbose) cat("Analysis completed successfully!\n")

  return(results_list)
}

# Helper function: Validate inputs
.validate_inputs <- function(data, adjust_p_method, alpha, min_group_size) {
  # Check data structure
  if (!is.list(data)) {
    stop("Input 'data' must be a list from perform_PreprocessingPeakData function")
  }

  required_elements <- c("data_scaledPCA_rsdFiltered_varFiltered", "Metadata")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    stop("Missing required data elements: ", paste(missing_elements, collapse = ", "))
  }

  # Check p-value adjustment method
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!adjust_p_method %in% valid_methods) {
    stop("Invalid adjust_p_method. Must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Check numeric parameters
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }

  if (!is.numeric(min_group_size) || min_group_size < 2) {
    stop("min_group_size must be a numeric value >= 2")
  }

  # Check if required packages are available for plotting
  if (!is.null(plot_metabolites) && !requireNamespace("ggstatsplot", quietly = TRUE)) {
    warning("ggstatsplot package not available. Plots will not be generated.")
  }
}

# Helper function: Initialize results structure
.initialize_results <- function() {
  list(
    results = NULL,
    plots = list(),
    summary = NULL,
    assumptions = NULL,
    metadata = NULL
  )
}

# Helper function: Extract and validate analysis data
.extract_analysis_data <- function(data, min_group_size, verbose) {
  # Identify QC samples
  qc_patterns <- c("SQC", "EQC", "QC", "BLANK", "blank", "Blank")
  qc_indices <- data$Metadata$Group %in% qc_patterns

  if (verbose && sum(qc_indices) > 0) {
    cat("Removing", sum(qc_indices), "QC samples from analysis\n")
  }

  # Extract non-QC data
  non_qc_indices <- !qc_indices
  df <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, , drop = FALSE]
  groups <- factor(data$Metadata$Group[non_qc_indices])

  # Validate data dimensions
  if (nrow(df) == 0) {
    stop("No samples remaining after removing QC samples")
  }

  if (ncol(df) == 0) {
    stop("No metabolites found in the data")
  }

  # Check group sizes
  group_counts <- table(groups)
  small_groups <- names(group_counts)[group_counts < min_group_size]

  if (length(small_groups) > 0) {
    stop("Groups with insufficient sample size (< ", min_group_size, "): ",
         paste(small_groups, collapse = ", "))
  }

  # Check for missing values
  if (any(is.na(df))) {
    warning("Missing values detected in data. Results may be unreliable.")
  }

  return(list(df = df, groups = groups))
}

# Helper function: Validate paired analysis requirements
.validate_paired_analysis <- function(groups, num_groups) {
  group_counts <- table(groups)

  if (num_groups > 2 && !all(group_counts == group_counts[1])) {
    stop("Paired analysis requires equal sample sizes across all groups")
  }

  if (num_groups == 2 && group_counts[1] != group_counts[2]) {
    stop("Paired analysis requires equal sample sizes in both groups")
  }
}

# Helper function: Perform two-group analysis
.perform_two_group_analysis <- function(df, groups, paired, alpha, verbose) {
  unique_groups <- levels(groups)
  idx_group1 <- which(groups == unique_groups[1])
  idx_group2 <- which(groups == unique_groups[2])

  # Pre-allocate results matrix for efficiency
  n_metabolites <- ncol(df)
  results_matrix <- matrix(NA, nrow = n_metabolites, ncol = 5,
                           dimnames = list(colnames(df),
                                           c("test_used", "p_value", "normality_p",
                                             "variance_p", "effect_size")))

  # Vectorized computation where possible
  for (i in seq_len(n_metabolites)) {
    x <- df[, i]
    x1 <- x[idx_group1]
    x2 <- x[idx_group2]

    # Test assumptions
    normality_test <- .test_normality(x, groups)
    variance_test <- if (!paired) .test_variance_equality(x1, x2) else list(p_value = 1)

    # Select and perform appropriate test
    test_result <- .select_and_perform_test_two_group(
      x1, x2, normality_test$is_normal, variance_test$p_value > alpha, paired
    )

    # Calculate effect size
    effect_size <- .calculate_effect_size_two_group(x1, x2, paired)

    results_matrix[i, ] <- c(
      test_result$test_used,
      test_result$p_value,
      normality_test$p_value,
      variance_test$p_value,
      effect_size
    )
  }

  return(list(
    results = results_matrix,
    group_info = list(groups = unique_groups, indices = list(idx_group1, idx_group2))
  ))
}

# Helper function: Perform multi-group analysis
.perform_multi_group_analysis <- function(df, groups, paired, alpha, verbose) {
  n_metabolites <- ncol(df)
  results_matrix <- matrix(NA, nrow = n_metabolites, ncol = 5,
                           dimnames = list(colnames(df),
                                           c("test_used", "p_value", "normality_p",
                                             "variance_p", "effect_size")))

  for (i in seq_len(n_metabolites)) {
    x <- df[, i]

    # Test assumptions
    normality_test <- .test_normality_multi_group(x, groups)
    variance_test <- if (!paired) .test_variance_homogeneity(x, groups) else list(p_value = 1)

    # Select and perform appropriate test
    test_result <- .select_and_perform_test_multi_group(
      x, groups, normality_test$all_normal, variance_test$p_value > alpha, paired
    )

    # Calculate effect size
    effect_size <- .calculate_effect_size_multi_group(x, groups)

    results_matrix[i, ] <- c(
      test_result$test_used,
      test_result$p_value,
      normality_test$p_value,
      variance_test$p_value,
      effect_size
    )
  }

  return(list(
    results = results_matrix,
    group_info = list(groups = levels(groups))
  ))
}

# Helper function: Test normality for two-group analysis
.test_normality <- function(x, groups) {
  # Test residuals from linear model for normality
  tryCatch({
    residuals <- residuals(lm(x ~ groups))
    if (length(residuals) < 3) {
      return(list(is_normal = FALSE, p_value = 0))
    }
    shapiro_result <- shapiro.test(residuals)
    list(is_normal = shapiro_result$p.value > 0.05, p_value = shapiro_result$p.value)
  }, error = function(e) {
    list(is_normal = FALSE, p_value = 0)
  })
}

# Helper function: Test normality for multi-group analysis
.test_normality_multi_group <- function(x, groups) {
  tryCatch({
    group_normality <- tapply(x, groups, function(g) {
      if (length(g) < 3) return(0)
      shapiro.test(g)$p.value
    })

    min_p <- min(group_normality, na.rm = TRUE)
    all_normal <- all(group_normality > 0.05, na.rm = TRUE)

    list(all_normal = all_normal, p_value = min_p)
  }, error = function(e) {
    list(all_normal = FALSE, p_value = 0)
  })
}

# Helper function: Test variance equality (two groups)
.test_variance_equality <- function(x1, x2) {
  tryCatch({
    var_test <- var.test(x1, x2)
    list(p_value = var_test$p.value)
  }, error = function(e) {
    list(p_value = 0)
  })
}

# Helper function: Test variance homogeneity (multi-group)
.test_variance_homogeneity <- function(x, groups) {
  tryCatch({
    bartlett_test <- bartlett.test(x ~ groups)
    list(p_value = bartlett_test$p.value)
  }, error = function(e) {
    list(p_value = 0)
  })
}

# Helper function: Select and perform appropriate test (two groups)
.select_and_perform_test_two_group <- function(x1, x2, is_normal, equal_variance, paired) {
  tryCatch({
    if (paired) {
      if (is_normal) {
        result <- t.test(x1, x2, paired = TRUE)
        list(test_used = "Paired t-test", p_value = result$p.value)
      } else {
        result <- wilcox.test(x1, x2, paired = TRUE)
        list(test_used = "Wilcoxon Signed-Rank test", p_value = result$p.value)
      }
    } else {
      if (is_normal) {
        if (equal_variance) {
          result <- t.test(x1, x2, var.equal = TRUE)
          list(test_used = "Independent Samples t-test", p_value = result$p.value)
        } else {
          result <- t.test(x1, x2, var.equal = FALSE)
          list(test_used = "Welch's t-test", p_value = result$p.value)
        }
      } else {
        result <- wilcox.test(x1, x2)
        list(test_used = "Mann-Whitney U test", p_value = result$p.value)
      }
    }
  }, error = function(e) {
    list(test_used = "Test Failed", p_value = 1)
  })
}

# Helper function: Select and perform appropriate test (multi-group)
.select_and_perform_test_multi_group <- function(x, groups, all_normal, equal_variance, paired) {
  tryCatch({
    if (paired) {
      # For paired multi-group, use repeated measures ANOVA
      subject_id <- rep(1:length(unique(groups)), length(levels(groups)))
      aov_result <- aov(x ~ groups + Error(factor(subject_id)))
      p_value <- summary(aov_result)$`Error: Within`[[1]]$`Pr(>F)`[1]
      list(test_used = "Repeated Measures ANOVA", p_value = p_value)
    } else {
      if (all_normal && equal_variance) {
        aov_result <- aov(x ~ groups)
        p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
        list(test_used = "One-way ANOVA", p_value = p_value)
      } else {
        kw_result <- kruskal.test(x ~ groups)
        list(test_used = "Kruskal-Wallis test", p_value = kw_result$p.value)
      }
    }
  }, error = function(e) {
    list(test_used = "Test Failed", p_value = 1)
  })
}

# Helper function: Calculate effect size (two groups)
.calculate_effect_size_two_group <- function(x1, x2, paired) {
  tryCatch({
    if (paired) {
      # Cohen's d for paired samples
      diff <- x1 - x2
      cohen_d <- mean(diff) / sd(diff)
    } else {
      # Cohen's d for independent samples
      pooled_sd <- sqrt(((length(x1) - 1) * var(x1) + (length(x2) - 1) * var(x2)) /
                          (length(x1) + length(x2) - 2))
      cohen_d <- (mean(x1) - mean(x2)) / pooled_sd
    }
    return(abs(cohen_d))
  }, error = function(e) {
    return(0)
  })
}

# Helper function: Calculate effect size (multi-group)
.calculate_effect_size_multi_group <- function(x, groups) {
  tryCatch({
    # Eta-squared
    aov_result <- aov(x ~ groups)
    ss_total <- sum((x - mean(x))^2)
    ss_between <- sum(summary(aov_result)[[1]][["Sum Sq"]])
    eta_squared <- ss_between / ss_total
    return(eta_squared)
  }, error = function(e) {
    return(0)
  })
}

# Helper function: Process results
.process_results <- function(stat_results, adjust_p_method, sort_p, verbose) {
  results_df <- as.data.frame(stat_results$results, stringsAsFactors = FALSE)

  # Convert numeric columns
  numeric_cols <- c("p_value", "normality_p", "variance_p", "effect_size")
  results_df[numeric_cols] <- lapply(results_df[numeric_cols], as.numeric)

  # Apply p-value adjustment
  results_df$adj_p_value <- p.adjust(results_df$p_value, method = adjust_p_method)

  # Sort if requested
  if (sort_p) {
    results_df <- results_df[order(results_df$adj_p_value), ]
  }

  # Format p-values for display
  results_df$p_value_formatted <- ifelse(results_df$p_value < 0.001, "<.001",
                                         round(results_df$p_value, 3))
  results_df$adj_p_value_formatted <- ifelse(results_df$adj_p_value < 0.001, "<.001",
                                             round(results_df$adj_p_value, 3))

  # Add significance indicators
  results_df$significance <- case_when(
    results_df$adj_p_value < 0.001 ~ "***",
    results_df$adj_p_value < 0.01  ~ "**",
    results_df$adj_p_value < 0.05  ~ "*",
    TRUE ~ ""
  )

  # Create assumptions summary
  assumptions_summary <- data.frame(
    metabolite = rownames(results_df),
    test_used = results_df$test_used,
    normality_assumption = results_df$normality_p > 0.05,
    variance_assumption = results_df$variance_p > 0.05,
    stringsAsFactors = FALSE
  )

  return(list(
    results = results_df,
    assumptions = assumptions_summary
  ))
}

# Helper function: Generate summary statistics
.generate_summary <- function(results, groups) {
  list(
    total_metabolites = nrow(results),
    significant_metabolites = sum(results$adj_p_value < 0.05, na.rm = TRUE),
    highly_significant = sum(results$adj_p_value < 0.01, na.rm = TRUE),
    very_highly_significant = sum(results$adj_p_value < 0.001, na.rm = TRUE),
    test_distribution = table(results$test_used),
    group_sizes = table(groups),
    median_effect_size = median(results$effect_size, na.rm = TRUE)
  )
}

# Helper function: Generate plots
.generate_plots <- function(plot_metabolites, df, groups, results, adjust_p_method, verbose) {
  if (!requireNamespace("ggstatsplot", quietly = TRUE)) {
    warning("ggstatsplot package not available. Skipping plot generation.")
    return(list())
  }

  plots_list <- list()

  for (metabolite in plot_metabolites) {
    if (!metabolite %in% rownames(results)) {
      if (verbose) {
        cat("Warning: Metabolite '", metabolite, "' not found in results. Skipping plot.\n")
      }
      next
    }

    tryCatch({
      # Get test type for the metabolite
      test_used <- results[metabolite, "test_used"]
      test_type <- if (grepl("Mann-Whitney|Kruskal-Wallis|Wilcoxon", test_used)) {
        "nonparametric"
      } else {
        "parametric"
      }

      # Create plot data
      plot_data <- data.frame(
        group = groups,
        value = df[[metabolite]],
        stringsAsFactors = FALSE
      )

      # Generate plot
      plot_obj <- ggstatsplot::ggbetweenstats(
        data = plot_data,
        x = group,
        y = value,
        type = test_type,
        p.adjust.method = adjust_p_method,
        title = paste("Comparative Analysis:", metabolite),
        package = "ggplot2",
        palette = "Set2"
      )

      plots_list[[metabolite]] <- plot_obj

      if (verbose) {
        cat("Generated plot for:", metabolite, "\n")
      }

    }, error = function(e) {
      if (verbose) {
        cat("Error generating plot for", metabolite, ":", e$message, "\n")
      }
    })
  }

  return(plots_list)
}

#' Perform Comparative Statistical Analysis on a Preprocessed Data
#'
#' @description
#' This function performs comparative analysis based on the characteristics of the data. This means that this function
#' checks for the assumptions first before proceeding with the statistical analysis. The tests include t-tests and
#' Analysis of Variance (ANOVA), and its types such as paired samples t-test, and their nonparametric counterparts such as
#' Mann-Whitney U test, Kruskal-Wallis test, etc.
#'
#' @param data  List. This list must be a result from the `perform_PreprocessingPeakData` function.
#' @param adjust_p_method String. The p-value correction method. Defaults to "BH". Read more at `?p.adjust`.
#'   \itemize{
#'     \item "holm": Holm (1979). Less conservative than "bonferroni".
#'     \item "hochberg": Hochberg (1988). Less conservative than "bonferroni".
#'     \item "hommel": Hommel (1988). Less conservative than "bonferroni".
#'     \item "bonferroni": Bonferroni correction. The p-values are multiplied by the number of comparisons.
#'     \item "BH": Benjamini & Hochberg (1995). Less conservative than "bonferroni".
#'     \item "BY": Benjamini & Yekutieli (2001). Less conservative than "bonferroni".
#'     \item "fdr": Same as "BH".
#'     \item "none":
#'     }
#'     Note: The first four methods are designed to give strong control of the family-wise error rate. There seems no reason to use the unmodified Bonferroni correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions. Hochberg's and Hommel's methods are valid when the hypothesis tests are independent or when they are non-negatively associated ( Sarkar, 1998; Sarkar and Chang, 1997). Hommel's method is more powerful than Hochberg's, but the difference is usually small and the Hochberg p-values are faster to compute. The "BH" (aka "fdr") and "BY" methods of Benjamini, Hochberg, and Yekutieli control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses. The false discovery rate is a less stringent condition than the family-wise error rate, so these methods are more powerful than the others.
#' @param sort_p Boolean. If `TRUE` (default), sorts the adjusted p-values in ascending order.
#' @param paired Boolean. If `TRUE`, this instructs the function that the data is paired, not independent. All p-values are adjusted using Bonferroni and Hochberg method.
#' @param plot_iden_met Boolean. Setting this to `TRUE` will output a plot of the identified metabolites using `ggstatsplot` packagge, which means that statistical tests and p-values are displayed in the plot.
#'
#' @returns A list of data frames and plots, if requested.
#' @export
#'
#' @examples
#' \dontrun{
#' perform_ComparativeAnalysis(data = results_from_perform_PreprocessingPeakData_function)
#'}
#'
perform_ComparativeAnalysis <- function(
    data,
    adjust_p_method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")[5],
    sort_p          = TRUE,
    paired          = FALSE,
    plot_iden_met   = NULL
) {

  comparativeAnalysisResults                <- base::list()
  comparativeAnalysisResults$FunctionOrigin <- "perform_ComparativeAnalysis"

  qc_indices     <- data$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices <- !qc_indices
  df             <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ]
  groups         <- data$Metadata$Group[non_qc_indices]

  num_groups <- base::length(base::unique(groups))

  # Two-group comparison
  if (num_groups == 2) {
    unique_groups     <- base::unique(groups)
    idx_group1        <- base::which(groups == unique_groups[1])
    idx_group2        <- base::which(groups == unique_groups[2])

    results           <- base::apply(df, 2, function(x) {
      residuals       <- stats::residuals(stats::lm(x ~ groups))
      # normality       <- shapiro.test(x)$p.value > 0.05
      normality       <- stats::shapiro.test(residuals)$p.value > 0.05
      variance        <- stats::var.test(x[idx_group1], x[idx_group2])$p.value > 0.05

      if (paired) {
        if (normality) {
          test_used   <- "Paired t-test"
          p_value     <- stats::t.test(x[idx_group1], x[idx_group2],      paired = TRUE)$p.value
        } else {
          test_used   <- "Wilcoxon Signed-Rank test"
          p_value     <- stats::wilcox.test(x[idx_group1], x[idx_group2], paired = TRUE)$p.value
        }
      } else {
        if (normality) {
          if (variance) {
            test_used <- "Independent Samples t-test"
            p_value   <- stats::t.test(x[idx_group1], x[idx_group2])$p.value
          } else {
            test_used <- "Welch's t-test"
            p_value   <- stats::t.test(x[idx_group1], x[idx_group2], var.equal = FALSE)$p.value
          }
        } else {
          test_used   <- "Mann-Whitney U test"
          p_value     <- stats::wilcox.test(x[idx_group1], x[idx_group2])$p.value
        }
      }
      return(base::c(test_used, p_value))
    })

  } else { # Multi-group comparison

    results <- base::apply(df, 2, function(x) {
      normality     <- base::all(base::by(x, groups, function(g) stats::shapiro.test(g)$p.value > 0.05))
      variance      <- stats::bartlett.test(x ~ groups)$p.value > 0.05

      if (paired) {
        test_used   <- "Repeated Measures ANOVA (MANOVA)"
        p_value     <- base::summary(stats::aov(x ~ groups + Error(groups)))$p.value
      } else {
        if (normality && variance) {
          test_used <- "ANOVA"
          p_value   <- base::summary(stats::aov(x ~ groups))[[1]][["Pr(>F)"]][1]
        } else {
          test_used <- "Kruskal-Wallis"
          p_value   <- stats::kruskal.test(x ~ groups)$p.value
        }
      }
      return(c(test_used, p_value))
    })
  }

  results                <- base::as.data.frame(base::t(results), stringsAsFactors = FALSE)
  colnames(results)      <- base::c("Test used", "p-value")
  results$`p-value`      <- base::as.numeric(results$`p-value`)
  results$`adj. p-value` <- stats::p.adjust(results$`p-value`, method = adjust_p_method)

  if(sort_p) {
    results <- results %>% dplyr::arrange(`adj. p-value`) # Sort p-values
  }

  results$`p-value2`      <- base::ifelse(results$`p-value`      < 0.001, "<.001", base::round(results$`p-value`,      3))
  results$`adj. p-value2` <- base::ifelse(results$`adj. p-value` < 0.001, "<.001", base::round(results$`adj. p-value`, 3))

  results$Significance <- dplyr::case_when(
    base::as.numeric(results$`adj. p-value`) < 0.001 ~ "***",
    base::as.numeric(results$`adj. p-value`) < 0.01  ~ "**",
    base::as.numeric(results$`adj. p-value`) < 0.05  ~ "*",
    TRUE                                             ~ ""
  )

  comparativeAnalysisResults$results <- results

  # Plotting identified metabolites using ggstatsplot
  if (!base::is.null(plot_iden_met)) {
    for (metabolite in plot_iden_met) {
      if (metabolite %in% base::rownames(results)) {  # Check if metabolite is in results

        # Get the test used for the metabolite
        test_used <- results[metabolite, "Test used"]

        # Set the test type dynamically based on the test used
        test_type <- base::ifelse(test_used %in% base::c("Mann-Whitney U test", "Kruskal-Wallis"),
                                  "nonparametric", "parametric")

        # Create the plot with dynamic test type
        plot_the_met <- NULL # Control, because if 'plot_the_met' will be returned as 'not found' otherwise
        plot_the_met <- ggstatsplot::ggbetweenstats(
          data = base::data.frame(x = groups,
                                  y = df[[metabolite]]),
          x    = "x",
          y    = "y",
          type = test_type,  # Use the dynamically determined test type
          # pairwise = TRUE,
          p.adjust.method = adjust_p_method,  # Apply p-value correction
          title = base::paste0("Comparative Analysis of ", metabolite)
        )

        if(!base::is.null(plot_the_met)) {

          comparativeAnalysisResults$plots[[metabolite]] <- plot_the_met

          print(plot_the_met)

        } else {
          next
        }
      } else {
        print(base::paste0("Metabolite '", metabolite, "' is not found in the preprocessed data (it has been removed in one of the data preprocessing steps). Boxplot is not generated."))

      }
    }
  }

  return(comparativeAnalysisResults)
}
