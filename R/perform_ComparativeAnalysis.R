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
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{perform_PreprocessingPeakData}}, \code{\link[stats]{t.test}},
#'   \code{\link[stats]{aov}}, \code{\link[ggstatsplot]{ggbetweenstats}}
#'
#' @importFrom stats aov bartlett.test kruskal.test lm p.adjust sd shapiro.test t.test var var.test wilcox.test
#'
#' @references Bartlett, M. S. (1937). Properties of sufficiency and statistical tests. Proceedings of the Royal Society of London Series A 160, 268–282. doi:10.1098/rspa.1937.0109. (for bartlett.test)
#' @references David F. Bauer (1972). Constructing confidence sets using rank statistics. Journal of the American Statistical Association 67, 687–690. doi:10.1080/01621459.1972.10481279. (for wilcox.test)
#' @references Myles Hollander and Douglas A. Wolfe (1973). Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 27–33 (one-sample), 68–75 (two-sample). Or second edition (1999). (for wilcox.test)
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). The New S Language. Wadsworth & Brooks/Cole. (for var)
#' @references Kendall, M. G. (1938). A new measure of rank correlation, Biometrika, 30, 81–93. doi:10.1093/biomet/30.1-2.81. (for var)
#' @references Kendall, M. G. (1945). The treatment of ties in rank problems. Biometrika, 33 239–251. doi:10.1093/biomet/33.3.239 (for var)
#' @references Chambers, J. M. (1992) Linear models. Chapter 4 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole. (for lm)
#' @references Wilkinson, G. N. and Rogers, C. E. (1973). Symbolic descriptions of factorial models for analysis of variance. Applied Statistics, 22, 392–399. doi:10.2307/2346786. (for lm)
#' @references Chambers, J. M., Freeny, A and Heiberger, R. M. (1992) Analysis of variance; designed experiments. Chapter 5 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole. (for aov)
#' @references Myles Hollander and Douglas A. Wolfe (1973), Nonparametric Statistical Methods. New York: John Wiley & Sons. Pages 115–120. (for kruskal.test)
#' @references Patrick Royston (1982). An extension of Shapiro and Wilk's W test for normality to large samples. Applied Statistics, 31, 115–124. doi:10.2307/2347973. (for shapiro.test)
#' @references Patrick Royston (1982). Algorithm AS 181: The W test for Normality. Applied Statistics, 31, 176–180. doi:10.2307/2347986. (for shapiro.test)
#' @references Patrick Royston (1995). Remark AS R94: A remark on Algorithm AS 181: The W test for normality. Applied Statistics, 44, 547–551. doi:10.2307/2986146. (for shapiro.test)
#' @references Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x. (for p.adjust)
#' @references Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165–1188. doi:10.1214/aos/1013699998. (for p.adjust)
#' @references Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70. https://www.jstor.org/stable/4615733. (for p.adjust)
#' @references Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383–386. doi:10.2307/2336190. (for p.adjust)
#' @references Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800–803. doi:10.2307/2336325. (for p.adjust)
#' @references Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561–584. doi:10.1146/annurev.ps.46.020195.003021. (An excellent review of the area.) (for p.adjust)
#' @references Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics, 26, 494–504. doi:10.1214/aos/1028144846. (for p.adjust)
#' @references Sarkar, S., and Chang, C. K. (1997). The Simes method for multiple hypothesis testing with positively dependent test statistics. Journal of the American Statistical Association, 92, 1601–1608. doi:10.2307/2965431. (for p.adjust)
#' @references Wright, S. P. (1992). Adjusted P-values for simultaneous inference. Biometrics, 48, 1005–1013. doi:10.2307/2532694. (Explains the adjusted P-value approach.) (for p.adjust)
#'
#' @export
#'
perform_ComparativeAnalysis <- function(
    data,
    adjust_p_method = "BH",
    sort_p = TRUE,
    paired = FALSE,
    plot_metabolites = NULL,
    alpha = 0.05,
    min_group_size = 3,
    verbose = FALSE
) {

  # Validate inputs - Fixed to pass all parameters correctly
  .validate_inputs_comparativeanalysis(
    data = data,
    adjust_p_method = adjust_p_method,
    alpha = alpha,
    min_group_size = min_group_size,
    plot_metabolites = plot_metabolites,
    paired = paired,
    sort_p = sort_p,
    verbose = verbose
  )

  # Initialize results structure
  results_list <- .initialize_results()

  if (verbose) cat("Starting comparative analysis...\n")

  # Extract and validate data
  analysis_data <- .extract_analysis_data_comparativeanalysis(data, min_group_size, verbose)
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
  if (!is.null(plot_metabolites) && length(plot_metabolites) > 0) {
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

# Helper function: Validate inputs - Fixed to include all parameters
.validate_inputs_comparativeanalysis <- function(data, adjust_p_method, alpha, min_group_size,
                                                 plot_metabolites = NULL, paired = FALSE,
                                                 sort_p = TRUE, verbose = FALSE) {
  # Check data structure
  if (!is.list(data)) {
    stop("Input 'data' must be a list from perform_PreprocessingPeakData function")
  }

  required_elements <- c("data_scaledPCA_rsdFiltered_varFiltered", "Metadata")
  missing_elements <- setdiff(required_elements, names(data))
  if (length(missing_elements) > 0) {
    stop("Missing required data elements: ", paste(missing_elements, collapse = ", "))
  }

  # Validate data types
  if (!is.data.frame(data$Metadata) && !is.matrix(data$Metadata)) {
    stop("Metadata must be a data frame or matrix")
  }

  if (!"Group" %in% colnames(data$Metadata)) {
    stop("Metadata must contain a 'Group' column")
  }

  if (!is.matrix(data$data_scaledPCA_rsdFiltered_varFiltered) &&
      !is.data.frame(data$data_scaledPCA_rsdFiltered_varFiltered)) {
    stop("data_scaledPCA_rsdFiltered_varFiltered must be a matrix or data frame")
  }

  # Check p-value adjustment method
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!adjust_p_method %in% valid_methods) {
    stop("Invalid adjust_p_method. Must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Check numeric parameters
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value between 0 and 1")
  }

  if (!is.numeric(min_group_size) || length(min_group_size) != 1 || min_group_size < 2) {
    stop("min_group_size must be a single numeric value >= 2")
  }

  # Check logical parameters
  if (!is.logical(paired) || length(paired) != 1) {
    stop("paired must be a single logical value")
  }

  if (!is.logical(sort_p) || length(sort_p) != 1) {
    stop("sort_p must be a single logical value")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Check plot_metabolites parameter
  if (!is.null(plot_metabolites)) {
    if (!is.character(plot_metabolites)) {
      stop("plot_metabolites must be a character vector or NULL")
    }

    if (length(plot_metabolites) == 0) {
      warning("plot_metabolites is an empty vector. No plots will be generated.")
    }

    # Check if required packages are available for plotting
    if (!requireNamespace("ggstatsplot", quietly = TRUE)) {
      warning("ggstatsplot package not available. Plots will not be generated.")
    }
  }

  # Check data dimensions compatibility
  if (nrow(data$data_scaledPCA_rsdFiltered_varFiltered) != nrow(data$Metadata)) {
    stop("Number of rows in data and metadata must match")
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
.extract_analysis_data_comparativeanalysis <- function(data, min_group_size, verbose) {
  # Identify QC samples more efficiently
  qc_patterns <- c("SQC", "EQC", "QC", "BLANK", "blank", "Blank")
  groups_col <- data$Metadata$Group

  # Convert to character for pattern matching if factor
  if (is.factor(groups_col)) {
    groups_char <- as.character(groups_col)
  } else {
    groups_char <- groups_col
  }

  qc_indices <- groups_char %in% qc_patterns

  if (verbose && sum(qc_indices) > 0) {
    cat("Removing", sum(qc_indices), "QC samples from analysis\n")
  }

  # Extract non-QC data
  non_qc_indices <- !qc_indices

  # Ensure we have data left after QC removal
  if (sum(non_qc_indices) == 0) {
    stop("No samples remaining after removing QC samples")
  }

  df <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, , drop = FALSE]
  groups <- factor(groups_char[non_qc_indices])

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

  # Check for at least 2 groups
  if (length(group_counts) < 2) {
    stop("At least 2 groups are required for comparative analysis")
  }

  # Check for missing values and handle them
  na_count <- sum(is.na(df))
  if (na_count > 0) {
    warning("Missing values detected in data (", na_count, " values). Results may be unreliable.")
    if (verbose) {
      cat("Missing values found in", sum(apply(is.na(df), 2, any)), "metabolites\n")
    }
  }

  return(list(df = as.data.frame(df), groups = groups))
}

# Helper function: Validate paired analysis requirements
.validate_paired_analysis <- function(groups, num_groups) {
  group_counts <- table(groups)

  if (num_groups > 2) {
    if (!all(group_counts == group_counts[1])) {
      stop("Paired analysis requires equal sample sizes across all groups. Current sizes: ",
           paste(names(group_counts), "=", group_counts, collapse = ", "))
    }
  } else if (num_groups == 2) {
    if (group_counts[1] != group_counts[2]) {
      stop("Paired analysis requires equal sample sizes in both groups. Current sizes: ",
           paste(names(group_counts), "=", group_counts, collapse = ", "))
    }
  }
}

# Helper function: Perform two-group analysis with better error handling
.perform_two_group_analysis <- function(df, groups, paired, alpha, verbose) {
  unique_groups <- levels(groups)
  idx_group1 <- which(groups == unique_groups[1])
  idx_group2 <- which(groups == unique_groups[2])

  # Pre-allocate results for better performance
  n_metabolites <- ncol(df)
  metabolite_names <- colnames(df)

  # Initialize result vectors
  test_used <- character(n_metabolites)
  p_values <- numeric(n_metabolites)
  normality_p <- numeric(n_metabolites)
  variance_p <- numeric(n_metabolites)
  effect_sizes <- numeric(n_metabolites)

  # Initialize vectors for group means/medians
  group1_stat <- numeric(n_metabolites)
  group2_stat <- numeric(n_metabolites)
  group1_stat_name <- character(n_metabolites)
  group2_stat_name <- character(n_metabolites)

  # Process metabolites with progress indication
  if (verbose && n_metabolites > 100) {
    cat("Processing", n_metabolites, "metabolites...\n")
    progress_points <- seq(1, n_metabolites, length.out = min(10, n_metabolites))
  }

  for (i in seq_len(n_metabolites)) {
    if (verbose && n_metabolites > 100 && i %in% progress_points) {
      cat("Progress:", round(100 * i / n_metabolites), "%\n")
    }

    tryCatch({
      x <- df[, i]
      x1 <- x[idx_group1]
      x2 <- x[idx_group2]

      # Skip if too many missing values
      if (sum(is.na(c(x1, x2))) > 0.5 * length(c(x1, x2))) {
        test_used[i] <- "Insufficient Data"
        p_values[i] <- NA
        normality_p[i] <- NA
        variance_p[i] <- NA
        effect_sizes[i] <- NA
        group1_stat[i] <- NA
        group2_stat[i] <- NA
        group1_stat_name[i] <- NA
        group2_stat_name[i] <- NA
        next
      }

      # Remove missing values
      if (paired) {
        complete_cases <- complete.cases(x1, x2)
        x1 <- x1[complete_cases]
        x2 <- x2[complete_cases]
        if (length(x1) < 3) {
          test_used[i] <- "Insufficient Data"
          p_values[i] <- NA
          normality_p[i] <- NA
          variance_p[i] <- NA
          effect_sizes[i] <- NA
          group1_stat[i] <- NA
          group2_stat[i] <- NA
          group1_stat_name[i] <- NA
          group2_stat_name[i] <- NA
          next
        }
      } else {
        x1 <- x1[!is.na(x1)]
        x2 <- x2[!is.na(x2)]
        if (length(x1) < 3 || length(x2) < 3) {
          test_used[i] <- "Insufficient Data"
          p_values[i] <- NA
          normality_p[i] <- NA
          variance_p[i] <- NA
          effect_sizes[i] <- NA
          group1_stat[i] <- NA
          group2_stat[i] <- NA
          group1_stat_name[i] <- NA
          group2_stat_name[i] <- NA
          next
        }
      }

      # Test assumptions
      normality_test <- .test_normality_safe(c(x1, x2), c(rep(unique_groups[1], length(x1)),
                                                          rep(unique_groups[2], length(x2))))
      variance_test <- if (!paired) .test_variance_equality_safe(x1, x2) else list(p_value = 1)

      # Select and perform appropriate test
      test_result <- .select_and_perform_test_two_group_safe(
        x1, x2, normality_test$is_normal, variance_test$p_value > alpha, paired
      )

      # Calculate effect size
      effect_size <- .calculate_effect_size_two_group_safe(x1, x2, paired)

      # Determine which statistic to use (mean or median)
      is_parametric <- grepl("t-test", test_result$test_used, ignore.case = TRUE)
      if (is_parametric) {
        group1_stat[i] <- mean(x1, na.rm = TRUE)
        group2_stat[i] <- mean(x2, na.rm = TRUE)
        group1_stat_name[i] <- paste0("mean_", unique_groups[1])
        group2_stat_name[i] <- paste0("mean_", unique_groups[2])
      } else {
        group1_stat[i] <- median(x1, na.rm = TRUE)
        group2_stat[i] <- median(x2, na.rm = TRUE)
        group1_stat_name[i] <- paste0("median_", unique_groups[1])
        group2_stat_name[i] <- paste0("median_", unique_groups[2])
      }

      # Store results
      test_used[i] <- test_result$test_used
      p_values[i] <- test_result$p_value
      normality_p[i] <- normality_test$p_value
      variance_p[i] <- variance_test$p.value
      effect_sizes[i] <- effect_size

    }, error = function(e) {
      if (verbose) {
        cat("Error processing metabolite", metabolite_names[i], ":", e$message, "\n")
      }
      test_used[i] <- "Test Failed"
      p_values[i] <- NA
      normality_p[i] <- NA
      variance_p[i] <- NA
      effect_sizes[i] <- NA
      group1_stat[i] <- NA
      group2_stat[i] <- NA
      group1_stat_name[i] <- NA
      group2_stat_name[i] <- NA
    })
  }

  # Create results matrix with group statistics
  results_matrix <- data.frame(
    test_used = test_used,
    p_value = p_values,
    normality_p = normality_p,
    variance_p = variance_p,
    effect_size = effect_sizes,
    row.names = metabolite_names,
    stringsAsFactors = FALSE
  )

  # Add group statistics, handling different stat names for each row
  unique_stat_names <- unique(c(group1_stat_name, group2_stat_name))
  for (name in unique_stat_names) {
    results_matrix[[name]] <- NA
  }

  for (i in seq_len(n_metabolites)) {
    if (!is.na(group1_stat_name[i])) {
      results_matrix[i, group1_stat_name[i]] <- group1_stat[i]
    }
    if (!is.na(group2_stat_name[i])) {
      results_matrix[i, group2_stat_name[i]] <- group2_stat[i]
    }
  }


  return(list(
    results = results_matrix,
    group_info = list(groups = unique_groups, indices = list(idx_group1, idx_group2))
  ))
}

# Helper function: Perform multi-group analysis
.perform_multi_group_analysis <- function(df, groups, paired, alpha, verbose) {
  n_metabolites <- ncol(df)
  metabolite_names <- colnames(df)
  unique_groups <- levels(groups)
  num_groups <- length(unique_groups)

  # Initialize result vectors
  test_used <- character(n_metabolites)
  p_values <- numeric(n_metabolites)
  normality_p <- numeric(n_metabolites)
  variance_p <- numeric(n_metabolites)
  effect_sizes <- numeric(n_metabolites)
  group_stat_type <- character(n_metabolites)

  # Initialize a list to hold the group stats for each metabolite
  group_stats_list <- vector("list", n_metabolites)

  # Process metabolites with progress indication
  if (verbose && n_metabolites > 100) {
    cat("Processing", n_metabolites, "metabolites...\n")
  }

  for (i in seq_len(n_metabolites)) {
    tryCatch({
      x <- df[, i]

      # Handle missing values
      complete_cases <- complete.cases(x)
      x_complete <- x[complete_cases]
      groups_complete <- groups[complete_cases]

      # Check if we have enough data
      group_sizes <- table(groups_complete)
      if (any(group_sizes < 3) || length(group_sizes) < 2) {
        test_used[i] <- "Insufficient Data"
        p_values[i] <- NA
        normality_p[i] <- NA
        variance_p[i] <- NA
        effect_sizes[i] <- NA
        group_stat_type[i] <- NA
        group_stats_list[[i]] <- NA
        next
      }

      # Test assumptions
      normality_test <- .test_normality_multi_group_safe(x_complete, groups_complete)
      variance_test <- if (!paired) .test_variance_homogeneity_safe(x_complete, groups_complete) else list(p_value = 1)

      # Select and perform appropriate test
      test_result <- .select_and_perform_test_multi_group_safe(
        x_complete, groups_complete, normality_test$all_normal, variance_test$p_value > alpha, paired
      )

      # Calculate effect size
      effect_size <- .calculate_effect_size_multi_group_safe(x_complete, groups_complete)

      # Determine which statistic to use (mean or median)
      is_parametric <- grepl("ANOVA", test_result$test_used, ignore.case = TRUE)
      if (is_parametric) {
        group_stats <- tapply(x_complete, groups_complete, mean, na.rm = TRUE)
        group_stat_type[i] <- "mean"
      } else {
        group_stats <- tapply(x_complete, groups_complete, median, na.rm = TRUE)
        group_stat_type[i] <- "median"
      }
      # Store the named vector of stats
      group_stats_list[[i]] <- group_stats

      # Store results
      test_used[i] <- test_result$test_used
      p_values[i] <- test_result$p_value
      normality_p[i] <- normality_test$p_value
      variance_p[i] <- variance_test$p_value
      effect_sizes[i] <- effect_size

    }, error = function(e) {
      if (verbose) {
        cat("Error processing metabolite", metabolite_names[i], ":", e$message, "\n")
      }
      test_used[i] <- "Test Failed"
      p_values[i] <- NA
      normality_p[i] <- NA
      variance_p[i] <- NA
      effect_sizes[i] <- NA
      group_stat_type[i] <- NA
      group_stats_list[[i]] <- NA
    })
  }

  # Prepare the group stats for merging
  final_group_stats_df <- data.frame(row.names = metabolite_names)
  group_names_to_add <- levels(groups)

  for (i in 1:length(metabolite_names)) {
    if (!is.na(group_stat_type[i])) {
      stat_type <- group_stat_type[i]
      current_stats <- group_stats_list[[i]]

      for (group_name in group_names_to_add) {
        col_name <- paste0(stat_type, "_", group_name)

        # Initialize column if it doesn't exist
        if (!col_name %in% names(final_group_stats_df)) {
          final_group_stats_df[[col_name]] <- NA
        }

        # Assign the value
        if (group_name %in% names(current_stats)) {
          final_group_stats_df[metabolite_names[i], col_name] <- current_stats[group_name]
        }
      }
    }
  }


  # Create the main results data frame
  results_matrix <- data.frame(
    test_used = test_used,
    p_value = p_values,
    normality_p = normality_p,
    variance_p = variance_p,
    effect_size = effect_sizes,
    row.names = metabolite_names,
    stringsAsFactors = FALSE
  )

  # Combine results with group statistics
  results_with_stats <- cbind(results_matrix, final_group_stats_df)


  return(list(
    results = results_with_stats,
    group_info = list(groups = levels(groups))
  ))
}

# Safe wrapper functions for statistical tests
.test_normality_safe <- function(x, groups) {
  tryCatch({
    if (length(x) < 3) {
      return(list(is_normal = FALSE, p_value = 0))
    }

    # Remove infinite values
    finite_indices <- is.finite(x)
    if (sum(finite_indices) < 3) {
      return(list(is_normal = FALSE, p_value = 0))
    }

    x_clean <- x[finite_indices]
    groups_clean <- groups[finite_indices]

    residuals <- residuals(lm(x_clean ~ groups_clean))

    if (length(residuals) < 3 || var(residuals) == 0) {
      return(list(is_normal = FALSE, p_value = 0))
    }

    shapiro_result <- shapiro.test(residuals)
    list(is_normal = shapiro_result$p.value > 0.05, p_value = shapiro_result$p.value)
  }, error = function(e) {
    list(is_normal = FALSE, p_value = 0)
  })
}

.test_normality_multi_group_safe <- function(x, groups) {
  tryCatch({
    group_normality <- tapply(x, groups, function(g) {
      if (length(g) < 3 || !any(is.finite(g)) || var(g, na.rm = TRUE) == 0) return(0)
      g_clean <- g[is.finite(g)]
      if (length(g_clean) < 3) return(0)
      shapiro.test(g_clean)$p.value
    })

    # Handle NULL results
    group_normality[sapply(group_normality, is.null)] <- 0
    group_normality <- as.numeric(group_normality)

    min_p <- min(group_normality, na.rm = TRUE)
    all_normal <- all(group_normality > 0.05, na.rm = TRUE) && !is.na(min_p)

    list(all_normal = all_normal, p_value = ifelse(is.finite(min_p), min_p, 0))
  }, error = function(e) {
    list(all_normal = FALSE, p_value = 0)
  })
}

.test_variance_equality_safe <- function(x1, x2) {
  tryCatch({
    # Remove non-finite values
    x1_clean <- x1[is.finite(x1)]
    x2_clean <- x2[is.finite(x2)]

    if (length(x1_clean) < 2 || length(x2_clean) < 2) {
      return(list(p_value = 0))
    }

    # Check for zero variance
    if (var(x1_clean) == 0 || var(x2_clean) == 0) {
      return(list(p_value = 0))
    }

    var_test <- var.test(x1_clean, x2_clean)
    list(p_value = var_test$p.value)
  }, error = function(e) {
    list(p_value = 0)
  })
}

.test_variance_homogeneity_safe <- function(x, groups) {
  tryCatch({
    # Remove non-finite values
    finite_indices <- is.finite(x)
    x_clean <- x[finite_indices]
    groups_clean <- groups[finite_indices]

    if (length(x_clean) < 5) {
      return(list(p_value = 0))
    }

    # Check group sizes
    group_sizes <- table(groups_clean)
    if (any(group_sizes < 2)) {
      return(list(p_value = 0))
    }

    bartlett_test <- bartlett.test(x_clean ~ groups_clean)
    list(p_value = bartlett_test$p.value)
  }, error = function(e) {
    list(p_value = 0)
  })
}

# Safe test selection and execution functions
.select_and_perform_test_two_group_safe <- function(x1, x2, is_normal, equal_variance, paired) {
  tryCatch({
    if (paired) {
      if (is_normal) {
        result <- t.test(x1, x2, paired = TRUE)
        list(test_used = "Paired t-test", p_value = result$p.value)
      } else {
        result <- wilcox.test(x1, x2, paired = TRUE, exact = FALSE)
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
        result <- wilcox.test(x1, x2, exact = FALSE)
        list(test_used = "Mann-Whitney U test", p_value = result$p.value)
      }
    }
  }, error = function(e) {
    list(test_used = "Test Failed", p_value = NA)
  })
}

.select_and_perform_test_multi_group_safe <- function(x, groups, all_normal, equal_variance, paired) {
  tryCatch({
    if (paired) {
      # For paired multi-group, create subject IDs
      group_levels <- levels(groups)
      n_per_group <- table(groups)[1]  # Should be equal for paired analysis
      subject_id <- rep(seq_len(n_per_group), length(group_levels))

      if (length(subject_id) != length(x)) {
        # Fallback to regular ANOVA if subject matching fails
        aov_result <- aov(x ~ groups)
        p_value <- summary(aov_result)[[1]][["Pr(>F)"]][1]
        list(test_used = "One-way ANOVA", p_value = p_value)
      } else {
        aov_result <- aov(x ~ groups + Error(factor(subject_id)))
        summary_result <- summary(aov_result)

        # Extract p-value safely
        if ("Error: Within" %in% names(summary_result)) {
          p_value <- summary_result$`Error: Within`[[1]]$`Pr(>F)`[1]
        } else {
          p_value <- summary_result[[1]][[1]][["Pr(>F)"]][1]
        }

        list(test_used = "Repeated Measures ANOVA", p_value = p_value)
      }
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
    list(test_used = "Test Failed", p_value = NA)
  })
}

# Safe effect size calculation functions
.calculate_effect_size_two_group_safe <- function(x1, x2, paired) {
  tryCatch({
    # Remove non-finite values
    if (paired) {
      complete_cases <- complete.cases(x1, x2)
      x1_clean <- x1[complete_cases]
      x2_clean <- x2[complete_cases]

      if (length(x1_clean) < 3) return(0)

      diff <- x1_clean - x2_clean
      if (sd(diff) == 0) return(0)

      cohen_d <- mean(diff) / sd(diff)
    } else {
      x1_clean <- x1[is.finite(x1)]
      x2_clean <- x2[is.finite(x2)]

      if (length(x1_clean) < 3 || length(x2_clean) < 3) return(0)

      # Check for zero variance
      if (var(x1_clean) == 0 && var(x2_clean) == 0) return(0)

      pooled_sd <- sqrt(((length(x1_clean) - 1) * var(x1_clean) +
                           (length(x2_clean) - 1) * var(x2_clean)) /
                          (length(x1_clean) + length(x2_clean) - 2))

      if (pooled_sd == 0) return(0)

      cohen_d <- (mean(x1_clean) - mean(x2_clean)) / pooled_sd
    }

    return(abs(cohen_d))
  }, error = function(e) {
    return(0)
  })
}

.calculate_effect_size_multi_group_safe <- function(x, groups) {
  tryCatch({
    # Remove non-finite values
    finite_indices <- is.finite(x)
    x_clean <- x[finite_indices]
    groups_clean <- groups[finite_indices]

    if (length(x_clean) < 5) return(0)

    # Check for sufficient group sizes
    group_sizes <- table(groups_clean)
    if (any(group_sizes < 2)) return(0)

    # Eta-squared calculation
    aov_result <- aov(x_clean ~ groups_clean)
    aov_summary <- summary(aov_result)[[1]]

    ss_between <- aov_summary[["Sum Sq"]][1]
    ss_total <- sum(aov_summary[["Sum Sq"]])

    if (ss_total == 0) return(0)

    eta_squared <- ss_between / ss_total
    return(eta_squared)
  }, error = function(e) {
    return(0)
  })
}

# Helper function: Process results with better error handling
.process_results <- function(stat_results, adjust_p_method, sort_p, verbose) {
  tryCatch({
    results_df <- stat_results$results
    group_info <- stat_results$group_info

    # Get group names
    group_names <- group_info$groups

    # Identify stat columns
    stat_cols_temp <- grep("^(mean|median)_", colnames(results_df), value = TRUE)
    stat_cols_stat <- grep("^stat_", colnames(results_df), value = TRUE)

    # Build the new results data frame from scratch
    if (length(group_names) == 2) {
      final_results_df <- data.frame(
        test_used = results_df$test_used,
        stringsAsFactors = FALSE
      )

      # Add stat columns in the desired order
      for (col_name in stat_cols_temp) {
        final_results_df[[col_name]] <- results_df[[col_name]]
      }

      # Add other columns
      other_cols <- setdiff(colnames(results_df), c("test_used", stat_cols_temp))
      final_results_df <- cbind(final_results_df, results_df[, other_cols, drop = FALSE])
      rownames(final_results_df) <- rownames(results_df)
      results_df <- final_results_df

    } else { # Multi-group case

      # Initialize final data frame
      final_results_df <- data.frame(
        test_used = results_df$test_used,
        row.names = rownames(results_df),
        stringsAsFactors = FALSE
      )

      # Dynamically build the group stat columns and populate them
      for (i in seq_along(results_df$test_used)) {
        test_type <- results_df$test_used[i]

        # Determine stat type based on test
        stat_type <- if (grepl("ANOVA", test_type)) "mean" else if (grepl("Kruskal-Wallis", test_type)) "median" else NA

        if (!is.na(stat_type)) {
          for (j in seq_along(group_names)) {
            col_name <- paste0(stat_type, "_", group_names[j])

            # Check if the temporary stat column exists in the original df
            temp_stat_col_name <- paste0("stat_", group_names[j])
            if (temp_stat_col_name %in% colnames(results_df)) {
              stat_value <- results_df[i, temp_stat_col_name]

              # Add the new column if it doesn't exist
              if (!(col_name %in% names(final_results_df))) {
                final_results_df[[col_name]] <- NA
              }

              final_results_df[i, col_name] <- stat_value
            }
          }
        }
      }

      # Add the remaining original columns
      other_cols <- setdiff(colnames(results_df), c("test_used", stat_cols_stat))
      final_results_df <- cbind(final_results_df, results_df[, other_cols, drop = FALSE])
      results_df <- final_results_df
    }

    # Handle missing p-values
    valid_p_indices <- !is.na(results_df$p_value) & is.finite(results_df$p_value)

    if (sum(valid_p_indices) == 0) {
      stop("No valid p-values found in results")
    }

    if (verbose && sum(!valid_p_indices) > 0) {
      cat("Warning:", sum(!valid_p_indices), "metabolites had invalid p-values\n")
    }

    # Apply p-value adjustment only to valid p-values
    results_df$adj_p_value <- NA
    if (sum(valid_p_indices) > 0) {
      results_df$adj_p_value[valid_p_indices] <- p.adjust(
        results_df$p_value[valid_p_indices],
        method = adjust_p_method
      )
    }

    # Sort if requested - handle NA values properly
    if (sort_p) {
      # Sort by adjusted p-values, putting NAs at the end
      results_df <- results_df[order(results_df$adj_p_value, na.last = TRUE), ]
    }

    # Format p-values for display - handle NA values
    results_df$p_value_formatted <- ifelse(
      is.na(results_df$p_value),
      "NA",
      ifelse(results_df$p_value < 0.001, "<.001", round(results_df$p_value, 3))
    )

    results_df$adj_p_value_formatted <- ifelse(
      is.na(results_df$adj_p_value),
      "NA",
      ifelse(results_df$adj_p_value < 0.001, "<.001", round(results_df$adj_p_value, 3))
    )

    # Add significance indicators - handle NA values
    results_df$significance <- ifelse(
      is.na(results_df$adj_p_value), "",
      ifelse(results_df$adj_p_value < 0.001, "***",
             ifelse(results_df$adj_p_value < 0.01, "**",
                    ifelse(results_df$adj_p_value < 0.05, "*", "")))
    )

    # Create assumptions summary
    assumptions_summary <- data.frame(
      metabolite = rownames(results_df),
      test_used = results_df$test_used,
      normality_assumption = ifelse(is.na(results_df$normality_p), FALSE,
                                    results_df$normality_p > 0.05),
      variance_assumption = ifelse(is.na(results_df$variance_p), FALSE,
                                   results_df$variance_p > 0.05),
      stringsAsFactors = FALSE
    )

    return(list(
      results = results_df,
      assumptions = assumptions_summary
    ))

  }, error = function(e) {
    stop("Error processing results: ", e$message)
  })
}

# Helper function: Generate summary statistics
.generate_summary <- function(results, groups) {
  tryCatch({
    # Count valid results
    valid_results <- !is.na(results$adj_p_value)

    summary_stats <- list(
      total_metabolites = nrow(results),
      valid_results = sum(valid_results),
      failed_tests = sum(results$test_used == "Test Failed", na.rm = TRUE),
      insufficient_data = sum(results$test_used == "Insufficient Data", na.rm = TRUE),
      significant_metabolites = sum(results$adj_p_value < 0.05, na.rm = TRUE),
      highly_significant = sum(results$adj_p_value < 0.01, na.rm = TRUE),
      very_highly_significant = sum(results$adj_p_value < 0.001, na.rm = TRUE),
      test_distribution = table(results$test_used),
      group_sizes = table(groups)
    )

    # Add effect size summary if available
    valid_effects <- !is.na(results$effect_size) & is.finite(results$effect_size)
    if (sum(valid_effects) > 0) {
      summary_stats$median_effect_size <- median(results$effect_size[valid_effects])
      summary_stats$mean_effect_size <- mean(results$effect_size[valid_effects])
    } else {
      summary_stats$median_effect_size <- NA
      summary_stats$mean_effect_size <- NA
    }

    return(summary_stats)

  }, error = function(e) {
    warning("Error generating summary statistics: ", e$message)
    return(list(
      total_metabolites = nrow(results),
      error = e$message
    ))
  })
}

# Helper function: Generate plots with better error handling
.generate_plots <- function(plot_metabolites, df, groups, results, adjust_p_method, verbose) {
  if (!requireNamespace("ggstatsplot", quietly = TRUE)) {
    warning("ggstatsplot package not available. Skipping plot generation.")
    return(list())
  }

  plots_list <- list()
  available_metabolites <- colnames(df)

  for (metabolite in plot_metabolites) {
    if (!metabolite %in% available_metabolites) {
      if (verbose) {
        cat("Warning: Metabolite '", metabolite, "' not found in data. Available metabolites: ",
            paste(head(available_metabolites), collapse = ", "), "...\n")
      }
      next
    }

    if (!metabolite %in% rownames(results)) {
      if (verbose) {
        cat("Warning: Metabolite '", metabolite, "' not found in results. Skipping plot.\n")
      }
      next
    }

    tryCatch({
      # Get test type for the metabolite
      test_used <- results[metabolite, "test_used"]

      # Skip if test failed
      if (is.na(test_used) || test_used %in% c("Test Failed", "Insufficient Data")) {
        if (verbose) {
          cat("Skipping plot for", metabolite, "- test failed or insufficient data\n")
        }
        next
      }

      test_type <- if (grepl("Mann-Whitney|Kruskal-Wallis|Wilcoxon", test_used)) {
        "nonparametric"
      } else {
        "parametric"
      }

      # Create plot data - handle missing values
      metabolite_data <- df[[metabolite]]
      complete_cases <- !is.na(metabolite_data) & is.finite(metabolite_data)

      if (sum(complete_cases) < 5) {
        if (verbose) {
          cat("Skipping plot for", metabolite, "- insufficient valid data points\n")
        }
        next
      }

      plot_data <- data.frame(
        group = groups[complete_cases],
        value = metabolite_data[complete_cases],
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

  if (length(plots_list) == 0 && length(plot_metabolites) > 0) {
    warning("No plots were successfully generated")
  }

  return(plots_list)
}
