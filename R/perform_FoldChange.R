#' Perform Fold Change Analysis on Preprocessed Peak Data
#'
#' @description
#' Performs pairwise fold change analysis on preprocessed metabolomics or proteomics
#' peak data. The function calculates fold changes and log2 fold changes between all
#' possible group pairs while excluding quality control (QC) samples. The analysis
#' uses group means for comparison and includes options for custom group ordering
#' and result sorting.
#'
#' @param data List. A preprocessed data object returned by `perform_PreprocessingPeakData()`
#'   or similar preprocessing functions. Must contain:
#'   - `$data_scaledPCA_rsdFiltered_varFiltered`: Numeric matrix/data.frame of peak data
#'   - `$Metadata`: Data.frame with at least a 'Group' column
#' @param arrangeLevels Character vector. Optional custom ordering for group levels.
#'   Must contain all unique group names present in the data (excluding QC samples).
#'   If `NULL` (default), groups are sorted alphabetically.
#' @param sortFC Logical. If `TRUE` (default), sorts results by fold change values
#'   in descending order within each comparison.
#' @param qc_patterns Character vector. Patterns to identify QC samples in the Group
#'   column. Default is `c("SQC", "EQC", "QC")`.
#' @param min_samples_per_group Integer. Minimum number of samples required per group
#'   for analysis. Default is 2.
#' @param epsilon Numeric. Small value added to prevent log2(0) issues when data
#'   contains zeros or negative values. Default is `1e-8`.
#'
#' @return A list containing:
#'   \item{FunctionOrigin}{Character. Function name for traceability}
#'   \item{data_shifted}{Data.frame. Adjusted data with minimum value shifted to 1}
#'   \item{group_summary}{Data.frame. Summary of groups and sample counts}
#'   \item{comparison_matrix}{Matrix. All pairwise comparisons performed}
#'   \item{data_combined_GROUP1 vs. GROUP2}{Data.frame. Fold change results for each comparison}
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input data structure and parameters
#' 2. Filters out QC samples based on specified patterns
#' 3. Checks for sufficient samples per group
#' 4. Handles zero/negative values by adding appropriate shift
#' 5. Performs all pairwise group comparisons
#' 6. Calculates fold changes as mean(group1)/mean(group2)
#' 7. Computes log2 fold changes
#' 8. Optionally sorts results by fold change magnitude
#'
#' @section Warning:
#' The function assumes that higher values indicate higher abundance/expression.
#' Ensure your data is appropriately normalized before analysis.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' fc_results <- perform_FoldChange(preprocessed_data)
#'
#' # Custom group ordering and no sorting
#' fc_results <- perform_FoldChange(
#'   data = preprocessed_data,
#'   arrangeLevels = c("Control", "Treatment1", "Treatment2"),
#'   sortFC = FALSE
#' )
#'
#' # Custom QC patterns
#' fc_results <- perform_FoldChange(
#'   data = preprocessed_data,
#'   qc_patterns = c("QC", "Blank", "POOL")
#' )
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{perform_PreprocessingPeakData}}
#'
#' @export
perform_FoldChange <- function(
    data,
    arrangeLevels = NULL,
    sortFC = TRUE,
    qc_patterns = c("SQC", "EQC", "QC"),
    min_samples_per_group = 2,
    epsilon = 1e-8
) {

  # Input validation
  .validate_foldchange_inputs(data, arrangeLevels, sortFC, qc_patterns,
                              min_samples_per_group, epsilon)

  # Initialize results list
  fc_results <- list(
    FunctionOrigin = "perform_FoldChange",
    analysis_timestamp = Sys.time()
  )

  # Extract and validate data components
  peak_data <- .extract_peak_data(data)
  metadata <- .extract_metadata(data)

  # Filter QC samples
  qc_filter_result <- .filter_qc_samples(metadata, qc_patterns)
  non_qc_indices <- qc_filter_result$non_qc_indices

  if (length(non_qc_indices) == 0) {
    stop("No non-QC samples found after filtering. Check your QC patterns.")
  }

  # Subset data and metadata
  df_filtered <- peak_data[non_qc_indices, , drop = FALSE]
  groups_filtered <- metadata$Group[non_qc_indices]

  # Validate group levels and arrange if specified
  groups_factor <- .process_group_levels(groups_filtered, arrangeLevels)

  # Check minimum samples per group
  .check_min_samples_per_group(groups_factor, min_samples_per_group)

  # Handle zero/negative values
  df_adjusted <- .handle_zero_negative_values(df_filtered, epsilon)
  fc_results$data_shifted <- as.data.frame(df_adjusted)

  # Add group summary to results
  fc_results$group_summary <- .create_group_summary(groups_factor, qc_filter_result)

  # Perform pairwise comparisons
  comparison_results <- .perform_pairwise_comparisons(df_adjusted, groups_factor, sortFC)

  # Add comparison results to main results list
  fc_results <- c(fc_results, comparison_results$results)
  fc_results$comparison_matrix <- comparison_results$comparison_matrix

  # Add analysis metadata
  fc_results$analysis_info <- list(
    n_comparisons = comparison_results$n_comparisons,
    n_features = ncol(df_adjusted),
    n_samples_analyzed = nrow(df_adjusted),
    groups_analyzed = levels(groups_factor)
  )

  return(fc_results)
}

#' Validate inputs for fold change analysis
#' @noRd
.validate_foldchange_inputs <- function(data, arrangeLevels, sortFC, qc_patterns,
                                        min_samples_per_group, epsilon) {

  # Check data structure
  if (!is.list(data)) {
    stop("'data' must be a list object from preprocessing functions.")
  }

  required_components <- c("data_scaledPCA_rsdFiltered_varFiltered", "Metadata")
  missing_components <- setdiff(required_components, names(data))
  if (length(missing_components) > 0) {
    stop(paste("Missing required components in 'data':",
               paste(missing_components, collapse = ", ")))
  }

  # Check arrangeLevels
  if (!is.null(arrangeLevels) && !is.character(arrangeLevels)) {
    stop("'arrangeLevels' must be a character vector or NULL.")
  }

  # Check logical parameters
  if (!is.logical(sortFC) || length(sortFC) != 1) {
    stop("'sortFC' must be a single logical value (TRUE or FALSE).")
  }

  # Check qc_patterns
  if (!is.character(qc_patterns) || length(qc_patterns) == 0) {
    stop("'qc_patterns' must be a non-empty character vector.")
  }

  # Check numeric parameters
  if (!is.numeric(min_samples_per_group) || min_samples_per_group < 1) {
    stop("'min_samples_per_group' must be a positive integer.")
  }

  if (!is.numeric(epsilon) || epsilon <= 0) {
    stop("'epsilon' must be a positive numeric value.")
  }
}

#' Extract peak data from input object
#' @noRd
.extract_peak_data <- function(data) {
  peak_data <- data$data_scaledPCA_rsdFiltered_varFiltered

  if (is.null(peak_data)) {
    stop("Peak data component is NULL.")
  }

  if (!is.data.frame(peak_data) && !is.matrix(peak_data)) {
    stop("Peak data must be a data.frame or matrix.")
  }

  # Convert to matrix for faster computation
  if (is.data.frame(peak_data)) {
    peak_data <- as.matrix(peak_data)
  }

  if (nrow(peak_data) == 0 || ncol(peak_data) == 0) {
    stop("Peak data is empty.")
  }

  if (!is.numeric(peak_data)) {
    stop("Peak data must contain only numeric values.")
  }

  return(peak_data)
}

#' Extract metadata from input object
#' @noRd
.extract_metadata <- function(data) {
  metadata <- data$Metadata

  if (is.null(metadata) || !is.data.frame(metadata)) {
    stop("Metadata must be a data.frame.")
  }

  if (!"Group" %in% colnames(metadata)) {
    stop("Metadata must contain a 'Group' column.")
  }

  if (nrow(metadata) == 0) {
    stop("Metadata is empty.")
  }

  return(metadata)
}

#' Filter QC samples from metadata
#' @noRd
.filter_qc_samples <- function(metadata, qc_patterns) {
  qc_indices <- metadata$Group %in% qc_patterns
  non_qc_indices <- which(!qc_indices)

  n_qc_removed <- sum(qc_indices)
  n_samples_remaining <- length(non_qc_indices)

  message(sprintf("Filtered out %d QC samples. %d samples remaining for analysis.",
                  n_qc_removed, n_samples_remaining))

  return(list(
    qc_indices = which(qc_indices),
    non_qc_indices = non_qc_indices,
    n_qc_removed = n_qc_removed,
    n_samples_remaining = n_samples_remaining
  ))
}

#' Process group levels and arrange if specified
#' @noRd
.process_group_levels <- function(groups_filtered, arrangeLevels) {
  unique_groups <- unique(groups_filtered)

  if (is.null(arrangeLevels)) {
    # Default alphabetical sorting
    groups_factor <- factor(groups_filtered, levels = sort(unique_groups))
  } else {
    # Validate arrangeLevels
    missing_levels <- setdiff(unique_groups, arrangeLevels)
    extra_levels <- setdiff(arrangeLevels, unique_groups)

    if (length(missing_levels) > 0) {
      stop(paste("Groups present in data but missing from 'arrangeLevels':",
                 paste(missing_levels, collapse = ", ")))
    }

    if (length(extra_levels) > 0) {
      warning(paste("Groups in 'arrangeLevels' but not in data:",
                    paste(extra_levels, collapse = ", ")))
    }

    # Filter arrangeLevels to only include groups present in data
    valid_levels <- intersect(arrangeLevels, unique_groups)
    groups_factor <- factor(groups_filtered, levels = valid_levels)
  }

  return(groups_factor)
}

#' Check minimum samples per group requirement
#' @noRd
.check_min_samples_per_group <- function(groups_factor, min_samples_per_group) {
  group_counts <- table(groups_factor)
  insufficient_groups <- group_counts < min_samples_per_group

  if (any(insufficient_groups)) {
    insufficient_names <- names(group_counts)[insufficient_groups]
    insufficient_counts <- group_counts[insufficient_groups]

    error_msg <- paste(
      sprintf("Groups with insufficient samples (< %d):", min_samples_per_group),
      paste(sprintf("%s: %d", insufficient_names, insufficient_counts), collapse = ", ")
    )
    stop(error_msg)
  }

  if (length(levels(groups_factor)) < 2) {
    stop("At least 2 groups are required for fold change analysis.")
  }
}

#' Handle zero and negative values in the data
#' @noRd
.handle_zero_negative_values <- function(df_filtered, epsilon) {
  min_value <- min(df_filtered, na.rm = TRUE)

  if (is.na(min_value)) {
    stop("Data contains only NA values.")
  }

  if (min_value <= 0) {
    # Calculate shift value to make minimum = epsilon
    shift_value <- epsilon - min_value
    df_adjusted <- df_filtered + shift_value

    message(sprintf("Added %.2e to all values to handle zero/negative values. New minimum: %.2e",
                    shift_value, min(df_adjusted, na.rm = TRUE)))
  } else {
    df_adjusted <- df_filtered
  }

  return(df_adjusted)
}

#' Create group summary information
#' @noRd
.create_group_summary <- function(groups_factor, qc_filter_result) {
  group_counts <- table(groups_factor)

  summary_df <- data.frame(
    Group = names(group_counts),
    Sample_Count = as.numeric(group_counts),
    stringsAsFactors = FALSE
  )

  # Add QC information
  summary_df$QC_Samples_Removed <- qc_filter_result$n_qc_removed
  summary_df$Total_Samples_Available <- qc_filter_result$n_samples_remaining + qc_filter_result$n_qc_removed

  return(summary_df)
}

#' Perform all pairwise group comparisons
#' @noRd
.perform_pairwise_comparisons <- function(df_adjusted, groups_factor, sortFC) {
  unique_groups <- levels(groups_factor)
  n_groups <- length(unique_groups)

  if (n_groups < 2) {
    stop("Need at least 2 groups for comparison.")
  }

  group_combinations <- combn(unique_groups, 2, simplify = FALSE)
  n_comparisons <- length(group_combinations)

  message(sprintf("Performing %d pairwise comparisons between %d groups.",
                  n_comparisons, n_groups))

  # Initialize results storage
  comparison_results <- vector("list", n_comparisons)
  comparison_names <- character(n_comparisons)
  comparison_matrix <- matrix(NA, nrow = n_comparisons, ncol = 3,
                              dimnames = list(NULL, c("Comparison", "Group1", "Group2")))

  # Perform comparisons
  for (i in seq_along(group_combinations)) {
    pair <- group_combinations[[i]]
    group1 <- pair[1]
    group2 <- pair[2]
    comparison_label <- paste(group1, "vs.", group2, sep = " ")

    # Store comparison info
    comparison_matrix[i, ] <- c(comparison_label, group1, group2)
    comparison_names[i] <- paste0("data_combined_", comparison_label)

    # Calculate fold changes
    fc_result <- .calculate_fold_change_pair(df_adjusted, groups_factor, group1, group2, sortFC)
    comparison_results[[i]] <- fc_result
  }

  # Name the results list
  names(comparison_results) <- comparison_names

  return(list(
    results = comparison_results,
    comparison_matrix = as.data.frame(comparison_matrix),
    n_comparisons = n_comparisons
  ))
}

#' Calculate fold change for a specific pair of groups
#' @noRd
.calculate_fold_change_pair <- function(df_adjusted, groups_factor, group1, group2, sortFC) {
  # Get indices for each group
  idx_group1 <- which(groups_factor == group1)
  idx_group2 <- which(groups_factor == group2)

  # Calculate group means efficiently using vectorized operations
  if (length(idx_group1) == 1) {
    mean_group1 <- df_adjusted[idx_group1, ]
  } else {
    mean_group1 <- colMeans(df_adjusted[idx_group1, , drop = FALSE], na.rm = TRUE)
  }

  if (length(idx_group2) == 1) {
    mean_group2 <- df_adjusted[idx_group2, ]
  } else {
    mean_group2 <- colMeans(df_adjusted[idx_group2, , drop = FALSE], na.rm = TRUE)
  }

  # Calculate fold change (avoid division by zero)
  fold_change <- ifelse(mean_group2 == 0, Inf, mean_group1 / mean_group2)

  # Calculate log2 fold change
  log2_fc <- log2(fold_change)

  # Create results data frame
  feature_names <- colnames(df_adjusted)
  if (is.null(feature_names)) {
    feature_names <- paste0("Feature_", seq_len(ncol(df_adjusted)))
  }

  combined_fc <- data.frame(
    Feature = feature_names,
    fold_change = fold_change,
    log2_fc = log2_fc,
    mean_group1 = mean_group1,
    mean_group2 = mean_group2,
    stringsAsFactors = FALSE,
    row.names = feature_names
  )

  # Rename so that the mean has the group name
  colnames(combined_fc)[4:5] <- c(paste0("mean_for_fc_", group1), paste0("mean_for_fc_", group2))

  # Sort if requested
  if (sortFC) {
    combined_fc <- combined_fc[order(combined_fc$fold_change, decreasing = TRUE), ]
  }

  return(combined_fc)
}
