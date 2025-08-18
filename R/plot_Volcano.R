#' Create Volcano Plots for Differential Expression Analysis
#'
#' @description
#' Generates volcano plots to visualize the results of differential expression analysis
#' by plotting log2 fold changes against negative log10 adjusted p-values. The function
#' creates plots for all pairwise group comparisons and applies significance thresholds
#' to highlight upregulated, downregulated, and non-significant features.
#'
#' @param PPData A list object returned by \code{perform_PreprocessingPeakData()}.
#'   Must contain a \code{Metadata} component with a \code{Group} column.
#' @param FCData A list object returned by \code{perform_FoldChange()}.
#'   Must contain fold change data for group comparisons.
#' @param CAData A list object returned by \code{perform_ComparativeAnalysis()}.
#'   Must contain a \code{results} component with p-values and adjusted p-values.
#' @param arrangeLevels A character vector specifying the order of group levels.
#'   Format: \code{c('group1', 'group2', ...)}. If \code{NULL} (default),
#'   groups are sorted alphabetically. Recommended order: control, case1, case2
#'   (e.g., Control, Non-severe dengue, Severe dengue).
#' @param fcUP A numeric value specifying the upper fold change threshold for
#'   significance. Features with fold change >= \code{fcUP} are considered
#'   upregulated. Default: 2.
#' @param fcDown A numeric value specifying the lower fold change threshold for
#'   significance. Features with fold change <= \code{fcDown} are considered
#'   downregulated. Default: 0.5.
#' @param adjpvalue A numeric value specifying the adjusted p-value threshold
#'   for significance. Default: 0.05.
#' @param show_plots A logical value indicating whether to display plots.
#'   Default: \code{TRUE}.
#' @param verbose A logical value indicating whether to print progress messages.
#'   Default: \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{FunctionOrigin}{Character string identifying the source function}
#'   \item{Parameters}{List of input parameters used}
#'   \item{VolcanoPlots}{Named list of ggplot2 objects, one for each comparison}
#'   \item{VolcanoData_[Comparison]}{Data frame with complete volcano plot data for each comparison}
#'   \item{VolcanoData_filtered_[Comparison]}{Data frame with significant features only for each comparison}
#'   \item{Summary}{Summary statistics for each comparison}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input data structures and parameters
#'   \item Extracts non-QC samples from the metadata
#'   \item Creates all possible pairwise group comparisons
#'   \item Merges fold change and statistical test results
#'   \item Applies significance thresholds to classify features
#'   \item Generates volcano plots with customizable thresholds
#'   \item Returns comprehensive results including plots and data
#' }
#'
#' Features are classified as:
#' \itemize{
#'   \item \strong{Upregulated}: Fold change >= \code{fcUP} AND adjusted p-value < \code{adjpvalue}
#'   \item \strong{Downregulated}: Fold change <= \code{fcDown} AND adjusted p-value < \code{adjpvalue}
#'   \item \strong{Not Significant}: All other features
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' volcano_results <- plot_Volcano(
#'   PPData = preprocess_results,
#'   FCData = foldchange_results,
#'   CAData = comparative_results
#' )
#'
#' # With custom thresholds and group ordering
#' volcano_results <- plot_Volcano(
#'   PPData = preprocess_results,
#'   FCData = foldchange_results,
#'   CAData = comparative_results,
#'   arrangeLevels = c("Control", "Treatment1", "Treatment2"),
#'   fcUP = 1.5,
#'   fcDown = 0.67,
#'   adjpvalue = 0.01
#' )
#'
#' # Access specific results
#' plot_obj <- volcano_results$VolcanoPlots[["Control vs. Treatment1"]]
#' sig_features <- volcano_results$VolcanoData_filtered_Control_vs._Treatment1
#' }
#'
#' @seealso
#' \code{\link{perform_PreprocessingPeakData}}, \code{\link{perform_FoldChange}},
#' \code{\link{perform_ComparativeAnalysis}}
#'
#' @author John Lennon L. Calorio
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom utils combn
#' @export
plot_Volcano <- function(PPData,
                         FCData,
                         CAData,
                         arrangeLevels = NULL,
                         fcUP = 2,
                         fcDown = 0.5,
                         adjpvalue = 0.05,
                         show_plots = TRUE,
                         verbose = TRUE) {

  # Validate inputs
  .validate_volcano_inputs(PPData, FCData, CAData, arrangeLevels, fcUP, fcDown, adjpvalue)

  if (verbose) message("Initializing volcano plot analysis...")

  # Initialize results structure
  volcano_results <- .initialize_volcano_results(FCData, CAData, fcUP, fcDown, adjpvalue)

  # Extract and process groups
  groups_info <- .extract_groups(PPData, arrangeLevels, verbose)

  # Generate all pairwise comparisons
  group_combinations <- utils::combn(groups_info$unique_groups, 2, simplify = FALSE)

  if (verbose) {
    message(sprintf("Processing %d group comparison(s):", length(group_combinations)))
    for (i in seq_along(group_combinations)) {
      message(sprintf("  %d. %s vs. %s", i, group_combinations[[i]][1], group_combinations[[i]][2]))
    }
  }

  # Process each comparison
  volcano_results$VolcanoPlots <- list()
  volcano_results$Summary <- list()

  for (i in seq_along(group_combinations)) {
    group_pair <- group_combinations[[i]]

    if (verbose) message(sprintf("Processing comparison %d/%d: %s vs. %s",
                                 i, length(group_combinations), group_pair[1], group_pair[2]))

    # Process single comparison
    comparison_result <- .process_single_comparison(
      group_pair, FCData, CAData, fcUP, fcDown, adjpvalue, show_plots, verbose
    )

    if (!is.null(comparison_result)) {
      group_name <- comparison_result$group_name

      # Store results
      volcano_results$VolcanoPlots[[group_name]] <- comparison_result$plot
      volcano_results[[paste0("VolcanoData_", gsub(" ", "_", group_name))]] <- comparison_result$data
      volcano_results[[paste0("VolcanoData_filtered_", gsub(" ", "_", group_name))]] <- comparison_result$filtered_data
      volcano_results$Summary[[group_name]] <- comparison_result$summary
    }
  }

  if (verbose) {
    message("Volcano plot analysis completed successfully!")
    message(sprintf("Generated %d plot(s) and datasets", length(volcano_results$VolcanoPlots)))
  }

  return(volcano_results)
}

# Helper function: Validate inputs
.validate_volcano_inputs <- function(PPData, FCData, CAData, arrangeLevels, fcUP, fcDown, adjpvalue) {
  # Check required data structures
  if (!is.list(PPData) || is.null(PPData$Metadata) || is.null(PPData$Metadata$Group)) {
    stop("PPData must be a list with Metadata$Group component", call. = FALSE)
  }

  if (!is.list(FCData) || length(FCData) == 0) {
    stop("FCData must be a non-empty list", call. = FALSE)
  }

  if (!is.list(CAData) || is.null(CAData$results)) {
    stop("CAData must be a list with a 'results' component", call. = FALSE)
  }

  # Validate numeric parameters
  if (!is.numeric(fcUP) || length(fcUP) != 1 || fcUP <= 0) {
    stop("fcUP must be a positive numeric value", call. = FALSE)
  }

  if (!is.numeric(fcDown) || length(fcDown) != 1 || fcDown <= 0) {
    stop("fcDown must be a positive numeric value", call. = FALSE)
  }

  if (!is.numeric(adjpvalue) || length(adjpvalue) != 1 || adjpvalue <= 0 || adjpvalue >= 1) {
    stop("adjpvalue must be a numeric value between 0 and 1", call. = FALSE)
  }

  # Check logical relationship between thresholds
  if (fcDown >= fcUP) {
    stop("fcDown must be less than fcUP", call. = FALSE)
  }

  # Validate arrangeLevels if provided
  if (!is.null(arrangeLevels) && !is.character(arrangeLevels)) {
    stop("arrangeLevels must be a character vector or NULL", call. = FALSE)
  }
}

# Helper function: Initialize results structure
.initialize_volcano_results <- function(FCData, CAData, fcUP, fcDown, adjpvalue) {
  list(
    FunctionOrigin = "plot_Volcano",
    Parameters = list(
      FoldChange_Upper = fcUP,
      FoldChange_Down = fcDown,
      AdjustedPValue = adjpvalue,
      Log2_FoldChange_Upper = log2(fcUP),
      Log2_FoldChange_Down = log2(fcDown)
    ),
    InputData = list(
      FoldChangeData = FCData,
      ComparativeAnalysisData = CAData
    )
  )
}

# Helper function: Extract and validate groups
.extract_groups <- function(PPData, arrangeLevels, verbose) {
  # Identify QC samples
  qc_patterns <- c("SQC", "EQC", "QC", "POOLED", "BLANK")
  qc_indices <- grepl(paste(qc_patterns, collapse = "|"), PPData$Metadata$Group, ignore.case = TRUE)
  non_qc_indices <- !qc_indices

  if (sum(non_qc_indices) == 0) {
    stop("No non-QC samples found in the data", call. = FALSE)
  }

  # Extract groups
  groups <- PPData$Metadata$Group[non_qc_indices]

  if (verbose) message(sprintf("Found %d non-QC samples across %d group(s)",
                               sum(non_qc_indices), length(unique(groups))))

  # Validate and arrange levels
  unique_groups <- unique(groups)

  if (!is.null(arrangeLevels)) {
    missing_levels <- setdiff(arrangeLevels, unique_groups)
    extra_levels <- setdiff(unique_groups, arrangeLevels)

    if (length(missing_levels) > 0) {
      stop(sprintf("arrangeLevels contains groups not found in data: %s",
                   paste(missing_levels, collapse = ", ")), call. = FALSE)
    }

    if (length(extra_levels) > 0) {
      warning(sprintf("Data contains groups not in arrangeLevels: %s",
                      paste(extra_levels, collapse = ", ")))
    }

    unique_groups <- arrangeLevels[arrangeLevels %in% unique_groups]
  } else {
    unique_groups <- sort(unique_groups)
  }

  if (length(unique_groups) < 2) {
    stop("At least 2 groups are required for comparison", call. = FALSE)
  }

  list(groups = groups, unique_groups = unique_groups)
}

# Helper function: Process single comparison
.process_single_comparison <- function(group_pair, FCData, CAData, fcUP, fcDown, adjpvalue, show_plots, verbose) {
  group1 <- group_pair[1]
  group2 <- group_pair[2]

  # Try both possible group name formats
  possible_names <- c(
    paste(group1, "vs.", group2, sep = " "),
    paste(group2, "vs.", group1, sep = " ")
  )

  # Find matching fold change data
  fc_data <- NULL
  group_name <- NULL

  for (name in possible_names) {
    fc_key <- paste0("data_combined_", name)
    if (fc_key %in% names(FCData) && !is.null(FCData[[fc_key]])) {
      fc_data <- FCData[[fc_key]]
      group_name <- name
      break
    }
  }

  # Remove 'Feature' column
  fc_data$Feature <- NULL

  if (is.null(fc_data)) {
    if (verbose) warning(sprintf("No fold change data found for %s vs. %s", group1, group2))
    return(NULL)
  }

  # Process fold change data
  fc_processed <- fc_data %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::arrange(.data$Feature)

  # Process comparative analysis data
  ca_processed <- CAData$results %>%
    tibble::rownames_to_column("Feature") %>%
    dplyr::arrange(.data$Feature)

  # Merge datasets with error handling
  tryCatch({
    volcano_data <- .merge_and_classify_data(fc_processed, ca_processed, fcUP, fcDown, adjpvalue)
  }, error = function(e) {
    if (verbose) warning(sprintf("Error merging data for %s: %s", group_name, e$message))
    return(NULL)
  })

  if (nrow(volcano_data) == 0) {
    if (verbose) message(sprintf("No features found for %s", group_name))
    return(NULL)
  }

  # Create filtered dataset
  filtered_data <- volcano_data %>%
    dplyr::filter(.data$Significance != "Not Significant")

  # Generate summary statistics
  summary_stats <- .generate_summary_stats(volcano_data, group_name)

  # Create plot
  plot_obj <- .create_volcano_plot(volcano_data, group_name, fcUP, fcDown, adjpvalue)

  if (show_plots) {
    print(plot_obj)
  }

  list(
    group_name = group_name,
    data = volcano_data,
    filtered_data = filtered_data,
    plot = plot_obj,
    summary = summary_stats
  )
}

# Helper function: Merge and classify data
.merge_and_classify_data <- function(fc_data, ca_data, fcUP, fcDown, adjpvalue) {
  # Standardize column names for merging
  required_fc_cols <- c("Feature", "fold_change", "log2_fc")
  required_ca_cols <- c("Feature", "p_value", "adj_p_value")

  if (!all(required_fc_cols %in% colnames(fc_data))) {
    missing_cols <- setdiff(required_fc_cols, colnames(fc_data))
    stop(sprintf("Missing required columns in fold change data: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  if (!all(required_ca_cols %in% colnames(ca_data))) {
    missing_cols <- setdiff(required_ca_cols, colnames(ca_data))
    stop(sprintf("Missing required columns in comparative analysis data: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Merge and classify
  merged_data <- dplyr::full_join(fc_data, ca_data, by = "Feature") %>%
    dplyr::rename(
      Fold_Change = .data$fold_change,
      Log2_Fold_Change = .data$log2_fc,
      P_Value = .data$p_value,
      Adjusted_P_Value = .data$adj_p_value
    ) %>%
    dplyr::mutate(
      # Handle potential NA values
      Fold_Change = ifelse(is.na(.data$Fold_Change), 1, .data$Fold_Change),
      Adjusted_P_Value = ifelse(is.na(.data$Adjusted_P_Value), 1, .data$Adjusted_P_Value),

      # Classify significance
      Significance = dplyr::case_when(
        .data$Fold_Change >= fcUP & .data$Adjusted_P_Value < adjpvalue ~ "Upregulated",
        .data$Fold_Change <= fcDown & .data$Adjusted_P_Value < adjpvalue ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),

      # Create log10 p-value for plotting
      Neg_Log10_Adj_P = -log10(pmax(.data$Adjusted_P_Value, .Machine$double.eps))
    ) %>%
    dplyr::arrange(dplyr::desc(.data$Fold_Change))

  return(merged_data)
}

# Helper function: Generate summary statistics
.generate_summary_stats <- function(data, group_name) {
  summary_table <- data %>%
    dplyr::count(.data$Significance, name = "Count") %>%
    dplyr::mutate(Percentage = round(.data$Count / sum(.data$Count) * 100, 2))

  list(
    Comparison = group_name,
    Total_Features = nrow(data),
    Significant_Features = sum(data$Significance != "Not Significant"),
    Upregulated_Count = sum(data$Significance == "Upregulated"),
    Downregulated_Count = sum(data$Significance == "Downregulated"),
    Not_Significant_Count = sum(data$Significance == "Not Significant"),
    Summary_Table = summary_table,
    Max_Log2_FC = max(abs(data$Log2_Fold_Change), na.rm = TRUE),
    Min_Adj_P_Value = min(data$Adjusted_P_Value, na.rm = TRUE)
  )
}

# Helper function: Create volcano plot
.create_volcano_plot <- function(data, group_name, fcUP, fcDown, adjpvalue) {
  # Count significant features for subtitle
  n_up <- sum(data$Significance == "Upregulated")
  n_down <- sum(data$Significance == "Downregulated")
  n_total <- nrow(data)
  fc_max <- max(data$Fold_Change) # max(data[sapply(data$Fold_Change, is.numeric)], na.rm = TRUE)
  fc_min <- min(data$Fold_Change) # min(data[sapply(data$Fold_Change, is.numeric)], na.rm = TRUE)

  # Create the plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$Log2_Fold_Change, y = .data$Neg_Log10_Adj_P)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$Significance),
      alpha = 0.6,
      size = 1.2
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Not Significant" = "#808080",
        "Upregulated" = "#E31A1C",
        "Downregulated" = "#1F78B4"
      ),
      name = "Regulation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = sprintf("Volcano Plot: %s", group_name),
      subtitle = sprintf(
        # "Features: %d total | %d upregulated | %d downregulated\nThresholds: |log2FC| >= %.2f, adj. p-value < %.3f",
        # n_total, n_up, n_down, log2(fcUP), adjpvalue
        "Features: %d total | %d upregulated | %d downregulated | Min FC = %.2f | Max FC = %.2f \nThresholds: %.2f <= log2FC <= %.2f, adj. p-value < %.3f",
        n_total, n_up, n_down, fc_min, fc_max, log2(fcUP), log2(fcDown), adjpvalue
      ),
      x = "Log2 Fold Change",
      y = "-Log10 (Adjusted P-value)",
      caption = sprintf("FC thresholds: >=%.2f (up) or <=%.2f (down)", fcUP, fcDown)
    ) +
    # Add threshold lines
    ggplot2::geom_hline(
      yintercept = -log10(adjpvalue),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.8,
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = c(log2(fcDown), log2(fcUP)),
      linetype = "dashed",
      color = "gray50",
      linewidth = 0.8,
      alpha = 0.8
    )

  return(p)
}
