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

  # ========== INPUT VALIDATION ==========

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

  if (fcDown >= fcUP) {
    stop("fcDown must be less than fcUP", call. = FALSE)
  }

  if (!is.null(arrangeLevels) && !is.character(arrangeLevels)) {
    stop("arrangeLevels must be a character vector or NULL", call. = FALSE)
  }

  if (verbose) message("Initializing volcano plot analysis...")

  # ========== INITIALIZE RESULTS STRUCTURE ==========

  volcano_results <- list(
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
    ),
    VolcanoPlots = list(),
    Summary = list()
  )

  # ========== EXTRACT AND VALIDATE GROUPS ==========

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

  # ========== GENERATE PAIRWISE COMPARISONS ==========

  group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)

  if (verbose) {
    message(sprintf("Processing %d group comparison(s):", length(group_combinations)))
    for (i in seq_along(group_combinations)) {
      message(sprintf("  %d. %s vs. %s", i, group_combinations[[i]][1], group_combinations[[i]][2]))
    }
  }

  # ========== PROCESS EACH COMPARISON ==========

  for (i in seq_along(group_combinations)) {
    group_pair <- group_combinations[[i]]
    group1 <- group_pair[1]
    group2 <- group_pair[2]

    if (verbose) message(sprintf("Processing comparison %d/%d: %s vs. %s",
                                 i, length(group_combinations), group1, group2))

    # --- Find matching fold change data ---
    possible_names <- c(
      paste(group1, "vs.", group2, sep = " "),
      paste(group2, "vs.", group1, sep = " ")
    )

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

    if (is.null(fc_data)) {
      if (verbose) warning(sprintf("No fold change data found for %s vs. %s", group1, group2))
      next
    }

    # Remove 'Feature' column if it exists
    if ("Feature" %in% colnames(fc_data)) {
      fc_data$Feature <- NULL
    }

    # --- Process fold change data ---
    fc_processed <- fc_data %>%
      tibble::rownames_to_column("Feature") %>%
      dplyr::arrange(.data$Feature)

    # --- Process comparative analysis data ---
    ca_processed <- CAData$results %>%
      tibble::rownames_to_column("Feature") %>%
      dplyr::arrange(.data$Feature)

    # # --- Merge and classify data ---
    # required_fc_cols <- c("Feature", "fold_change", "log2_fc")
    # required_ca_cols <- c("Feature", "omnibus_p_value")
    #
    # if (!all(required_fc_cols %in% colnames(fc_processed))) {
    #   missing_cols <- setdiff(required_fc_cols, colnames(fc_processed))
    #   if (verbose) warning(sprintf("Missing required columns in fold change data for %s: %s",
    #                                group_name, paste(missing_cols, collapse = ", ")))
    #   next
    # }
    #
    # if (!all(required_ca_cols %in% colnames(ca_processed))) {
    #   missing_cols <- setdiff(required_ca_cols, colnames(ca_processed))
    #   if (verbose) warning(sprintf("Missing required columns in comparative analysis data for %s: %s",
    #                                group_name, paste(missing_cols, collapse = ", ")))
    #   next
    # }
    #
    # volcano_data <- dplyr::full_join(fc_processed, ca_processed, by = "Feature") %>%
    #   dplyr::rename(
    #     Fold_Change = fold_change,
    #     Log2_Fold_Change = log2_fc,
    #     P_Value = `omnibus_p_value`,
    #     Adjusted_P_Value = `omnibus_p_value`
    #   ) %>%
    #   dplyr::rowwise() %>%
    #   dplyr::mutate(
    #     # Handle potential NA values
    #     Fold_Change = ifelse(is.na(.data$Fold_Change), 1, .data$Fold_Change),
    #     Adjusted_P_Value = min(dplyr::c_across(dplyr::starts_with("posthoc_p_"))),
    #
    #     # Classify significance
    #     Significance = dplyr::case_when(
    #       round(.data$Fold_Change, 2) >= fcUP & round(.data$Adjusted_P_Value, 3) < adjpvalue ~ "Upregulated",
    #       round(.data$Fold_Change, 2) <= fcDown & round(.data$Adjusted_P_Value, 3) < adjpvalue ~ "Downregulated",
    #       TRUE ~ "Not Significant"
    #     ),
    #
    #     # Create log10 p-value for plotting
    #     Neg_Log10_Adj_P = -log10(pmax(.data$Adjusted_P_Value, .Machine$double.eps))
    #   ) %>%
    #   dplyr::arrange(dplyr::desc(.data$Fold_Change))

    # --- Merge and classify data ---
    required_fc_cols <- c("Feature", "fold_change", "log2_fc")
    required_ca_cols <- c("Feature", "omnibus_p_value")

    if (!all(required_fc_cols %in% colnames(fc_processed))) {
      missing_cols <- setdiff(required_fc_cols, colnames(fc_processed))
      if (verbose) warning(sprintf("Missing required columns in fold change data for %s: %s",
                                   group_name, paste(missing_cols, collapse = ", ")))
      next
    }

    if (!all(required_ca_cols %in% colnames(ca_processed))) {
      missing_cols <- setdiff(required_ca_cols, colnames(ca_processed))
      if (verbose) warning(sprintf("Missing required columns in comparative analysis data for %s: %s",
                                   group_name, paste(missing_cols, collapse = ", ")))
      next
    }

    # volcano_data <- dplyr::full_join(fc_processed, ca_processed, by = "Feature") %>%
    #   dplyr::rename(
    #     Fold_Change = fold_change,
    #     Log2_Fold_Change = log2_fc,
    #     P_Value = `omnibus_p_value`
    #   ) %>%
    #   dplyr::rowwise() %>%
    #   dplyr::mutate(
    #     # Handle potential NA values
    #     Fold_Change = ifelse(is.na(.data$Fold_Change), 1, .data$Fold_Change),
    #
    #     # For 2 groups: use omnibus p-value directly
    #     # For 3+ groups: use minimum of post-hoc p-values
    #     Adjusted_P_Value = {
    #       posthoc_cols <- dplyr::c_across(dplyr::starts_with("posthoc_p_"))
    #       if (length(posthoc_cols) > 0 && any(!is.na(posthoc_cols))) {
    #         min(posthoc_cols, na.rm = TRUE)
    #       } else {
    #         .data$P_Value
    #       }
    #     },
    #
    #     # Classify significance
    #     Significance = dplyr::case_when(
    #       round(.data$Fold_Change, 2) >= fcUP & round(.data$Adjusted_P_Value, 3) < adjpvalue ~ "Upregulated",
    #       round(.data$Fold_Change, 2) <= fcDown & round(.data$Adjusted_P_Value, 3) < adjpvalue ~ "Downregulated",
    #       TRUE ~ "Not Significant"
    #     ),
    #
    #     # Create log10 p-value for plotting
    #     Neg_Log10_Adj_P = -log10(pmax(.data$Adjusted_P_Value, .Machine$double.eps))
    #   ) %>%
    #   dplyr::arrange(dplyr::desc(.data$Fold_Change))

    volcano_data <- dplyr::full_join(fc_processed, ca_processed, by = "Feature") %>%
      dplyr::rename(
        Fold_Change = fold_change,
        Log2_Fold_Change = log2_fc,
        P_Value = `omnibus_p_value`
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        # Handle potential NA values
        Fold_Change = ifelse(is.na(.data$Fold_Change), 1, .data$Fold_Change),

        # For 2 groups: use omnibus p-value directly
        # For 3+ groups: use minimum of post-hoc p-values
        Adjusted_P_Value = {
          posthoc_cols <- dplyr::c_across(dplyr::starts_with("posthoc_p_"))
          if (length(posthoc_cols) > 0 && any(!is.na(posthoc_cols))) {
            min(posthoc_cols, na.rm = TRUE)
          } else {
            .data$P_Value
          }
        },

        # Create rounded versions for consistent classification and display
        FC_rounded = round(.data$Fold_Change, 2),
        AdjP_rounded = round(.data$Adjusted_P_Value, 3),

        # Classify significance using rounded values
        Significance = dplyr::case_when(
          .data$FC_rounded >= fcUP & .data$AdjP_rounded < adjpvalue ~ "Upregulated",
          .data$FC_rounded <= fcDown & .data$AdjP_rounded < adjpvalue ~ "Downregulated",
          TRUE ~ "Not Significant"
        ),

        # Create log10 p-value for plotting
        Neg_Log10_Adj_P = -log10(pmax(.data$Adjusted_P_Value, .Machine$double.eps))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dplyr::desc(.data$Fold_Change))

    if (nrow(volcano_data) == 0) {
      if (verbose) message(sprintf("No features found for %s", group_name))
      next
    }

    # --- Create filtered dataset ---
    filtered_data <- volcano_data %>%
      dplyr::filter(.data$Significance != "Not Significant")

    # # --- Generate summary statistics ---
    # summary_table <- volcano_data %>%
    #   dplyr::count(.data$Significance, name = "Count") %>%
    #   dplyr::mutate(Percentage = round(.data$Count / sum(.data$Count) * 100, 2))
    #
    # n_up <- sum(volcano_data$Significance == "Upregulated")
    # n_down <- sum(volcano_data$Significance == "Downregulated")
    #
    # summary_stats <- list(
    #   Comparison = group_name,
    #   Total_Features = nrow(volcano_data),
    #   Significant_Features = sum(volcano_data$Significance != "Not Significant"),
    #   Upregulated_Count = n_up,
    #   Downregulated_Count = n_down,
    #   Not_Significant_Count = sum(volcano_data$Significance == "Not Significant"),
    #   Summary_Table = summary_table,
    #   Max_Log2_FC = max(abs(volcano_data$Log2_Fold_Change), na.rm = TRUE),
    #   Min_Adj_P_Value = min(volcano_data$Adjusted_P_Value, na.rm = TRUE)
    # )
    #
    # # --- Create volcano plot ---
    # n_total <- nrow(volcano_data)
    # fc_max <- max(volcano_data$Fold_Change, na.rm = TRUE)
    # fc_min <- min(volcano_data$Fold_Change, na.rm = TRUE)
    #
    # plot_obj <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = .data$Log2_Fold_Change, y = .data$Neg_Log10_Adj_P)) +
    #   ggplot2::geom_point(
    #     ggplot2::aes(color = .data$Significance),
    #     alpha = 0.6,
    #     size = 1.2
    #   ) +
    #   ggplot2::scale_color_manual(
    #     values = c(
    #       "Not Significant" = "#808080",
    #       "Upregulated" = "#E31A1C",
    #       "Downregulated" = "#1F78B4"
    #     ),
    #     name = "Regulation"
    #   ) +
    #   ggplot2::theme_minimal() +
    #   ggplot2::theme(
    #     plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    #     plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, color = "gray40"),
    #     axis.title = ggplot2::element_text(size = 12),
    #     axis.text = ggplot2::element_text(size = 10),
    #     legend.title = ggplot2::element_text(size = 11),
    #     legend.text = ggplot2::element_text(size = 10),
    #     panel.grid.minor = ggplot2::element_blank()
    #   ) +
    #   ggplot2::labs(
    #     title = sprintf("Volcano Plot: %s", group_name),
    #     subtitle = sprintf(
    #       "Features: %d total | %d upregulated | %d downregulated | Min FC = %.2f | Max FC = %.2f \nThresholds: %.2f >= log2FC >= %.2f, adj. p-value < %.3f",
    #       n_total, n_up, n_down, fc_min, fc_max, log2(fcDown), log2(fcUP), adjpvalue
    #     ),
    #     x = "Log2 Fold Change",
    #     y = "-Log10 (Adjusted P-value)",
    #     caption = sprintf("FC thresholds: >=%.2f (up) or <=%.2f (down)", fcUP, fcDown)
    #   ) +
    #   ggplot2::geom_hline(
    #     yintercept = -log10(adjpvalue),
    #     linetype = "dashed",
    #     color = "gray50",
    #     linewidth = 0.8,
    #     alpha = 0.8
    #   ) +
    #   ggplot2::geom_vline(
    #     xintercept = c(log2(fcDown), log2(fcUP)),
    #     linetype = "dashed",
    #     color = "gray50",
    #     linewidth = 0.8,
    #     alpha = 0.8
    #   )

    # --- Generate summary statistics ---
    summary_table <- volcano_data %>%
      dplyr::count(.data$Significance, name = "Count") %>%
      dplyr::mutate(Percentage = round(.data$Count / sum(.data$Count) * 100, 2))

    n_up <- sum(volcano_data$Significance == "Upregulated")
    n_down <- sum(volcano_data$Significance == "Downregulated")

    summary_stats <- list(
      Comparison = group_name,
      Total_Features = nrow(volcano_data),
      Significant_Features = sum(volcano_data$Significance != "Not Significant"),
      Upregulated_Count = n_up,
      Downregulated_Count = n_down,
      Not_Significant_Count = sum(volcano_data$Significance == "Not Significant"),
      Summary_Table = summary_table,
      Max_Log2_FC = max(abs(volcano_data$Log2_Fold_Change), na.rm = TRUE),
      Min_Adj_P_Value = min(volcano_data$Adjusted_P_Value, na.rm = TRUE)
    )

    # --- Create volcano plot ---
    n_total <- nrow(volcano_data)
    # Use rounded values for min/max display to match classification logic
    fc_max <- max(volcano_data$FC_rounded, na.rm = TRUE)
    fc_min <- min(volcano_data$FC_rounded, na.rm = TRUE)

    plot_obj <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = .data$Log2_Fold_Change, y = .data$Neg_Log10_Adj_P)) +
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
          "Features: %d total | %d upregulated | %d downregulated | Min FC = %.2f | Max FC = %.2f \nThresholds: %.2f >= log2FC >= %.2f, adj. p-value < %.3f",
          n_total, n_up, n_down, fc_min, fc_max, log2(fcDown), log2(fcUP), adjpvalue
        ),
        x = "Log2 Fold Change",
        y = "-Log10 (Adjusted P-value)",
        caption = sprintf("FC thresholds: >=%.2f (up) or <=%.2f (down)", fcUP, fcDown)
      ) +
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

    if (show_plots) {
      print(plot_obj)
    }

    # --- Store results ---
    volcano_results$VolcanoPlots[[group_name]] <- plot_obj
    volcano_results[[paste0("VolcanoData_", gsub(" ", "_", group_name))]] <- volcano_data
    volcano_results[[paste0("VolcanoData_filtered_", gsub(" ", "_", group_name))]] <- filtered_data
    volcano_results$Summary[[group_name]] <- summary_stats
  }

  # ========== FINALIZE ==========

  if (verbose) {
    message("Volcano plot analysis completed successfully!")
    message(sprintf("Generated %d plot(s) and datasets", length(volcano_results$VolcanoPlots)))
  }

  return(volcano_results)
}
