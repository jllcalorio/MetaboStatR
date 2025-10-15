#' Perform Area Under the Receiver Operating Characteristic Curve Analysis
#'
#' @description
#' Performs comprehensive Area Under the Receiver Operating Characteristic (AUROC)
#' analysis for metabolomics data. This function integrates results from data
#' preprocessing, dimension reduction (OPLS-DA), fold change analysis, and
#' comparative analysis to identify and visualize discriminative metabolic features
#' between different groups. The function automatically filters features based on
#' Variable Importance in Projection (VIP) scores, fold change thresholds, and
#' statistical significance, then generates ROC curves for the most discriminative features.
#'
#' @param data_PP List. Results from the \code{performPreprocessingPeakData} function.
#'   Must contain \code{$Metadata} with \code{$Group} column and
#'   \code{$data_scaledPCA_varFiltered} or \code{$data_scaledPCA_merged} matrix.
#'   The latter is for the case when replicates were merged.
#' @param data_DR List. Results from the \code{performDimensionReduction} function
#'   (OPLS-DA). Must contain VIP score data frames with naming pattern
#'   \code{data_VIPScores_[group1] vs. [group2]}.
#' @param data_FC List. Results from the \code{performFoldChange} function. Must
#'   contain fold change data frames with naming pattern
#'   \code{data_combined_[group1] vs. [group2]}.
#' @param data_CA List. Results from the \code{performComparativeAnalysis} function.
#'   Must contain \code{$results} data frame with adjusted p-values.
#' @param arrangeLevels Character vector or NULL. Specifies the order of groups
#'   for analysis (e.g., \code{c("Control", "Case1", "Case2")}). When NULL,
#'   groups are sorted alphabetically. It is recommended to arrange from least
#'   severe to most severe condition (e.g., "Control", "Mild", "Severe").
#' @param VIPmin Numeric. Minimum Variable Importance in Projection (VIP) score
#'   threshold for feature selection. Features with VIP < VIPmin are excluded.
#'   Default: 1.0.
#' @param fcUP Numeric. Minimum fold change threshold for up-regulated features.
#'   Features with fold change >= fcUP are considered up-regulated. Default: 2.0.
#' @param fcDown Numeric. Maximum fold change threshold for down-regulated features.
#'   Features with fold change <= fcDown are considered down-regulated. Default: 0.5.
#' @param adjpvalue Numeric. Adjusted p-value threshold for statistical significance.
#'   Features with adjusted p-value >= adjpvalue are excluded. Default: 0.05.
#' @param direction Character. Direction parameter for ROC analysis. Options:
#'   \itemize{
#'     \item \code{"auto"}: Automatically determines optimal direction
#'     \item \code{">"}: Use when predictor values for control group > case group
#'     \item \code{"<"}: Use when predictor values for control group < case group
#'   }
#'   See \code{\link[pROC]{roc}} for details. Default: "auto".
#' @param top_n Integer. Number of top features (by AUC) to display in ROC plots.
#'   Must be positive. Default: 5.
#' @param plot_iden_met Character vector or NULL. Specific metabolite names to plot
#'   individually. When provided, generates separate ROC plots for each specified
#'   metabolite regardless of filtering criteria. Default: NULL.
#' @param confidence_level Numeric. Confidence level for AUC confidence intervals.
#'   Must be between 0 and 1. Default: 0.95.
#' @param min_group_size Integer. Minimum number of samples required per group
#'   for ROC analysis. Default: 3.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{FunctionOrigin}: Character indicating the source function
#'     \item \code{data_Merged_[comparison]}: Merged data frames for each group comparison
#'     \item \code{data_Filtered_[comparison]}: Filtered data matrices for each comparison
#'     \item \code{data_results_[comparison]}: AUROC results data frames
#'     \item \code{plots}: List of combined ROC plots for top features
#'     \item \code{plots_identifiedMetabolites}: Individual plots for specified metabolites (if requested)
#'     \item \code{summary}: Overall summary statistics
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input data and parameters
#'   \item Excludes quality control (QC) samples
#'   \item Generates all possible pairwise group combinations
#'   \item For each comparison, merges VIP, fold change, and statistical data
#'   \item Filters features based on specified thresholds
#'   \item Computes ROC curves and AUC values with confidence intervals
#'   \item Generates publication-ready plots
#'   \item Optionally creates individual plots for specified metabolites
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' auroc_results <- perform_AUROC(
#'   data_PP = preprocessing_results,
#'   data_DR = dimension_reduction_results,
#'   data_FC = fold_change_results,
#'   data_CA = comparative_analysis_results
#' )
#'
#' # Advanced usage with custom parameters
#' auroc_results <- perform_AUROC(
#'   data_PP = preprocessing_results,
#'   data_DR = dimension_reduction_results,
#'   data_FC = fold_change_results,
#'   data_CA = comparative_analysis_results,
#'   arrangeLevels = c("Control", "Mild", "Severe"),
#'   VIPmin = 1.5,
#'   fcUP = 1.5,
#'   fcDown = 0.67,
#'   adjpvalue = 0.01,
#'   top_n = 10,
#'   plot_iden_met = c("Metabolite1", "Metabolite2")
#' )
#' }
#'
#' @seealso
#' \code{\link[pROC]{roc}}, \code{\link[pROC]{auc}}, \code{\link[pROC]{ci}}
#'
#' @importFrom dplyr arrange filter select any_of full_join desc
#' @importFrom ggplot2 ggplot geom_line aes labs theme_minimal theme scale_color_manual geom_abline
#' @importFrom pROC roc auc ci
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices rainbow
#' @importFrom utils combn head
#' @importFrom stats median
#'
#' @references Tom Fawcett (2006) “An introduction to ROC analysis”. Pattern Recognition Letters 27, 861–874. DOI: doi: 10.1016/j.patrec.2005.10.010. (for roc, auc)
#' @references Xavier Robin, Natacha Turck, Alexandre Hainard, et al. (2011) “pROC: an open-source package for R and S+ to analyze and compare ROC curves”. BMC Bioinformatics, 7, 77. DOI: doi: 10.1186/1471-2105-12-77. (for roc, auc, ci)
#' @references David J. Hand and Robert J. Till (2001). A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. Machine Learning 45(2), p. 171–186. DOI: doi: 10.1023/A:1010920819831. (for auc)
#' @references Donna Katzman McClish (1989) “Analyzing a Portion of the ROC Curve”. Medical Decision Making 9(3), 190–195. DOI: doi: 10.1177/0272989X8900900307. (for auc)
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. (for median)
#'
#' @export
perform_AUROC <- function(
    data_PP,
    data_DR,
    data_FC,
    data_CA,
    arrangeLevels = NULL,
    VIPmin = 1.0,
    fcUP = 2.0,
    fcDown = 0.5,
    adjpvalue = 0.05,
    direction = "auto",
    top_n = 5L,
    plot_iden_met = NULL,
    confidence_level = 0.95,
    min_group_size = 3L
) {

  # Validate inputs
  .validate_auroc_inputs(data_PP, data_DR, data_FC, data_CA, arrangeLevels,
                         VIPmin, fcUP, fcDown, adjpvalue, direction, top_n,
                         confidence_level, min_group_size)

  # Initialize results list
  auroc_results <- list(
    FunctionOrigin = "perform_AUROC",
    timestamp = Sys.time(),
    parameters = list(
      arrangeLevels = arrangeLevels,
      VIPmin = VIPmin,
      fcUP = fcUP,
      fcDown = fcDown,
      adjpvalue = adjpvalue,
      direction = direction,
      top_n = top_n,
      confidence_level = confidence_level,
      min_group_size = min_group_size
    ),
    plots = list(),
    plots_identifiedMetabolites = list(),
    summary = list()
  )

  # Process group data
  group_info <- .process_group_data_auroc(data_PP, arrangeLevels, min_group_size)
  groups <- group_info$groups
  unique_groups <- group_info$unique_groups
  non_qc_indices <- group_info$non_qc_indices

  # Generate group combinations
  group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)

  if (length(group_combinations) == 0) {
    warning("No valid group combinations found for analysis.")
    return(auroc_results)
  }

  # Initialize summary statistics
  processed_comparisons <- 0
  failed_comparisons <- character(0)
  total_features_analyzed <- 0

  # Process each group combination
  for (group_pair in group_combinations) {

    tryCatch({
      # Note: In "arrangeLevels" parameter, the input should be c(Control, Case1, Case2, ...) as required by AUROC
      # Set up group comparison (reverse for ROC: control vs case)
      group1 <- group_pair[2]  # case (the 2nd item in the pair)
      group2 <- group_pair[1]  # control (the 1st item in the pair)
      group_name <- paste(group1, "vs.", group2, sep = " ")
      # group_name <- paste(group2, "vs.", group1, sep = " ") # Interchanged group1 and group2 positions

      message("\nProcessing comparison: ", group_name)

      # Check for required data
      vip_key <- paste0("data_VIPScores_", group_name)
      if (is.null(data_DR[[vip_key]])) {
        message("Skipping '", group_name, "' - missing VIP data (no OPLS-DA model)")
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      # Extract and merge datasets
      merged_data <- .extract_and_merge_data(data_DR, data_FC, data_CA, group_name)

      # Save the merged data in a list
      auroc_results[[paste0("merged_data_", group_name)]] <- merged_data

      if (is.null(merged_data)) {
        message("Skipping '", group_name, "' - failed to merge data")
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      # Filter features based on criteria
      filtered_data <- .filter_features(merged_data, VIPmin, fcUP, fcDown, adjpvalue)

      # Save the filtered data in a list
      auroc_results[[paste0("filtered_data_", group_name)]] <- filtered_data

      if (nrow(filtered_data) == 0) {
        message("No features passed filtering for ", group_name)
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      # Store merged and filtered data
      auroc_results[[paste0("data_Merged_", group_name)]] <- merged_data

      # Prepare data matrix
      feature_data <- .prepare_feature_data(data_PP, filtered_data$Feature,
                                            groups, non_qc_indices)
      auroc_results[[paste0("data_Filtered_", group_name)]] <- feature_data

      # Perform AUROC analysis
      auroc_df <- .perform_auroc_analysis(feature_data, group1, group2,
                                          direction, confidence_level)

      if (nrow(auroc_df) == 0) {
        message("No valid AUROC results for ", group_name)
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      auroc_results[[paste0("data_results_", group_name)]] <- auroc_df

      # Create main ROC plot
      main_plot <- .create_roc_plot(feature_data, auroc_df, group1, group2,
                                    group_name, direction, top_n)
      auroc_results$plots[[group_name]] <- main_plot

      # Display plot
      print(main_plot)

      # Create individual metabolite plots if requested
      if (!is.null(plot_iden_met)) {
        individual_plots <- .create_individual_plots(data_PP, data_DR, data_FC,
                                                     data_CA, plot_iden_met,
                                                     group1, group2, group_name,
                                                     groups, non_qc_indices,
                                                     direction, confidence_level)
        if (length(individual_plots) > 0) {
          auroc_results$plots_identifiedMetabolites[[group_name]] <- individual_plots
        }
      }

      processed_comparisons <- processed_comparisons + 1
      total_features_analyzed <- total_features_analyzed + nrow(auroc_df)

    }, error = function(e) {
      warning("Error processing ", group_name, ": ", e$message)
      failed_comparisons <<- c(failed_comparisons, group_name)
    })
  }

  # Add summary information
  auroc_results$summary <- list(
    total_comparisons = length(group_combinations),
    processed_comparisons = processed_comparisons,
    failed_comparisons = failed_comparisons,
    total_features_analyzed = total_features_analyzed,
    unique_groups = unique_groups,
    samples_per_group = table(groups)
  )

  message("\nAUROC analysis completed: ", processed_comparisons, "/",
          length(group_combinations), " comparisons processed successfully")

  return(auroc_results)
}

# Helper Functions --------------------------------------------------------

#' Validate AUROC function inputs
#' @keywords internal
.validate_auroc_inputs <- function(data_PP, data_DR, data_FC, data_CA, arrangeLevels,
                                   VIPmin, fcUP, fcDown, adjpvalue, direction, top_n,
                                   confidence_level, min_group_size) {

  # Check required list inputs
  required_lists <- list(data_PP = data_PP, data_DR = data_DR,
                         data_FC = data_FC, data_CA = data_CA)

  for (name in names(required_lists)) {
    if (!is.list(required_lists[[name]])) {
      stop(paste(name, "must be a list"))
    }
  }

  # Check data_PP structure
  if (is.null(data_PP$Metadata) || is.null(data_PP$Metadata$Group)) {
    stop("data_PP must contain $Metadata$Group")
  }

  # if (is.null(data_PP$data_scaledPCA_varFiltered)) {
  #   stop("data_PP must contain $data_scaledPCA_varFiltered")
  # }

  if(is.null(data_PP$data_scaledPLS_merged)) {
    if(is.null(data_PP$data_scaledPCA_varFiltered)) {
      stop("data_PP must contain $data_scaledPCA_varFiltered")
    }
  }

  # Check data_CA structure
  if (is.null(data_CA$results)) {
    stop("data_CA must contain $results")
  }

  # Validate numeric parameters
  numeric_params <- list(
    VIPmin = VIPmin, fcUP = fcUP, fcDown = fcDown, adjpvalue = adjpvalue,
    confidence_level = confidence_level, min_group_size = min_group_size
  )

  for (name in names(numeric_params)) {
    if (!is.numeric(numeric_params[[name]]) || length(numeric_params[[name]]) != 1) {
      stop(paste(name, "must be a single numeric value"))
    }
  }

  # Validate specific numeric ranges
  if (VIPmin < 0) stop("VIPmin must be non-negative")
  if (fcUP <= 1) stop("fcUP must be greater than 1")
  if (fcDown >= 1) stop("fcDown must be less than 1")
  if (fcDown <= 0) stop("fcDown must be positive")
  if (adjpvalue <= 0 || adjpvalue >= 1) stop("adjpvalue must be between 0 and 1")
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1")
  }
  if (min_group_size < 2) stop("min_group_size must be at least 2")

  # Validate direction parameter
  if (!direction %in% c("auto", ">", "<")) {
    stop("direction must be 'auto', '>', or '<'")
  }

  # Validate top_n
  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1 || top_n != round(top_n)) {
    stop("top_n must be a positive integer")
  }

  # Validate arrangeLevels if provided
  if (!is.null(arrangeLevels) && !is.character(arrangeLevels)) {
    stop("arrangeLevels must be a character vector or NULL")
  }

  # Check for data matrix - handle both merged and non-merged cases
  has_varFiltered <- !is.null(data_PP$data_scaledPCA_varFiltered)
  has_merged <- !is.null(data_PP$data_scaledPCA_merged)
  auto_merge_replicates <- isTRUE(data_PP$Parameters$auto_merge_replicates)

  if (auto_merge_replicates) {
    if (!has_merged) {
      stop("data_PP$Parameters$auto_merge_replicates is TRUE but data_PP$data_scaledPCA_merged is missing")
    }
  } else {
    if (!has_varFiltered) {
      stop("data_PP$Parameters$auto_merge_replicates is FALSE but data_PP$data_scaledPCA_varFiltered is missing")
    }
  }
}

#' Process group data and handle QC samples
#' @keywords internal
.process_group_data_auroc <- function(data_PP, arrangeLevels, min_group_size) {

  # Identify QC and non-QC samples
  qc_indices <- data_PP$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices <- !qc_indices

  if (sum(non_qc_indices) == 0) {
    stop("No non-QC samples found in data")
  }

  # Extract non-QC groups
  groups <- data_PP$Metadata$Group[non_qc_indices]

  # Handle group arrangement
  if (!is.null(arrangeLevels)) {
    missing_groups <- setdiff(arrangeLevels, unique(groups))
    extra_groups <- setdiff(unique(groups), arrangeLevels)

    if (length(missing_groups) > 0) {
      stop("arrangeLevels contains groups not found in data: ",
           paste(missing_groups, collapse = ", "))
    }

    if (length(extra_groups) > 0) {
      warning("Data contains groups not in arrangeLevels: ",
              paste(extra_groups, collapse = ", "))
    }

    groups <- factor(groups, levels = arrangeLevels)
    unique_groups <- arrangeLevels
  } else {
    unique_groups <- sort(unique(groups))
  }

  # Check minimum group sizes
  group_counts <- table(groups)
  small_groups <- names(group_counts)[group_counts < min_group_size]

  if (length(small_groups) > 0) {
    warning("Groups with fewer than ", min_group_size, " samples: ",
            paste(small_groups, collapse = ", "))
  }

  # Filter to groups with sufficient size
  valid_groups <- names(group_counts)[group_counts >= min_group_size]
  unique_groups <- intersect(unique_groups, valid_groups)

  if (length(unique_groups) < 2) {
    stop("Need at least 2 groups with minimum ", min_group_size, " samples each")
  }

  return(list(
    groups = groups,
    unique_groups = unique_groups,
    non_qc_indices = non_qc_indices,
    group_counts = group_counts
  ))
}

#' Extract and merge VIP, fold change, and comparative analysis data
#' @keywords internal
.extract_and_merge_data <- function(data_DR, data_FC, data_CA, group_name) {

  tryCatch({
    # Extract datasets with error checking
    vip_data <- data_DR[[paste0("data_VIPScores_", group_name)]]
    fc_data <- data_FC[[paste0("data_combined_", group_name)]]
    ca_results <- data_CA$results

    # Remove the "Feature" column in Fold Change results
    fc_data$Feature <- NULL

    # Validate extracted data
    if (is.null(vip_data) || nrow(vip_data) == 0) {
      stop("VIP data is empty or missing")
    }

    if (is.null(fc_data) || nrow(fc_data) == 0) {
      stop("Fold change data is empty or missing")
    }

    if (is.null(ca_results) || nrow(ca_results) == 0) {
      stop("Comparative analysis results are empty or missing")
    }

    # Prepare data frames with Feature column
    vip_data <- vip_data %>%
      dplyr::arrange(Feature)

    fc_data <- fc_data %>%
      cbind(Feature = rownames(.), .) %>%
      `rownames<-`(NULL) %>%
      dplyr::arrange(Feature)

    ca_results <- ca_results %>%
      cbind(Feature = rownames(.), .) %>%
      `rownames<-`(NULL) %>%
      dplyr::arrange(Feature)

    # Merge datasets
    merged_data <- vip_data %>%
      dplyr::full_join(fc_data, by = "Feature") %>%
      dplyr::full_join(ca_results, by = "Feature")

    # Remove rows with missing critical data
    merged_data <- merged_data %>%
      filter(!is.na(VIP), !is.na(fold_change), !is.na(adj_p_value))

    return(merged_data)

  }, error = function(e) {
    message("Error merging data for ", group_name, ": ", e$message)
    return(NULL)
  })
}

#' Filter features based on specified criteria
#' @keywords internal
.filter_features <- function(merged_data, VIPmin, fcUP, fcDown, adjpvalue) {

  filtered_data <- merged_data %>%
    dplyr::filter(
      VIP > VIPmin,
      (fold_change >= fcUP | fold_change <= fcDown),
      adj_p_value < adjpvalue
    )

  return(filtered_data)
}

# REPLACED WITH UPDATED CODES BELOW

# #' Prepare feature data matrix for AUROC analysis
# #' @keywords internal
# .prepare_feature_data <- function(data_PP, feature_names, groups, non_qc_indices) {
#
#   # Get the data matrix
#   data_matrix <- data_PP$data_scaledPCA_varFiltered[non_qc_indices, , drop = FALSE]
#
#   # Check which features are available
#   available_features <- intersect(feature_names, colnames(data_matrix))
#
#   if (length(available_features) == 0) {
#     stop("No requested features found in data matrix")
#   }
#
#  if (length(available_features) < length(feature_names)) {
#    missing_features <- setdiff(feature_names, available_features)
#     warning("Some features not found in data matrix: ",
#             paste(missing_features, collapse = ", "))
#   }
#
#   # Create feature data with groups
#   feature_data <- data_matrix %>%
#     dplyr::select(dplyr::any_of(available_features)) %>%
#     cbind(Groups = groups, .)
#
#   return(feature_data)
# }

#' Prepare feature data matrix for AUROC analysis
#' @keywords internal
.prepare_feature_data <- function(data_PP, feature_names, groups, non_qc_indices) {

  tryCatch({
    # Get the appropriate data matrix based on auto_merge_replicates setting
    data_matrix <- .get_data_matrix_AUROC(data_PP)

    # Apply non-QC indices filtering
    data_matrix <- data_matrix[non_qc_indices, , drop = FALSE]

    # Check which features are available
    available_features <- intersect(feature_names, colnames(data_matrix))

    if (length(available_features) == 0) {
      stop("No requested features found in data matrix")
    }

    if (length(available_features) < length(feature_names)) {
      missing_features <- setdiff(feature_names, available_features)
      warning("Some features not found in data matrix: ",
              paste(missing_features, collapse = ", "))
    }

    # Create feature data with groups
    feature_data <- data_matrix %>%
      dplyr::select(dplyr::any_of(available_features)) %>%
      cbind(Groups = groups, .)

    return(feature_data)

  }, error = function(e) {
    stop("Error preparing feature data: ", e$message)
  })
}

#' Perform AUROC analysis for filtered features
#' @keywords internal
.perform_auroc_analysis <- function(feature_data, group1, group2, direction, confidence_level) {

  # Get feature names (exclude Groups column)
  feature_names <- setdiff(colnames(feature_data), "Groups")

  # Initialize results data frame
  auroc_df <- data.frame(
    Feature = character(),
    Group_Comparison = character(),
    AUROC = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    stringsAsFactors = FALSE
  )

  group_name <- paste(group1, "vs.", group2, sep = " ")

  for (feature_name in feature_names) {

    tryCatch({
      # Compute ROC curve
      roc_obj <- pROC::roc(
        response = feature_data$Groups,
        predictor = feature_data[[feature_name]],
        levels = c(group1, group2),
        direction = direction,
        quiet = TRUE
      )

      # Get AUC and confidence interval
      auc_val <- as.numeric(pROC::auc(roc_obj))
      ci_vals <- as.numeric(pROC::ci(roc_obj, conf.level = confidence_level))

      # Add to results
      auroc_df <- rbind(auroc_df, data.frame(
        Feature = feature_name,
        Group_Comparison = group_name,
        AUROC = auc_val,
        CI_Lower = ci_vals[1],
        CI_Upper = ci_vals[3],
        stringsAsFactors = FALSE
      ))

    }, error = function(e) {
      warning("Failed to compute ROC for feature ", feature_name, ": ", e$message)
    })
  }

  # Sort by AUC (descending)
  auroc_df <- auroc_df %>%
    dplyr::arrange(dplyr::desc(AUROC))

  return(auroc_df)
}

#' Create ROC plot for top features - FIXED AXES
#' @keywords internal
.create_roc_plot <- function(feature_data, auroc_df, group1, group2, group_name,
                             direction, top_n) {

  # Get top features
  top_features <- head(auroc_df, min(top_n, nrow(auroc_df)))

  if (nrow(top_features) == 0) {
    warning("No features available for plotting")
    return(ggplot2::ggplot())
  }

  # Determine colors
  n_features <- nrow(top_features)
  if (n_features <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n_features), "Dark2")[1:n_features]
  } else {
    colors <- grDevices::rainbow(n_features)
  }

  # Initialize plot and data for legend ordering
  gg_roc <- ggplot2::ggplot()
  legend_data <- data.frame()

  # Add ROC curves
  for (i in seq_len(nrow(top_features))) {
    feature <- top_features$Feature[i]
    auc_val <- top_features$AUROC[i]
    ci_low <- top_features$CI_Lower[i]
    ci_up <- top_features$CI_Upper[i]

    # Compute ROC
    roc_obj <- pROC::roc(
      response = feature_data$Groups,
      predictor = feature_data[[feature]],
      levels = c(group1, group2),
      direction = direction,
      quiet = TRUE
    )

    # Create ROC data frame - FIXED: Correct axes assignment and sorting by x-axis
    roc_df <- data.frame(
      x = 1 - roc_obj$specificities,  # 1-Specificity on x-axis
      y = roc_obj$sensitivities,      # Sensitivity on y-axis
      Feature = sprintf("%s (AUC = %.3f) (95%% CI: %.3f-%.3f)",
                        feature, auc_val, ci_low, ci_up)
    ) %>%
      dplyr::arrange(x, y)  # Sort by x-axis (1-Specificity) then y-axis for smooth curves

    # Add to plot
    gg_roc <- gg_roc +
      ggplot2::geom_line(
        data = roc_df,
        ggplot2::aes(x = x, y = y, color = Feature),
        size = 1.2
      )

    # Store for legend ordering
    legend_data <- rbind(legend_data, data.frame(
      Feature = sprintf("%s (AUC = %.3f) (95%% CI: %.3f-%.3f)",
                        feature, auc_val, ci_low, ci_up),
      AUC = auc_val
    ))
  }

  # Order legend by AUC (descending)
  legend_data <- legend_data %>%
    dplyr::arrange(dplyr::desc(AUC))

  # Finalize plot - FIXED: Correct axis labels and reference line
  final_plot <- gg_roc +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray50", size = 0.8) +  # Reference line
    ggplot2::scale_color_manual(
      values = colors,
      breaks = legend_data$Feature
    ) +
    ggplot2::labs(
      title = sprintf("Top %d Features' ROC Curves - %s", min(top_n, nrow(top_features)), group_name),
      x = "1 - Specificity (False Positive Rate)",  # FIXED: Correct x-axis label
      y = "Sensitivity (True Positive Rate)",       # FIXED: Correct y-axis label
      color = "Feature (AUC, 95% CI)"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = c(0.8, 0.2),
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.width = ggplot2::unit(1.5, "lines")
    )

  return(final_plot)
}

#' Create individual plots for specified metabolites - FIXED AXES
#' @keywords internal
.create_individual_plots <- function(data_PP, data_DR, data_FC, data_CA, plot_iden_met,
                                     group1, group2, group_name, groups, non_qc_indices,
                                     direction, confidence_level) {

  individual_plots <- list()

  tryCatch({
    # Get merged data (without filtering)
    merged_data_full <- .extract_and_merge_data(data_DR, data_FC, data_CA, group_name)

    if (is.null(merged_data_full)) {
      return(individual_plots)
    }

    # Find available metabolites
    available_metabolites <- intersect(plot_iden_met, merged_data_full$Feature)

    if (length(available_metabolites) == 0) {
      warning("None of the specified metabolites found in data for ", group_name)
      return(individual_plots)
    }

    # Prepare data matrix - UPDATED CALL
    feature_data <- .prepare_feature_data(data_PP, available_metabolites,
                                          groups, non_qc_indices)

    # Create individual plots
    for (metabolite in available_metabolites) {

      tryCatch({
        # Compute ROC
        roc_obj <- pROC::roc(
          response = feature_data$Groups,
          predictor = feature_data[[metabolite]],
          levels = c(group1, group2),
          direction = direction,
          quiet = TRUE
        )

        # Get AUC and CI
        auc_val <- as.numeric(pROC::auc(roc_obj))
        ci_vals <- as.numeric(pROC::ci(roc_obj, conf.level = confidence_level))

        # Create plot data - FIXED: Correct axes assignment and sorting by x-axis
        roc_df <- data.frame(
          x = 1 - roc_obj$specificities,  # 1-Specificity on x-axis
          y = roc_obj$sensitivities,      # Sensitivity on y-axis
          Label = sprintf("%s (AUC = %.3f, 95%% CI: %.3f-%.3f)",
                          metabolite, auc_val, ci_vals[1], ci_vals[3])
        ) %>%
          dplyr::arrange(x, y)  # Sort by x-axis (1-Specificity) then y-axis for smooth curves

        # Create plot - FIXED: Correct axis labels and reference line
        individual_plot <- ggplot2::ggplot() +
          ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                               color = "gray50", size = 0.8) +  # Reference line
          ggplot2::geom_line(
            data = roc_df,
            ggplot2::aes(x = x, y = y, color = Label),
            size = 1.2
          ) +
          ggplot2::labs(
            title = sprintf("ROC Curve - %s", group_name),
            x = "1 - Specificity (False Positive Rate)",  # FIXED: Correct x-axis label
            y = "Sensitivity (True Positive Rate)",       # FIXED: Correct y-axis label
            color = "Metabolite (AUC, 95% CI)"
          ) +
          ggplot2::xlim(0, 1) +
          ggplot2::ylim(0, 1) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
            legend.text = ggplot2::element_text(size = 10)
          )

        # Display and store plot
        print(individual_plot)
        individual_plots[[metabolite]] <- individual_plot

      }, error = function(e) {
        warning("Failed to create plot for metabolite ", metabolite, ": ", e$message)
      })
    }

  }, error = function(e) {
    warning("Error creating individual plots for ", group_name, ": ", e$message)
  })

  return(individual_plots)
}

# # Additional Utility Functions --------------------------------------------
#
# #' Extract AUROC Results
# #'
# #' @description
# #' Extract comprehensive AUROC results from perform_AUROC output including
# #' feature names, AUC values, confidence intervals, and additional metrics.
# #' This function provides flexible options to extract results for specific
# #' comparisons or all comparisons, with optional filtering and sorting.
# #'
# #' @param auroc_results List. Results from the \code{perform_AUROC} function.
# #' @param comparison Character vector or NULL. Specific group comparison(s) to extract
# #'   (e.g., \code{c("Case1 vs. Control", "Case2 vs. Control")}). If NULL, extracts
# #'   all comparisons. Default: NULL.
# #' @param min_auc Numeric. Minimum AUC threshold for filtering results.
# #'   Features with AUC < min_auc are excluded. Default: 0 (no filtering).
# #' @param max_auc Numeric. Maximum AUC threshold for filtering results.
# #'   Features with AUC > max_auc are excluded. Default: 1 (no filtering).
# #' @param top_n Integer or NULL. Number of top features (by AUC) to return per
# #'   comparison. If NULL, returns all features. Default: NULL.
# #' @param include_merged_data Logical. Whether to include merged data (VIP, fold
# #'   change, statistical results) for each comparison. Default: FALSE.
# #' @param sort_by Character. Sorting criteria for results. Options:
# #'   \itemize{
# #'     \item \code{"auc_desc"}: Sort by AUC descending (highest first)
# #'     \item \code{"auc_asc"}: Sort by AUC ascending (lowest first)
# #'     \item \code{"feature"}: Sort by feature name alphabetically
# #'     \item \code{"comparison"}: Sort by comparison name alphabetically
# #'   }
# #'   Default: "auc_desc".
# #'
# #' @return A list containing the following components:
# #'   \itemize{
# #'     \item \code{results}: Data frame with AUROC results including Feature,
# #'       Group_Comparison, AUROC, CI_Lower, CI_Upper, and additional metrics
# #'     \item \code{summary}: Summary statistics for extracted results
# #'     \item \code{parameters}: Parameters used for extraction
# #'     \item \code{merged_data}: List of merged data frames for each comparison
# #'       (if \code{include_merged_data = TRUE})
# #'   }
# #'
# #' @details
# #' The extracted results data frame includes:
# #' \enumerate{
# #'   \item Feature: Metabolite/feature name
# #'   \item Group_Comparison: Comparison groups (e.g., "Case vs. Control")
# #'   \item AUROC: Area Under the ROC Curve value
# #'   \item CI_Lower: Lower bound of confidence interval
# #'   \item CI_Upper: Upper bound of confidence interval
# #'   \item CI_Width: Width of confidence interval (CI_Upper - CI_Lower)
# #'   \item Performance_Category: Categorical assessment of AUC performance
# #' }
# #'
# #' Performance categories:
# #' \itemize{
# #'   \item "Excellent" (AUC ≥ 0.9)
# #'   \item "Good" (0.8 ≤ AUC < 0.9)
# #'   \item "Fair" (0.7 ≤ AUC < 0.8)
# #'   \item "Poor" (0.6 ≤ AUC < 0.7)
# #'   \item "Fail" (AUC < 0.6)
# #' }
# #'
# #' @author John Lennon L. Calorio
# #'
# #' @examples
# #' \dontrun{
# #' # Extract all AUROC results
# #' all_results <- extract_auroc_results(auroc_results)
# #'
# #' # Extract results for specific comparison
# #' case_results <- extract_auroc_results(
# #'   auroc_results,
# #'   comparison = "Case vs. Control"
# #' )
# #'
# #' # Extract top 10 features with AUC > 0.7
# #' top_results <- extract_auroc_results(
# #'   auroc_results,
# #'   min_auc = 0.7,
# #'   top_n = 10,
# #'   include_merged_data = TRUE
# #' )
# #'
# #' # Extract results sorted by feature name
# #' sorted_results <- extract_auroc_results(
# #'   auroc_results,
# #'   sort_by = "feature"
# #' )
# #' }
# #'
# #' @seealso \code{\link{perform_AUROC}}, \code{\link{get_auroc_summary}}
# #'
# #' @export
extract_auroc_results <- function(
    auroc_results,
    comparison = NULL,
    min_auc = 0,
    max_auc = 1,
    top_n = NULL,
    include_merged_data = FALSE,
    sort_by = "auc_desc"
) {

  # Validate inputs
  if (!"perform_AUROC" %in% auroc_results$FunctionOrigin) {
    stop("Input must be results from perform_AUROC function")
  }

  # Validate parameters
  if (!is.numeric(min_auc) || !is.numeric(max_auc) || min_auc < 0 || max_auc > 1 || min_auc >= max_auc) {
    stop("min_auc and max_auc must be numeric values between 0 and 1, with min_auc < max_auc")
  }

  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1 || top_n != round(top_n))) {
    stop("top_n must be a positive integer or NULL")
  }

  if (!sort_by %in% c("auc_desc", "auc_asc", "feature", "comparison")) {
    stop("sort_by must be one of: 'auc_desc', 'auc_asc', 'feature', 'comparison'")
  }

  # Find all AUROC result data frames
  result_names <- names(auroc_results)[grepl("^data_results_", names(auroc_results))]

  if (length(result_names) == 0) {
    warning("No AUROC results found")
    return(list(
      results = data.frame(),
      summary = list(),
      parameters = list(),
      merged_data = list()
    ))
  }

  # Filter by specific comparisons if requested
  if (!is.null(comparison)) {
    comparison_names <- paste0("data_results_", comparison)
    valid_comparisons <- intersect(comparison_names, result_names)

    if (length(valid_comparisons) == 0) {
      stop("None of the specified comparisons found in results. Available comparisons: ",
           paste(gsub("^data_results_", "", result_names), collapse = ", "))
    }

    result_names <- valid_comparisons

    # Check for invalid comparisons
    invalid_comparisons <- setdiff(comparison_names, result_names)
    if (length(invalid_comparisons) > 0) {
      warning("Some specified comparisons not found: ",
              paste(gsub("^data_results_", "", invalid_comparisons), collapse = ", "))
    }
  }

  # Extract and combine all results
  all_results <- data.frame()
  merged_data_list <- list()

  for (result_name in result_names) {
    comparison_name <- gsub("^data_results_", "", result_name)
    result_df <- auroc_results[[result_name]]

    if (!is.null(result_df) && nrow(result_df) > 0) {
      all_results <- rbind(all_results, result_df)

      # Include merged data if requested
      if (include_merged_data) {
        merged_name <- paste0("data_Merged_", comparison_name)
        if (!is.null(auroc_results[[merged_name]])) {
          merged_data_list[[comparison_name]] <- auroc_results[[merged_name]]
        }
      }
    }
  }

  if (nrow(all_results) == 0) {
    warning("No valid AUROC results found for specified comparisons")
    return(list(
      results = data.frame(),
      summary = list(),
      parameters = list(),
      merged_data = list()
    ))
  }

  # Apply AUC filtering
  filtered_results <- all_results[all_results$AUROC >= min_auc & all_results$AUROC <= max_auc, ]

  if (nrow(filtered_results) == 0) {
    warning("No results pass the AUC filtering criteria (min_auc: ", min_auc, ", max_auc: ", max_auc, ")")
    return(list(
      results = data.frame(),
      summary = list(),
      parameters = list(),
      merged_data = merged_data_list
    ))
  }

  # Add additional calculated columns
  filtered_results$CI_Width <- filtered_results$CI_Upper - filtered_results$CI_Lower

  # Add performance category
  filtered_results$Performance_Category <- cut(
    filtered_results$AUROC,
    breaks = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c("Fail", "Poor", "Fair", "Good", "Excellent"),
    include.lowest = TRUE,
    right = FALSE
  )

  # Apply sorting
  filtered_results <- switch(sort_by,
                             "auc_desc" = filtered_results[order(-filtered_results$AUROC), ],
                             "auc_asc" = filtered_results[order(filtered_results$AUROC), ],
                             "feature" = filtered_results[order(filtered_results$Feature), ],
                             "comparison" = filtered_results[order(filtered_results$Group_Comparison, -filtered_results$AUROC), ]
  )

  # Apply top_n filtering per comparison if specified
  if (!is.null(top_n)) {
    if (sort_by %in% c("feature", "comparison")) {
      # For non-AUC sorting, need to sort by AUC within each comparison first
      filtered_results <- do.call(rbind, lapply(split(filtered_results, filtered_results$Group_Comparison), function(df) {
        df_sorted <- df[order(-df$AUROC), ]
        head(df_sorted, min(top_n, nrow(df_sorted)))
      }))

      # Re-apply original sorting
      filtered_results <- switch(sort_by,
                                 "feature" = filtered_results[order(filtered_results$Feature), ],
                                 "comparison" = filtered_results[order(filtered_results$Group_Comparison, -filtered_results$AUROC), ]
      )
    } else {
      # For AUC-based sorting, apply top_n per comparison
      filtered_results <- do.call(rbind, lapply(split(filtered_results, filtered_results$Group_Comparison), function(df) {
        head(df, min(top_n, nrow(df)))
      }))
    }
  }

  # Reset row names
  rownames(filtered_results) <- NULL

  # Calculate summary statistics
  summary_stats <- list(
    total_features = nrow(filtered_results),
    unique_comparisons = length(unique(filtered_results$Group_Comparison)),
    comparisons = unique(filtered_results$Group_Comparison),
    auc_statistics = list(
      mean = mean(filtered_results$AUROC, na.rm = TRUE),
      median = median(filtered_results$AUROC, na.rm = TRUE),
      min = min(filtered_results$AUROC, na.rm = TRUE),
      max = max(filtered_results$AUROC, na.rm = TRUE),
      sd = sd(filtered_results$AUROC, na.rm = TRUE)
    ),
    performance_distribution = table(filtered_results$Performance_Category),
    ci_width_statistics = list(
      mean = mean(filtered_results$CI_Width, na.rm = TRUE),
      median = median(filtered_results$CI_Width, na.rm = TRUE),
      min = min(filtered_results$CI_Width, na.rm = TRUE),
      max = max(filtered_results$CI_Width, na.rm = TRUE)
    )
  )

  # Store extraction parameters
  extraction_params <- list(
    comparison = comparison,
    min_auc = min_auc,
    max_auc = max_auc,
    top_n = top_n,
    include_merged_data = include_merged_data,
    sort_by = sort_by,
    extraction_timestamp = Sys.time()
  )

  return(list(
    results = filtered_results,
    summary = summary_stats,
    parameters = extraction_params,
    merged_data = merged_data_list
  ))
}

#' Get AUROC summary statistics
#'
#' @description Extract summary statistics from AUROC results
#' @param auroc_results List returned by perform_AUROC
#' @return Data frame with summary statistics
#' @export
get_auroc_summary <- function(auroc_results) {

  if (!"perform_AUROC" %in% auroc_results$FunctionOrigin) {
    stop("Input must be results from perform_AUROC function")
  }

  # Extract all results data frames
  result_names <- names(auroc_results)[grepl("^data_results_", names(auroc_results))]

  if (length(result_names) == 0) {
    warning("No AUROC results found")
    return(data.frame())
  }

  # Combine all results
  all_results <- do.call(rbind, lapply(result_names, function(name) {
    auroc_results[[name]]
  }))

  # Calculate summary statistics
  summary_stats <- data.frame(
    Total_Comparisons = auroc_results$summary$total_comparisons,
    Processed_Comparisons = auroc_results$summary$processed_comparisons,
    Total_Features = nrow(all_results),
    Mean_AUC = mean(all_results$AUROC, na.rm = TRUE),
    Median_AUC = median(all_results$AUROC, na.rm = TRUE),
    Min_AUC = min(all_results$AUROC, na.rm = TRUE),
    Max_AUC = max(all_results$AUROC, na.rm = TRUE),
    Features_AUC_Above_0.7 = sum(all_results$AUROC > 0.7, na.rm = TRUE),
    Features_AUC_Above_0.8 = sum(all_results$AUROC > 0.8, na.rm = TRUE),
    Features_AUC_Above_0.9 = sum(all_results$AUROC > 0.9, na.rm = TRUE)
  )

  return(summary_stats)
}

#' Plot AUROC distribution
#'
#' @description Create histogram of AUC values across all comparisons
#' @param auroc_results List returned by perform_AUROC
#' @param bins Integer. Number of histogram bins
#' @return ggplot2 object
#' @export
plot_auroc_distribution <- function(auroc_results, bins = 20) {

  if (!"perform_AUROC" %in% auroc_results$FunctionOrigin) {
    stop("Input must be results from perform_AUROC function")
  }

  # Extract all AUC values
  result_names <- names(auroc_results)[grepl("^data_results_", names(auroc_results))]

  if (length(result_names) == 0) {
    warning("No AUROC results found")
    return(ggplot2::ggplot())
  }

  all_aucs <- do.call(c, lapply(result_names, function(name) {
    auroc_results[[name]]$AUROC
  }))

  # Create histogram
  auc_df <- data.frame(AUC = all_aucs)

  p <- ggplot2::ggplot(auc_df, ggplot2::aes(x = AUC)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue", alpha = 0.7, color = "white") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
    ggplot2::geom_vline(xintercept = 0.7, linetype = "dashed", color = "orange", size = 1) +
    ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed", color = "green", size = 1) +
    ggplot2::annotate("text", x = 0.5, y = Inf, label = "Random (0.5)",
                      vjust = 1.5, hjust = -0.1, color = "red") +
    ggplot2::annotate("text", x = 0.7, y = Inf, label = "Good (0.7)",
                      vjust = 1.5, hjust = -0.1, color = "orange") +
    ggplot2::annotate("text", x = 0.8, y = Inf, label = "Excellent (0.8)",
                      vjust = 1.5, hjust = -0.1, color = "green") +
    ggplot2::labs(
      title = "Distribution of AUC Values",
      subtitle = paste("Total features:", length(all_aucs)),
      x = "Area Under Curve (AUC)",
      y = "Frequency"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}

#' Get appropriate data matrix based on auto_merge_replicates setting
#' @keywords internal
.get_data_matrix_AUROC <- function(data_PP) {

  auto_merge_replicates <- isTRUE(data_PP$Parameters$auto_merge_replicates)

  if (auto_merge_replicates) {
    if (is.null(data_PP$data_scaledPCA_merged)) {
      stop("auto_merge_replicates is TRUE but data_scaledPCA_merged is not available")
    }
    return(data_PP$data_scaledPCA_merged)
  } else {
    if (is.null(data_PP$data_scaledPCA_varFiltered)) {
      stop("auto_merge_replicates is FALSE but data_scaledPCA_varFiltered is not available")
    }
    return(data_PP$data_scaledPCA_varFiltered)
  }
}
