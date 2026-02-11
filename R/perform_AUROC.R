#' Perform Area Under the Receiver Operating Characteristic Curve Analysis
#'
#' @description
#' Performs comprehensive Area Under the Receiver Operating Characteristic (AUROC)
#' analysis for metabolomics data. This function integrates results from data
#' preprocessing, dimension reduction (OPLS-DA), fold change analysis, and
#' comparative analysis to identify and visualize discriminative metabolic features
#' between different groups.
#'
#' @param data_PP List. Results from the \code{performPreprocessingPeakData} function.
#' @param data_DR List. Results from the \code{performDimensionReduction} function.
#' @param data_FC List. Results from the \code{performFoldChange} function.
#' @param data_CA List. Results from the \code{performComparativeAnalysis} function.
#' @param arrangeLevels Character vector or NULL. Specifies the order of groups.
#' Should be arranged from control to case.
#' @param VIPmin Numeric. Minimum VIP score threshold. Default: 1.0.
#' @param fcUP Numeric. Minimum fold change threshold for up-regulation. Default: 2.0.
#' Checks after rounding off to 2 decimal places.
#' @param fcDown Numeric. Maximum fold change threshold for down-regulation. Default: 0.5.
#' Checks after rounding off to 2 decimal places.
#' @param adjpvalue Numeric. Adjusted p-value threshold. Default: 0.05.
#' @param direction Character. Direction parameter for ROC analysis. Default: "auto".
#' @param top_n Integer. Number of top features to display. Default: 5.
#' @param minimum_auc Numeric or NULL. Minimum AUC threshold for feature selection.
#' When specified, overrides 'top_n' and plots all features with AUC >= minimum_auc.
#' Must be between 0 and 1. Default: 0.8.
#' @param plot_iden_met Character vector or NULL. Specific metabolites to plot. Default: NULL.
#' @param confidence_level Numeric. Confidence level for AUC intervals. Default: 0.95.
#' @param min_group_size Integer. Minimum samples required per group. Default: 3.
#' @param skip_tests Character vector or NULL. Specifies which tests to skip when filtering features.
#' Options: "DR" (dimension reduction/VIP), "FC" (fold change), "CA" (comparative analysis).
#' When a test is skipped, its filtering criteria are ignored. For example, skip_tests = c("DR", "FC")
#' will only filter by comparative analysis p-value. Use skip_tests = c("DR", "FC", "CA") to
#' perform AUROC on all available features without any filtering. Default: NULL (use all tests).
#'
#' @return A list containing AUROC analysis results
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom dplyr arrange filter select any_of full_join desc
#' @importFrom ggplot2 ggplot geom_line aes labs theme_minimal theme scale_color_manual geom_abline
#' @importFrom pROC roc auc ci
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices rainbow
#' @importFrom utils combn head
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#' # Standard analysis using all tests
#' auroc_full <- perform_AUROC(
#'   data_PP = preprocessed_data,
#'   data_DR = dimension_reduction_data,
#'   data_FC = fold_change_data,
#'   data_CA = comparative_analysis_data
#' )
#'
#' # Skip dimension reduction when OPLS-DA component is not significant
#' auroc_no_dr <- perform_AUROC(
#'   data_PP = preprocessed_data,
#'   data_DR = dimension_reduction_data,
#'   data_FC = fold_change_data,
#'   data_CA = comparative_analysis_data,
#'   skip_tests = "DR"
#' )
#'
#' # Skip both DR and FC, use only statistical significance
#' auroc_stats_only <- perform_AUROC(
#'   data_PP = preprocessed_data,
#'   data_DR = dimension_reduction_data,
#'   data_FC = fold_change_data,
#'   data_CA = comparative_analysis_data,
#'   skip_tests = c("DR", "FC")
#' )
#'
#' # Perform AUROC on all features without filtering
#' auroc_all_features <- perform_AUROC(
#'   data_PP = preprocessed_data,
#'   data_DR = dimension_reduction_data,
#'   data_FC = fold_change_data,
#'   data_CA = comparative_analysis_data,
#'   skip_tests = c("DR", "FC", "CA")
#' )
#' }
#'
#' @export
perform_AUROC <- function(
    data_PP,
    data_DR,
    data_FC,
    data_CA,
    arrangeLevels    = NULL,
    VIPmin           = 1.0,
    fcUP             = 2.0,
    fcDown           = 0.5,
    adjpvalue        = 0.05,
    direction        = "auto",
    top_n            = 5L,
    minimum_auc      = 0.8,
    plot_iden_met    = NULL,
    confidence_level = 0.95,
    min_group_size   = 3L,
    skip_tests       = NULL
) {

  # ========== INPUT VALIDATION ==========
  # Check required list inputs
  if (!is.list(data_PP) || !is.list(data_DR) || !is.list(data_FC) || !is.list(data_CA)) {
    stop("data_PP, data_DR, data_FC, and data_CA must be lists")
  }

  if (is.null(data_PP$Metadata) || is.null(data_PP$Metadata$Group)) {
    stop("data_PP must contain $Metadata$Group")
  }

  if (is.null(data_PP$data_scaledPLS_merged)) {
    if (is.null(data_PP$data_scaledNONPLS_varFiltered)) {
      stop("data_PP must contain $data_scaledNONPLS_varFiltered")
    }
  }

  if (is.null(data_CA$results)) {
    stop("data_CA must contain $results")
  }

  # Validate numeric parameters
  if (!is.numeric(VIPmin) || VIPmin < 0) stop("VIPmin must be non-negative")
  if (!is.numeric(fcUP) || fcUP <= 1) stop("fcUP must be greater than 1")
  if (!is.numeric(fcDown) || fcDown >= 1 || fcDown <= 0) stop("fcDown must be between 0 and 1")
  if (!is.numeric(adjpvalue) || adjpvalue <= 0 || adjpvalue >= 1) {
    stop("adjpvalue must be between 0 and 1")
  }
  if (!is.numeric(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1")
  }
  if (!is.numeric(min_group_size) || min_group_size < 2) {
    stop("min_group_size must be at least 2")
  }

  # Validate skip_tests parameter
  if (!is.null(skip_tests)) {
    if (!is.character(skip_tests)) {
      stop("skip_tests must be a character vector or NULL")
    }
    valid_skip_options <- c("DR", "FC", "CA")
    invalid_options <- setdiff(skip_tests, valid_skip_options)
    if (length(invalid_options) > 0) {
      stop("Invalid skip_tests options: ", paste(invalid_options, collapse = ", "),
           ". Valid options are: ", paste(valid_skip_options, collapse = ", "))
    }
    skip_tests <- toupper(skip_tests)  # Ensure uppercase
  }

  if (!direction %in% c("auto", ">", "<")) {
    stop("direction must be 'auto', '>', or '<'")
  }

  if (!is.numeric(top_n) || top_n < 1 || top_n != round(top_n)) {
    stop("top_n must be a positive integer")
  }

  # Validate minimum_auc parameter
  if (!is.null(minimum_auc)) {
    if (!is.numeric(minimum_auc) || length(minimum_auc) != 1) {
      stop("minimum_auc must be a single numeric value or NULL")
    }
    if (minimum_auc < 0 || minimum_auc > 1) {
      stop("minimum_auc must be between 0 and 1")
    }
    if (minimum_auc < 0.5) {
      warning("minimum_auc < 0.5 may include features with poor discriminative ability")
    }
  }

  if (!is.null(arrangeLevels) && !is.character(arrangeLevels)) {
    stop("arrangeLevels must be a character vector or NULL")
  }

  # ========== INITIALIZE RESULTS ==========
  auroc_results <- list(
    FunctionOrigin              = "perform_AUROC",
    timestamp                   = Sys.time(),
    parameters                  = list(
      arrangeLevels             = arrangeLevels, VIPmin = VIPmin, fcUP = fcUP,
      fcDown                    = fcDown, adjpvalue = adjpvalue, direction = direction,
      top_n                     = top_n, minimum_auc = minimum_auc, confidence_level = confidence_level,
      min_group_size            = min_group_size, skip_tests = skip_tests
    ),
    plots                       = list(),
    plots_identifiedMetabolites = list(),
    summary                     = list()
  )

  # ========== PROCESS GROUP DATA ==========
  # Identify QC and non-QC samples (vectorized)
  qc_groups      <- c("SQC", "EQC", "QC")
  qc_indices     <- data_PP$Metadata$Group %in% qc_groups
  non_qc_indices <- !qc_indices

  if (sum(non_qc_indices) == 0) {
    stop("No non-QC samples found in data")
  }

  groups <- data_PP$Metadata$Group[non_qc_indices]

  # Handle group arrangement
  if (!is.null(arrangeLevels)) {
    unique_groups_data <- unique(groups)
    missing_groups     <- setdiff(arrangeLevels, unique_groups_data)
    extra_groups       <- setdiff(unique_groups_data, arrangeLevels)

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

  # Check minimum group sizes (vectorized)
  group_counts <- table(groups)
  small_groups <- names(group_counts)[group_counts < min_group_size]

  if (length(small_groups) > 0) {
    warning("Groups with fewer than ", min_group_size, " samples: ",
            paste(small_groups, collapse = ", "))
  }

  valid_groups  <- names(group_counts)[group_counts >= min_group_size]
  unique_groups <- intersect(unique_groups, valid_groups)

  if (length(unique_groups) < 2) {
    stop("Need at least 2 groups with minimum ", min_group_size, " samples each")
  }

  # ========== GENERATE GROUP COMBINATIONS ==========
  group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)

  if (length(group_combinations) == 0) {
    warning("No valid group combinations found for analysis.")
    return(auroc_results)
  }

  # Initialize summary statistics
  processed_comparisons   <- 0
  failed_comparisons      <- character(0)
  total_features_analyzed <- 0

  # ========== PROCESS EACH GROUP COMBINATION ==========
  for (i in seq_along(group_combinations)) {

    group_pair <- group_combinations[[i]]

    tryCatch({
      # Set up group comparison (reverse for ROC: control vs case)
      group1     <- group_pair[2]  # case
      group2     <- group_pair[1]  # control
      group_name <- paste(group1, "vs.", group2, sep = " ")

      message("\nProcessing comparison: ", group_name)

      # Check for required data (with skip_tests logic)
      vip_key <- paste0("data_VIPScores_", group_name)
      skip_DR <- "DR" %in% skip_tests

      if (!skip_DR && is.null(data_DR[[vip_key]])) {
        message("Skipping '", group_name, "' - missing VIP data (no OPLS-DA model). ",
                "Consider using skip_tests = 'DR' to bypass this requirement.")
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      # ========== EXTRACT AND MERGE DATA ==========
      skip_DR <- "DR" %in% skip_tests
      skip_FC <- "FC" %in% skip_tests
      skip_CA <- "CA" %in% skip_tests

      # Initialize merged_data with Feature column
      merged_data <- data.frame(Feature = character(), stringsAsFactors = FALSE)

      # Add VIP data if not skipped
      if (!skip_DR) {
        vip_data <- data_DR[[vip_key]]
        if (!is.null(vip_data) && nrow(vip_data) > 0) {
          vip_data <- vip_data[order(vip_data$Feature), ]
          if (nrow(merged_data) == 0) {
            merged_data <- vip_data
          } else {
            merged_data <- dplyr::full_join(merged_data, vip_data, by = "Feature")
          }
        }
      }

      # Add FC data if not skipped
      if (!skip_FC) {
        fc_data <- data_FC[[paste0("data_combined_", group_name)]]
        if (!is.null(fc_data) && nrow(fc_data) > 0) {
          fc_data$Feature <- NULL
          fc_data <- cbind(Feature = rownames(fc_data), fc_data)
          rownames(fc_data) <- NULL
          fc_data <- fc_data[order(fc_data$Feature), ]
          if (nrow(merged_data) == 0) {
            merged_data <- fc_data
          } else {
            merged_data <- dplyr::full_join(merged_data, fc_data, by = "Feature")
          }
        }
      }

      # Add CA data if not skipped
      if (!skip_CA) {
        ca_results <- data_CA$results
        if (!is.null(ca_results) && nrow(ca_results) > 0) {
          ca_results <- cbind(Feature = rownames(ca_results), ca_results)
          rownames(ca_results) <- NULL
          ca_results <- ca_results[order(ca_results$Feature), ]
          if (nrow(merged_data) == 0) {
            merged_data <- ca_results
          } else {
            merged_data <- dplyr::full_join(merged_data, ca_results, by = "Feature")
          }
        }
      }

      # If all tests are skipped, get all features from data matrix
      if (skip_DR && skip_FC && skip_CA) {
        message("All tests skipped - performing AUROC on all available features for ", group_name)
        data_matrix_temp <- if (merge_replicates) {
          data_PP$data_scaledNONPLS_merged[non_qc_indices, , drop = FALSE]
        } else {
          data_PP$data_scaledNONPLS_varFiltered[non_qc_indices, , drop = FALSE]
        }
        merged_data <- data.frame(Feature = colnames(data_matrix_temp), stringsAsFactors = FALSE)
      }

      auroc_results[[paste0("merged_", group_name)]] <- merged_data

      if (nrow(merged_data) == 0) {
        message("Skipping '", group_name, "' - no data to merge")
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      # ========== FILTER FEATURES ==========
      # Build filter condition based on skip_tests
      filter_conditions <- rep(TRUE, nrow(merged_data))

      if (!skip_DR && "VIP" %in% colnames(merged_data)) {
        filter_conditions <- filter_conditions & (merged_data$VIP > VIPmin)
      }

      if (!skip_FC && "fold_change" %in% colnames(merged_data)) {
        filter_conditions <- filter_conditions &
          ((round(merged_data$fold_change, 2) >= fcUP) | (round(merged_data$fold_change, 2) <= fcDown))
      }

      if (!skip_CA && "omnibus_p_value" %in% colnames(merged_data)) {
        filter_conditions <- filter_conditions & (merged_data$omnibus_p_value < adjpvalue)
      }

      filtered_data <- merged_data[filter_conditions, ]

      auroc_results[[paste0("filtered_", group_name)]] <- filtered_data

      if (nrow(filtered_data) == 0) {
        skip_msg <- if (length(skip_tests) > 0) {
          paste0(" (skipped tests: ", paste(skip_tests, collapse = ", "), ")")
        } else {
          ""
        }
        message("No features passed filtering for ", group_name, skip_msg)
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      # ========== PREPARE FEATURE DATA ==========
      # Get appropriate data matrix
      merge_replicates <- isTRUE(data_PP$Parameters$merge_replicates)

      if (merge_replicates) {
        if (is.null(data_PP$data_scaledNONPLS_merged)) {
          stop("merge_replicates is TRUE but data_scaledNONPLS_merged is not available")
        }
        data_matrix <- data_PP$data_scaledNONPLS_merged
      } else {
        if (is.null(data_PP$data_scaledNONPLS_varFiltered)) {
          stop("merge_replicates is FALSE but data_scaledNONPLS_varFiltered is not available")
        }
        data_matrix <- data_PP$data_scaledNONPLS_varFiltered
      }

      data_matrix        <- data_matrix[non_qc_indices, , drop = FALSE]
      available_features <- intersect(filtered_data$Feature, colnames(data_matrix))

      if (length(available_features) == 0) {
        stop("No requested features found in data matrix")
      }

      if (length(available_features) < nrow(filtered_data)) {
        missing_features <- setdiff(filtered_data$Feature, available_features)
        warning("Some features not found in data matrix: ",
                paste(missing_features, collapse = ", "))
      }

      feature_data <- data_matrix[, available_features, drop = FALSE]
      feature_data <- cbind(Groups = groups, feature_data)

      auroc_results[[paste0("samples.x.features_", group_name)]] <- feature_data

      # ========== PERFORM AUROC ANALYSIS ==========
      feature_names <- setdiff(colnames(feature_data), "Groups")

      # Pre-allocate results data frame
      n_features <- length(feature_names)
      auroc_df <- data.frame(
        Feature          = character(n_features),
        Group_Comparison = character(n_features),
        AUROC            = numeric(n_features),
        CI_Lower         = numeric(n_features),
        CI_Upper         = numeric(n_features),
        stringsAsFactors = FALSE
      )

      valid_idx <- 0

      for (j in seq_along(feature_names)) {
        feature_name <- feature_names[j]
        tryCatch({
          roc_obj <- pROC::roc(
            response  = feature_data$Groups,
            predictor = feature_data[[feature_name]],
            levels    = c(group1, group2),
            direction = direction,
            quiet     = TRUE
          )

          auc_val <- as.numeric(pROC::auc(roc_obj))
          ci_vals <- as.numeric(pROC::ci(roc_obj, conf.level = confidence_level))

          valid_idx                            <- valid_idx + 1
          auroc_df$Feature[valid_idx]          <- feature_name
          auroc_df$Group_Comparison[valid_idx] <- group_name
          auroc_df$AUROC[valid_idx]            <- auc_val
          auroc_df$CI_Lower[valid_idx]         <- ci_vals[1]
          auroc_df$CI_Upper[valid_idx]         <- ci_vals[3]

        }, error = function(e) {
          warning("Failed to compute ROC for feature ", feature_name, ": ", e$message)
        })
      }

      # Remove unused rows
      auroc_df <- auroc_df[1:valid_idx, ]
      auroc_df <- auroc_df[order(-auroc_df$AUROC), ]

      if (nrow(auroc_df) == 0) {
        message("No valid AUROC results for ", group_name)
        failed_comparisons <- c(failed_comparisons, group_name)
        next
      }

      auroc_results[[paste0("auroc_results_", group_name)]] <- auroc_df

      # # ========== DETERMINE FEATURES TO PLOT ==========
      # # Use minimum_auc if specified, otherwise use top_n
      # if (!is.null(minimum_auc)) {
      #   plot_features   <- auroc_df[round(auroc_df$AUROC, 2) >= minimum_auc, ]
      #   n_plot_features <- nrow(plot_features)
      #
      #   if (n_plot_features == 0) {
      #     message("No features meet minimum_auc threshold (", minimum_auc, ") for ", group_name)
      #     failed_comparisons <- c(failed_comparisons, group_name)
      #     next
      #   }
      #
      #   message("Plotting ", n_plot_features, " feature(s) with AUC >= ", minimum_auc)
      # } else {
      #   plot_features   <- head(auroc_df, min(top_n, nrow(auroc_df)))
      #   n_plot_features <- nrow(plot_features)
      # }
      #
      # # ========== CREATE MAIN ROC PLOT ==========
      # if (n_plot_features > 0) {
      #   # Generate colors
      #   colors <- if (n_plot_features <= 8) {
      #     RColorBrewer::brewer.pal(max(3, n_plot_features), "Dark2")[1:n_plot_features]
      #   } else {
      #     grDevices::rainbow(n_plot_features)
      #   }
      #
      #   gg_roc <- ggplot2::ggplot()
      #
      #   # Prepare legend data (vectorized)
      #   legend_data <- data.frame(
      #     Feature          = sprintf("%s (AUC = %.2f) (95%% CI: %.2f-%.2f)",  # Changed %.3f to %.2f
      #                                plot_features$Feature,
      #                                round(plot_features$AUROC, 2),
      #                                round(plot_features$CI_Lower, 2),
      #                                round(plot_features$CI_Upper, 2)),
      #     AUC              = plot_features$AUROC,
      #     stringsAsFactors = FALSE
      #   )
      #   legend_data <- legend_data[order(-legend_data$AUC), ]
      #
      #   # Add ROC curves
      #   for (j in seq_len(n_plot_features)) {
      #     feature <- plot_features$Feature[j]
      #
      #     roc_obj <- pROC::roc(
      #       response  = feature_data$Groups,
      #       predictor = feature_data[[feature]],
      #       levels    = c(group1, group2),
      #       direction = direction,
      #       quiet     = TRUE
      #     )
      #
      #     roc_df <- data.frame(
      #       x       = 1 - roc_obj$specificities,
      #       y       = roc_obj$sensitivities,
      #       Feature = legend_data$Feature[j]
      #     )
      #     roc_df <- roc_df[order(roc_df$x, roc_df$y), ]
      #
      #     gg_roc <- gg_roc +
      #       ggplot2::geom_line(
      #         data           = roc_df,
      #         ggplot2::aes(x = x, y = y, color = Feature),
      #         size           = 1.2
      #       )
      #   }
      #
      #   # Finalize plot
      #   plot_title <- if (!is.null(minimum_auc)) {
      #     sprintf("ROC Curves (AUC >= %.2f) - %s", minimum_auc, group_name)
      #   } else {
      #     sprintf("Top %d Features' ROC Curves - %s", min(top_n, nrow(auroc_df)), group_name)
      #   }
      #
      #   main_plot <- gg_roc +
      #     ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
      #                          color = "gray50", size = 0.8) +
      #     ggplot2::scale_color_manual(values = colors, breaks = legend_data$Feature) +
      #     ggplot2::labs(
      #       title = plot_title,
      #       x = "1 - Specificity (False Positive Rate)",
      #       y = "Sensitivity (True Positive Rate)",
      #       color = "Feature (AUC, 95% CI)"
      #     ) +
      #     ggplot2::xlim(0, 1) +
      #     ggplot2::ylim(0, 1) +
      #     ggplot2::theme_minimal() +
      #     ggplot2::theme(
      #       legend.position = c(0.8, 0.2),
      #       plot.title       = ggplot2::element_text(size = 12, hjust = 0.5),
      #       legend.text      = ggplot2::element_text(size = 8),
      #       legend.key.width = ggplot2::unit(1.5, "lines")
      #     )
      #
      #   auroc_results$plots[[group_name]] <- main_plot
      #   print(main_plot)
      # }

      # ========== DETERMINE FEATURES TO PLOT ==========
      # Use minimum_auc if specified, otherwise use top_n
      if (!is.null(minimum_auc)) {
        plot_features   <- auroc_df[round(auroc_df$AUROC, 2) >= minimum_auc, ]
        n_plot_features <- nrow(plot_features)
        if (n_plot_features == 0) {
          message("No features meet minimum_auc threshold (", minimum_auc, ") for ", group_name)
          failed_comparisons <- c(failed_comparisons, group_name)
          next
        }
        message("Plotting ", n_plot_features, " feature(s) with AUC >= ", minimum_auc)
      } else {
        plot_features   <- head(auroc_df, min(top_n, nrow(auroc_df)))
        n_plot_features <- nrow(plot_features)
      }

      # ========== CREATE MAIN ROC PLOT ==========
      if (n_plot_features > 0) {
        # Generate colors
        colors <- if (n_plot_features <= 8) {
          RColorBrewer::brewer.pal(max(3, n_plot_features), "Dark2")[1:n_plot_features]
        } else {
          grDevices::rainbow(n_plot_features)
        }

        gg_roc <- ggplot2::ggplot()

        # Helper function to determine appropriate decimal places for upper CI
        get_upper_ci_decimals <- function(ci_upper) {
          if (ci_upper == 1) {
            return(2)  # If truly 1, use 2 decimals
          }

          # Test increasing decimal places until rounded value is < 1
          for (decimals in 2:10) {
            rounded_val <- round(ci_upper, decimals)
            if (rounded_val < 1) {
              return(decimals)
            }
          }
          return(10)  # Maximum fallback
        }

        # Prepare legend data with adaptive decimal precision
        legend_labels <- character(nrow(plot_features))

        for (i in seq_len(nrow(plot_features))) {
          ci_upper <- plot_features$CI_Upper[i]
          upper_decimals <- get_upper_ci_decimals(ci_upper)

          legend_labels[i] <- sprintf(
            "%s (AUC = %.2f) (95%% CI: %.2f-%.*f)",
            plot_features$Feature[i],
            round(plot_features$AUROC[i], 2),
            round(plot_features$CI_Lower[i], 2),
            upper_decimals,
            round(ci_upper, upper_decimals)
          )
        }

        legend_data <- data.frame(
          Feature          = legend_labels,
          AUC              = plot_features$AUROC,
          stringsAsFactors = FALSE
        )
        legend_data <- legend_data[order(-legend_data$AUC), ]

        # Add ROC curves
        for (j in seq_len(n_plot_features)) {
          feature <- plot_features$Feature[j]
          roc_obj <- pROC::roc(
            response  = feature_data$Groups,
            predictor = feature_data[[feature]],
            levels    = c(group1, group2),
            direction = direction,
            quiet     = TRUE
          )
          roc_df <- data.frame(
            x       = 1 - roc_obj$specificities,
            y       = roc_obj$sensitivities,
            Feature = legend_data$Feature[j]
          )
          roc_df <- roc_df[order(roc_df$x, roc_df$y), ]
          gg_roc <- gg_roc +
            ggplot2::geom_line(
              data           = roc_df,
              ggplot2::aes(x = x, y = y, color = Feature),
              size           = 1.2
            )
        }

        # Finalize plot
        plot_title <- if (!is.null(minimum_auc)) {
          sprintf("ROC Curves (AUC >= %.2f) - %s", minimum_auc, group_name)
        } else {
          sprintf("Top %d Features' ROC Curves - %s", min(top_n, nrow(auroc_df)), group_name)
        }

        main_plot <- gg_roc +
          ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                               color = "gray50", size = 0.8) +
          ggplot2::scale_color_manual(values = colors, breaks = legend_data$Feature) +
          ggplot2::labs(
            title = plot_title,
            x = "1 - Specificity (False Positive Rate)",
            y = "Sensitivity (True Positive Rate)",
            color = "Feature (AUC, 95% CI)"
          ) +
          ggplot2::xlim(0, 1) +
          ggplot2::ylim(0, 1) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            legend.position = c(0.8, 0.2),
            plot.title       = ggplot2::element_text(size = 12, hjust = 0.5),
            legend.text      = ggplot2::element_text(size = 8),
            legend.key.width = ggplot2::unit(1.5, "lines")
          )

        auroc_results$plots[[group_name]] <- main_plot
        print(main_plot)
      }

      # ========== CREATE INDIVIDUAL METABOLITE PLOTS ==========
      if (!is.null(plot_iden_met)) {
        individual_plots <- list()

        # Get full merged data (without filtering)
        vip_data_full   <- data_DR[[vip_key]]
        fc_data_full    <- data_FC[[paste0("data_combined_", group_name)]]
        ca_results_full <- data_CA$results

        fc_data_full$Feature <- NULL

        # Vectorized data preparation
        vip_data_full             <- vip_data_full[order(vip_data_full$Feature), ]
        fc_data_full              <- cbind(Feature = rownames(fc_data_full), fc_data_full)
        rownames(fc_data_full)    <- NULL
        fc_data_full              <- fc_data_full[order(fc_data_full$Feature), ]

        ca_results_full           <- cbind(Feature = rownames(ca_results_full), ca_results_full)
        rownames(ca_results_full) <- NULL
        ca_results_full           <- ca_results_full[order(ca_results_full$Feature), ]

        merged_data_full          <- vip_data_full %>%
          dplyr::full_join(fc_data_full, by = "Feature") %>%
          dplyr::full_join(ca_results_full, by = "Feature")

        # Vectorized filtering for complete cases
        complete_mask         <- !is.na(merged_data_full$VIP) &
          !is.na(merged_data_full$fold_change) &
          !is.na(merged_data_full$omnibus_adj_p_value)
        merged_data_full      <- merged_data_full[complete_mask, ]

        available_metabolites <- intersect(plot_iden_met, merged_data_full$Feature)

        if (length(available_metabolites) > 0) {
          # Prepare feature data for individual plots
          data_matrix_full <- if (merge_replicates) {
            data_PP$data_scaledNONPLS_merged[non_qc_indices, , drop = FALSE]
          } else {
            data_PP$data_scaledNONPLS_varFiltered[non_qc_indices, , drop = FALSE]
          }

          available_in_matrix <- intersect(available_metabolites, colnames(data_matrix_full))

          if (length(available_in_matrix) > 0) {
            feature_data_ind <- data_matrix_full[, available_in_matrix, drop = FALSE]
            feature_data_ind <- cbind(Groups = groups, feature_data_ind)

            for (metabolite in available_in_matrix) {
              tryCatch({
                roc_obj <- pROC::roc(
                  response  = feature_data_ind$Groups,
                  predictor = feature_data_ind[[metabolite]],
                  levels    = c(group1, group2),
                  direction = direction,
                  quiet     = TRUE
                )

                auc_val <- as.numeric(pROC::auc(roc_obj))
                ci_vals <- as.numeric(pROC::ci(roc_obj, conf.level = confidence_level))

                roc_df <- data.frame(
                  x     = 1 - roc_obj$specificities,
                  y     = roc_obj$sensitivities,
                  Label = sprintf("%s (AUC = %.2f, 95%% CI: %.2f-%.2f)",
                                  metabolite, round(auc_val, 2), round(ci_vals[1], 2), round(ci_vals[3], 2))
                )
                roc_df <- roc_df[order(roc_df$x, roc_df$y), ]

                individual_plot <- ggplot2::ggplot() +
                  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                                       color = "gray50", size = 0.8) +
                  ggplot2::geom_line(
                    data           = roc_df,
                    ggplot2::aes(x = x, y = y, color = Label),
                    size           = 1.2
                  ) +
                  ggplot2::labs(
                    title = sprintf("ROC Curve - %s", group_name),
                    x     = "1 - Specificity (False Positive Rate)",
                    y     = "Sensitivity (True Positive Rate)",
                    color = "Metabolite (AUC, 95% CI)"
                  ) +
                  ggplot2::xlim(0, 1) +
                  ggplot2::ylim(0, 1) +
                  ggplot2::theme_minimal() +
                  ggplot2::theme(
                    legend.position = "bottom",
                    plot.title      = ggplot2::element_text(size = 12, hjust = 0.5),
                    legend.text     = ggplot2::element_text(size = 10)
                  )

                print(individual_plot)
                individual_plots[[metabolite]] <- individual_plot

              }, error = function(e) {
                warning("Failed to create plot for metabolite ", metabolite, ": ", e$message)
              })
            }
          }
        }

        if (length(individual_plots) > 0) {
          auroc_results$plots_identifiedMetabolites[[group_name]] <- individual_plots
        }
      }

      processed_comparisons   <- processed_comparisons + 1
      total_features_analyzed <- total_features_analyzed + nrow(auroc_df)

    }, error = function(e) {
      warning("Error processing ", group_name, ": ", e$message)
      failed_comparisons <<- c(failed_comparisons, group_name)
    })
  }

  # ========== ADD SUMMARY INFORMATION ==========
  auroc_results$summary <- list(
    total_comparisons       = length(group_combinations),
    processed_comparisons   = processed_comparisons,
    failed_comparisons      = failed_comparisons,
    total_features_analyzed = total_features_analyzed,
    unique_groups           = unique_groups,
    samples_per_group       = as.list(table(groups)),
    minimum_auc_used        = !is.null(minimum_auc),
    skip_tests              = skip_tests
  )

  message("\nAUROC analysis completed: ", processed_comparisons, "/",
          length(group_combinations), " comparisons processed successfully")

  class(auroc_results) <- c("perform_AUROC", "list")
  return(auroc_results)
}

# ========== UTILITY FUNCTIONS ==========

#' Extract AUROC Results
#' @export
extract_auroc_results <- function(
    auroc_results,
    comparison          = NULL,
    min_auc             = 0,
    max_auc             = 1,
    top_n               = NULL,
    include_merged_data = FALSE,
    sort_by             = "auc_desc"
) {

  if (!"perform_AUROC" %in% auroc_results$FunctionOrigin) {
    stop("Input must be results from perform_AUROC function")
  }

  if (!is.numeric(min_auc) || !is.numeric(max_auc) || min_auc < 0 ||
      max_auc > 1 || min_auc >= max_auc) {
    stop("min_auc and max_auc must be numeric values between 0 and 1, with min_auc < max_auc")
  }

  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n) != 1 ||
                          top_n < 1 || top_n != round(top_n))) {
    stop("top_n must be a positive integer or NULL")
  }

  if (!sort_by %in% c("auc_desc", "auc_asc", "feature", "comparison")) {
    stop("sort_by must be one of: 'auc_desc', 'auc_asc', 'feature', 'comparison'")
  }

  result_names <- names(auroc_results)[grepl("^data_results_", names(auroc_results))]

  if (length(result_names) == 0) {
    warning("No AUROC results found")
    return(list(results = data.frame(), summary = list(), parameters = list(),
                merged_data = list()))
  }

  if (!is.null(comparison)) {
    comparison_names  <- paste0("data_results_", comparison)
    valid_comparisons <- intersect(comparison_names, result_names)

    if (length(valid_comparisons) == 0) {
      stop("None of the specified comparisons found in results. Available comparisons: ",
           paste(gsub("^data_results_", "", result_names), collapse = ", "))
    }

    result_names <- valid_comparisons

    invalid_comparisons <- setdiff(comparison_names, result_names)
    if (length(invalid_comparisons) > 0) {
      warning("Some specified comparisons not found: ",
              paste(gsub("^data_results_", "", invalid_comparisons), collapse = ", "))
    }
  }

  all_results <- data.frame()
  merged_data_list <- list()

  for (result_name in result_names) {
    comparison_name <- gsub("^data_results_", "", result_name)
    result_df       <- auroc_results[[result_name]]

    if (!is.null(result_df) && nrow(result_df) > 0) {
      all_results <- rbind(all_results, result_df)

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
    return(list(results = data.frame(), summary = list(), parameters = list(),
                merged_data = list()))
  }

  filtered_results <- all_results[round(all_results$AUROC, 2) >= min_auc &
                                    round(all_results$AUROC, 2) <= max_auc, ]

  if (nrow(filtered_results) == 0) {
    warning("No results pass the AUC filtering criteria (min_auc: ", min_auc,
            ", max_auc: ", max_auc, ")")
    return(list(results = data.frame(), summary = list(), parameters = list(),
                merged_data = merged_data_list))
  }

  filtered_results$CI_Width <- filtered_results$CI_Upper - filtered_results$CI_Lower
  filtered_results$Performance_Category <- cut(
    filtered_results$AUROC,
    breaks         = c(0, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels         = c("Fail", "Poor", "Fair", "Good", "Excellent"),
    include.lowest = TRUE, right = FALSE
  )

  filtered_results <- switch(sort_by,
                             "auc_desc"   = filtered_results[order(-filtered_results$AUROC), ],
                             "auc_asc"    = filtered_results[order(filtered_results$AUROC), ],
                             "feature"    = filtered_results[order(filtered_results$Feature), ],
                             "comparison" = filtered_results[order(filtered_results$Group_Comparison,
                                                                   -filtered_results$AUROC), ]
  )

  if (!is.null(top_n)) {
    if (sort_by %in% c("feature", "comparison")) {
      filtered_results <- do.call(rbind, lapply(
        split(filtered_results, filtered_results$Group_Comparison),
        function(df) {
          df_sorted <- df[order(-df$AUROC), ]
          head(df_sorted, min(top_n, nrow(df_sorted)))
        }
      ))

      filtered_results <- switch(sort_by,
                                 "feature"    = filtered_results[order(filtered_results$Feature), ],
                                 "comparison" = filtered_results[order(filtered_results$Group_Comparison,
                                                                       -filtered_results$AUROC), ]
      )
    } else {
      filtered_results <- do.call(rbind, lapply(
        split(filtered_results, filtered_results$Group_Comparison),
        function(df) head(df, min(top_n, nrow(df)))
      ))
    }
  }

  rownames(filtered_results) <- NULL

  summary_stats <- list(
    total_features           = nrow(filtered_results),
    unique_comparisons       = length(unique(filtered_results$Group_Comparison)),
    comparisons              = unique(filtered_results$Group_Comparison),
    auc_statistics           = list(
      mean                   = mean(filtered_results$AUROC, na.rm = TRUE),
      median                 = median(filtered_results$AUROC, na.rm = TRUE),
      min                    = min(filtered_results$AUROC, na.rm = TRUE),
      max                    = max(filtered_results$AUROC, na.rm = TRUE),
      sd                     = sd(filtered_results$AUROC, na.rm = TRUE)
    ),
    performance_distribution = table(filtered_results$Performance_Category),
    ci_width_statistics      = list(
      mean                   = mean(filtered_results$CI_Width, na.rm = TRUE),
      median                 = median(filtered_results$CI_Width, na.rm = TRUE),
      min                    = min(filtered_results$CI_Width, na.rm = TRUE),
      max                    = max(filtered_results$CI_Width, na.rm = TRUE)
    )
  )

  extraction_params <- list(
    comparison           = comparison,
    min_auc              = min_auc,
    max_auc              = max_auc,
    top_n                = top_n,
    include_merged_data  = include_merged_data,
    sort_by              = sort_by,
    extraction_timestamp = Sys.time()
  )

  return(list(
    results     = filtered_results,
    summary     = summary_stats,
    parameters  = extraction_params,
    merged_data = merged_data_list
  ))
}

#' Get AUROC summary statistics
#' @param auroc_results List returned by perform_AUROC
#' @return Data frame with summary statistics
#' @export
get_auroc_summary <- function(auroc_results) {

  if (!"perform_AUROC" %in% auroc_results$FunctionOrigin) {
    stop("Input must be results from perform_AUROC function")
  }

  result_names <- names(auroc_results)[grepl("^data_results_", names(auroc_results))]

  if (length(result_names) == 0) {
    warning("No AUROC results found")
    return(data.frame())
  }

  all_results <- do.call(rbind, lapply(result_names, function(name) {
    auroc_results[[name]]
  }))

  summary_stats <- data.frame(
    Total_Comparisons      = auroc_results$summary$total_comparisons,
    Processed_Comparisons  = auroc_results$summary$processed_comparisons,
    Total_Features         = nrow(all_results),
    Mean_AUC               = mean(all_results$AUROC, na.rm = TRUE),
    Median_AUC             = median(all_results$AUROC, na.rm = TRUE),
    Min_AUC                = min(all_results$AUROC, na.rm = TRUE),
    Max_AUC                = max(all_results$AUROC, na.rm = TRUE),
    Features_AUC_Above_0.7 = sum(round(all_results$AUROC, 2) > 0.7, na.rm = TRUE),
    Features_AUC_Above_0.8 = sum(round(all_results$AUROC, 2) > 0.8, na.rm = TRUE),
    Features_AUC_Above_0.9 = sum(round(all_results$AUROC, 2) > 0.9, na.rm = TRUE)
  )

  return(summary_stats)
}

#' Plot AUROC distribution
#' @param auroc_results List returned by perform_AUROC
#' @param bins Integer. Number of histogram bins
#' @return ggplot2 object
#' @export
plot_auroc_distribution <- function(auroc_results, bins = 20) {

  if (!"perform_AUROC" %in% auroc_results$FunctionOrigin) {
    stop("Input must be results from perform_AUROC function")
  }

  result_names <- names(auroc_results)[grepl("^data_results_", names(auroc_results))]

  if (length(result_names) == 0) {
    warning("No AUROC results found")
    return(ggplot2::ggplot())
  }

  all_aucs <- do.call(c, lapply(result_names, function(name) {
    auroc_results[[name]]$AUROC
  }))

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
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}

# S3 Methods
#' @export
print.perform_AUROC <- function(x, ...) {
  cat("=== AUROC Analysis ===\n")
  cat("Comparisons:     ", x$summary$processed_comparisons, "\n")
  cat("Features Tested: ", x$summary$total_features_analyzed, "\n")

  # Check if we have results to calculate average AUC
  res_names <- names(x)[grepl("^data_results_", names(x))]
  if(length(res_names) > 0) {
    all_aucs <- unlist(lapply(res_names, function(n) x[[n]]$AUROC))
    cat("Mean AUC:        ", round(mean(all_aucs, na.rm=TRUE), 3), "\n")
  }
  invisible(x)
}

#' @export
summary.perform_AUROC <- function(object, ...) {
  stats <- get_auroc_summary(object) # Using the helper function inside your code
  ans <- list(
    stats  = stats,
    failed = object$summary$failed_comparisons
  )
  class(ans) <- "summary.perform_AUROC"
  return(ans)
}

#' @export
print.summary.perform_AUROC <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("AUROC Performance Summary\n")
  cat("---------------------------------------\n")
  if(nrow(x$stats) > 0) {
    print(t(x$stats)) # Transpose for better readability if single row
  } else {
    cat("No significant results found.\n")
  }

  if(length(x$failed) > 0) {
    cat("\nFailed Comparisons:", paste(x$failed, collapse=", "), "\n")
  }
  invisible(x)
}
