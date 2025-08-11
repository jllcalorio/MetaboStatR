#' Perform Area Under the Receiver Operating Characteristic (AUROC)
#'
#' @description
#' This function performs the AUROC analysis using the data from data processing, dimension-reduction,
#' fold change analysis, and comparative analysis. Data will be automatically selected from those lists.
#'
#' @param data_PP List. This list must be a result from the `performPreprocessingPeakData` function.
#' @param data_DR List. This list must be a result from the `performDimensionReduction` function, specifically from OPLS-DA.
#' @param data_FC List. This list must be a result from the `performFoldChange` function.
#' @param data_CA List. This list must be a result from the `performComparativeAnalysis` function.
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order. Suggests to be inputted (control, case1, case2) where case2 is worse than case1 e.g., Severe dengue than Non-severe dengue.
#' @param VIPmin Numeric. Minimum VIP score to consider.
#' @param fcUP Numeric. Minimum fold change to consider feature as up-regulated.
#' @param fcDown Numeric. Maximum fold change to consider feature as down-regulated.
#' @param adjpvalue Numeric. Defaults to 0.05. Defines the adjusted p-value threshold.
#' @param direction Defines the direction in AUROC.
#'   \itemize{
#'     \item "auto": (This is not applied.)
#'     \item ">": median of control < median of case
#'     \item "<": median of control > median of case
#'     \item "auto2": Dynamically change the direction. If `median of control < median of case`, then `>` will be applied, ottherwise, `<`.
#'     }
#'     Defaults to ">".
#' @param top_n Numeric. Defines the features with the highest AUC. Defaults to plotting the top 5.
#' @param plot_iden_met Boolean. If `TRUE`, plots the identified metabolites, one-by-one. Defaults to `NULL`.
#'
#' @returns Returns a list of results and plots if requested.
#' @export
#'
#' @examples
#' \dontrun{
#'   data_PP = results_from_performPreprocessingPeakData_function,
#'   data_DR = results_from_performDimensionReduction_function,
#'   data_FC = results_from_performFoldChange_function,
#'   data_CA = results_from_performComparativeAnalysis_function
#' }
performAUROC <- function(
    data_PP,
    data_DR,
    data_FC,
    data_CA,
    arrangeLevels  = NULL,
    VIPmin         = 1,
    fcUP           = 2,
    fcDown         = 0.5,
    adjpvalue      = 0.05,
    direction      = "auto2",
    top_n          = 5,
    plot_iden_met = NULL
) {

  auroc_results                <- list() # List to store all relevant results
  auroc_results$FunctionOrigin <- "performAUROC"

  non_qc_indices <- data_PP$Metadata$Groups != "QC"
  groups         <- data_PP$Metadata$Groups[non_qc_indices]

  # Check if arrangeLevels have the same actual values of groups
  if (!is.null(arrangeLevels)) {
    if (length(setdiff(arrangeLevels, unique(groups))) > 0) {
      stop("Invalid arrangeLevels: check for typos or missing groups.")
    } else {
      groups        <- factor(groups, levels = arrangeLevels)
      unique_groups <- arrangeLevels
    }
  } else {
    unique_groups <- unique(groups)
  }

  # unique_groups <- ifelse(is.null(arrangeLevels), unique(groups), levels(groups))

  group_combinations <- combn(unique_groups, 2, simplify = FALSE)

  for (group_pair in group_combinations) {

    group2     <- group_pair[1] # The group here is originally group1, but since ROC requires (control, case), thus
    group1     <- group_pair[2] # The group here is originally group2, but since ROC requires (control, case), thus
    group_name <- paste(group1, "vs.", group2, sep = " ")

    if (is.null(data_DR[[paste0("data_VIPScores_", group_name)]])) {
      message(paste0("Skipping '", group_name, "' due to missing VIP data. This is due to no model being created in OPLS-DA."))
      next
    }

    # Extract matching datasets from DR, FC, and CA
    vip_data   <- data_DR[[paste0("data_VIPScores_", group_name)]] %>%
      dplyr::arrange(Feature)
    fc_data    <- data_FC[[paste0("data_combined_", group_name)]] %>%
      cbind(Feature = rownames(.), .) %>%
      `rownames<-`(NULL) %>%
      dplyr::arrange(Feature)
    ca_results <- data_CA$results %>%
      cbind(Feature = rownames(.), .) %>%
      `rownames<-`(NULL) %>%
      dplyr::arrange(Feature)

    if (is.null(vip_data) || is.null(fc_data) || is.null(ca_results)) {
      message("Skipping ", group_name, " due to missing data.")
      next
    }

    # Merge the datasets
    merged_data <- dplyr::full_join(vip_data, fc_data, by = "Feature") %>%
      dplyr::full_join(ca_results, by = "Feature") %>%
      dplyr::filter(VIP > VIPmin,
                    fold_change >= fcUP | fold_change <= fcDown,
                    `adj. p-value` < adjpvalue)

    if (nrow(merged_data) == 0) {
      message(paste0("No features passed filtering for the ", group_name, " comparison."))
      next
    }

    # auroc_results$data_merged[[group_name]] <- merged_data
    auroc_results[[paste0("data_Merged_", group_name)]] <- merged_data

    filtered_features <- merged_data$Feature

    new_data <- data_PP$data_scaledOPLSDA[non_qc_indices, ] %>%
      dplyr::select(any_of(filtered_features)) %>% # Filter features that passed the set thresholds
      cbind(Groups = groups, .)

    # auroc_results$data_filtered[[group_name]] <- new_data
    auroc_results[[paste0("data_Filtered_", group_name)]] <- new_data

    # Initialize per-group result container
    auroc_results_df <- data.frame(
      Feature          = character(),
      Group_Comparison = character(),
      AUROC            = numeric(),
      CI_Lower         = numeric(),
      CI_Upper         = numeric(),
      stringsAsFactors = FALSE
    )

    for (feature_name in filtered_features) {
      if (!(feature_name %in% colnames(new_data))) next

      # Manually specify direction
      if (direction == "auto2") {
        direction = ifelse(
          median(new_data[[feature_name]][new_data$Groups == group1]) < median(new_data[[feature_name]][new_data$Groups == group2]),
          ">", # if median of control < median of case
          "<"  # if otherwise
        )
      }

      # Compute ROC curve
      roc_obj <- pROC::roc(
        response  = new_data$Groups,
        predictor = new_data[[feature_name]],
        levels    = c(group1, group2),
        direction = direction
      )

      auroc_results_df <- rbind(auroc_results_df, data.frame(
        Feature          = feature_name,
        Group_Comparison = group_name,
        AUROC            = as.numeric(pROC::auc(roc_obj)),
        CI_Lower         = as.numeric(pROC::ci(roc_obj)[1]),
        CI_Upper         = as.numeric(pROC::ci(roc_obj)[3])
      )) %>% as.data.frame() %>% dplyr::arrange(desc(AUROC)) # Sort: Highest AUC first
    }


    ########################
    # PLOT
    ########################

    # Store top N features
    top_features <- head(auroc_results_df, min(top_n, nrow(auroc_results_df)))

    # Generate ROC objects for each feature
    roc_list <- lapply(seq_len(nrow(top_features)), function(i) {
      feature <- top_features$Feature[i]
      pROC::roc(
        response  = new_data$Groups,
        predictor = new_data[[feature]],
        levels    = c(group1, group2),
        direction = direction
      )
    })

    # Create a single ggplot object with all ROC curves
    gg_roc   <- ggplot2::ggplot()
    auc_data <- data.frame()

    for (i in base::seq_along(roc_list)) {
      feature <- top_features$Feature[i]
      roc_obj <- roc_list[[i]]
      auc_val <- top_features$AUROC[i]
      ci_low  <- top_features$CI_Lower[i]
      ci_up   <- top_features$CI_Upper[i]

      # Get ROC curve data as data frame
      roc_df <- data.frame(
        x       = 1 - roc_obj$specificities,   # Correct the X-axis, i.e., plot must go up from left to right, not from bottom-up
        y       = roc_obj$sensitivities,
        Feature = paste0(
          feature,
          " (AUC = ", sprintf("%.3f", auc_val), ")",
          " (95% CI = ", sprintf("%.3f", ci_low), "-", sprintf("%.3f", ci_up), ")"
        )
      ) %>% dplyr::arrange(y)

      gg_roc <- gg_roc +
        geom_line(
          data = roc_df,
          aes(x = x, y = y, color = Feature),
          size = 1.2
        )

      # This is for sorting the legend in AUC plot
      # Create a data frame with feature names and AUC values
      auc_data <- rbind(auc_data, data.frame(
        Feature = paste0(
          feature,
          " (AUC = ", sprintf("%.3f", auc_val), ")",
          " (95% CI = ", sprintf("%.3f", ci_low), "-", sprintf("%.3f", ci_up), ")"
        ),
        AUC = auc_val
      ))
    }

    auc_data       <- auc_data %>% dplyr::arrange(desc(AUC)) # Sort by AUC in ascending order
    feature_levels <- auc_data$Feature # Use the sorted Feature names

    # Define colors for plotting
    # colors <- c("blue", "red", "green", "purple", "orange")
    # Dynamically generate distinct colors based on top_n
    # top_n <- nrow(top_features)
    if (top_n <= 8) {
      colors <- RColorBrewer::brewer.pal(top_n, "Dark2")
    } else {
      colors <- grDevices::rainbow(top_n)
    }

    # Finalize plot
    combined_plot <- gg_roc +
      # scale_color_manual(values = colors[seq_len(nrow(top_features))]) +
      scale_color_manual(
        values = colors,
        breaks = feature_levels  # Ensures the legend follows the sorted order
      ) +
      labs(
        title = paste0("Top ", top_n, " features' ROC Curves - ", group_name),
        x     = "1 - Specificity",
        y     = "Sensitivity",
        color = "Feature (AUC)"
      ) +
      theme_minimal() +
      # theme(legend.position = "bottom")
      theme(legend.position = c(0.8, 0.2))  # Adjust the coordinates

    ########################
    ########################

    # Store results specific to the group pair
    auroc_results$results_auroc[[group_name]] <- roc_obj
    # auroc_results$results[[group_name]]       <- auroc_results_df
    auroc_results[[paste0("data_results_", group_name)]]       <- auroc_results_df
    # Save to results
    auroc_results$plots[[group_name]]         <- combined_plot

    print(combined_plot)

    # Plot identified metabolites one-by-one
    if (!is.null(plot_iden_met)) {

      # Merged data regardless if it has VIP, fcUP, fcDown, or adjusted p-value according to set thresholds
      merged_data2 <- dplyr::full_join(vip_data, fc_data, by = "Feature") %>%
        dplyr::full_join(ca_results, by = "Feature")

      merged_data2_metabolites <- merged_data2$Feature # Get metabolite names

      metabolites_to_plot  <- base::intersect(plot_iden_met, merged_data2_metabolites)

      new_data2 <- data_PP$data_scaledOPLSDA[non_qc_indices, ] %>%
        dplyr::select(any_of(metabolites_to_plot)) %>% # Filter metabolites that passed the set thresholds
        cbind(Groups = groups, .)

      if (length(metabolites_to_plot) > 0) {
        auroc_results$plots_identifiedMetabolites[[group_name]] <- list()
        for (metabolite in metabolites_to_plot) {

          # Initialize per-group result container
          auroc_results_df2 <- data.frame(
            Feature          = character(),
            Group_Comparison = character(),
            AUROC            = numeric(),
            CI_Lower         = numeric(),
            CI_Upper         = numeric(),
            stringsAsFactors = FALSE
          )

          # Determine direction
          if (direction == "auto2") {
            direction <- ifelse(
              median(new_data2[[metabolite]][new_data2$Groups == group1]) < median(new_data2[[metabolite]][new_data2$Groups == group2]),
              ">",
              "<"
            )
          }

          # Compute ROC
          roc_obj2 <- pROC::roc(
            response  = new_data2$Groups,
            predictor = new_data2[[metabolite]],
            levels    = c(group1, group2),
            direction = direction
          )

          auc_val <- as.numeric(pROC::auc(roc_obj2))
          ci_vals <- as.numeric(pROC::ci(roc_obj2))
          auc_label <- paste0(
            metabolite,
            " (AUC = ", sprintf("%.3f", auc_val),
            ", 95% CI: ", sprintf("%.3f", ci_vals[1]), "-", sprintf("%.3f", ci_vals[3]), ")"
          )

          auroc_results_df2 <- rbind(auroc_results_df2, data.frame(
            Feature          = metabolite,
            Group_Comparison = group_name,
            AUROC            = auc_val,
            CI_Lower         = ci_vals[1],
            CI_Upper         = ci_vals[3]
          ))

          # Create ROC dataframe
          roc_df <- data.frame(
            x          = 1 - roc_obj2$specificities,
            y          = roc_obj2$sensitivities,
            Label      = auc_label
          ) %>% dplyr::arrange(y)

          # Plot
          now_plot_iden_met <- ggplot() +
            geom_line(
              data = roc_df,
              aes(x = x, y = y, color = Label),
              size = 1.2
            ) +
            labs(
              title = paste("ROC Curves -", group_name),
              x     = "1 - Specificity",
              y     = "Sensitivity",
              color = "Metabolite (AUC, 95% CI)"
            ) +
            theme_minimal() +
            theme(legend.position = "bottom")

          print(now_plot_iden_met)

          auroc_results$plots_identifiedMetabolites[[group_name]][[metabolite]] <- now_plot_iden_met
        }
      }
    }
  }

  return(auroc_results)
}
