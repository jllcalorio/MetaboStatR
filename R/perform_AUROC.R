#' Perform Area Under the Receiver Operating Characteristic (AUROC)
#'
#' @description
#' This function performs the AUROC analysis using the data from data processing, dimension-reduction,
#' fold change analysis, and comparative analysis. Data will be automatically selected from those lists.
#'
#' @param data_PP List. This list must be a result from the `perform_PreprocessingPeakData` function.
#' @param data_DR List. This list must be a result from the `perform_DimensionReduction` function, specifically from OPLS-DA.
#' @param data_FC List. This list must be a result from the `perform_FoldChange` function.
#' @param data_CA List. This list must be a result from the `perform_ComparativeAnalysis` function.
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order. Suggests to be inputted (control, case1, case2) where case2 is worse than case1 e.g., Severe dengue than Non-severe dengue.
#' @param VIPmin Numeric. Minimum VIP score to consider.
#' @param fcUP Numeric. Minimum fold change to consider feature as up-regulated.
#' @param fcDown Numeric. Maximum fold change to consider feature as down-regulated.
#' @param adjpvalue Numeric. Defaults to 0.05. Defines the adjusted p-value threshold.
#' @param direction Defines the direction in AUROC. See more details in `?pROC::roc`. The `>` and `<` are used when resampling or randomizing the data.
#'   \itemize{
#'     \item "auto": Dynamically change the direction of which group is higher.
#'     \item ">": if the predictor values for the control group are higher than the values of the case group will be applied
#'     \item "<": otherwise of ">".
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
#'   data_PP = results_from_perform_PreprocessingPeakData_function,
#'   data_DR = results_from_perform_DimensionReduction_function,
#'   data_FC = results_from_perform_FoldChange_function,
#'   data_CA = results_from_perform_ComparativeAnalysis_function
#' }
perform_AUROC <- function(
    data_PP,
    data_DR,
    data_FC,
    data_CA,
    arrangeLevels  = NULL,
    VIPmin         = 1,
    fcUP           = 2,
    fcDown         = 0.5,
    adjpvalue      = 0.05,
    direction      = c("auto", ">", "<")[1],
    top_n          = 5,
    plot_iden_met  = NULL
) {

  auroc_results                <- base::list() # List to store all relevant results
  auroc_results$FunctionOrigin <- "perform_AUROC"

  if (!(direction %in% c("auto", ">", "<"))) {
    stop("The parameter 'direction' must be either >, <, or auto. Choose '>' or '<' only when 'resampling or randomizing the data'. See more in ?pROC::roc.")
  }

  qc_indices     <- data_PP$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices <- !qc_indices

  # Define groups
  groups         <- data_PP$Metadata$Group[non_qc_indices]

  # Check if arrangeLevels have the same actual values of groups
  if (!base::is.null(arrangeLevels)) {
    if (base::length(base::setdiff(arrangeLevels, base::unique(groups))) > 0) {
      stop("Invalid arrangeLevels: check for typos or missing groups.")
    } else {
      groups        <- base::factor(groups, levels = arrangeLevels)
      unique_groups <- arrangeLevels
    }
  } else {
    unique_groups <- base::unique(groups)
  }

  # unique_groups <- ifelse(is.null(arrangeLevels), unique(groups), levels(groups))

  group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)

  for (group_pair in group_combinations) {

    group2     <- group_pair[1] # The group here is originally group1, but since ROC requires (control, case), thus
    group1     <- group_pair[2] # The group here is originally group2, but since ROC requires (control, case), thus
    group_name <- base::paste(group1, "vs.", group2, sep = " ")

    if (is.null(data_DR[[base::paste0("data_VIPScores_", group_name)]])) {
      message(base::paste0("Skipping '", group_name, "' due to missing VIP data. This is due to no model being created in OPLS-DA."))
      next
    }

    # Extract matching datasets from DR, FC, and CA
    vip_data   <- data_DR[[base::paste0("data_VIPScores_", group_name)]] %>%
      dplyr::arrange(Feature)
    fc_data    <- data_FC[[base::paste0("data_combined_",  group_name)]] %>%
      base::cbind(Feature = base::rownames(.), .) %>%
      base::`rownames<-`(NULL) %>%
      dplyr::arrange(Feature)
    ca_results <- data_CA$results %>%
      base::cbind(Feature = base::rownames(.), .) %>%
      base::`rownames<-`(NULL) %>%
      dplyr::arrange(Feature)

    if (base::is.null(vip_data) || base::is.null(fc_data) || base::is.null(ca_results)) {
      message("Skipping ", group_name, " due to missing data.")
      next
    }

    # Merge the datasets
    merged_data <- dplyr::full_join(vip_data, fc_data, by = "Feature") %>%
      dplyr::full_join(ca_results, by = "Feature") %>%
      dplyr::filter(VIP > VIPmin,
                    fold_change >= fcUP | fold_change <= fcDown,
                    `adj. p-value` < adjpvalue)

    if (base::nrow(merged_data) == 0) {
      message(base::paste0("No features passed filtering for the ", group_name, " comparison."))
      next
    }

    # auroc_results$data_merged[[group_name]] <- merged_data
    auroc_results[[base::paste0("data_Merged_", group_name)]] <- merged_data

    filtered_features <- merged_data$Feature

    new_data <- data_PP$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ] %>%
      dplyr::select(dplyr::any_of(filtered_features)) %>% # Filter features that passed the set thresholds
      base::cbind(Groups = groups, .)

    # auroc_results$data_filtered[[group_name]] <- new_data
    auroc_results[[base::paste0("data_Filtered_", group_name)]] <- new_data

    # Initialize per-group result container
    auroc_results_df <- base::data.frame(
      Feature          = base::character(),
      Group_Comparison = base::character(),
      AUROC            = base::numeric(),
      CI_Lower         = base::numeric(),
      CI_Upper         = base::numeric(),
      stringsAsFactors = FALSE
    )

    for (feature_name in filtered_features) {
      if (!(feature_name %in% base::colnames(new_data))) next

      # # Manually specify direction
      # if (direction == "auto2") {
      #   direction = base::ifelse(
      #     stats::median(new_data[[feature_name]][new_data$Groups == group1]) < stats::median(new_data[[feature_name]][new_data$Groups == group2]),
      #     ">", # if median of control < median of case
      #     "<"  # if otherwise
      #   )
      # }

      # Compute ROC curve
      roc_obj <- pROC::roc(
        response  = new_data$Groups,
        predictor = new_data[[feature_name]],
        levels    = base::c(group1, group2),
        direction = direction
      )

      auroc_results_df <- base::rbind(auroc_results_df, base::data.frame(
        Feature          = feature_name,
        Group_Comparison = group_name,
        AUROC            = as.numeric(pROC::auc(roc_obj)),
        CI_Lower         = as.numeric(pROC::ci(roc_obj)[1]),
        CI_Upper         = as.numeric(pROC::ci(roc_obj)[3])
      )) %>% base::as.data.frame() %>% dplyr::arrange(dplyr::desc(AUROC)) # Sort: Highest AUC first
    }


    ########################
    # PLOT
    ########################

    # Store top N features
    top_features <- utils::head(auroc_results_df, base::min(top_n, base::nrow(auroc_results_df)))

    # Generate ROC objects for each feature
    roc_list <- base::lapply(base::seq_len(base::nrow(top_features)), function(i) {
      feature <- top_features$Feature[i]
      pROC::roc(
        response  = new_data$Groups,
        predictor = new_data[[feature]],
        levels    = base::c(group1, group2),
        direction = direction
      )
    })

    # Create a single ggplot object with all ROC curves
    gg_roc   <- ggplot2::ggplot()
    auc_data <- base::data.frame()

    for (i in base::seq_along(roc_list)) {
      feature <- top_features$Feature[i]
      roc_obj <- roc_list[[i]]
      auc_val <- top_features$AUROC[i]
      ci_low  <- top_features$CI_Lower[i]
      ci_up   <- top_features$CI_Upper[i]

      # Get ROC curve data as data frame
      roc_df <- base::data.frame(
        x       = 1 - roc_obj$specificities,   # Correct the X-axis, i.e., plot must go up from left to right, not from bottom-up
        y       = roc_obj$sensitivities,
        Feature = base::paste0(
          feature,
          " (AUC = ", base::sprintf("%.3f", auc_val), ")",
          " (95% CI = ", base::sprintf("%.3f", ci_low), "-", base::sprintf("%.3f", ci_up), ")"
        )
      ) %>% dplyr::arrange(y)

      gg_roc <- gg_roc +
        ggplot2::geom_line(
          data = roc_df,
          ggplot2::aes(x = x, y = y, color = Feature),
          size = 1.2
        )

      # This is for sorting the legend in AUC plot
      # Create a data frame with feature names and AUC values
      auc_data <- base::rbind(auc_data, data.frame(
        Feature = base::paste0(
          feature,
          " (AUC = ", base::sprintf("%.3f", auc_val), ")",
          " (95% CI = ", base::sprintf("%.3f", ci_low), "-", base::sprintf("%.3f", ci_up), ")"
        ),
        AUC = auc_val
      ))
    }

    auc_data       <- auc_data %>% dplyr::arrange(dplyr::desc(AUC)) # Sort by AUC in ascending order
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
      ggplot2::scale_color_manual(
        values = colors,
        breaks = feature_levels  # Ensures the legend follows the sorted order
      ) +
      ggplot2::labs(
        title = base::paste0("Top ", top_n, " features' ROC Curves - ", group_name),
        x     = "1 - Specificity",
        y     = "Sensitivity",
        color = "Feature (AUC)"
      ) +
      ggplot2::theme_minimal() +
      # ggplot2::theme(legend.position = "bottom")
      ggplot2::theme(legend.position = base::c(0.8, 0.2))  # Adjust the coordinates

    ########################
    ########################

    # Store results specific to the group pair
    auroc_results$results_auroc[[group_name]] <- roc_obj
    # auroc_results$results[[group_name]]       <- auroc_results_df
    auroc_results[[base::paste0("data_results_", group_name)]]       <- auroc_results_df
    # Save to results
    auroc_results$plots[[group_name]]         <- combined_plot

    print(combined_plot)

    # Plot identified metabolites one-by-one
    if (!base::is.null(plot_iden_met)) {

      # Merged data regardless if it has VIP, fcUP, fcDown, or adjusted p-value according to set thresholds
      merged_data2 <- dplyr::full_join(vip_data, fc_data, by = "Feature") %>%
        dplyr::full_join(ca_results, by = "Feature")

      merged_data2_metabolites <- merged_data2$Feature # Get metabolite names

      metabolites_to_plot      <- base::intersect(plot_iden_met, merged_data2_metabolites)

      new_data2 <- data_PP$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ] %>%
        dplyr::select(dplyr::any_of(metabolites_to_plot)) %>% # Filter metabolites that passed the set thresholds
        base::cbind(Groups = groups, .)

      if (base::length(metabolites_to_plot) > 0) {
        auroc_results$plots_identifiedMetabolites[[group_name]] <- base::list()
        for (metabolite in metabolites_to_plot) {

          # Initialize per-group result container
          auroc_results_df2 <- base::data.frame(
            Feature          = base::character(),
            Group_Comparison = base::character(),
            AUROC            = base::numeric(),
            CI_Lower         = base::numeric(),
            CI_Upper         = base::numeric(),
            stringsAsFactors = FALSE
          )

          # # Determine direction
          # if (direction == "auto2") {
          #   direction <- base::ifelse(
          #     stats::median(new_data2[[metabolite]][new_data2$Groups == group1]) < stats::median(new_data2[[metabolite]][new_data2$Groups == group2]),
          #     ">",
          #     "<"
          #   )
          # }

          # Compute ROC
          roc_obj2 <- pROC::roc(
            response  = new_data2$Groups,
            predictor = new_data2[[metabolite]],
            levels    = base::c(group1, group2),
            direction = direction
          )

          auc_val   <- base::as.numeric(pROC::auc(roc_obj2))
          ci_vals   <- base::as.numeric(pROC::ci(roc_obj2))
          auc_label <- base::paste0(
            metabolite,
            " (AUC = ", base::sprintf("%.3f", auc_val),
            ", 95% CI: ", base::sprintf("%.3f", ci_vals[1]), "-", base::sprintf("%.3f", ci_vals[3]), ")"
          )

          auroc_results_df2 <- base::rbind(auroc_results_df2, base::data.frame(
            Feature          = metabolite,
            Group_Comparison = group_name,
            AUROC            = auc_val,
            CI_Lower         = ci_vals[1],
            CI_Upper         = ci_vals[3]
          ))

          # Create ROC dataframe
          roc_df <- base::data.frame(
            x          = 1 - roc_obj2$specificities,
            y          = roc_obj2$sensitivities,
            Label      = auc_label
          ) %>% dplyr::arrange(y)

          # Plot
          now_plot_iden_met <- ggplot2::ggplot() +
            ggplot2::geom_line(
              data = roc_df,
              ggplot2::aes(x = x, y = y, color = Label),
              size = 1.2
            ) +
            ggplot2::labs(
              title = base::paste("ROC Curves -", group_name),
              x     = "1 - Specificity",
              y     = "Sensitivity",
              color = "Metabolite (AUC, 95% CI)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "bottom")

          base::print(now_plot_iden_met)

          auroc_results$plots_identifiedMetabolites[[group_name]][[metabolite]] <- now_plot_iden_met
        }
      }
    }
  }

  return(auroc_results)
}
