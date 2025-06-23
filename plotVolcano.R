# Function for Volcano Plot

plotVolcano <- function(
    PPData,
    FCData,
    CAData,
    arrangeLevels = NULL,        # Vector. A user-input data. Defaults to NULL = unique groups. Suggests to be inputted (control, case1, case2) where case2 is worse than case1 e.g., Severe dengue than Non-severe dengue
    fcUP          = 2,
    fcDown        = 0.5,
    adjpvalue     = 0.05
) {

  volcanoPlotResults                         <- list()
  volcanoPlotResults$Class                   <- "plotVolcano"
  volcanoPlotResults$FoldChangeData          <- FCData
  volcanoPlotResults$ComparativeAnalysisData <- CAData
  volcanoPlotResults$FoldChange_Upper        <- fcUP
  volcanoPlotResults$FoldChange_Down         <- fcDown
  volcanoPlotResults$AdjustedPValue          <- adjpvalue

  non_qc_indices <- PPData$Metadata$Groups != "QC"
  groups         <- PPData$Metadata$Groups[non_qc_indices]

  # Check if arrangeLevels have the same actual values of groups
  if (!is.null(arrangeLevels)) {
    if (length(setdiff(arrangeLevels, unique(groups))) > 0) {
      stop("Invalid arrangeLevels: check for typos or missing groups.")
    } else {
      groups <- factor(groups, levels = arrangeLevels) %>% levels()
    }
  }

  unique_groups <- unique(groups)
  group_combinations <- combn(unique_groups, 2, simplify = FALSE)


  for (group_pair in group_combinations) {

    group1 <- group_pair[1] # The group here is originally group1, but since ROC requires (control, case), thus
    group2 <- group_pair[2] # The group here is originally group2, but since ROC requires (control, case), thus
    group_name <- paste(group1, "vs.", group2, sep = " ")

    ## NOT INCLUDED SINCE RESULTS FROM OPLS-DA DOES NOT MATTER HERE
    ## KEEPING THE SCRIPT JSUT CAUSE

    # if (is.null(data_DR[[paste0("data_VIPScores_", group_name)]])) {
    #   message(paste0("Skipping '", group_name, "' due to missing VIP data. This is due to no model being created in OPLS-DA."))
    #   next
    # }


    # Extract matching datasets from DR, FC, and CA
    # vip_data   <- data_DR[[paste0("data_VIPScores_", group_name)]] %>%
    #   dplyr::arrange(Feature)
    fc_data    <- FCData[[paste0("data_combined_", group_name)]] %>%
      cbind(Feature = rownames(.), .) %>%
      `rownames<-`(NULL) %>%
      dplyr::arrange(Feature)
    ca_results <- CAData$results %>%
      cbind(Feature = rownames(.), .) %>%
      `rownames<-`(NULL) %>%
      dplyr::arrange(Feature)

    # Merge the datasets
    volcanoPlotData <- dplyr::full_join(fc_data, ca_results, by = "Feature") %>%
      # dplyr::full_join(ca_results, by = "Feature") %>%
      # dplyr::filter(fold_change >= fcUP | fold_change <= fcDown,
      #               `adj. p-value` < adjpvalue) %>%
      dplyr::rename(
        c("Fold_Change" = fold_change,
          "Log2_Fold_Change" = log2_fc,
          "P_Value" = `p-value`,
          "Adjusted_P_Value" = `adj. p-value`)
      ) %>%
      mutate(
        Significance = case_when(
          Fold_Change >= fcUP   & Adjusted_P_Value < adjpvalue ~ "Upregulated",
          Fold_Change <= fcDown & Adjusted_P_Value < adjpvalue ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      ) %>%
      dplyr::arrange(desc(Fold_Change)) # Sort by highest FC first

    # volcanoPlotResults$VolcanoData[[group_name]] <- volcanoPlotData # UPdate list()
    volcanoPlotResults[[paste0("VolcanoData_", group_name)]] <- volcanoPlotData # UPdate list()


    # Filter the merged the datasets
    volcanoPlotData_filtered <- dplyr::full_join(fc_data, ca_results, by = "Feature") %>%
      # dplyr::full_join(ca_results, by = "Feature") %>%
      dplyr::filter(fold_change >= fcUP | fold_change <= fcDown,
                    `adj. p-value` < adjpvalue) %>%
      dplyr::rename(
        c("Fold_Change" = fold_change,
          "Log2_Fold_Change" = log2_fc,
          "P_Value" = `p-value`,
          "Adjusted_P_Value" = `adj. p-value`)
      ) %>%
      mutate(
        Sig = case_when(
          Fold_Change >= fcUP   & Adjusted_P_Value < adjpvalue ~ "Up",
          Fold_Change <= fcDown & Adjusted_P_Value < adjpvalue ~ "Down",
          TRUE ~ "Not Sig"
        )
      ) %>%
      dplyr::arrange(desc(Fold_Change)) # Sort by highest FC first

    # volcanoPlotResults$VolcanoData_filtered[[group_name]] <- volcanoPlotData_filtered # UPdate list()
    volcanoPlotResults[[paste0("VolcanoData_filtered_", group_name)]] <- volcanoPlotData_filtered # UPdate list()


    if (nrow(volcanoPlotData) == 0) {
      message("No features passed filtering for ", group_name)
      next
    }

    # Plot
    p <- ggplot(volcanoPlotData, aes(x = Log2_Fold_Change,
                                     y = -log10(Adjusted_P_Value))) +
      geom_point(aes(color = Significance), alpha = 0.5) + # alpha is how strong the color is
      scale_color_manual(values = c("Not Significant" = "grey",
                                    "Upregulated"     = "red",
                                    "Downregulated"   = "blue")) +
      theme_minimal() +
      labs(
        title = paste0("Volcano Plot (", group_name,") "),
        subtitle = paste0("Log2 Fold Change ≥ ", log2(fcUP) %>% round(2),
                          " or ≤ ", log2(fcDown) %>% round(2),
                          "; adjusted p-value < ", adjpvalue),
        x = "Log2 Fold Change",
        y = "Log10 adjusted p-value",
        caption = paste0("Fold Change ≥ ", fcUP,
                         " or ≤ ", fcDown)
      ) +
      geom_hline(yintercept = -log10(adjpvalue), linetype = "dashed", color = "grey", linewidth = 1) +
      geom_vline(xintercept = c(log2(fcDown), log2(fcUP)), linetype = "dashed", color = "grey", linewidth = 1)

    volcanoPlotResults$VolcanoPlot[[group_name]] <- p

    print(p) # print volcano plot

    # print(volcanoPlotData_filtered %>%
    #         dplyr::select(
    #           Feature, Fold_Change, Log2_Fold_Change, `adj. p-value2`, Sig
    #         )) # Print the filtered data



    # '#####################################'
    # '#####################################'
    #
    #
    #
    # filtered_features <- volcanoPlotData$Feature
    #
    # new_data <- PPData$data_scaledOPLSDA[non_qc_indices, ] %>%
    #   dplyr::select(any_of(filtered_features)) %>%
    #   cbind(Groups = groups, .)
    #
    # volcanoPlotResults$data_filtered[[group_name]] <- new_data # Update list()


  }

  return(volcanoPlotResults)


}


# Usage
# myvolcano <- plotVolcano(
#   FCData    = myfoldchange,
#   CAData    = myComparative,
#   fcUP      = 1.01,
#   fcDown    = 0.97,
#   adjpvalue = 0.05
# )

