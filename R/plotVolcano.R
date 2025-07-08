#' Plot Volcano Plot
#'
#' @description
#' This function plots a volcano plot.
#'
#' @param PPData List. This list must be a result from the `performPreprocessingPeakData` function.
#' @param FCData List. This list must be a result from the `performFoldChange` function.
#' @param CAData List. This list must be a result from the `performComparativeAnalysis` function.
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order. Suggests to be inputted (control, case1, case2) where case2 is worse than case1 e.g., Severe dengue than Non-severe dengue.
#' @param fcUP Numeric. The threshold for the upper fold change to be considered.
#' @param fcDown Numeric. The threshold for the lower fold change to be considered.
#' @param adjpvalue Numeric. The adjusted p-value threshold.
#'
#' @returns Resurns a list of data used in the analysis, volcano plot data, and the plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plotVolcano(
#'   PPData    = results_from_performPreprocessingPeakData_function,
#'   FCData    = results_from_performFoldChange_function,
#'   CAData    = results_from_performComparativeAnalysis_function
#' )
#' }
#'
plotVolcano <- function(
    PPData,
    FCData,
    CAData,
    arrangeLevels = NULL,
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

    # Extract matching datasets from DR, FC, and CA
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
    volcanoPlotResults[[paste0("VolcanoData_", group_name)]] <- volcanoPlotData # Update list()

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
  }

  return(volcanoPlotResults)
}
