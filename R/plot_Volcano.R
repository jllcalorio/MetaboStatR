#' Plot Volcano Plot
#'
#' @description
#' This function plots a volcano plot.
#'
#' @param PPData List. This list must be a result from the `perform_PreprocessingPeakData` function.
#' @param FCData List. This list must be a result from the `perform_FoldChange` function.
#' @param CAData List. This list must be a result from the `perform_ComparativeAnalysis` function.
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order. Suggests to be inputted (control, case1, case2) where case2 is worse than case1 e.g., Severe dengue than Non-severe dengue.
#' @param fcUP Numeric. The threshold for the upper fold change to be considered.
#' @param fcDown Numeric. The threshold for the lower fold change to be considered.
#' @param adjpvalue Numeric. The adjusted p-value threshold.
#'
#' @returns Returns a list of data used in the analysis, volcano plot data, and the plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_Volcano(
#'   PPData    = results_from_perform_PreprocessingPeakData_function,
#'   FCData    = results_from_perform_FoldChange_function,
#'   CAData    = results_from_perform_ComparativeAnalysis_function
#' )
#' }
#'
plot_Volcano <- function(
    PPData,
    FCData,
    CAData,
    arrangeLevels = NULL,
    fcUP          = 2,
    fcDown        = 0.5,
    adjpvalue     = 0.05
) {
  volcanoPlotResults                         <- list()
  volcanoPlotResults$FunctionOrigin          <- "plot_Volcano"
  volcanoPlotResults$FoldChangeData          <- FCData
  volcanoPlotResults$ComparativeAnalysisData <- CAData
  volcanoPlotResults$FoldChange_Upper        <- fcUP
  volcanoPlotResults$FoldChange_Down         <- fcDown
  volcanoPlotResults$AdjustedPValue          <- adjpvalue
  qc_indices                                 <- PPData$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices                             <- !qc_indices
  # Define the groups
  groups                                     <- PPData$Metadata$Group[non_qc_indices]
  # Check if arrangeLevels have the same actual values of groups
  if (!base::is.null(arrangeLevels)) {
    if (base::length(base::setdiff(arrangeLevels, base::unique(groups))) > 0) {
      stop("Invalid arrangeLevels: check for typos or missing groups.")
    } else {
      groups <- base::factor(groups, levels = arrangeLevels) %>% base::levels()
    }
  }
  unique_groups <- base::unique(groups)
  group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)

  # Initialize VolcanoPlot as a list
  volcanoPlotResults$VolcanoPlot <- list()

  for (group_pair in group_combinations) {
    group1     <- group_pair[1] # The group here is originally group1, but since ROC requires (control, case), thus
    group2     <- group_pair[2] # The group here is originally group2, but since ROC requires (control, case), thus
    group_name <- base::paste(group1, "vs.", group2, sep = " ")
    # Extract matching datasets from DR, FC, and CA
    # From FC
    if (!base::is.null(FCData[[base::paste0("data_combined_", group_name)]])) {
      fc_data    <- FCData[[base::paste0("data_combined_", group_name)]] %>%
        base::cbind(Feature = base::rownames(.), .) %>%
        base::`rownames<-`(NULL) %>%
        dplyr::arrange(Feature)
    } else {
      group_name <- base::paste(group2, "vs.", group1, sep = " ") # When NULL, change the order of groups
      fc_data    <- FCData[[base::paste0("data_combined_", group_name)]] %>%
        base::cbind(Feature = base::rownames(.), .) %>%
        base::`rownames<-`(NULL) %>%
        dplyr::arrange(Feature)
    }
    # From CA
    ca_results <- CAData$results %>%
      base::cbind(Feature = base::rownames(.), .) %>%
      base::`rownames<-`(NULL) %>%
      dplyr::arrange(Feature)
    # Merge the datasets
    volcanoPlotData <- dplyr::full_join(fc_data, ca_results, by = "Feature") %>%
      dplyr::rename(
        Fold_Change      = fold_change,
        Log2_Fold_Change = log2_fc,
        P_Value          = `p-value`,
        Adjusted_P_Value = `adj. p-value`
      ) %>%
      dplyr::mutate(
        Significance = dplyr::case_when(
          Fold_Change >= fcUP   & Adjusted_P_Value < adjpvalue ~ "Upregulated",
          Fold_Change <= fcDown & Adjusted_P_Value < adjpvalue ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      ) %>%
      dplyr::arrange(desc(Fold_Change))                                               # Sort by highest FC first
    # volcanoPlotResults$VolcanoData[[group_name]] <- volcanoPlotData # UPdate list()
    volcanoPlotResults[[base::paste0("VolcanoData_", group_name)]] <- volcanoPlotData # Update list()
    # Filter the merged the datasets
    volcanoPlotData_filtered <- dplyr::full_join(fc_data, ca_results, by = "Feature") %>%
      # dplyr::full_join(ca_results, by = "Feature") %>%
      dplyr::filter(fold_change >= fcUP | fold_change <= fcDown,
                    `adj. p-value` < adjpvalue) %>%
      dplyr::rename(
        Fold_Change      = fold_change,
        Log2_Fold_Change = log2_fc,
        P_Value          = `p-value`,
        Adjusted_P_Value = `adj. p-value`
      ) %>%
      dplyr::mutate(
        Sig = dplyr::case_when(
          Fold_Change >= fcUP   & Adjusted_P_Value < adjpvalue ~ "Up",
          Fold_Change <= fcDown & Adjusted_P_Value < adjpvalue ~ "Down",
          TRUE ~ "Not Sig"
        )
      ) %>%
      dplyr::arrange(dplyr::desc(Fold_Change))                                                          # Sort by highest FC first
    volcanoPlotResults[[base::paste0("VolcanoData_filtered_", group_name)]] <- volcanoPlotData_filtered # Update list()
    if (nrow(volcanoPlotData) == 0) {
      message("No features passed filtering for ", group_name)
      next
    }
    # Plot
    p <- ggplot2::ggplot(volcanoPlotData, ggplot2::aes(x = Log2_Fold_Change,
                                                       y = -base::log10(Adjusted_P_Value))) +
      ggplot2::geom_point(ggplot2::aes(color = Significance), alpha = 0.5) +                           # alpha is how strong the color is
      ggplot2::scale_color_manual(values = base::c("Not Significant" = "grey",
                                                   "Upregulated"     = "red",
                                                   "Downregulated"   = "blue")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title    = base::paste0("Volcano Plot (", group_name,") "),
        subtitle = base::paste0("Log2 Fold Change >= ", base::log2(fcUP) %>% base::round(2),
                                " or <= ", base::log2(fcDown) %>% base::round(2),
                                "; adjusted p-value < ", adjpvalue),
        x        = "Log2 Fold Change",
        y        = "Log10 adjusted p-value",
        caption  = base::paste0("Fold Change >= ", fcUP,
                               " or <= ", fcDown)
      ) +
      ggplot2::geom_hline(yintercept = -base::log10(adjpvalue),
                          linetype   = "dashed",
                          color      = "grey",
                          linewidth  = 1) +
      ggplot2::geom_vline(xintercept = base::c(base::log2(fcDown),
                                               base::log2(fcUP)),
                          linetype   = "dashed",
                          color      = "grey",
                          linewidth  = 1)
    volcanoPlotResults$VolcanoPlot[[group_name]] <- p
    print(p) # print volcano plot
  }
  return(volcanoPlotResults)
}
