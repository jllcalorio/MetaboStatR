#' Perform Fold Change Analysis on a Preprocessed Data
#'
#' @description
#' This function performs fold change analysis on a preprocessed Data.
#'
#' @param data List. This list must be a result from the `performPreprocessingPeakData` function. It can also be from a data scaling technique.
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order.
#' @param sortFC Boolean. If `TRUE` (default), sorts the fold changes in descending order.
#'
#' @returns A list of data frames.
#' @export
#'
#' @examples
#' \dontrun{
#' performFoldChange(data = results_from_performPreprocessingPeakData_function)
#'}
#'
performFoldChange <- function(
    data,
    arrangeLevels = NULL,
    sortFC        = TRUE
) {

  # Add some data checks
  foldChangeAnalysisResults       <- list() # Empty list to store all results
  foldChangeAnalysisResults$Class <- "performFoldChange"  # Update list

  non_qc_indices <- data$Metadata$Groups != "QC"

  # Filter out QC samples
  df     <- data$data_scaledOPLSDA[non_qc_indices, ]

  # Parameter-check for arrangeLevels
  if (is.null(arrangeLevels)) {
    groups <- data$Metadata$Groups[non_qc_indices]
  } else if (!is.null(arrangeLevels)) {
    if (length(setdiff(arrangeLevels, data$Metadata$Groups[non_qc_indices])) > 0) {
      stop("Check 'arrangeLevels' values. There might be (1) typos (2) Group not in the data or (3) Group is missing.")
    } else {
      groups <- data$Metadata$Groups[non_qc_indices] %>% factor(levels = arrangeLevels)
    }
  }

  # Add a small constant to avoid division by zero or log(0)
  if (any(df <= 0)) {
    shift_value <- 1 - min(df) # find minimum then subtract it to find the max value to add to the smallest value
    df <- df + shift_value     # Add the shift value so the minimum value = 1, since log2 of <= 0 is -Inf and NaNs
  }

  foldChangeAnalysisResults$data_Min_is_1 <- df %>% as.data.frame()  # Update list

  # Get unique group combinations
  unique_groups      <- levels(groups)
  group_combinations <- combn(unique_groups, 2, simplify = FALSE)

  for (pair in group_combinations) {

    group1 <- pair[1]
    group2 <- pair[2]
    comparison_label <- paste(group1, "vs.", group2, sep = " ")

    # Subset data for the current pair
    idx_group1 <- which(groups == group1)
    idx_group2 <- which(groups == group2)

    # Compute Fold Change and Log2 Fold Change
    if (sortFC == TRUE) {
      fold_change <- apply(df, 2, function(x) mean(x[idx_group1]) / mean(x[idx_group2])) %>%
        as.data.frame() %>%
        `colnames<-`("fold_change") %>%
        dplyr::arrange(desc(fold_change))
    } else if (sortFC == FALSE) {
      fold_change <- apply(df, 2, function(x) mean(x[idx_group1]) / mean(x[idx_group2])) %>%
        as.data.frame() %>%
        `colnames<-`("fold_change")
    }

    log2_fc <- log2(fold_change) %>%
      as.data.frame() %>%
      `colnames<-`("log2_fc")

    # Combine results
    combined_fc <- cbind(fold_change, log2_fc)

    # Store results in the list using the comparison label
    # foldChangeAnalysisResults[[paste0("data_foldChange_", comparison_label)]] <- fold_change
    # foldChangeAnalysisResults[[paste0("data_log2foldChange_", comparison_label)]] <- log2_fc
    foldChangeAnalysisResults[[paste0("data_combined_", comparison_label)]] <- combined_fc
  }

  return(foldChangeAnalysisResults)
}
