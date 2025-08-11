#' Perform Fold Change Analysis on a Preprocessed Data
#'
#' @description
#' This function performs fold change analysis on a preprocessed Data.
#'
#' @param data List. This list must be a result from the `perform_PreprocessingPeakData` function. It can also be from a data scaling technique.
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order.
#' @param sortFC Boolean. If `TRUE` (default), sorts the fold changes in descending order.
#'
#' @returns A list of data frames.
#' @export
#'
#' @examples
#' \dontrun{
#' perform_FoldChange(data = results_from_perform_PreprocessingPeakData_function)
#'}
#'
perform_FoldChange <- function(
    data,
    arrangeLevels = NULL,
    sortFC        = TRUE
) {

  # Add some data checks
  foldChangeAnalysisResults                <- base::list() # Empty list to store all results
  foldChangeAnalysisResults$FunctionOrigin <- "perform_FoldChange"  # Update list

  qc_indices                               <- data$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices                           <- !qc_indices

  # Filter out QC samples
  df                                       <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ]

  # Parameter-check for arrangeLevels
  if (base::is.null(arrangeLevels)) {
    groups <- data$Metadata$Group[non_qc_indices]
  } else if (!base::is.null(arrangeLevels)) {
    if (base::length(base::setdiff(arrangeLevels, data$Metadata$Group[non_qc_indices])) > 0) {
      stop("Check 'arrangeLevels' values. There might be (1) typos (2) Group not in the data or (3) Group is missing.")
    } else {
      groups <- data$Metadata$Group[non_qc_indices] %>% factor(levels = arrangeLevels)
    }
  }

  # Add a small constant to avoid division by zero or log(0)
  if (any(df <= 0)) {
    shift_value <- 1 - base::min(df) # find minimum then subtract it to find the max value to add to the smallest value
    df <- df + shift_value     # Add the shift value so the minimum value = 1, since log2 of <= 0 is -Inf and NaNs
  }

  foldChangeAnalysisResults$data_Min_is_1 <- df %>% base::as.data.frame()  # Update list

  # Get unique group combinations
  unique_groups      <- base::levels(base::factor(groups))
  group_combinations <- utils::combn(unique_groups, 2, simplify = FALSE)

  for (pair in group_combinations) {

    group1           <- pair[1]
    group2           <- pair[2]
    comparison_label <- base::paste(group1, "vs.", group2, sep = " ")

    # Subset data for the current pair
    idx_group1 <- base::which(groups == group1)
    idx_group2 <- base::which(groups == group2)

    # Compute Fold Change and Log2 Fold Change
    fold_change <- base::apply(df, 2, function(x) base::mean(x[idx_group1]) / base::mean(x[idx_group2])) %>%
      base::as.data.frame() %>%
      base::`colnames<-`("fold_change")

    if (sortFC == TRUE) {
      fold_change <- dplyr::arrange(fold_change, desc(fold_change))
    }

    # Calculate log2 of fold changes
    log2_fc <- base::log2(fold_change) %>%
      base::as.data.frame() %>%
      base::`colnames<-`("log2_fc")

    # Combine results
    combined_fc <- base::cbind(fold_change, log2_fc)

    # Store results in the list using the comparison label
    foldChangeAnalysisResults[[base::paste0("data_combined_", comparison_label)]] <- combined_fc
  }

  return(foldChangeAnalysisResults)
}
