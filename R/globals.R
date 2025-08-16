#' @importFrom magrittr %>%
NULL

#' Global variables to avoid R CMD check notes
#'
#' These variables are used in dplyr/ggplot2 code where they refer to column names
#' in data frames, not actual global variables.
#' @noRd
#'
globalVariables(c(
  ".", "Feature", "VIP", "fold_change", "adj. p-value",
  "AUROC", "AUC", "Label", "PC1", "PC2", "Batch", "Group",
  "Groups", "Injection", "Sample", "Value", "Abundance",
  "Covariance", "Correlation", "before", "after", "correction",
  "sample_id", "intensity", "log2_fc", "p-value", "Fold_Change",
  "Log2_Fold_Change", "Adjusted_P_Value", "Significance",
  "scalePCA", "base", "purrr", "desc", "y", "x", "labeller",

  "group", "value", "Loading", "Comp1", "Comp2", "PC_num", "batch", "var", "sd"
))
