#' Plot Random Data Before and After Preprocessing
#'
#' @description
#' This function will plot random samples or features the before and after data preprocessing.
#'
#'
#' @param data List. This list must be a result from the `performPreprocessingPeakData` function.
#' @param scaled String. Choose one below to plot one of the final results from the data preprocessing.
#'   \itemize{
#'     \item 'PCA: Plots results from PCA.
#'     \item 'OPLS-DA': Plots results from OPLS-DA.
#'     }
#'     Defaults to "OPLS-DA".
#' @param group_by String. Whether to plot the "Samples" or the "Features".
#'   \itemize{
#'     \item 'Sample: Plots the samples.
#'     \item 'Feature': Plots the features.
#'     }
#'     Defaults to "Sample".
#' @param n_random_samples Numeric. The number of random 'samples' to be created with box plots. If `NULL`, plots all. Defaults to 30, or the samples exceeds 30.
#' @param n_random_features Numeric. he number of random 'features' to be created with box plots. If `NULL`, plots all. Defaults to 30, or the features exceeds 30.
#'
#' @returns A list of plots of the before and after data preprocessing.
#' @export
#'
#' @examples
#' \dontrun{
#' mplotBeforeAfter(
#'   data     = results_from_performPreprocessingPeakData_function,
#'   scaled   = "OPLS-DA",
#'   group_by = "Sample"
#' )
#' }

plotBeforeAfter <- function(
    data,
    scaled = "OPLS-DA",
    group_by = "Sample",
    n_random_samples = 30,
    n_random_features = 30
) {

  # NOTE
  # THIS IS NOT YET WORKING

  # Visualize before and after normalization

  # Combined Plot Function

  # Install/Load required package
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' is required but not installed.")
    message("Installing the package...")
    install.packages('patchwork')
    message("Loading the package...")
    library(patchwork)
  } else {
    library(patchwork)
  }

  # Parameter checks
  if (data$Metadata$FunctionOrigin != "performPreprocessingPeakData") {
    stop("Data did not come from 'performProcessingPeakData' function.")
  }

  if (!is.null(scaled)) {
    if (scaled != "PCA" && scaled != "OPLS-DA") {
      stop("The 'scaled' parameter must be either PCA or OPLS-DA.")
    }
  }

  beforeAfterResults <- list() # Empty list to store results
  beforeAfterResults$Class <- "plotBeforeAfter" # Update list

  # Prepare Data for Visualization

  non_qc_indices <- data$Metadata$Groups != "QC" # Non-QC indices

  original <- data$data_no_NA[non_qc_indices, ] # Select this as original data

  if (scaled == "PCA") {
    transformed <- data$data_scaledPCA[non_qc_indices, ]
  } else if (scaled == "OPLS-DA") {
    transformed <- data$data_scaledOPLSDA[non_qc_indices, ]
  } else {
    stop("scaled: Must be either PCA or OPLS-DA")
  }

  beforeAfterResults$data_before <- original # Update list
  beforeAfterResults$data_after <- transformed # Update list

  # Convert data to long format for ggplot2
  long_format <- function(data) {
    data %>%
      as.data.frame() %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(-Sample, names_to = "Feature", values_to = "Value")
  }

  # Convert data to long format
  orig_long  <- long_format(original)
  trans_long <- long_format(transformed)

  # Density plots

  p1 <- ggplot(orig_long, aes(x = Value)) +
    geom_density(alpha = 0.5) +
    # ggtitle(paste0(title_prefix, " - Before Scaling")) +
    theme_minimal() +
    theme(legend.position = "none")

  p2 <- ggplot(trans_long, aes(x = Value)) +
    geom_density(alpha = 0.5) +
    # ggtitle(paste0(title_prefix, " - After Scaling")) +
    theme_minimal() +
    theme(legend.position = "none")

  # Boxplots

  random_samples <- sample(rownames(transformed), n_random_samples)  # Select random features

  if (group_by == "Sample") {

    # update random_features if n_random_features is NULL or not
    if (!is.null(n_random_samples)) {
      random_samples <- sample(colnames(transformed), n_random_samples)  # Select random features
    } else {
      random_samples <- dim(transformed)[1]
    }

    p3 <- ggplot(orig_long %>% filter(Sample %in% random_samples),
                 aes(x = Sample, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot Before Scaling")) +
      theme_minimal() +
      coord_flip()

    p4 <- ggplot(trans_long %>% filter(Sample %in% random_samples),
                 aes(x = Sample, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot After Scaling")) +
      theme_minimal() +
      coord_flip()
  } else if (group_by == "Feature") {

    # update random_features if n_random_features is NULL or not
    if (!is.null(n_random_features)) {
      random_features <- sample(colnames(transformed), n_random_features)  # Select random features
    } else {
      random_features <- dim(transformed)[2]
    }

    p3 <- ggplot(orig_long %>% filter(Feature %in% random_features),
                 aes(x = Feature, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot Before Scaling")) +
      theme_minimal() +
      coord_flip()

    p4 <- ggplot(trans_long %>% filter(Feature %in% random_features),
                 aes(x = Feature, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot After Scaling")) +
      theme_minimal() +
      coord_flip()
  } else {
    stop("The 'group_by' parameter must be 'Sample' or 'Feature'")
  }

  beforeAfterResults$plot_density_before <- p1 # Update list
  beforeAfterResults$plot_density_after  <- p2 # Update list
  beforeAfterResults$plot_box_before     <- p3 # Update list
  beforeAfterResults$plot_box_after      <- p4 # Update list

  # Combine plots
  all_plots <- (p1 | p2) / (p3 | p4)

  beforeAfterResults$plot_4in1           <- all_plots # Update list

  return(beforeAfterResults)
}
