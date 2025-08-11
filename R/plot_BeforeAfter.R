#' Plot Random Data Before and After Preprocessing
#'
#' @description
#' This function will plot random samples or features the before and after data preprocessing.
#'
#'
#' @param data List. This list must be a result from the `perform_PreprocessingPeakData` function.
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
#' @param n_random_features Numeric. The number of random 'features' to be created with box plots. If `NULL`, plots all. Defaults to 30, or the features exceeds 30.
#'
#' @returns A list of plots of the before and after data preprocessing.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_BeforeAfter(
#'   data     = results_from_perform_PreprocessingPeakData_function,
#'   scaled   = "OPLS-DA",
#'   group_by = "Sample"
#' )
#' }

plot_BeforeAfter <- function(
    data,
    scaled = "OPLS-DA",
    group_by = "Sample",
    n_random_samples = 30,
    n_random_features = 30
) {

  # Load required libraries
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required but not installed.")
  }

  # Parameter checks
  if (data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("Data did not come from 'performProcessingPeakData' function.")
  }

  if (!base::is.null(scaled)) {
    if (scaled != "PCA" && scaled != "OPLS-DA") {
      stop("The 'scaled' parameter must be either PCA or OPLS-DA.")
    }
  }

  if (!group_by %in% c("Sample", "Feature")) {
    stop("The 'group_by' parameter must be 'Sample' or 'Feature'")
  }

  beforeAfterResults                <- base::list() # Empty list to store results
  beforeAfterResults$FunctionOrigin <- "plot_BeforeAfter" # Update list

  # Prepare Data for Visualization
  non_qc_indices <- data$Metadata$Group_ != "QC" # Non-QC indices
  original <- data$data_no_NA[non_qc_indices, ] # Select this as original data

  if (scaled == "PCA") {
    transformed <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ]
  } else if (scaled == "OPLS-DA") {
    transformed <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, ]
  } else {
    stop("scaled: Must be either PCA or OPLS-DA")
  }

  beforeAfterResults$data_before <- original # Update list
  beforeAfterResults$data_after  <- transformed # Update list

  # Convert data to long format for ggplot2
  long_format <- function(data) {
    data %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "Sample") %>%
      tidyr::pivot_longer(-Sample, names_to = "Feature", values_to = "Value")
  }

  # Convert data to long format
  orig_long  <- long_format(original)
  trans_long <- long_format(transformed)

  # Create plots based on group_by parameter
  if (group_by == "Sample") {

    # Density plots for all data (sample-wise perspective)
    p1 <- ggplot2::ggplot(orig_long, ggplot2::aes(x = Value)) +
      ggplot2::geom_density(alpha = 0.5,
                            fill  = "blue") +
      ggplot2::ggtitle("Density - Before Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    p2 <- ggplot2::ggplot(trans_long, ggplot2::aes(x = Value)) +
      ggplot2::geom_density(alpha = 0.5,
                            fill  = "red") +
      ggplot2::ggtitle("Density - After Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    # Boxplots for random samples
    n_samples <- base::nrow(transformed)
    if (base::is.null(n_random_samples)) {
      selected_samples <- base::rownames(transformed)
    } else {
      min_samples <- base::min(n_samples, n_random_samples)
      selected_samples <- base::sample(base::rownames(transformed), min_samples)
    }

    p3 <- ggplot2::ggplot(orig_long %>% dplyr::filter(Sample %in% selected_samples),
                          ggplot2::aes(x = Sample, y = Value)) +
      ggplot2::geom_boxplot() +
      ggplot2::ggtitle("Boxplot (Samples) - Before Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

    p4 <- ggplot2::ggplot(trans_long %>% dplyr::filter(Sample %in% selected_samples),
                          ggplot2::aes(x = Sample,
                                       y = Value)) +
      ggplot2::geom_boxplot() +
      ggplot2::ggtitle("Boxplot (Samples) - After Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

  } else if (group_by == "Feature") {

    # For feature-wise analysis, we want to see distribution of each feature across samples
    n_features <- base::ncol(transformed)
    if (base::is.null(n_random_features)) {
      selected_features <- base::colnames(transformed)
    } else {
      min_features      <- base::min(n_features, n_random_features)
      selected_features <- base::sample(base::colnames(transformed), min_features)
    }

    # Density plots for selected features (feature-wise perspective)
    orig_long_subset  <- orig_long  %>% dplyr::filter(Feature %in% selected_features)
    trans_long_subset <- trans_long %>% dplyr::filter(Feature %in% selected_features)

    p1 <- ggplot2::ggplot(orig_long_subset, ggplot2::aes(x = Value)) +
      ggplot2::geom_density(alpha = 0.5,
                            fill  = "blue") +
      ggplot2::ggtitle("Density (Features) - Before Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    p2 <- ggplot2::ggplot(trans_long_subset, ggplot2::aes(x = Value)) +
      ggplot2::geom_density(alpha = 0.5,
                            fill  = "red") +
      ggplot2::ggtitle("Density (Features) - After Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")

    # Boxplots for selected features
    p3 <- ggplot2::ggplot(orig_long_subset,
                          ggplot2::aes(x = Feature,
                                       y = Value)) +
      ggplot2::geom_boxplot() +
      ggplot2::ggtitle("Boxplot (Features) - Before Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

    p4 <- ggplot2::ggplot(trans_long_subset,
                          ggplot2::aes(x = Feature,
                                       y = Value)) +
      ggplot2::geom_boxplot() +
      ggplot2::ggtitle("Boxplot (Features) - After Preprocessing") +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  }

  # Store individual plots
  beforeAfterResults$plot_density_before <- p1 # Update list
  beforeAfterResults$plot_density_after  <- p2 # Update list
  beforeAfterResults$plot_box_before     <- p3 # Update list
  beforeAfterResults$plot_box_after      <- p4 # Update list

  # Combine plots
  all_plots <- gridExtra::grid.arrange(grobs = list(p1, p2, p3, p4),
                                       nrow  = 2,
                                       ncol  = 2,
                                       top   = paste("Data Visualization:", group_by, "Perspective"))

  beforeAfterResults$plot_4in1 <- all_plots # Update list

  return(beforeAfterResults)
}
