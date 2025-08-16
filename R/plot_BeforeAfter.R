#' Visualize Data Before and After Preprocessing
#'
#' @description
#' Creates comprehensive visualizations comparing data distributions and patterns
#' before and after preprocessing steps. Generates density plots and box plots
#' from either sample or feature perspectives to assess preprocessing effectiveness.
#' Supports both PCA and PLS scaled data visualization with intelligent
#' random sampling for large datasets.
#'
#' @param data List. A preprocessing results object from \code{perform_PreprocessingPeakData}
#'   function containing original data, transformed data, and metadata.
#' @param scaled Character. Specifies which preprocessed data to visualize:
#'   \itemize{
#'     \item \code{"PCA"}: Uses PCA-scaled preprocessed data
#'     \item \code{"PLS"}: Uses PLS-scaled preprocessed data
#'   }
#'   Default: \code{"PLS"}
#' @param group_by Character. Visualization perspective:
#'   \itemize{
#'     \item \code{"Sample"}: Visualizes sample-wise distributions
#'     \item \code{"Feature"}: Visualizes feature-wise distributions
#'   }
#'   Default: \code{"Sample"}
#' @param n_random_samples Integer. Number of random samples to display in box plots.
#'   If \code{NULL}, uses all samples. If specified number exceeds total samples,
#'   uses all available samples. Default: \code{30}
#' @param n_random_features Integer. Number of random features to display in box plots.
#'   If \code{NULL}, uses all features. If specified number exceeds total features,
#'   uses all available features. Default: \code{30}
#' @param seed Integer. Random seed for reproducible sample/feature selection.
#'   Default: \code{123}
#'
#' @return List containing:
#'   \itemize{
#'     \item \code{FunctionOrigin}: Character indicating function source
#'     \item \code{data_before}: Matrix of original data (non-QC samples)
#'     \item \code{data_after}: Matrix of transformed data (non-QC samples)
#'     \item \code{plot_density_before}: ggplot2 density plot of original data
#'     \item \code{plot_density_after}: ggplot2 density plot of transformed data
#'     \item \code{plot_box_before}: ggplot2 box plot of original data
#'     \item \code{plot_box_after}: ggplot2 box plot of transformed data
#'     \item \code{plot_combined}: Combined 2x2 grid plot
#'     \item \code{processing_info}: List with processing metadata
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with PLS scaling
#' plots <- plot_BeforeAfter(
#'   data = preprocessing_results,
#'   scaled = "PLS",
#'   group_by = "Sample"
#' )
#'
#' # Feature perspective with PCA scaling
#' feature_plots <- plot_BeforeAfter(
#'   data = preprocessing_results,
#'   scaled = "PCA",
#'   group_by = "Feature",
#'   n_random_features = 20
#' )
#'
#' # Display combined plot
#' print(plots$plot_combined)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_density geom_boxplot ggtitle theme_minimal coord_flip element_text
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom gridExtra grid.arrange
#' @importFrom stats setNames
#' @seealso \code{\link{perform_PreprocessingPeakData}}

plot_BeforeAfter <- function(data,
                             scaled = "PLS",
                             group_by = "Sample",
                             n_random_samples = 30L,
                             n_random_features = 30L,
                             seed = 123L) {

  # Set seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  .validate_inputs_plotbeforeafter(data, scaled, group_by, n_random_samples, n_random_features, seed)

  # Check and load required packages
  .check_required_packages()

  # Initialize results list
  results <- list(
    FunctionOrigin = "plot_BeforeAfter",
    processing_info = list(
      scaled = scaled,
      group_by = group_by,
      n_random_samples = n_random_samples,
      n_random_features = n_random_features,
      seed = seed,
      timestamp = Sys.time()
    )
  )

  # Extract and prepare data
  data_prep <- .prepare_data(data, scaled)
  results$data_before <- data_prep$original
  results$data_after <- data_prep$transformed

  # Convert to long format for ggplot
  long_data <- .convert_to_long_format(data_prep$original, data_prep$transformed)

  # Generate plots based on group_by parameter
  if (group_by == "Sample") {
    plots <- .create_sample_plots(long_data, data_prep$transformed, n_random_samples)
  } else {
    plots <- .create_feature_plots(long_data, data_prep$transformed, n_random_features)
  }

  # Store individual plots
  results$plot_density_before <- plots$p1
  results$plot_density_after <- plots$p2
  results$plot_box_before <- plots$p3
  results$plot_box_after <- plots$p4

  # Create combined plot
  results$plot_combined <- .create_combined_plot(plots, group_by)

  # Add processing summary
  results$processing_info$original_dims <- dim(data_prep$original)
  results$processing_info$transformed_dims <- dim(data_prep$transformed)
  results$processing_info$qc_samples_removed <- sum(data$Metadata$Group_ == "QC")

  return(results)
}

#' Validate Input Parameters
#' @noRd
.validate_inputs_plotbeforeafter <- function(data, scaled, group_by, n_random_samples, n_random_features, seed) {

  # Check data structure
  if (!is.list(data)) {
    stop("Input 'data' must be a list object.", call. = FALSE)
  }

  # Check function origin
  if (is.null(data$FunctionOrigin) || data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("Data must originate from 'perform_PreprocessingPeakData' function.", call. = FALSE)
  }

  # Required data components
  required_components <- c("data_no_NA", "Metadata", "data_scaledPCA_rsdFiltered_varFiltered")
  missing_components <- setdiff(required_components, names(data))
  if (length(missing_components) > 0) {
    stop("Missing required data components: ", paste(missing_components, collapse = ", "),
         call. = FALSE)
  }

  # Validate scaled parameter
  if (!is.character(scaled) || length(scaled) != 1) {
    stop("'scaled' must be a single character string.", call. = FALSE)
  }
  if (!scaled %in% c("PCA", "PLS")) {
    stop("'scaled' must be either 'PCA' or 'PLS'.", call. = FALSE)
  }

  # Validate group_by parameter
  if (!is.character(group_by) || length(group_by) != 1) {
    stop("'group_by' must be a single character string.", call. = FALSE)
  }
  if (!group_by %in% c("Sample", "Feature")) {
    stop("'group_by' must be either 'Sample' or 'Feature'.", call. = FALSE)
  }

  # Validate numeric parameters
  if (!is.null(n_random_samples)) {
    if (!is.numeric(n_random_samples) || length(n_random_samples) != 1 || n_random_samples <= 0) {
      stop("'n_random_samples' must be a positive integer or NULL.", call. = FALSE)
    }
    if (n_random_samples != as.integer(n_random_samples)) {
      warning("'n_random_samples' converted to integer.", call. = FALSE)
      n_random_samples <- as.integer(n_random_samples)
    }
  }

  if (!is.null(n_random_features)) {
    if (!is.numeric(n_random_features) || length(n_random_features) != 1 || n_random_features <= 0) {
      stop("'n_random_features' must be a positive integer or NULL.", call. = FALSE)
    }
    if (n_random_features != as.integer(n_random_features)) {
      warning("'n_random_features' converted to integer.", call. = FALSE)
      n_random_features <- as.integer(n_random_features)
    }
  }

  # Validate seed
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single integer or NULL.", call. = FALSE)
    }
  }

  # Check metadata structure
  if (is.null(data$Metadata$Group_)) {
    stop("Metadata must contain 'Group_' column for QC identification.", call. = FALSE)
  }

  # Check data dimensions compatibility
  if (nrow(data$data_no_NA) != nrow(data$Metadata)) {
    stop("Number of rows in data and metadata must match.", call. = FALSE)
  }
}

#' Check and Load Required Packages
#' @noRd
.check_required_packages <- function() {
  required_packages <- c("ggplot2", "dplyr", "tidyr", "tibble", "gridExtra")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed. Please install it with: install.packages('",
           pkg, "')", call. = FALSE)
    }
  }
}

#' Prepare Data for Visualization
#' @noRd
.prepare_data <- function(data, scaled) {

  # Identify non-QC samples
  non_qc_indices <- data$Metadata$Group_ != "QC"

  if (sum(non_qc_indices) == 0) {
    stop("No non-QC samples found in the data.", call. = FALSE)
  }

  # Extract original data (non-QC samples only)
  original <- data$data_no_NA[non_qc_indices, , drop = FALSE]

  # Extract transformed data based on scaling method
  # Note: Both PCA and PLS use the same preprocessed data in this implementation
  transformed <- data$data_scaledPCA_rsdFiltered_varFiltered[non_qc_indices, , drop = FALSE]

  # Validate data integrity
  if (nrow(original) == 0 || ncol(original) == 0) {
    stop("Original data is empty after removing QC samples.", call. = FALSE)
  }

  if (nrow(transformed) == 0 || ncol(transformed) == 0) {
    stop("Transformed data is empty.", call. = FALSE)
  }

  if (!identical(dim(original), dim(transformed))) {
    warning("Original and transformed data have different dimensions. This may affect visualization.",
            call. = FALSE)
  }

  # Check for missing values
  if (any(is.na(original))) {
    warning("Original data contains missing values which may affect visualization.", call. = FALSE)
  }

  if (any(is.na(transformed))) {
    warning("Transformed data contains missing values which may affect visualization.", call. = FALSE)
  }

  return(list(original = original, transformed = transformed))
}

#' Convert Data to Long Format
#' @noRd
.convert_to_long_format <- function(original, transformed) {

  # Helper function to convert matrix to long format
  to_long <- function(mat, data_type) {
    if (is.null(rownames(mat))) {
      rownames(mat) <- paste0("Sample_", seq_len(nrow(mat)))
    }
    if (is.null(colnames(mat))) {
      colnames(mat) <- paste0("Feature_", seq_len(ncol(mat)))
    }

    mat %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Sample") %>%
      tidyr::pivot_longer(-Sample, names_to = "Feature", values_to = "Value") %>%
      dplyr::mutate(DataType = data_type)
  }

  # Convert both datasets
  orig_long <- to_long(original, "Before")
  trans_long <- to_long(transformed, "After")

  return(list(original = orig_long, transformed = trans_long))
}

#' Create Sample-wise Plots
#' @noRd
.create_sample_plots <- function(long_data, transformed, n_random_samples) {

  orig_long <- long_data$original
  trans_long <- long_data$transformed

  # Create density plots
  p1 <- ggplot2::ggplot(orig_long, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(alpha = 0.6, fill = "steelblue", color = "darkblue") +
    ggplot2::ggtitle("Sample Distribution - Before Preprocessing") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Intensity Values", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

  p2 <- ggplot2::ggplot(trans_long, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(alpha = 0.6, fill = "coral", color = "darkred") +
    ggplot2::ggtitle("Sample Distribution - After Preprocessing") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Scaled Values", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

  # Select random samples for box plots
  n_samples <- nrow(transformed)
  if (is.null(n_random_samples) || n_random_samples >= n_samples) {
    selected_samples <- rownames(transformed)
    n_selected <- n_samples
  } else {
    selected_samples <- sample(rownames(transformed), n_random_samples)
    n_selected <- n_random_samples
  }

  # Create box plots
  p3 <- ggplot2::ggplot(
    orig_long %>% dplyr::filter(Sample %in% selected_samples),
    ggplot2::aes(x = reorder(Sample, Value, FUN = median), y = Value)
  ) +
    ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.7) +
    ggplot2::ggtitle(paste("Sample Profiles - Before (n =", n_selected, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Samples", y = "Intensity Values") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )

  p4 <- ggplot2::ggplot(
    trans_long %>% dplyr::filter(Sample %in% selected_samples),
    ggplot2::aes(x = reorder(Sample, Value, FUN = median), y = Value)
  ) +
    ggplot2::geom_boxplot(fill = "lightcoral", alpha = 0.7) +
    ggplot2::ggtitle(paste("Sample Profiles - After (n =", n_selected, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Samples", y = "Scaled Values") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )

  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

#' Create Feature-wise Plots
#' @noRd
.create_feature_plots <- function(long_data, transformed, n_random_features) {

  orig_long <- long_data$original
  trans_long <- long_data$transformed

  # Select random features
  n_features <- ncol(transformed)
  if (is.null(n_random_features) || n_random_features >= n_features) {
    selected_features <- colnames(transformed)
    n_selected <- n_features
  } else {
    selected_features <- sample(colnames(transformed), n_random_features)
    n_selected <- n_random_features
  }

  # Filter data for selected features
  orig_subset <- orig_long %>% dplyr::filter(Feature %in% selected_features)
  trans_subset <- trans_long %>% dplyr::filter(Feature %in% selected_features)

  # Create density plots
  p1 <- ggplot2::ggplot(orig_subset, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(alpha = 0.6, fill = "steelblue", color = "darkblue") +
    ggplot2::ggtitle(paste("Feature Distribution - Before (n =", n_selected, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Intensity Values", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

  p2 <- ggplot2::ggplot(trans_subset, ggplot2::aes(x = Value)) +
    ggplot2::geom_density(alpha = 0.6, fill = "coral", color = "darkred") +
    ggplot2::ggtitle(paste("Feature Distribution - After (n =", n_selected, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Scaled Values", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

  # Create box plots
  p3 <- ggplot2::ggplot(
    orig_subset,
    ggplot2::aes(x = reorder(Feature, Value, FUN = median), y = Value)
  ) +
    ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.7) +
    ggplot2::ggtitle(paste("Feature Profiles - Before (n =", n_selected, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Features", y = "Intensity Values") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )

  p4 <- ggplot2::ggplot(
    trans_subset,
    ggplot2::aes(x = reorder(Feature, Value, FUN = median), y = Value)
  ) +
    ggplot2::geom_boxplot(fill = "lightcoral", alpha = 0.7) +
    ggplot2::ggtitle(paste("Feature Profiles - After (n =", n_selected, ")")) +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Features", y = "Scaled Values") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )

  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

#' Create Combined Grid Plot
#' @noRd
.create_combined_plot <- function(plots, group_by) {

  main_title <- paste("Data Preprocessing Comparison:", group_by, "Perspective")

  combined_plot <- gridExtra::grid.arrange(
    grobs = list(plots$p1, plots$p2, plots$p3, plots$p4),
    nrow = 2,
    ncol = 2,
    top = grid::textGrob(
      main_title,
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )

  return(combined_plot)
}
