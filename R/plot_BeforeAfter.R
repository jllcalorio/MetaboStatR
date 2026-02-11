#' Visualize Data Before and After Preprocessing
#'
#' @description
#' Creates comprehensive visualizations comparing data distributions and patterns
#' before and after preprocessing steps. Generates density plots and box plots
#' from either sample or feature perspectives to assess preprocessing effectiveness.
#' Supports both NON-PLS and PLS scaled data visualization with intelligent
#' random sampling for large datasets.
#'
#' @param data List. A preprocessing results object from \code{perform_PreprocessingPeakData}
#'   function containing original data, transformed data, and metadata.
#' @param scaled Character. Specifies which preprocessed data to visualize:
#'   \itemize{
#'     \item \code{"NONPLS"}: Uses NONPLS-scaled preprocessed data
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
#' # Feature perspective with NONPLS scaling
#' feature_plots <- plot_BeforeAfter(
#'   data = preprocessing_results,
#'   scaled = "NONPLS",
#'   group_by = "Feature",
#'   n_random_features = 20
#' )
#'
#' # Display combined plot
#' print(plots$plot_combined)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_density geom_boxplot ggtitle theme_minimal coord_flip element_text labs theme
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob gpar
#' @importFrom stats setNames
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{perform_PreprocessingPeakData}}

plot_BeforeAfter <- function(data,
                             scaled = "PLS",
                             group_by = "Sample",
                             n_random_samples = 30L,
                             n_random_features = 30L,
                             seed = 123L) {

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # ===== Input Validation =====
  if (!is.list(data)) stop("Input 'data' must be a list object.", call. = FALSE)

  if (is.null(data$FunctionOrigin) || data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("Data must originate from 'perform_PreprocessingPeakData' function.", call. = FALSE)
  }

  required_components <- c("data_no_NA", "Metadata", "data_scaledNONPLS_varFiltered")
  missing <- setdiff(required_components, names(data))
  if (length(missing) > 0) {
    stop("Missing required data components: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!scaled %in% c("NONPLS", "PLS")) stop("'scaled' must be either 'NONPLS' or 'PLS'.", call. = FALSE)
  if (!group_by %in% c("Sample", "Feature")) stop("'group_by' must be either 'Sample' or 'Feature'.", call. = FALSE)

  if (!is.null(n_random_samples) && (n_random_samples <= 0 || !is.numeric(n_random_samples))) {
    stop("'n_random_samples' must be a positive integer or NULL.", call. = FALSE)
  }

  if (!is.null(n_random_features) && (n_random_features <= 0 || !is.numeric(n_random_features))) {
    stop("'n_random_features' must be a positive integer or NULL.", call. = FALSE)
  }

  if (is.null(data$Metadata$Group_)) {
    stop("Metadata must contain 'Group_' column for QC identification.", call. = FALSE)
  }

  # Check required packages
  required_packages <- c("ggplot2", "dplyr", "tidyr", "tibble", "gridExtra", "grid")
  missing_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Required packages not installed: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall with: install.packages(c('", paste(missing_pkgs, collapse = "', '"), "'))",
         call. = FALSE)
  }

  # ===== Data Preparation =====
  non_qc_idx <- data$Metadata$Group_ != "QC"
  if (!any(non_qc_idx)) stop("No non-QC samples found in the data.", call. = FALSE)

  # Extract matrices using vectorized subsetting
  original <- data$data_no_NA[non_qc_idx, , drop = FALSE]
  transformed <- data$data_scaledNONPLS_varFiltered[non_qc_idx, , drop = FALSE]

  if (nrow(original) == 0 || ncol(original) == 0) {
    stop("Original data is empty after removing QC samples.", call. = FALSE)
  }

  # Ensure consistent naming
  if (is.null(rownames(original))) rownames(original) <- rownames(transformed)
  if (is.null(colnames(original))) colnames(original) <- colnames(transformed)

  # Match features between datasets
  common_features <- intersect(colnames(original), colnames(transformed))
  if (length(common_features) < ncol(transformed)) {
    original <- original[, common_features, drop = FALSE]
    transformed <- transformed[, common_features, drop = FALSE]
  }

  # ===== Convert to Long Format =====
  # Vectorized conversion avoiding loops
  orig_df <- as.data.frame(original)
  trans_df <- as.data.frame(transformed)

  orig_df$Sample <- rownames(original)
  trans_df$Sample <- rownames(transformed)

  # Use data.table-style melting for speed if available, otherwise tidyr
  orig_long <- tidyr::pivot_longer(orig_df, -Sample, names_to = "Feature", values_to = "Value")
  trans_long <- tidyr::pivot_longer(trans_df, -Sample, names_to = "Feature", values_to = "Value")

  orig_long$DataType <- "Before"
  trans_long$DataType <- "After"

  # ===== Generate Plots Based on Group By =====
  if (group_by == "Sample") {
    # Density plots
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

    # Box plots - vectorized sample selection
    common_samples <- intersect(unique(orig_long$Sample), unique(trans_long$Sample))
    n_samples <- length(common_samples)

    if (is.null(n_random_samples) || n_random_samples >= n_samples) {
      selected_samples <- common_samples
      n_selected <- n_samples
    } else {
      selected_samples <- sample(common_samples, min(n_random_samples, n_samples))
      n_selected <- length(selected_samples)
    }

    # Vectorized filtering
    orig_subset <- orig_long[orig_long$Sample %in% selected_samples, ]
    trans_subset <- trans_long[trans_long$Sample %in% selected_samples, ]

    p3 <- ggplot2::ggplot(orig_subset, ggplot2::aes(x = reorder(Sample, Value, FUN = median), y = Value)) +
      ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.7) +
      ggplot2::ggtitle(paste("Sample Profiles - Before (n =", n_selected, ")")) +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Samples", y = "Intensity Values") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

    p4 <- ggplot2::ggplot(trans_subset, ggplot2::aes(x = reorder(Sample, Value, FUN = median), y = Value)) +
      ggplot2::geom_boxplot(fill = "lightcoral", alpha = 0.7) +
      ggplot2::ggtitle(paste("Sample Profiles - After (n =", n_selected, ")")) +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Samples", y = "Scaled Values") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

  } else {  # Feature perspective
    # Vectorized feature selection
    common_features <- intersect(unique(orig_long$Feature), unique(trans_long$Feature))
    n_features <- length(common_features)

    if (is.null(n_random_features) || n_random_features >= n_features) {
      selected_features <- common_features
      n_selected <- n_features
    } else {
      selected_features <- sample(common_features, min(n_random_features, n_features))
      n_selected <- length(selected_features)
    }

    # Vectorized filtering
    orig_subset <- orig_long[orig_long$Feature %in% selected_features, ]
    trans_subset <- trans_long[trans_long$Feature %in% selected_features, ]

    # Density plots
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

    # Box plots
    p3 <- ggplot2::ggplot(orig_subset, ggplot2::aes(x = reorder(Feature, Value, FUN = median), y = Value)) +
      ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.7) +
      ggplot2::ggtitle(paste("Feature Profiles - Before (n =", n_selected, ")")) +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Features", y = "Intensity Values") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))

    p4 <- ggplot2::ggplot(trans_subset, ggplot2::aes(x = reorder(Feature, Value, FUN = median), y = Value)) +
      ggplot2::geom_boxplot(fill = "lightcoral", alpha = 0.7) +
      ggplot2::ggtitle(paste("Feature Profiles - After (n =", n_selected, ")")) +
      ggplot2::theme_minimal() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Features", y = "Scaled Values") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))
  }

  # ===== Create Combined Plot =====
  main_title <- paste("Data Preprocessing Comparison:", group_by, "Perspective")
  combined_plot <- gridExtra::grid.arrange(
    grobs = list(p1, p2, p3, p4),
    nrow = 2, ncol = 2,
    top = grid::textGrob(main_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )

  # ===== Return Results =====
  list(
    FunctionOrigin = "plot_BeforeAfter",
    data_before = original,
    data_after = transformed,
    plot_density_before = p1,
    plot_density_after = p2,
    plot_box_before = p3,
    plot_box_after = p4,
    plot_combined = combined_plot,
    processing_info = list(
      scaled = scaled,
      group_by = group_by,
      n_random_samples = n_random_samples,
      n_random_features = n_random_features,
      seed = seed,
      original_dims = dim(original),
      transformed_dims = dim(transformed),
      qc_samples_removed = sum(data$Metadata$Group_ == "QC"),
      timestamp = Sys.time()
    )
  )
}
