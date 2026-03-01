#' Perform Random Forest Classification on Preprocessed Metabolomics Data
#'
#' @param prepped_data Output from perform_PreprocessingPeakData function
#' @param outcome_var Character string specifying the outcome variable column name in metadata (default: "Group")
#' @param features Character vector of feature/metabolite names to include in analysis.
#'        If NULL, uses all available features (default: NULL)
#' @param data_type Character string specifying which data matrix to use: "NONPLS" or "PLS" (default: "NONPLS")
#' @param use_merged Logical indicating whether to use merged replicates if available (default: TRUE)
#' @param ntree Integer specifying number of trees to build (default: 1000)
#' @param mtry Integer or NULL specifying number of variables to consider at each split.
#'        If NULL, will optimize automatically (default: NULL)
#' @param optimize_mtry Logical indicating whether to optimize mtry parameter (default: TRUE)
#' @param mtry_range Integer vector specifying range of mtry values to test (default: 1:10)
#' @param seed Integer for reproducibility (default: 123)
#' @param proximity Logical indicating whether to compute proximity matrix (default: TRUE)
#'
#' @return List containing:
#'   - model: The final randomForest model object
#'   - oob_error: OOB error rate
#'   - confusion_matrix: Confusion matrix
#'   - error_trajectory: Data frame of error rates across trees
#'   - optimal_mtry: Optimal mtry value (if optimized)
#'   - mtry_optimization: Vector of OOB errors for different mtry values (if optimized)
#'   - metadata: Metadata used for the model
#'   - data_matrix: Data matrix used for the model
#'   - outcome_var: Name of outcome variable
#'   - features_used: Character vector of features used in the model
#'
#' @examples
#' rf_results <- perform_RandomForest(df_prepped, outcome_var = "Group")
#' rf_results <- perform_RandomForest(df_prepped, outcome_var = "Group2",
#'                                     data_type = "PLS", ntree = 500)
#' # Use specific features only
#' rf_results <- perform_RandomForest(df_prepped,
#'                                     features = c("Met001", "Met002", "Met010"))
#' @export
perform_RandomForest <- function(prepped_data,
                                 outcome_var = "Group",
                                 features = NULL,
                                 data_type = "NONPLS",
                                 use_merged = TRUE,
                                 ntree = 1000,
                                 mtry = NULL,
                                 optimize_mtry = TRUE,
                                 mtry_range = 1:10,
                                 seed = 123,
                                 proximity = TRUE) {

  # Load required library
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required but not installed.")
  }

  # Validate inputs
  if (!inherits(prepped_data, "perform_PreprocessingPeakData")) {
    stop("Input must be output from perform_PreprocessingPeakData function")
  }

  if (!data_type %in% c("NONPLS", "PLS")) {
    stop("data_type must be either 'NONPLS' or 'PLS'")
  }

  # Determine which data and metadata to use
  if (use_merged && !is.null(prepped_data$data_scaledNONPLS_merged)) {
    # Use merged data if available
    data_matrix <- if (data_type == "NONPLS") {
      prepped_data$data_scaledNONPLS_merged
    } else {
      prepped_data$data_scaledPLS_merged
    }
    metadata <- prepped_data$Metadata_merged
    message("Using merged replicate data")
  } else {
    # Use non-merged data
    if (data_type == "NONPLS") {
      # Check for variance filtered version first
      if (!is.null(prepped_data$data_scaledNONPLS_varFiltered)) {
        data_matrix <- prepped_data$data_scaledNONPLS_varFiltered
      } else {
        data_matrix <- prepped_data$data_scaledNONPLS
      }
    } else {
      if (!is.null(prepped_data$data_scaledPLS_varFiltered)) {
        data_matrix <- prepped_data$data_scaledPLS_varFiltered
      } else {
        data_matrix <- prepped_data$data_scaledPLS
      }
    }
    metadata <- prepped_data$Metadata
    message("Using non-merged data")
  }

  # Check if outcome variable exists in metadata
  if (!outcome_var %in% colnames(metadata)) {
    stop(paste0("Outcome variable '", outcome_var, "' not found in metadata. ",
                "Available columns: ", paste(colnames(metadata), collapse = ", ")))
  }

  # Remove QC samples before analysis
  if ("Group_" %in% colnames(metadata)) {
    qc_indices <- which(metadata$Group_ == "QC")
    if (length(qc_indices) > 0) {
      message(paste0("Removing ", length(qc_indices), " QC samples from analysis"))
      data_matrix <- data_matrix[-qc_indices, , drop = FALSE]
      metadata <- metadata[-qc_indices, , drop = FALSE]
    }
  }

  # Filter features if specified
  if (!is.null(features)) {
    # Check which features are available
    available_features <- colnames(data_matrix)
    missing_features <- setdiff(features, available_features)

    if (length(missing_features) > 0) {
      warning(paste0("The following features were not found and will be excluded: ",
                     paste(missing_features, collapse = ", ")))
    }

    features_to_use <- intersect(features, available_features)

    if (length(features_to_use) == 0) {
      stop("None of the specified features were found in the data matrix.")
    }

    message(paste0("Using ", length(features_to_use), " out of ",
                   length(features), " requested features"))

    # Subset data matrix to selected features
    data_matrix <- data_matrix[, features_to_use, drop = FALSE]
  } else {
    features_to_use <- colnames(data_matrix)
  }

  # Prepare the data for modeling
  outcome <- as.factor(metadata[[outcome_var]])

  # Check for missing values in outcome
  if (any(is.na(outcome))) {
    warning("Missing values found in outcome variable. Removing those samples.")
    valid_idx <- !is.na(outcome)
    data_matrix <- data_matrix[valid_idx, ]
    outcome <- outcome[valid_idx]
    metadata <- metadata[valid_idx, ]
  }

  # Combine data for modeling
  model_data <- data.frame(outcome = outcome, data_matrix)

  # Set seed for reproducibility
  set.seed(seed)

  # Optimize mtry if requested
  optimal_mtry_value <- NULL
  mtry_optimization_results <- NULL

  if (optimize_mtry && is.null(mtry)) {
    message("Optimizing mtry parameter...")
    n_features <- ncol(data_matrix)

    # Ensure mtry_range doesn't exceed number of features
    mtry_test_range <- mtry_range[mtry_range <= n_features]

    if (length(mtry_test_range) == 0) {
      warning(paste0("All values in mtry_range exceed the number of features (",
                     n_features, "). Using default mtry."))
      optimal_mtry_value <- floor(sqrt(n_features))
    } else {
      oob.values <- vector(length = length(mtry_test_range))
      for (i in seq_along(mtry_test_range)) {
        temp_model <- randomForest::randomForest(outcome ~ .,
                                                 data = model_data,
                                                 mtry = mtry_test_range[i],
                                                 ntree = ntree)
        oob.values[i] <- temp_model$err.rate[nrow(temp_model$err.rate), 1]
      }

      names(oob.values) <- paste0("mtry_", mtry_test_range)
      mtry_optimization_results <- oob.values
      optimal_mtry_value <- mtry_test_range[which.min(oob.values)]
      message(paste0("Optimal mtry: ", optimal_mtry_value,
                     " (OOB error: ", round(min(oob.values) * 100, 2), "%)"))
    }
  } else if (!is.null(mtry)) {
    optimal_mtry_value <- mtry
    message(paste0("Using specified mtry: ", mtry))
  } else {
    # Use default (sqrt of features for classification)
    optimal_mtry_value <- floor(sqrt(ncol(data_matrix)))
    message(paste0("Using default mtry: ", optimal_mtry_value))
  }

  # Build final model
  message("Building Random Forest model...")
  final_model <- randomForest::randomForest(outcome ~ .,
                                            data = model_data,
                                            ntree = ntree,
                                            mtry = optimal_mtry_value,
                                            proximity = proximity,
                                            importance = TRUE)

  # Extract error trajectory
  error_trajectory <- data.frame(
    Trees = rep(1:nrow(final_model$err.rate), times = (ncol(final_model$err.rate))),
    Type = rep(c("OOB", levels(outcome)), each = nrow(final_model$err.rate)),
    Error = c(final_model$err.rate)
  )

  # Create output list
  results <- list(
    model = final_model,
    oob_error = final_model$err.rate[nrow(final_model$err.rate), "OOB"],
    confusion_matrix = final_model$confusion,
    error_trajectory = error_trajectory,
    optimal_mtry = optimal_mtry_value,
    mtry_optimization = mtry_optimization_results,
    metadata = metadata,
    data_matrix = data_matrix,
    outcome_var = outcome_var,
    features_used = features_to_use,
    data_type = data_type,
    parameters = list(
      ntree = ntree,
      seed = seed,
      proximity = proximity,
      use_merged = use_merged
    )
  )

  class(results) <- c("perform_RandomForest", "list")

  message("Random Forest analysis complete!")
  message(paste0("OOB error rate: ", round(results$oob_error * 100, 2), "%"))

  return(results)
}


#' Create Plots from Random Forest Results
#'
#' @param rf_results Output from perform_RandomForest function
#' @param plot_type Character string specifying which plot to create:
#'                  "mds", "error", or "mtry" (default: "mds")
#' @param color_by Character string specifying metadata column to color points by
#'                 (for MDS plot, default: uses outcome_var from RF)
#' @param point_size Numeric specifying size of points (for MDS plot, default: 3)
#' @param add_labels Logical indicating whether to add sample labels (for MDS plot, default: FALSE)
#' @param label_size Numeric specifying size of labels (for MDS plot, default: 3)
#' @param line_size Numeric specifying line width (for error plot, default: 1)
#'
#' @return ggplot2 object
#'
#' @examples
#' rf_results <- perform_RandomForest(df_prepped, outcome_var = "Group")
#'
#' # MDS plot
#' mds_plot <- plot_RandomForest(rf_results, plot_type = "mds")
#' print(mds_plot)
#'
#' # MDS with labels and different color
#' mds_plot2 <- plot_RandomForest(rf_results, plot_type = "mds",
#'                                 color_by = "Group2", add_labels = TRUE)
#'
#' # Error trajectory plot
#' error_plot <- plot_RandomForest(rf_results, plot_type = "error")
#' print(error_plot)
#'
#' # mtry optimization plot
#' mtry_plot <- plot_RandomForest(rf_results, plot_type = "mtry")
#' print(mtry_plot)
#' @export
plot_RandomForest <- function(rf_results,
                              plot_type = "mds",
                              color_by = NULL,
                              point_size = 3,
                              add_labels = FALSE,
                              label_size = 3,
                              line_size = 1) {

  # Load required library
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }

  # Validate input
  if (!inherits(rf_results, "perform_RandomForest")) {
    stop("Input must be output from perform_RandomForest function")
  }

  # Validate plot_type
  if (!plot_type %in% c("mds", "error", "mtry")) {
    stop("plot_type must be one of: 'mds', 'error', or 'mtry'")
  }

  # Create the requested plot
  if (plot_type == "mds") {
    p <- .plot_rf_mds(rf_results, color_by, point_size, add_labels, label_size)
  } else if (plot_type == "error") {
    p <- .plot_rf_error(rf_results, line_size)
  } else if (plot_type == "mtry") {
    p <- .plot_rf_mtry(rf_results)
  }

  return(p)
}


# Internal function for MDS plot
.plot_rf_mds <- function(rf_results, color_by, point_size, add_labels, label_size) {

  if (is.null(rf_results$model$proximity)) {
    stop("Random Forest model does not contain proximity matrix. ",
         "Re-run perform_RandomForest with proximity = TRUE")
  }

  # Determine color variable
  if (is.null(color_by)) {
    color_by <- rf_results$outcome_var
  }

  if (!color_by %in% colnames(rf_results$metadata)) {
    stop(paste0("Color variable '", color_by, "' not found in metadata. ",
                "Available columns: ", paste(colnames(rf_results$metadata), collapse = ", ")))
  }

  # Convert proximity to distance matrix
  distance_matrix <- as.dist(1 - rf_results$model$proximity)

  # Perform MDS
  mds_result <- cmdscale(distance_matrix, eig = TRUE, x.ret = TRUE)

  # Calculate percentage of variation explained
  mds_var_per <- round(mds_result$eig / sum(mds_result$eig) * 100, 1)

  # Create data frame for plotting
  mds_values <- mds_result$points
  mds_data <- data.frame(
    Sample = rownames(mds_values),
    X = mds_values[, 1],
    Y = mds_values[, 2],
    ColorVar = rf_results$metadata[[color_by]]
  )

  # Create plot
  p <- ggplot2::ggplot(data = mds_data, ggplot2::aes(x = X, y = Y, color = ColorVar))

  if (add_labels) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), size = label_size)
  } else {
    p <- p + ggplot2::geom_point(size = point_size)
  }

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste0("MDS1 - ", mds_var_per[1], "%"),
      y = paste0("MDS2 - ", mds_var_per[2], "%"),
      color = color_by,
      title = "MDS Plot using (1 - Random Forest Proximities)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  return(p)
}


# Internal function for error trajectory plot
.plot_rf_error <- function(rf_results, line_size) {

  p <- ggplot2::ggplot(data = rf_results$error_trajectory,
                       ggplot2::aes(x = Trees, y = Error, color = Type)) +
    ggplot2::geom_line(size = line_size) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Number of Trees",
      y = "Error Rate",
      color = "Type",
      title = "Random Forest Out-of-Bag Error Rate"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  return(p)
}


# Internal function for mtry optimization plot
.plot_rf_mtry <- function(rf_results) {

  if (is.null(rf_results$mtry_optimization)) {
    stop("No mtry optimization results found. Re-run perform_RandomForest with optimize_mtry = TRUE")
  }

  # Prepare data
  mtry_data <- data.frame(
    mtry = as.numeric(gsub("mtry_", "", names(rf_results$mtry_optimization))),
    OOB_Error = rf_results$mtry_optimization
  )

  # Create plot
  p <- ggplot2::ggplot(data = mtry_data, ggplot2::aes(x = mtry, y = OOB_Error)) +
    ggplot2::geom_line(color = "blue", size = 1) +
    ggplot2::geom_point(color = "blue", size = 3) +
    ggplot2::geom_vline(xintercept = rf_results$optimal_mtry,
                        linetype = "dashed", color = "red", size = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "mtry (Number of Variables at Each Split)",
      y = "Out-of-Bag Error Rate",
      title = "mtry Parameter Optimization"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::annotate("text",
                      x = rf_results$optimal_mtry,
                      y = max(mtry_data$OOB_Error),
                      label = paste0("Optimal mtry = ", rf_results$optimal_mtry),
                      hjust = -0.1, color = "red")

  return(p)
}
