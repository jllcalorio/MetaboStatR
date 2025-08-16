#' Perform Principal Component Analysis (PCA)
#'
#' @description
#' This function performs PCA and generates multiple scores plots at once based on a list of
#' principal component combinations, plus an optional scree plot. Each combination in the list
#' will produce a separate scores plot.
#'
#' @param data List. This list must be a result from the `perform_PreprocessingPeakData` function.
#' @param scorePC List. A list of 2-element numeric vectors, where each vector specifies
#'   the principal components to plot. For example: list(c(1,2), c(1,3), c(2,3))
#'   will generate 3 plots: PC1 vs PC2, PC1 vs PC3, and PC2 vs PC3.
#' @param includeQC Boolean. If `TRUE`, includes QC (Quality Control) samples in the analysis and plots.
#'   If `FALSE`, uses only biological samples (BS). Defaults to `FALSE`.
#' @param arrangeLevels Vector. Determines how the groups will be arranged.
#'   The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which
#'   sorts the groups in alphabetical order.
#' @param scoresEllipse Boolean. If `TRUE` (default), adds an ellipse in the scores plot.
#' @param scoresTitle String or Vector. The scores plot title(s). Can be a single string
#'   (applied to all plots) or a vector of strings (one for each plot). If `NULL`,
#'   defaults to "PCA Scores Plot (PC\{i\} vs PC\{j\})".
#' @param scoresLegend String. The title in the legend section of the scores plot.
#'   Defaults to `NULL` which means no legend title.
#' @param includeScree Boolean. If `TRUE`, generates a scree plot in addition to scores plots.
#'   Defaults to `TRUE`.
#' @param screeWhat Character. What to plot in the scree plot: "variance" or "eigenvalue".
#'   Defaults to "variance".
#' @param screeType Character. Type of scree plot: "bar", "line", or "both".
#'   Defaults to "both".
#' @param screeTitle Character. Title for the scree plot. Defaults to "Scree Plot".
#' @param screeMaxPC Integer. Maximum number of PCs to show in scree plot.
#'   If NULL, shows all available PCs up to 20. Defaults to `NULL`.
#' @param screeShowValues Boolean. If `TRUE`, shows variance/eigenvalue values on top of bars
#'   and removes y-axis labels. Defaults to `TRUE`.
#' @param screeShowCumulative Boolean. If `TRUE`, shows cumulative variance explained
#'   on top of the scree plot. Only works when screeWhat="variance". Defaults to `FALSE`.
#'
#' @returns Returns a list containing:
#'   \itemize{
#'     \item plots: A list of ggplot objects, one for each PC combination
#'     \item scree_plot: A ggplot object for the scree plot (if includeScree is TRUE)
#'     \item plot_info: A data frame with information about each plot (PC combinations, variance explained)
#'     \item pca_results: The PCA results object
#'     \item data_used: The data matrix used for PCA
#'     \item variance_explained: Vector of variance explained by each PC
#'     \item eigenvalues: Vector of eigenvalues for each PC
#'   }
#' @export
#'
#' @importFrom stats complete.cases
#'
#' @examples
#' \dontrun{
#' # Generate 3 different scores plots plus scree plot
#' multi_plots <- perform_PCA(
#'   data = data_from_perform_PreprocessingPeakData_function,
#'   scorePC = list(c(1,2), c(1,3), c(2,3)),
#'   includeQC = FALSE,
#'   scoresEllipse = TRUE,
#'   includeScree = TRUE,
#'   screeWhat = "variance",
#'   screeType = "both"
#' )
#'
#' # Access individual plots
#' plot1 <- multi_plots$plots[[1]]  # PC1 vs PC2
#' plot2 <- multi_plots$plots[[2]]  # PC1 vs PC3
#' scree <- multi_plots$scree_plot  # Scree plot
#'
#' # Display all plots
#' display_AllPlots(multi_plots, include_scree = TRUE)
#' }
perform_PCA <- function(
    data,
    scorePC = list(c(1, 2), c(1, 3), c(2, 3)),
    includeQC = FALSE,
    arrangeLevels = NULL,
    scoresEllipse = TRUE,
    scoresTitle = NULL,
    scoresLegend = NULL,
    includeScree = TRUE,
    screeWhat = c("variance", "eigenvalue")[1],
    screeType = c("bar", "line", "both")[3],
    screeTitle = "Scree Plot",
    screeMaxPC = NULL,
    screeShowValues = TRUE,
    screeShowCumulative = FALSE
) {

  # ============================================================================
  # PARAMETER VALIDATION
  # ============================================================================

  # Check if data is from the correct function
  if (data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("The parameter 'data' must be from the 'perform_PreprocessingPeakData' function.")
  }

  # Check scorePC parameter
  if (!is.list(scorePC)) {
    stop("scorePC: Must be a list of numeric vectors. Example: list(c(1,2), c(1,3), c(2,3))")
  }

  if (length(scorePC) == 0) {
    stop("scorePC: List cannot be empty. Provide at least one PC combination.")
  }

  # Validate each element in scorePC list
  for (i in seq_along(scorePC)) {
    pc_combo <- scorePC[[i]]

    # Check if each element is numeric
    if (!is.numeric(pc_combo)) {
      stop(paste0("scorePC element ", i, ": Must be numeric. Found: ", class(pc_combo)[1]))
    }

    # Check if each element has exactly 2 values
    if (length(pc_combo) != 2) {
      stop(paste0("scorePC element ", i, ": Must contain exactly 2 values. Found: ", length(pc_combo), " values"))
    }

    # Check if values are positive integers
    if (any(pc_combo <= 0) || any(pc_combo != floor(pc_combo))) {
      stop(paste0("scorePC element ", i, ": Values must be positive integers. Found: c(",
                  paste(pc_combo, collapse = ", "), ")"))
    }

    # Check for duplicate values within the same combination
    if (pc_combo[1] == pc_combo[2]) {
      stop(paste0("scorePC element ", i, ": Cannot plot the same PC against itself. Found: c(",
                  paste(pc_combo, collapse = ", "), ")"))
    }
  }

  # Check boolean parameters
  if (!is.logical(includeQC) || length(includeQC) != 1) {
    stop("includeQC: Must be a single logical value (TRUE or FALSE).")
  }

  if (!is.logical(scoresEllipse) || length(scoresEllipse) != 1) {
    stop("scoresEllipse: Must be a single logical value (TRUE or FALSE).")
  }

  if (!is.logical(includeScree) || length(includeScree) != 1) {
    stop("includeScree: Must be a single logical value (TRUE or FALSE).")
  }

  # Check arrangeLevels
  if (!is.null(arrangeLevels)) {
    if (!is.vector(arrangeLevels) || !is.character(arrangeLevels)) {
      stop("arrangeLevels: Must be NULL or a character vector of group names.")
    }
  }

  # Check scoresTitle
  if (!is.null(scoresTitle)) {
    if (!is.character(scoresTitle)) {
      stop("scoresTitle: Must be NULL, a single character string, or a character vector.")
    }
    if (length(scoresTitle) > 1 && length(scoresTitle) != length(scorePC)) {
      stop(paste0("scoresTitle: If providing multiple titles, must provide exactly one for each plot. ",
                  "Expected: ", length(scorePC), " titles, Found: ", length(scoresTitle), " titles"))
    }
  }

  # Check scoresLegend
  if (!is.null(scoresLegend)) {
    if (!is.character(scoresLegend) || length(scoresLegend) != 1) {
      stop("scoresLegend: Must be NULL or a single character string.")
    }
  }

  # Check scree plot parameters
  if (!is.character(screeWhat) || length(screeWhat) != 1 || !screeWhat %in% c("variance", "eigenvalue")) {
    stop("screeWhat: Must be either 'variance' or 'eigenvalue'.")
  }

  if (!is.character(screeType) || length(screeType) != 1 || !screeType %in% c("bar", "line", "both")) {
    stop("screeType: Must be 'bar', 'line', or 'both'.")
  }

  if (!is.character(screeTitle) || length(screeTitle) != 1) {
    stop("screeTitle: Must be a single character string.")
  }

  if (!is.null(screeMaxPC)) {
    if (!is.numeric(screeMaxPC) || length(screeMaxPC) != 1 || screeMaxPC <= 0 || screeMaxPC != floor(screeMaxPC)) {
      stop("screeMaxPC: Must be NULL or a single positive integer.")
    }
  }

  if (!is.logical(screeShowValues) || length(screeShowValues) != 1) {
    stop("screeShowValues: Must be a single logical value (TRUE or FALSE).")
  }

  if (!is.logical(screeShowCumulative) || length(screeShowCumulative) != 1) {
    stop("screeShowCumulative: Must be a single logical value (TRUE or FALSE).")
  }

  # Check cumulative variance parameter compatibility
  if (screeShowCumulative && screeWhat != "variance") {
    warning("screeShowCumulative only works with screeWhat='variance'. Setting screeShowCumulative to FALSE.")
    screeShowCumulative <- FALSE
  }

  # ============================================================================
  # DATA PREPARATION AND PCA
  # ============================================================================

  # Determine QC and biological sample indices
  qc_indices <- data$Metadata$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices <- !qc_indices

  # Determine which samples to use based on includeQC parameter
  if (includeQC) {
    sample_indices <- rep(TRUE, nrow(data$data_scaledPCA_rsdFiltered_varFiltered))
    sample_type <- "All samples (QC + Biological)"
  } else {
    sample_indices <- non_qc_indices
    sample_type <- "Biological samples only"
  }

  # Validate that we have samples to analyze
  if (sum(sample_indices) == 0) {
    stop("No samples available for analysis with the current settings.")
  }

  # Select appropriate filtered data
  if (is.null(data$data_scaledPCA_rsdFiltered_varFiltered) ||
      length(data$data_scaledPCA_rsdFiltered_varFiltered) == 0) {
    stop("The preprocessed data does not exist or does not have any data on it.")
  } else {
    data_pca <- data$data_scaledPCA_rsdFiltered_varFiltered[sample_indices, , drop = FALSE]
  }

  # Check if we have enough samples and variables for PCA
  if (nrow(data_pca) < 2) {
    stop("Need at least 2 samples to perform PCA.")
  }

  if (ncol(data_pca) < 2) {
    stop("Need at least 2 variables to perform PCA.")
  }

  # Check for missing values
  if (any(is.na(data_pca))) {
    warning("Missing values detected in data. Consider additional preprocessing.")
    # Remove columns with any missing values for PCA
    data_pca <- data_pca[, complete.cases(t(data_pca)), drop = FALSE]
    if (ncol(data_pca) < 2) {
      stop("After removing missing values, insufficient variables remain for PCA.")
    }
  }

  # Perform PCA with error handling
  pca_res <- tryCatch({
    stats::prcomp(data_pca, center = FALSE, scale. = FALSE)
  }, error = function(e) {
    stop(paste("PCA failed:", e$message))
  })

  # Check if requested PCs exist
  max_pc <- max(unlist(scorePC))
  if (max_pc > ncol(pca_res$x)) {
    stop(paste0("Requested PC", max_pc, " but only ", ncol(pca_res$x),
                " principal components are available."))
  }

  # Calculate variance explained and eigenvalues
  eigenvalues <- pca_res$sdev^2
  variance_explained <- round(100 * (eigenvalues / sum(eigenvalues)), 2)

  # ============================================================================
  # GENERATE SCREE PLOT (if requested)
  # ============================================================================

  scree_plot <- NULL
  if (includeScree) {
    # Determine number of PCs to show in scree plot
    n_pcs <- ncol(pca_res$x)
    max_show <- if (is.null(screeMaxPC)) min(n_pcs, 20) else min(screeMaxPC, n_pcs)

    # Prepare data for scree plot
    if (screeWhat == "variance") {
      y_values <- variance_explained[1:max_show]
      y_label <- "Variance Explained (%)"
      if (screeShowCumulative) {
        y_values <- variance_explained[1:max_show]
        y_label <- "Variance Explained (%)\nCumulative Explained (%)"
      }
    } else {
      y_values <- eigenvalues[1:max_show]
      y_label <- "Eigenvalue"
    }

    scree_data <- data.frame(
      PC = factor(paste0("PC", 1:max_show), levels = paste0("PC", 1:max_show)),
      PC_num = 1:max_show,
      Value = y_values,
      CumValue = if (screeWhat == "variance") cumsum(y_values) else NA,
      stringsAsFactors = FALSE
    )

    # Create scree plot based on type
    scree_plot <- ggplot2::ggplot(scree_data, ggplot2::aes(x = PC_num, y = Value))

    if (screeType == "bar") {
      scree_plot <- scree_plot +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
        ggplot2::scale_x_continuous(breaks = 1:max_show, labels = paste0("PC", 1:max_show))
    } else if (screeType == "line") {
      scree_plot <- scree_plot +
        ggplot2::geom_line(color = "steelblue", size = 1, group = 1) +
        ggplot2::geom_point(color = "steelblue", size = 3) +
        ggplot2::scale_x_continuous(breaks = 1:max_show, labels = paste0("PC", 1:max_show))
    } else { # "both"
      scree_plot <- scree_plot +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.5, width = 0.6) +
        ggplot2::geom_line(color = "darkblue", size = 1, group = 1) +
        ggplot2::geom_point(color = "darkblue", size = 3) +
        ggplot2::scale_x_continuous(breaks = 1:max_show, labels = paste0("PC", 1:max_show))
    }

    # Add value labels on bars/points if requested
    if (screeShowValues) {
      # Format the values for display
      value_labels <- if (screeWhat == "variance") {
        paste0(round(scree_data$Value, 1), "%")
      } else {
        round(scree_data$Value, 2)
      }

      # Position labels to the right of bars for bar/both plots, on top for line plots
      if (screeType == "line") {
        label_hjust <- 0.5
        label_vjust <- -0.5
      } else {
        label_hjust <- -0.2  # Position to the right
        label_vjust <- -0.5   # Center vertically
      }

      scree_plot <- scree_plot +
        ggplot2::geom_text(ggplot2::aes(label = value_labels),
                           vjust = label_vjust,
                           hjust = label_hjust,
                           size = 3,
                           color = "black")
    }

    # Add cumulative variance labels if requested and applicable
    if (screeShowCumulative && screeWhat == "variance") {
      cum_labels <- paste0(round(scree_data$CumValue, 1), "%")

      # Position cumulative labels higher and more to the right than individual values
      if (screeType == "line") {
        cum_hjust <- 0.5
        cum_vjust <- if (screeShowValues) -1.2 else -0.5
      } else {
        cum_hjust <- -0.2  # Same horizontal position as individual values
        cum_vjust <- if (screeShowValues) -1.8 else -1.8  # Higher than individual values
      }

      scree_plot <- scree_plot +
        ggplot2::geom_text(ggplot2::aes(label = cum_labels),
                           vjust = cum_vjust,
                           hjust = cum_hjust,
                           size = 2.8,
                           color = "darkred",
                           fontface = "italic")
    }

    # Determine y-axis label and whether to show y-axis values
    if (screeShowValues) {
      # Remove y-axis values but keep the title
      y_axis_text <- ggplot2::element_blank()
      y_axis_ticks <- ggplot2::element_blank()
      y_title <- y_label  # Keep the y-axis title
    } else {
      y_axis_text <- ggplot2::element_text(size = 10)
      y_axis_ticks <- ggplot2::element_line()
      y_title <- y_label
    }

    scree_plot <- scree_plot +
      ggplot2::labs(x = "Principal Component",
                    y = y_title,
                    title = screeTitle) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.title = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.5, size = 10),
        axis.text.y = y_axis_text,
        axis.ticks.y = y_axis_ticks,
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  # ============================================================================
  # GENERATE MULTIPLE SCORES PLOTS
  # ============================================================================

  plots_list <- list()
  plot_info <- data.frame(
    plot_number = integer(),
    PC1 = integer(),
    PC2 = integer(),
    PC1_variance = numeric(),
    PC2_variance = numeric(),
    title = character(),
    stringsAsFactors = FALSE
  )

  # Prepare scores data based on sample selection
  metadata_subset <- data$Metadata[sample_indices, , drop = FALSE]

  if (includeQC) {
    qc_subset <- qc_indices[sample_indices]
    combined_groups <- ifelse(qc_subset,
                              paste0("QC_", metadata_subset$Group),
                              as.character(metadata_subset$Group))
    base_scores_df <- data.frame(
      Group = combined_groups,
      Batch = metadata_subset$Batch,
      Label = rownames(data_pca),
      SampleType = ifelse(qc_subset, "QC", "Biological"),
      stringsAsFactors = FALSE
    )
  } else {
    base_scores_df <- data.frame(
      Group = metadata_subset$Group,
      Batch = metadata_subset$Batch,
      Label = rownames(data_pca),
      SampleType = "Biological",
      stringsAsFactors = FALSE
    )
  }

  # Apply arrangeLevels if specified
  if (!is.null(arrangeLevels)) {
    # Check if all specified levels exist in the data
    existing_groups <- unique(base_scores_df$Group)
    missing_groups <- arrangeLevels[!arrangeLevels %in% existing_groups]
    if (length(missing_groups) > 0) {
      warning(paste("The following groups in arrangeLevels were not found in data:",
                    paste(missing_groups, collapse = ", ")))
    }
    # Only use levels that exist in the data
    valid_levels <- arrangeLevels[arrangeLevels %in% existing_groups]
    # Add any groups that weren't specified in arrangeLevels
    remaining_groups <- existing_groups[!existing_groups %in% valid_levels]
    final_levels <- c(valid_levels, sort(remaining_groups))
    base_scores_df$Group <- factor(base_scores_df$Group, levels = final_levels)
  } else {
    base_scores_df$Group <- factor(base_scores_df$Group)
  }

  # Generate each plot
  for (i in seq_along(scorePC)) {
    pc_combo <- scorePC[[i]]
    pc1 <- pc_combo[1]
    pc2 <- pc_combo[2]

    # Create scores dataframe for this PC combination
    df_scores <- base_scores_df
    df_scores$PC1 <- pca_res$x[, pc1]
    df_scores$PC2 <- pca_res$x[, pc2]

    # Create axis labels
    x_label <- paste0("PC", pc1, " (", variance_explained[pc1], "%)")
    y_label <- paste0("PC", pc2, " (", variance_explained[pc2], "%)")

    # Determine plot title
    if (is.null(scoresTitle)) {
      plot_title <- paste0("PCA Scores Plot (PC", pc1, " vs PC", pc2, ")")
    } else if (length(scoresTitle) == 1) {
      plot_title <- scoresTitle
    } else {
      plot_title <- scoresTitle[i]
    }

    # Create the plot with optimized aesthetics
    scores_plot <- ggplot2::ggplot(df_scores,
                                   ggplot2::aes(x = PC1,
                                                y = PC2,
                                                color = Group,
                                                shape = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = Label),
                         vjust = -0.7,
                         hjust = 0.5,
                         size = 2.8,
                         alpha = 0.7)

    # Add ellipses if requested
    if (scoresEllipse && nlevels(df_scores$Group) > 1) {
      # Check if each group has enough points for ellipse
      group_counts <- table(df_scores$Group)
      groups_for_ellipse <- names(group_counts)[group_counts >= 3]

      if (length(groups_for_ellipse) > 0) {
        df_ellipse <- df_scores[df_scores$Group %in% groups_for_ellipse, ]
        scores_plot <- scores_plot +
          ggplot2::stat_ellipse(data = df_ellipse,
                                ggplot2::aes(group = Group, fill = Group),
                                level = 0.95,
                                geom = "polygon",
                                alpha = 0.2,
                                show.legend = FALSE)
      }
    }

    # Complete the plot with improved theme
    scores_plot <- scores_plot +
      ggplot2::labs(x = x_label,
                    y = y_label,
                    title = plot_title,
                    color = scoresLegend,
                    shape = scoresLegend,
                    fill = scoresLegend) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 11),
        legend.text = ggplot2::element_text(size = 10)
      )

    # Store the plot
    plots_list[[i]] <- scores_plot

    # Store plot information
    plot_info <- rbind(plot_info, data.frame(
      plot_number = i,
      PC1 = pc1,
      PC2 = pc2,
      PC1_variance = variance_explained[pc1],
      PC2_variance = variance_explained[pc2],
      title = plot_title,
      stringsAsFactors = FALSE
    ))
  }

  # Name the plots list
  names(plots_list) <- paste0("PC", plot_info$PC1, "_vs_PC", plot_info$PC2)

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  results <- list(
    plots = plots_list,
    scree_plot = scree_plot,
    plot_info = plot_info,
    pca_results = pca_res,
    data_used = data_pca,
    sample_type = sample_type,
    variance_explained = as.data.frame(variance_explained),
    eigenvalues = as.data.frame(eigenvalues),
    function_call = match.call(),
    scree_included = includeScree
  )

  class(results) <- c("multipleScoresPlots", "list")

  return(results)
}

# Enhanced helper function to print all plots
#' @export
print.multipleScoresPlots <- function(x, ...) {
  cat("Multiple PCA Scores Plots Object\n")
  cat("================================\n")
  cat("Number of scores plots:", length(x$plots), "\n")
  cat("Scree plot included:", ifelse(x$scree_included && !is.null(x$scree_plot), "Yes", "No"), "\n")
  cat("Sample type:", x$sample_type, "\n")
  cat("Total variance explained by all PCs:", round(sum(x$variance_explained), 1), "%\n")
  cat("PC combinations:\n")
  for (i in 1:nrow(x$plot_info)) {
    cat(sprintf("  Plot %d: PC%d vs PC%d (%.1f%% vs %.1f%% variance)\n",
                i, x$plot_info$PC1[i], x$plot_info$PC2[i],
                x$plot_info$PC1_variance[i], x$plot_info$PC2_variance[i]))
  }
  cat("\nUse $plots[[i]] to access individual scores plots\n")
  if (x$scree_included && !is.null(x$scree_plot)) {
    cat("Use $scree_plot to access the scree plot\n")
  }
  cat("Use $plot_info for detailed information about each plot\n")
  cat("Use display_AllPlots() to show all plots at once\n")
}

# Enhanced helper function to display all plots at once
#' Display All Plots from Multiple Scores Plots Object
#'
#' @param multi_plots_obj Object returned from perform_PCA()
#' @param arrange Logical. If TRUE, attempts to arrange plots in a grid using patchwork (if available).
#'   If FALSE or patchwork unavailable, prints plots individually. Default: TRUE.
#' @param include_scree Logical. Whether to include the scree plot in the display. Default: TRUE.
#' @param ncol Number of columns in the plot arrangement (default: 2)
#' @param nrow Number of rows in the plot arrangement (default: NULL, auto-calculated)
#' @export
display_AllPlots <- function(multi_plots_obj, arrange = TRUE, include_scree = TRUE, ncol = 2, nrow = NULL) {
  if (!inherits(multi_plots_obj, "multipleScoresPlots")) {
    stop("Object must be created by perform_PCA() function")
  }

  # Determine which plots to display
  plots_to_show <- multi_plots_obj$plots

  if (include_scree && multi_plots_obj$scree_included && !is.null(multi_plots_obj$scree_plot)) {
    plots_to_show <- c(list(scree = multi_plots_obj$scree_plot), plots_to_show)
  }

  # If arrange is TRUE and patchwork is available, use it
  if (arrange && requireNamespace("patchwork", quietly = TRUE)) {
    # Use patchwork to combine plots
    wrap_plots_fn <- getFromNamespace("wrap_plots", "patchwork")
    combined_plot <- wrap_plots_fn(plots_to_show, ncol = ncol, nrow = nrow)
    print(combined_plot)
  } else {
    # Print plots individually
    if (arrange) {
      cat("Note: 'patchwork' package not available. Displaying plots individually.\n")
      cat("To display plots in a grid, install patchwork: install.packages('patchwork')\n\n")
    }

    for (i in seq_along(plots_to_show)) {
      plot_name <- names(plots_to_show)[i]
      if (is.null(plot_name) || plot_name == "") {
        plot_name <- paste("Plot", i)
      }
      cat("Plot:", plot_name, "\n")
      print(plots_to_show[[i]])
      if (i < length(plots_to_show)) {
        cat("\n", paste(rep("-", 50), collapse = ""), "\n\n")
      }
    }
  }
}

# Additional utility function to extract scree plot data
#' Extract Scree Plot Data
#'
#' @param multi_plots_obj Object returned from perform_PCA()
#' @param what Character. "variance" or "eigenvalue" to specify what values to return
#' @param max_pc Integer. Maximum number of PCs to include. If NULL, includes all.
#' @return Data frame with PC information
#' @export
getScreeData <- function(multi_plots_obj, what = "variance", max_pc = NULL) {
  if (!inherits(multi_plots_obj, "multipleScoresPlots")) {
    stop("Object must be created by perform_PCA() function")
  }

  n_pcs <- length(multi_plots_obj$variance_explained)
  max_show <- if (is.null(max_pc)) n_pcs else min(max_pc, n_pcs)

  if (what == "variance") {
    values <- multi_plots_obj$variance_explained[1:max_show]
    value_name <- "Variance_Explained_Percent"
  } else if (what == "eigenvalue") {
    values <- multi_plots_obj$eigenvalues[1:max_show]
    value_name <- "Eigenvalue"
  } else {
    stop("'what' must be either 'variance' or 'eigenvalue'")
  }

  result <- data.frame(
    PC = paste0("PC", 1:max_show),
    PC_Number = 1:max_show,
    stringsAsFactors = FALSE
  )
  result[[value_name]] <- values

  if (what == "variance") {
    result$Cumulative_Variance <- cumsum(values)
  }

  return(result)
}
