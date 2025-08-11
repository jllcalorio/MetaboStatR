#' Generate Multiple PCA Scores Plots
#'
#' @description
#' This function performs PCA generates multiple scores plots at once based on a list of
#' principal component combinations. Each combination in the list will produce
#' a separate scores plot.
#'
#' @param data List. This list must be a result from the `performPreprocessingPeakData` function.
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
#'   defaults to "PCA Scores Plot (PC{i} vs PC{j})".
#' @param scoresLegend String. The title in the legend section of the scores plot.
#'   Defaults to `NULL` which means no legend title.
#'
#' @returns Returns a list containing:
#'   \itemize{
#'     \item plots: A list of ggplot objects, one for each PC combination
#'     \item plot_info: A data frame with information about each plot (PC combinations, variance explained)
#'     \item pca_results: The PCA results object
#'     \item data_used: The data matrix used for PCA
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate 3 different scores plots
#' multi_plots <- generate_MultipleScoresPlots(
#'   data = data_from_performPreprocessingPeakData_function,
#'   scorePC = list(c(1,2), c(1,3), c(2,3)),
#'   includeQC = FALSE,
#'   scoresEllipse = TRUE
#' )
#'
#' # Access individual plots
#' plot1 <- multi_plots$plots[[1]]  # PC1 vs PC2
#' plot2 <- multi_plots$plots[[2]]  # PC1 vs PC3
#' plot3 <- multi_plots$plots[[3]]  # PC2 vs PC3
#'
#' # Display all plots
#' for(i in seq_along(multi_plots$plots)) {
#'   print(multi_plots$plots[[i]])
#' }
#' }
generate_MultipleScoresPlots <- function(
    data,
    scorePC = list(c(1, 2), c(1, 3), c(2, 3)),
    includeQC = FALSE,
    arrangeLevels = NULL,
    scoresEllipse = TRUE,
    scoresTitle = NULL,
    scoresLegend = NULL
) {

  # ============================================================================
  # PARAMETER VALIDATION
  # ============================================================================

  # Check if data is from the correct function
  if (!exists("FunctionOrigin", where = data) || data$FunctionOrigin != "performPreprocessingPeakData") {
    stop("The parameter 'data' must be from the 'performPreprocessingPeakData' function.")
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

  # Select appropriate filtered data
  if (is.null(data$data_scaledPCA_rsdFiltered_varFiltered) ||
      length(data$data_scaledPCA_rsdFiltered_varFiltered) == 0) {
    # data_pca <- data$data_preprocessed[sample_indices, ]
    stop("The preprocessed data does not exist or does not have any data on it.")
  } else {
    data_pca <- data$data_scaledPCA_rsdFiltered_varFiltered[sample_indices, ]
  }

  # Perform PCA
  pca_res <- stats::prcomp(data_pca, center = FALSE, scale. = FALSE)

  # Check if requested PCs exist
  max_pc <- max(unlist(scorePC))
  if (max_pc > ncol(pca_res$x)) {
    stop(paste0("Requested PC", max_pc, " but only ", ncol(pca_res$x),
                " principal components are available."))
  }

  # Calculate variance explained
  variance_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)

  # ============================================================================
  # GENERATE MULTIPLE PLOTS
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
  if (includeQC) {
    combined_groups <- ifelse(qc_indices[sample_indices],
                              paste0("QC_", data$Metadata$Group[sample_indices]),
                              as.character(data$Metadata$Group[sample_indices]))
    base_scores_df <- data.frame(
      Group = combined_groups,
      Batch = data$Metadata$Batch[sample_indices],
      Label = rownames(data_pca),
      SampleType = ifelse(qc_indices[sample_indices], "QC", "Biological"),
      stringsAsFactors = FALSE
    )
  } else {
    base_scores_df <- data.frame(
      Group = data$Metadata$Group[sample_indices],
      Batch = data$Metadata$Batch[sample_indices],
      Label = rownames(data_pca),
      SampleType = "Biological",
      stringsAsFactors = FALSE
    )
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

    # Create the plot
    scores_plot <- ggplot2::ggplot(df_scores,
                                   ggplot2::aes(x = PC1,
                                                y = PC2,
                                                color = Group,
                                                shape = Group)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_text(ggplot2::aes(label = Label),
                         vjust = -0.5,
                         hjust = 0.5,
                         size = 3)

    # Add ellipses if requested
    if (scoresEllipse) {
      scores_plot <- scores_plot +
        ggplot2::stat_ellipse(ggplot2::aes(group = Group, fill = Group),
                              level = 0.95,
                              geom = "polygon",
                              alpha = 0.3)
    }

    # Complete the plot
    scores_plot <- scores_plot +
      ggplot2::labs(x = x_label,
                    y = y_label,
                    title = plot_title,
                    color = scoresLegend,
                    shape = scoresLegend,
                    fill = scoresLegend) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(hjust = 0.5))

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
    plot_info = plot_info,
    pca_results = pca_res,
    data_used = data_pca,
    sample_type = sample_type,
    variance_explained = variance_explained,
    function_call = match.call()
  )

  class(results) <- c("multipleScoresPlots", "list")

  return(results)
}

# Helper function to print all plots
#' @export
print.multipleScoresPlots <- function(x, ...) {
  cat("Multiple PCA Scores Plots Object\n")
  cat("================================\n")
  cat("Number of plots:", length(x$plots), "\n")
  cat("Sample type:", x$sample_type, "\n")
  cat("PC combinations:\n")
  for (i in 1:nrow(x$plot_info)) {
    cat(sprintf("  Plot %d: PC%d vs PC%d (%.1f%% vs %.1f%% variance)\n",
                i, x$plot_info$PC1[i], x$plot_info$PC2[i],
                x$plot_info$PC1_variance[i], x$plot_info$PC2_variance[i]))
  }
  cat("\nUse $plots[[i]] to access individual plots\n")
  cat("Use $plot_info for detailed information about each plot\n")
}

# Helper function to display all plots at once
#' Display All Plots from Multiple Scores Plots Object
#'
#' @param multi_plots_obj Object returned from generate_MultipleScoresPlots()
#' @param arrange Logical. If TRUE, attempts to arrange plots in a grid using patchwork (if available).
#'   If FALSE or patchwork unavailable, prints plots individually. Default: TRUE.
#' @param ncol Number of columns in the plot arrangement (default: 2)
#' @param nrow Number of rows in the plot arrangement (default: NULL, auto-calculated)
#' @export
display_AllPlots <- function(multi_plots_obj, arrange = TRUE, ncol = 2, nrow = NULL) {
  if (!inherits(multi_plots_obj, "multipleScoresPlots")) {
    stop("Object must be created by generate_MultipleScoresPlots() function")
  }

  # If arrange is TRUE and patchwork is available, use it
  if (arrange && requireNamespace("patchwork", quietly = TRUE)) {
    # Use patchwork to combine plots (accessing functions properly)
    wrap_plots_fn <- getFromNamespace("wrap_plots", "patchwork")
    combined_plot <- wrap_plots_fn(multi_plots_obj$plots, ncol = ncol, nrow = nrow)
    print(combined_plot)
  } else {
    # Print plots individually
    if (arrange) {
      cat("Note: 'patchwork' package not available. Displaying plots individually.\n")
      cat("To display plots in a grid, install patchwork: install.packages('patchwork')\n\n")
    }

    for (i in seq_along(multi_plots_obj$plots)) {
      cat("Plot", i, ":", names(multi_plots_obj$plots)[i], "\n")
      print(multi_plots_obj$plots[[i]])
      if (i < length(multi_plots_obj$plots)) {
        cat("\n", paste(rep("-", 50), collapse = ""), "\n\n")
      }
    }
  }
}
