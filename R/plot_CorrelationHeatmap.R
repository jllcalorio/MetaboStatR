#' Plot Correlation and Hierarchical Clustering Heatmaps
#'
#' @description
#' Generates correlation or hierarchical clustering heatmaps with advanced filtering
#' capabilities. The function creates publication-ready visualizations by selecting
#' the most variable features/samples based on interquartile range (IQR) and applies
#' statistical significance masking for correlation plots. Automatically selects the
#' appropriate data source based on replicate merging settings.
#'
#' @param data A list object containing preprocessed data. Must include a component
#'   named \code{data_scaledPCA_varFiltered} or \code{data_scaledPCA_merged} (typically output from
#'   \code{perform_PreprocessingPeakData} function). The function automatically selects
#'   the appropriate data source based on \code{data$Parameters$auto_merge_replicates}.
#' @param method Character string specifying the correlation method. One of:
#'   \itemize{
#'     \item \code{"pearson"}: Pearson product-moment correlation (default)
#'     \item \code{"spearman"}: Spearman's rank correlation
#'     \item \code{"kendall"}: Kendall's tau correlation
#'   }
#' @param plot_top_n Positive integer specifying the number of top variable
#'   features/samples to include in the plot. Must be > 0. Default is 1000.
#' @param plot_what Character string specifying what to plot. One of:
#'   \itemize{
#'     \item \code{"Features"}: Plot correlations between features (default)
#'     \item \code{"Samples"}: Plot correlations between samples
#'   }
#' @param plot_type Character string specifying the plot type. One of:
#'   \itemize{
#'     \item \code{"correlation"}: Correlation heatmap with significance masking (default)
#'     \item \code{"hierarchical"}: Hierarchical clustering heatmap of raw values
#'   }
#' @param show_rownames Logical. Whether to display row names. Default is \code{TRUE}.
#' @param show_colnames Logical. Whether to display column names. Default is \code{TRUE}.
#' @param clustering_distance_rows Character string specifying the distance metric
#'   for row clustering. See \code{\link[stats]{dist}} for details. One of:
#'   \code{"euclidean"} (default), \code{"maximum"}, \code{"manhattan"},
#'   \code{"canberra"}, \code{"binary"}, \code{"minkowski"}.
#' @param clustering_distance_cols Character string specifying the distance metric
#'   for column clustering. Same options as \code{clustering_distance_rows}.
#'   Default is \code{"euclidean"}.
#' @param clustering_method Character string specifying the clustering algorithm.
#'   See \code{\link[stats]{hclust}} for details. One of:
#'   \itemize{
#'     \item \code{"ward.D"}: Ward's minimum variance method (default)
#'     \item \code{"ward.D2"}: Implements Ward's (1963) criterion
#'     \item \code{"single"}: Single linkage clustering
#'     \item \code{"complete"}: Complete linkage clustering
#'     \item \code{"average"}: UPGMA clustering
#'     \item \code{"mcquitty"}: WPGMA clustering
#'     \item \code{"median"}: WPGMC clustering
#'     \item \code{"centroid"}: UPGMC clustering
#'   }
#' @param significance_threshold Numeric value between 0 and 1 specifying the
#'   p-value threshold for significance masking in correlation plots.
#'   Default is 0.05.
#' @param color_palette Character vector of length 3 specifying colors for
#'   negative correlations, neutral/non-significant, and positive correlations.
#'   Default is \code{c("blue", "white", "red")}.
#' @param fontsize_main Numeric value for main title font size. Default is 12.
#' @param fontsize_labels Numeric value for axis label font size. Default is 8.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{plot}}{The generated heatmap plot object}
#'     \item{\code{filtered_data}}{Data frame of filtered data used for plotting}
#'     \item{\code{correlation_matrix}}{Correlation matrix (correlation plots only)}
#'     \item{\code{p_values}}{Matrix of p-values (correlation plots only)}
#'     \item{\code{top_features}}{Names of selected top variable features/samples}
#'     \item{\code{iqr_values}}{Data frame of IQR values for all features/samples}
#'     \item{\code{parameters}}{List of all function parameters used}
#'     \item{\code{summary_stats}}{Summary statistics of the analysis}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data structure
#'   \item Automatically selects between \code{data_scaledPCA_varFiltered} and \code{data_scaledPCA_merged} based on \code{auto_merge_replicates} parameter
#'   \item Calculates interquartile ranges (IQR) to identify most variable features/samples
#'   \item Selects top N most variable features/samples based on IQR
#'   \item For correlation plots: computes correlation matrix with p-values and applies significance masking
#'   \item For hierarchical plots: uses raw filtered data with clustering
#'   \item Generates publication-ready heatmap with customizable aesthetics
#' }
#'
#' For correlation plots, non-significant correlations (p >= significance_threshold)
#' are displayed as neutral color (white by default) to highlight statistically
#' significant relationships.
#'
#' @examples
#' \dontrun{
#' # Basic correlation heatmap
#' result <- plot_CorrelationHeatmap(
#'   data = preprocessed_data,
#'   method = "pearson",
#'   plot_top_n = 500
#' )
#'
#' # Hierarchical clustering of samples
#' result <- plot_CorrelationHeatmap(
#'   data = preprocessed_data,
#'   plot_what = "Samples",
#'   plot_type = "hierarchical",
#'   clustering_method = "ward.D2"
#' )
#'
#' # Custom correlation plot with Spearman correlation
#' result <- plot_CorrelationHeatmap(
#'   data = preprocessed_data,
#'   method = "spearman",
#'   significance_threshold = 0.01,
#'   color_palette = c("darkblue", "grey90", "darkred")
#' )
#' }
#'
#' @seealso
#' \code{\link[stats]{cor}}, \code{\link[Hmisc]{rcorr}}, \code{\link[pheatmap]{pheatmap}}
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom stats IQR cor
#' @importFrom Hmisc rcorr
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom grid grid.newpage
#'
#' @references Tukey, J. W. (1977). Exploratory Data Analysis. Reading: Addison-Wesley. (for IQR)
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988). The New S Language. Wadsworth & Brooks/Cole. (for cor)
#' @references Kendall, M. G. (1938). A new measure of rank correlation, Biometrika, 30, 81â€"93. doi:10.1093/biomet/30.1-2.81. (for cor)
#' @references Kendall, M. G. (1945). The treatment of ties in rank problems. Biometrika, 33 239â€"251. doi:10.1093/biomet/33.3.239 (for cor)
#' @references Hollander M. and Wolfe D.A. (1973). Nonparametric Statistical Methods. New York: Wiley. (for Hmisc::rcorr)
#' @references Press WH, Flannery BP, Teukolsky SA, Vetterling, WT (1988): Numerical Recipes in C. Cambridge: Cambridge University Press. (for Hmisc::rcorr)
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. (for stats::dist and stats::hclust in pheatmap)
#' @references Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979) Multivariate Analysis. Academic Press. (for stats::dist in pheatmap)
#' @references Borg, I. and Groenen, P. (1997) Modern Multidimensional Scaling. Theory and Applications. Springer. (for stats::dist in pheatmap)
#' @references Everitt, B. (1974). Cluster Analysis. London: Heinemann Educ. Books. (for stats::hclust in pheatmap)
#' @references Hartigan, J.A. (1975). Clustering Algorithms. New York: Wiley. (for stats::hclust in pheatmap)
#' @references Sneath, P. H. A. and R. R. Sokal (1973). Numerical Taxonomy. San Francisco: Freeman. (for stats::hclust in pheatmap)
#' @references Anderberg, M. R. (1973). Cluster Analysis for Applications. Academic Press: New York. (for stats::hclust in pheatmap)
#' @references Gordon, A. D. (1999). Classification. Second Edition. London: Chapman and Hall / CRC (for stats::hclust in pheatmap)
#' @references Murtagh, F. (1985). â€œMultidimensional Clustering Algorithmsâ€, in COMPSTAT Lectures 4. Wuerzburg: Physica-Verlag (for algorithmic details of algorithms used). (for stats::hclust in pheatmap)
#' @references McQuitty, L.L. (1966). Similarity Analysis by Reciprocal Pairs for Discrete and Continuous Data. Educational and Psychological Measurement, 26, 825â€"831. doi:10.1177/001316446602600402. (for stats::hclust in pheatmap)
#' @references Legendre, P. and L. Legendre (2012). Numerical Ecology, 3rd English ed. Amsterdam: Elsevier Science BV. (for stats::hclust in pheatmap)
#' @references Murtagh, Fionn and Legendre, Pierre (2014). Ward's hierarchical agglomerative clustering method: which algorithms implement Ward's criterion? Journal of Classification, 31, 274â€"295. doi:10.1007/s00357-014-9161-z. (for stats::hclust in pheatmap)
#'
#' @export
plot_CorrelationHeatmap <- function(
    data,
    method = "pearson",
    plot_top_n = 1000,
    plot_what = "Features",
    plot_type = "correlation",
    show_rownames = TRUE,
    show_colnames = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D",
    significance_threshold = 0.05,
    color_palette = c("blue", "white", "red"),
    fontsize_main = 12,
    fontsize_labels = 8
) {

  # Helper function for input validation
  validate_inputs_corrheatmap <- function() {
    # Check required packages
    required_packages <- c("stats", "Hmisc", "pheatmap", "grDevices", "grid")
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
    if (length(missing_packages) > 0) {
      stop(sprintf("Required packages missing: %s. Please install them.",
                   paste(missing_packages, collapse = ", ")))
    }

    # Validate data structure and determine correct data source
    if (!is.list(data)) {
      stop("'data' must be a list object")
    }

    # Determine which data component to use based on auto_merge_replicates parameter
    use_merged_data <- FALSE
    if ("Parameters" %in% names(data) &&
        is.list(data$Parameters) &&
        "auto_merge_replicates" %in% names(data$Parameters)) {
      use_merged_data <- isTRUE(data$Parameters$auto_merge_replicates)
    }

    # Check for appropriate data component
    if (use_merged_data) {
      if (!"data_scaledPCA_merged" %in% names(data)) {
        stop("'data' must contain 'data_scaledPCA_merged' component when auto_merge_replicates is TRUE")
      }
      data_component_name <- "data_scaledPCA_merged"
      data_matrix <- data$data_scaledPCA_merged
    } else {
      if (!"data_scaledPCA_varFiltered" %in% names(data)) {
        stop("'data' must contain 'data_scaledPCA_varFiltered' component when auto_merge_replicates is FALSE or not specified")
      }
      data_component_name <- "data_scaledPCA_varFiltered"
      data_matrix <- data$data_scaledPCA_varFiltered
    }

    if (!is.matrix(data_matrix) && !is.data.frame(data_matrix)) {
      stop(sprintf("'%s' must be a matrix or data frame", data_component_name))
    }

    if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
      stop("Data matrix cannot be empty")
    }

    # Validate method
    valid_methods <- c("pearson", "spearman", "kendall")
    if (!method %in% valid_methods) {
      stop(sprintf("'method' must be one of: %s", paste(valid_methods, collapse = ", ")))
    }

    # Validate plot_top_n
    if (!is.numeric(plot_top_n) || length(plot_top_n) != 1 || plot_top_n <= 0 || !is.finite(plot_top_n)) {
      stop("'plot_top_n' must be a positive finite number")
    }
    plot_top_n <- as.integer(plot_top_n)

    # Validate plot_what
    valid_plot_what <- c("Samples", "Features")
    if (!plot_what %in% valid_plot_what) {
      stop(sprintf("'plot_what' must be one of: %s", paste(valid_plot_what, collapse = ", ")))
    }

    # Validate plot_type
    valid_plot_types <- c("correlation", "hierarchical")
    if (!plot_type %in% valid_plot_types) {
      stop(sprintf("'plot_type' must be one of: %s", paste(valid_plot_types, collapse = ", ")))
    }

    # Validate logical parameters
    if (!is.logical(show_rownames) || length(show_rownames) != 1) {
      stop("'show_rownames' must be a single logical value")
    }
    if (!is.logical(show_colnames) || length(show_colnames) != 1) {
      stop("'show_colnames' must be a single logical value")
    }

    # Validate distance methods
    valid_distances <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
    if (!clustering_distance_rows %in% valid_distances) {
      stop(sprintf("'clustering_distance_rows' must be one of: %s", paste(valid_distances, collapse = ", ")))
    }
    if (!clustering_distance_cols %in% valid_distances) {
      stop(sprintf("'clustering_distance_cols' must be one of: %s", paste(valid_distances, collapse = ", ")))
    }

    # Validate clustering method
    valid_clustering <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
    if (!clustering_method %in% valid_clustering) {
      stop(sprintf("'clustering_method' must be one of: %s", paste(valid_clustering, collapse = ", ")))
    }

    # Validate significance threshold
    if (!is.numeric(significance_threshold) || length(significance_threshold) != 1 ||
        significance_threshold < 0 || significance_threshold > 1) {
      stop("'significance_threshold' must be a number between 0 and 1")
    }

    # Validate color palette
    if (!is.character(color_palette) || length(color_palette) != 3) {
      stop("'color_palette' must be a character vector of length 3")
    }

    # Validate font sizes
    if (!is.numeric(fontsize_main) || length(fontsize_main) != 1 || fontsize_main <= 0) {
      stop("'fontsize_main' must be a positive number")
    }
    if (!is.numeric(fontsize_labels) || length(fontsize_labels) != 1 || fontsize_labels <= 0) {
      stop("'fontsize_labels' must be a positive number")
    }

    # Return both the validated plot_top_n and the determined data info
    return(list(plot_top_n = plot_top_n,
                data_matrix = data_matrix,
                data_component_name = data_component_name,
                use_merged_data = use_merged_data))
  }

  # Helper function to prepare data matrix
  prepare_data_matrix_corrheatmap <- function(data_matrix, plot_what) {
    # Convert to matrix if data.frame
    if (is.data.frame(data_matrix)) {
      data_matrix <- as.matrix(data_matrix)
    }

    # Remove rows/columns with all missing values
    if (plot_what == "Samples") {
      data_matrix <- t(data_matrix)
      # Remove samples (rows) with all NA
      complete_samples <- rowSums(!is.na(data_matrix)) > 0
      if (!any(complete_samples)) {
        stop("All samples contain only missing values")
      }
      data_matrix <- data_matrix[complete_samples, , drop = FALSE]
    } else {
      # Remove features (columns) with all NA
      complete_features <- colSums(!is.na(data_matrix)) > 0
      if (!any(complete_features)) {
        stop("All features contain only missing values")
      }
      data_matrix <- data_matrix[, complete_features, drop = FALSE]
    }

    return(data_matrix)
  }

  # Helper function to calculate IQR values efficiently
  calculate_iqr_values <- function(data_matrix) {
    iqr_values <- apply(data_matrix, 2, function(x) {
      if (sum(!is.na(x)) < 2) return(NA_real_)
      stats::IQR(x, na.rm = TRUE)
    })

    # Remove features with NA IQR (insufficient data)
    valid_iqr <- !is.na(iqr_values)
    if (!any(valid_iqr)) {
      stop("Cannot calculate IQR for any features (insufficient non-missing data)")
    }

    iqr_df <- data.frame(
      Feature = names(iqr_values)[valid_iqr],
      IQR = iqr_values[valid_iqr],
      stringsAsFactors = FALSE
    )
    iqr_df <- iqr_df[order(iqr_df$IQR, decreasing = TRUE), ]
    rownames(iqr_df) <- iqr_df$Feature

    return(iqr_df)
  }

  # Helper function to create correlation plot
  create_correlation_plot <- function(filtered_df, method, significance_threshold,
                                      color_palette, clustering_params, font_params, plot_what) {

    # Check if we have enough data for correlation
    if (ncol(filtered_df) < 2) {
      stop("Need at least 2 features/samples for correlation analysis")
    }

    # Compute correlation with error handling
    tryCatch({
      cor_results <- Hmisc::rcorr(as.matrix(filtered_df), type = method)
    }, error = function(e) {
      stop(sprintf("Error computing correlations: %s", e$message))
    })

    cor_matrix <- cor_results$r
    p_values <- cor_results$P

    # Handle case where p-values might be NULL (single variable case)
    if (is.null(p_values)) {
      p_values <- matrix(0, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
      dimnames(p_values) <- dimnames(cor_matrix)
    }

    # Create masked correlation matrix
    cor_matrix_masked <- cor_matrix
    cor_matrix_masked[p_values >= significance_threshold] <- 0

    # Create color palette
    my_palette <- grDevices::colorRampPalette(color_palette)(50)

    # Clear graphics device
    grid::grid.newpage()

    # Create plot with error handling
    plot_obj <- tryCatch({
      pheatmap::pheatmap(
        cor_matrix_masked,
        color = my_palette,
        breaks = seq(-1, 1, length.out = 51),
        display_numbers = FALSE,
        show_rownames = if (plot_what == "Samples") TRUE else clustering_params$show_rownames,
        show_colnames = if (plot_what == "Samples") TRUE else clustering_params$show_colnames,
        fontsize = font_params$fontsize_main,
        fontsize_row = font_params$fontsize_labels,
        fontsize_col = font_params$fontsize_labels,
        clustering_distance_rows = clustering_params$distance_rows,
        clustering_distance_cols = clustering_params$distance_cols,
        clustering_method = clustering_params$method,
        main = sprintf("Correlation Heatmap of %s (p >= %.3f shown as neutral)",
                       plot_what, significance_threshold),
        silent = TRUE
      )
    }, error = function(e) {
      stop(sprintf("Error creating correlation heatmap: %s", e$message))
    })

    return(list(
      plot = plot_obj,
      correlation_matrix = cor_matrix,
      p_values = p_values,
      masked_matrix = cor_matrix_masked
    ))
  }

  # Helper function to create hierarchical plot
  create_hierarchical_plot <- function(filtered_df, color_palette, clustering_params,
                                       font_params, plot_what) {

    # Create color palette for raw values
    my_palette <- grDevices::colorRampPalette(color_palette)(50)

    # Clear graphics device
    grid::grid.newpage()

    # Create plot with error handling
    plot_obj <- tryCatch({
      pheatmap::pheatmap(
        filtered_df,
        color = my_palette,
        show_rownames = if (plot_what == "Samples") TRUE else clustering_params$show_rownames,
        show_colnames = if (plot_what == "Samples") TRUE else clustering_params$show_colnames,
        fontsize = font_params$fontsize_main,
        fontsize_row = font_params$fontsize_labels,
        fontsize_col = font_params$fontsize_labels,
        clustering_distance_rows = clustering_params$distance_rows,
        clustering_distance_cols = clustering_params$distance_cols,
        clustering_method = clustering_params$method,
        main = sprintf("Hierarchical Clustering Heatmap of %s", plot_what),
        silent = TRUE
      )
    }, error = function(e) {
      stop(sprintf("Error creating hierarchical heatmap: %s", e$message))
    })

    return(list(plot = plot_obj))
  }

  # Main function execution starts here

  # Validate all inputs and get data matrix info
  validation_results <- validate_inputs_corrheatmap()
  plot_top_n <- validation_results$plot_top_n
  selected_data_matrix <- validation_results$data_matrix
  data_component_name <- validation_results$data_component_name

  # Prepare data matrix
  data_matrix <- prepare_data_matrix_corrheatmap(selected_data_matrix, plot_what)

  # Calculate IQR values
  iqr_results <- calculate_iqr_values(data_matrix)

  # Select top features
  n_available <- nrow(iqr_results)
  n_select <- min(plot_top_n, n_available)
  top_features <- iqr_results$Feature[1:n_select]

  # Filter data to top features
  filtered_df <- data_matrix[, top_features, drop = FALSE]

  # Prepare clustering and font parameters
  clustering_params <- list(
    distance_rows = clustering_distance_rows,
    distance_cols = clustering_distance_cols,
    method = clustering_method,
    show_rownames = show_rownames,
    show_colnames = show_colnames
  )

  font_params <- list(
    fontsize_main = fontsize_main,
    fontsize_labels = fontsize_labels
  )

  # Create appropriate plot
  if (plot_type == "correlation") {
    plot_results <- create_correlation_plot(
      filtered_df, method, significance_threshold, color_palette,
      clustering_params, font_params, plot_what
    )
  } else {
    plot_results <- create_hierarchical_plot(
      filtered_df, color_palette, clustering_params, font_params, plot_what
    )
  }

  # Compile summary statistics
  summary_stats <- list(
    n_original_features = ifelse(plot_what == "Samples", nrow(selected_data_matrix),
                                 ncol(selected_data_matrix)),
    n_selected_features = n_select,
    n_samples = ifelse(plot_what == "Samples", ncol(filtered_df), nrow(filtered_df)),
    missing_data_percent = round(100 * sum(is.na(filtered_df)) / (nrow(filtered_df) * ncol(filtered_df)), 2),
    iqr_range = c(min = min(iqr_results$IQR), max = max(iqr_results$IQR)),
    data_source = data_component_name,
    auto_merge_replicates = validation_results$use_merged_data
  )

  if (plot_type == "correlation") {
    n_comparisons <- ncol(filtered_df) * (ncol(filtered_df) - 1) / 2
    n_significant <- sum(plot_results$p_values[upper.tri(plot_results$p_values)] < significance_threshold, na.rm = TRUE)
    summary_stats$n_correlations_tested <- n_comparisons
    summary_stats$n_significant_correlations <- n_significant
    summary_stats$percent_significant <- round(100 * n_significant / n_comparisons, 2)
  }

  # Compile results
  results <- list(
    plot = plot_results$plot,
    filtered_data = as.data.frame(filtered_df),
    top_features = top_features,
    iqr_values = iqr_results,
    parameters = list(
      method = method,
      plot_top_n = plot_top_n,
      plot_what = plot_what,
      plot_type = plot_type,
      significance_threshold = significance_threshold,
      clustering_distance_rows = clustering_distance_rows,
      clustering_distance_cols = clustering_distance_cols,
      clustering_method = clustering_method,
      color_palette = color_palette
    ),
    summary_stats = summary_stats,
    function_origin = "plot_CorrelationHeatmap"
  )

  # Add plot-specific results
  if (plot_type == "correlation") {
    results$correlation_matrix <- as.data.frame(plot_results$correlation_matrix)
    results$p_values <- as.data.frame(plot_results$p_values)
    results$masked_correlation_matrix <- as.data.frame(plot_results$masked_matrix)
  }

  return(results)
}
