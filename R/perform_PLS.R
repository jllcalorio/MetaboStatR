#' Perform Partial Least Squares (PLS) Related Methods
#'
#' This function performs various Partial Least Squares (PLS) methods including PLS, PLS-DA,
#' OPLS-DA, and sparse PLS-DA (sPLS-DA) on metabolomics or other omics data. It provides
#' comprehensive analysis with visualization and feature importance assessment.
#'
#' @param data A list containing processed data with the following required components:
#'   \itemize{
#'     \item \code{data_scaledPLS_rsdFiltered_varFiltered}: Matrix or data.frame of processed data
#'     \item \code{Metadata}: Data.frame containing sample metadata with 'Group' and 'Group_' columns
#'   }
#' @param method Character string specifying the PLS method to use. Options are:
#'   \itemize{
#'     \item \code{"pls"}: Standard PLS regression (experimental)
#'     \item \code{"plsda"}: PLS Discriminant Analysis (experimental)
#'     \item \code{"oplsda"}: Orthogonal PLS Discriminant Analysis (default)
#'     \item \code{"splsda"}: sparse PLS Discriminant Analysis (experimental)
#'   }
#' @param arrangeLevels Character vector specifying the order of group levels for analysis.
#'   If NULL (default), all unique groups will be used in alphabetical order.
#' @param includeQC Logical indicating whether to include QC samples in visualizations.
#'   Note: QC samples are never included in model building. Default is FALSE.
#' @param predI Integer specifying the number of predictive components. Default is 1.
#' @param orthoI Integer specifying the number of orthogonal components for OPLS-DA.
#'   If NA (default), optimal number is determined automatically.
#' @param crossvalI Integer specifying the number of cross-validation folds. Default is 10.
#' @param permI Integer specifying the number of permutations for validation. Default is 20.
#' @param scaleC Character string specifying scaling method. Options are "none", "center",
#'   "pareto", "unit". Default is "none".
#' @param top_features Integer specifying the number of top VIP features to display in plots.
#'   Default is 20.
#' @param keepX Integer vector for sPLS-DA specifying the number of variables to keep
#'   on each component. Default is NULL (automatic selection).
#' @param ncomp Integer specifying the number of components for sPLS-DA. Default is 2.
#' @param validation Character string for sPLS-DA validation method. Default is "Mfold".
#' @param folds Integer for sPLS-DA cross-validation folds. Default is 10.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return A list containing analysis results with the following components:
#'   \itemize{
#'     \item \code{method}: The PLS method used
#'     \item \code{data_used}: The data matrix used for analysis
#'     \item \code{results_[method]_[comparison]}: Model results for each pairwise comparison
#'     \item \code{data_VIPScores_[comparison]}: VIP scores for each comparison (PLS-DA/OPLS-DA)
#'     \item \code{data_Abundance_[comparison]}: Abundance data for top features
#'     \item \code{plot_VIPAbundance_[comparison]}: Combined VIP and abundance plots
#'     \item \code{data_SPlot_[comparison]}: S-plot data (OPLS-DA only)
#'     \item \code{plot_SPlot_[comparison]}: S-plots (OPLS-DA only)
#'     \item \code{plot_Scores_[comparison]}: Score plots for each comparison
#'     \item \code{summary}: Summary statistics and model performance metrics
#'   }
#'
#' @examples
#' \dontrun{
#' # Example data structure
#' data <- list(
#'   data_scaledPLS_rsdFiltered_varFiltered = matrix(rnorm(1000), nrow = 50, ncol = 20),
#'   Metadata = data.frame(
#'     Group = c(rep(c("Control", "Treatment"), c(20, 20)), rep(c("SQC", "EQC"), c(5, 5))),
#'     Group_ = c(rep(c("Control", "Treatment"), c(20, 20)), rep("QC", 10))
#'   )
#' )
#' colnames(data$data_scaledPLS_rsdFiltered_varFiltered) <- paste0("Feature_", 1:20)
#'
#' # Perform OPLS-DA
#' results <- perform_PLS(data, method = "oplsda")
#'
#' # Perform PLS-DA with specific group arrangement
#' results <- perform_PLS(data, method = "plsda",
#'                       arrangeLevels = c("Control", "Treatment"))
#'
#' # Perform sparse PLS-DA
#' results <- perform_PLS(data, method = "splsda", ncomp = 3, keepX = c(10, 5, 5))
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @import ggplot2
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import purrr
#' @importFrom ropls opls getVipVn
#' @importFrom mixOmics splsda perf tune.splsda plotIndiv plotVar
#' @importFrom patchwork wrap_plots
#' @importFrom stats cov cor reorder
#' @importFrom utils combn
#' @importFrom grDevices colorRampPalette
#'
#' @references Eriksson et al. (2006). Multi- and Megarvariate Data Analysis. Umetrics Academy. Rosipal and Kramer (2006). Overview and recent advances in partial least squares Tenenhaus (1990). La regression PLS : theorie et pratique. Technip. Wehrens (2011). Chemometrics with R. Springer. Wold et al. (2001). PLS-regression: a basic tool of chemometrics (for PLS)
#' @references Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics feature selection and multiple data integration. PLoS Comput Biol 13(11): e1005752 (for sPLS)
#' @references Galindo-Prieto B., Eriksson L. and Trygg J. (2014). Variable influence on projection (VIP) for orthogonal projections to latent structures (OPLS). Journal of Chemometrics 28, 623-632. (for getVipVn)
#'
#' @export
perform_PLS <- function(data,
                        method = "oplsda",
                        arrangeLevels = NULL,
                        includeQC = FALSE,
                        predI = 1,
                        orthoI = NA,
                        crossvalI = 10,
                        permI = 20,
                        scaleC = "none",
                        top_features = 20,
                        keepX = NULL,
                        ncomp = 2,
                        validation = "Mfold",
                        folds = 10,
                        verbose = TRUE) {

  # # Helper function for memory-efficient operations
  # gc_if_needed <- function() {
  #   if (gc.time()[1] %% 10 == 0) invisible(gc())
  # }

  # Input validation
  validate_inputs(data, method, arrangeLevels, predI, orthoI, crossvalI,
                  permI, scaleC, top_features, ncomp, verbose)

  # Initialize results list
  results <- list(
    method = method,
    parameters = list(
      predI = predI, orthoI = orthoI, crossvalI = crossvalI,
      permI = permI, scaleC = scaleC, top_features = top_features,
      ncomp = ncomp, keepX = keepX
    ),
    summary = list()
  )

  if (verbose) cat("Starting", toupper(method), "analysis...\n")

  # Prepare data
  pls_data <- prepare_pls_data(data, verbose)
  results$data_used <- pls_data$data_matrix

  # Get non-QC indices
  non_qc_indices <- get_non_qc_indices(data, verbose)

  # Handle group arrangements
  group_combinations <- setup_group_combinations(data, arrangeLevels, non_qc_indices, verbose)

  # Perform analysis based on method
  if (method == "splsda") {
    results <- perform_splsda_analysis(data, pls_data, non_qc_indices,
                                       group_combinations, results,
                                       ncomp, keepX, validation, folds,
                                       top_features, verbose)
  } else {
    results <- perform_pls_analysis(data, pls_data, non_qc_indices,
                                    group_combinations, results, method,
                                    predI, orthoI, crossvalI, permI,
                                    scaleC, top_features, verbose)
  }

  # Generate summary
  results$summary <- generate_analysis_summary(results, method, verbose)

  # Assign class to results and return
  class(results) <- c("pls_results", "list")

  if (verbose) cat("Analysis completed successfully!\n")
  # gc_if_needed()

  return(results)
}

# Helper function: Input validation
validate_inputs <- function(data, method, arrangeLevels, predI, orthoI,
                            crossvalI, permI, scaleC, top_features, ncomp, verbose) {

  # Check required packages
  required_packages <- c("ggplot2", "dplyr", "tibble", "tidyr", "purrr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }

  # Validate data structure
  if (!is.list(data)) {
    stop("'data' must be a list")
  }

  required_components <- c("data_scaledPLS_rsdFiltered_varFiltered", "Metadata")
  missing_components <- required_components[!required_components %in% names(data)]

  if (length(missing_components) > 0) {
    stop("Missing required data components: ", paste(missing_components, collapse = ", "))
  }

  # Validate method
  valid_methods <- c("pls", "plsda", "oplsda", "splsda")
  if (!method %in% valid_methods) {
    stop("Method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Method-specific package checks
  if (method %in% c("pls", "plsda", "oplsda")) {
    if (!requireNamespace("ropls", quietly = TRUE)) {
      stop("Package 'ropls' is required for ", method)
    }
  }

  if (method == "splsda") {
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
      stop("Package 'mixOmics' is required for sPLS-DA")
    }
  }

  # Validate metadata
  required_metadata_cols <- c("Group", "Group_")
  missing_metadata_cols <- required_metadata_cols[!required_metadata_cols %in% colnames(data$Metadata)]

  if (length(missing_metadata_cols) > 0) {
    stop("Metadata must contain the following columns: ", paste(missing_metadata_cols, collapse = ", "))
  }

  # Validate numeric parameters
  if (!is.numeric(predI) || predI < 1) {
    stop("predI must be a positive integer")
  }

  if (!is.na(orthoI) && (!is.numeric(orthoI) || orthoI < 0)) {
    stop("orthoI must be a non-negative integer or NA")
  }

  if (!is.numeric(crossvalI) || crossvalI < 2) {
    stop("crossvalI must be an integer >= 2")
  }

  if (!is.numeric(top_features) || top_features < 1) {
    stop("top_features must be a positive integer")
  }

  # Validate scaling method
  valid_scales <- c("none", "center", "pareto", "unit")
  if (!scaleC %in% valid_scales) {
    stop("scaleC must be one of: ", paste(valid_scales, collapse = ", "))
  }

  if (verbose) cat("Input validation passed.\n")
}

# Helper function: Prepare PLS data
prepare_pls_data <- function(data, verbose) {

  if (verbose) cat("Preparing data for analysis...\n")

  # Check if data exists and is not empty
  if (is.null(data$data_scaledPLS_rsdFiltered_varFiltered)) {
    stop("data_scaledPLS_rsdFiltered_varFiltered is NULL. Please ensure data processing is completed.")
  }

  data_matrix <- data$data_scaledPLS_rsdFiltered_varFiltered

  # Check if data is empty
  if (length(data_matrix) == 0 || nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    stop("data_scaledPLS_rsdFiltered_varFiltered is empty. Please ensure data processing is completed.")
  }

  if (verbose) cat("Using processed and filtered data (", nrow(data_matrix), " samples x ", ncol(data_matrix), " features).\n")

  # Ensure data is in proper format
  data_matrix <- as.data.frame(data_matrix)

  # Check for missing values
  if (any(is.na(data_matrix))) {
    warning("Missing values detected in data. Consider imputation before analysis.")
  }

  # Check for zero variance features
  zero_var_features <- names(data_matrix)[apply(data_matrix, 2, var, na.rm = TRUE) == 0]
  if (length(zero_var_features) > 0) {
    warning("Zero variance features detected: ", paste(head(zero_var_features, 5), collapse = ", "))
  }

  return(list(data_matrix = data_matrix))
}

# Helper function: Get non-QC indices
get_non_qc_indices <- function(data, verbose) {

  # Use Group_ column to identify QC samples (where Group_ == "QC")
  non_qc_indices <- data$Metadata$Group_ != "QC"
  n_bio <- sum(non_qc_indices)
  n_qc <- sum(!non_qc_indices)

  if (verbose) {
    cat("Found", n_bio, "biological samples and", n_qc, "QC samples.\n")

    # Show the specific QC types from the Group column
    qc_types <- unique(data$Metadata$Group[!non_qc_indices])
    if (length(qc_types) > 0) {
      cat("QC types found:", paste(qc_types, collapse = ", "), "\n")
    }
  }

  if (n_bio < 6) {
    warning("Very few biological samples (n=", n_bio, "). Results may be unreliable.")
  }

  return(non_qc_indices)
}

# Helper function: Setup group combinations
setup_group_combinations <- function(data, arrangeLevels, non_qc_indices, verbose) {

  # Use Group column for biological samples (exclude QC samples)
  bio_groups <- data$Metadata$Group[non_qc_indices]
  unique_groups <- unique(bio_groups)

  if (is.null(arrangeLevels)) {
    groups_to_use <- sort(unique_groups)
    if (verbose) cat("Using all biological groups in alphabetical order:", paste(groups_to_use, collapse = ", "), "\n")
  } else {
    # Validate arrangeLevels
    missing_groups <- setdiff(arrangeLevels, unique_groups)
    if (length(missing_groups) > 0) {
      stop("Groups not found in biological samples: ", paste(missing_groups, collapse = ", "),
           "\nAvailable biological groups: ", paste(unique_groups, collapse = ", "))
    }
    groups_to_use <- arrangeLevels
    if (verbose) cat("Using specified biological groups:", paste(groups_to_use, collapse = ", "), "\n")
  }

  if (length(groups_to_use) < 2) {
    stop("At least 2 biological groups are required for analysis. Found: ", length(groups_to_use))
  }

  # Generate pairwise combinations
  group_combinations <- combn(groups_to_use, 2, simplify = FALSE)

  if (verbose) {
    cat("Will perform", length(group_combinations), "pairwise comparisons:\n")
    for (i in seq_along(group_combinations)) {
      cat("  ", group_combinations[[i]][1], "vs", group_combinations[[i]][2], "\n")
    }
  }

  return(group_combinations)
}

# Helper function: Perform PLS/PLS-DA/OPLS-DA analysis
perform_pls_analysis <- function(data, pls_data, non_qc_indices, group_combinations,
                                 results, method, predI, orthoI, crossvalI, permI,
                                 scaleC, top_features, verbose) {

  data_matrix <- pls_data$data_matrix

  # Loop through each pairwise comparison
  for (i in seq_along(group_combinations)) {

    pair <- group_combinations[[i]]
    group1 <- pair[1]
    group2 <- pair[2]

    if (verbose) cat("\nAnalyzing", group1, "vs", group2, "...\n")

    # Subset data for current pair
    current_groups <- data$Metadata$Group %in% c(group1, group2) & non_qc_indices

    if (sum(current_groups) < 6) {
      warning("Too few samples for ", group1, " vs ", group2, " comparison (n=",
              sum(current_groups), "). Skipping.")
      next
    }

    data_pair <- data_matrix[current_groups, ]
    y_pair <- factor(data$Metadata$Group[current_groups])

    # comparison_label <- paste0(group1, "vs", group2, sep = "_")
    comparison_label <- paste0(group1, " vs. ", group2)

    # Perform analysis based on method
    model_results <- tryCatch({

      if (method == "oplsda") {
        ropls::opls(x = data_pair, y = y_pair, predI = predI, orthoI = orthoI,
                    algoC = "nipals", crossvalI = crossvalI, log10L = FALSE,
                    permI = permI, scaleC = scaleC, subset = NULL,
                    plotSubC = paste0(group1, " vs ", group2),
                    fig.pdfC = "interactive")

      } else if (method == "plsda") {
        ropls::opls(x = data_pair, y = y_pair, predI = predI, orthoI = 0,
                    algoC = "nipals", crossvalI = crossvalI, log10L = FALSE,
                    permI = permI, scaleC = scaleC, subset = NULL,
                    plotSubC = paste0(group1, " vs ", group2))

      } else if (method == "pls") {
        # For PLS regression, we need continuous Y
        y_numeric <- as.numeric(y_pair) - 1  # Convert to 0/1
        ropls::opls(x = data_pair, y = y_numeric, predI = predI, orthoI = 0,
                    algoC = "nipals", crossvalI = crossvalI, log10L = FALSE,
                    permI = permI, scaleC = scaleC, subset = NULL,
                    plotSubC = paste0(group1, " vs ", group2))
      }

    }, error = function(e) {
      warning("Failed to build model for ", group1, " vs ", group2, ": ", e$message)
      return(NULL)
    })

    if (is.null(model_results)) next

    # Store model results
    results[[paste0("results_", toupper(method), "_", comparison_label)]] <- model_results

    # Generate plots and additional analyses
    if (inherits(model_results, "opls") && length(model_results@summaryDF) > 0) {

      results <- generate_feature_analysis(results, model_results, data_pair, y_pair,
                                           group1, group2, method, top_features, verbose)

      results <- generate_score_plots(results, model_results, y_pair, group1, group2, method)

      if (method == "oplsda") {
        results <- generate_splot(results, model_results, data_pair, group1, group2)
      }

      # gc_if_needed()
    }
  }

  return(results)
}

# Helper function: Perform sPLS-DA analysis
perform_splsda_analysis <- function(data, pls_data, non_qc_indices, group_combinations,
                                    results, ncomp, keepX, validation, folds,
                                    top_features, verbose) {

  data_matrix <- pls_data$data_matrix

  for (i in seq_along(group_combinations)) {

    pair <- group_combinations[[i]]
    group1 <- pair[1]
    group2 <- pair[2]

    if (verbose) cat("\nAnalyzing", group1, "vs", group2, "with sPLS-DA...\n")

    current_groups <- data$Metadata$Group %in% c(group1, group2) & non_qc_indices
    data_pair <- data_matrix[current_groups, ]
    y_pair <- factor(data$Metadata$Group[current_groups])

    # comparison_label <- paste0(group1, "vs", group2, sep = "_")
    comparison_label <- paste0(group1, " vs. ", group2)

    model_results <- tryCatch({

      # Tune keepX if not provided
      if (is.null(keepX)) {
        if (verbose) cat("  Tuning keepX parameters...\n")
        tune_results <- mixOmics::tune.splsda(X = data_pair, Y = y_pair,
                                              ncomp = ncomp, validation = validation,
                                              folds = folds, test.keepX = c(seq(5, 50, 5)))
        optimal_keepX <- tune_results$choice.keepX
      } else {
        optimal_keepX <- keepX[1:min(length(keepX), ncomp)]
      }

      # Fit sPLS-DA model
      mixOmics::splsda(X = data_pair, Y = y_pair, ncomp = ncomp, keepX = optimal_keepX)

    }, error = function(e) {
      warning("Failed to build sPLS-DA model for ", group1, " vs ", group2, ": ", e$message)
      return(NULL)
    })

    if (is.null(model_results)) next

    results[[paste0("results_SPLSDA_", comparison_label)]] <- model_results

    # Generate sPLS-DA specific plots and analyses
    results <- generate_splsda_analysis(results, model_results, data_pair, y_pair,
                                        group1, group2, top_features, verbose)

    # gc_if_needed()
  }

  return(results)
}

# Helper function: Generate feature analysis (VIP, abundance)
generate_feature_analysis <- function(results, model_results, data_pair, y_pair,
                                      group1, group2, method, top_features, verbose) {

  # comparison_label <- paste0(group1, "vs", group2, sep = "_")
  comparison_label <- paste0(group1, " vs. ", group2)

  if (method %in% c("plsda", "oplsda")) {

    # Extract VIP scores
    all_vip <- tryCatch({
      ropls::getVipVn(model_results, orthoL = FALSE) %>%
        tibble::enframe(name = "Feature", value = "VIP") %>%
        dplyr::arrange(desc(VIP)) %>%
        as.data.frame()
    }, error = function(e) {
      warning("Could not extract VIP scores: ", e$message)
      return(NULL)
    })

    if (!is.null(all_vip)) {

      results[[paste0("data_VIPScores_", group1, " vs. ", group2)]] <- all_vip

      top_vip <- all_vip %>%
        dplyr::slice_head(n = min(top_features, nrow(all_vip)))

      # Prepare abundance data
      abundance_data <- data_pair %>%
        as.data.frame() %>%
        dplyr::select(dplyr::all_of(top_vip$Feature)) %>%
        dplyr::mutate(Group = factor(y_pair)) %>%
        tidyr::pivot_longer(-Group, names_to = "Feature", values_to = "Abundance") %>%
        dplyr::mutate(Feature = factor(Feature, levels = rev(top_vip$Feature)))

      results[[paste0("data_Abundance_", group1, " vs. ", group2)]] <- abundance_data

      # Create combined VIP and abundance plot
      results <- create_vip_abundance_plot(results, top_vip, abundance_data,
                                           group1, group2, method)
    }
  }

  return(results)
}

# Helper function: Create VIP and abundance plot
create_vip_abundance_plot <- function(results, top_vip, abundance_data,
                                      group1, group2, method) {

  # Check if patchwork is available
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' not available. Creating separate plots.")
    create_separate_plots <- TRUE
  } else {
    create_separate_plots <- FALSE
  }

  # VIP Score Plot
  vip_plot <- ggplot2::ggplot(top_vip, ggplot2::aes(x = VIP, y = stats::reorder(Feature, VIP))) +
    ggplot2::geom_point(size = 4, color = "steelblue") +
    # ggplot2::geom_segment(ggplot2::aes(x = 0, xend = VIP, y = Feature, yend = Feature),
    #                       color = "steelblue", alpha = 0.6) +
    ggplot2::labs(title = paste0("Top ", nrow(top_vip), " Features by VIP Score"),
                  subtitle = paste0(group1, " vs. ", group2),
                  x = "VIP Score", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 12, face = "bold"))

  # Abundance Heatmap
  abundance_plot <- ggplot2::ggplot(abundance_data,
                                    ggplot2::aes(x = Group, y = Feature, fill = Abundance)) +
    ggplot2::geom_tile() + # color = "white", size = 0.1
    ggplot2::scale_fill_gradient(low = "darkgreen", high = "red" #mid = "white",
                                 # ,midpoint = median(abundance_data$Abundance, na.rm = TRUE),
                                 # name = "Abundance"
    ) +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   aspect.ratio = nrow(top_vip) / 3)

  # Combine plots
  if (!create_separate_plots) {
    combined_plot <- vip_plot + abundance_plot +
      patchwork::plot_layout() # widths = c(2, 1)
  } else {
    combined_plot <- list(vip_plot = vip_plot, abundance_plot = abundance_plot)
  }

  results[[paste0("plot_VIPAbundance_", group1, " vs. ", group2)]] <- combined_plot

  return(results)
}

# Helper function: Generate score plots
generate_score_plots <- function(results, model_results, y_pair, group1, group2, method) {

  score_data <- data.frame(
    PC1 = model_results@scoreMN[, 1],
    Group = y_pair
  )

  if (ncol(model_results@scoreMN) > 1) {
    score_data$PC2 <- model_results@scoreMN[, 2]

    score_plot <- ggplot2::ggplot(score_data, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::stat_ellipse(level = 0.95, linetype = "dashed") +
      ggplot2::labs(title = paste0("Score Plot - ", toupper(method)),
                    subtitle = paste0(group1, " vs. ", group2),
                    x = paste0("t[1] (", round(model_results@modelDF$R2X[1] * 100, 1), "%)"),
                    y = paste0("t[2] (", round(model_results@modelDF$R2X[2] * 100, 1), "%)")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  } else {
    score_plot <- ggplot2::ggplot(score_data, ggplot2::aes(x = PC1, y = 0, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.7, position = ggplot2::position_jitter(height = 0.1)) +
      ggplot2::labs(title = paste0("Score Plot - ", toupper(method)),
                    subtitle = paste0(group1, " vs. ", group2),
                    x = paste0("t[1] (", round(model_results@modelDF$R2X[1] * 100, 1), "%)"),
                    y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  }

  results[[paste0("plot_Scores_", group1, " vs. ", group2)]] <- score_plot

  return(results)
}

# Helper function: Generate S-plot (OPLS-DA only)
generate_splot <- function(results, model_results, data_pair, group1, group2) {

  splot_data <- tryCatch({
    tibble::tibble(
      Variable = colnames(data_pair),
      Covariance = purrr::map_dbl(as.data.frame(data_pair),
                                  ~ stats::cov(.x, model_results@scoreMN[, 1], use = "complete.obs")),
      Correlation = purrr::map_dbl(as.data.frame(data_pair),
                                   ~ stats::cor(.x, model_results@scoreMN[, 1], use = "complete.obs"))
    )
  }, error = function(e) {
    warning("Could not generate S-plot data: ", e$message)
    return(NULL)
  })

  if (!is.null(splot_data)) {

    results[[paste0("data_SPlot_", group1, " vs. ", group2)]] <- splot_data

    # Create S-plot
    s_plot <- ggplot2::ggplot(splot_data, ggplot2::aes(x = Covariance, y = Correlation)) +
      ggplot2::geom_point(ggplot2::aes(color = abs(Correlation)), size = 2, alpha = 0.7) +
      ggplot2::scale_color_gradient(low = "lightblue", high = "red", name = "|Correlation|") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste0("S-Plot - OPLS-DA"),
        subtitle = paste0(group1, " vs. ", group2),
        x = "Covariance [p1]",
        y = "Correlation [p(corr)1]"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 10)
      )

    results[[paste0("plot_SPlot_", group1, " vs. ", group2)]] <- s_plot
  }

  return(results)
}

# Helper function: Generate sPLS-DA specific analysis
generate_splsda_analysis <- function(results, model_results, data_pair, y_pair,
                                     group1, group2, top_features, verbose) {

  # comparison_label <- paste0(group1, "vs", group2, sep = "_")
  comparison_label <- paste0(group1, " vs. ", group2)

  # Extract selected features
  selected_features <- tryCatch({
    selected_vars <- mixOmics::selectVar(model_results, comp = 1)
    if (length(selected_vars$name) > 0) {
      data.frame(
        Feature = selected_vars$name,
        Loading = selected_vars$value[, 1],
        stringsAsFactors = FALSE
      ) %>%
        dplyr::arrange(desc(abs(Loading))) %>%
        dplyr::slice_head(n = min(top_features, nrow(.)))
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Could not extract selected features for sPLS-DA: ", e$message)
    return(NULL)
  })

  if (!is.null(selected_features)) {

    results[[paste0("data_SelectedFeatures_", group1, " vs. ", group2)]] <- selected_features

    # Create loading plot
    loading_plot <- ggplot2::ggplot(selected_features,
                                    ggplot2::aes(x = Loading, y = stats::reorder(Feature, abs(Loading)))) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
      ggplot2::labs(title = paste0("Selected Features - sPLS-DA"),
                    subtitle = paste0(group1, " vs. ", group2),
                    x = "Loading Value", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))

    results[[paste0("plot_Loadings_", group1, " vs. ", group2)]] <- loading_plot

    # Create abundance heatmap for selected features
    if (nrow(selected_features) > 0) {
      abundance_data <- data_pair %>%
        as.data.frame() %>%
        dplyr::select(dplyr::all_of(selected_features$Feature)) %>%
        dplyr::mutate(Group = factor(y_pair)) %>%
        tidyr::pivot_longer(-Group, names_to = "Feature", values_to = "Abundance") %>%
        dplyr::mutate(Feature = factor(Feature, levels = rev(selected_features$Feature)))

      abundance_plot <- ggplot2::ggplot(abundance_data,
                                        ggplot2::aes(x = Group, y = Feature, fill = Abundance)) +
        ggplot2::geom_tile(color = "white", size = 0.1) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                      midpoint = median(abundance_data$Abundance, na.rm = TRUE),
                                      name = "Abundance") +
        ggplot2::labs(title = "Feature Abundance Heatmap",
                      subtitle = paste0(group1, " vs. ", group2),
                      x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5))

      results[[paste0("plot_Abundance_", group1, " vs. ", group2)]] <- abundance_plot
    }
  }

  # Generate score plot for sPLS-DA
  score_data <- data.frame(
    Comp1 = model_results$variates$X[, 1],
    Group = y_pair
  )

  if (ncol(model_results$variates$X) > 1) {
    score_data$Comp2 <- model_results$variates$X[, 2]

    score_plot <- ggplot2::ggplot(score_data, ggplot2::aes(x = Comp1, y = Comp2, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::stat_ellipse(level = 0.95, linetype = "dashed") +
      ggplot2::labs(title = "sPLS-DA Score Plot",
                    subtitle = paste0(group1, " vs. ", group2),
                    x = "Component 1", y = "Component 2") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  } else {
    score_plot <- ggplot2::ggplot(score_data, ggplot2::aes(x = Comp1, y = 0, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.7, position = ggplot2::position_jitter(height = 0.1)) +
      ggplot2::labs(title = "sPLS-DA Score Plot",
                    subtitle = paste0(group1, " vs. ", group2),
                    x = "Component 1", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
  }

  results[[paste0("plot_Scores_", group1, " vs. ", group2)]] <- score_plot

  return(results)
}

# Helper function: Generate analysis summary
generate_analysis_summary <- function(results, method, verbose) {

  if (verbose) cat("\nGenerating analysis summary...\n")

  # Count successful analyses
  model_keys <- names(results)[grepl("^results_", names(results))]
  n_comparisons <- length(model_keys)

  # Extract model performance metrics
  performance_metrics <- list()

  for (key in model_keys) {
    model <- results[[key]]
    comparison_name <- gsub("^results_[A-Z]+_", "", key)

    if (method %in% c("pls", "plsda", "oplsda") && inherits(model, "opls")) {
      performance_metrics[[comparison_name]] <- list(
        R2X = if(length(model@modelDF$R2X) > 0) sum(model@modelDF$R2X, na.rm = TRUE) else NA,
        R2Y = if(length(model@modelDF$R2Y) > 0) sum(model@modelDF$R2Y, na.rm = TRUE) else NA,
        Q2 = if(length(model@modelDF$Q2) > 0) sum(model@modelDF$Q2, na.rm = TRUE) else NA,
        n_components = if(length(model@modelDF$R2X) > 0) length(model@modelDF$R2X) else 0,
        n_features = nrow(model@loadingMN)
      )
    } else if (method == "splsda" && inherits(model, "mixo_splsda")) {
      # Extract performance metrics for sPLS-DA if available
      performance_metrics[[comparison_name]] <- list(
        n_components = model$ncomp,
        n_selected_features = if(!is.null(model$keepX)) sum(model$keepX) else NA,
        keepX = model$keepX
      )
    }
  }

  summary_info <- list(
    method = method,
    n_comparisons = n_comparisons,
    comparison_names = gsub("^results_[A-Z]+_", "", model_keys),
    performance_metrics = performance_metrics,
    timestamp = Sys.time()
  )

  if (verbose) {
    cat("Summary:\n")
    cat("  Method:", method, "\n")
    cat("  Comparisons analyzed:", n_comparisons, "\n")
    cat("  Comparison names:", paste(summary_info$comparison_names, collapse = ", "), "\n")
  }

  return(summary_info)
}

#' Print Method for PLS Analysis Results
#'
#' @param x A pls_results object returned by perform_PLS
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @export
print.pls_results <- function(x, ...) {
  cat("PLS Analysis Results\n")
  cat("====================\n")
  cat("Method:", x$method, "\n")
  cat("Number of comparisons:", x$summary$n_comparisons, "\n")
  cat("Timestamp:", as.character(x$summary$timestamp), "\n\n")

  if (length(x$summary$comparison_names) > 0) {
    cat("Comparisons:\n")
    for (i in seq_along(x$summary$comparison_names)) {
      cat("  ", i, ".", x$summary$comparison_names[i], "\n")
    }
  }

  cat("\nAvailable result components:\n")
  result_types <- c("results_", "data_", "plot_")
  for (type in result_types) {
    matching_keys <- names(x)[grepl(paste0("^", type), names(x))]
    if (length(matching_keys) > 0) {
      cat("  ", gsub("_$", "", type), ":", length(matching_keys), "items\n")
    }
  }

  invisible(x)
}

# #' Plot Method for PLS Analysis Results
# #'
# #' @param x A pls_results object returned by perform_PLS
# #' @param comparison Integer or character specifying which comparison to plot (default: 1)
# #' @param type Character specifying plot type: "scores", "vipabundance", "splot", "loadings", "abundance"
# #' @param ... Additional arguments (not used)
# #' @return A ggplot object
# #' @export
# plot.pls_results <- function(x, comparison = 1, type = "scores", ...) {
#
#   if (length(x$summary$comparison_names) == 0) {
#     stop("No comparisons found in results")
#   }
#
#   if (is.numeric(comparison)) {
#     if (comparison > length(x$summary$comparison_names)) {
#       stop("Comparison index out of range")
#     }
#     comp_name <- x$summary$comparison_names[comparison]
#   } else {
#     comp_name <- comparison
#   }
#
#   # Find the appropriate plot
#   plot_key <- paste0("plot_", stringr::str_to_title(type), "_", comp_name)
#
#   if (!plot_key %in% names(x)) {
#     available_plots <- names(x)[grepl(paste0("_", comp_name, "$"), names(x)) & grepl("^plot_", names(x))]
#     stop("Plot not found. Available plots for this comparison: ",
#          paste(gsub("^plot_|_.*$", "", available_plots), collapse = ", "))
#   }
#
#   return(x[[plot_key]])
# }
