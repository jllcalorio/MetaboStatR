#' Perform Principal Component Analysis (PCA)
#'
#' @description
#' This function performs PCA and generates multiple scores plots at once based on a list of
#' principal component combinations, plus an optional scree plot. Each combination in the list
#' will produce a separate scores plot.
#'
#' @param data List. This list must be a result from the `perform_PreprocessingPeakData` function.
#'   The function will automatically detect if replicates were merged during preprocessing
#'   and use the appropriate data matrix.
#' @param scorePC List. A list of 2-element numeric vectors, where each vector specifies
#'   the principal components to plot. For example: list(c(1,2), c(1,3), c(2,3))
#'   will generate 3 plots: PC1 vs PC2, PC1 vs PC3, and PC2 vs PC3.
#' @param includeQC Boolean. If `TRUE`, includes QC (Quality Control) samples in the analysis and plots.
#'   If `FALSE`, uses only biological samples (BS). Defaults to `FALSE`.
#' @param arrangeLevels Vector. Determines how the groups will be arranged.
#'   The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which
#'   sorts the groups in alphabetical order.
#' @param scoresEllipse Boolean. If `TRUE` (default), adds an ellipse in the scores plot.
#' @param scoresColorVar String. Optional. Column name from metadata to use for
#'   continuous gradient coloring in scores plots (e.g., "Time", "Concentration").
#'   If `NULL` (default), uses discrete `Group` colors instead.
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
#' @param showOutliers Boolean. If `TRUE`, labels samples that fall outside the 95\%
#'   confidence ellipse of their respective group directly on the scores plot(s).
#'   Outlier detection exactly mirrors \code{ggplot2::stat_ellipse}: the ellipse is
#'   computed from the group covariance matrix scaled by a t-distribution multiplier
#'   (\code{qt(0.975, df = n - 1)^2}), and a point is flagged as an outlier when its
#'   Mahalanobis distance (scaled by the same factor) exceeds 1. Requires at least
#'   3 samples per group (the same minimum needed to draw an ellipse). Defaults to
#'   \code{FALSE}. Note: has no effect when \code{scoresColorVar} is supplied, because
#'   ellipses are not drawn in gradient-color mode.
#'   Labels are rendered with \pkg{ggrepel} to prevent overlap; connecting segments
#'   are drawn from each label to its corresponding point.
#' @param outlierGroups Character vector controlling which groups are checked for
#'   outliers and labelled on the plot. Accepted values are:
#'   \itemize{
#'     \item \code{"all"} — every group present in the data (default).
#'     \item Any group name(s) found in \code{Metadata$Group} (or
#'       \code{Metadata_merged$Group} when replicates were merged), e.g.
#'       \code{c("GroupA", "GroupB")}.
#'     \item \code{"QC"} — a convenience alias that selects all QC-type groups
#'       (\code{"SQC"}, \code{"EQC"}, \code{"QC"}) that are present in the data.
#'       Note that QC samples must also be included via \code{includeQC = TRUE}
#'       for this to have any effect.
#'   }
#'   Multiple values can be combined, e.g. \code{c("GroupA", "QC")}. Any value not
#'   recognised as a group name present in the data will trigger a warning.
#'   Defaults to \code{"all"}.
#'
#' @returns Returns a list containing:
#'   \itemize{
#'     \item \code{plots}: A list of ggplot objects, one for each PC combination.
#'     \item \code{scree_plot}: A ggplot object for the scree plot (if \code{includeScree}
#'       is \code{TRUE}).
#'     \item \code{plot_info}: A data frame with information about each plot (PC
#'       combinations, variance explained).
#'     \item \code{pca_results}: The PCA results object from \code{stats::prcomp}.
#'     \item \code{data_used}: The data matrix used for PCA.
#'     \item \code{variance_explained}: Data frame of variance explained (\%) per PC.
#'     \item \code{eigenvalues}: Data frame of eigenvalues per PC.
#'     \item \code{scores_data}: A list of data frames, one per PC combination (named
#'       identically to \code{plots}). Each data frame contains:
#'       \itemize{
#'         \item \code{Label} — sample name / row name.
#'         \item \code{Group} — group assignment.
#'         \item \code{Batch} — batch assignment.
#'         \item \code{SampleType} — \code{"QC"} or \code{"Biological"}.
#'         \item \code{PC<i>}, \code{PC<j>} — scores on the two plotted PCs.
#'         \item \code{mahal_dist} — Mahalanobis distance from the group centroid
#'           in the 2-D PC space (groups with < 3 samples receive \code{NA}).
#'         \item \code{ellipse_threshold} — the t-based scaling factor used as the
#'           95\% CI boundary, matching \code{ggplot2::stat_ellipse} exactly
#'           (\code{NA} for groups with < 3 samples).
#'         \item \code{is_outlier} — logical; \code{TRUE} when the point lies
#'           outside the 95\% ellipse as drawn by \code{stat_ellipse}.
#'       }
#'   }
#' @export
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom stats complete.cases mahalanobis cov qt
#' @importFrom utils getFromNamespace
#'
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. (for prcomp)
#' @references Mardia, K. V., J. T. Kent, and J. M. Bibby (1979) Multivariate Analysis, London: Academic Press. (for prcomp)
#' @references Venables, W. N. and B. D. Ripley (2002) Modern Applied Statistics with S, Springer-Verlag. (for prcomp)
#'
#' @examples
#' \dontrun{
#' # Standard PCA with discrete groups
#' multi_plots <- perform_PCA(
#'   data = data_from_perform_PreprocessingPeakData_function,
#'   scorePC = list(c(1,2), c(1,3)),
#'   includeQC = FALSE
#' )
#'
#' # PCA with outlier labeling for all groups
#' multi_plots_outliers <- perform_PCA(
#'   data = data_from_perform_PreprocessingPeakData_function,
#'   scorePC = list(c(1,2)),
#'   showOutliers = TRUE,
#'   outlierGroups = "all"
#' )
#'
#' # Label outliers only in specific groups
#' multi_plots_sel <- perform_PCA(
#'   data = data_from_perform_PreprocessingPeakData_function,
#'   scorePC = list(c(1,2)),
#'   showOutliers = TRUE,
#'   outlierGroups = c("GroupA", "GroupB")
#' )
#'
#' # PCA with gradient coloring by a continuous variable (e.g., "Time")
#' multi_plots_time <- perform_PCA(
#'   data = data_from_perform_PreprocessingPeakData_function,
#'   scorePC = list(c(1,2)),
#'   includeQC = TRUE,
#'   scoresColorVar = "Time",
#'   scoresLegend = "Collection Time"
#' )
#'
#' plot_PCAScreeScores(multi_plots_time)
#'
#' # Inspect per-sample scores and outlier flags
#' head(multi_plots_outliers$scores_data[["PC1_vs_PC2"]])
#' }

perform_PCA <- function(
    data,
    scorePC = list(c(1, 2), c(1, 3), c(2, 3)),
    includeQC = FALSE,
    arrangeLevels = NULL,
    scoresEllipse = TRUE,
    scoresColorVar = NULL,
    scoresTitle = NULL,
    scoresLegend = NULL,
    includeScree = TRUE,
    screeWhat = c("variance", "eigenvalue")[1],
    screeType = c("bar", "line", "both")[3],
    screeTitle = "Scree Plot",
    screeMaxPC = NULL,
    screeShowValues = TRUE,
    screeShowCumulative = FALSE,
    showOutliers = FALSE,
    outlierGroups = "all"
) {

  # ============================================================================
  # PARAMETER VALIDATION
  # ============================================================================

  if (data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("The parameter 'data' must be from the 'perform_PreprocessingPeakData' function.")
  }

  if (!is.list(scorePC)) {
    stop("scorePC: Must be a list of numeric vectors. Example: list(c(1,2), c(1,3), c(2,3))")
  }
  if (length(scorePC) == 0) {
    stop("scorePC: List cannot be empty. Provide at least one PC combination.")
  }

  for (i in seq_along(scorePC)) {
    pc_combo <- scorePC[[i]]
    if (!is.numeric(pc_combo))
      stop(paste0("scorePC element ", i, ": Must be numeric. Found: ", class(pc_combo)[1]))
    if (length(pc_combo) != 2)
      stop(paste0("scorePC element ", i, ": Must contain exactly 2 values. Found: ", length(pc_combo)))
    if (any(pc_combo <= 0) || any(pc_combo != floor(pc_combo)))
      stop(paste0("scorePC element ", i, ": Values must be positive integers. Found: c(",
                  paste(pc_combo, collapse = ", "), ")"))
    if (pc_combo[1] == pc_combo[2])
      stop(paste0("scorePC element ", i, ": Cannot plot the same PC against itself. Found: c(",
                  paste(pc_combo, collapse = ", "), ")"))
  }

  for (nm in c("includeQC", "scoresEllipse", "includeScree", "showOutliers")) {
    val <- get(nm)
    if (!is.logical(val) || length(val) != 1)
      stop(paste0(nm, ": Must be a single logical value (TRUE or FALSE)."))
  }

  if (!is.null(arrangeLevels) && (!is.vector(arrangeLevels) || !is.character(arrangeLevels)))
    stop("arrangeLevels: Must be NULL or a character vector of group names.")

  if (!is.null(scoresTitle)) {
    if (!is.character(scoresTitle))
      stop("scoresTitle: Must be NULL, a single character string, or a character vector.")
    if (length(scoresTitle) > 1 && length(scoresTitle) != length(scorePC))
      stop(paste0("scoresTitle: Expected ", length(scorePC), " titles, found ", length(scoresTitle), "."))
  }

  if (!is.null(scoresLegend) && (!is.character(scoresLegend) || length(scoresLegend) != 1))
    stop("scoresLegend: Must be NULL or a single character string.")

  if (!is.character(screeWhat) || length(screeWhat) != 1 || !screeWhat %in% c("variance", "eigenvalue"))
    stop("screeWhat: Must be either 'variance' or 'eigenvalue'.")
  if (!is.character(screeType) || length(screeType) != 1 || !screeType %in% c("bar", "line", "both"))
    stop("screeType: Must be 'bar', 'line', or 'both'.")
  if (!is.character(screeTitle) || length(screeTitle) != 1)
    stop("screeTitle: Must be a single character string.")
  if (!is.null(screeMaxPC) &&
      (!is.numeric(screeMaxPC) || length(screeMaxPC) != 1 || screeMaxPC <= 0 || screeMaxPC != floor(screeMaxPC)))
    stop("screeMaxPC: Must be NULL or a single positive integer.")
  if (!is.logical(screeShowValues)    || length(screeShowValues) != 1)
    stop("screeShowValues: Must be a single logical value (TRUE or FALSE).")
  if (!is.logical(screeShowCumulative) || length(screeShowCumulative) != 1)
    stop("screeShowCumulative: Must be a single logical value (TRUE or FALSE).")

  if (!is.character(outlierGroups) || length(outlierGroups) == 0)
    stop("outlierGroups: Must be a non-empty character vector. ",
         "Use \"all\", \"QC\", or specific group names from Metadata$Group.")

  if (screeShowCumulative && screeWhat != "variance") {
    warning("screeShowCumulative only works with screeWhat='variance'. Setting screeShowCumulative to FALSE.")
    screeShowCumulative <- FALSE
  }

  if (showOutliers && !is.null(scoresColorVar))
    warning("showOutliers = TRUE has no effect when scoresColorVar is supplied, ",
            "because ellipses are not drawn in gradient-color mode.")

  # Check for ggrepel when showOutliers is TRUE
  if (showOutliers && is.null(scoresColorVar) &&
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for showOutliers = TRUE. ",
         "Install it with: install.packages('ggrepel')")
  }

  if (is.null(data$Parameters$merge_replicates)) {
    merge_replicates <- FALSE
    warning("merge_replicates parameter not found in data. Using default FALSE for backward compatibility.")
  } else {
    merge_replicates <- data$Parameters$merge_replicates
    if (!is.logical(merge_replicates) || length(merge_replicates) != 1)
      stop("data$Parameters$merge_replicates must be a single logical value (TRUE or FALSE).")
  }

  # ============================================================================
  # DATA PREPARATION AND PCA
  # ============================================================================

  if (merge_replicates) {
    if (is.null(data$data_scaledNONPLS_merged) || length(data$data_scaledNONPLS_merged) == 0)
      stop("merge_replicates is TRUE but data$data_scaledNONPLS_merged does not exist or is empty.")
    if (is.null(data$Metadata_merged) || nrow(data$Metadata_merged) == 0)
      stop("merge_replicates is TRUE but data$Metadata_merged does not exist or is empty.")
    if (nrow(data$data_scaledNONPLS_merged) != nrow(data$Metadata_merged))
      stop("Dimension mismatch: data_scaledNONPLS_merged and Metadata_merged row counts differ.")

    full_data     <- data$data_scaledNONPLS_merged
    metadata_full <- data$Metadata_merged
    data_source   <- "merged"
  } else {
    if (is.null(data$data_scaledNONPLS_varFiltered) || length(data$data_scaledNONPLS_varFiltered) == 0)
      stop("The preprocessed data does not exist or does not have any data on it.")

    full_data     <- data$data_scaledNONPLS_varFiltered
    metadata_full <- data$Metadata
    data_source   <- "filtered"
  }

  qc_indices     <- metadata_full$Group %in% c("SQC", "EQC", "QC")
  non_qc_indices <- !qc_indices

  if (includeQC) {
    sample_indices <- rep(TRUE, nrow(full_data))
    sample_type    <- paste0("All samples (QC + Biological) - ", data_source, " data")
  } else {
    sample_indices <- non_qc_indices
    sample_type    <- paste0("Biological samples only - ", data_source, " data")
  }

  if (sum(sample_indices) == 0)
    stop("No samples available for analysis with the current settings.")

  data_pca <- full_data[sample_indices, , drop = FALSE]

  pca_res <- tryCatch(
    stats::prcomp(data_pca, center = FALSE, scale. = FALSE),
    error = function(e) stop(paste("PCA failed:", e$message))
  )

  max_pc <- max(unlist(scorePC))
  if (max_pc > ncol(pca_res$x))
    stop(paste0("Requested PC", max_pc, " but only ", ncol(pca_res$x), " PCs are available."))

  eigenvalues        <- pca_res$sdev^2
  variance_explained <- round(100 * (eigenvalues / sum(eigenvalues)), 1)

  # ============================================================================
  # BASE SCORES DATA FRAME
  # ============================================================================

  metadata_subset <- metadata_full[sample_indices, , drop = FALSE]

  if (includeQC) {
    qc_subset       <- qc_indices[sample_indices]
    combined_groups <- ifelse(qc_subset,
                              paste0("QC_", metadata_subset$Group),
                              as.character(metadata_subset$Group))
    base_scores_df <- data.frame(
      Group      = combined_groups,
      Batch      = metadata_subset$Batch,
      Label      = rownames(data_pca),
      SampleType = ifelse(qc_subset, "QC", "Biological"),
      stringsAsFactors = FALSE
    )
  } else {
    base_scores_df <- data.frame(
      Group      = metadata_subset$Group,
      Batch      = metadata_subset$Batch,
      Label      = rownames(data_pca),
      SampleType = "Biological",
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(scoresColorVar)) {
    if (!scoresColorVar %in% colnames(metadata_subset))
      stop(paste0("scoresColorVar '", scoresColorVar, "' not found in metadata."))
    base_scores_df$ColorVar <- metadata_subset[[scoresColorVar]]
  }

  all_groups <- unique(base_scores_df$Group)

  if (!is.null(arrangeLevels)) {
    missing_grps <- arrangeLevels[!arrangeLevels %in% all_groups]
    if (length(missing_grps) > 0)
      warning("The following groups in arrangeLevels were not found in data: ",
              paste(missing_grps, collapse = ", "))
    valid_levels     <- arrangeLevels[arrangeLevels %in% all_groups]
    remaining_groups <- all_groups[!all_groups %in% valid_levels]
    base_scores_df$Group <- factor(base_scores_df$Group,
                                   levels = c(valid_levels, sort(remaining_groups)))
  } else {
    base_scores_df$Group <- factor(base_scores_df$Group)
  }

  # ============================================================================
  # VALIDATE outlierGroups
  # ============================================================================

  qc_type_groups <- c("SQC", "EQC", "QC")

  if (showOutliers && is.null(scoresColorVar)) {
    qc_alias_used <- "QC" %in% outlierGroups
    non_alias     <- outlierGroups[outlierGroups != "QC" & outlierGroups != "all"]
    resolved_outlier_groups <- character(0)

    if ("all" %in% outlierGroups) {
      resolved_outlier_groups <- levels(base_scores_df$Group)
    } else {
      if (qc_alias_used) {
        qc_matched <- levels(base_scores_df$Group)[
          grepl(paste(qc_type_groups, collapse = "|"), levels(base_scores_df$Group))
        ]
        resolved_outlier_groups <- c(resolved_outlier_groups, qc_matched)
      }
      if (length(non_alias) > 0) {
        unrecognised <- non_alias[!non_alias %in% levels(base_scores_df$Group)]
        if (length(unrecognised) > 0)
          warning("outlierGroups: The following values were not found and will be ignored: ",
                  paste(unrecognised, collapse = ", "),
                  ".\nAvailable groups: ", paste(levels(base_scores_df$Group), collapse = ", "))
        resolved_outlier_groups <- c(resolved_outlier_groups,
                                     non_alias[non_alias %in% levels(base_scores_df$Group)])
      }
      resolved_outlier_groups <- unique(resolved_outlier_groups)
    }

    if (length(resolved_outlier_groups) == 0) {
      warning("outlierGroups: No valid groups resolved. Outlier labeling will be skipped.")
      showOutliers <- FALSE
    }
  } else {
    resolved_outlier_groups <- character(0)
  }

  # ============================================================================
  # HELPER: replicate ggplot2::stat_ellipse outlier detection (t-distribution method)
  #
  # ggplot2 computes the confidence ellipse as:
  #   centre + sqrt(t_crit) * chol(cov(x)) * unit_circle
  # where t_crit = qt(level/2 + 0.5, df = n - 1)^2  (two-sided, df = n-1)
  #
  # A point p is INSIDE the ellipse when its scaled Mahalanobis distance
  #   d^2 = (p - mu) %*% solve(cov(x)) %*% (p - mu)
  # satisfies d^2 <= t_crit.
  # Points where d^2 > t_crit are OUTSIDE = outliers.
  # ============================================================================

  .compute_outliers <- function(df, pc1_col, pc2_col) {
    mahal_d   <- rep(NA_real_, nrow(df))
    threshold <- rep(NA_real_, nrow(df))
    is_out    <- rep(NA, nrow(df))

    for (grp in levels(df$Group)) {
      idx <- which(df$Group == grp)
      n   <- length(idx)
      if (n < 3) next   # need >= 3 to compute a 2x2 covariance matrix

      pts     <- as.matrix(df[idx, c(pc1_col, pc2_col), drop = FALSE])
      cov_mat <- tryCatch(stats::cov(pts), error = function(e) NULL)
      if (is.null(cov_mat) || det(cov_mat) <= .Machine$double.eps) next

      # t-based threshold — matches ggplot2::stat_ellipse(level = 0.95, type = "t")
      t_crit <- stats::qt(0.975, df = n - 1)^2   # = qt((0.95/2 + 0.5), df = n-1)^2

      d2 <- tryCatch(
        stats::mahalanobis(pts, colMeans(pts), cov_mat),
        error = function(e) rep(NA_real_, n)
      )

      mahal_d[idx]   <- d2
      threshold[idx] <- t_crit
      is_out[idx]    <- !is.na(d2) & (d2 > t_crit)
    }

    df$mahal_dist        <- mahal_d
    df$ellipse_threshold <- threshold
    df$is_outlier        <- is_out
    df
  }

  # ============================================================================
  # SCREE PLOT
  # ============================================================================

  scree_plot <- NULL
  if (includeScree) {
    n_pcs    <- ncol(pca_res$x)
    max_show <- if (is.null(screeMaxPC)) min(n_pcs, 20) else min(screeMaxPC, n_pcs)

    y_values <- if (screeWhat == "variance") variance_explained[1:max_show] else eigenvalues[1:max_show]
    y_label  <- if (screeWhat == "variance") "Variance Explained (%)" else "Eigenvalue"

    scree_data <- data.frame(
      PC       = factor(paste0("PC", 1:max_show), levels = paste0("PC", 1:max_show)),
      PC_num   = 1:max_show,
      Value    = y_values,
      CumValue = if (screeWhat == "variance") cumsum(y_values) else NA,
      stringsAsFactors = FALSE
    )

    scree_plot <- ggplot2::ggplot(scree_data, ggplot2::aes(x = PC_num, y = Value))

    if (screeType == "bar") {
      scree_plot <- scree_plot +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
        ggplot2::scale_x_continuous(breaks = 1:max_show, labels = paste0("PC", 1:max_show))
    } else if (screeType == "line") {
      scree_plot <- scree_plot +
        ggplot2::geom_line(color = "steelblue", linewidth = 1, group = 1) +
        ggplot2::geom_point(color = "steelblue", size = 3) +
        ggplot2::scale_x_continuous(breaks = 1:max_show, labels = paste0("PC", 1:max_show))
    } else {
      scree_plot <- scree_plot +
        ggplot2::geom_col(fill = "steelblue", alpha = 0.5, width = 0.6) +
        ggplot2::geom_line(color = "darkblue", linewidth = 1, group = 1) +
        ggplot2::geom_point(color = "darkblue", size = 3) +
        ggplot2::scale_x_continuous(breaks = 1:max_show, labels = paste0("PC", 1:max_show))
    }

    if (screeShowValues) {
      val_labs    <- if (screeWhat == "variance") paste0(round(scree_data$Value, 1), "%") else round(scree_data$Value, 1)
      lbl_hjust   <- if (screeType == "line") 0.5  else -0.2
      lbl_vjust   <- if (screeType == "line") -0.5 else -0.5
      scree_plot  <- scree_plot +
        ggplot2::geom_text(ggplot2::aes(label = val_labs),
                           vjust = lbl_vjust, hjust = lbl_hjust, size = 3, color = "black")
    }

    if (screeShowCumulative && screeWhat == "variance") {
      cum_labs  <- paste0(round(scree_data$CumValue, 1), "%")
      cum_hjust <- if (screeType == "line") 0.5  else -0.2
      cum_vjust <- if (screeType == "line") { if (screeShowValues) -1.2 else -0.5 } else -1.8
      scree_plot <- scree_plot +
        ggplot2::geom_text(ggplot2::aes(label = cum_labs),
                           vjust = cum_vjust, hjust = cum_hjust,
                           size = 2.8, color = "darkred", fontface = "italic")
    }

    y_ax_text  <- if (screeShowValues) ggplot2::element_blank() else ggplot2::element_text(size = 10)
    y_ax_ticks <- if (screeShowValues) ggplot2::element_blank() else ggplot2::element_line()

    scree_plot <- scree_plot +
      ggplot2::labs(x = "Principal Component", y = y_label, title = screeTitle) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.title       = ggplot2::element_text(size = 12),
        axis.text.x      = ggplot2::element_text(angle = 90, hjust = 0.5, size = 10),
        axis.text.y      = y_ax_text,
        axis.ticks.y     = y_ax_ticks,
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  # ============================================================================
  # SCORES PLOTS
  # ============================================================================

  plots_list <- list()
  scores_data <- list()
  plot_info <- data.frame(
    plot_number = integer(), PC1 = integer(), PC2 = integer(),
    PC1_variance = numeric(), PC2_variance = numeric(), title = character(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(scorePC)) {
    pc_combo <- scorePC[[i]]
    pc1 <- pc_combo[1];  pc2 <- pc_combo[2]
    pc1_col <- paste0("PC", pc1);  pc2_col <- paste0("PC", pc2)

    df_scores            <- base_scores_df
    df_scores[[pc1_col]] <- pca_res$x[, pc1]
    df_scores[[pc2_col]] <- pca_res$x[, pc2]

    # Always compute outlier stats (available in scores_data regardless of showOutliers)
    df_scores <- .compute_outliers(df_scores, pc1_col, pc2_col)

    x_label    <- paste0("PC", pc1, " (", variance_explained[pc1], "%)")
    y_label    <- paste0("PC", pc2, " (", variance_explained[pc2], "%)")
    plot_title <- if (is.null(scoresTitle)) {
      paste0("PCA Scores Plot (PC", pc1, " vs PC", pc2, ")")
    } else if (length(scoresTitle) == 1) {
      scoresTitle
    } else {
      scoresTitle[i]
    }

    # Base plot
    if (is.null(scoresColorVar)) {
      scores_plot <- ggplot2::ggplot(
        df_scores,
        ggplot2::aes(x = .data[[pc1_col]], y = .data[[pc2_col]],
                     color = Group, shape = Group)
      ) + ggplot2::geom_point(size = 3, alpha = 0.8)
    } else {
      scores_plot <- ggplot2::ggplot(
        df_scores,
        ggplot2::aes(x = .data[[pc1_col]], y = .data[[pc2_col]],
                     color = ColorVar, shape = Group)
      ) +
        ggplot2::geom_point(size = 3, alpha = 0.9) +
        viridis::scale_color_viridis(option = "plasma", name = scoresLegend)
    }

    # Ellipses
    if (scoresEllipse && is.null(scoresColorVar) && nlevels(df_scores$Group) > 1) {
      grp_counts      <- table(df_scores$Group)
      grps_for_ellipse <- names(grp_counts)[grp_counts >= 3]
      if (length(grps_for_ellipse) > 0) {
        df_ell <- df_scores[df_scores$Group %in% grps_for_ellipse, ]
        scores_plot <- scores_plot +
          ggplot2::stat_ellipse(
            data    = df_ell,
            ggplot2::aes(group = Group, color = Group, fill = Group),
            level   = 0.95,
            geom    = "polygon",
            alpha   = 0.2,
            show.legend = FALSE
          )
      }
    }

    # Outlier labels via ggrepel
    if (showOutliers && is.null(scoresColorVar) && length(resolved_outlier_groups) > 0) {
      df_out_label <- df_scores[
        !is.na(df_scores$is_outlier) &
          df_scores$is_outlier &
          as.character(df_scores$Group) %in% resolved_outlier_groups,
      ]

      if (nrow(df_out_label) > 0) {
        scores_plot <- scores_plot +
          ggrepel::geom_label_repel(
            data              = df_out_label,
            ggplot2::aes(
              x     = .data[[pc1_col]],
              y     = .data[[pc2_col]],
              label = Label,
              color = Group          # match group colour for the label border
            ),
            fill              = "white",
            size              = 3,
            label.padding     = ggplot2::unit(0.15, "lines"),
            box.padding       = ggplot2::unit(0.5,  "lines"),
            point.padding     = ggplot2::unit(0.3,  "lines"),
            segment.color     = "grey40",
            # segment.linewidth = 0.5,
            segment.linetype  = "solid",
            min.segment.length = 0,          # always draw the segment
            max.overlaps      = Inf,         # never silently drop labels
            force             = 3,           # repulsion strength
            force_pull        = 0.5,
            show.legend       = FALSE
          )
      }
    }

    scores_plot <- scores_plot +
      ggplot2::labs(x = x_label, y = y_label, title = plot_title) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position  = "bottom",
        plot.title       = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.title       = ggplot2::element_text(size = 12),
        axis.text        = ggplot2::element_text(size = 10),
        legend.title     = ggplot2::element_text(size = 11),
        legend.text      = ggplot2::element_text(size = 10)
      )

    plots_list[[i]]  <- scores_plot
    scores_data[[i]] <- df_scores

    plot_info <- rbind(plot_info, data.frame(
      plot_number  = i, PC1 = pc1, PC2 = pc2,
      PC1_variance = variance_explained[pc1],
      PC2_variance = variance_explained[pc2],
      title        = plot_title,
      stringsAsFactors = FALSE
    ))
  }

  plot_names         <- paste0("PC", plot_info$PC1, "_vs_PC", plot_info$PC2)
  names(plots_list)  <- plot_names
  names(scores_data) <- plot_names

  # ============================================================================
  # RETURN
  # ============================================================================

  results <- list(
    data_source        = data_source,
    plots              = plots_list,
    scree_plot         = scree_plot,
    plot_info          = plot_info,
    pca_results        = pca_res,
    data_used          = data_pca,
    sample_type        = sample_type,
    variance_explained = as.data.frame(variance_explained),
    eigenvalues        = as.data.frame(eigenvalues),
    scores_data        = scores_data,
    function_call      = match.call(),
    scree_included     = includeScree
  )

  class(results) <- c("perform_PCA", "list")
  return(results)
}

# ==============================================================================
# S3 METHODS
# ==============================================================================

#' @export
print.perform_PCA <- function(x, ...) {
  cat("=== Principal Component Analysis ===\n")
  cat("Source Data:      ", x$data_source, "\n")
  cat("Sample Type:      ", x$sample_type, "\n")
  cat("Samples:          ", nrow(x$data_used), "\n")
  cat("Total Variance:   ", sum(x$variance_explained$variance_explained), "%\n")
  cat("Plots Generated:  ", length(x$plots), "\n")
  cat("\nPC combinations:\n")
  for (i in seq_len(nrow(x$plot_info))) {
    cat(sprintf("  Plot %d: PC%d vs PC%d (%.1f%% vs %.1f%% variance)\n",
                i, x$plot_info$PC1[i], x$plot_info$PC2[i],
                x$plot_info$PC1_variance[i], x$plot_info$PC2_variance[i]))
  }
  cat("\nUse $plots[[i]] to access individual scores plots\n")
  if (x$scree_included && !is.null(x$scree_plot))
    cat("Use $scree_plot to access the scree plot\n")
  cat("Use $scores_data[[i]] for per-sample scores, Mahalanobis distances, and outlier flags\n")
  cat("Use $plot_info for detailed plot information\n")
  cat("Use plot_PCAScreeScores() to display all plots at once\n")
  invisible(x)
}

#' @export
summary.perform_PCA <- function(object, ...) {
  var_summ <- head(object$variance_explained, 5)
  ans <- list(
    variance_table      = var_summ,
    plot_info           = object$plot_info,
    cumulative_variance = cumsum(object$variance_explained)[1:5, ]
  )
  class(ans) <- "summary.perform_PCA"
  return(ans)
}

#' @export
print.summary.perform_PCA <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("PCA Variance Explained (Top 5 PCs)\n")
  cat("---------------------------------------\n")
  print(x$variance_table)
  cat("\n-- Plot Configurations --\n")
  print(x$plot_info[, c("plot_number", "title")])
  invisible(x)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Display All Plots from a perform_PCA Object
#'
#' @param multi_plots_obj Object returned from \code{perform_PCA()}.
#' @param arrange Logical. If \code{TRUE} (default), arranges plots in a grid via
#'   \pkg{patchwork}. Falls back to individual printing if unavailable.
#' @param include_scree Logical. Whether to include the scree plot. Default: \code{TRUE}.
#' @param ncol Number of columns in the plot grid. Default: \code{2}.
#' @param nrow Number of rows (default: \code{NULL}, auto-calculated).
#' @export
plot_PCAScreeScores <- function(multi_plots_obj, arrange = TRUE,
                                include_scree = TRUE, ncol = 2, nrow = NULL) {
  if (!inherits(multi_plots_obj, "perform_PCA"))
    stop("Object must be created by perform_PCA() function")

  plots_to_show <- multi_plots_obj$plots
  if (include_scree && multi_plots_obj$scree_included && !is.null(multi_plots_obj$scree_plot))
    plots_to_show <- c(list(scree = multi_plots_obj$scree_plot), plots_to_show)

  if (arrange && requireNamespace("patchwork", quietly = TRUE)) {
    wrap_plots_fn <- getFromNamespace("wrap_plots", "patchwork")
    print(wrap_plots_fn(plots_to_show, ncol = ncol, nrow = nrow))
  } else {
    if (arrange)
      cat("Note: 'patchwork' not available. Displaying plots individually.\n",
          "Install with: install.packages('patchwork')\n\n")
    for (i in seq_along(plots_to_show)) {
      nm <- names(plots_to_show)[i]
      if (is.null(nm) || nm == "") nm <- paste("Plot", i)
      cat("Plot:", nm, "\n")
      print(plots_to_show[[i]])
      if (i < length(plots_to_show)) cat("\n", strrep("-", 50), "\n\n")
    }
  }
}

#' Extract Scree Plot Data
#'
#' @param multi_plots_obj Object returned from \code{perform_PCA()}.
#' @param what Character. \code{"variance"} or \code{"eigenvalue"}.
#' @param max_pc Integer. Maximum number of PCs. If \code{NULL}, includes all.
#' @return Data frame with PC-level summary statistics.
#' @export
get_ScreeData <- function(multi_plots_obj, what = "variance", max_pc = NULL) {
  if (!inherits(multi_plots_obj, "perform_PCA"))
    stop("Object must be created by perform_PCA() function")

  n_pcs    <- nrow(multi_plots_obj$variance_explained)
  max_show <- if (is.null(max_pc)) n_pcs else min(max_pc, n_pcs)

  if (what == "variance") {
    values     <- multi_plots_obj$variance_explained[1:max_show, ]
    value_name <- "Variance_Explained_Percent"
  } else if (what == "eigenvalue") {
    values     <- multi_plots_obj$eigenvalues[1:max_show, ]
    value_name <- "Eigenvalue"
  } else {
    stop("'what' must be either 'variance' or 'eigenvalue'")
  }

  result <- data.frame(PC = paste0("PC", 1:max_show), PC_Number = 1:max_show,
                       stringsAsFactors = FALSE)
  result[[value_name]] <- values
  if (what == "variance") result$Cumulative_Variance <- cumsum(values)
  return(result)
}
