#' Comprehensive Data Preprocessing Pipeline
#'
#' @description
#' Performs a complete data preprocessing workflow to prepare raw data for
#' downstream analysis. This function applies preprocessing steps sequentially
#' in the order specified by the parameters to ensure optimal data quality
#' and analytical readiness. This function also automatically merges
#' technical replicates.
#'
#' @param raw_data List. Quality-checked data from the `perform_DataQualityCheck` function.
#' @param outliers Vector. Biological samples and/or QC samples considered as outliers. Example format: `c('Sample1', 'Sample2', 'QC1', 'QC2', ...)`. Defaults to `NULL`.
#' @param filterMissing Numeric. Minimum percentage of missing values across all groups required to remove a feature. Set to `NULL` to skip this step.
#' @param filterMissing_by_group Boolean. Determines whether `filterMissing` should assess group-specific missingness before feature removal.
#' @param filterMissing_includeQC Boolean. If `FALSE` (default), QC samples are excluded when implementing `filterMissing`.
#' @param denMissing Numeric. Denominator value used in the fraction `1/denMissing` to replace missing values.
#' @param driftBatchCorrection Boolean. If `TRUE` (default), perform Quality Control-Robust Spline Correction (QC-RSC) algorithm for signal drift and batch effect correction.
#' @param spline_smooth_param Numeric. Spline smoothing parameter ranging from 0 to 1.
#' @param spline_smooth_param_limit Vector. A vector of format `c(min, max)` for spline parameter limits.
#' @param log_scale Boolean. If `TRUE` (default), performs signal correction fit on log-scaled data.
#' @param min_QC Numeric. Minimum number of QC samples required for signal correction per batch.
#' @param removeUncorrectedFeatures Boolean. If `TRUE` (default), removes features that were not corrected by QCRSC due to insufficient QC samples meeting the min_QC threshold.
#' @param dataNormalize String. Data normalization method. Options:
#'  \itemize{
#'    \item \code{"none"}: No normalization
#'    \item \code{"Normalization"}: Using the values from "Normalization" row
#'    \item \code{"sum"}: By sum
#'    \item \code{"median"}: By median
#'    \item \code{"PQN1"}: By median of reference spectrum
#'    \item \code{"PQN2"}: By reference sample supplied in `refSample`
#'    \item \code{"groupPQN"}: By group in `c("SQC", "EQC", "QC")`, if both (default) then all QCs are considered as QC
#'    \item \code{"quantile"}: By quantile
#'    }
#'    Default: "Normalization" (if present, otherwise, "sum")
#' @param refSample String. Reference sample for `dataNormalize = "PQN2"`.
#' @param groupSample String. Used only if `dataNormalize = "groupPQN"`.
#' @param reference_method String. Method for computing reference from QC samples in `dataNormalize = "quantile"`. Options:
#'  \itemize{
#'    \item \code{"mean"}: Default
#'    \item \code{"median"}
#'  }
#' @param dataTransform String. Data transformation method applied after `dataNormalize`. Options:
#'  \itemize{
#'    \item \code{"none"}: No transformation
#'    \item \code{"log2"}: log base 2
#'    \item \code{"log10"}: log base 10
#'    \item \code{"sqrt"}: Square-root
#'    \item \code{"cbrt"}: Cube-root
#'    \item \code{"vsn"}: Variance Stabilizing Normalization
#'    \item \code{"glog"}: Variance stabilising generalised logarithm (glog) transformation (added August 24, 2025)
#'    }
#' @param dataScalePCA String. Data scaling for PCA analysis. This scaled data will be used in all analysis except PLS-type analysis. Options:
#'  \itemize{
#'    \item \code{"none"}: No data scaling
#'    \item \code{"mean"}: Scale by mean (average)
#'    \item \code{"meanSD"}: Auto-scaling. Scale by mean divided by standard deviation (SD)
#'    \item \code{"mean2SD"}: Pareto-scaling. Scale by mean divided by square-root of SD. Always use this for PLS-type analysis.
#'    }
#' @param dataScalePLS String. Data scaling for PLS analysis only. Same options as `dataScalePCA`.
#' @param filterMaxRSD Numeric. Threshold for Relative Standard Deviation (RSD) filtering. Set to `NULL` to skip this step.
#' @param filterMaxRSD_by String. Which QC samples to use for RSD filtering. Options:
#'  \itemize{
#'    \item \code{"SQC"}: Filter by sample QC
#'    \item \code{"EQC"}: Filter by extract QC (default)
#'    \item \code{"both"}: Filter by by both SQC and EQC (or QC altogether). Use this when there are no SQC and EQC in the "Group" row.
#'    }
#' @param filterMaxVarSD Numeric. Remove `nth percentile` of features with lowest variability. Set to `NULL` to skip this step.
#' @param verbose Logical. Whether to print detailed progress messages. Default TRUE.
#' @param auto_merge_replicates Logical. If `TRUE` (default), automatically merge technical replicates based on SubjectID and Replicate columns. Only applies when both columns contain meaningful values.
#'
#' @returns A list containing results from all preprocessing steps and a plot of six random features/metabolites before and after drift- and batch correction (if `TRUE`).
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{perform_DataQualityCheck}} for the data quality check and \code{\link{plot_before_after_correction}} to plot before and after drift- and batch correction
#'
#' @importFrom stats quantile
#'
#' @references Kirwan, J.A., Broadhurst, D.I., Davidson, R.L. et al. Characterising and correcting batch variation in an automated direct infusion mass spectrometry (DIMS) metabolomics workflow. Anal Bioanal Chem 405, 5147–5157 (2013). <https://doi.org/10.1007/s00216-013-6856-7> (pmp::QCRSC)
#' @references Parsons, H.M., Ludwig, C., Günther, U.L. et al. Improved classification accuracy in 1- and 2-dimensional NMR metabolomics data using the variance stabilising generalised logarithm transformation. BMC Bioinformatics 8, 234 (2007). <https://doi.org/10.1186/1471-2105-8-234> (pmp::glog_transformation)
#' @references Frank Dieterle, Alfred Ross, Götz Schlotterbeck, and Hans Senn. Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. Application in 1H NMR Metabonomics. Analytical Chemistry 2006 78 (13), 4281-4290. DOI: 10.1021/ac051632c <https://pubs.acs.org/doi/10.1021/ac051632c> (pmp::pqn_normalisation)
#' @references Jankevics A, Lloyd GR, Weber RJM (2025). pmp: Peak Matrix Processing and signal batch correction for metabolomics datasets. doi:10.18129/B9.bioc.pmp, R package version 1.20.0, <https://bioconductor.org/packages/pmp>. (pmp package)
#' @references Broadhurst, D.I. (2025). QC:MXP Repeat Injection based Quality Control, Batch Correction, Exploration & Data Cleaning (Version 2.1) Zendono. <https://doi.org/10.5281/zenodo.16824822>. Retrieved from <https://github.com/broadhurstdavid/QC-MXP>.
#' @references Variance stabilization applied to microarray data calibration and to the quantification of differential expression, Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie Poustka, Martin Vingron; Bioinformatics (2002) 18 Suppl.1 S96-S104. <https://doi.org/10.1093/bioinformatics/18.suppl_1.s96> (vsn::vsn2)
#' @references Huber W, von Heydebreck A, Sueltmann H, Poustka A, Vingron M. Parameter estimation for the calibration and variance stabilization of microarray data. Stat Appl Genet Mol Biol. 2003;2:Article3. doi: 10.2202/1544-6115.1008. Epub 2003 Apr 5. PMID: 16646781. <https://pubmed.ncbi.nlm.nih.gov/16646781/> (vsn::vsn2)
#' @references L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained Optimization, C. Zhu, R.H. Byrd, P. Lu and J. Nocedal, Technical Report, Northwestern University (1996). (vsn::vsn2)
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole. <https://doi.org/10.1201/9781351074988> (for scale)
#'
#' @export

perform_PreprocessingPeakData <- function(
    raw_data,
    outliers                    = NULL,
    filterMissing               = 20,
    filterMissing_by_group      = TRUE,
    filterMissing_includeQC     = FALSE,
    denMissing                  = 5,
    driftBatchCorrection        = TRUE,
    spline_smooth_param         = 0,
    spline_smooth_param_limit   = c(-1.5, 1.5),
    log_scale                   = TRUE,
    min_QC                      = 5,
    removeUncorrectedFeatures   = TRUE,
    dataNormalize               = "Normalization",
    refSample                   = NULL,
    groupSample                 = NULL,
    reference_method            = "mean",
    dataTransform               = "vsn",
    dataScalePCA                = "meanSD",
    dataScalePLS                = "mean2SD",
    filterMaxRSD                = 30,
    filterMaxRSD_by             = "EQC",
    filterMaxVarSD              = 10,
    verbose                     = TRUE,
    auto_merge_replicates       = TRUE
) {

  # Empty list for results
  listPreprocessed <- list()

  # Add timestamp start
  listPreprocessed$ProcessingTimestampStart <- Sys.time()

  # Helper function for conditional messaging
  msg <- function(...) if (verbose) message(...)

  # Helper function for safe numeric conversion
  safe_numeric <- function(x) {
    tryCatch(as.numeric(x),
             warning = function(w) {
               warning(paste("Warning in numeric conversion:", w$message))
               as.numeric(x)
             },
             error = function(e) {
               stop(paste("Error in numeric conversion:", e$message))
             })
  }

  # Helper function for memory-efficient operations
  gc_if_needed <- function() {
    if (gc.time()[1] %% 10 == 0) invisible(gc())
  }

  utils::flush.console()

  msg("Starting data preprocessing pipeline...")
  msg("Validating input data structure...")

  # Enhanced input validation
  if (!is.list(raw_data)) {
    stop("'raw_data' must be a list object from perform_DataQualityCheck function.")
  }

  if (!"FunctionOrigin" %in% names(raw_data)) {
    stop("'raw_data' is missing 'FunctionOrigin' field. Ensure data comes from perform_DataQualityCheck function.")
  }

  if (raw_data$FunctionOrigin != "perform_DataQualityCheck") {
    stop("The 'raw_data' must be from the 'perform_DataQualityCheck' function for quality check.")
  }

  if (!"raw_data" %in% names(raw_data)) {
    stop("'raw_data' is missing the actual data field 'raw_data'.")
  }

  msg("Input validation passed.")

  # Optimized parameter validation function
  validate_parameters <- function() {
    errors <- character()

    # Numeric parameter validation with optimized checks
    numeric_params <- list(
      filterMissing = list(value = filterMissing, range = c(1, 100)),
      denMissing = list(value = denMissing, range = c(1, 100)),
      spline_smooth_param = list(value = spline_smooth_param, range = c(0, 1)),
      min_QC = list(value = min_QC, range = c(1, Inf)),
      filterMaxRSD = list(value = filterMaxRSD, range = c(1, 100)),
      filterMaxVarSD = list(value = filterMaxVarSD, range = c(1, 100))
    )

    for (param_name in names(numeric_params)) {
      param_info <- numeric_params[[param_name]]
      value <- param_info$value

      # Skip validation for NULL values where appropriate
      if (is.null(value) && param_name %in% c("filterMissing", "filterMaxRSD", "filterMaxVarSD")) {
        next
      }

      if (!is.numeric(value) || length(value) != 1) {
        errors <- c(errors, paste(param_name, ": Must be a single numeric value."))
      } else if (value < param_info$range[1] || value > param_info$range[2]) {
        errors <- c(errors, paste(param_name, ": Out of valid range [",
                                  param_info$range[1], ",", param_info$range[2], "]."))
      }
    }

    # Logical parameter validation
    logical_params <- c("filterMissing_by_group", "filterMissing_includeQC",
                        "driftBatchCorrection", "log_scale", "removeUncorrectedFeatures",
                        "auto_merge_replicates")
    for (param_name in logical_params) {
      if (!is.logical(get(param_name))) {
        errors <- c(errors, paste(param_name, ": Must be logical (TRUE/FALSE)."))
      }
    }

    # Choice parameter validation
    choice_params <- list(
      dataNormalize = c("none", "Normalization", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile"),
      dataTransform = c("none", "log2", "log10", "sqrt", "cbrt", "vsn", "glog"),
      dataScalePCA = c("none", "mean", "meanSD", "mean2SD"),
      dataScalePLS = c("none", "mean", "meanSD", "mean2SD"),
      reference_method = c("mean", "median"),
      filterMaxRSD_by = c("SQC", "EQC", "both")
    )

    for (param_name in names(choice_params)) {
      value <- get(param_name)
      if (is.null(value) || length(value) != 1 || !value %in% choice_params[[param_name]]) {
        errors <- c(errors, paste(param_name, ": Must be one of:",
                                  paste(choice_params[[param_name]], collapse = ", ")))
      }
    }

    # Vector parameter validation
    if (!is.null(outliers) && !is.vector(outliers)) {
      errors <- c(errors, "outliers: Must be a vector of sample identifiers.")
    }

    if (!is.numeric(spline_smooth_param_limit) || length(spline_smooth_param_limit) != 2) {
      errors <- c(errors, "spline_smooth_param_limit: Must be a numeric vector of length 2.")
    }

    if (length(errors) > 0) {
      stop("Parameter validation failed:\n", paste(errors, collapse = "\n"))
    }

    msg("Parameter validation passed.")
  }

  validate_parameters()

  msg("Initializing preprocessing pipeline...")

  # Initialize results list with better structure
  listPreprocessed <- list(
    FunctionOrigin = "perform_PreprocessingPeakData",
    Parameters = list(
      outliers = outliers,
      filterMissing = filterMissing,
      filterMissing_by_group = filterMissing_by_group,
      filterMissing_includeQC = filterMissing_includeQC,
      denMissing = denMissing,
      driftBatchCorrection = driftBatchCorrection,
      removeUncorrectedFeatures = removeUncorrectedFeatures,
      dataNormalize = dataNormalize,
      dataTransform = dataTransform,
      dataScalePCA = dataScalePCA,
      dataScalePLS = dataScalePLS,
      filterMaxRSD = filterMaxRSD,
      filterMaxVarSD = filterMaxVarSD,
      auto_merge_replicates = auto_merge_replicates
    ),
    ProcessingSteps = character(),
    Dimensions = data.frame(
      Step = character(),
      Samples = numeric(),
      Features = numeric(),
      stringsAsFactors = FALSE
    )
  )

  # Optimized data transposition and initial processing
  msg("Transposing and preparing data...")

  tryCatch({
    # More efficient transposition
    data_transposed <- as.data.frame(t(raw_data$quality_checked_data))

    if (nrow(data_transposed) == 0) {
      stop("Transposed data has zero rows. Check input data structure.")
    }

    # Set column names more efficiently
    colnames(data_transposed) <- as.character(data_transposed[1, ])
    data_transposed <- data_transposed[-1, ]

    # Safer numeric conversion with error handling
    if ("Injection" %in% colnames(data_transposed)) {
      data_transposed$Injection <- safe_numeric(data_transposed$Injection)
      data_transposed <- data_transposed[order(data_transposed$Injection), ]
    } else {
      warning("'Injection' column not found. Data will not be sorted by injection sequence.")
    }

  }, error = function(e) {
    stop("Error in data transposition and preparation: ", e$message)
  })

  n_samples <- nrow(data_transposed)
  n_features_metabolites <- ncol(data_transposed) - 8

  msg(sprintf("Dataset contains %d samples and %d features/metabolites.",
              n_samples, n_features_metabolites))

  # Record initial dimensions
  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "Original", Samples = n_samples, Features = n_features_metabolites)
  )

  # Optimized outlier removal
  if (!is.null(outliers) && length(outliers) > 0) {
    msg("Removing specified outliers...")

    # More efficient outlier identification
    outliers_present <- outliers[outliers %in% data_transposed$Sample]
    outliers_missing <- setdiff(outliers, outliers_present)

    if (length(outliers_present) > 0) {
      data_transposed <- data_transposed[!data_transposed$Sample %in% outliers_present, ]
      msg(sprintf("Removed %d outliers: %s",
                  length(outliers_present),
                  paste(outliers_present, collapse = ", ")))
    }

    if (length(outliers_missing) > 0) {
      msg(sprintf("Outliers not found in data: %s",
                  paste(outliers_missing, collapse = ", ")))
    }

    listPreprocessed$outliers_removed <- outliers_present
    listPreprocessed$outliers_not_found <- outliers_missing

    listPreprocessed$Dimensions <- rbind(
      listPreprocessed$Dimensions,
      data.frame(Step = "After outlier removal",
                 Samples = nrow(data_transposed),
                 Features = ncol(data_transposed) - 8)
    )
  }

  # More efficient metadata extraction
  msg("Extracting metadata...")

  required_cols <- c("Sample", "SubjectID", "Replicate", "Group", "Group2", "Batch", "Injection", "Normalization", "Response")
  missing_cols <- setdiff(required_cols, colnames(data_transposed))

  if (length(missing_cols) > 0) {
    stop("Missing required columns in data: ", paste(missing_cols, collapse = ", "))
  }

  # Optimized metadata creation
  listPreprocessed$Metadata <- data.frame(
    Samples = data_transposed$Sample,
    SubjectID = safe_numeric(data_transposed$SubjectID),
    TechnicalReplicates = data_transposed$Replicate,
    Group = data_transposed$Group,
    Group2 = data_transposed$Group2,
    Group_ = gsub("SQC|EQC", "QC", data_transposed$Group),
    Batches = safe_numeric(data_transposed$Batch),
    InjectionSequence = safe_numeric(data_transposed$Injection),
    Normalization = safe_numeric(data_transposed$Normalization),
    Response = safe_numeric(data_transposed$Response),
    stringsAsFactors = FALSE
  )

  # More efficient QC indices
  indices_qc <- listPreprocessed$Metadata$Group %in% c("SQC", "EQC", "QC")
  indices_non_qc <- !indices_qc

  # Extract feature data more efficiently
  msg("Extracting feature data...")

  metadata_columns <- c("Sample", "SubjectID", "Replicate", "Group", "Group2", "Batch", "Injection", "Normalization", "Response")
  df <- data_transposed[, !colnames(data_transposed) %in% metadata_columns, drop = FALSE]

  # Store feature names
  listPreprocessed$All_Features_Metabolites <- colnames(df)

  # Optimized missing value handling
  msg("Processing missing values...")

  # More efficient zero-to-NA conversion
  df[df == 0] <- NA
  num_na <- sum(is.na(df))
  msg(sprintf("Found %d missing values in the data.", num_na))

  # Store number of missing data
  listPreprocessed$n_missing <- num_na

  # Set rownames efficiently
  rownames(df) <- data_transposed$Sample

  # Convert to numeric more safely
  df[] <- lapply(df, function(x) {
    tryCatch(as.numeric(x),
             warning = function(w) as.numeric(x),
             error = function(e) {
               warning("Could not convert column to numeric: ", e$message)
               x
             })
  })

  # Remove all-NA columns efficiently
  all_na_cols <- colSums(is.na(df)) == nrow(df)
  if (any(all_na_cols)) {
    n_removed <- sum(all_na_cols)
    df <- df[, !all_na_cols, drop = FALSE]
    msg(sprintf("Removed %d features with all missing values.", n_removed))
  }

  listPreprocessed$SamplesXFeatures <- df
  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After removing all-NA features",
               Samples = nrow(df), Features = ncol(df))
  )

  if (!is.null(filterMissing)) {

    # Optimized missing value filtering
    msg(sprintf("Applying missing value filter (threshold: %g%%)", filterMissing))

    filter_missing_features <- function(data, threshold, by_group, include_qc, metadata) {
      if (by_group) {
        if (include_qc) {
          groups <- unique(metadata$Group)
        } else {
          groups <- unique(metadata$Group[!metadata$Group %in% c("SQC", "EQC", "QC")])
        }

        # Calculate missing percentages per group more efficiently
        missing_by_group <- sapply(groups, function(g) {
          group_indices <- metadata$Group == g
          colMeans(is.na(data[group_indices, , drop = FALSE]))
        })

        # Remove features with high missingness in ALL groups
        remove_features <- apply(missing_by_group, 1, function(x) all(x >= threshold / 100))
      } else {
        if (include_qc) {
          missing_overall <- colMeans(is.na(data))
        } else {
          non_qc_indices <- !metadata$Group %in% c("SQC", "EQC", "QC")
          missing_overall <- colMeans(is.na(data[non_qc_indices, , drop = FALSE]))
        }
        remove_features <- missing_overall >= threshold / 100
      }

      return(!remove_features)  # Return which features to keep
    }

    # listPreprocessed$data_filteredMissing <- df
    # listPreprocessed$Dimensions <- rbind(
    #   listPreprocessed$Dimensions,
    #   data.frame(Step = "After missing value filter",
    #              Samples = nrow(df), Features = ncol(df))
    # )
  } else {
    msg("Skipping missing value filtering (filterMissing = NULL)")
    listPreprocessed$data_filteredMissing <- df
  }

  features_to_keep <- filter_missing_features(df, filterMissing, filterMissing_by_group,
                                              filterMissing_includeQC, listPreprocessed$Metadata)

  n_features_before <- ncol(df)
  df <- df[, features_to_keep, drop = FALSE]
  n_features_after <- ncol(df)

  msg(sprintf("Removed %d features due to missing value threshold.",
              n_features_before - n_features_after))

  listPreprocessed$data_filteredMissing <- df
  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After missing value filter",
               Samples = nrow(df), Features = ncol(df))
  )

  # Optimized missing value imputation
  msg("Performing missing value imputation...")

  # More efficient imputation
  df[] <- lapply(df, function(x) {
    if (any(is.na(x))) {
      positive_values <- x[x > 0 & !is.na(x)]
      if (length(positive_values) > 0) {
        min_val <- min(positive_values)
        x[is.na(x)] <- min_val / denMissing
      }
    }
    return(x)
  })

  listPreprocessed$data_no_NA <- df
  msg(sprintf("Imputed missing values using 1/%d of minimum positive value per feature.", denMissing))

  # Optimized drift and batch correction with uncorrected feature detection
  msg("Applying drift and batch correction...")

  if (driftBatchCorrection) {
    # Check for required packages
    if (!requireNamespace("pmp", quietly = TRUE)) {
      warning("Package 'pmp' not available. Skipping drift/batch correction.")
      driftBatchCorrection <- FALSE
    }
  }

  # Initialize variables for tracking uncorrected features
  uncorrected_features <- character(0)
  features_removed_uncorrected <- 0

  if (driftBatchCorrection) {
    tryCatch({
      # Store data before correction for comparison
      df_before_correction <- df

      df_corrected <- pmp::QCRSC(
        df = df,
        order = listPreprocessed$Metadata$InjectionSequence,
        batch = listPreprocessed$Metadata$Batches,
        classes = as.vector(listPreprocessed$Metadata$Group_),
        spar = spline_smooth_param,
        log = log_scale,
        minQC = min_QC,
        qc_label = "QC",
        spar_lim = spline_smooth_param_limit
      )

      df_corrected <- as.data.frame(t(df_corrected))
      msg("Applied QC-RSC drift and batch correction.")

      # Robust detection of uncorrected features
      tryCatch({
        msg("Identifying uncorrected features...")

        # Ensure both datasets have same dimensions and column names
        if (ncol(df_before_correction) != ncol(df_corrected) ||
            nrow(df_before_correction) != nrow(df_corrected)) {
          warning("Dimension mismatch between before/after correction data. Cannot identify uncorrected features.")
        } else if (!identical(colnames(df_before_correction), colnames(df_corrected))) {
          warning("Column name mismatch between before/after correction data. Cannot identify uncorrected features.")
        } else {
          # Compare each feature with tolerance for floating point precision
          tolerance <- .Machine$double.eps^0.5

          uncorrected_indices <- logical(ncol(df_before_correction))

          for (i in seq_len(ncol(df_before_correction))) {
            # Calculate maximum absolute difference for this feature
            max_diff <- max(abs(df_before_correction[, i] - df_corrected[, i]), na.rm = TRUE)

            # If difference is within tolerance, feature was not corrected
            uncorrected_indices[i] <- max_diff <= tolerance
          }

          uncorrected_features <- colnames(df_before_correction)[uncorrected_indices]

          msg(sprintf("Identified %d features that were not corrected due to insufficient QC samples.",
                      length(uncorrected_features)))

          # Store information about uncorrected features
          listPreprocessed$UncorrectedFeatures <- list(
            feature_names = uncorrected_features,
            count = length(uncorrected_features),
            reason = paste("Insufficient QC samples (< min_QC =", min_QC, ") for reliable correction")
          )

          # Remove uncorrected features if requested
          if (removeUncorrectedFeatures && length(uncorrected_features) > 0) {
            features_to_keep_corrected <- !colnames(df_corrected) %in% uncorrected_features
            df_corrected <- df_corrected[, features_to_keep_corrected, drop = FALSE]
            features_removed_uncorrected <- length(uncorrected_features)

            msg(sprintf("Removed %d uncorrected features from dataset.", features_removed_uncorrected))

            # Also remove from the before-correction data for consistency
            df_before_correction <- df_before_correction[, features_to_keep_corrected, drop = FALSE]
          } else if (!removeUncorrectedFeatures && length(uncorrected_features) > 0) {
            msg(sprintf("Keeping %d uncorrected features in dataset as requested.", length(uncorrected_features)))
          }
        }

      }, error = function(e) {
        warning("Error in identifying uncorrected features: ", e$message, ". Proceeding without removal.")
        listPreprocessed$UncorrectedFeatures <- list(
          feature_names = character(0),
          count = 0,
          reason = paste("Could not identify uncorrected features due to error:", e$message)
        )
      })

      # Update df with the corrected (and potentially filtered) data
      df <- df_corrected

      # Handle any new NAs introduced by correction
      new_na_count <- sum(is.na(df)) - sum(is.na(listPreprocessed$data_no_NA))
      if (new_na_count > 0) {
        msg(sprintf("Correction introduced %d new missing values. Re-imputing...", new_na_count))
        df[] <- lapply(df, function(x) {
          if (any(is.na(x))) {
            positive_values <- x[x > 0 & !is.na(x)]
            if (length(positive_values) > 0) {
              min_val <- min(positive_values)
              x[is.na(x)] <- min_val / denMissing
            }
          }
          return(x)
        })

        # Store new number of missing data
        listPreprocessed$new_n_missing <- new_na_count
      }

    }, error = function(e) {
      warning("Drift/batch correction failed: ", e$message, ". Proceeding without correction.")
      driftBatchCorrection <- FALSE
      # Initialize empty uncorrected features info for failed correction
      listPreprocessed$UncorrectedFeatures <- list(
        feature_names = character(0),
        count = 0,
        reason = paste("Drift/batch correction failed:", e$message)
      )
    })
  }

  # Add summary message about uncorrected features
  if (driftBatchCorrection) {
    if (length(uncorrected_features) > 0) {
      if (removeUncorrectedFeatures) {
        msg(sprintf("Removed %d uncorrected features from the dataset.", features_removed_uncorrected))
      } else {
        msg(sprintf("Kept %d uncorrected features in the dataset.", length(uncorrected_features)))
      }
    } else {
      msg("All features were successfully corrected by QCRSC.")
    }
  } else {
    msg("Skipping drift and batch correction.")
    # Initialize empty uncorrected features info when correction is skipped
    listPreprocessed$UncorrectedFeatures <- list(
      feature_names = character(0),
      count = 0,
      reason = "Drift/batch correction was not performed"
    )
  }

  listPreprocessed$data_driftBatchCorrected <- df

  # Add dimension entry for uncorrected feature removal (even if no features were removed)
  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After removing uncorrected features",
               Samples = nrow(df),
               Features = ncol(df))
  )

  if (features_removed_uncorrected > 0) {
    msg(sprintf("Dataset after removing uncorrected features: %d samples X %d features",
                nrow(df), ncol(df)))
  }

  gc_if_needed()

  # Optimized normalization function
  normalize_data <- function(data, method, metadata, indices_qc, indices_non_qc) {
    msg(sprintf("Applying %s normalization...", method))

    switch(method,
           "none" = {
             msg("No normalization applied.")
             return(data)
           },

           "Normalization" = {
             sg_values <- as.numeric(metadata$Normalization[indices_non_qc])
             if (all(is.na(sg_values))) {
               msg("Normalization values are all NA. Using sum normalization instead.")
               sample_sums <- rowSums(data, na.rm = TRUE)
               return(data / sample_sums)
             } else {
               data[indices_non_qc, ] <- data[indices_non_qc, ] / sg_values
               data[indices_qc, ] <- data[indices_qc, ] / rowSums(data[indices_qc, ], na.rm = TRUE)
               msg("Applied normalization using provided values.")
               return(data)
             }
           },

           "sum" = {
             sample_sums <- rowSums(data, na.rm = TRUE)
             msg("Applied sum normalization.")
             return(data / sample_sums)
           },

           "median" = {
             sample_medians <- apply(data, 1, median, na.rm = TRUE)
             msg("Applied median normalization.")
             return(data / sample_medians)
           },

           "PQN1" = {
             reference_spectrum <- apply(data, 2, median, na.rm = TRUE)
             quotients <- sweep(data, 2, reference_spectrum, "/")
             median_quotients <- apply(quotients, 1, median, na.rm = TRUE)
             msg("Applied PQN (global median approach).")
             return(data / median_quotients)
           },

           "PQN2" = {
             if (is.null(refSample) || !refSample %in% rownames(data)) {
               stop("Reference sample not specified or not found in data for PQN2 normalization.")
             }
             reference_spectrum <- as.numeric(data[refSample, ])
             quotients <- sweep(data, 2, reference_spectrum, "/")
             median_quotients <- apply(quotients, 1, median, na.rm = TRUE)
             msg(sprintf("Applied PQN using %s as reference.", refSample))
             return(data / median_quotients)
           },

           "groupPQN" = {
             if (is.null(groupSample)) {
               stop("Group sample not specified for groupPQN normalization.")
             }

             pooled_indices <- if (groupSample == "both") {
               metadata$Group %in% c("SQC", "EQC", "QC")
             } else {
               metadata$Group == groupSample
             }

             if (sum(pooled_indices) == 0) {
               stop("No samples found for specified group sample.")
             }

             reference_spectrum <- apply(data[pooled_indices, , drop = FALSE], 2, median, na.rm = TRUE)
             quotients <- sweep(data, 2, reference_spectrum, "/")
             median_quotients <- apply(quotients, 1, median, na.rm = TRUE)
             msg(sprintf("Applied group PQN using %s samples.", groupSample))
             return(data / median_quotients)
           },

           "quantile" = {
             if (!requireNamespace("pmp", quietly = TRUE)) {
               stop("Package 'pmp' required for quantile normalization.")
             }

             tryCatch({
               normalized <- pmp::pqn_normalisation(
                 t(data),
                 classes = metadata$Group_,
                 qc_label = "QC",
                 ref_mean = NULL,
                 qc_frac = 0,
                 sample_frac = 0,
                 ref_method = reference_method
               )
               msg("Applied quantile normalization.")
               return(as.data.frame(t(normalized)))
             }, error = function(e) {
               warning("Quantile normalization failed: ", e$message, ". Using sum normalization.")
               sample_sums <- rowSums(data, na.rm = TRUE)
               return(data / sample_sums)
             })
           },

           {
             warning("Unknown normalization method. No normalization applied.")
             return(data)
           }
    )
  }

  # Apply normalization
  df <- normalize_data(df, dataNormalize, listPreprocessed$Metadata, indices_qc, indices_non_qc)
  listPreprocessed$data_normalized <- df
  gc_if_needed()

  # Optimized transformation function
  transform_data <- function(data, method) {
    msg(sprintf("Applying %s transformation...", method))

    switch(method,
           "none" = {
             msg("No transformation applied.")
             return(data)
           },

           "log2" = {
             # Ensure positive values
             min_val <- min(data, na.rm = TRUE)
             if (min_val <= 0) {
               data <- data - min_val + 1
             }
             msg("Applied log2 transformation.")
             return(log2(data))
           },

           "log10" = {
             min_val <- min(data, na.rm = TRUE)
             if (min_val <= 0) {
               data <- data - min_val + 1
             }
             msg("Applied log10 transformation.")
             return(log10(data))
           },

           "sqrt" = {
             min_val <- min(data, na.rm = TRUE)
             if (min_val < 0) {
               data <- data - min_val
             }
             msg("Applied square root transformation.")
             return(sqrt(data))
           },

           "cbrt" = {
             min_val <- min(data, na.rm = TRUE)
             if (min_val < 0) {
               data <- data - min_val + 1
             }
             msg("Applied cube root transformation.")
             return(data^(1/3))
           },

           "vsn" = {
             if (!requireNamespace("vsn", quietly = TRUE)) {
               warning("Package 'vsn' not available. Using log10 transformation instead.")
               min_val <- min(data, na.rm = TRUE)
               if (min_val <= 0) {
                 data <- data - min_val + 1
               }
               return(log10(data))
             }

             tryCatch({
               data_matrix <- as.matrix(data)
               vsn_fit <- vsn::vsn2(data_matrix)
               result <- vsn::predict(vsn_fit, newdata = data_matrix) # Added "vsn::" to specify "predict" function origin
               msg("Applied VSN transformation.")
               return(as.data.frame(result))
             }, error = function(e) {
               warning("VSN transformation failed: ", e$message, ". Using log10 instead.")
               min_val <- min(data, na.rm = TRUE)
               if (min_val <= 0) {
                 data <- data - min_val + 1
               }
               return(log10(data))
             })
           },

           # Added August 24, 2025
           "glog" = {
             if (!requireNamespace("pmp", quietly = TRUE)) {
               warning("Package 'pmp' not available. Using log10 transformation instead.")
               min_val <- min(data, na.rm = TRUE)
               if (min_val <= 0) {
                 data <- data - min_val + 1
               }
               return(log10(data))
             }

             tryCatch({
               # pmp::glog_transformation expects Features x Samples, so transpose input
               data_matrix <- as.matrix(t(data))
               glog_result <- pmp::glog_transformation(data_matrix, classes = listPreprocessed$Metadata$Group_, qc_label = "QC")
               # Transpose back to Samples x Features
               result <- as.data.frame(t(glog_result))
               msg("Applied glog transformation.")
               return(result)
             }, error = function(e) {
               warning("glog transformation failed: ", e$message, ". Using log10 instead.")
               min_val <- min(data, na.rm = TRUE)
               if (min_val <= 0) {
                 data <- data - min_val + 1
               }
               return(log10(data))
             })
           },

           {
             warning("Unknown transformation method. No transformation applied.")
             return(data)
           }
    )
  }

  # Apply transformation
  df <- transform_data(df, dataTransform)
  listPreprocessed$data_transformed <- df
  gc_if_needed()

  # Optimized scaling function
  scale_data_optimized <- function(data, method, purpose) {
    msg(sprintf("Applying %s scaling for %s analysis...", method, purpose))

    switch(method,
           "none" = {
             msg(sprintf("No scaling applied for %s data.", purpose))
             return(data)
           },

           "mean" = {
             scaled_data <- scale(data, center = TRUE, scale = FALSE)
             msg(sprintf("Applied mean centering for %s analysis.", purpose))
             return(as.data.frame(scaled_data))
           },

           "meanSD" = {
             scaled_data <- scale(data, center = TRUE, scale = TRUE)
             msg(sprintf("Applied mean centering and unit variance scaling for %s analysis.", purpose))
             return(as.data.frame(scaled_data))
           },

           "mean2SD" = {
             # Pareto scaling
             sds <- apply(data, 2, sd, na.rm = TRUE)
             scaled_data <- scale(data, center = TRUE, scale = sqrt(sds))
             msg(sprintf("Applied Pareto scaling for %s analysis.", purpose))
             return(as.data.frame(scaled_data))
           },

           {
             warning("Unknown scaling method. No scaling applied.")
             return(data)
           }
    )
  }

  # Apply scaling for both PCA and PLS
  listPreprocessed$data_scaledPCA <- scale_data_optimized(df, dataScalePCA, "PCA")
  listPreprocessed$data_scaledPLS <- scale_data_optimized(df, dataScalePLS, "PLS")
  gc_if_needed()

  # Optimized RSD filtering function
  apply_rsd_filtering <- function(data, threshold, filter_by, metadata, data_for) {
    if (is.null(threshold)) {
      msg(sprintf("No RSD filtering applied for %s data.", data_for))
      return(data)
    }

    msg(sprintf("Applying RSD filtering (threshold: %g%%) for %s data...", threshold, data_for))

    if (!requireNamespace("pmp", quietly = TRUE)) {
      warning("Package 'pmp' not available. Skipping RSD filtering.")
      return(data)
    }

    tryCatch({
      # Determine QC label and classes based on filter_by parameter
      qc_label <- ifelse(filter_by == "both", "QC", filter_by)

      # Get the correct classes vector based on filter_by
      if (filter_by == "both") {
        classes <- metadata$Group_
      } else {
        classes <- metadata$Group
      }

      # Ensure we have the right number of samples by matching rownames
      data_sample_names <- rownames(data)
      metadata_sample_names <- metadata$Samples

      # Find matching samples between data and metadata
      matching_indices <- match(data_sample_names, metadata_sample_names)

      # Check for any non-matching samples
      if (any(is.na(matching_indices))) {
        warning("Some samples in data are not found in metadata. This may cause issues.")
        # Remove NAs and use only matching samples
        valid_indices <- !is.na(matching_indices)
        data <- data[valid_indices, , drop = FALSE]
        matching_indices <- matching_indices[valid_indices]
      }

      # Subset classes to match the data samples
      classes_matched <- classes[matching_indices]

      # Final dimension check
      if (length(classes_matched) != nrow(data)) {
        warning(sprintf("Dimension mismatch persists: classes (%d) vs data rows (%d). Skipping RSD filtering.",
                        length(classes_matched), nrow(data)))
        return(data)
      }

      # Check for NAs in classes
      if (any(is.na(classes_matched))) {
        warning("NA values found in sample classes. Skipping RSD filtering.")
        return(data)
      }

      # Check if we have the required QC samples
      if (!qc_label %in% classes_matched) {
        warning(sprintf("No %s samples found for RSD filtering. Skipping RSD filtering.", qc_label))
        return(data)
      }

      filtered_data <- pmp::filter_peaks_by_rsd(
        t(data),
        max_rsd = threshold,
        class = as.vector(classes_matched),
        qc_label = qc_label
      )

      result <- as.data.frame(t(filtered_data))
      n_removed <- ncol(data) - ncol(result)
      msg(sprintf("Removed %d features with RSD > %g%% for %s analysis.", n_removed, threshold, data_for))

      return(result)

    }, error = function(e) {
      warning("RSD filtering failed: ", e$message, ". Returning original data.")
      return(data)
    })
  }

  # Apply RSD filtering
  listPreprocessed$data_scaledPCA_rsdFiltered <- apply_rsd_filtering(
    listPreprocessed$data_scaledPCA, filterMaxRSD, filterMaxRSD_by,
    listPreprocessed$Metadata, "PCA"
  )

  listPreprocessed$data_scaledPLS_rsdFiltered <- apply_rsd_filtering(
    listPreprocessed$data_scaledPLS, filterMaxRSD, filterMaxRSD_by,
    listPreprocessed$Metadata, "PLS"
  )

  # Update dimensions for RSD filtered data
  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After RSD filtering (PCA)",
               Samples = nrow(listPreprocessed$data_scaledPCA_rsdFiltered),
               Features = ncol(listPreprocessed$data_scaledPCA_rsdFiltered))
  )

  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After RSD filtering (PLS)",
               Samples = nrow(listPreprocessed$data_scaledPLS_rsdFiltered),
               Features = ncol(listPreprocessed$data_scaledPLS_rsdFiltered))
  )

  gc_if_needed()

  # Optimized variance filtering function
  apply_variance_filtering <- function(data, threshold, data_for) {
    if (is.null(threshold)) {
      msg(sprintf("No variance filtering applied for %s data.", data_for))
      return(data)
    }

    msg(sprintf("Applying variance filtering (%gth percentile) for %s data...", threshold, data_for))

    tryCatch({
      # Calculate standard deviations more efficiently
      sds <- apply(data, 2, sd, na.rm = TRUE)

      # Handle invalid SDs
      valid_sds <- is.finite(sds) & !is.na(sds) & sds > 0

      if (sum(valid_sds) == 0) {
        warning(sprintf("No valid standard deviations for %s data. Returning original data.", data_for))
        return(data)
      }

      # Calculate threshold using only valid SDs
      threshold_value <- quantile(sds[valid_sds], threshold/100, na.rm = TRUE)

      # Keep features with SD > threshold AND valid SDs
      keep_features <- (sds > threshold_value) & valid_sds

      result <- data[, keep_features, drop = FALSE]
      n_removed <- ncol(data) - ncol(result)
      msg(sprintf("Removed %d features with low variance for %s analysis.", n_removed, data_for))

      return(result)

    }, error = function(e) {
      warning("Variance filtering failed: ", e$message, ". Returning original data.")
      return(data)
    })
  }

  # Apply variance filtering
  listPreprocessed$data_scaledPCA_varFiltered <- apply_variance_filtering(
    listPreprocessed$data_scaledPCA_rsdFiltered, filterMaxVarSD, "PCA"
  )

  listPreprocessed$data_scaledPLS_varFiltered <- apply_variance_filtering(
    listPreprocessed$data_scaledPLS_rsdFiltered, filterMaxVarSD, "PLS"
  )

  # Update dimensions for low variance-filtered data
  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After low variance filtering (PCA)",
               Samples = nrow(listPreprocessed$data_scaledPCA_varFiltered),
               Features = ncol(listPreprocessed$data_scaledPCA_varFiltered))
  )

  listPreprocessed$Dimensions <- rbind(
    listPreprocessed$Dimensions,
    data.frame(Step = "After low variance filtering (PLS)",
               Samples = nrow(listPreprocessed$data_scaledPLS_varFiltered),
               Features = ncol(listPreprocessed$data_scaledPLS_varFiltered))
  )

  # Helper function for merging technical replicates
  merge_technical_replicates <- function(data_pca, data_pls, metadata, verbose = TRUE) {
    msg <- function(...) if (verbose) message(...)

    tryCatch({
      msg("Checking for technical replicates to merge...")

      # Check if SubjectID and TechnicalReplicates columns have meaningful values
      has_subject_id <- !all(is.na(metadata$SubjectID)) &&
        length(unique(metadata$SubjectID[!is.na(metadata$SubjectID)])) > 1

      has_replicates <- !all(is.na(metadata$TechnicalReplicates)) &&
        length(unique(metadata$TechnicalReplicates[!is.na(metadata$TechnicalReplicates)])) > 1

      if (!has_subject_id || !has_replicates) {
        msg("No meaningful replicates found. Skipping replicate merging.")
        return(list(
          data_pca = data_pca,
          data_pls = data_pls,
          metadata = metadata,
          merged = FALSE,
          n_samples_before = nrow(metadata),
          n_samples_after = nrow(metadata)
        ))
      }

      msg("Technical replicates detected. Proceeding with merging...")

      # Create a comprehensive data frame for processing
      comprehensive_data_pca <- cbind(metadata, data_pca)
      comprehensive_data_pls <- cbind(metadata, data_pls)

      # Separate QC and non-QC samples
      qc_indices <- metadata$Group %in% c("SQC", "EQC", "QC")

      # Process QC samples (no averaging needed, just restructure)
      if (any(qc_indices)) {
        qc_data_pca <- comprehensive_data_pca[qc_indices, ]
        qc_data_pls <- comprehensive_data_pls[qc_indices, ]

        # Create consistent identifiers for QC samples (SIMPLIFIED)
        qc_data_pca$MergeID <- qc_data_pca$Samples
        qc_data_pls$MergeID <- qc_data_pls$Samples
      } else {
        qc_data_pca <- NULL
        qc_data_pls <- NULL
      }

      # Process non-QC samples (these need averaging)
      if (any(!qc_indices)) {
        non_qc_data_pca <- comprehensive_data_pca[!qc_indices, ]
        non_qc_data_pls <- comprehensive_data_pls[!qc_indices, ]

        # Create merge identifier for non-QC samples (SIMPLIFIED)
        non_qc_data_pca$MergeID <- paste0(non_qc_data_pca$TechnicalReplicates)
        non_qc_data_pls$MergeID <- paste0(non_qc_data_pls$TechnicalReplicates)

        # Get feature column indices for PCA data
        feature_cols_pca <- which(!colnames(non_qc_data_pca) %in%
                                    c("Samples", "SubjectID", "TechnicalReplicates", "Group", "Group2",
                                      "Group_", "Batches", "InjectionSequence", "Normalization",
                                      "Response", "MergeID"))

        # Get feature column indices for PLS data (SEPARATE!)
        feature_cols_pls <- which(!colnames(non_qc_data_pls) %in%
                                    c("Samples", "SubjectID", "TechnicalReplicates", "Group", "Group2",
                                      "Group_", "Batches", "InjectionSequence", "Normalization",
                                      "Response", "MergeID"))

        # Average technical replicates for PCA data
        merged_non_qc_pca <- non_qc_data_pca %>%
          group_by(MergeID) %>%
          summarise(
            # Keep first occurrence of metadata
            TechnicalReplicates = first(TechnicalReplicates),
            SubjectID = first(SubjectID),
            Group = first(Group),
            Group2 = first(Group2),
            Batches = first(Batches),
            # Add to check if this works
            InjectionSequence = first(InjectionSequence),
            Normalization = first(Normalization),
            Response = first(Response),
            # Average the feature data
            across(all_of(feature_cols_pca), ~mean(., na.rm = TRUE)),
            .groups = 'drop'
          ) %>%
          arrange(SubjectID)

        # Average technical replicates for PLS data
        merged_non_qc_pls <- non_qc_data_pls %>%
          group_by(MergeID) %>%
          summarise(
            # Keep first occurrence of metadata
            TechnicalReplicates = first(TechnicalReplicates),
            SubjectID = first(SubjectID),
            Group = first(Group),
            Group2 = first(Group2),
            Batches = first(Batches),
            # Add to check if this works
            InjectionSequence = first(InjectionSequence),
            Normalization = first(Normalization),
            Response = first(Response),
            # Average the feature data
            across(all_of(feature_cols_pls), ~mean(., na.rm = TRUE)),
            .groups = 'drop'
          ) %>%
          arrange(SubjectID)
      } else {
        merged_non_qc_pca <- NULL
        merged_non_qc_pls <- NULL
      }

      # Combine QC and merged non-QC data
      if (!is.null(qc_data_pca) && !is.null(merged_non_qc_pca)) {
        # Ensure QC data has the same structure as merged non-QC data
        qc_subset_pca <- qc_data_pca[, colnames(merged_non_qc_pca)]
        qc_subset_pls <- qc_data_pls[, colnames(merged_non_qc_pls)]

        final_data_pca <- rbind(merged_non_qc_pca, qc_subset_pca)
        final_data_pls <- rbind(merged_non_qc_pls, qc_subset_pls)
      } else if (!is.null(merged_non_qc_pca)) {
        final_data_pca <- merged_non_qc_pca
        final_data_pls <- merged_non_qc_pls
      } else if (!is.null(qc_data_pca)) {
        feature_cols_qc_pca <- which(!colnames(qc_data_pca) %in%
                                       c("Samples", "SubjectID", "TechnicalReplicates", "Group", "Group2",
                                         "Group_", "Batches", "InjectionSequence", "Normalization", "Response", "MergeID"))
        feature_cols_qc_pls <- which(!colnames(qc_data_pls) %in%
                                       c("Samples", "SubjectID", "TechnicalReplicates", "Group", "Group2",
                                         "Group_", "Batches", "InjectionSequence", "Normalization", "Response", "MergeID"))

        final_data_pca <- qc_data_pca[, c("MergeID", "TechnicalReplicates", "SubjectID", "Group", "Group2",
                                          "Batches", "InjectionSequence", "Normalization", "Response",
                                          colnames(qc_data_pca)[feature_cols_qc_pca])]
        final_data_pls <- qc_data_pls[, c("MergeID", "TechnicalReplicates", "SubjectID", "Group", "Group2",
                                          "Batches", "InjectionSequence", "Normalization", "Response",
                                          colnames(qc_data_pls)[feature_cols_qc_pls])]
      } else {
        stop("No valid data found for merging.")
      }

      # Arrange by SubjectID
      final_data_pca <- final_data_pca %>% arrange(SubjectID)
      final_data_pls <- final_data_pls %>% arrange(SubjectID)

      # Prepare final outputs
      # Extract feature data (remove metadata columns)
      feature_cols_final_pca <- which(!colnames(final_data_pca) %in%
                                        c("MergeID", "TechnicalReplicates", "SubjectID",
                                          "Group", "Group2", "Batches", "InjectionSequence", "Normalization", "Response"))
      feature_cols_final_pls <- which(!colnames(final_data_pls) %in%
                                        c("MergeID", "TechnicalReplicates", "SubjectID",
                                          "Group", "Group2", "Batches", "InjectionSequence", "Normalization", "Response"))

      merged_features_pca <- final_data_pca[, feature_cols_final_pca, drop = FALSE]
      merged_features_pls <- final_data_pls[, feature_cols_final_pls, drop = FALSE]

      # Convert to data.frame and set row names to MergeID (FIXED)
      merged_features_pca <- as.data.frame(merged_features_pca)
      merged_features_pls <- as.data.frame(merged_features_pls)
      rownames(merged_features_pca) <- final_data_pca$MergeID
      rownames(merged_features_pls) <- final_data_pls$MergeID

      # Create new metadata
      merged_metadata <- data.frame(
        Samples = final_data_pca$MergeID,  # Use MergeID as sample identifier
        SubjectID = final_data_pca$SubjectID,
        TechnicalReplicates = final_data_pca$TechnicalReplicates,
        Group = final_data_pca$Group,
        Group2 = final_data_pca$Group2,
        Group_ = gsub("SQC|EQC", "QC", final_data_pca$Group),
        Batches = final_data_pca$Batches,
        # Removed if code above works
        # InjectionSequence = seq_len(nrow(final_data_pca)),  # Create new sequence
        # Normalization = rep(1, nrow(final_data_pca)),       # Default normalization
        # Response = rep(NA, nrow(final_data_pca)),           # Default response
        # Another code
        # InjectionSequence = first(metadata$InjectionSequence),
        # Normalization = first(metadata$Normalization),
        # Response = first(metadata$Response),
        # Yet another code to see if this works
        InjectionSequence = final_data_pca$InjectionSequence,
        Normalization = final_data_pca$Normalization,
        Response = final_data_pca$Response,
        stringsAsFactors = FALSE
      )

      n_samples_before <- nrow(metadata)
      n_samples_after <- nrow(merged_metadata)
      n_qc <- sum(merged_metadata$Group_ == "QC")

      msg(sprintf("Successfully merged technical replicates: %d samples -> %d samples (including %d QC samples)",
                  n_samples_before, n_samples_after, n_qc))

      return(list(
        data_pca = merged_features_pca,
        data_pls = merged_features_pls,
        metadata = merged_metadata,
        merged = TRUE,
        n_samples_before = n_samples_before,
        n_samples_after = n_samples_after
      ))

    }, error = function(e) {
      warning("Error in merging technical replicates: ", e$message, ". Returning original data.")
      return(list(
        data_pca = data_pca,
        data_pls = data_pls,
        metadata = metadata,
        merged = FALSE,
        n_samples_before = nrow(metadata),
        n_samples_after = nrow(metadata),
        error = e$message
      ))
    })
  }

  # Technical replicate merging
  if (auto_merge_replicates) {
    msg("Attempting to merge technical replicates...")

    merge_result <- merge_technical_replicates(
      data_pca = listPreprocessed$data_scaledPCA_varFiltered,
      data_pls = listPreprocessed$data_scaledPLS_varFiltered,
      metadata = listPreprocessed$Metadata,
      verbose = verbose
    )

    if (merge_result$merged) {
      # Store merged data
      listPreprocessed$data_scaledPCA_merged <- merge_result$data_pca
      listPreprocessed$data_scaledPLS_merged <- merge_result$data_pls
      listPreprocessed$Metadata_merged <- merge_result$metadata

      # Update dimensions for merged data
      listPreprocessed$Dimensions <- rbind(
        listPreprocessed$Dimensions,
        data.frame(Step = "After replicate merging (PCA)",
                   Samples = nrow(merge_result$data_pca),
                   Features = ncol(merge_result$data_pca))
      )

      listPreprocessed$Dimensions <- rbind(
        listPreprocessed$Dimensions,
        data.frame(Step = "After replicate merging (PLS)",
                   Samples = nrow(merge_result$data_pls),
                   Features = ncol(merge_result$data_pls))
      )

      # Store merging information
      listPreprocessed$ReplicateMerging <- list(
        merged = TRUE,
        n_samples_before = merge_result$n_samples_before,
        n_samples_after = merge_result$n_samples_after,
        reduction_count = merge_result$n_samples_before - merge_result$n_samples_after
      )

      # msg(sprintf("Replicate merging completed: %d -> %d samples",
      #             merge_result$n_samples_before, merge_result$n_samples_after))

    } else {
      listPreprocessed$ReplicateMerging <- list(
        merged = FALSE,
        reason = "No meaningful replicates found or merging failed"
      )

      if (!is.null(merge_result$error)) {
        listPreprocessed$ReplicateMerging$error <- merge_result$error
      }
    }

    # # Update final dimensions
    # listPreprocessed$Dimensions <- rbind(
    #   listPreprocessed$Dimensions,
    #   data.frame(Step = "Final (PCA ready)",
    #              Samples = nrow(listPreprocessed$data_scaledPCA_varFiltered),
    #              Features = ncol(listPreprocessed$data_scaledPCA_varFiltered))
    # )
    # listPreprocessed$Dimensions <- rbind(
    #   listPreprocessed$Dimensions,
    #   data.frame(Step = "Final (PLS ready)",
    #              Samples = nrow(listPreprocessed$data_scaledPLS_varFiltered),
    #              Features = ncol(listPreprocessed$data_scaledPLS_varFiltered))
    # )

  } else {
    msg("Automatic replicate merging disabled.")
    listPreprocessed$ReplicateMerging <- list(
      merged = FALSE,
      reason = "Disabled by user (auto_merge_replicates = FALSE)"
    )

    # # Update final dimensions
    # listPreprocessed$Dimensions <- rbind(
    #   listPreprocessed$Dimensions,
    #   data.frame(Step = "Final (PCA ready)",
    #              Samples = nrow(listPreprocessed$data_scaledPCA_varFiltered),
    #              Features = ncol(listPreprocessed$data_scaledPCA_varFiltered))
    # )
    # listPreprocessed$Dimensions <- rbind(
    #   listPreprocessed$Dimensions,
    #   data.frame(Step = "Final (PLS ready)",
    #              Samples = nrow(listPreprocessed$data_scaledPLS_varFiltered),
    #              Features = ncol(listPreprocessed$data_scaledPLS_varFiltered))
    # )
  }

  gc_if_needed()

  # # Update final dimensions
  # listPreprocessed$Dimensions <- rbind(
  #   listPreprocessed$Dimensions,
  #   data.frame(Step = "Final (PCA ready)",
  #              Samples = nrow(listPreprocessed$data_scaledPCA_varFiltered),
  #              Features = ncol(listPreprocessed$data_scaledPCA_varFiltered))
  # )
  #
  # listPreprocessed$Dimensions <- rbind(
  #   listPreprocessed$Dimensions,
  #   data.frame(Step = "Final (PLS ready)",
  #              Samples = nrow(listPreprocessed$data_scaledPLS_varFiltered),
  #              Features = ncol(listPreprocessed$data_scaledPLS_varFiltered))
  # )

  # Add processing summary (updated to include uncorrected features info)
  listPreprocessed$ProcessingSummary <- list(
    total_samples_processed = nrow(listPreprocessed$Metadata),
    total_features_original = length(listPreprocessed$All_Features_Metabolites),
    total_features_final_PCA = ncol(listPreprocessed$data_scaledPCA_varFiltered),
    total_features_final_PLS = ncol(listPreprocessed$data_scaledPLS_varFiltered),
    outliers_removed = length(listPreprocessed$outliers_removed %||% character(0)),
    missing_values_imputed = sum(is.na(listPreprocessed$SamplesXFeatures)),
    drift_batch_correction_applied = driftBatchCorrection,
    uncorrected_features_count = length(uncorrected_features),
    uncorrected_features_removed = features_removed_uncorrected,
    normalization_method = dataNormalize,
    transformation_method = dataTransform,
    scaling_method_PCA = dataScalePCA,
    scaling_method_PLS = dataScalePLS,
    rsd_filtering_applied = !is.null(filterMaxRSD),
    variance_filtering_applied = !is.null(filterMaxVarSD),
    replicates_merged = listPreprocessed$ReplicateMerging$merged
  )

  # Add data quality metrics (updated)
  listPreprocessed$QualityMetrics <- list(
    data_completeness = 1 - sum(is.na(listPreprocessed$data_transformed)) /
      (nrow(listPreprocessed$data_transformed) * ncol(listPreprocessed$data_transformed)),
    feature_retention_rate_PCA = ncol(listPreprocessed$data_scaledPCA_varFiltered) /
      length(listPreprocessed$All_Features_Metabolites),
    feature_retention_rate_PLS = ncol(listPreprocessed$data_scaledPLS_varFiltered) /
      length(listPreprocessed$All_Features_Metabolites),
    sample_retention_rate = nrow(listPreprocessed$Metadata) /
      (nrow(listPreprocessed$Metadata) + length(listPreprocessed$outliers_removed %||% character(0))),
    uncorrected_features_percentage = ifelse(driftBatchCorrection,
                                             length(uncorrected_features) /
                                               length(listPreprocessed$All_Features_Metabolites) * 100,
                                             0)
  )

  # Final cleanup
  gc_if_needed()

  if (auto_merge_replicates) {
    # Summary message about replicate merging
    if (listPreprocessed$ReplicateMerging$merged) {
      # msg(sprintf("Technical replicates merged: %d -> %d samples",
      #             listPreprocessed$ReplicateMerging$n_samples_before,
      #             listPreprocessed$ReplicateMerging$n_samples_after))
      msg(sprintf("Final merged dataset dimensions: %d samples X %d features (PCA-ready)",
                  nrow(listPreprocessed$data_scaledPCA_merged),
                  ncol(listPreprocessed$data_scaledPCA_merged)))
      msg(sprintf("Final merged dataset dimensions: %d samples X %d features (PLS-ready)",
                  nrow(listPreprocessed$data_scaledPLS_merged),
                  ncol(listPreprocessed$data_scaledPLS_merged)))
    }
  } else {
    msg(sprintf("Final dataset dimensions: %d samples X %d features (PCA-ready)",
                nrow(listPreprocessed$data_scaledPCA_varFiltered),
                ncol(listPreprocessed$data_scaledPCA_varFiltered)))
    msg(sprintf("Final dataset dimensions: %d samples X %d features (PLS-ready)",
                nrow(listPreprocessed$data_scaledPLS_varFiltered),
                ncol(listPreprocessed$data_scaledPLS_varFiltered)))
  }

  msg("Data preprocessing pipeline completed successfully!")

  # Add timestamp end
  listPreprocessed$ProcessingTimestampEnd <- Sys.time()

  # Add warnings if any steps were skipped (updated)
  warnings_list <- character(0)
  if (!driftBatchCorrection) {
    warnings_list <- c(warnings_list, "Drift/batch correction was skipped")
  }
  if (length(uncorrected_features) > 0 && !removeUncorrectedFeatures) {
    warnings_list <- c(warnings_list,
                       paste("Dataset contains", length(uncorrected_features),
                             "uncorrected features that were retained"))
  }

  if (length(warnings_list) > 0) {
    listPreprocessed$Warnings <- warnings_list
  }

  return(listPreprocessed)
}

# Added August 24, 2025

# Function to create before and after batch correction plots:

#' Create Before and After Batch Correction Plots
#'
#' @description
#' Creates visualization plots comparing data before and after drift/batch correction.
#' This function generates plots showing the correction effects for selected features.
#'
#' @param preprocessed_data List. Results from the `perform_PreprocessingPeakData` function.
#' @param n_features Numeric. Number of random features to plot. Default is 6.
#' @param feature_names Vector. Specific feature names to plot. If provided, overrides n_features.
#'
#' @returns A combined plot showing before and after correction for selected features.
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap scale_color_manual labs theme_minimal theme element_text
#' @importFrom tidyr pivot_longer
#' @importFrom gridExtra grid.arrange
#'
#' @export

plot_before_after_correction <- function(preprocessed_data, n_features = 6, feature_names = NULL) {

  # Input validation
  if (!is.list(preprocessed_data)) {
    stop("'preprocessed_data' must be a list object from perform_PreprocessingPeakData function.")
  }

  if(!preprocessed_data$Parameters$driftBatchCorrection) {
    stop("The 'driftBatchCorrection' in the 'preprocessed_data' provided was set to FALSE, hence, no data from drift and batch correction was generated.")
  }

  if (!"FunctionOrigin" %in% names(preprocessed_data)) {
    stop("'preprocessed_data' is missing 'FunctionOrigin' field.")
  }

  if (preprocessed_data$FunctionOrigin != "perform_PreprocessingPeakData") {
    stop("The 'preprocessed_data' must be from the 'perform_PreprocessingPeakData' function.")
  }

  # Check if required data fields exist
  required_fields <- c("data_no_NA", "data_driftBatchCorrected", "Metadata")
  missing_fields <- setdiff(required_fields, names(preprocessed_data))

  if (length(missing_fields) > 0) {
    stop("Missing required fields in preprocessed_data: ",
         paste(missing_fields, collapse = ", "))
  }

  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Required packages for plotting not available. Please install: ggplot2, tidyr, gridExtra")
  }

  # Extract data
  df_before <- preprocessed_data$data_no_NA
  df_after <- preprocessed_data$data_driftBatchCorrected
  metadata <- preprocessed_data$Metadata

  # Check if correction was actually applied
  if (identical(df_before, df_after)) {
    warning("Data before and after correction are identical. No correction was applied.")
  }

  # Feature selection
  if (!is.null(feature_names)) {
    # Use specified features
    available_features <- intersect(feature_names, colnames(df_before))
    if (length(available_features) == 0) {
      stop("None of the specified feature names found in the data.")
    }
    if (length(available_features) < length(feature_names)) {
      warning(paste("Some specified features not found:",
                    paste(setdiff(feature_names, available_features), collapse = ", ")))
    }
    selected_features <- match(available_features, colnames(df_before))
    n_features <- length(selected_features)
  } else {
    # Random feature selection
    n_features <- min(n_features, ncol(df_before), ncol(df_after))
    if (n_features < 1) {
      stop("No features available for plotting.")
    }
    selected_features <- sample(1:ncol(df_before), n_features, replace = FALSE)
  }

  # Create individual plots
  plot_list <- list()

  for (i in seq_along(selected_features)) {
    feature_idx <- selected_features[i]
    feature_name <- colnames(df_before)[feature_idx]

    # Prepare plot data with better error handling
    plot_data <- data.frame(
      sample_id = 1:nrow(df_before),
      before = log10(pmax(df_before[, feature_idx], 1e-10)),  # Prevent log of zero
      after = log10(pmax(df_after[, feature_idx], 1e-10)),
      class = metadata$Group_,
      batch = as.factor(metadata$Batches),
      stringsAsFactors = FALSE
    )

    # Transform to long format
    plot_data_long <- tidyr::pivot_longer(
      plot_data,
      cols = c("before", "after"),
      names_to = "correction",
      values_to = "intensity"
    )

    plot_data_long$correction <- factor(plot_data_long$correction,
                                        levels = c("before", "after"))

    # Create plot
    p <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = sample_id, y = intensity)) +
      ggplot2::geom_point(ggplot2::aes(color = class, shape = batch),
                          size = 2, alpha = 0.7) +
      ggplot2::geom_smooth(data = subset(plot_data_long, class == "QC"),
                           formula = y ~ x,
                           method = "loess", se = TRUE, alpha = 0.3,
                           color = "darkred", linetype = "dashed") +
      ggplot2::facet_wrap(~ correction, scales = "free_y",
                          labeller = ggplot2::labeller(
                            correction = c(before = "Before Correction",
                                           after = "After Correction"))) +
      ggplot2::scale_color_manual(values = c("QC" = "red", "Sample" = "blue",
                                             "Blank" = "green")) +
      ggplot2::labs(title = paste("Feature:", feature_name),
                    x = "Sample Index",
                    y = "log10(Intensity)",
                    color = "Class", shape = "Batch") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 11),
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 9),
        strip.text = ggplot2::element_text(size = 10, face = "bold")
      )

    plot_list[[i]] <- p
  }

  # Create combined plot
  n_batches <- length(unique(metadata$Batches))
  main_title <- ifelse(n_batches == 1,
                       "Random Features Before and After Drift-Correction",
                       "Random Features Before and After Drift- and Batch-Correction")

  combined_plot <- do.call(gridExtra::grid.arrange,
                           c(plot_list, list(ncol = 2, top = main_title)))

  message("Before/after correction plots created successfully.")
  return(combined_plot)
}
