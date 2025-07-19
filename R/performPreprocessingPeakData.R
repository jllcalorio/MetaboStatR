#' Perform Data Preprocessing on a Quality-Checked Data
#'
#' @description
#' This function performs data preprocessing techniques to prepare the data for downstream analysis.
#' The data preprocessing includes data filtration, missing value imputation, drift and batch correction,
#' transformation, scaling, data normalization, and removal of known, identified, or perceived outliers.
#'
#' @param raw_data List. A quality-checked data from the `performDataQualityCheck` function.
#' @param filterMissing Numeric. Minimum %missing in all groups to remove feature.
#' @param filterMissing_by_group Boolean. Determines if the `filterMissing` should detect group missingness before removal of features/metabolites. Example, if there are 2 groups (Group1, and Group2), and say `filterMissing`% is present in Group1 and the same percentage also applies to Group2 in Feature1, then Feature1 is removed. Set to `FALSE` to ignore grouping.
#' @param filterMissing_includeQC Boolean. If `TRUE` (default), when implementing `filterMissing`, includes QC samples. Example: Feature1 is 20% missing in QC, Group1, then Feature1 will be removed. Set to `FALSE` if only the biological samples are required.
#' @param denMissing Numeric. A value to be used in the denominator `1/denMissing` which will be used to replace missing values. The missing values in each feature will be replaced with the `1/denMissing` of the smallest positive value in a feature.
#' @param driftBatchCorrection Boolean. If `TRUE` (default), perform Quality Control-Robust Spline Correction (QC-RSC) algorithm for signal drift and batch effect correction within/across a multi-batch direct infusion mass spectrometry (DIMS) and liquid chromatography mass spectrometry (LCMS) datasets. Read more on `?pmp::QCRSC`.
#' @param spline_smooth_param Numeric. Only used in `driftBatchCorrection`. Spline smoothing parameter. Should be in the range 0 to 1. If set to 0 it will be estimated using leave-one-out cross-validation (to avoid overfitting).
#' @param spline_smooth_param_limit Vector. Only used in `driftBatchCorrection`. A vector of format `c(num1, num2)` signifying the minimum and maximum values of `spline_smooth_param` when searching for an optimum.
#' @param log_scale Boolean. Only used in `driftBatchCorrection`. If `TRUE` (default), performs the signal correction fit on the log scaled data.
#' @param min_QC Numeric. Only used in `driftBatchCorrection`. The minimum number of measured quality control (QC) samples required for signal correction within feature per batch. For features where signal was measured in less QC samples than threshold signal correction won't be applied.
#' @param dataNormalize String. Perform data normalization. Options are:
#'   \itemize{
#'     \item "none": No normalization.
#'     \item "DilutionMarker": Normalization using values provided in the "DilutionMarker" row. Otherwise, normalized by 'sum'.
#'     \item "sum": Normalization by sum.
#'     \item "median": Normalization by median.
#'     \item "PQN1": Probabilistic Quotient Normalization (PQN) according to global median approach.
#'     \item "PQN2": Probabilistic Quotient Normalization (PQN) using a reference sample indicated in `refSample` as the reference sample in the normalization.
#'     \item "groupPQN": Group Probabilistic Quotient Normalization using pooled group indicated in `groupSample`.
#'     \item "quantile": Normalization by quantile. Read more on `?pmp::pqn_normalisation`.
#'     }
#'     Defaults to 'DilutionMarker' if present, otherwise, normalization by 'sum' will be performed.
#' @param refSample String. The reference sample in the case of `dataNormalize = "PQN2"`. Must be in the samples, and is not part of "outliers" vector c("SQC", "EQC", "both").
#' @param groupSample String. Used only if `dataNormalize = "groupPQN"`. Ignored if not otherwise. Default to "EQC" if "groupPQN". Other choice is "SQC".
#' @param reference_method String. Only used if `dataNormalize = "quantile"`. The method used to compute the reference from the QC samples. Choices are below. Read more on `?pmp::pqn_normalisation`.
#'   \itemize{
#'     \item "mean"
#'     \item "median"
#'     }
#'     Defaults to "mean".
#' @param dataTransform String. A transformation method to transform the data after `dataNormalize`.
#'   \itemize{
#'     \item "none": No data transformation is done.
#'     \item "log2": Perform log2 transformation.
#'     \item "log10": Perform log10 transformation.
#'     \item "sqrt": Perform square root transformation.
#'     \item "cbrt": Perform cube root transformation.
#'     \item "vsn": Perform Variance Stabilizing Normalization. Read more on `?vsn::justvsn`. In its description, "The method uses a robust variant of the maximum-likelihood estimator for an additive-multiplicative error model and affine calibration. The model incorporates data calibration step (a.k.a. normalization), a model for the dependence of the variance on the mean intensity and a variance stabilizing data transformation. Differences between transformed intensities are analogous to "normalized log-ratios". However, in contrast to the latter, their variance is independent of the mean, and they are usually more sensitive and specific in detecting differential transcription."
#'     }
#'     Defaults to "log10".
#' @param dataScalePCA String. A data scaling done to the data after `dataTransform`. This data will be used later for Principal Component Analysis (PCA) "only".
#'   \itemize{
#'     \item "none": No data scaling is performed.
#'     \item "mean": Mean-centered only
#'     \item "meanSD": Mean-centered and divided by SD of each feature
#'     \item "meanSD2": Mean-centered and divided by the square root of SD of each feature. Also called pareto-scaling.
#'     }
#'     Defaults to "meanSD".
#' @param dataScaleOPLSDA String. A data scaling done to the data after `dataTransform`. This data will be used later for analyses "other than PCA", e.g., Orthogonal Partial Least Squares-Discriminant Analysis (OPLS-DA), among others.
#'   \itemize{
#'     \item "none": No data scaling is performed.
#'     \item "mean": Mean-centered only.
#'     \item "meanSD": Mean-centered and divided by SD of each feature.
#'     \item "meanSD2": Mean-centered and divided by the square root of SD of each feature. Also called pareto-scaling.
#'     }
#'     Defaults to "meanSD2".
#' @param filterMaxRSD Numeric. The threshold for the Relative Standard Deviation (RSD). Suggestions below. RSD is used to assess and filter out unreliable features, ensuring that only high-quality, reproducible data are used for statistical and biological interpretation. RSD is calculated for each feature across Quality Control (QC) samples. Features where QCs have RSD greater than the desired `filterMaxRSD` threshold are removed.
#'   \itemize{
#'     \item 20: for LC-MS analyzed data. Remove feature if RSD of QC >= 20%.
#'     \item 30: for GC-MS analyzed data. Remove feature if RSD of QC >= 30%.
#'     }
#'     Defaults to 30.
#' @param filterMaxRSD_by String. Choose below where to apply the RSD filtering in `filterMaxRSD`.
#'   \itemize{
#'     \item "SQC": Apply on Sample QCs only.
#'     \item "EQC": Apply on Extract QCs only.
#'     \item "both": Apply on both, which means all SQCs and EQCs are taken as one data. This is also the option when SQC and EQC are not present.
#'     }
#'     Defaults to "EQC".
#' @param filterMaxVarSD Numeric. Remove `nth` percentile (e.g., 10 for 10%) of features with the lowest variability. This removes features where variation between the groups is very low. Set to `NULL` to skip this filtering.
#' @param outliers A vector of biological samples and/or QC that are considered as outliers. Example format is "c('Sample1', 'Sample2', 'QC1', 'QC2', ...)". Defaults to `NULL`.
#'
#' @returns A list of results from all of the tests done (e.g., batch-corrected data, normalized data, transformed data, scaled data, etc.).
#' @export
#'
#' @examples
#' \dontrun{
#' # Using the default parameters
#' performPreprocessingPeakData(
#'  raw_data = a_csv_file
#' )
#' }

performPreprocessingPeakData <- function(
    raw_data,
    filterMissing             = 20,
    filterMissing_by_group    = TRUE,
    filterMissing_includeQC   = TRUE,
    denMissing                = 5,
    driftBatchCorrection      = TRUE,
    spline_smooth_param       = 0,
    spline_smooth_param_limit = c(-1.5, 1.5),
    log_scale                 = TRUE,
    min_QC                    = 5,
    dataTransform             = c("none", "log2", "log10", "sqrt", "cbrt", "vsn")[3],
    dataScalePCA              = c("none", "mean", "meanSD", "meanSD2")[3],
    dataScaleOPLSDA           = c("none", "mean", "meanSD", "meanSD2")[4],
    dataNormalize             = c("none", "DilutionMarker", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile")[2],
    refSample                 = NULL,
    groupSample               = NULL,
    reference_method          = c("mean", "median")[1],
    filterMaxRSD              = 30,
    filterMaxRSD_by           = c("SQC", "EQC", "both")[2],
    filterMaxVarSD            = 10,
    outliers                  = NULL
) {

  # Check if data is from "performDataQualityCheck" function
  if (raw_data$FunctionOrigin != "performDataQualityCheck") {
    stop("The 'raw_data' must be from the 'performDataQualityCheck' function for quality check.")
  }

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # --------------------------------------- Parameter Validation (start) ----------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  check_parameters <- function(filterMissing,
                               filterMissing_by_group,
                               filterMissing_includeQC,
                               denMissing,
                               driftBatchCorrection,
                               spline_smooth_param,
                               spline_smooth_param_limit,
                               log_scale,
                               min_QC,
                               dataNormalize,
                               refSample,
                               groupSample,
                               reference_method,
                               dataTransform,
                               dataScalePCA,
                               dataScaleOPLSDA,
                               filterMaxRSD,
                               filterMaxRSD_by,
                               filterMaxVarSD,
                               outliers) {
    errors <- character()

    # Check numeric parameters
    if (!is.numeric(filterMissing)) {
      errors <- c(errors, "filterMissing: Not numeric.")
    } else if (filterMissing < 1 || filterMissing > 100) {
      errors <- c(errors, "filterMissing: Out of range [1, 100].")
    }

    if (!is.logical(filterMissing_by_group)) {
      errors <- c(errors, "filterMissing_by_group: Must be logical.")
    }

    if (!is.logical(filterMissing_includeQC)) {
      errors <- c(errors, "filterMissing_includeQC: Must be logical.")
    }

    if (!is.numeric(denMissing)) {
      errors <- c(errors, "denMissing: Not numeric.")
    } else if (denMissing < 1 || denMissing > 100) {
      errors <- c(errors, "denMissing: Out of range [1, 100].")
    }

    if (!is.logical(driftBatchCorrection)) {
      errors <- c(errors, "driftBatchCorrection: Must be logical.")
    }

    if (!is.numeric(spline_smooth_param)) {
      errors <- c(errors, "spline_smooth_param: Not numeric.")
    } else if (spline_smooth_param < 0 || spline_smooth_param > 1) {
      errors <- c(errors, "spline_smooth_param: Out of range [0, 1].")
    }

    if (!is.numeric(spline_smooth_param_limit)) {
      errors <- c(errors, "spline_smooth_param_limit: Not numeric. Must be a numeric vector of 2 values (e.g., `c(-1.5, 1.5)` (default)).")
    }

    if (!is.logical(log_scale)) {
      errors <- c(errors, "log_scale: Must be logical.")
    }

    if (!is.numeric(min_QC)) {
      errors <- c(errors, "min_QC: Not numeric.")
    }

    if (!is.null(filterMaxRSD)) {
      if (!is.numeric(filterMaxRSD)) {
        errors <- c(errors, "filterMaxRSD: Not numeric. Suggestions are 20 (20% for LC-MS analyzed) and 30 (30% for GC-MS analyzed).")
      } else if (filterMaxRSD < 1 || filterMaxRSD > 100) {
        errors <- c(errors, "filterMaxRSD: Out of range [1, 100].")
      }
    }

    if (!is.null(filterMaxVarSD)) {
      if (!is.numeric(filterMaxVarSD)) {
        errors <- c(errors, "filterMaxVarSD: Not numeric. Must be between 1 and 100.")
      } else if (filterMaxVarSD < 1 || filterMaxVarSD > 100) {
        errors <- c(errors, "filterMaxVarSD: Out of range [1, 100].")
      }
    }

    # Check choices
    allowed_scales <- c("none", "DilutionMarker", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile")
    if (!is.logical(dataNormalize)) {
      if (!(dataNormalize %in% allowed_scales)) {
        errors <- c(errors, "dataNormalize: Must either be none, DilutionMarker, sum, median, PQN1, PQN2, groupPQN, quantile.")
      }
    }

    allowed_scales <- c("none", "log2", "log10", "sqrt", "cbrt", "vsn")
    if (!is.logical(dataTransform)) {
      if (!(dataTransform %in% allowed_scales)) {
        errors <- c(errors, "dataTransform: Must either be none, log2, log10, sqrt, cbrt, vsn. log2 and log10 applies log2 and log10 transformation, sqrt and cbrt applies square root and cube root transformation, and vsn applies variance stabilizing normalization.")
      }
    }

    allowed_scales <- c("none", "mean", "meanSD", "meanSD2")
    if (!is.null(dataScalePCA)) {
      if (!(dataScalePCA %in% allowed_scales)) {
        errors <- c(errors, "dataScalePCA: Must be either NULL, mean, meanSD, or meanSD2.")
      }
    }
    if (!is.null(dataScaleOPLSDA)) {
      if (!(dataScaleOPLSDA %in% allowed_scales)) {
        errors <- c(errors, "dataScaleOPLSDA: Must be either NULL, mean, meanSD, or meanSD2.")
      }
    }

    allowed_scales <- c("mean", "median")
    if (!is.null(reference_method)) {
      if (!(reference_method %in% allowed_scales)) {
        errors <- c(errors, "reference_method: Must be either mean or median.")
      }
    }

    allowed_scales <- c("SQC", "EQC", "both")
    if (!is.null(filterMaxRSD_by)) {
      if (!(filterMaxRSD_by %in% allowed_scales)) {
        errors <- c(errors, "filterMaxRSD_by: Must be either SQC, EQC, or both.")
      }
    }

    # Check vector
    if (!is.null(outliers)) {
      if (!is.vector(outliers)) {
        errors <- c(errors, "outliers: Must be a vector of biological samples and QC samples identified as outliers.")
      }
    }

    return(errors)
  }

  # Validate parameters
  validation_errors <- check_parameters(filterMissing,
                                        filterMissing_by_group,
                                        filterMissing_includeQC,
                                        denMissing,
                                        driftBatchCorrection,
                                        spline_smooth_param,
                                        spline_smooth_param_limit,
                                        log_scale,
                                        min_QC,
                                        dataNormalize,
                                        refSample,
                                        groupSample,
                                        reference_method,
                                        dataTransform,
                                        dataScalePCA,
                                        dataScaleOPLSDA,
                                        filterMaxRSD,
                                        filterMaxRSD_by,
                                        filterMaxVarSD,
                                        outliers)

  if (length(validation_errors) > 0) {
    stop(paste("Invalid parameters:", paste(validation_errors,
                                            collapse = "\n")))
  }

  # If all parameters are okay, proceed with code.

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ---------------------------------------- Parameter Validation (end) -----------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  listPreprocessed            <- base::list() # Create empty list to contain all results. A list can store any data (vector, list, data frame, matrix, etc.)
  # This is to make sure that the other analysis uses results only from this function
  listPreprocessed$FunctionOrigin <- "performPreprocessingPeakData"
  df_dimensions               <- base::data.frame()
  listPreprocessed$Dimensions <- base::rbind(df_dimensions, c("preprocessing_step", "n_samples", "n_features")) %>% `colnames<-`(NULL) # Update list

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------------ Load the data (start) --------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # # Transpose the data
  data_transposed <- as.data.frame(t(raw_data$raw_data)) %>%
    setNames(as.character(.[1, ])) %>% .[-1, ] %>% # Set column names from 1st row (features, etc.)
    # Arrange by Injection Sequence
    dplyr::mutate(Injection = as.numeric(Injection)) %>% # Convert to numeric so it will be arranged correctly in the line below
    dplyr::arrange(Injection) # Sort by injection, needed in drift-correction

  message("This might take more or less 5 minutes due to a large number of data preprocessing methods included in this function. This is especially true when there is more than 1 modelName in driftCorrection.")
  message("While processing, RStudio maybe be used (unless maybe there is another separate session/window open). Feel free to sip a coffee!")

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------------ Remove outliers (start) ------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Remove user-defined outliers
  if (!is.null(outliers)) {

    message("The 'outliers' parameter is not empty, removing biological and/or QC samples.")

    # Identify which outliers are present in the 'Sample' column
    removed_outliers   <- outliers[ outliers %in% data_transposed$Sample]
    not_found_outliers <- outliers[!outliers %in% data_transposed$Sample]

    # Filter outliers
    data_transposed <- data_transposed %>%
      dplyr::filter(!(Sample %in% outliers))

    # Print them
    cat("✅ Removed outliers:\n")
    print(removed_outliers %>% data.frame() %>% `colnames<-`(NULL))

    cat("\n⚠️ Outliers NOT found in the 'Sample' column (probably already removed):\n")
    print(not_found_outliers %>% data.frame() %>% `colnames<-`(NULL))

    listPreprocessed$outliers             <- outliers # Update list
    listPreprocessed$outliers_removed     <- removed_outliers # Update list
    listPreprocessed$outliers_not_removed <- not_found_outliers # Update list

  } else {
    removed_outliers <- data.frame() # Quickest way to solve as of the moment
  }

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------------ Remove outliers (end) --------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # List of metadata indices
  indices_qc      <- data_transposed$Group %in% c("SQC", "EQC", "QC")    # Indices for all QCs
  indices_non_qc  <- !indices_qc                                         # Indices for Biological Samples, i.e., non-QC

  # Saving metadata in 1 data frame
  listPreprocessed$Metadata <- data.frame(
    Samples             = data_transposed$Sample,
    SubjectID           = data_transposed$SubjectID %>% as.numeric(),
    TechnicalReplicates = data_transposed$Replicate,
    Group               = data_transposed$Group,
    Group_              = gsub("SQC|EQC", "QC", data_transposed$Group), # This is another set of groups where SQC and EQC are converted to QC (for plot purposes only)
    Batches             = data_transposed$Batch %>% as.numeric(),
    InjectionSequence   = data_transposed$Injection %>% as.numeric(),
    DilutionMarker      = data_transposed$DilutionMarker %>% as.numeric(),
    Response            = data_transposed$Response %>% as.numeric()
  )

  # Extract Samples x Features data only
  metadata_columns <- c("Sample", "SubjectID", "Replicate", "Group", "Batch", "Injection", "DilutionMarker", "Response") # The columns we don't need anymore
  df <- data_transposed %>% dplyr::select(-one_of(metadata_columns)) # remove metadata columns columns

  listPreprocessed$Dimensions <- base::rbind(listPreprocessed$Dimensions,
                                             # Add the removed outliers. I did this to avoid creating multiple variables
                                             c("Original", dim(df)[1] + length(removed_outliers), dim(df)[2])) # Update list
  listPreprocessed$Dimensions <- base::rbind(listPreprocessed$Dimensions, c("Removed outliers",
                                                                            dim(df)[1], dim(df)[2])) # Update list

  listPreprocessed$SamplesXFeatures <- df # Update list

  df[df == 0] <- NA # Replace 0 with NA
  num_na      <- sum(is.na(df), na.rm = TRUE) # Count NAs
  df          <- `rownames<-`(df, data_transposed$Sample) # set row names to sample IDs
  df[]        <- lapply(df, as.numeric) # Convert to numeric
  df          <- df[, colSums(is.na(df)) < nrow(df)] # Remove features if all NAs

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------- Filtering Out Features/Metabolites with Large RSD (start) -------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # In filterMissing
  if (filterMissing_by_group) {
    # Group-based filtering (your existing logic)
    if (filterMissing_includeQC) {
      # Compute missing percentages per group (including QC)
      missing_by_group <- base::sapply(listPreprocessed$Metadata$Group, function(g) {
        rows <- listPreprocessed$Metadata$Group == g
        base::colMeans(is.na(df[rows, , drop = FALSE]))
      }) %>% t() %>% .[!duplicated(.), ]
      # Remove feature if missingness is >= threshold in all groups
      remove_features <- base::apply(missing_by_group, 2, function(x) all(x >= filterMissing / 100))
    } else {
      # Compute missing percentages per group (excluding SQC, EQC, and QC)
      missing_by_group <- base::sapply(listPreprocessed$Metadata$Group[!listPreprocessed$Metadata$Group %in% c("SQC", "EQC", "QC")], function(g) {
        rows <- listPreprocessed$Metadata$Group == g
        base::colMeans(is.na(df[rows, , drop = FALSE]))
      }) %>% t() %>% .[!duplicated(.), ]
      # Remove feature if missingness is >= threshold in all groups
      remove_features <- base::apply(missing_by_group, 2, function(x) all(x >= filterMissing / 100))
    }
  } else {
    # Overall filtering (ignoring groups)
    if (filterMissing_includeQC) {
      # Calculate missing percentages across all samples (including QC)
      missing_overall <- base::colMeans(is.na(df))
      # Remove features if missingness is >= threshold
      remove_features <- missing_overall >= filterMissing / 100
    } else {
      # Calculate missing percentages excluding SQC, EQC, and QC samples
      rows_to_include <- !listPreprocessed$Metadata$Group %in% c("SQC", "EQC", "QC")
      missing_overall <- base::colMeans(is.na(df[rows_to_include, , drop = FALSE]))
      # Remove features if missingness is >= threshold
      remove_features <- missing_overall >= filterMissing / 100
    }

    # Save removed features as data frame
    listPreprocessed$data_Removed_Columns_in_filteredMissing <- df[, remove_features, drop = FALSE]

    # Filter the dataset
    df <- df[, !remove_features]
  }

  # Filter the dataset
  df <- df[, !remove_features]

  listPreprocessed$data_filteredMissing <- df # Update list

  # Update results
  listPreprocessed$Dimensions            <- base::rbind(listPreprocessed$Dimensions, c("Missing Percent Greater", dim(df)[1], dim(df)[2])) # Update list
  listPreprocessed$NumberOfMissingValues <- num_na

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------- Filtering Out Features/Metabolites with Large RSD (end) --------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------- Missing Value Imputation (start) -------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Replace missing values with a fraction of the minimum positive value
  df[]        <- base::lapply(df, function(x) {
    min_val   <- base::min(x[x > 0], na.rm = TRUE)
    base::replace(
      x      = x,                   # vector
      list   = is.na(x),            # index vector
      values = min_val / denMissing # replacement values
    )
  })

  listPreprocessed$data_no_NA <- df # Update list

  message("✅ Employed: All 0's have been replaced with NA's.")
  message(paste0("✅ Employed: Removal of features where at least ", filterMissing, "% is missing in each Group."))
  message(paste0("✅ Employed: Missing values replaced with 1/", denMissing, " of the smallest positive peak."))

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------- Missing Value Imputation (end) ---------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------ Drift- and Batch-Correction (start) ------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("Employing: Drift correction then batch correction. This only outputs 1 data frame.")

  if (isTRUE(driftBatchCorrection)) {
    df <- pmp::QCRSC(
      df       = df,
      order    = listPreprocessed$Metadata$InjectionSequence,
      batch    = listPreprocessed$Metadata$Batches,
      # Replace SQC and EQC with QC and back to vector
      classes  = data.frame(Groups = listPreprocessed$Metadata$Group) %>% dplyr::mutate(Groups = ifelse(Groups %in% c("SQC", "EQC"), "QC", Groups)) %>% as.vector() %>% .$Groups,
      spar     = spline_smooth_param, # 0 Spline smoothing parameter. Should be in the range 0 to 1. If set to 0 it will be estimated using leave-one-out cross-validation.
      log      = log_scale,
      minQC    = min_QC,
      qc_label = "QC",
      spar_lim = spline_smooth_param_limit
    ) %>% t() %>% as.data.frame()

    listPreprocessed$data_driftBatchCorrected <- df # Update list

    if (spline_smooth_param == 0) {
      message("✅ Employed: Drift correction then batch correction with Spline smoothing parameter (leave-one-out cross-validation).")
    } else {
      message("✅ Employed: Drift correction then batch correction.")
    }

  } else { # if driftBatchCorrection is not TRUE
    df <- df
    listPreprocessed$data_driftBatchCorrected <- df # Update list

    message("❌ NOT Employed: Drift correction then batch correction with Spline smoothing parameter (leave-one-out cross-validation).")
  }

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------ Drift- and Batch-Correction (end) --------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("Employing: 2nd Replacing of missing values using 1/denMissing, because drift correction and batch correction likely produced NAs.")


  # Replace missing values 2, because drift correction and batch correction produced NAs
  df[]        <- base::lapply(df, function(x) { # Replace missing values with a fraction of the minimum positive value
    min_val   <- base::min(x[x > 0], na.rm = TRUE)
    base::replace(
      x      = x,                     # vector
      list   = is.na(x),              # index vector
      values = min_val / denMissing   # replacement values ## Replace with minimum value only since only some values were NAs
    )
  })

  listPreprocessed$data_driftBatchCorrected_noNa <- df # Update list

  message("✅ Employed: 2nd Replacing of missing values introduced by drift- and batch-correction.")

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------- Plot Before and After Drift- and Batch-Correction (start) -------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Before and After Drift- and Batch-correction with trend lines on QCs
  create_batch_drift_plots_with_trends <- function(df_before, df_after, classes, batch, indexes, qc_label = "QC") {

    create_single_plot_with_trend <- function(feature_idx, feature_name) {

      plot_data <- base::data.frame(
        sample_id = 1:nrow(df_before),
        before    = base::log10(df_before[, feature_idx]),   # log10 transformation
        after     = base::log10(df_after[, feature_idx]),    # log10 transformation
        class     = classes,
        batch     = base::as.factor(batch)
      )

      plot_data_long <- plot_data %>%
        tidyr::pivot_longer(cols      = c(before, after),
                            names_to  = "correction",
                            values_to = "intensity") %>%
        dplyr::mutate(correction = base::factor(correction, levels = c("before", "after")))

      # Add trend lines for QC samples
      qc_data <- plot_data_long %>% dplyr::filter(class == qc_label)

      p <- ggplot2::ggplot(plot_data_long, aes(x = sample_id, y = intensity)) +
        ggplot2::geom_point(aes(color = class, shape = batch), size = 2, alpha = 0.7) +
        # Add trend line for QC samples
        ggplot2::geom_smooth(data = qc_data, aes(x = sample_id, y = intensity),
                             method = "loess", se = TRUE, alpha = 0.3,
                             color = "darkred", linetype = "dashed") +
        ggplot2::facet_wrap(~ correction, scales = "free_y",
                            labeller = labeller(correction = c(before = "Before Correction",
                                                               after = "After Correction"))) +
        ggplot2::scale_color_manual(values = c("QC" = "red",
                                               "Sample" = "blue",
                                               "Blank" = "green")) +
        ggplot2:: scale_shape_manual(values = c("1" = 16, "2" = 17, "3" = 15, "4" = 18,
                                                "5" = 8, "6" = 9, "7" = 10, "8" = 11)) +
        ggplot2::labs(title = paste("Feature:", feature_name),
                      x     = "Sample Index",
                      y     = "log10(Intensity)",
                      color = "Class",
                      shape = "Batch") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title   = element_text(size = 12, hjust = 0.5),
          axis.text    = element_text(size = 10),
          axis.title   = element_text(size = 11),
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 9),
          strip.text   = element_text(size = 10, face = "bold")
        )

      return(p)
    }

    plot_list     <- base::list()
    feature_names <- base::colnames(df_before)[indexes]

    for (i in seq_along(indexes)) {
      plot_list[[i]] <- create_single_plot_with_trend(indexes[i], feature_names[i])
    }

    return(plot_list)
  }

  # Select 6 random features
  selected_indexes <- base::sample(x = 1:ncol(df), size = 6, replace = FALSE)

  # Create the enhanced plots with trend lines and log10 transformation
  beforeAfterBatchDrift <- create_batch_drift_plots_with_trends(
    df_before = listPreprocessed$data_no_NA,
    df_after  = df,
    classes   = listPreprocessed$Metadata$Group_,
    batch     = listPreprocessed$Metadata$Batches,
    indexes   = selected_indexes,
    qc_label  = "QC"
  )

  # Determine the main title based on batch information
  main_title <- if (length(unique(listPreprocessed$Metadata$Batches)) == 1) {
    "Plot of 6 Random Features/Metabolites Before and After Drift-Correction"
  } else {
    "Plot of 6 Random Features/Metabolites Before and After Drift- and Batch-Correction"
  }

  # Arrange plots in a grid with main title and suppress all messages from this block
  plot_beforeAfterBatchDrift <- base::suppressMessages({
    gridExtra::grid.arrange(
      beforeAfterBatchDrift[[1]],
      beforeAfterBatchDrift[[2]],
      beforeAfterBatchDrift[[3]],
      beforeAfterBatchDrift[[4]],
      beforeAfterBatchDrift[[5]],
      beforeAfterBatchDrift[[6]],
      ncol = 2, nrow = 3,
      top = main_title
    )
  })

  # Store in the list
  listPreprocessed$plot_beforeAfterDriftBatchCorrection <- plot_beforeAfterBatchDrift

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------- Plot Before and After Drift- and Batch-Correction (end) --------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ---------------------------------------- Data Transformation (start) ----------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("Employing: Data transformation.")

  # Assume `dataTransform` is one of: "none", "log2", "log10", "sqrt", "cbrt"
  if (dataTransform == "log2") {
    df <- df - base::min(df) + 1  # Ensure all values are positive for log2
    df <- base::log2(df) %>% as.data.frame()
    message("✅ Employed: Log2 transformation.")

  } else if (dataTransform == "log10") {
    df <- df - min(df) + 1  # Ensure all values are positive for log10
    df <- base::log10(df) %>% as.data.frame()
    message("✅ Employed: Log10 transformation.")

  } else if (dataTransform == "sqrt") {
    df <- df - min(df) + 1  # Ensure all values are non-negative for sqrt
    df <- base::sqrt(df) %>% as.data.frame()
    message("✅ Employed: Square root transformation.")

  } else if (dataTransform == "cbrt") {
    # df <- sign(df) * abs(df)^(1/3) %>% as.data.frame()  # Cube root handles negative values
    df <- df - base::min(df) + 1  # Ensure all values are non-negative
    df <- df^(1/3) %>% as.data.frame()
    message("✅ Employed: Cube root transformation.")

  } else if (dataTransform == "vsn") {
    df <- t(df) %>% # Transpose to Features x Samples format as requirement of vsn
      vsn::justvsn(df) %>%
      t() %>% as.data.frame() # Transpose back to Samples x Features format
    message("✅ Employed: Variance Stabilizing Normalization.")

  } else {
    message("❌  NOT Employed: 'none' was selected, hence, no data transformation was carried out.")
  }

  listPreprocessed$data_transformed <- df # Update list

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------------ Data Transformation (end) ----------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------------- Data Scaling (start) --------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Scaling function
  scale_data <- function(data, method, data_type) {
    message(paste0("Employing: Data scaling for ", data_type, "."))

    if (is.null(method)) {
      message(paste0("❌ NOT Employed: Data scaling for ", data_type, " data."))
      return(data)
    }

    scaled_data <- base::switch(method,
                                "mean" = scale(data, center = TRUE, scale = FALSE),
                                "meanSD" = scale(data, center = TRUE, scale = TRUE),
                                "meanSD2" = scale(data, center = TRUE, scale = sqrt(apply(data, 2, sd))),
                                data  # fallback if method not recognized
    ) %>% as.data.frame()

    message(paste0("✅ Employed: ", method, " Data scaling for ", data_type, " data."))
    return(scaled_data)
  }

  # Apply scaling for both datasets
  listPreprocessed$data_scaledPCA <- scale_data(df, dataScalePCA, "PCA")
  listPreprocessed$data_scaledOPLSDA <- scale_data(df, dataScaleOPLSDA, "OPLS-DA")

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------- Data Scaling (end) ---------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ---------------------------------------- Data Normalization (start) -----------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("Employing: Data Normalization.")

  # Function to apply normalization to a dataset
  normalize_data <- function(df_scaled) {
    if (dataNormalize == "none") {
      message("❌ NOT Employed: Data Normalization.")
      return(df_scaled)
    }

    if (dataNormalize == "DilutionMarker") {
      sg_values <- base::as.numeric(data_transposed$DilutionMarker)
      if (all(is.na(sg_values))) {
        message("⚠️ Dilution Marker values are all NA. Normalization by Sum will be performed instead.")
        sample_sums <- base::rowSums(df_scaled, na.rm = TRUE)
        df_scaled <- df_scaled / sample_sums
      } else {
        df_scaled[indices_non_qc, ] <- base::sweep(df_scaled[indices_non_qc, ], 1, sg_values, "/")
        df_scaled[indices_qc, ] <- base::sweep(df_scaled[indices_qc, ], 1, base::rowSums(df_scaled[indices_qc, ], na.rm = TRUE), "/")
        message("✅ Employed: Normalization using Dilution Marker values.")
      }
    } else if (dataNormalize == "sum") {
      sample_sums <- base::rowSums(df_scaled, na.rm = TRUE)
      df_scaled <- df_scaled / sample_sums
      message("✅ Employed: Normalization by Sum.")
    } else if (dataNormalize == "median") {
      sample_medians <- base::apply(df_scaled, 1, median, na.rm = TRUE)
      df_scaled <- df_scaled / sample_medians
      message("✅ Employed: Normalization by Median.")
    } else if (dataNormalize == "PQN1") {
      reference_spectrum <- base::apply(df_scaled, 2, median, na.rm = TRUE)
      quotients <- df_scaled / reference_spectrum
      median_quotients <- base::apply(quotients, 1, median, na.rm = TRUE)
      df_scaled <- df_scaled / median_quotients
      message("✅ Employed: Probabilistic Quotient Normalization (PQN) according to global median approach.")
    } else if (dataNormalize == "PQN2") {
      reference_spectrum <- df_scaled[refSample, ] %>% t()
      quotients <- base::sweep(df_scaled, 2, reference_spectrum, "/")
      median_quotients <- apply(quotients, 1, median, na.rm = TRUE)
      df_scaled <- base::sweep(df_scaled, 1, median_quotients, "/")
      message(paste0("✅ Employed: Probabilistic Quotient Normalization (PQN) using '", refSample, "' as the reference sample in the normalization."))
    } else if (dataNormalize == "groupPQN") {
      if (is.null(groupSample)) {
        stop("Please specify 'groupSample' for groupPQN normalization.")
      }
      pooled_indices <- if (groupSample == "both") {
        data_transposed$Group %in% c("SQC", "EQC", "QC")
      } else {
        data_transposed$Group == groupSample
      }
      if (sum(pooled_indices) == 0) {
        stop("No samples found for the specified 'groupSample'.")
      }
      reference_spectrum <- base::apply(df_scaled[pooled_indices, , drop = FALSE], 2, median, na.rm = TRUE)
      quotients <- df_scaled / reference_spectrum
      median_quotients <- base::apply(quotients, 1, median, na.rm = TRUE)
      df_scaled <- df_scaled / median_quotients
      message("✅ Employed: Group Probabilistic Quotient Normalization (groupPQN) using pooled group: ", groupSample)
    } else if (dataNormalize == "quantile") {
      df_scaled <- pmp::pqn_normalisation(
        df_scaled,
        classes = base::data.frame(Groups = listPreprocessed$Metadata$Group) %>%
          dplyr::mutate(Groups = ifelse(Groups %in% c("SQC", "EQC"), "QC", Groups)) %>%
          dplyr::pull(Groups),
        qc_label = "QC",
        ref_mean = NULL,
        qc_frac = 0,
        sample_frac = 0,
        ref_method = reference_method
      ) %>% t() %>% as.data.frame()
      message("✅ Employed: Quantile Normalization.")
    }

    return(df_scaled)
  }

  # Apply normalization to both datasets
  listPreprocessed$data_normalizedPCA <- normalize_data(listPreprocessed$data_scaledPCA)
  listPreprocessed$data_normalizedOPLSDA <- normalize_data(listPreprocessed$data_scaledOPLSDA)

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ----------------------------------------- Data Normalization (end) ------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ---------------------- Data Filtering Using Relative Standard Deviation (RSD) (start) -----------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("Employing: Filtering of uninformative features using RSD (Relative Standard Deviation).")

  # Function to apply RSD filtering
  filter_by_rsd <- function(data, data_type) {
    if (!is.null(filterMaxRSD)) {
      df <- pmp::filter_peaks_by_rsd(
        data,
        max_rsd = filterMaxRSD,
        class = if (filterMaxRSD_by == "both") {
          base::data.frame(Groups = listPreprocessed$Metadata$Group) %>%
            dplyr::mutate(Groups = ifelse(Groups %in% c("SQC", "EQC"), "QC", Groups)) %>%
            dplyr::pull(Groups)
        } else {
          listPreprocessed$Metadata$Group
        },
        qc_label = ifelse(filterMaxRSD_by == "both", "QC", filterMaxRSD_by)
      ) %>% t() %>% as.data.frame()

      # Update dimensions
      listPreprocessed$Dimensions <<- base::rbind(listPreprocessed$Dimensions,
                                                  c(paste0("Relative SD (for ", data_type, "data)"),
                                                    dim(df)[1], dim(df)[2])) %>% `colnames<-`(NULL)

      message(paste0("✅ Employed: Filtering of uninformative features using RSD = ",
                     filterMaxRSD, "%. (for ", data_type, " data)"))
      return(df)
    } else {
      message(paste0("❌ NOT Employed: Filtering of uninformative features using RSD = ",
                     filterMaxRSD, "%. (for ", data_type, " data)"))
      return(data)
    }
  }

  # Apply RSD filtering to both datasets
  listPreprocessed$data_filteredMaxRSD_PCA <- filter_by_rsd(listPreprocessed$data_normalizedPCA, "PCA")
  listPreprocessed$data_filteredMaxRSD_OPLSDA <- filter_by_rsd(listPreprocessed$data_normalizedOPLSDA, "OPLS-DA")

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ---------------------- Data Filtering Using Relative Standard Deviation (RSD) (end) -------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------ Data Filtering on Features with Low Variability Across Groups (start) --------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("Employing: Filtering features with low variability across samples.")

  # Function to apply variance filtering
  filter_by_variance <- function(data, data_type) {
    if (!is.null(filterMaxVarSD)) {
      df <- data %>%
        { .[, apply(., 2, sd) > quantile(apply(., 2, sd), filterMaxVarSD/100)] }

      # Update dimensions
      listPreprocessed$Dimensions <<- base::rbind(listPreprocessed$Dimensions,
                                                  c(paste0("Low Variance (for ", data_type, "data)"),
                                                    dim(df)[1], dim(df)[2])) %>%
        `colnames<-`(NULL)

      message(paste0("✅ Employed: Filtered uninformative features with low variability across samples. Used threshold = ",
                     filterMaxVarSD, "th percentile (for ", data_type, " data)."))
      return(df)
    } else {
      message(paste0("❌ NOT Employed: Filtering of uninformative features with low variability across samples (for ",
                     data_type, " data)."))
      return(data)
    }
  }

  # Apply variance filtering to both datasets
  listPreprocessed$data_filteredMaxVarSD_PCA <- filter_by_variance(listPreprocessed$data_filteredMaxRSD_PCA, "PCA")
  listPreprocessed$data_filteredMaxVarSD_OPLSDA <- filter_by_variance(listPreprocessed$data_filteredMaxRSD_OPLSDA, "OPLS-DA")

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------- Data Filtering on Features with Low Variability Across Groups (end) --------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  message("✅✅✅ Data preprocessing is complete. The data is now ready for any statistical analysis based on selected data preprocessing techniques!")

  message("⚠️⚠️⚠️ In the results, despite choosing 'none' for some methods (for example no normalization or no filtering vis RSD), there will still be results as such. However, this is for recording purposes only, and if chosen to have no normalization or filtering or any other step, there will be performing of such data preprocessing step. The data is the same as the previous step in this process.")
  flush.console()

  return(listPreprocessed)
}
