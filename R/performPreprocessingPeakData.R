#' Perform Data Preprocessing on a Raw Data (csv)
#'
#' @description
#' This function performs data preprocessing techniques in bioinformatics to prepare the data for downstream analysis.
#' The data preprocessing includes data filtration, drift correction, data normalization, transformation, scaling,
#' and removal of known or identified outliers.
#'
#' @param raw_data A csv file having the following characteristics:
#'   \itemize{
#'     \item 1st row: Contains QC and the groups (with/out spaces).
#'     \item 2nd row: Contains SQC, EQC, and the same groups as in the 1st row.
#'     \item 3rd row: The sample names, ideally shortened versions to make it easier to view in visualizations, without spaces.
#'     \item 4th row: The batch numbers.
#'     \item 5th row: The injection sequence.
#'     \item 6th row: The osmolality values or specific gravity.
#'     \item Others below:
#'     \item Missing values are best to be left or encoded as 0s.
#'     }
#' @param filterMissing Numeric. Minimum %missing in all sample groups required to remove feature
#' @param denMissing Numeric. A value to be used in missing value imputation. A denominator in 1/denMissing.
#' @param driftBatchCorrection Boolean. If `TRUE` (default), perform multi-batch alignment to merge features artificially split between batches.
#' @param filterMaxRSD Numeric. The threshold for the Relative Standard Deviation (RSD). Suggestions below. RSD is used to assess and filter out unreliable metabolites, ensuring that only high-quality, reproducible data are used for statistical and biological interpretation. In metabolomics, RSD is commonly calculated for each metabolite across Quality Control (QC) samples. QC samples are repeated pooled samples meant to assess the analytical precision and reproducibility of the instrument and the experimental setup.
#'   \itemize{
#'     \item 20: for LC-MS analyzed data
#'     \item 30: for GC-MS analyzed data
#'     }
#' @param filterMaxRSD_by String. Choose below where to apply the RSD filtering in `filterMaxRSD`.
#'   \itemize{
#'     \item "SQC": Apply on Sample QCs only.
#'     \item "EQC": Apply on Extract QCs only.
#'     \item "both": Apply on both.
#'     }
#'     Defaults to "EQC".
#' @param filterMaxVarSD Numeric. Remove nth percentile of features with the lowest variability (e.g., 10 for 10%). Set to `NULL` to skip this filtering.
#' @param dataNormalize String. Perform data normalization. Options are:
#'   \itemize{
#'     \item "none": No normalization.
#'     \item "SpecificGravity": Normalization using specific gravity, if provided. Otherwise, normalized by 'sum'.
#'     \item "sum": Normalization by sum.
#'     \item "median": Normalization by median.
#'     \item "PQN1": Probabilistic Quotient Normalization (PQN) according to global median approach."
#'     \item "PQN2": Probabilistic Quotient Normalization (PQN) using a reference sample indicated in `refSample` as the reference sample in the normalization."
#'     \item "groupPQN": Group Probabilistic Quotient Normalization using pooled group indicated in `groupSample`.
#'     \item "quantile": Normalization by quantile,
#'     }
#'     Defaults to 'SpecificGravity' is present, otherwise, normalization by 'sum' will be performed.
#' @param refSample String. The reference sample in the case of normalization method "PQN2". Must be in the samples, and is not part of "outliers" vector c("SQC", "EQC", "both")
#' @param groupSample String. Used in the `groupPQN.` NULL if not "groupPQN". Default to "EQC" if "groupPQN" Other choice is "SQC".
#' @param dataTransform String. A transformation method to transform the data after `dataNormalize`.
#'   \itemize{
#'     \item "none": No data transformation is done.
#'     \item "log10": Perform log10 transformation.
#'     \item "sqrt": Perform square root transformation.
#'     \item "cbrt": Perform cube root transformation.
#'     }
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
#'     \item "mean": Mean-centered only
#'     \item "meanSD": Mean-centered and divided by SD of each feature
#'     \item "meanSD2": Mean-centered and divided by the square root of SD of each feature. Also called pareto-scaling.
#'     }
#'     Defaults to "meanSD2".
#' @param outliers A vector of biological samples and/or QC that are considered as outliers. Example format is "c('Sample1', 'Sample2', 'QC1', 'QC2', ...)".
#'
#' @returns A list of results from all of the tests done (e.g., batch-corrected data, normalized data, transformed data, scaled data, etc.)
#' @export
#'
#' @examples
#' \dontrun{
#' # Using the default parameters
#' mydata <- performPreprocessingPeakData(
#'  raw_data = a_csv_file
#' )
#' }

performPreprocessingPeakData <- function(
    raw_data,
    filterMissing        = 20,
    denMissing           = 5,
    driftBatchCorrection = TRUE,
    filterMaxRSD         = 30,
    filterMaxRSD_by      = c("SQC", "EQC", "both")[2],
    filterMaxVarSD       = 10,
    dataNormalize        = c("none", "SpecificGravity", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile")[2],
    refSample            = NULL,
    groupSample          = NULL,
    dataTransform        = c("none", "log10", "sqrt", "cbrt")[2],
    dataScalePCA         = c("none", "mean", "meanSD", "meanSD2")[3],
    dataScaleOPLSDA      = c("none", "mean", "meanSD", "meanSD2")[4],
    outliers             = NULL
) {


  # Parameter Validation Function
  check_parameters <- function(filterMissing,
                               denMissing,
                               driftBatchCorrection,
                               filterMaxRSD,
                               filterMaxVarSD,
                               dataNormalize,
                               dataTransform,
                               dataScalePCA,
                               dataScaleOPLSDA,
                               outliers) {
    errors <- character()

    # Check numeric parameters
    if (!is.numeric(filterMissing)) {
      errors <- c(errors, "filterMissing: Not numeric")
    } else if (filterMissing < 1 || filterMissing > 100) {
      errors <- c(errors, "filterMissing: Out of range [1, 100].")
    }

    if (!is.numeric(denMissing)) {
      errors <- c(errors, "denMissing: Not numeric")
    } else if (denMissing < 1 || denMissing > 100) {
      errors <- c(errors, "denMissing: Out of range [1, 100]")
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
    allowed_scales <- c("none", "SpecificGravity", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile")
    if (!is.logical(dataNormalize)) {
      if (!(dataNormalize %in% allowed_scales)) {
        errors <- c(errors, "dataNormalize: Must either be none, SpecificGravity, sum, median, PQN1, PQN2, groupPQN, quantile.")
      }
    }

    allowed_scales <- c("none", "log10", "sqrt", "cbrt")
    if (!is.logical(dataTransform)) {
      if (!(dataTransform %in% allowed_scales)) {
        errors <- c(errors, "dataTransform: Must either be none, log10, sqrt, cbrt. log10 applies log10 transformation, sqrt and cbrt applies square root and cube root transformation.")
      }
    }

    # Check choices
    allowed_scales <- c("none", "mean", "meanSD", "meanSD2")
    if (!is.null(dataScalePCA)) {
      if (!(dataScalePCA %in% allowed_scales)) {
        errors <- c(errors, "dataScalePCA: Must be either NULL, mean, meanS', or meanSD2.")
      }
    }

    if (!is.null(dataScaleOPLSDA)) {
      if (!(dataScaleOPLSDA %in% allowed_scales)) {
        errors <- c(errors, "dataScaleOPLSDA: Must be either NULL, mean, meanS', or meanSD2.")
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

  library(tidyverse)

  # Validate parameters
  validation_errors <- check_parameters(filterMissing,
                                        denMissing,
                                        driftBatchCorrection,
                                        filterMaxRSD,
                                        filterMaxVarSD,
                                        dataNormalize,
                                        dataTransform,
                                        dataScalePCA,
                                        dataScaleOPLSDA,
                                        outliers)

  if (length(validation_errors) > 0) {
    stop(paste("Invalid parameters:", paste(validation_errors,
                                            collapse = "\n")))
  }

  # If all parameters are okay, proceed with code.

  message("This might take more or less 5 minutes due to a large number of data preprocessing methods included in this function. This is especially true when there is more than 1 modelName in driftCorrection.")
  message("While processing, R cannot be used (unless maybe there is another separate session). Feel free to do other tasks.")
  flush.console() # Force message to be displayed

  listPreprocessed            <- list() # Create empty list to contain all results. A list can store any data (vector, list, data frame, matrix, etc.)
  # This is to make sure that the other analysis uses results only from this function
  listPreprocessed$Metadata$FunctionOrigin <- "performPreprocessingPeakData"
  df_dimensions               <- data.frame()
  listPreprocessed$Dimensions <- rbind(df_dimensions, c("step", "n_samples", "n_features")) %>% `colnames<-`(NULL) # Update list

  # Transpose the data
  data_transposed <- as.data.frame(t(raw_data)) %>%
    setNames(as.character(.[1, ])) %>% .[-1, ] # Set column names from 1st row (features, etc.)

  # Arrange by Injection Sequence
  data_transposed <- data_transposed %>%
    dplyr::mutate(Injection = as.numeric(Injection)) %>% # Convert to numeric so it will be arranged correctly in the line below
    dplyr::arrange(Injection) # Sort by injection, needed in drift-correction

  # Data checking steps. Stop if these columns to not exist in the data
  required_headers <- c("Group", "Group2", "Sample", "Batch", "Injection", "Osmolality")
  if (!all(required_headers %in% colnames(data_transposed))) {
    stop("⚠️ Data preprocessing is halted. Data does not contain the expected column names (Group, Group2, Sample, Batch, Injection, Osmolality). Check for typo or the 1st 6 rows in your csv file.")
  }

  # Remove user-defined outliers
  if (!is.null(outliers)) {

    message("outliers is not empty, removing biological and/or QC samples.")
    flush.console() # Force message to be displayed

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

    # message("outliers is not empty, removing biological and/or QC samples.")
    flush.console() # Force message to be displayed

  } else {
    removed_outliers <- data.frame() # Quickest way to solve as of the moment
  }

  # List of metadata indices
  indices_sqc     <- data_transposed$Group2 == "SQC"                # Indices for Sample QC
  indices_eqc     <- data_transposed$Group2 == "EQC"                # Indices for Extract QC
  indices_qc      <- data_transposed$Group2 %in% c("SQC", "EQC")    # Indices for both QCs
  indices_non_qc  <- !indices_qc                                   # Indices for Biological Samples, i.e., non-QC

  # Saving metadata

  # Compile all metadata in 1 list
  listPreprocessed$Metadata$Groups            <- data_transposed$Group
  listPreprocessed$Metadata$Groups2           <- data_transposed$Group2
  listPreprocessed$Metadata$Samples           <- data_transposed$Sample
  listPreprocessed$Metadata$Batches           <- data_transposed$Batch %>% as.numeric()
  listPreprocessed$Metadata$InjectionSequence <- data_transposed$Injection %>% as.numeric()
  listPreprocessed$Metadata$OsmolalityValues  <- data_transposed$Osmolality[indices_non_qc]

  # Extract Samples x Features data only
  df <- data_transposed %>% dplyr::select(-one_of(required_headers)) # remove required_headers columns

  listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                       # Add the removed outliers. I did this to avoid creating multiple variables
                                       c("Original", dim(df)[1] + length(removed_outliers), dim(df)[2])) # Update list
  listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions, c("Removed outliers", dim(df)[1], dim(df)[2])) # Update list

  listPreprocessed$SamplesXFeatures <- df # Update list

  df[df == 0] <- NA # Replace 0 with NA
  num_na      <- sum(is.na(df), na.rm = TRUE) # Count NAs
  df          <- `rownames<-`(df, data_transposed$Sample) # set row names to sample IDs
  df[]        <- lapply(df, as.numeric) # Convert to numeric
  df          <- df[, colSums(is.na(df)) < nrow(df)] # Remove features if all NAs
  df          <- df[, colMeans(is.na(df)) <= filterMissing/100] # Remove columns with more than a specified percentage of missing values

  listPreprocessed$Dimensions        <- rbind(listPreprocessed$Dimensions, c("Missing Percent Greater", dim(df)[1], dim(df)[2])) # Update list
  listPreprocessed$MissingValues     <- num_na

  df[]        <- lapply(df, function(x) { # Replace missing values with a fraction of the minimum positive value
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
  flush.console() # Force message to be displayed

  message("Employing: Drift correction then batch correction. This only outputs 1 data frame.")
  flush.console()

  if (isTRUE(driftBatchCorrection)) {


    df <- pmp::QCRSC(
      df       = df,
      order    = listPreprocessed$Metadata$InjectionSequence,
      batch    = listPreprocessed$Metadata$Batches,
      # Replace SQC and EQC with QC and back to vector
      classes  = data.frame(Groups = listPreprocessed$Metadata$Groups) %>% dplyr::mutate(Groups = ifelse(Groups %in% c("SQC", "EQC"), "QC", Groups)) %>% as.vector() %>% .$Groups,
      spar     = 0, # 0 Spline smoothing parameter. Should be in the range 0 to 1. If set to 0 it will be estimated using leave-one-out cross-validation.
      log      = TRUE,
      minQC    = 5,
      qc_label = "QC",
      spar_lim = c(-1.5, 1.5)
    ) %>% t() %>% as.data.frame()

    listPreprocessed$data_driftBatchCorrected <- df # Update list

    message("✅ Employed: Drift correction then batch correction with Spline smoothing parameter (leave-one-out cross-validation).")
    flush.console()
  } else { # if driftBatchCorrection is not TRUE
    df <- df
    listPreprocessed$data_driftBatchCorrected <- df # Update list

    message("❌ NOT Employed: Drift correction then batch correction with Spline smoothing parameter (leave-one-out cross-validation).")
    flush.console()

  }

  message("Employing: 2nd Replacing of missing values using 1/denMissing, because drift correction and batch correction likely produced NAs.")
  flush.console()

  # Replace missing values 2, because drift correction and batch correction produced NAs
  df[]        <- lapply(df, function(x) { # Replace missing values with a fraction of the minimum positive value
    min_val   <- base::min(x[x > 0], na.rm = TRUE)
    base::replace(
      x      = x,                     # vector
      list   = is.na(x),              # index vector
      values = min_val / denMissing   # replacement values ## Replace with minimum value only since only some values were NAs
    )
  })

  listPreprocessed$data_driftBatchCorrected_noNa <- df # Update list

  message("✅ Employed: 2nd Replacing of missing values introduced by drift- and batch-correction.")
  flush.console()


  # Plot before and after drift- and batch-correction using 6 random features
  beforeAfterBatchDrift <- pmp::sbc_plot(
    df           = listPreprocessed$data_no_NA,
    corrected_df = df,
    classes      = listPreprocessed$Metadata$Groups,
    batch        = listPreprocessed$Metadata$Batches,
    indexes      = sample(x = 1:dim(df)[2], size = 6, replace = FALSE), # select up 'size' number of plots, choices are from 1 to the number of features in the latest df
    qc_label     = "QC",
    output       = NULL
  )

  plot_beforeAfterBatchDrift <- gridExtra::grid.arrange(
    beforeAfterBatchDrift[[1]],
    beforeAfterBatchDrift[[2]],
    beforeAfterBatchDrift[[3]],
    beforeAfterBatchDrift[[4]],
    beforeAfterBatchDrift[[5]],
    beforeAfterBatchDrift[[6]]
  )

  listPreprocessed$plot_beforeAfterDriftBatchCorrection <- plot_beforeAfterBatchDrift

  message("Employing: Filtering of uninformative features using RSD (Relative Standard Deviation).")
  flush.console()

  # Filter uninformative features using RSD (Relative Standard Deviation)
  if (!is.null(filterMaxRSD)) {
    df <- pmp::filter_peaks_by_rsd(
      df, max_rsd = filterMaxRSD, # 20 = 20% for LC-MS, 30 = 30% for GC-MS
      # class = data_transposed$Group2, qc_label = filterMaxRSD_by,

      class    = if (filterMaxRSD_by == "both") {data_transposed$Group} else {data_transposed$Group2},
      qc_label = ifelse(filterMaxRSD_by == "both", "QC", filterMaxRSD_by)

    ) %>% t() %>% as.data.frame()

    listPreprocessed$data_filteredMaxRSD <- df # Update list

    listPreprocessed$Dimensions        <- rbind(listPreprocessed$Dimensions, c("Relative SD", dim(df)[1], dim(df)[2])) %>% `colnames<-`(NULL) # Update list

    message(paste0("✅ Employed: Filtering of uninformative features using RSD (Relative Standard Deviation) = ", filterMaxRSD, "."))
    flush.console()
  } else {
    listPreprocessed$data_filteredMaxRSD <- df # Update list

    message(paste0("❌  NOT Employed: Filtering of uninformative features using RSD (Relative Standard Deviation) = ", filterMaxRSD, "."))
    flush.console()
  }

  message("Employing: Filtering features with low variability across samples.")
  flush.console()

  # Filter features with low variability across samples
  if (!is.null(filterMaxVarSD)) {
    df <- df %>%
      { .[, apply(., 2, sd) > quantile(apply(., 2, sd), filterMaxVarSD/100)] } # 0.10 = 10% standard deviation

    listPreprocessed$data_filteredMaxVarSD <- df # Update list

    listPreprocessed$Dimensions        <- rbind(listPreprocessed$Dimensions, c("Low Variance", dim(df)[1], dim(df)[2])) %>% `colnames<-`(NULL) # Update list

    message(paste0("✅ Employed: Filtered uninformative features with low variability across samples. Used threshold = ", filterMaxVarSD, "th percentile."))
    flush.console()
  } else {
    listPreprocessed$data_filteredMaxVarSD <- df # Update list

    message("❌  NOT Employed: Filtering of uninformative features with low variability across samples.")
    flush.console()
  }

  message("Employing: Data Normalization.")
  flush.console()

  # Apply normalization based on the specified method
  if (dataNormalize == "none") {
    message("❌ NOT Employed: Data Normalization.")
    # return(df)

  } else if (dataNormalize == "SpecificGravity") {
    sg_values <- as.numeric(data_transposed$Osmolality)
    if (all(is.na(sg_values))) {
      message("⚠️ SpecificGravity values are all NA. Normalization by Sum will be performed instead.")
      sample_sums <- rowSums(df, na.rm = TRUE)
      df <- df / sample_sums
    } else {
      # df <- df / sg_values
      df[indices_non_qc, ] <- sweep(df[indices_non_qc, ], 1, sg_values, "/") # normalize for those with osmolality values
      df[indices_qc, ]     <- sweep(df[indices_qc, ]    , 1, rowSums(df(df[indices_qc, ]), na.rm = TRUE)) # normalize by sum the QCs

      message("✅ Employed: Normalization using SpecificGravity values.")
    }
    # return(df)

  } else if (dataNormalize == "sum") {
    sample_sums <- rowSums(df, na.rm = TRUE)
    df <- df / sample_sums
    message("✅ Employed: Normalization by Sum.")
    # return(df)

  } else if (dataNormalize == "median") {
    sample_medians <- apply(df, 1, median, na.rm = TRUE)
    df <- df / sample_medians
    message("✅ Employed: Normalization by Median.")
    # return(df)

  } else if (dataNormalize == "PQN1") {

    reference_spectrum <- apply(df, 2, median, na.rm = TRUE)        # Calculate reference spectrum as median across samples
    quotients          <- df / reference_spectrum                   # Compute quotients
    median_quotients   <- apply(quotients, 1, median, na.rm = TRUE) # Calculate median quotient for each sample
    df                 <- df / median_quotients                     # Normalize samples

    message("✅ Employed: Probabilistic Quotient Normalization (PQN) according to global median approach.")
    # return(df)

  } else if (dataNormalize == "PQN2") {
    reference_spectrum <- df[refSample, ] %>% t()                   # Extract the reference sample's metabolite intensities
    quotients          <- sweep(df, 2, reference_spectrum, "/")     # Compute quotients by dividing each sample by the reference spectrum
    median_quotients   <- apply(quotients, 1, median, na.rm = TRUE) # Calculate median quotients for each sample
    df_normalized      <- sweep(df, 1, median_quotients, "/")       # Normalize by dividing each sample by its corresponding median quotient

    message(paste0("✅ Employed:  Probabilistic Quotient Normalization (PQN) using '", refSample, "' as the reference sample in the normalization."))
    # return(df)

  } else if (dataNormalize == "groupPQN") {
    if (is.null(groupSample)) {
      stop("Please specify 'groupSample' for groupPQN normalization.")
    }
    # Identify pooled group samples
    # pooled_indices <- data_transposed$Group == groupSample
    if (groupSample == "both") {
      pooled_indices <- data_transposed$Group == "QC"
    } else {
      pooled_indices <- data_transposed$Group2 == groupSample
    }

    if (sum(pooled_indices) == 0) {
      stop("No samples found for the specified 'groupSample'.")
    }
    # Calculate reference spectrum as median of pooled group
    reference_spectrum <- apply(df[pooled_indices, , drop = FALSE], 2, median, na.rm = TRUE)
    # Compute quotients
    quotients <- df / reference_spectrum
    # Calculate median quotient for each sample
    median_quotients <- apply(quotients, 1, median, na.rm = TRUE)
    # Normalize samples
    df <- df / median_quotients
    message("✅ Employed: Group Probabilistic Quotient Normalization (groupPQN) using pooled group: ", groupSample)

  } else if (dataNormalize == "quantile") {

    df <- pmp::pqn_normalisation(
      df          = df,
      classes     = listPreprocessed$Metadata$Groups,
      qc_label    = "QC",
      ref_mean    = NULL,
      qc_frac     = 0,
      sample_frac = 0,
      ref_method  = "mean"
    ) %>% t() %>% as.data.frame()



    message("✅ Employed: Quantile Normalization.")

  }

  listPreprocessed$data_normalized <- df # Update list

  message("Employing: Data transformation.")
  flush.console()

  # Assume `dataTransform` is one of: "none", "log10", "sqrt", "cbrt"
  if (dataTransform == "log10") {
    df <- df - min(df) + 1  # Ensure all values are positive for log10
    df <- log10(df) %>% as.data.frame()
    message("✅ Employed: Log10 transformation.")

  } else if (dataTransform == "sqrt") {
    df <- df - min(df) + 1  # Ensure all values are non-negative for sqrt
    df <- sqrt(df) %>% as.data.frame()
    message("✅ Employed: Square root transformation.")

  } else if (dataTransform == "cbrt") {
    # df <- sign(df) * abs(df)^(1/3) %>% as.data.frame()  # Cube root handles negative values
    df <- df - min(df) + 1  # Ensure all values are non-negative
    df <- df^(1/3) %>% as.data.frame()
    message("✅ Employed: Cube root transformation.")

  } else {
    message("❌  NOT Employed: Any transformation ('none' was selected).")
  }

  flush.console()
  listPreprocessed$data_transformed <- df # Update list

  # Data Scaling

  # Scaling function
  scale_data <- function(data, method) {
    if (is.null(method)) {
      return(data)
    } else if (method == "mean") {
      return(scale(data, center = TRUE, scale = FALSE) %>% as.data.frame())
    } else if (method == "meanSD") {
      return(scale(data, center = TRUE, scale = TRUE) %>% as.data.frame())
    } else if (method == "meanSD2") {
      return(scale(data, center = TRUE, scale = sqrt(apply(data, 2, sd))) %>% as.data.frame())
    }
  }

  message("Employing: Data scaling for PCA.")
  flush.console()

  # Apply scaling for PCA
  scaled_pca_data <- scale_data(data = df, method = dataScalePCA) #%>% t() %>% as.data.frame()

  # Update list
  if (is.null(dataScalePCA)) {
    listPreprocessed$data_scaledPCA <- df
  } else {
    listPreprocessed$data_scaledPCA <- scaled_pca_data
  }

  message(paste0(ifelse(is.null(dataScalePCA), "❌  NOT Employed: ", "✅ Employed: "), "Data scaling for PCA."))
  flush.console()

  message("Employing: Data scaling for OPLS-DA.")
  flush.console()

  # Apply scaling for OPLS-DA
  scaled_oplsda_data <- scale_data(data = df, method = dataScaleOPLSDA) #%>% t() %>% as.data.frame()

  # Update list
  if (is.null(dataScaleOPLSDA)) {
    listPreprocessed$data_scaledOPLSDA <- df
  } else {
    listPreprocessed$data_scaledOPLSDA <- scaled_oplsda_data
  }

  message(paste0(ifelse(is.null(dataScaleOPLSDA), "❌  NOT Employed: ", "✅ Employed: "), "Data scaling for OPLS-DA."))
  flush.console()

  message("✅✅✅ Data preprocessing is complete. The data is now ready for any statistical analysis based on selected data preprocessing techniques.")
  flush.console()

  return(listPreprocessed)
}
