performPreprocessingPeakData <- function(
    raw_data,
    filterMissing   = 20, # Minimum %missing in all sample groups required to remove feature
    denMissing      = 5, # Missing value imputation. A denominator in 1/denMissing
    # batchAlign      = FALSE, # Perform multi-batch alignment to merge features artificially split between batches
    # tolerance_mz    = 0.001, # Tolerance for difference in m/z
    # tolerance_rt    = 15, # Tolerance for difference in retention time
    # driftCorrection = TRUE, # Perform drift correction
    # modelName       = c("VVV", "VVE", "VEV", "VEE", "VEI", "VVI", "VII"), # After run 1, users have the option to re-run again using the results of clustering. Look for "MClust final model with 5 clusters and VVV geometry." and take the 3 letters before the word "geometry."
    # nClusters       = seq(5, 35, by = 10), # After run 1, users have the option to re-run again using the results of clustering. Look for the number of clusters.
    # batchCorrection = TRUE, # c(0, 1, 2) 0 = Do not perform batch correction; 1 = Perform batch correction for 1 batch, 2 = for more 2 or more batches

    driftBatchCorrection = TRUE, # Logical. Defaults to TRUE to perform between-batch correction.

    filterMaxRSD    = 30, # NULL to skip this filtering. LC-MS = 20 for 20%; GC-MS = 30 for 30%
    filterMaxRSD_by = c("SQC", "EQC", "both")[2], # c("SQC", "EQC", "both"). If "both" then all QCs will be used. Defaults to "EQC"


    filterMaxVarSD  = 10, # NULL to skip this filtering. Remove nth percentile of features with the lowest variability

    # choices are TRUE (automatic checking of osmolality values, if none, then by sum)
    #             FALSE (no normalization)
    #             SUM, MEDIAN  (normalize by sum or median)
    #             PQN, groupPQN (by reference sample, pooled sample from group)
    #             QUANTILE (by quantile, for > 1,000 features)
    dataNormalize       = c("none", "SpecificGravity", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile")[2],
    refSample           = NULL, # Provide the reference sample in the case of normalization method "PQN2" Must be in the samples, and is not part of "outliers" vector
    # c("SQC", "EQC", "both")
    groupSample         = NULL, #ifelse(dataNormalize == "groupPQN", "EQC", NULL), # Used in the groupPQN. NULL if not "groupPQN" Default to "EQC" if "groupPQN" Other choice is "SQC"

    dataTransform       = c("none", "log10", "sqrt", "cbrt")[2], # Transform the data
    dataScalePCA        = c("none", "mean", "meanSD", "meanSD2")[3], # Defaults to "meanSD". c("none", "mean", "meanSD", "meanSD2"). "mean" = mean-centered only; meanSD = mean-centered and divided by SD of each feature; meanSD2 = mean-centered and divided by the square root of SD of each feature
    dataScaleOPLSDA     = c("none", "mean", "meanSD", "meanSD2")[4], # Same choices as in dataScalePCA. Defaults to "meanSD2".
    outliers            = NULL # A vector of biological samples and/or QC that are considered as outliers
) {


  # Parameter Validation Function
  check_parameters <- function(filterMissing,
                               denMissing,
                               # batchAlign,
                               # tolerance_mz,
                               # tolerance_rt,
                               # driftCorrection,
                               # modelName,
                               # nClusters,
                               # batchCorrection,

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

    # if (!is.numeric(tolerance_mz)) {
    #   errors <- c(errors, "tolerance_mz: Not numeric.")
    # } else if (tolerance_mz < 0) {
    #   errors <- c(errors, "tolerance_mz: Is a negative number.")
    # }
    #
    # if (!is.numeric(tolerance_rt)) {
    #   errors <- c(errors, "tolerance_rt: Not numeric.")
    # } else if (tolerance_rt < 0) {
    #   errors <- c(errors, "tolerance_rt: Is a negative number.")
    # }

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

    # if (!is.logical(batchAlign)) {
    #   errors <- c(errors, "batchAlign: Not logical. Must be either TRUE or FALSE.")
    # }
    #
    # if (!is.logical(driftCorrection)) {
    #   errors <- c(errors, "driftCorrection: Not logical. Must be either TRUE or FALSE.")
    # }
    #
    # if (!is.logical(batchCorrection)) {
    #   errors <- c(errors, "batchCorrection: Not logical. Must be either TRUE or FALSE.")
    # }


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

    # allowed_scales <- c("VVV", "VVE", "VEV", "VEE", "VEI", "VVI", "VII")
    # if (!is.null(modelName)) {
    #   if (!any(modelName %in% allowed_scales)) {
    #     errors <- c(errors, "modelName: Must be either VVV, VVE, VEV, VEE, VEI, VVI, VII (see mclust documentation).")
    #   }
    # }

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
                                        # batchAlign,
                                        # tolerance_mz,
                                        # tolerance_rt,
                                        # driftCorrection,
                                        # modelName,
                                        # nClusters,
                                        # batchCorrection,

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



  # listPreprocessed$data_transposed   <- data_transposed # Update list

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

  # df_with_NA  <- df # For drift correction purposes, specifically in the alignBatches function, see ?alignBatches

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


  # df_not_yet_driftBatchCorrected_withNA <- df


  {
    # message("Employing: Multi-batch alignment to merge features artificially split between batches.")
    # flush.console() # Force message to be displayed
    #
    #
    # if (batchAlign == TRUE) {
    #   if (length(unique(listPreprocessed$Metadata$Batches)) == 1) {
    #
    #     listPreprocessed$data_batchAligned <- df
    #
    #     message("NOT Employed: Multi-batch alignment to merge features artificially split between batches.")
    #     message(paste0("Reason: There is only 1 batch, i.e., Batch ", unique(listPreprocessed$Metadata$Batches), "."))
    #     flush.console() # Force message to be displayed
    #
    #
    #   } else { #if (length(unique(listPreprocessed$Metadata$Batches) == 1)) {
    #
    #     df <- batchCorr::alignBatches(
    #       peakInfo      = batchCorr::peakInfo(df, sep = "@", 2, 0),
    #       PeakTabNoFill = df_with_NA,
    #       PeakTabFilled = df,
    #       batches       = listPreprocessed$Metadata$Batches,
    #       sampleGroups  = listPreprocessed$Metadata$Groups,
    #       selectGroup   = "QC",
    #       mzdiff        = tolerance_mz,
    #       rtdiff        = tolerance_rt,
    #       report        = FALSE # export diagnostic plots into your work directory
    #     )
    #
    #     listPreprocessed$data_batchAligned <- df
    #
    #
    #     message("Employed: Multi-batch alignment to merge features artificially split between batches.")
    #     flush.console() # Force message to be displayed
    #
    #
    #   }
    # } else if (batchAlign == FALSE) {
    #   message("NOT Employed: batchAlign = FALSE. Multi-batch alignment to merge features artificially split between batches.")
    #   flush.console() # Force message to be displayed
    # }
    #
    #
    #
    #
    #
    # message("Employing: Drift correction.")
    # flush.console()
    #
    # if (driftCorrection == TRUE) {
    #
    # # Drift Correction
    # df <- batchCorr::correctDrift(
    #   peakTable    = listPreprocessed$data_no_NA,
    #   injections   = listPreprocessed$Metadata$InjectionSequence,
    #   sampleGroups = listPreprocessed$Metadata$Groups,
    #   QCID         = "QC",
    #   modelNames   = modelName, # c('VVV', 'VVE', 'VEV', 'VEE', 'VEI', 'VVI', 'VII'),
    #   G            = nClusters, # seq(5, 35, by = 10),
    #   report       = FALSE # TRUE if to print pdf reports of drift models
    # )
    #
    # listPreprocessed$data_allResultsDriftCorrection <- df # Update list
    #
    # df <- as.data.frame(df$TestFeatsCorr) # Extract drift-corrected data
    #
    # listPreprocessed$data_driftCorrected <- df # Update list
    #
    # message("Employed: Drift correction.")
    # flush.console()
    #
    # } else if (driftCorrection == FALSE) {
    #   listPreprocessed$data_allResultsDriftCorrection <- df # Update list
    #
    #   message("NOT Employed: Drift correction. driftCorrection = FALSE")
    #   flush.console()
    # } else {
    #   message("NOT Employed: Drift correction. driftCorrection must be either TRUE or FALSE.")
    #   flush.console()
    # }
    #
    #
    #
    #
    #
    #
    #
    #
    # # message("Employing: Batch correction.")
    # # flush.console()
    # #
    # # # Batch correction using sva
    # # if (length(unique(data_transposed$Batch)) > 1) {
    # #   if (batchCorrection == FALSE) {
    # #     message("NOT Employed: Batch correction.")
    # #     flush.console()
    # #     # df remains unchanged
    # #   } else if (batchCorrection == TRUE) {
    # #     df <- t(df) %>%
    # #       sva::ComBat(batch = as.factor(data_transposed$Batch),
    # #                   mod   = model.matrix(~ as.factor(data_transposed$Group))) %>%
    # #       t() %>% `rownames<-`(data_transposed$Sample)
    # #
    # #     listPreprocessed$data_batchCorrected <- df # Update list
    # #
    # #     message("Employed: Batch correction has been completed using ComBat.")
    # #     flush.console()
    # #   }
    # # } else {
    # #   listPreprocessed$data_batchCorrected <- df # Update list
    # #
    # #   message(paste0("NOT Employed: Batch correction since there is only 1 batch (Batch ", unique(data_transposed$Batch) ,")."))
    # #   flush.console()
    # # }
    #
    #
    #
    # message("Employing: Batch correction.")
    # flush.console()
    #
    # # Batch correction using sva
    # if (batchCorrection == TRUE) {
    #   if (length(unique(data_transposed$Batch)) == 1) {
    #     message(paste0("NOT Employed: Batch correction since there is only 1 batch (Batch ", unique(data_transposed$Batch) ,")."))
    #     flush.console()
    #   } else if (length(unique(data_transposed$Batch)) >= 1) {
    #     df <- t(df) %>%
    #       sva::ComBat(batch = as.factor(data_transposed$Batch),
    #                   mod   = model.matrix(~ as.factor(data_transposed$Group))) %>%
    #       t() %>% `rownames<-`(data_transposed$Sample) %>% as.data.frame()
    #
    #     listPreprocessed$data_batchCorrected <- df # Update list
    #
    #     message("Employed: Batch correction has been completed using ComBat.")
    #     flush.console()
    #   } else {
    #     message("NOT Employed: Batch correction. batchCorrection might not be TRUE/FALSE or there are no batches.")
    #     flush.console()
    #   }
    #
    # } else if (batchCorrection == FALSE) {
    #   message("NOT Employed: Batch correction.")
    #   flush.console()
    #   # df remains unchanged
    # } else {
    #   message("NOT Employed: Batch correction. batchCorrection might not be TRUE/FALSE or there are no batches.")
    #   flush.console()
    # }
  }

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

  # message("Employing: Normalization.")
  # flush.console()
  # # Normalize using osmolality values if present, otherwise by sum. No normalization is set to FALSE
  # if (dataNormalize == TRUE) {
  #   osmolality_values <- as.numeric(data_transposed$Osmolality)
  #   non_qc_indices <- data_transposed$Group != "QC"
  #   osmolality_values_non_qc <- osmolality_values[non_qc_indices]
  #
  #   num_na <- sum(is.na(osmolality_values_non_qc))
  #   total_samples <- length(osmolality_values_non_qc)
  #
  #   if (num_na == total_samples) { # Checks if all osmolality values are NAs then normalize by sum
  #     sample_sums <- rowSums(df, na.rm = TRUE)
  #     df <- sweep(df, 1, sample_sums, "/")
  #
  #     message("✅ Employed: All osmolality values are NA. Normalization by Sum was performed instead.")
  #     flush.console()
  #
  #     listPreprocessed$data_normalized <- df # Update list
  #   } else if (num_na > 0) { # Checks if any osmolality values are NAs, then normalize by sum
  #     sample_sums <- rowSums(df, na.rm = TRUE)
  #     df <- sweep(df, 1, sample_sums, "/")
  #
  #     message("✅ Employed: At least one osmolality value is NA. Normalization by Sum was performed instead.")
  #     flush.console()
  #
  #     listPreprocessed$data_normalized <- df # Update list
  #   } else { # No NA values, normalize by osmolality
  #     df[non_qc_indices, ] <- sweep(df[non_qc_indices, ], 1, osmolality_values_non_qc, "/")
  #
  #     qc_sums <- rowSums(df[!non_qc_indices, ], na.rm = TRUE) # Filter QC rows
  #     df[!non_qc_indices, ] <- sweep(df[!non_qc_indices, ], 1, qc_sums, "/") # Normalize by sum the QC rows
  #
  #
  #     message("✅ Employed: Normalization using osmolality values was performed to biological samples, while normalization by Sum was perform to QC samples.")
  #     flush.console()
  #
  #     listPreprocessed$data_normalized <- df # Update list
  #   }
  # } else { # When dataNormalize = FALSE
  #   listPreprocessed$data_normalized <- df # Update list
  #
  #   message("❌  NOT Employed: Data Normalization.")
  #   flush.console()
  # }



  #   message("Employing: Normalization.")
  #   flush.console()
  #
  #   # Normalize data based on the specified method
  #   if (dataNormalize == TRUE) {
  #     # Ensure 'dataNormalize' is defined; default to 'sum' if not
  #     dataNormalize <- tolower(dataNormalize)
  #
  #     # Extract osmolality values and identify non-QC samples
  #     osmolality_values        <- as.numeric(data_transposed$Osmolality)
  #     non_qc_indices           <- data_transposed$Group != "QC"
  #     osmolality_values_non_qc <- osmolality_values[non_qc_indices]
  #
  #     # Count NA values in osmolality
  #     num_na        <- sum(is.na(osmolality_values_non_qc))
  #     total_samples <- length(osmolality_values_non_qc)
  #
  #     # Define normalization functions
  #     normalize_by_sum              <- function(x) x / sum(x, na.rm = TRUE)
  #     normalize_by_median           <- function(x) x / median(x, na.rm = TRUE)
  #     normalize_by_reference_sample <- function(df, refSample) {
  #       ref_values <- as.numeric(df[refSample, ])
  #       t(apply(df, 1, function(x) x / median(x / ref_values, na.rm = TRUE)))
  #     }
  #     normalize_by_pooled_group     <- function(df, group_vector, group_name) {
  #       pooled_sample <- colMeans(df[group_vector == group_name, , drop = FALSE], na.rm = TRUE)
  #       t(apply(df, 1, function(x) x / median(x / pooled_sample, na.rm = TRUE)))
  #     }
  #     quantile_normalization        <- function(df) {
  #       if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  #         stop("The 'preprocessCore' package is required for quantile normalization. Please install it.")
  #       }
  #       preprocessCore::normalize.quantiles(as.matrix(df))
  #     }
  #
  #     # Apply the selected normalization method
  #     if (dataNormalize == "osmolality") {
  #       if (num_na == total_samples) {
  #         df <- t(apply(df, 1, normalize_by_sum))
  #         message("✅ Employed: All osmolality values are NA. Normalization by Sum was performed instead.")
  #       } else if (num_na > 0) {
  #         df <- t(apply(df, 1, normalize_by_sum))
  #         message("✅ Employed: At least one osmolality value is NA. Normalization by Sum was performed instead.")
  #       } else {
  #         df[non_qc_indices, ] <- df[non_qc_indices, ] / osmolality_values_non_qc
  #         qc_indices           <- !non_qc_indices
  #         df[qc_indices, ]     <- t(apply(df[qc_indices, ], 1, normalize_by_sum))
  #         message("✅ Employed: Normalization using osmolality values for biological samples; normalization by Sum for QC samples.")
  #       }
  #     } else if (dataNormalize == "sum") {
  #       df <- t(apply(df, 1, normalize_by_sum))
  #       message("✅ Employed: Normalization by Sum.")
  #     } else if (dataNormalize == "median") {
  #       df <- t(apply(df, 1, normalize_by_median))
  #       message("✅ Employed: Normalization by Median.")
  #     } else if (dataNormalize == "reference_sample") {
  #       # Ensure 'reference_sample' is defined
  #       if (!exists("reference_sample") || !(reference_sample %in% rownames(df))) {
  #         stop("Please specify a valid 'reference_sample' present in the data.")
  #       }
  #       df <- normalize_by_reference_sample(df, reference_sample)
  #       message("✅ Employed: Normalization by Reference Sample (PQN).")
  #     } else if (dataNormalize == "pooled_group") {
  #       # Ensure 'pooled_group_name' is defined
  #       if (!exists("pooled_group_name") || !(pooled_group_name %in% data_transposed$Group)) {
  #         stop("Please specify a valid 'pooled_group_name' present in the Group column.")
  #       }
  #       group_vector <- data_transposed$Group
  #       df           <- normalize_by_pooled_group(df, group_vector, pooled_group_name)
  #       message("✅ Employed: Normalization by Pooled Sample Group (Group PQN).")
  #     } else if (dataNormalize == "quantile") {
  #       df           <- quantile_normalization(df)
  #       message("✅ Employed: Quantile Normalization.")
  #     } else {
  #       stop("Invalid normalization method specified.")
  #     }
  #
  #     # Update the preprocessed data
  #     listPreprocessed$data_normalized <- df
  #   } else {
  #     listPreprocessed$data_normalized <- df
  #     message("❌ NOT Employed: Data Normalization.")
  #
  #
  #
  #   }
  #
  #
  #
  # }













  message("Employing: Data Normalization.")
  flush.console()


  # Define normalization function
  # normalize_data <- function(df, data_transposed, dataNormalize = "none", pooled_group_name = NULL) {
  # Ensure df is a matrix for efficient computation
  # df <- as.matrix(df)

  # Identify non-QC samples
  # non_qc_indices <- data_transposed$Group != "QC"

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


      # df[non_qc_indices, ] <- sweep(df[non_qc_indices, ], 1, osmolality_values_non_qc, "/")
      #
      # qc_sums <- rowSums(df[!non_qc_indices, ], na.rm = TRUE) # Filter QC rows
      # df[!non_qc_indices, ] <- sweep(df[!non_qc_indices, ], 1, qc_sums, "/") # Normalize by sum the QC rows


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
    # return(df)

  } else if (dataNormalize == "quantile") {
    # Quantile normalization
    # quantile_normalize <- function(mat) {
    #   mat_rank <- apply(mat, 2, rank, ties.method = "min")
    #   mat_sorted <- apply(mat, 2, sort)
    #   row_means <- rowMeans(mat_sorted)
    #   index_to_mean <- function(ranks) row_means[ranks]
    #   mat_normalized <- apply(mat_rank, 2, index_to_mean)
    #   return(mat_normalized)
    # }
    # df <- quantile_normalize(df)


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
    # return(df)

  }
  #   else {
  #     stop("Invalid normalization method specified.")
  #   }
  # }

  listPreprocessed$data_normalized <- df # Update list




















  message("Employing: Data transformation.")
  flush.console()
  # # Log10 transformation
  # if (dataTransform == TRUE) {
  #   # Ensure the minimum value is positive before applying log10
  #   df <- df - min(df) + 1 # Shifts the data to make all values positive
  #
  #   log_transformed_data <- log10(df) %>% as.data.frame()
  #
  #   message("✅ Employed: Log10 transformation.")
  #   flush.console()
  #
  #   df <- log_transformed_data
  #
  #   listPreprocessed$data_log10transformed <- df # Update list
  # } else { # When dataTransform = FALSE
  #   listPreprocessed$data_log10transformed <- df # Update list
  #
  #   message("❌  NOT Employed: Log10 transformation.")
  #   flush.console()
  # }


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


## Usage
# mydata <- performPreprocessingPeakData(
#   raw_data = org_data,
#   filterMissing   = 20, # Minimum %missing in all sample groups required to remove feature
#   denMissing      = 5, # Missing value imputation. A denominator in 1/denMissing
#   batchCorrection = 0, # c(0, 1, 2) 0 = Do not perform batch correction; 1 = Perform batch correction for 1 batch, 2 = for more 2 or more batches
#   filterMaxRSD    = NULL, # LC-MS = 20 for 20%; GC-MS = 30 for 30%; NULL to skip this filtering
#   filterMaxVarSD  = NULL, # Remove nth percentile of features with the lowest variability; NULL to skip this filtering
#   dataNormalize       = FALSE, # Logical. FALSE to not normalize. Defaults to TRUE to normalize data using osmolality values if given, otherwise normalizes by sum
#   dataTransform  = FALSE, # Logical. Defaults to TRUE to log10 transform the data
#   dataScalePCA        = NULL, # Defaults to "meanSD". c(NULL, "mean", "meanSD", "meanSD2"). "mean" = mean-centered only; meanSD = mean-centered and divided by SD of each feature; meanSD2 = mean-centered and divided by the square root of SD of each feature
#   dataScaleOPLSDA     = NULL # Same choices as in dataScalePCA. Defaults to "meanSD2".
# )
