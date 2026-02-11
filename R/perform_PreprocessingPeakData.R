#' Comprehensive Data Preprocessing Pipeline
#'
#' @description
#' Performs a complete data preprocessing workflow. This function
#' is designed to work seamlessly with the output of `perform_DataQualityCheck`,
#' directly using the pre-processed `raw_data`. It applies
#' vectorized operations for speed and includes a robust,
#' technical replicate merging step.
#'
#' @param raw_data List. Quality-checked data *list* from the `perform_DataQualityCheck` function.
#' @param outliers Vector. Biological samples and/or QC samples considered as outliers.
#' @param filterMissing Numeric. Minimum percentage of missing values (0-100).
#' @param filterMissing_by_group Boolean. Assess missingness group-wise.
#' @param filterMissing_includeQC Boolean. Include QCs in missingness filter.
#' @param denMissing Numeric. Denominator for missing value imputation (e.g., 5 -> 1/5th of min).
#' @param driftBatchCorrection Boolean. If `TRUE`, perform QC-RSC.
#' @param spline_smooth_param Numeric. Spline smoothing parameter (0-1).
#' @param spline_smooth_param_limit Vector. `c(min, max)` for spline limits.
#' @param log_scale Boolean. Fit correction on log-scaled data.
#' @param min_QC Numeric. Minimum QCs required per batch for correction.
#' @param removeUncorrectedFeatures Boolean. If `TRUE`, remove features QCRSC couldn't correct. Defaults to `FALSE` and we let RSD and low-variance filtering remove uninformative features.
#' @param dataNormalize String. Normalization method.
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
#'    Default: "Normalization" (if present, otherwise, defaults to "sum"). QC samples are normalized by `median` (see https://pmc.ncbi.nlm.nih.gov/articles/PMC6835889/).
#' @param refSample String. Reference sample for `dataNormalize = "PQN2"`.
#' @param groupSample String. Group for `dataNormalize = "groupPQN"`.
#' @param reference_method String. Method (`"mean"`, `"median"`) for quantile ref.
#' @param dataTransform String. Data transformation method. Options:
#'  \itemize{
#'    \item \code{"none"}: No transformation
#'    \item \code{"log2"}: log base 2
#'    \item \code{"log10"}: log base 10
#'    \item \code{"sqrt"}: Square-root
#'    \item \code{"cbrt"}: Cube-root
#'    \item \code{"vsn"}: Variance Stabilizing Normalization (vsn)
#'    \item \code{"glog"}: Variance stabilising generalised logarithm (glog) transformation
#'    }
#' @param dataScaleNONPLS String. Data scaling for non-PLS type analysis (PCA, t-test/ANOVA, fold change, correlation, etc.). Options:
#'  \itemize{
#'    \item \code{"none"}: No data scaling
#'    \item \code{"mean"}: Scale by mean (average)
#'    \item \code{"auto"}: Default. Auto-scaling. Scale by mean divided by standard deviation (SD)
#'    \item \code{"pareto"}: Pareto-scaling. Scale by mean divided by square-root of SD. Suggested for PLS-type analysis (PLS, PLS-DA, OPLS-DA, sPLS-DA)
#'    }
#' @param dataScalePLS String. Data scaling for PLS-type analysis. Same options as in `dataScaleNONPLS`. Defaults to `"pareto"`.
#' @param qcNormalize String. QC normalization method. Options:
#'  \itemize{
#'    \item \code{"none"}: No QC normalization
#'    \item \code{"mean"}: Scale by mean (average) sample `Normalization` values
#'    \item \code{"median"}: Scale by median sample `Normalization` values (default)
#'    }
#' @param filterReference String. Which scaled data to use for filtering so both `NONPLS` and `PLS` scaled data have the same features. Options:
#'  \itemize{
#'    \item \code{"NONPLS"}: Use auto-scaled (or dataScaleNONPLS) features as reference for filtering
#'    \item \code{"PLS"}: Use pareto-scaled (or dataScalePLS) features as reference for filtering
#'    \item \code{"auto"}: Apply filters independently to both datasets and keep only features passing in BOTH (intersection, most robust - default)
#'    }
#' @param filterMaxRSD Numeric. Relative Standard Deviation (RSD) threshold (0-100). Default: 30 (remove features with RSD >= .3)
#' @param filterMaxRSD_by String. QCs to use for RSD. Options:
#'  \itemize{
#'    \item \code{"SQC"}: Filter by sample QC
#'    \item \code{"EQC"}: Filter by extract QC (default)
#'    \item \code{"both"}: Filter by both SQC and EQC (or QC altogether). Use this when there are no SQC and EQC in the "Group" row
#'    }
#' @param filterLowVariability Numeric. Percentile (0-100) for variance filtering. Default: 10 (remove '10th' percentile of features with the lowest variability in the 'Group' row.)
#' @param merge_replicates Logical. If `TRUE`, merge technical replicates by SubjectID.
#' @param num_cores Integer. Number of cores for `dataTransform = "vsn"`. Can be `max` where auto detection of
#' the total number of cores `max` will be performed, and `num_cores <- max - 2`. Defaults to 1.
#' @param verbose Logical. Print progress messages.
#'
#' @returns A list containing all data generated from all preprocessing steps.
#'
#' @details
#'
#' ## Overview
#'
#' This function performs comprehensive data preprocessing for metabolomics data
#' in a fixed, sequential workflow optimized for LC-MS and GC-MS datasets. The
#' preprocessing steps are executed in a predefined order regardless of how
#' parameters are arranged in the function call. This ensures reproducibility
#' and follows established best practices in metabolomics data processing.
#'
#' ## Preprocessing Workflow (Fixed Order)
#'
#' The function applies the following steps sequentially:
#'
#' 1. **Input Validation & Data Preparation**
#'    - Validates input from `perform_DataQualityCheck`
#'    - Extracts metadata and feature data
#'    - Converts matrix to samples × features format
#'
#' 2. **Outlier Removal**
#'    - Removes specified biological and/or QC samples
#'    - Updates sample count accordingly
#'
#' 3. **Missing Value Filtering**
#'    - Removes features exceeding `filterMissing` threshold
#'    - Can be assessed group-wise (`filterMissing_by_group = TRUE`)
#'    - Option to include/exclude QCs in assessment
#'
#' 4. **Missing Value Imputation**
#'    - Replaces zeros and NAs with 1/`denMissing` of column minimum
#'    - Applied to prevent downstream numerical issues
#'
#' 5. **Signal Drift and Batch Correction (QCRSC)**
#'    - **Only applied if `driftBatchCorrection = TRUE`**
#'    - Uses Quality Control-based Robust Spline Correction
#'    - See 'Signal Drift and Batch Correction Details' section below
#'
#' 6. **Normalization**
#'    - Accounts for dilution and sample-to-sample variation
#'    - Multiple methods available (see `dataNormalize` parameter)
#'    - QC samples normalized separately (see `qcNormalize` parameter)
#'
#' 7. **Transformation**
#'    - Stabilizes variance and reduces heteroscedasticity
#'    - Methods include log, sqrt, VSN, glog (see `dataTransform` parameter)
#'
#' 8. **Scaling**
#'    - Two separate scaled datasets generated:
#'      - **NONPLS**: For PCA, t-test, ANOVA, correlation (typically auto-scaling)
#'      - **PLS**: For PLS-DA, OPLS-DA, sPLS-DA (typically Pareto-scaling)
#'
#' 9. **Quality Filtering**
#'    - **RSD Filtering**: Removes features with high QC variability
#'    - **Variance Filtering**: Removes low-variance features
#'    - Applied based on `filterReference` (NONPLS, PLS, or both)
#'
#' 10. **Technical Replicate Merging**
#'     - **Only if `merge_replicates = TRUE`**
#'     - Merges replicates by SubjectID (averages feature intensities)
#'     - QC samples and samples without SubjectID are never merged
#'
#' ## Signal Drift and Batch Correction Details
#'
#' ### What is QCRSC?
#'
#' Quality Control-based Robust Spline Correction (QCRSC) is a widely-used
#' method for correcting systematic signal drift and batch effects in
#' metabolomics data. It was originally developed by Kirwan et al. (2013)
#' and is implemented in the `pmp` Bioconductor package.
#'
#' ### How QCRSC Works
#'
#' 1. **QC Samples as Reference**
#'    - QC samples (pooled quality control samples) are injected regularly
#'      throughout the analytical run
#'    - These should theoretically show constant signal since they represent
#'      the same biological matrix
#'    - Any variation in QC signal indicates systematic drift/batch effects
#'
#' 2. **Feature-wise Correction**
#'    - For each feature (metabolite) independently:
#'      - Extracts QC sample intensities across injection order
#'      - Fits a smoothing spline curve through QC points
#'      - This spline represents the drift/batch effect pattern
#'
#' 3. **Spline Fitting Parameters**
#'    - `spline_smooth_param`: Controls spline smoothness (0 = more flexible)
#'    - `spline_smooth_param_limit`: Constrains correction magnitude
#'    - `log_scale`: Whether to fit spline on log-transformed data
#'    - `min_QC`: Minimum QCs required per batch for reliable fitting
#'
#' 4. **Correction Application**
#'    - The fitted drift curve is used to correct **all samples** (biological + QC)
#'    - Correction formula: `corrected = original / fitted_curve`
#'    - Normalizes all samples to remove systematic drift
#'
#' 5. **Batch Effect Handling**
#'    - If multiple batches exist, QCRSC fits separate splines per batch
#'    - Ensures batch-specific drift patterns are corrected independently
#'
#' ### Uncorrected Features
#'
#' Some features may remain uncorrected if:
#' - **Insufficient QCs**: Fewer than `min_QC` samples in a batch
#' - **Poor spline fit**: Algorithm cannot fit a reliable correction curve
#' - **Numerical issues**: Division by zero or infinite values
#'
#' For uncorrected features, QCRSC returns the **original values unchanged**.
#'
#' This function identifies uncorrected features by comparing pre- and
#' post-correction data using machine precision tolerance. Features with
#' maximum absolute change ≤ `.Machine$double.eps^0.5` are flagged as
#' uncorrected.
#'
#' **Handling Uncorrected Features:**
#' - If `removeUncorrectedFeatures = TRUE`: Remove these features entirely
#' - If `removeUncorrectedFeatures = FALSE`: Flag but retain for downstream analysis
#'
#' The default (`FALSE`) is recommended as these features may still be biologically
#' informative, and subsequent RSD filtering will remove high-variability features
#' regardless of correction status.
#'
#' ### When to Use QCRSC
#'
#' **Use QCRSC when:**
#' - Long analytical runs with visible signal drift
#' - Multiple batches analyzed at different times
#' - QC samples show systematic trends across injection order
#'
#' **Skip QCRSC when:**
#' - Short runs with minimal drift (< 50 samples)
#' - Insufficient QC samples (< 5 per batch)
#' - QC samples already show stable signal
#'
#' ## Normalization Methods
#'
#' - **"none"**: No normalization (not recommended for most datasets)
#' - **"Normalization"**: Uses values from "Normalization" metadata row
#'   (e.g., specific gravity, osmolality for urine; total protein for plasma)
#' - **"sum"**: Total sum normalization (assumes similar total abundance)
#' - **"median"**: Median normalization (robust to outliers)
#' - **"PQN1"**: Probabilistic Quotient Normalization using global median reference
#' - **"PQN2"**: PQN using specific reference sample (`refSample`)
#' - **"groupPQN"**: PQN using pooled QC samples as reference
#' - **"quantile"**: Quantile normalization (forces identical distributions)
#'
#' **QC Normalization:** When using "Normalization" method, QC samples are
#' normalized separately using the mean or median of biological sample
#' normalization values (controlled by `qcNormalize`). This prevents
#' inappropriate correction of QC samples with biological sample-specific values.
#'
#' ## Scaling Methods
#'
#' - **"none"**: No scaling (preserves original magnitude relationships)
#' - **"mean"**: Mean-centering only (centers distribution at zero)
#' - **"auto"**: Auto-scaling = mean-centering + unit variance scaling
#'   (recommended for PCA, correlation, general multivariate analysis)
#' - **"pareto"**: Pareto-scaling = mean-centering + sqrt(SD) scaling
#'   (recommended for PLS-DA, OPLS-DA; balances between auto and none)
#'
#' ## Filter Reference Strategy
#'
#' The `filterReference` parameter determines which scaled dataset is used
#' as the reference for RSD and variance filtering:
#'
#' - **"NONPLS"** (default): Apply filters to NONPLS-scaled data, then use
#'   the resulting feature set for both NONPLS and PLS datasets. More stringent
#'   as auto-scaling amplifies small variations.
#'
#' - **"PLS"**: Apply filters to PLS-scaled data, then use the resulting
#'   feature set for both datasets. Less stringent as Pareto-scaling
#'   preserves more of the original structure.
#'
#' - **"auto"**: Apply filters independently to each scaled dataset and keep
#'   only features passing in BOTH (intersection). Most conservative approach.
#'
#' ## Technical Replicate Merging
#'
#' When `merge_replicates = TRUE`, the function automatically detects and
#' merges technical replicates based on the `SubjectID` column:
#'
#' - **What gets merged:** Only biological samples with matching SubjectID
#' - **What stays separate:** QC samples, blanks, and samples without SubjectID
#' - **Merging method:** Averages feature intensities across replicates
#' - **Metadata handling:** Uses earliest injection order, first batch encountered
#'
#' Merging is applied **after all filtering** to ensure only high-quality
#' features are included in the final averaged values.
#'
#' ## Parallelization (VSN transformation)
#'
#' The Variance Stabilizing Normalization (VSN) transformation can be
#' computationally intensive for large datasets. The `num_cores` parameter
#' enables parallel processing:
#'
#' - Set to integer (1 to max available cores) for specific core count
#' - Set to `"max"` for automatic detection (uses max cores - 2)
#' - Only affects VSN; other steps run sequentially
#'
#' ## Output Structure
#'
#' The function returns a comprehensive list containing:
#'
#' - **Metadata**: Sample information, groups, batches, injection order
#' - **Raw & intermediate data**: Data at each preprocessing step
#' - **Final scaled data**:
#'   - `data_scaledNONPLS_varFiltered`: For PCA, t-test, ANOVA, correlation
#'   - `data_scaledPLS_varFiltered`: For PLS-DA, OPLS-DA, sPLS-DA
#'   - `data_scaledNONPLS_merged`: Merged replicates (if applicable)
#'   - `data_scaledPLS_merged`: Merged replicates (if applicable)
#' - **Processing summary**: Dimensions at each step, parameters used
#' - **QC metrics**: Uncorrected features, outliers removed, time elapsed
#'
#' ## Important Notes
#'
#' 1. **Order is fixed**: Preprocessing steps execute in the order documented
#'    above, regardless of parameter order in the function call
#'
#' 2. **Use appropriate downstream data**:
#'    - For PCA, heatmaps, clustering: Use `data_scaledNONPLS_*`
#'    - For PLS-DA, OPLS-DA: Use `data_scaledPLS_*`
#'    - If replicates merged: Use `*_merged` datasets
#'
#' 3. **QC samples**: Retained throughout for quality assessment but should
#'    be removed before statistical analysis of biological groups
#'
#' 4. **Missing values**: After imputation, no NAs should remain. If new NAs
#'    appear after correction, they are automatically re-imputed.
#'
#' 5. **Feature count**: Expect substantial feature reduction through filtering
#'    steps. The `Dimensions` data frame tracks feature count at each stage.
#'
#' @author John Lennon L. Calorio
#'
#' @seealso \code{\link{perform_DataQualityCheck}}
#'
#' @importFrom stats quantile median sd
#' @importFrom parallel detectCores
#' @importFrom BiocParallel register SnowParam SerialParam
#' @importFrom matrixStats rowMedians colMedians colSds colMins colMaxs
#' @importFrom dplyr group_by summarise across all_of first
#' @importFrom parallelly availableCores
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
    removeUncorrectedFeatures   = FALSE,
    dataNormalize               = "Normalization",
    refSample                   = NULL,
    groupSample                 = NULL,
    reference_method            = "mean",
    dataTransform               = "vsn",
    dataScaleNONPLS             = "auto",
    dataScalePLS                = "pareto",
    qcNormalize                 = "median",
    filterReference             = "auto",
    filterMaxRSD                = 30,
    filterMaxRSD_by             = "EQC",
    filterLowVariability        = 10,
    merge_replicates            = TRUE,
    num_cores                   = 1,
    verbose                     = TRUE
) {

  start_time <- Sys.time()

  # Conditional messaging function
  msg <- function(...) if (verbose) message(...)
  utils::flush.console()
  msg("Starting data preprocessing pipeline...")

  # ============================================================================
  # INITIALIZE RESULTS LIST
  # ============================================================================

  listPreprocessed <- list(
    ProcessingTimestampStart = start_time,
    FunctionOrigin = "perform_PreprocessingPeakData",
    version = "2.0.0",
    Parameters = list(
      outliers = outliers, filterMissing = filterMissing, filterMissing_by_group = filterMissing_by_group,
      filterMissing_includeQC = filterMissing_includeQC, denMissing = denMissing,
      driftBatchCorrection = driftBatchCorrection, removeUncorrectedFeatures = removeUncorrectedFeatures,
      dataNormalize = dataNormalize, dataTransform = dataTransform, dataScaleNONPLS = dataScaleNONPLS,
      dataScalePLS = dataScalePLS, qcNormalize = qcNormalize, filterReference = filterReference, filterMaxRSD = filterMaxRSD,
      filterLowVariability = filterLowVariability, merge_replicates = merge_replicates,
      num_cores = num_cores
    ),
    Dimensions = data.frame(Step = character(), Samples = numeric(), Features = numeric(), stringsAsFactors = FALSE),
    ProcessingTimestampEnd = NULL,
    TotalProcessingTime_in_seconds = NULL,
    Error = NULL
  )

  tryCatch({

    # ============================================================================
    # INPUT VALIDATION
    # ============================================================================

    msg("Validating input data structure...")

    if (!is.list(raw_data) ||
        !identical(raw_data$FunctionOrigin, "perform_DataQualityCheck") ||
        !"raw_data" %in% names(raw_data) ||
        !"metadata_summary" %in% names(raw_data)) {
      stop("'raw_data' must be a list output from perform_DataQualityCheck function ",
           "containing 'raw_data' and 'metadata_summary'.")
    }

    # Get max number of cores available
    max_available_cores <- parallelly::availableCores(omit = 2)

    # --- Parameter validation ---
    validate_params <- function() {
      errors <- character()

      # Numeric parameters
      num_checks <- list(
        list(filterMissing, "filterMissing", 0, 100, TRUE),
        list(denMissing, "denMissing", 1, Inf, FALSE),
        list(spline_smooth_param, "spline_smooth_param", 0, 1, FALSE),
        list(min_QC, "min_QC", 1, Inf, FALSE),
        list(filterMaxRSD, "filterMaxRSD", 0, 100, TRUE),
        list(filterLowVariability, "filterLowVariability", 0, 100, TRUE)
        # list(num_cores, "num_cores", 1, max_available_cores, FALSE)
      )

      for (check in num_checks) {
        val <- check[[1]]; name <- check[[2]]; min_v <- check[[3]]; max_v <- check[[4]]; allow_null <- check[[5]]
        if (is.null(val) && allow_null) next
        if (!is.numeric(val) || length(val) != 1 || val < min_v || val > max_v) {
          errors <- c(errors, sprintf("%s: Must be numeric in [%g, %g]", name, min_v, max_v))
        }
      }

      # Logical parameters
      for (param in c("filterMissing_by_group", "filterMissing_includeQC", "driftBatchCorrection",
                      "log_scale", "removeUncorrectedFeatures", "merge_replicates")) {
        val <- get(param)
        if (!is.logical(val) || length(val) != 1) errors <- c(errors, paste(param, "must be a single logical (TRUE/FALSE)"))
      }

      # Choice parameters
      choices <- list(
        dataNormalize    = c("none", "Normalization", "sum", "median", "PQN1", "PQN2", "groupPQN", "quantile"),
        dataTransform    = c("none", "log2", "log10", "sqrt", "cbrt", "vsn", "glog"),
        dataScaleNONPLS  = c("none", "mean", "auto", "pareto"),
        dataScalePLS     = c("none", "mean", "auto", "pareto"),
        qcNormalize      = c("none", "mean", "median"),
        filterReference  = c("NONPLS", "PLS", "auto"),
        reference_method = c("mean", "median"),
        filterMaxRSD_by  = c("SQC", "EQC", "both")
      )

      for (param in names(choices)) {
        if (!get(param) %in% choices[[param]]) {
          errors <- c(errors, sprintf("%s: Must be one of: %s", param, paste(choices[[param]], collapse = ", ")))
        }
      }

      # Vector validations
      if (!is.null(outliers) && !is.vector(outliers)) errors <- c(errors, "outliers must be a vector")
      if (!is.numeric(spline_smooth_param_limit) || length(spline_smooth_param_limit) != 2) {
        errors <- c(errors, "spline_smooth_param_limit must be numeric vector of length 2")
      }

      if (length(errors) > 0) stop("Parameter validation failed:\n", paste(errors, collapse = "\n"))
      msg("Parameter validation passed.")
    }

    validate_params()

    # Check if num_cores is either max or integer
    if (num_cores != "max") {
      # Convert to integer if numeric and check range
      if (is.numeric(num_cores)) {
        num_cores <- as.integer(num_cores)
        if (num_cores < 1 || num_cores > max_available_cores) {
          stop(sprintf("'num_cores' must be between 1 and %d (available cores)", max_available_cores))
        }
      } else {
        stop("'num_cores' must either be 'max' or a positive integer")
      }
    }

    # Once all is validated, set final number of cores
    if (num_cores == "max") {
      num_cores <- max_available_cores
    }

    # Parallel setup
    BPPARAM_to_use <- BiocParallel::SnowParam(workers = num_cores)
    BiocParallel::register(BPPARAM_to_use)

    # Ensure SerialParam is restored even if the function errors
    on.exit({
      BiocParallel::register(BiocParallel::SerialParam())
    }, add = TRUE)

    # ============================================================================
    # DATA PREPARATION
    # ============================================================================

    msg("Loading pre-checked data...")

    data_processed <- data.frame(raw_data$raw_data)

    # Transpose and set row/column names
    data_processed <- data_processed %>% t() %>% `colnames<-`(., .[1, ]) %>% .[-1, ]
    data_processed <- data_processed %>% `rownames<-`(., .[, 1]) %>% as.data.frame()

    metadata_cols <- raw_data$metadata_summary$metadata_columns
    feature_cols  <- setdiff(colnames(data_processed), metadata_cols)

    # Handle merge_replicates check
    if (merge_replicates &&
        (!"Replicate" %in% colnames(data_processed) ||
         all(is.na(data_processed$Replicate) | data_processed$Replicate == "" | data_processed$Replicate == "0"))) {
      merge_replicates <- FALSE
      listPreprocessed$Parameters$merge_replicates <- FALSE
      msg("merge_replicates set to FALSE: no technical replicates found in 'Replicate' column.")
    }

    n_samples  <- nrow(data_processed)
    n_features <- length(feature_cols)
    msg(sprintf("Dataset contains %d samples and %d features/metabolites.", n_samples, n_features))

    listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                         data.frame(Step = "Original", Samples = n_samples, Features = n_features))

    # ============================================================================
    # OUTLIER REMOVAL
    # ============================================================================

    if (!is.null(outliers) && length(outliers) > 0) {
      msg("Removing specified outliers...")
      outliers_present <- outliers[outliers %in% data_processed$Sample]
      outliers_missing <- setdiff(outliers, outliers_present)

      if (length(outliers_present) > 0) {
        data_processed <- data_processed[!data_processed$Sample %in% outliers_present, ]
        msg(sprintf("Removed %d outliers: %s", length(outliers_present), paste(outliers_present, collapse = ", ")))
      }
      if (length(outliers_missing) > 0) {
        msg(sprintf("Outliers not found: %s", paste(outliers_missing, collapse = ", ")))
      }

      listPreprocessed$outliers_removed   <- outliers_present
      listPreprocessed$outliers_not_found <- outliers_missing
      listPreprocessed$Dimensions         <- rbind(listPreprocessed$Dimensions,
                                                   data.frame(Step = "After outlier removal", Samples = nrow(data_processed), Features = length(feature_cols)))
    }

    # ============================================================================
    # METADATA EXTRACTION
    # ============================================================================

    msg("Extracting metadata...")

    # Helper to safely get a column or return NA if it doesn't exist
    get_col_or_na <- function(df, colname) {
      if (colname %in% colnames(df)) df[[colname]] else rep(NA, nrow(df))
    }

    listPreprocessed$Metadata <- data.frame(
      Samples             = get_col_or_na(data_processed, "Sample"),
      SubjectID           = suppressWarnings(as.character(get_col_or_na(data_processed, "SubjectID"))), # Keep as char for merging
      TechnicalReplicates = get_col_or_na(data_processed, "Replicate"), # Renamed
      Group               = get_col_or_na(data_processed, "Group"),
      Group2              = get_col_or_na(data_processed, "Group2"),
      Group_              = gsub("SQC|EQC", "QC", get_col_or_na(data_processed, "Group")),
      Batches             = suppressWarnings(as.numeric(get_col_or_na(data_processed, "Batch"))),
      InjectionSequence   = suppressWarnings(as.numeric(get_col_or_na(data_processed, "Injection"))),
      Normalization       = suppressWarnings(as.numeric(get_col_or_na(data_processed, "Normalization"))),
      Response            = suppressWarnings(as.numeric(get_col_or_na(data_processed, "Response"))),
      stringsAsFactors    = FALSE
    )
    rownames(listPreprocessed$Metadata) <- listPreprocessed$Metadata$Samples

    indices_qc <- listPreprocessed$Metadata$Group %in% c("SQC", "EQC", "QC")
    indices_non_qc <- !indices_qc

    # ============================================================================
    # FEATURE DATA EXTRACTION
    # ============================================================================

    msg("Extracting feature data...")

    df                                        <- data_processed[, feature_cols, drop = FALSE]
    rownames(df)                              <- data_processed$Sample

    listPreprocessed$All_Features_Metabolites <- colnames(df)

    # Process missing values (DataQualityCheck converts to numeric, 0s are 0s)
    msg("Processing missing values (0 -> NA)...")
    df[df == 0] <- NA

    # Ensure all columns are numeric
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))

    listPreprocessed$n_missing <- sum(is.na(df))
    msg(sprintf("Found %d missing values (0s or NAs) in the data.", listPreprocessed$n_missing))

    # Remove all-NA columns
    all_na_cols <- colSums(is.na(df), na.rm = TRUE) == nrow(df)
    if (any(all_na_cols)) {
      n_removed <- sum(all_na_cols)
      df        <- df[, !all_na_cols, drop = FALSE]
      msg(sprintf("Removed %d features with all missing values.", n_removed))
    }

    listPreprocessed$SamplesXFeatures <- df
    listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                         data.frame(Step = "After removing all-NA features", Samples = nrow(df), Features = ncol(df)))

    # ============================================================================
    # MISSING VALUE FILTERING
    # ============================================================================

    if (!is.null(filterMissing)) {
      msg(sprintf("Applying missing value filter (threshold: %g%%)", filterMissing))

      if (filterMissing_by_group) {
        groups <- if (filterMissing_includeQC) {
          unique(listPreprocessed$Metadata$Group)
        } else {
          unique(listPreprocessed$Metadata$Group[!indices_qc])
        }

        missing_by_group <- sapply(groups, function(g) {
          group_indices <- listPreprocessed$Metadata$Group == g
          if (sum(group_indices) == 0) return(rep(0, ncol(df)))
          colMeans(is.na(df[group_indices, , drop = FALSE]))
        })

        if (is.vector(missing_by_group)) { # Only one group
          features_to_keep <- missing_by_group < (filterMissing / 100)
        } else { # Multiple groups
          features_to_keep <- !apply(missing_by_group, 1, function(x) all(x >= filterMissing / 100))
        }

      } else {
        if (filterMissing_includeQC) {
          missing_overall <- colMeans(is.na(df))
        } else {
          missing_overall <- colMeans(is.na(df[indices_non_qc, , drop = FALSE]))
        }
        features_to_keep <- missing_overall < (filterMissing / 100)
      }

      n_before <- ncol(df)
      df <- df[, features_to_keep, drop = FALSE]
      msg(sprintf("Removed %d features due to missing value threshold.", n_before - ncol(df)))
    } else {
      msg("Skipping missing value filtering (filterMissing = NULL)")
    }

    listPreprocessed$data_filteredMissing <- df
    listPreprocessed$Dimensions           <- rbind(listPreprocessed$Dimensions,
                                                   data.frame(Step = "After missing value filter", Samples = nrow(df), Features = ncol(df)))

    # ============================================================================
    # MISSING VALUE IMPUTATION
    # ============================================================================

    msg("Performing missing value imputation...")

    if (ncol(df) > 0) {
      # Use matrixStats for a fast, vectorized column minimum
      min_vals                       <- matrixStats::colMins(as.matrix(df), na.rm = TRUE) / denMissing

      # Handle features that were all NA (min_val will be Inf)
      min_vals[!is.finite(min_vals)] <- 1e-9 # Impute with a small number

      # Fast matrix-based NA replacement
      na_indices                     <- is.na(df)
      df[na_indices]                 <- min_vals[col(df)][na_indices]
    }

    listPreprocessed$data_no_NA <- df
    msg(sprintf("Imputed missing values using 1/%d of minimum positive value per feature.", denMissing))

    # ============================================================================
    # SIGNAL DRIFT AND BATCH CORRECTION
    # ============================================================================

    uncorrected_features <- character(0)
    features_removed_uncorrected <- 0

    if (driftBatchCorrection) {
      if (!requireNamespace("pmp", quietly = TRUE)) {
        warning("Package 'pmp' not available. Skipping drift/batch correction.")
        driftBatchCorrection <- FALSE
        listPreprocessed$Parameters$driftBatchCorrection <- FALSE
      }
    }

    if (driftBatchCorrection) {
      tryCatch({
        df_before_correction <- df

        df_corrected <- pmp::QCRSC(
          df       = df,
          order    = listPreprocessed$Metadata$InjectionSequence,
          batch    = listPreprocessed$Metadata$Batches,
          classes  = as.vector(listPreprocessed$Metadata$Group_),
          spar     = spline_smooth_param,
          log      = log_scale,
          minQC    = min_QC,
          qc_label = "QC",
          spar_lim = spline_smooth_param_limit
        ) %>% t() %>% as.data.frame()

        msg("Applied QC-RSC drift and batch correction.")

        # Identify uncorrected features
        tryCatch({
          msg("Identifying uncorrected features...")

          if (identical(dim(df_before_correction), dim(df_corrected)) &&
              identical(colnames(df_before_correction), colnames(df_corrected))) {

            # Use a tolerance that captures "essentially unchanged" values
            # This catches features where QCRSC returned original values
            tolerance            <- .Machine$double.eps^0.5

            # Check if values are essentially identical (unchanged by QCRSC)
            diff_matrix          <- abs(as.matrix(df_before_correction) - as.matrix(df_corrected))
            col_max_diffs        <- matrixStats::colMaxs(diff_matrix, na.rm = TRUE)

            # A feature is uncorrected if maximum change is below tolerance
            # This indicates QCRSC left it unchanged (likely due to insufficient QCs)
            uncorrected_indices  <- (col_max_diffs <= tolerance)

            uncorrected_features <- colnames(df_before_correction)[uncorrected_indices]
            msg(sprintf("Identified %d uncorrected features.", length(uncorrected_features)))

            listPreprocessed$UncorrectedFeatures <- list(
              feature_names       = uncorrected_features,
              count               = length(uncorrected_features),
              max_absolute_change = col_max_diffs,
              reason              = paste("Insufficient QC samples (< min_QC =", min_QC,
                                          ") for reliable spline fitting")
            )

            if (removeUncorrectedFeatures && length(uncorrected_features) > 0) {
              features_to_keep             <- !uncorrected_indices
              df_corrected                 <- df_corrected[, features_to_keep, drop = FALSE]
              df_before_correction         <- df_before_correction[, features_to_keep, drop = FALSE]
              features_removed_uncorrected <- length(uncorrected_features)
              msg(sprintf("Removed %d uncorrected features.", features_removed_uncorrected))
            } else if (!removeUncorrectedFeatures && length(uncorrected_features) > 0) {
              msg(sprintf("Keeping %d uncorrected features.", length(uncorrected_features)))
            }
          } else {
            warning("Dimension/column mismatch. Cannot identify uncorrected features.")
          }
        }, error = function(e) {
          warning("Error identifying uncorrected features: ", e$message)
        })

        df <- df_corrected

        # Re-impute any new NAs
        new_na_count <- sum(is.na(df))
        if (new_na_count > 0) {
          msg(sprintf("Correction introduced %d new missing values. Re-imputing...", new_na_count))
          min_vals                       <- matrixStats::colMins(as.matrix(df), na.rm = TRUE) / denMissing
          min_vals[!is.finite(min_vals)] <- 1e-9
          na_indices                     <- is.na(df)
          df[na_indices]                 <- min_vals[col(df)][na_indices]
          listPreprocessed$new_n_missing <- new_na_count
        }

      }, error = function(e) {
        warning("Drift/batch correction failed: ", e$message, ". Proceeding without correction.")
        driftBatchCorrection <- FALSE
        listPreprocessed$Parameters$driftBatchCorrection <- FALSE
      })
    }

    if (!driftBatchCorrection) {
      msg("Skipping drift and batch correction.")
      listPreprocessed$UncorrectedFeatures <- list(
        feature_names = character(0), count = 0,
        reason        = "Drift/batch correction was not performed"
      )
    }

    listPreprocessed$data_driftBatchCorrected <- df
    listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                         data.frame(Step = "After drift correction", Samples = nrow(df), Features = ncol(df)))

    # ============================================================================
    # NORMALIZATION
    # ============================================================================

    msg(sprintf("Applying %s normalization...", dataNormalize))

    # Convert to matrix for fast row/col operations
    df_matrix <- as.matrix(df)

    df <- switch(dataNormalize,
                 "none" = {
                   msg("No normalization applied.")
                   df
                 },
                 "Normalization" = {
                   sg_values <- as.numeric(listPreprocessed$Metadata$Normalization[indices_non_qc])
                   if (all(is.na(sg_values) | sg_values == 0)) {
                     msg("Normalization values are all NA/0. Using normalization by 'sum' instead.")
                     df_matrix / rowSums(df_matrix, na.rm = TRUE)
                   } else {
                     # Normalize biological samples
                     df_matrix[indices_non_qc, ] <- df_matrix[indices_non_qc, ] / sg_values

                     # For QC samples
                     if (qcNormalize == "mean") {
                       qc_norm_value <- mean(sg_values, na.rm = TRUE)
                       msg("Applied 'normalization using provided values'. QCs normalized using 'mean' of 'Normalization' values.")
                     } else if (qcNormalize == "median") {
                       qc_norm_value <- median(sg_values, na.rm = TRUE)
                       msg("Applied 'normalization using provided values'. QCs normalized using 'median' of 'Normalization' values.")
                     } else if (qcNormalize == "none") {
                       qc_norm_value <- 1
                       msg("Applied 'normalization using provided values'. QCs are not normalized.")
                     }
                     df_matrix[indices_qc, ] <- df_matrix[indices_qc, ] / qc_norm_value
                     df_matrix
                   }
                 },
                 "sum" = {
                   msg("Applied 'sum' normalization.")
                   df_matrix / rowSums(df_matrix, na.rm = TRUE)
                 },
                 "median" = {
                   msg("Applied 'median' normalization.")
                   df_matrix / matrixStats::rowMedians(df_matrix, na.rm = TRUE)
                 },
                 "PQN1" = {
                   msg("Applied 'PQN' (global median).")
                   reference_spectrum <- matrixStats::colMedians(df_matrix, na.rm = TRUE)
                   quotients          <- sweep(df_matrix, 2, reference_spectrum, "/")
                   median_quotients   <- matrixStats::rowMedians(quotients, na.rm = TRUE)
                   df_matrix / median_quotients
                 },
                 "PQN2" = {
                   if (is.null(refSample) || !refSample %in% rownames(df_matrix)) {
                     stop("Reference sample not specified or not found for PQN2 normalization.")
                   }
                   msg(sprintf("Applied 'PQN' (ref sample %s).", refSample))
                   reference_spectrum <- as.numeric(df_matrix[refSample, ])
                   quotients          <- sweep(df_matrix, 2, reference_spectrum, "/")
                   median_quotients   <- matrixStats::rowMedians(quotients, na.rm = TRUE)
                   df_matrix / median_quotients
                 },
                 "groupPQN" = {
                   if (is.null(groupSample)) groupSample <- "both"
                   pooled_indices <- if (groupSample == "both") {
                     listPreprocessed$Metadata$Group %in% c("SQC", "EQC", "QC")
                   } else {
                     listPreprocessed$Metadata$Group == groupSample
                   }
                   if (sum(pooled_indices) == 0) stop("No samples found for specified group sample.")
                   msg(sprintf("Applied 'group PQN' (group %s).", groupSample))
                   reference_spectrum <- matrixStats::colMedians(df_matrix[pooled_indices, , drop = FALSE], na.rm = TRUE)
                   quotients          <- sweep(df_matrix, 2, reference_spectrum, "/")
                   median_quotients   <- matrixStats::rowMedians(quotients, na.rm = TRUE)
                   df_matrix / median_quotients
                 },
                 "quantile" = {
                   if (!requireNamespace("pmp", quietly = TRUE)) {
                     stop("Package 'pmp' required for quantile normalization.")
                   }
                   tryCatch({
                     normalized <- pmp::pqn_normalisation(
                       t(df_matrix), classes = listPreprocessed$Metadata$Group_,
                       qc_label = "QC", ref_mean = NULL, qc_frac = 0,
                       sample_frac = 0, ref_method = reference_method
                     )
                     msg("Applied 'quantile' normalization.")
                     t(normalized)
                   }, error = function(e) {
                     warning("'Quantile' normalization failed: ", e$message, ". Using 'sum' normalization.")
                     df_matrix / rowSums(df_matrix, na.rm = TRUE)
                   })
                 },
                 {
                   warning("Unknown normalization method. No normalization applied.")
                   df_matrix
                 }
    )

    # Convert back to data.frame
    df <- as.data.frame(df)

    # Handle Infs produced by division by zero (if any)
    if (any(is.infinite(as.matrix(df)))) {
      msg("Replacing Inf values with NA (and re-imputing).")
      df[is.infinite(as.matrix(df))] <- NA
      # Re-impute NAs
      min_vals                       <- matrixStats::colMins(as.matrix(df), na.rm = TRUE) / denMissing
      min_vals[!is.finite(min_vals)] <- 1e-9
      na_indices                     <- is.na(df)
      df[na_indices]                 <- min_vals[col(df)][na_indices]
    }

    listPreprocessed$data_normalized <- df

    # ============================================================================
    # TRANSFORMATION
    # ============================================================================

    msg(sprintf("Applying %s transformation...", dataTransform))

    df <- switch(dataTransform,
                 "none" = {
                   msg("No transformation applied.")
                   df
                 },
                 "log2" = {
                   min_val              <- min(df, na.rm = TRUE)
                   if (min_val <= 0) df <- df - min_val + 1
                   msg("Applied 'log2' transformation.")
                   log2(df)
                 },
                 "log10" = {
                   min_val              <- min(df, na.rm = TRUE)
                   if (min_val <= 0) df <- df - min_val + 1
                   msg("Applied 'log10' transformation.")
                   log10(df)
                 },
                 "sqrt" = {
                   min_val             <- min(df, na.rm = TRUE)
                   if (min_val < 0) df <- df - min_val
                   msg("Applied 'square root' transformation.")
                   sqrt(df)
                 },
                 "cbrt" = {
                   min_val             <- min(df, na.rm = TRUE)
                   if (min_val < 0) df <- df - min_val + 1
                   msg("Applied 'cube root' transformation.")
                   df^(1/3)
                 },
                 "vsn" = {
                   if (!requireNamespace("vsn", quietly = TRUE)) {
                     warning("Package 'vsn' not available. Using 'log10' transformation instead.")
                     min_val <- min(df, na.rm = TRUE); if (min_val <= 0) df <- df - min_val + 1; log10(df)
                   } else {
                     # tryCatch({
                     #   BPPARAM_to_use <- BiocParallel::SnowParam(workers = num_cores)
                     #   BiocParallel::register(BPPARAM_to_use)
                     #
                     #   vsn_fit <- vsn::vsn2(as.matrix(df))
                     #   result <- vsn::predict(vsn_fit, newdata = as.matrix(df))
                     #   msg("Applied VSN transformation.")
                     #   as.data.frame(result)
                     #
                     # }, error = function(e) {
                     #   warning("VSN transformation failed: ", e$message, ". Using log10 instead.")
                     #   min_val <- min(df, na.rm = TRUE); if (min_val <= 0) df <- df - min_val + 1; log10(df)
                     # }, finally = {
                     #   BiocParallel::register(BiocParallel::SerialParam()) # Always reset
                     # })
                     tryCatch({
                       vsn_fit <- vsn::vsn2(as.matrix(df))
                       result  <- vsn::predict(vsn_fit, newdata = as.matrix(df))
                       msg("Applied 'VSN' transformation.")
                       result  <- as.data.frame(result)

                     }, error = function(e) {
                       warning("'VSN' transformation failed: ", e$message, ". Using 'log10' instead.")
                       min_val              <- min(df, na.rm = TRUE)
                       if (min_val <= 0) df <- df - min_val + 1
                       result               <- log10(df)
                     })
                   }
                 },
                 "glog" = {
                   if (!requireNamespace("pmp", quietly = TRUE)) {
                     warning("Package 'pmp' not available. Using 'log10' transformation instead.")
                     min_val <- min(df, na.rm = TRUE); if (min_val <= 0) df <- df - min_val + 1; log10(df)
                   } else {
                     tryCatch({
                       glog_result <- pmp::glog_transformation(
                         t(as.matrix(df)),
                         classes  = listPreprocessed$Metadata$Group_,
                         qc_label = "QC"
                       )
                       msg("Applied 'glog' transformation.")
                       as.data.frame(t(glog_result))
                     }, error = function(e) {
                       warning("'glog' transformation failed: ", e$message, ". Using 'log10' instead.")
                       min_val <- min(df, na.rm = TRUE); if (min_val <= 0) df <- df - min_val + 1; log10(df)
                     })
                   }
                 },
                 {
                   warning("Unknown transformation method. No transformation applied.")
                   df
                 }
    )

    listPreprocessed$data_transformed <- df

    # ============================================================================
    # SCALING
    # ============================================================================

    scale_data <- function(data, method, purpose) {
      msg(sprintf("Applying %s scaling for %s analysis...", method, purpose))
      data_matrix <- as.matrix(data)

      switch(method,
             "none" = {
               msg(sprintf("No scaling applied for %s data.", purpose))
               as.data.frame(data)
             },
             "mean" = {
               msg(sprintf("Applied 'mean centering' for %s analysis.", purpose))
               as.data.frame(scale(data_matrix, center = TRUE, scale = FALSE))
             },
             "auto" = {
               msg(sprintf("Applied 'mean centering and unit variance' scaling for %s analysis.", purpose))
               as.data.frame(scale(data_matrix, center = TRUE, scale = TRUE))
             },
             "pareto" = {
               msg(sprintf("Applied 'Pareto scaling' for %s analysis.", purpose))
               sds           <- matrixStats::colSds(data_matrix, na.rm = TRUE)
               sds[sds == 0] <- 1 # Avoid division by zero
               as.data.frame(scale(data_matrix, center = TRUE, scale = sqrt(sds)))
             },
             {
               warning("Unknown scaling method. No scaling applied.")
               as.data.frame(data)
             }
      )
    }

    listPreprocessed$data_scaledNONPLS <- scale_data(df, dataScaleNONPLS, "NONPLS")
    listPreprocessed$data_scaledPLS    <- scale_data(df, dataScalePLS, "PLS")

    # ============================================================================
    # RSD FILTERING
    # ============================================================================

    apply_rsd_filter <- function(data, threshold, filter_by, metadata, data_for) {
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
        # Ensure sample alignment
        data_samples   <- rownames(data)
        meta_samples   <- metadata$Samples

        common_samples <- intersect(data_samples, meta_samples)
        if (length(common_samples) == 0) {
          warning("No matching samples found. Skipping RSD filtering.")
          return(data)
        }

        data_aligned <- data[common_samples, , drop = FALSE]
        meta_aligned <- metadata[match(common_samples, meta_samples), , drop = FALSE]

        # Determine QC label and classes
        if (filter_by == "both") {
          qc_label <- "QC"
          classes  <- meta_aligned$Group_
        } else {
          qc_label <- filter_by
          classes  <- meta_aligned$Group
        }

        if (!qc_label %in% classes) {
          warning(sprintf("No %s samples found for RSD filtering. Skipping.", qc_label))
          return(data)
        }

        filtered_data <- pmp::filter_peaks_by_rsd(
          t(data_aligned), max_rsd = threshold,
          class                    = as.vector(classes),
          qc_label                 = qc_label
        )

        result <- as.data.frame(t(filtered_data))
        n_removed <- ncol(data) - ncol(result)
        msg(sprintf("Removed %d features with RSD > %g%% for %s analysis.", n_removed, threshold, data_for))

        # Return data with original sample set, but filtered features
        return(data[, colnames(result), drop = FALSE])

      }, error = function(e) {
        warning("RSD filtering failed: ", e$message, ". Returning original data.")
        data
      })
    }

    # ============================================================================
    # LOW-VARIANCE FILTERING
    # ============================================================================

    apply_var_filter <- function(data, threshold, data_for) {
      if (is.null(threshold)) {
        msg(sprintf("No variance filtering applied for %s data.", data_for))
        return(data)
      }

      msg(sprintf("Applying variance filtering (%gth percentile) for %s data...", threshold, data_for))

      tryCatch({
        sds <- matrixStats::colSds(as.matrix(data), na.rm = TRUE)
        valid_sds <- is.finite(sds) & !is.na(sds) & sds > 0

        if (sum(valid_sds) == 0) {
          warning(sprintf("No valid standard deviations for %s data. Returning original.", data_for))
          return(data)
        }

        threshold_value <- quantile(sds[valid_sds], threshold/100, na.rm = TRUE)
        keep_features   <- (sds > threshold_value) & valid_sds

        result          <- data[, keep_features, drop = FALSE]
        n_removed       <- ncol(data) - ncol(result)
        msg(sprintf("Removed %d features with low variance for %s analysis.", n_removed, data_for))
        result

      }, error = function(e) {
        warning("Variance filtering failed: ", e$message, ". Returning original data.")
        data
      })
    }

    # ============================================================================
    # APPLY FILTERING BASED ON REFERENCE
    # ============================================================================

    if (filterReference == "NONPLS") {
      msg("Using NONPLS-scaled data as reference for filtering...")

      listPreprocessed$data_scaledNONPLS_rsdFiltered <- apply_rsd_filter(
        listPreprocessed$data_scaledNONPLS, filterMaxRSD, filterMaxRSD_by,
        listPreprocessed$Metadata, "NONPLS"
      )

      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After RSD filtering (NONPLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledNONPLS_rsdFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledNONPLS_rsdFiltered)))

      listPreprocessed$data_scaledNONPLS_varFiltered <- apply_var_filter(
        listPreprocessed$data_scaledNONPLS_rsdFiltered, filterLowVariability, "NONPLS"
      )

      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After variance filtering (NONPLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledNONPLS_varFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledNONPLS_varFiltered)))

      # Apply NONPLS-filtered feature set to PLS (no actual filtering on PLS)
      features_to_keep                            <- colnames(listPreprocessed$data_scaledNONPLS_varFiltered)
      listPreprocessed$data_scaledPLS_rsdFiltered <- listPreprocessed$data_scaledPLS[, features_to_keep, drop = FALSE]
      listPreprocessed$data_scaledPLS_varFiltered <- listPreprocessed$data_scaledPLS_rsdFiltered

      msg(sprintf("Applied NONPLS feature set (%d features) to PLS data (no independent PLS filtering).",
                  length(features_to_keep)))

      # Record PLS dimensions (same feature set as NONPLS, not independently filtered)
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After applying NONPLS feature set (PLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledPLS_varFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledPLS_varFiltered)))

    } else if (filterReference == "PLS") {
      msg("Using PLS-scaled data as reference for filtering...")

      listPreprocessed$data_scaledPLS_rsdFiltered <- apply_rsd_filter(
        listPreprocessed$data_scaledPLS, filterMaxRSD, filterMaxRSD_by,
        listPreprocessed$Metadata, "PLS"
      )

      # Record dimensions after RSD filtering
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After RSD filtering (PLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledPLS_rsdFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledPLS_rsdFiltered)))

      listPreprocessed$data_scaledPLS_varFiltered    <- apply_var_filter(
        listPreprocessed$data_scaledPLS_rsdFiltered, filterLowVariability, "PLS"
      )

      # Record dimensions after variance filtering
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After variance filtering (PLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledPLS_varFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledPLS_varFiltered)))

      features_to_keep                               <- colnames(listPreprocessed$data_scaledPLS_varFiltered)
      listPreprocessed$data_scaledNONPLS_rsdFiltered <- listPreprocessed$data_scaledNONPLS[, features_to_keep, drop = FALSE]
      listPreprocessed$data_scaledNONPLS_varFiltered <- listPreprocessed$data_scaledNONPLS_rsdFiltered

      # Record NONPLS dimensions (same feature set as PLS)
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After RSD filtering (NONPLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledNONPLS_rsdFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledNONPLS_rsdFiltered)))
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After variance filtering (NONPLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledNONPLS_varFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledNONPLS_varFiltered)))

    } else if (filterReference == "auto") {
      msg("Applying independent filtering and keeping intersection...")

      data_nonpls_rsd  <- apply_rsd_filter(listPreprocessed$data_scaledNONPLS, filterMaxRSD,
                                           filterMaxRSD_by, listPreprocessed$Metadata, "NONPLS")
      data_nonpls_var  <- apply_var_filter(data_nonpls_rsd, filterLowVariability, "NONPLS")

      data_pls_rsd     <- apply_rsd_filter(listPreprocessed$data_scaledPLS, filterMaxRSD,
                                           filterMaxRSD_by, listPreprocessed$Metadata, "PLS")
      data_pls_var     <- apply_var_filter(data_pls_rsd, filterLowVariability, "PLS")

      # Record dimensions after RSD filtering (before intersection)
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After RSD filtering (NONPLS, pre-intersection)",
                                                      Samples  = nrow(data_nonpls_rsd),
                                                      Features = ncol(data_nonpls_rsd)))
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After RSD filtering (PLS, pre-intersection)",
                                                      Samples  = nrow(data_pls_rsd),
                                                      Features = ncol(data_pls_rsd)))

      # Record dimensions after variance filtering (before intersection)
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After variance filtering (NONPLS, pre-intersection)",
                                                      Samples  = nrow(data_nonpls_var),
                                                      Features = ncol(data_nonpls_var)))
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After variance filtering (PLS, pre-intersection)",
                                                      Samples  = nrow(data_pls_var),
                                                      Features = ncol(data_pls_var)))

      features_to_keep <- intersect(colnames(data_nonpls_var), colnames(data_pls_var))

      if (length(features_to_keep) == 0) {
        stop("No features passed filtering in both datasets. Consider less stringent parameters.")
      }

      listPreprocessed$data_scaledNONPLS_rsdFiltered <- listPreprocessed$data_scaledNONPLS[, features_to_keep, drop = FALSE]
      listPreprocessed$data_scaledNONPLS_varFiltered <- listPreprocessed$data_scaledNONPLS_rsdFiltered
      listPreprocessed$data_scaledPLS_rsdFiltered    <- listPreprocessed$data_scaledPLS[, features_to_keep, drop = FALSE]
      listPreprocessed$data_scaledPLS_varFiltered    <- listPreprocessed$data_scaledPLS_rsdFiltered

      # Record final dimensions after intersection
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After intersection (NONPLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledNONPLS_varFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledNONPLS_varFiltered)))
      listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                           data.frame(Step     = "After intersection (PLS)",
                                                      Samples  = nrow(listPreprocessed$data_scaledPLS_varFiltered),
                                                      Features = ncol(listPreprocessed$data_scaledPLS_varFiltered)))
    }

    # Create filtered transformed data for downstream use
    features_kept                                 <- colnames(listPreprocessed$data_scaledNONPLS_varFiltered)
    listPreprocessed$data_normalized_varFiltered  <- listPreprocessed$data_normalized[, features_kept, drop = FALSE]
    listPreprocessed$data_transformed_varFiltered <- listPreprocessed$data_transformed[, features_kept, drop = FALSE]

    # Update final filtering dimensions summary (keep these for backward compatibility)
    listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                         data.frame(Step     = "After final filtering (NONPLS)",
                                                    Samples  = nrow(listPreprocessed$data_scaledNONPLS_varFiltered),
                                                    Features = ncol(listPreprocessed$data_scaledNONPLS_varFiltered)))
    listPreprocessed$Dimensions <- rbind(listPreprocessed$Dimensions,
                                         data.frame(Step     = "After final filtering (PLS)",
                                                    Samples  = nrow(listPreprocessed$data_scaledPLS_varFiltered),
                                                    Features = ncol(listPreprocessed$data_scaledPLS_varFiltered)))

    # ============================================================================
    # TECHNICAL REPLICATE MERGING
    # ============================================================================

    merge_technical_replicates <- function(data_nonpls, data_pls, metadata, verbose = TRUE) {
      msg <- function(...) if (verbose) message(...)

      tryCatch({
        msg("Checking for technical replicates to merge...")

        has_subject_id <- !all(is.na(metadata$SubjectID) | metadata$SubjectID == "")
        has_replicates <- !all(is.na(metadata$TechnicalReplicates) | metadata$TechnicalReplicates == "")

        if (!has_subject_id || !has_replicates) {
          msg("No meaningful SubjectID/Replicate columns found. Skipping replicate merging.")
          return(list(merged = FALSE, reason = "No SubjectID/Replicate data"))
        }

        # Identify the 3 groups of samples:
        # 1. QCs (never merged)
        # 2. Non-QCs with valid SubjectID (to be merged)
        # 3. Non-QCs *without* valid SubjectID (e.g., Blanks) (never merged)

        qc_indices                   <- metadata$Group %in% c("SQC", "EQC", "QC")
        mergeable_indices            <- !qc_indices & !is.na(metadata$SubjectID) & metadata$SubjectID != ""
        non_mergeable_non_qc_indices <- !qc_indices & (is.na(metadata$SubjectID) | metadata$SubjectID == "")

        n_to_merge <- sum(mergeable_indices)
        if (n_to_merge == 0 || length(unique(metadata$SubjectID[mergeable_indices])) == n_to_merge) {
          msg("No technical replicates found to merge. Skipping.")
          return(list(merged = FALSE, reason = "No actual replicates to merge"))
        }

        msg("Technical replicates detected. Proceeding with merging...")

        # Combine metadata and feature data for easier processing
        comp_data_nonpls <- cbind(metadata, data_nonpls)
        comp_data_pls    <- cbind(metadata, data_pls)

        # --- 1. Get data for mergeable samples ---
        data_to_merge_nonpls <- comp_data_nonpls[mergeable_indices, , drop = FALSE]
        data_to_merge_pls    <- comp_data_pls[mergeable_indices, , drop = FALSE]

        # --- 2. Get data for samples to keep as-is (QCs and non-mergeable) ---
        data_to_keep_nonpls <- comp_data_nonpls[!mergeable_indices, , drop = FALSE]
        data_to_keep_pls    <- comp_data_pls[!mergeable_indices, , drop = FALSE]

        # Set 'MergeID' for kept samples to their original Sample name
        data_to_keep_nonpls$MergeID <- data_to_keep_nonpls$Samples
        data_to_keep_pls$MergeID    <- data_to_keep_pls$Samples

        # --- 3. Perform the merge using dplyr ---
        if (!requireNamespace("dplyr", quietly = TRUE)) {
          stop("Package 'dplyr' is required for merging replicates.")
        }

        # Get feature column names (not metadata)
        feature_cols_nonpls <- colnames(data_to_merge_nonpls)[!colnames(data_to_merge_nonpls) %in% colnames(metadata)]
        feature_cols_pls    <- colnames(data_to_merge_pls)[!colnames(data_to_merge_pls) %in% colnames(metadata)]

        # IMPORTANT: Group by BOTH SubjectID AND Group to preserve group identity
        merged_non_qc_nonpls <- data_to_merge_nonpls %>%
          dplyr::group_by(SubjectID, Group) %>%
          dplyr::summarise(
            TechnicalReplicates   = dplyr::first(TechnicalReplicates),
            Group2                = dplyr::first(Group2),
            Batches               = dplyr::first(Batches),
            InjectionSequence     = min(InjectionSequence, na.rm = TRUE),
            Normalization         = dplyr::first(Normalization),
            Response              = dplyr::first(Response),
            dplyr::across(dplyr::all_of(feature_cols_nonpls), ~mean(., na.rm = TRUE)),
            .groups               = 'drop'
          )

        merged_non_qc_pls <- data_to_merge_pls %>%
          dplyr::group_by(SubjectID, Group) %>%
          dplyr::summarise(
            TechnicalReplicates   = dplyr::first(TechnicalReplicates),
            Group2                = dplyr::first(Group2),
            Batches               = dplyr::first(Batches),
            InjectionSequence     = min(InjectionSequence, na.rm = TRUE),
            Normalization         = dplyr::first(Normalization),
            Response              = dplyr::first(Response),
            dplyr::across(dplyr::all_of(feature_cols_pls), ~mean(., na.rm = TRUE)),
            .groups               = 'drop'
          )

        # Set 'MergeID' for merged samples to their TechnicalReplicates value
        merged_non_qc_nonpls$MergeID <- merged_non_qc_nonpls$TechnicalReplicates
        merged_non_qc_pls$MergeID    <- merged_non_qc_pls$TechnicalReplicates

        # --- 4. Combine merged and kept data ---
        # Align columns before rbind
        cols_to_keep      <- colnames(merged_non_qc_nonpls)

        final_data_nonpls <- rbind(merged_non_qc_nonpls, data_to_keep_nonpls[, cols_to_keep, drop = FALSE])
        final_data_pls    <- rbind(merged_non_qc_pls, data_to_keep_pls[, cols_to_keep, drop = FALSE])

        # --- 5. Separate final metadata and features ---
        # Identify metadata columns (including MergeID)
        meta_col_names <- c(colnames(metadata), "MergeID")

        # Get feature columns (everything that's not metadata)
        final_feature_cols_nonpls <- colnames(final_data_nonpls)[!colnames(final_data_nonpls) %in% meta_col_names]
        final_feature_cols_pls    <- colnames(final_data_pls)[!colnames(final_data_pls) %in% meta_col_names]

        merged_features_nonpls           <- as.data.frame(final_data_nonpls[, final_feature_cols_nonpls, drop = FALSE])
        merged_features_pls              <- as.data.frame(final_data_pls[, final_feature_cols_pls, drop = FALSE])
        rownames(merged_features_nonpls) <- final_data_nonpls$MergeID
        rownames(merged_features_pls)    <- final_data_pls$MergeID

        merged_metadata <- data.frame(
          Samples             = final_data_nonpls$MergeID,
          SubjectID           = final_data_nonpls$SubjectID,
          TechnicalReplicates = final_data_nonpls$TechnicalReplicates,
          Group               = final_data_nonpls$Group,
          Group2              = final_data_nonpls$Group2,
          Group_              = gsub("SQC|EQC", "QC", final_data_nonpls$Group),
          Batches             = final_data_nonpls$Batches,
          InjectionSequence   = final_data_nonpls$InjectionSequence,
          Normalization       = final_data_nonpls$Normalization,
          Response            = final_data_nonpls$Response,
          stringsAsFactors    = FALSE
        )
        rownames(merged_metadata) <- merged_metadata$Samples

        n_samples_before <- nrow(metadata)
        n_samples_after  <- nrow(merged_metadata)

        msg(sprintf("Successfully merged technical replicates: %d samples -> %d samples",
                    n_samples_before, n_samples_after))

        return(list(
          data_nonpls      = merged_features_nonpls,
          data_pls         = merged_features_pls,
          metadata         = merged_metadata,
          merged           = TRUE,
          n_samples_before = n_samples_before,
          n_samples_after  = n_samples_after
        ))

      }, error = function(e) {
        warning("Error merging technical replicates: ", e$message, ". Returning original data.")
        return(list(merged = FALSE, reason = e$message))
      })
    }

    if (merge_replicates) {
      msg("Attempting to merge technical replicates...")

      merge_result <- merge_technical_replicates(
        data_nonpls = listPreprocessed$data_scaledNONPLS_varFiltered,
        data_pls    = listPreprocessed$data_scaledPLS_varFiltered,
        metadata    = listPreprocessed$Metadata,
        verbose     = verbose
      )

      if (merge_result$merged) {
        listPreprocessed$data_scaledNONPLS_merged <- merge_result$data_nonpls
        listPreprocessed$data_scaledPLS_merged    <- merge_result$data_pls
        listPreprocessed$Metadata_merged          <- merge_result$metadata

        listPreprocessed$Dimensions               <- rbind(listPreprocessed$Dimensions,
                                                           data.frame(Step = "After replicate merging (NONPLS)",
                                                                      Samples = nrow(merge_result$data_nonpls), Features = ncol(merge_result$data_nonpls)))
        listPreprocessed$Dimensions               <- rbind(listPreprocessed$Dimensions,
                                                           data.frame(Step = "After replicate merging (PLS)",
                                                                      Samples = nrow(merge_result$data_pls),    Features = ncol(merge_result$data_pls)))

        listPreprocessed$ReplicateMerging <- list(
          merged = TRUE,
          n_samples_before = merge_result$n_samples_before,
          n_samples_after  = merge_result$n_samples_after,
          reduction_count  = merge_result$n_samples_before - merge_result$n_samples_after
        )

      } else {
        listPreprocessed$ReplicateMerging <- list(merged = FALSE, reason = merge_result$reason)
      }
    } else {
      msg("Automatic replicate merging disabled.")
      listPreprocessed$ReplicateMerging <- list(
        merged = FALSE,
        reason = "Disabled by user (merge_replicates = FALSE)"
      )
    }

    # ============================================================================
    # PROCESSING SUMMARY
    # ============================================================================

    # Helper for NULL-safe length
    `%||%` <- function(a, b) if (is.null(a)) b else a

    listPreprocessed$ProcessingSummary <- list(
      total_samples_processed        = nrow(listPreprocessed$Metadata),
      total_features_original        = length(listPreprocessed$All_Features_Metabolites),
      total_features_final_NONPLS    = ncol(listPreprocessed$data_scaledNONPLS_varFiltered),
      total_features_final_PLS       = ncol(listPreprocessed$data_scaledPLS_varFiltered),
      outliers_removed               = length(listPreprocessed$outliers_removed %||% character(0)),
      missing_values_imputed         = listPreprocessed$n_missing,
      drift_batch_correction_applied = driftBatchCorrection,
      uncorrected_features_count     = length(uncorrected_features),
      uncorrected_features_removed   = features_removed_uncorrected,
      normalization_method           = dataNormalize,
      transformation_method          = dataTransform,
      scaling_method_NONPLS          = dataScaleNONPLS,
      scaling_method_PLS             = dataScalePLS,
      filtering_reference            = filterReference,
      rsd_filtering_applied          = !is.null(filterMaxRSD),
      variance_filtering_applied     = !is.null(filterLowVariability),
      replicates_merged              = listPreprocessed$ReplicateMerging$merged
    )

    # Final summary messages
    if (merge_replicates && listPreprocessed$ReplicateMerging$merged) {
      msg(sprintf("Final merged dataset dimensions: %d samples X %d features (NONPLS-ready)",
                  nrow(listPreprocessed$data_scaledNONPLS_merged), ncol(listPreprocessed$data_scaledNONPLS_merged)))
      msg(sprintf("Final merged dataset dimensions: %d samples X %d features (PLS-ready)",
                  nrow(listPreprocessed$data_scaledPLS_merged), ncol(listPreprocessed$data_scaledPLS_merged)))
    } else {
      msg(sprintf("Final dataset dimensions: %d samples X %d features (NONPLS-ready)",
                  nrow(listPreprocessed$data_scaledNONPLS_varFiltered), ncol(listPreprocessed$data_scaledNONPLS_varFiltered)))
      msg(sprintf("Final dataset dimensions: %d samples X %d features (PLS-ready)",
                  nrow(listPreprocessed$data_scaledPLS_varFiltered), ncol(listPreprocessed$data_scaledPLS_varFiltered)))
    }

    msg("Data preprocessing pipeline completed successfully!")

  }, error = function(e) {

    # --- CATCH ERROR ---
    msg("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    msg(sprintf("ERROR during preprocessing: %s", e$message))
    msg("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    listPreprocessed$Error <- e$message

  }, finally = {

    # Calculate total preprocessing time in seconds
    listPreprocessed$ProcessingTimestampEnd         <- Sys.time()
    listPreprocessed$TotalProcessingTime_in_seconds <- as.numeric(
      difftime(listPreprocessed$ProcessingTimestampEnd, listPreprocessed$ProcessingTimestampStart, units = "secs")
    )

    if (!is.null(listPreprocessed$Error)) {
      msg(sprintf("Pipeline FAILED after %g seconds.", listPreprocessed$TotalProcessingTime_in_seconds))
    } else {
      msg(sprintf("Pipeline finished in %g seconds.", listPreprocessed$TotalProcessingTime_in_seconds))
    }
  })

  class(listPreprocessed) <- c("perform_PreprocessingPeakData", "list")
  return(listPreprocessed)
}

# ============================================================================
# PLOTTING FUNCTION
# ============================================================================

#' Create Before and After Signal Drift and Batch Correction Plots
#'
#' @description
#' Creates visualization plots comparing data before and after drift/batch correction.
#'
#' @param preprocessed_data List. Results from the `perform_PreprocessingPeakData` function.
#' @param n_features Numeric. Number of random features to plot. Default is 6. Suggested to be an even number.
#' @param feature_names Vector. Specific feature names to plot.
#' @param seed Numeric. Random seed for reproducible feature selection. Default is 123.
#'
#' @returns A combined plot grob object (requires ggplot2, gridExtra, grid).
#'
#' @author John Lennon L. Calorio
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap scale_color_manual labs theme_minimal theme element_text
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.newpage grid.draw
#'
#' @export
plot_before_after_correction <- function(
    preprocessed_data,
    n_features = 6,
    feature_names = NULL,
    seed = 123
) {

  # Input validation
  if (!is.list(preprocessed_data) ||
      !grepl("perform_PreprocessingPeakData", preprocessed_data$FunctionOrigin)) {
    stop("'preprocessed_data' must be from 'perform_PreprocessingPeakData' function.")
  }

  if (!preprocessed_data$Parameters$driftBatchCorrection) {
    message("Cannot plot: 'driftBatchCorrection' was set to FALSE.")
    return(invisible(NULL))
  }

  if (!is.null(preprocessed_data$Error)) {
    message("Cannot plot: Preprocessing pipeline failed.")
    return(invisible(NULL))
  }

  required_fields <- c("data_no_NA", "data_driftBatchCorrected", "Metadata")
  missing_fields  <- setdiff(required_fields, names(preprocessed_data))
  if (length(missing_fields) > 0) {
    stop("Missing required fields: ", paste(missing_fields, collapse = ", "))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE) ||
      !requireNamespace("grid", quietly = TRUE)) {
    stop("Required packages not available. Install: ggplot2, gridExtra, grid")
  }

  df_before <- preprocessed_data$data_no_NA
  df_after  <- preprocessed_data$data_driftBatchCorrected
  metadata  <- preprocessed_data$Metadata

  if (identical(df_before, df_after)) {
    warning("Data before and after correction are identical. No correction may have been applied.")
  }

  # Align data (in case outliers were removed, etc.)
  common_samples <- Reduce(intersect, list(rownames(df_before), rownames(df_after), metadata$Samples))

  df_before <- df_before[common_samples, , drop = FALSE]
  df_after  <- df_after[common_samples, , drop = FALSE]
  metadata  <- metadata[match(common_samples, metadata$Samples), , drop = FALSE]

  # Feature selection
  if (!is.null(feature_names)) {
    available_features        <- intersect(feature_names, colnames(df_before))
    if (length(available_features) == 0) {
      stop("None of the specified feature names found in the data.")
    }
    selected_features         <- available_features
  } else {
    n_features                <- min(n_features, ncol(df_before), ncol(df_after))
    if (n_features < 1) stop("No features available for plotting.")

    # Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)

    selected_features_indices <- sample(ncol(df_before), n_features, replace = FALSE)
    selected_features         <- colnames(df_before)[selected_features_indices]
  }

  # Create plots
  plot_list <- lapply(selected_features, function(feature_name) {

    # Vectorized log transformation
    before_vals <- log10(pmax(df_before[, feature_name], 1e-10))
    after_vals  <- log10(pmax(df_after[, feature_name], 1e-10))

    # Create long format data without tidyr
    plot_data_long <- data.frame(
      sample_id        = rep(seq_len(nrow(df_before)), 2),
      injection        = rep(metadata$InjectionSequence, 2),
      intensity        = c(before_vals, after_vals),
      correction       = rep(c("before", "after"), each = nrow(df_before)),
      class            = rep(metadata$Group_, 2),
      batch            = rep(as.factor(metadata$Batches), 2),
      stringsAsFactors = FALSE
    )

    plot_data_long$correction <- factor(plot_data_long$correction, levels = c("before", "after"))

    ggplot2::ggplot(plot_data_long, ggplot2::aes(x = injection, y = intensity)) +
      ggplot2::geom_point(ggplot2::aes(color = class, shape = batch), size = 2, alpha = 0.7) +
      ggplot2::geom_smooth(data = subset(plot_data_long, class == "QC"),
                           formula = y ~ x, method = "loess", se = TRUE, alpha = 0.3,
                           color = "darkred", linetype = "dashed", fullrange = TRUE) +
      ggplot2::facet_wrap(~ correction, scales = "free_y",
                          labeller = ggplot2::labeller(
                            correction = c(before = "Before Correction", after = "After Correction"))) +
      ggplot2::scale_color_manual(values = c("QC" = "red", "Sample" = "blue", "Blank" = "green"),
                                  limits = c("QC", "Sample", "Blank"), drop = FALSE) +
      ggplot2::labs(title = paste("Feature:", feature_name),
                    x = "Injection Sequence", y = "log10(Intensity)",
                    color = "Class", shape = "Batch") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title   = ggplot2::element_text(size = 12, hjust = 0.5),
        axis.text    = ggplot2::element_text(size = 10),
        axis.title   = ggplot2::element_text(size = 11),
        legend.title = ggplot2::element_text(size = 10),
        legend.text  = ggplot2::element_text(size = 9),
        strip.text   = ggplot2::element_text(size = 10, face = "bold")
      )
  })

  n_batches <- length(unique(metadata$Batches))
  main_title <- ifelse(n_batches <= 1,
                       "Features Before and After Signal Drift-Correction",
                       "Features Before and After Signal Drift and Batch Correction")

  # Determine number of columns for plot layout
  plot_cols <- if (length(plot_list) == 1) 1 else if (length(plot_list) <= 6) 2 else 3

  # Use arrangeGrob to create the plot object without displaying it
  combined_plot <- gridExtra::arrangeGrob(
    grobs = plot_list,
    ncol  = plot_cols,
    top   = main_title
  )

  # Display the plot
  grid::grid.newpage()
  grid::grid.draw(combined_plot)

  message("Before/after correction plots created successfully.")

  # Return invisibly so it can be stored and redrawn later
  return(invisible(combined_plot))
}

# S3 Methods
#' @export
print.perform_PreprocessingPeakData <- function(x, ...) {
  cat("=== Metabolomics Preprocessing Pipeline ===\n")

  # Determine which final dataset is primary (Merged or VarFiltered)
  if (isTRUE(x$Parameters$merge_replicates) && x$ReplicateMerging$merged) {
    final_dim <- dim(x$data_scaledNONPLS_merged)
    suffix <- "(Merged)"
  } else {
    final_dim <- dim(x$data_scaledNONPLS_varFiltered)
    suffix <- "(Filtered)"
  }

  cat("Original Data: ", x$Dimensions$Samples[1], "samples x", x$Dimensions$Features[1], "features\n")
  cat("Final Data:    ", final_dim[1], "samples x", final_dim[2], "features", suffix, "\n")
  cat("Normalization: ", x$Parameters$dataNormalize, "\n")
  cat("Transform:     ", x$Parameters$dataTransform, "\n")
  cat("Scaling:       ", x$Parameters$dataScaleNONPLS, "(NON-PLS) /", x$Parameters$dataScalePLS, "(PLS)\n")
  invisible(x)
}

#' @export
summary.perform_PreprocessingPeakData <- function(object, ...) {
  ans <- list(
    params = object$Parameters,
    dim_history = object$Dimensions,
    missing_imputed = object$n_missing,
    uncorrected = object$ProcessingSummary$uncorrected_features_count,
    replicates = object$ReplicateMerging
  )
  class(ans) <- "summary.perform_PreprocessingPeakData"
  return(ans)
}

#' @export
print.summary.perform_PreprocessingPeakData <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("Preprocessing Step-by-Step Summary\n")
  cat("---------------------------------------\n")
  print(x$dim_history, row.names = FALSE)

  cat("\n-- Key Statistics --\n")
  cat("Missing Values Imputed:   ", x$missing_imputed, "\n")
  if(x$params$driftBatchCorrection) {
    cat("Uncorrected Features:     ", x$uncorrected, "\n")
  }

  if (x$replicates$merged) {
    cat("\n-- Technical Replicates --\n")
    cat("Status: Merged\n")
    cat("Reduction:", x$replicates$reduction_count, "samples averaged\n")
  }
  invisible(x)
}
