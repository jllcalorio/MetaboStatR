#' Perform Comprehensive Data Quality Check for Metabolomics Data
#'
#' @description
#' This function imports and validates metabolomics data from Excel, CSV, TSV, or text files,
#' performing comprehensive quality controls to ensure data integrity and compatibility with
#' downstream analysis functions. The function validates data structure, checks for required
#' metadata rows, ensures uniqueness constraints, cleans special characters from identifiers,
#' and prepares data for preprocessing pipelines.
#'
#' Key validation checks include:
#' \itemize{
#'   \item Presence of required metadata rows (Sample, SubjectID, Replicate, Group, Batch, Injection, Normalization, Response)
#'   \item Uniqueness of sample names, injection sequences, and feature/metabolite identifiers
#'   \item Proper QC sample annotation (empty values in SubjectID, Replicate, Normalization, Response for QC samples)
#'   \item Numeric validation of feature/metabolite data
#'   \item Character cleaning and standardization of identifiers
#' }
#'
#' @param file_location Character string specifying the file path. If NULL, an interactive
#'   file selection dialog will open. Supported formats: .xlsx, .csv, .tsv, .txt
#' @param sheet_name Character string specifying the Excel worksheet name. Ignored for
#'   non-Excel files. If NULL for Excel files, the first sheet is used.
#' @param skip_rows Integer specifying the number of rows to skip when reading the file.
#'   Default is 0.
#' @param separator Character string specifying the field separator for delimited files.
#' Common values: "," (comma), "\\t" (tab). Default is ",". Ignored for Excel files.
#' @param validate_qc Logical indicating whether to enforce QC validation rules.
#'   Default is TRUE.
#' @param allow_missing_optional Logical indicating whether to allow missing values in
#'   optional metadata rows (SubjectID, Replicate, Normalization, Response). Default is TRUE.
#' @param clean_names Logical indicating whether to clean special characters from names.
#'   Default is TRUE.
#' @param verbose Logical indicating whether to display progress messages. Default is TRUE.
#'
#' @return A list containing:
#'   \item{raw_data}{Data frame with the original data as loaded from the file}
#'   \item{quality_checked_data}{Data frame with validated and cleaned data, sorted by injection sequence}
#'   \item{metadata_summary}{Summary statistics of the metadata}
#'   \item{validation_report}{Detailed validation results}
#'   \item{file_info}{Information about the source file}
#'   \item{processing_log}{Log of all processing steps performed}
#'
#' @details
#' The input data must follow a specific structure:
#' \itemize{
#'   \item Row 1: "Sample" - Unique sample identifiers (no spaces recommended)
#'   \item Row 2: "SubjectID" - Numeric subject identifiers (can be non-unique)
#'   \item Row 3: "Replicate" - Replicate identifiers (can be non-unique)
#'   \item Row 4: "Group" - Group assignments including QC samples
#'   \item Row 5: "Batch" - Batch numbers
#'   \item Row 6: "Injection" - Unique injection sequence numbers
#'   \item Row 7: "Normalization" - Concentration markers (e.g., osmolality)
#'   \item Row 8: "Response" - Response variable values
#'   \item Rows 9+: Feature/metabolite data (e.g., m/z@retention_time format)
#' }
#'
#' Missing values should be left blank or encoded as 0. QC samples in the Group row
#' must have empty values in SubjectID, Replicate, Normalization, and Response rows.
#'
#' @examples
#' \dontrun{
#' # Basic usage with file selection dialog
#' result <- perform_DataQualityCheck()
#'
#' # Specify file location directly
#' result <- perform_DataQualityCheck(
#'   file_location = "path/to/metabolomics_data.xlsx",
#'   sheet_name = "Sheet1"
#' )
#'
#' # CSV file with custom separator
#' result <- perform_DataQualityCheck(
#'   file_location = "path/to/data.csv",
#'   separator = ";",
#'   skip_rows = 1
#' )
#'
#' # Access results
#' clean_data <- result$raw_data
#' validation_summary <- result$validation_report
#' }
#'
#' @seealso
#' \code{\link{perform_PreprocessingPeakData}} for the next step in the analysis pipeline
#'
#' @export
#' @importFrom readxl read_excel
#' @importFrom dplyr mutate arrange select all_of
#' @importFrom stringr str_replace_all str_trim
#' @importFrom utils read.csv read.delim
#' @importFrom tools file_ext
#' @importFrom stats setNames
perform_DataQualityCheck <- function(
    file_location = NULL,
    sheet_name = NULL,
    skip_rows = 0,
    separator = ",",
    validate_qc = TRUE,
    allow_missing_optional = TRUE,
    clean_names = TRUE,
    verbose = TRUE
) {

  # Initialize results object with processing metadata
  results <- initialize_results_object(file_location, sheet_name, validate_qc,
                                       allow_missing_optional, clean_names)

  tryCatch({
    # Step 1: Load and validate file
    if (verbose) message("Step 1/6: Reading and validating file...")
    original_data <- load_and_validate_file(file_location, sheet_name, skip_rows,
                                            separator, verbose)
    results$raw_data <- original_data  # Store original data
    results$file_info <- extract_file_info(file_location, original_data)

    # Create working copy for quality checking
    raw_data <- original_data

    # Step 2: Clean names if requested
    if (clean_names) {
      if (verbose) message("Step 2/6: Cleaning special characters in identifiers...")
      raw_data <- clean_data_identifiers(raw_data, verbose)
    } else {
      if (verbose) message("Step 2/6: Skipping name cleaning (clean_names = FALSE)")
    }

    # Step 3: Validate required structure
    if (verbose) message("Step 3/6: Validating data structure...")
    structure_validation <- validate_data_structure(raw_data, allow_missing_optional)
    results$validation_report$structure <- structure_validation

    # Step 4: Check for duplicates
    if (verbose) message("Step 4/6: Checking for duplicate identifiers...")
    raw_data <- transpose_and_prepare_data(raw_data)
    duplicate_validation <- validate_uniqueness_constraints(raw_data)
    results$validation_report$duplicates <- duplicate_validation

    # Step 5: Validate QC rules if requested
    if (validate_qc) {
      if (verbose) message("Step 5/6: Validating QC sample rules...")
      qc_validation <- validate_qc_rules(raw_data)
      results$validation_report$qc_rules <- qc_validation
    } else {
      if (verbose) message("Step 5/6: Skipping QC validation (validate_qc = FALSE)")
    }

    # Step 6: Validate numeric data and finalize
    if (verbose) message("Step 6/6: Validating numeric data and finalizing...")
    raw_data <- sort_by_injection_sequence(raw_data)
    numeric_validation <- validate_numeric_features(raw_data)
    results$validation_report$numeric_data <- numeric_validation

    # Finalize results
    results$quality_checked_data <- finalize_data_structure(raw_data)
    results$metadata_summary <- generate_metadata_summary(results$quality_checked_data)
    results$processing_log$completion_time <- Sys.time()
    results$processing_log$success <- TRUE

    if (verbose) message("✓ Data quality check completed successfully!")

  }, error = function(e) {
    results$processing_log$error <- e$message
    results$processing_log$success <- FALSE
    results$processing_log$completion_time <- Sys.time()

    stop("Data quality check failed: ", e$message, call. = FALSE)
  })

  return(results)
}

# ================================================================================
# HELPER FUNCTIONS
# ================================================================================

#' Initialize Results Object
#' @keywords internal
initialize_results_object <- function(file_location, sheet_name, validate_qc,
                                      allow_missing_optional, clean_names) {
  list(
    FunctionOrigin = "perform_DataQualityCheck",
    version = "2.0.0",
    file_location = file_location,
    sheet_name = sheet_name,
    validate_qc = validate_qc,
    allow_missing_optional = allow_missing_optional,
    clean_names = clean_names,
    processing_log = list(
      start_time = Sys.time(),
      completion_time = NULL,
      success = FALSE,
      error = NULL
    ),
    file_info = list(),
    validation_report = list(),
    metadata_summary = list(),
    raw_data = NULL,
    quality_checked_data = NULL
  )
}

#' Load and Validate File
#' @keywords internal
load_and_validate_file <- function(file_location, sheet_name, skip_rows,
                                   separator, verbose) {

  # Handle file selection
  if (is.null(file_location)) {
    file_location <- select_file_interactively()
  }

  # Validate file existence
  if (!file.exists(file_location)) {
    stop("File does not exist: ", file_location, call. = FALSE)
  }

  # Validate file size (optional safety check)
  file_size_mb <- file.size(file_location) / (1024^2)
  if (file_size_mb > 500) {  # 500MB threshold
    warning("Large file detected (", round(file_size_mb, 1), "MB). Processing may be slow.")
  }

  # Load data based on file type
  file_ext <- tolower(tools::file_ext(file_location))

  raw_data <- suppressMessages(
    switch(file_ext,
           "xlsx" = load_excel_file(file_location, sheet_name, skip_rows),
           "csv"  = load_csv_file(file_location, separator, skip_rows),
           "txt"  = load_delimited_file(file_location, separator, skip_rows),
           "tsv"  = load_delimited_file(file_location, "\t", skip_rows),
           stop("Unsupported file type: '", file_ext, "'. Supported types: xlsx, csv, tsv, txt",
                call. = FALSE)
    )
  )

  # Validate loaded data
  if (is.null(raw_data) || nrow(raw_data) == 0 || ncol(raw_data) == 0) {
    stop("Failed to load data or file is empty", call. = FALSE)
  }

  if (verbose) {
    message("✓ File loaded successfully: ", nrow(raw_data), " rows × ", ncol(raw_data), " columns")
  }

  return(raw_data)
}

#' Select File Interactively
#' @keywords internal
select_file_interactively <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(rstudioapi::selectFile(
      caption = "Select metabolomics data file",
      filter = "Data Files (*.xlsx, *.csv, *.tsv, *.txt)"
    ))
  } else {
    return(file.choose())
  }
}

#' Load Excel File
#' @keywords internal
load_excel_file <- function(file_location, sheet_name, skip_rows) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required for Excel files but not available", call. = FALSE)
  }

  tryCatch({
    readxl::read_excel(
      path = file_location,
      sheet = sheet_name,
      col_names = FALSE,
      skip = skip_rows,
      .name_repair = "minimal"
    )
  }, error = function(e) {
    stop("Failed to read Excel file: ", e$message, call. = FALSE)
  })
}

#' Load CSV File
#' @keywords internal
load_csv_file <- function(file_location, separator, skip_rows) {
  tryCatch({
    utils::read.csv(
      file = file_location,
      header = FALSE,
      sep = separator,
      skip = skip_rows,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }, error = function(e) {
    stop("Failed to read CSV file: ", e$message, call. = FALSE)
  })
}

#' Load Delimited File
#' @keywords internal
load_delimited_file <- function(file_location, separator, skip_rows) {
  tryCatch({
    utils::read.delim(
      file = file_location,
      header = FALSE,
      sep = separator,
      skip = skip_rows,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }, error = function(e) {
    stop("Failed to read delimited file: ", e$message, call. = FALSE)
  })
}

#' Extract File Information
#' @keywords internal
extract_file_info <- function(file_location, raw_data) {
  list(
    file_path = basename(file_location),  # Just the filename, not full path
    file_name = basename(file_location),
    file_size_mb = round(file.size(file_location) / (1024^2), 2),
    file_extension = tolower(tools::file_ext(file_location)),
    dimensions = list(rows = nrow(raw_data), columns = ncol(raw_data)),
    modification_time = file.mtime(file_location)
  )
}

#' Clean Data Identifiers
#' @keywords internal
clean_data_identifiers <- function(raw_data, verbose) {

  # Optimized text cleaning function
  clean_text_vectorized <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "[^a-zA-Z0-9@._]", "_")
    x <- stringr::str_replace_all(x, "_{2,}", "_")
    x <- stringr::str_replace_all(x, "^_+|_+$", "")
    x <- stringr::str_trim(x)
    return(x)
  }

  # Clean first column (row identifiers)
  raw_data[[1]] <- clean_text_vectorized(raw_data[[1]])

  # Clean specific rows (Sample and Replicate only - skip Group row)
  target_rows <- c(1, 3)  # Sample and Replicate rows only
  for (row_idx in target_rows) {
    if (row_idx <= nrow(raw_data)) {
      raw_data[row_idx, ] <- lapply(raw_data[row_idx, ], clean_text_vectorized)
    }
  }

  if (verbose) {
    message("✓ Special characters cleaned from identifiers (excluding Group row)")
  }

  return(raw_data)
}

#' Validate Data Structure
#' @keywords internal
validate_data_structure <- function(raw_data, allow_missing_optional) {
  required_headers <- c("Sample", "SubjectID", "Replicate", "Group",
                        "Batch", "Injection", "Normalization", "Response")

  # Check minimum dimensions
  if (nrow(raw_data) < 8) {
    stop("Data must have at least 8 rows (required metadata rows)", call. = FALSE)
  }

  if (ncol(raw_data) < 2) {
    stop("Data must have at least 2 columns", call. = FALSE)
  }

  # Extract first column values
  first_col_values <- as.character(unlist(raw_data[1:8, 1]))

  # Check for required headers
  missing_headers <- required_headers[!required_headers %in% first_col_values]

  if (length(missing_headers) > 0) {
    stop("Missing required metadata rows: ", paste(missing_headers, collapse = ", "),
         call. = FALSE)
  }

  # Validate header positions
  header_positions <- match(required_headers, first_col_values)
  expected_positions <- 1:8

  if (!all(header_positions == expected_positions)) {
    stop("Metadata rows are not in the correct order. Expected order: ",
         paste(required_headers, collapse = ", "), call. = FALSE)
  }

  return(list(
    required_headers_present = TRUE,
    header_positions = setNames(header_positions, required_headers),
    total_features = ncol(raw_data) - 1,
    total_samples = nrow(raw_data) - 8
  ))
}

#' Transpose and Prepare Data
#' @keywords internal
transpose_and_prepare_data <- function(raw_data) {
  # Transpose so that samples are rows and features are columns
  raw_data_t <- as.data.frame(t(raw_data))

  # Set column names from first row and remove it
  colnames(raw_data_t) <- as.character(raw_data_t[1, ])
  raw_data_t <- raw_data_t[-1, ]

  return(raw_data_t)
}

#' Validate Uniqueness Constraints
#' @keywords internal
validate_uniqueness_constraints <- function(raw_data) {

  validation_results <- list()

  # Check sample names uniqueness
  sample_names <- as.character(raw_data$Sample)
  sample_duplicates <- find_duplicates(sample_names)
  if (length(sample_duplicates) > 0) {
    stop("Duplicate sample names found: ", paste(sample_duplicates, collapse = ", "),
         call. = FALSE)
  }
  validation_results$sample_names <- "✓ All sample names are unique"

  # Check injection sequence uniqueness
  injection_seq <- as.numeric(raw_data$Injection)
  if (any(is.na(injection_seq))) {
    stop("Non-numeric values found in Injection sequence", call. = FALSE)
  }

  injection_duplicates <- find_duplicates(injection_seq)
  if (length(injection_duplicates) > 0) {
    stop("Duplicate injection sequences found: ", paste(injection_duplicates, collapse = ", "),
         call. = FALSE)
  }
  validation_results$injection_sequence <- "✓ All injection sequences are unique"

  # Check feature names uniqueness
  required_headers <- c("Sample", "SubjectID", "Replicate", "Group",
                        "Batch", "Injection", "Normalization", "Response")
  feature_names <- colnames(raw_data)[!colnames(raw_data) %in% required_headers]
  feature_duplicates <- find_duplicates(feature_names)
  if (length(feature_duplicates) > 0) {
    stop("Duplicate feature/metabolite names found: ", paste(feature_duplicates, collapse = ", "),
         call. = FALSE)
  }
  validation_results$feature_names <- "✓ All feature/metabolite names are unique"

  return(validation_results)
}

#' Find Duplicates in Vector
#' @keywords internal
find_duplicates <- function(x) {
  unique(x[duplicated(x) | duplicated(x, fromLast = TRUE)])
}

#' Sort Data by Injection Sequence
#' @keywords internal
sort_by_injection_sequence <- function(raw_data) {
  # Convert injection to numeric and sort
  raw_data$Injection <- as.numeric(raw_data$Injection)
  raw_data <- raw_data[order(raw_data$Injection), ]
  return(raw_data)
}

#' Validate QC Rules
#' @keywords internal
validate_qc_rules <- function(raw_data) {

  # Transpose back for row-based operations
  raw_data_t <- as.data.frame(t(raw_data))

  # Find QC columns (case insensitive)
  group_row <- as.character(raw_data_t[4, ])  # Group row
  qc_cols <- which(grepl("SQC|EQC|QC", group_row, ignore.case = TRUE))

  if (length(qc_cols) == 0) {
    return(list(message = "No QC columns found", qc_columns = 0))
  }

  # Check rows that should be empty for QC samples
  rows_to_check <- c(2, 3, 7, 8)  # SubjectID, Replicate, Normalization, Response
  violations <- c()

  for (col in qc_cols) {
    for (row in rows_to_check) {
      if (nrow(raw_data_t) >= row && ncol(raw_data_t) >= col) {
        cell_value <- raw_data_t[row, col]
        if (!is.na(cell_value) && cell_value != "" && cell_value != "0") {
          violations <- c(violations, paste0("Row ", row, " Col ", col))
        }
      }
    }
  }

  if (length(violations) > 0) {
    stop("QC rule violations: QC samples must have empty values in SubjectID, Replicate, Normalization, and Response rows. Violations found at: ",
         paste(violations, collapse = ", "), call. = FALSE)
  }

  return(list(
    message = "✓ QC rules validated successfully",
    qc_columns = length(qc_cols),
    qc_samples = sum(grepl("SQC|EQC|QC", raw_data$Group, ignore.case = TRUE))
  ))
}

#' Validate Numeric Features
#' @keywords internal
validate_numeric_features <- function(raw_data) {

  required_headers <- c("Sample", "SubjectID", "Replicate", "Group",
                        "Batch", "Injection", "Normalization", "Response")

  # Extract feature columns
  feature_cols <- raw_data[, !colnames(raw_data) %in% required_headers, drop = FALSE]

  if (ncol(feature_cols) == 0) {
    stop("No feature/metabolite columns found", call. = FALSE)
  }

  # Test numeric conversion
  conversion_issues <- c()

  for (col_name in colnames(feature_cols)) {
    col_data <- feature_cols[[col_name]]

    # Convert to numeric, checking for issues
    suppressWarnings({
      numeric_data <- as.numeric(col_data)
    })

    # Check for conversion issues (excluding originally missing values)
    original_na <- is.na(col_data) | col_data == "" | col_data == "0"
    conversion_na <- is.na(numeric_data)

    problematic <- !original_na & conversion_na

    if (any(problematic)) {
      problematic_values <- unique(col_data[problematic])
      conversion_issues <- c(conversion_issues,
                             paste0(col_name, ": ", paste(problematic_values, collapse = ", ")))
    }
  }

  if (length(conversion_issues) > 0) {
    stop("Non-numeric values found in feature columns:\n",
         paste(conversion_issues, collapse = "\n"), call. = FALSE)
  }

  return(list(
    message = "✓ All feature data is numeric or convertible to numeric",
    feature_count = ncol(feature_cols),
    sample_count = nrow(feature_cols)
  ))
}

#' Finalize Data Structure
#' @keywords internal
finalize_data_structure <- function(raw_data) {

  # Transpose back to Features x Samples format for output
  raw_data_final <- t(raw_data)
  raw_data_final <- cbind(rownames(raw_data_final), raw_data_final)
  rownames(raw_data_final) <- NULL
  raw_data_final <- as.data.frame(raw_data_final)

  return(raw_data_final)
}

#' Generate Metadata Summary
#' @keywords internal
generate_metadata_summary <- function(raw_data) {

  # Extract metadata from transposed format
  raw_data_samples <- as.data.frame(t(raw_data))
  colnames(raw_data_samples) <- raw_data_samples[1, ]
  raw_data_samples <- raw_data_samples[-1, ]

  summary_stats <- list(
    total_samples = nrow(raw_data_samples),
    total_features = ncol(raw_data_samples) - 8,
    groups = table(raw_data_samples$Group),
    batches = length(unique(raw_data_samples$Batch)),
    qc_samples = sum(grepl("QC", raw_data_samples$Group, ignore.case = TRUE)),
    injection_range = range(as.numeric(raw_data_samples$Injection), na.rm = TRUE),
    missing_data_summary = list(
      subjectid = sum(is.na(raw_data_samples$SubjectID) | raw_data_samples$SubjectID == ""),
      replicate = sum(is.na(raw_data_samples$Replicate) | raw_data_samples$Replicate == ""),
      normalization = sum(is.na(raw_data_samples$Normalization) | raw_data_samples$Normalization == ""),
      response = sum(is.na(raw_data_samples$Response) | raw_data_samples$Response == "")
    )
  )

  return(summary_stats)
}






#' #' Perform Data Quality Check
#' #'
#' #' @description
#' #' This function imports Excel or text files and performs quality controls prior to usage in the
#' #' succeeding functions. This is ideally used to check data requirements such as uniqueness of
#' #' inputs in the "Sample" row, uniqueness of Features/Metabolites, data types, and the
#' #' overall characteristics of the data. This also cleans the sample names, replicates, groups,
#' #' and feature/metabolite names from its special characters and replaces with underscores, except @ and . (period).
#' #' This is done to avoid errors later on during manual data extraction, plotting, etc, where R might
#' #' encounter an error from column/row names.
#' #'
#' #' @param file_location String. A string of the location of the data in the computer. Must have the following characteristics. Other files may be loaded as long as they can be read by `utils::read.delim`.
#' #'   \itemize{
#' #'     \item 1st row name is "Sample": The sample names, ideally shortened versions to make it easier to view in visualizations; without spaces.
#' #'     \item 2nd row name is "SubjectID": Must be numeric to be sorted correctly. Can contain non-unique identifiers of "Sample" names, usually used when "Replicate" is present. Also referred to as the individuals in the data.
#' #'     \item 3rd row name is "Replicate": The replicate IDs of Samples. Ideally non-unique identifiers.
#' #'     \item 4th row name is "Group": Contains QC and the groups (with/out spaces).
#' #'     \item 5th row name is "Batch": The batch numbers.
#' #'     \item 6th row name is "Injection": The injection sequence.
#' #'     \item 7th row name is "Normalization": Concentration markers to be used during data normalization to correct for variations in e.g., urine, concentration between samples. Can be osmolality values, specific gravity, etc.
#' #'     \item 8th row name is "Response": Numeric. Ideally used when the the response variable (dependent variable) to be used is not "Group".
#' #'     \item Others requirement/s are:
#' #'     \itemize{
#' #'       \item The 1st 8 rows above must be present.
#' #'       \item After the 8th row (Response), are the features possibly in a mass-to-charge ratio (m/z) and retention time format. Example `199.101@0.111` assuming you rounded off the decimal places to 3 decimal places (suggested for a cleaner output in visualizations).
#' #'       \item If the Normalization values, SubjectID, Replicate, and Response are not present, leave the rows blank except for the words in the 1st column (e.g., "Normalization' ).
#' #'       \item Missing values must be left blank or encoded as 0. Blanks and 0s are identified as NAs or missing values.
#' #'       }
#' #'     }
#' #' @param sheet_name String. The name of the worksheet where the data is located. Ignored when data type is ".csv".
#' #' @param skip_rows Numeric. The number of rows in the data to skip reading.
#' #' @param separator String. The separator in the ".csv" or ".txt" file. If this fails, it's likely that the file is supposed to be delimited, which should be saved as ".csv", ".tsv", or ".txt".
#' #'   \itemize{
#' #'     \item ",": A comma. From a comma-separated file (csv).
#' #'     \item "/t": A tab From a tab-separated file (.txt or .tsv).
#' #'     \item Also accepts other delimiters, as long as it is consistent in your data.
#' #'     }
#' #'     Defaults to ",". Ignored when the file is ".xlsx".
#' #' @returns A data frame ready to be used in the `perform_PreprocessingPeakData` function in the `MetaboStatR` package. The data frame is sorted in ascending order according to the injection sequence, as required during the "drift- and batch-correction" stage.
#' #' @export
#' #'
#' #' @examples
#' #' \dontrun{
#' #' # Using the default parameters
#' #' # The location typically looks like
#' #' # "C:/Users/Documents/Metabolomics Data Analysis Folder/file_name_of_data.xlsx"
#' #' perform_DataQualityCheck(
#' #'  file_location = the_complete_location_of_your_data
#' #' )
#' #' }
#' #'
#' perform_DataQualityCheck <- function(
#'     file_location = NULL,
#'     sheet_name    = NULL,
#'     skip_rows     = 0,
#'     separator     = ","
#' ) {
#'
#'   qualityCheckResults                <- base::list()
#'   qualityCheckResults$FunctionOrigin <- "perform_DataQualityCheck"
#'   qualityCheckResults$FileLocation   <- file_location
#'   qualityCheckResults$SheetName      <- sheet_name
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # ----------------------------------------------- Load the data -----------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   message("Reading the file...")
#'
#'   # If file_location is NULL, then choose a file.
#'   if (base::is.null(file_location)) {
#'     if (base::requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
#'       file_location <- rstudioapi::selectFile(
#'         caption = "Select a file"
#'         # ,filter = "CSV Files (*.csv)"
#'       )
#'     } else {
#'       # This part ideally should not be used, but just in case the above fails.
#'       file_location <- base::file.choose()
#'     }
#'   }
#'
#'   # Check if file exists
#'   if (!base::file.exists(file_location)) {
#'     stop("File does not exist: ", file_location)
#'   }
#'
#'   # Extract file extension
#'   file_ext <- base::tolower(tools::file_ext(file_location))
#'
#'   # Read data based on file type
#'   raw_data <- base::suppressMessages(
#'     switch(file_ext,
#'            "xlsx" = {
#'              # Read Excel file
#'              readxl::read_excel(
#'                file_location,
#'                sheet = sheet_name,
#'                range = NULL,
#'                col_names = FALSE,
#'                skip = skip_rows
#'              )
#'            },
#'            "csv" = {
#'              # Read CSV file
#'              utils::read.csv(
#'                file_location,
#'                header = FALSE,
#'                sep = separator,
#'                skip = skip_rows
#'              )
#'            },
#'            "txt" = {
#'              # Read text file (assuming tab-delimited by default)
#'              utils::read.delim(
#'                file_location,
#'                header = FALSE,
#'                sep = separator,
#'                skip = skip_rows
#'              )
#'            },
#'            "tsv" = {
#'              # Read text file (assuming tab-delimited by default)
#'              utils::read.delim(
#'                file_location,
#'                header = FALSE,
#'                sep = separator,
#'                skip = skip_rows
#'              )
#'            },
#'            {
#'              # Unsupported file type
#'              stop("Unsupported file type: ", file_ext,
#'                   ". Supported types are: xlsx, csv, tsv, txt")
#'            }
#'     )
#'   )
#'
#'   message(paste0("The file '", file_location, "' has been read and loaded successfully."))
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # ----------- Function to clean problematic characters in Sample names and Feature/Metabolite names -----------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   message("Cleaning the 'Sample', 'Replicate', 'Group', and 'Feature/Metabolite names' with its special characters. Removing colons, semi-colons, parentheses, brackets, braces, dashes, stars, etc...")
#'
#'   # Helper function to clean character strings
#'   clean_text <- function(x) {
#'     x <- base::as.character(x)
#'     x <- stringr::str_replace_all(x, "[^a-zA-Z0-9@._]", "_") # @, . (period), underscore. and numbers from 0-9
#'     x <- stringr::str_replace_all(x, " {2,}", " ") # multiple spaces
#'     x <- stringr::str_replace_all(x, " ", "_") # a space
#'     x <- stringr::str_trim(x) # trailing spaces
#'     x <- stringr::str_replace_all(x, "_{2,}", "_") # multiple underscores
#'     x <- stringr::str_replace_all(x, "^_+|_+$", "")
#'     return(x)
#'   }
#'
#'   # Clean values in a specific column
#'   clean_problematic_chars_col <- function(data, column = 1) {
#'     data[[column]] <- clean_text(data[[column]])
#'     return(data)
#'   }
#'
#'   # Clean values in a specific row
#'   clean_problematic_chars_row <- function(data, row = 1) {
#'     cleaned_row <- base::lapply(data[row, ], clean_text)  # keep as list to preserve structure
#'     data[row, ] <- cleaned_row
#'     return(data)
#'   }
#'
#'   # Apply cleaning
#'   raw_data <- raw_data %>%
#'     clean_problematic_chars_col(column = 1) %>%
#'     clean_problematic_chars_row(row    = 1) %>%
#'     clean_problematic_chars_row(row    = 3) %>%
#'     clean_problematic_chars_row(row    = 4)
#'
#'   message("Cleaning complete.")
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # ------------------------------------ Check if required rows are present -------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   message("Checking for required rows...")
#'
#'   # Feature/Metabolite names uniqueness, remove first not needed columns
#'   required_headers <- base::c("Sample", "SubjectID", "Replicate", "Group", "Batch", "Injection", "Normalization", "Response")
#'   if (!all(required_headers %in% base::as.vector(base::unlist(raw_data[, 1])))) {
#'     stop("Data preprocessing is halted. The raw data does not contain the expected rows (Sample, SubjectID, Replicate, Group, Batch, Injection, Normalization, Response). Check for typo or the 1st 8 rows in your file.")
#'   }
#'
#'   message("The required rows are present.")
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------- Check if there are duplicates in some rows ---------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # ------------------------------------- Function to check for duplicates --------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   check_duplicates <- function(data, data_name, action = "stop") {
#'     if (base::length(unique(data)) != base::length(data)) {
#'       duplicated_values <- base::unique(data[base::duplicated(data) | base::duplicated(data, fromLast = TRUE)])
#'       message_text_stop <- base::paste0("There are duplicate ", data_name, " in your data, which are required to be unique to proceed. Here are the duplicated ", data_name, ": ",
#'                                         base::paste(duplicated_values, collapse = ", "))
#'       message_text_warning <- base::paste0("There are duplicate ", data_name, " in your data, which you may want to be unique. Here are the duplicated ", data_name, ": ",
#'                                            base::paste(duplicated_values, collapse = ", "))
#'
#'       if (action == "stop") {
#'         stop(message_text_stop)
#'       } else if (action == "warning") {
#'         warning(message_text_warning)
#'       } else {
#'         stop("Invalid action parameter. Use 'stop' or 'warning'.")
#'       }
#'     }
#'   }
#'
#'   message("Checking for duplicate 'Sample names' and 'Features/metabolites names'...")
#'
#'   # Transpose the data
#'   raw_data <- base::as.data.frame(base::t(raw_data)) # The data is now Samples x Features (row x column)
#'
#'   # Apply checks for sample names and injection orders
#'   check_duplicates(base::as.vector(base::unlist(raw_data[, 1])), "Sample names",     action = "stop")
#'   check_duplicates(base::as.vector(base::unlist(raw_data[, 6])), "Injection orders", action = "stop")
#'
#'   # Apply checks for features/metabolites
#'   all_column_names       <- base::as.vector(base::unlist(raw_data[1, ])) # Get column names excluding required_headers and check for duplicates
#'   remaining_column_names <- all_column_names[!all_column_names %in% required_headers]
#'   check_duplicates(remaining_column_names, "Features/metabolites")
#'
#'   message("The 'Sample names' and the 'Features/metabolites' are unique, as required.")
#'
#'   # If all is good, proceed with the remaining tests
#'   raw_data <- raw_data %>%
#'     stats::setNames(base::as.character(.[1, ])) %>% .[-1, ] %>% # Set column names from 1st row (features, etc.)
#'     # Arrange by Injection Sequence
#'     dplyr::mutate(Injection = base::as.numeric(Injection)) %>% # Convert to numeric so it will be arranged correctly in the line below
#'     dplyr::arrange(Injection) # Sort by injection, needed in drift-correction
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------- Check if there are no values in the columns on the same row as QC ----------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   message("Checking if there are no values in the columns (SubjectID, Replicate, Normalization, and Response) on the same row as QC...")
#'
#'   # Transpose the data
#'   raw_data <- base::as.data.frame(base::t(raw_data)) # The data is now Features x Samples (row x column)
#'
#'   check_qc_rule_rows <- function(df) {
#'     # Find columns where row 4 contains SQC, EQC, or QC
#'     qc_cols <- base::which(base::grepl("SQC|EQC|QC", df[4, ], ignore.case = TRUE))
#'
#'     if (base::length(qc_cols) == 0) {
#'       base::message("No QC columns found.")
#'       return(base::invisible(NULL))
#'     }
#'
#'     # Check if rows 2, 3, 7, 8 are empty (NA or empty string) for QC columns
#'     rows_to_check <- base::c(2, 3, 7, 8)
#'     violations <- base::c()
#'
#'     for (col in qc_cols) {
#'       for (row in rows_to_check) {
#'         if (base::nrow(df) >= row && base::ncol(df) >= col) {
#'           if (!base::is.na(df[row, col]) && df[row, col] != "") {
#'             violations <- base::c(violations, col)
#'             break
#'           }
#'         }
#'       }
#'     }
#'
#'     if (base::length(violations) > 0) {
#'       base::stop("QC rule violations found. QC columns in 'Group' row must have empty values in rows 'SubjectID', 'Replicate', 'Normalization', and 'Response'.")
#'       # No specific column because the data frame has been sorted according to 'Injection' row, hence, actual column number is not the same as the original column number.
#'     }
#'
#'     base::message("No QC rule violations found. The rows 'SubjectID', 'Replicate', 'Normalization', and 'Response' do not have values in their QC columns in the 'Group' row, as required.")
#'     return(base::invisible(NULL))
#'   }
#'
#'   # Apply function
#'   check_qc_rule_rows(raw_data)
#'
#'   # Transpose back the data
#'   raw_data <- base::as.data.frame(base::t(raw_data)) # The data is now Samples x Features (row x column)
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------- Check if there are not numbers in features/metabolites ---------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   message("Checking data in the features/metabolites...")
#'
#'   # Metadata
#'   metadata <- raw_data %>%
#'     dplyr::select(dplyr::all_of(required_headers))
#'
#'   # Remove columns that are not part of the features/metabolites
#'   raw_data <- raw_data %>%
#'     dplyr::select(-dplyr::all_of(required_headers))
#'
#'   # Check for numeric-ness of data (can handle conversion)
#'   check_numeric <- function(data, data_name) {
#'     numeric_data   <- data %>% base::as.data.frame()
#'     numeric_data[] <- base::lapply(data, base::as.numeric) # Convert to numeric
#'
#'     # Check if conversion introduced NAs (meaning non-convertible values)
#'     if (base::any(base::is.na(numeric_data) & !base::is.na(data))) {
#'       non_convertible <- base::unique(data[base::is.na(numeric_data) & !base::is.na(data)])
#'       stop(base::paste0("There are non-numeric values in your '", data_name, "' data that cannot be converted to numbers. Here are the non-convertible values: ",
#'                         base::paste(non_convertible, collapse = ", ")))
#'     }
#'
#'     return(numeric_data)  # Return the converted data
#'   }
#'
#'   message("All data in the features/metabolites are either numbers, 0s, or NAs, as required.")
#'
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # ------------------- Combine the metadata and features/metabolites after all preprocessing -------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'   # -------------------------------------------------------------------------------------------------------------
#'
#'   raw_data <- base::cbind(metadata, raw_data) %>%
#'     base::t() %>% # The data is now Features x Samples (row x column)
#'     base::cbind(base::rownames(.), .) %>%
#'     base::`rownames<-`(NULL) %>%
#'     base::as.data.frame()
#'
#'   qualityCheckResults$raw_data <- raw_data
#'
#'   message("Data quality check is complete. The data is ready for analysis.")
#'
#'   return(qualityCheckResults)
#' }
