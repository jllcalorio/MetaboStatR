#' Perform Comprehensive Data Quality Check for Metabolomics Data
#'
#' @description
#' This function imports and validates metabolomics data, performing comprehensive
#' quality controls with strict structural requirements.
#'
#' Key validation checks include:
#' \itemize{
#'   \item Strict Metadata Order: Row 1=Sample, Row 2=Group, Row 3=Batch, Row 4=Injection.
#'   \item Variable Metadata: Allows any number of rows between Injection and Response.
#'   \item Metadata Uniqueness: Ensures no duplicate metadata row names.
#'   \item Numeric validation of all feature/metabolite data.
#'   \item Merging of duplicate samples by mean or median.
#' }
#'
#' @param file_location Character string specifying the file path.
#' @param sheet_name Character string specifying the Excel worksheet name.
#' @param skip_rows Integer specifying the number of rows to skip. Default is 0.
#' @param separator Character string specifying the field separator. Default is ",".
#' @param validate_qc Logical indicating whether to enforce QC validation rules. Default is TRUE.
#' @param allow_missing_optional Logical (currently unused but kept for compatibility). Default is TRUE.
#' @param clean_names Logical indicating whether to clean special characters. Default is TRUE.
#' @param merge_duplicate_samples_by Character string: NULL, "mean", or "median". Default is NULL.
#' @param verbose Logical indicating whether to display progress messages. Default is TRUE.
#'
#' @return A list containing:
#'   \item{raw_data}{Original data as loaded}
#'   \item{metadata_summary}{Summary statistics (Groups as dataframe)}
#'   \item{validation_report}{Detailed validation results}
#'   \item{file_info}{Information about the source file}
#'   \item{processing_log}{Log of all processing steps including execution time}
#'
#' @export
#' @importFrom readxl read_excel excel_sheets
#' @importFrom dplyr mutate arrange select all_of group_by summarise across first left_join ungroup n
#' @importFrom stringr str_replace_all str_trim
#' @importFrom utils read.csv read.delim
#' @importFrom tools file_ext
#' @importFrom stats setNames median
perform_DataQualityCheck <- function(
    file_location = NULL,
    sheet_name = NULL,
    skip_rows = 0,
    separator = ",",
    validate_qc = TRUE,
    allow_missing_optional = TRUE,
    clean_names = TRUE,
    merge_duplicate_samples_by = NULL,
    verbose = TRUE
) {

  # Record start time
  start_time <- Sys.time()

  # Validate parameters upfront
  validate_parameters(merge_duplicate_samples_by, file_location)

  # Initialize results object (Updated structure)
  results <- initialize_results_object(
    file_location, sheet_name, skip_rows, separator,
    validate_qc, allow_missing_optional, clean_names,
    merge_duplicate_samples_by, start_time
  )

  tryCatch({
    # Step 1: Load and validate file
    if (verbose) message("Step 1/9: Reading and validating file...")
    original_data <- load_and_validate_file(file_location, sheet_name, skip_rows,
                                            separator, verbose)
    results$raw_data <- original_data
    file_location <- attr(original_data, "file_location") # Get path if interactively chosen
    results$file_info <- extract_file_info(file_location, original_data)

    # Create working copy
    working_data <- original_data

    # Step 2: Clean names if requested
    if (clean_names) {
      if (verbose) message("Step 2/9: Cleaning special characters in identifiers...")
      working_data <- clean_data_identifiers(working_data, verbose)
    } else {
      if (verbose) message("Step 2/9: Skipping name cleaning (clean_names = FALSE)")
    }

    # Step 3: Validate strict structure and find metadata rows
    if (verbose) message("Step 3/9: Validating strict data structure and locating metadata rows...")
    structure_info <- validate_and_locate_metadata(working_data)
    results$validation_report$structure <- structure_info

    # Step 4: Transpose data for processing (samples as rows)
    if (verbose) message("Step 4/9: Transposing data for processing...")
    processed_data <- transpose_for_processing(working_data)

    # Step 5: Identify Metadata vs Feature columns (CRITICAL STEP)
    if (verbose) message("Step 5/9: Identifying metadata and feature columns...")
    all_processed_colnames <- colnames(processed_data)

    # Define metadata indices based on structure info
    # Includes 1:Sample, 2:Group, 3:Batch, 4:Injection, ... :Response
    metadata_indices <- 1:(structure_info$response_row)
    feature_indices <- (structure_info$feature_start_row):nrow(working_data)

    # Sanity check
    if (any(metadata_indices %in% feature_indices)) {
      stop("Internal error: Metadata and feature rows overlap.", call. = FALSE)
    }

    # These are the *actual* (mangled) column names in processed_data
    metadata_cols_processed <- all_processed_colnames[metadata_indices]
    feature_cols_processed <- all_processed_colnames[feature_indices]

    results$validation_report$column_identification <- list(
      total_metadata_cols = length(metadata_cols_processed),
      total_feature_cols = length(feature_cols_processed),
      metadata_colnames = metadata_cols_processed,
      feature_colnames = feature_cols_processed
    )

    # Step 6: Validate uniqueness constraints
    if (verbose) message("Step 6/9: Checking for duplicate identifiers...")
    duplicate_validation <- validate_uniqueness_constraints(
      processed_data, merge_duplicate_samples_by, feature_cols_processed, verbose
    )
    results$validation_report$duplicates <- duplicate_validation

    # Step 7: Merge duplicate samples if requested
    if (!is.null(merge_duplicate_samples_by)) {
      sample_names <- as.character(processed_data$Sample)
      sample_duplicates <- find_duplicates(sample_names)

      if (length(sample_duplicates) > 0) {
        if (verbose) message("Step 7/9: Merging duplicate samples by ", merge_duplicate_samples_by)

        processed_data <- merge_duplicate_samples(
          processed_data,
          merge_duplicate_samples_by,
          metadata_cols_processed,
          feature_cols_processed,
          verbose
        )

        results$validation_report$merge_info <- list(
          method = merge_duplicate_samples_by,
          duplicates_merged = sample_duplicates,
          samples_after_merge = nrow(processed_data)
        )
      } else {
        if (verbose) message("Step 7/9: No duplicate samples found, skipping merge...")
        results$validation_report$merge_info <- list(
          method = merge_duplicate_samples_by,
          duplicates_merged = character(0),
          message = "No duplicates found"
        )
      }
    } else {
      if (verbose) message("Step 7/9: Skipping duplicate merge (merge_duplicate_samples_by = NULL)")
    }

    # Step 8: Validate QC rules if requested
    if (validate_qc) {
      if (verbose) message("Step 8/9: Validating QC sample rules...")
      qc_validation <- validate_qc_rules(processed_data)
      results$validation_report$qc_rules <- qc_validation
    } else {
      if (verbose) message("Step 8/9: Skipping QC validation (validate_qc = FALSE)")
    }

    # Step 9: Validate numeric data, sort, and finalize
    if (verbose) message("Step 9/9: Validating numeric data and finalizing")

    # Sort by injection first
    processed_data <- sort_by_injection_sequence(processed_data)

    # Use numeric validator
    numeric_validation <- validate_numeric_features(
      processed_data, feature_cols_processed
    )
    results$validation_report$numeric_data <- numeric_validation

    # Finalize data (applies numeric conversion)
    # We create a temporary object for summary generation, but we DO NOT return it
    final_data_for_summary <- finalize_data_structure(
      processed_data, feature_cols_processed
    )

    # Generate summary
    results$metadata_summary <- generate_metadata_summary(
      final_data_for_summary, metadata_cols_processed, feature_cols_processed
    )

    results$processing_log$completion_time <- Sys.time()
    results$processing_log$success <- TRUE

    if (verbose) message("... Data quality check completed successfully!")

  }, error = function(e) {
    results$processing_log$error <- e$message
    results$processing_log$success <- FALSE
    results$processing_log$completion_time <- Sys.time()

    # Calculate final time even on error
    end_time <- results$processing_log$completion_time
    start_time <- results$processing_log$start_time
    if (!is.null(end_time) && !is.null(start_time)) {
      results$processing_log$TotalProcessingTime_in_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
    }

    stop("Data quality check failed: ", e$message, call. = FALSE)
  })

  # Calculate final processing time and put it INSIDE processing_log
  end_time <- results$processing_log$completion_time
  start_time <- results$processing_log$start_time
  if (!is.null(end_time) && !is.null(start_time)) {
    results$processing_log$TotalProcessingTime_in_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }

  class(results) <- c("perform_DataQualityCheck", "list")

  return(results)
}

# ================================================================================
# HELPER FUNCTIONS
# ================================================================================

#' Validate Function Parameters
#' @keywords internal
validate_parameters <- function(merge_method, file_location) {
  if (!is.null(merge_method) && !merge_method %in% c("mean", "median")) {
    stop("merge_duplicate_samples_by must be NULL, 'mean', or 'median'", call. = FALSE)
  }
  if (!is.null(file_location) && !is.character(file_location)) {
    stop("file_location must be a character string or NULL", call. = FALSE)
  }
}

#' Initialize Results Object
#' @keywords internal
initialize_results_object <- function(file_location, sheet_name, skip_rows, separator,
                                      validate_qc, allow_missing_optional,
                                      clean_names, merge_duplicate_samples_by,
                                      start_time) {
  # Removed: quality_checked_data, file_location (top level), version
  # Added: TotalProcessingTime_in_seconds inside processing_log
  list(
    FunctionOrigin = "perform_DataQualityCheck",
    parameters = list(
      file_location = file_location,
      sheet_name = sheet_name,
      skip_rows = skip_rows,
      separator = separator,
      validate_qc = validate_qc,
      allow_missing_optional = allow_missing_optional,
      clean_names = clean_names,
      merge_duplicate_samples_by = merge_duplicate_samples_by
    ),
    processing_log = list(
      start_time = start_time,
      completion_time = NULL,
      success = FALSE,
      error = NULL,
      TotalProcessingTime_in_seconds = NULL
    ),
    file_info = list(),
    validation_report = list(),
    metadata_summary = list(),
    raw_data = NULL
  )
}

#' Load and Validate File
#' @keywords internal
load_and_validate_file <- function(file_location, sheet_name, skip_rows,
                                   separator, verbose) {
  if (is.null(file_location)) {
    file_location <- select_file_interactively()
  }
  if (!file.exists(file_location)) {
    stop("File does not exist: ", file_location, call. = FALSE)
  }

  file_size_mb <- file.size(file_location) / (1024^2)
  if (file_size_mb > 500) {
    warning("Large file detected (", round(file_size_mb, 1), "MB). Processing may be slow.")
  }

  file_ext <- tolower(tools::file_ext(file_location))

  raw_data <- suppressMessages(
    switch(file_ext,
           "xlsx" = load_excel_file(file_location, sheet_name, skip_rows),
           "csv"  = load_csv_file(file_location, separator, skip_rows, verbose),
           "txt"  = load_delimited_file(file_location, separator, skip_rows, verbose),
           "tsv"  = load_delimited_file(file_location, "\t", skip_rows, verbose),
           stop("Unsupported file type: '", file_ext, "'. Supported: xlsx, csv, tsv, txt",
                call. = FALSE)
    )
  )

  if (is.null(raw_data) || nrow(raw_data) == 0 || ncol(raw_data) == 0) {
    stop("Failed to load data or file is empty", call. = FALSE)
  }

  if (verbose) {
    message("... File loaded successfully: ", nrow(raw_data), " rows X ", ncol(raw_data), " columns")
  }

  # Store the chosen file_location in the data object
  attr(raw_data, "file_location") <- file_location
  return(raw_data)
}

#' Select File Interactively
#' @keywords internal
select_file_interactively <- function() {
  message("file_location is NULL. Opening interactive file selection dialog...")
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    file_path <- rstudioapi::selectFile(
      caption = "Select metabolomics data file",
      filter = "Data Files (*.xlsx *.csv *.tsv *.txt)"
    )
  } else {
    file_path <- file.choose()
  }
  if (is.null(file_path) || file_path == "") {
    stop("No file selected. Aborting.", call. = FALSE)
  }
  message("... File selected: ", file_path)
  return(file_path)
}

#' Load Excel File
#' @keywords internal
load_excel_file <- function(file_location, sheet_name, skip_rows) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required for Excel files. Install with: install.packages('readxl')", call. = FALSE)
  }
  tryCatch({
    if (is.null(sheet_name)) {
      available_sheets <- readxl::excel_sheets(file_location)
      sheet_name <- available_sheets[1]
      message("... Using first sheet: ", sheet_name)
    }
    readxl::read_excel(
      path = file_location,
      sheet = sheet_name,
      col_names = FALSE,
      skip = skip_rows,
      .name_repair = "minimal"
    )
  }, error = function(e) {
    stop("Failed to read Excel file '", file_location, "' (sheet: ", sheet_name, "): ", e$message, call. = FALSE)
  })
}

#' Load CSV File
#' @keywords internal
load_csv_file <- function(file_location, separator, skip_rows, verbose) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (verbose) message("... Using data.table::fread for CSV loading")
    tryCatch({
      return(data.table::fread(
        file = file_location,
        header = FALSE,
        sep = separator,
        skip = skip_rows,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        data.table = FALSE,
        showProgress = FALSE
      ))
    }, error = function(e) {
      warning("data.table::fread failed: ", e$message, ". Falling back to utils::read.csv")
    })
  }

  if (verbose) message("... Using utils::read.csv (install 'data.table' for faster loading)")
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
load_delimited_file <- function(file_location, separator, skip_rows, verbose) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (verbose) message("... Using data.table::fread for delimited file loading")
    tryCatch({
      return(data.table::fread(
        file = file_location,
        header = FALSE,
        sep = separator,
        skip = skip_rows,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        data.table = FALSE,
        showProgress = FALSE
      ))
    }, error = function(e) {
      warning("data.table::fread failed: ", e$message, ". Falling back to utils::read.delim")
    })
  }

  if (verbose) message("... Using utils::read.delim (install 'data.table' for faster loading)")
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
    file_path = file_location,
    file_name = basename(file_location),
    file_size_mb = round(file.size(file_location) / (1024^2), 2),
    file_extension = tolower(tools::file_ext(file_location)),
    dimensions_raw = list(rows = nrow(raw_data), columns = ncol(raw_data)),
    modification_time = file.mtime(file_location)
  )
}

#' Clean Data Identifiers
#' @keywords internal
clean_data_identifiers <- function(raw_data, verbose) {
  clean_text_vectorized <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "[^a-zA-Z0-9@._-]", "_")
    x <- stringr::str_replace_all(x, "_{2,}", "_")
    x <- stringr::str_replace_all(x, "^_+|_+$", "")
    x <- stringr::str_trim(x)
    return(x)
  }

  # Clean first column (Metadata/Feature names)
  raw_data[, 1] <- clean_text_vectorized(raw_data[, 1])
  # Clean first row (Sample names)
  raw_data[1, ] <- clean_text_vectorized(raw_data[1, ])

  return(raw_data)
}

#' Validate and Locate Metadata Rows
#' @keywords internal
validate_and_locate_metadata <- function(raw_data) {

  if (nrow(raw_data) < 5) {
    stop("Data must have at least 5 rows (Sample, Group, Batch, Injection, Response)", call. = FALSE)
  }
  if (ncol(raw_data) < 2) {
    stop("Data must have at least 2 columns (Identifier column and at least one sample)", call. = FALSE)
  }

  first_col_values <- as.character(raw_data[, 1])

  # --- STRICT ORDER CHECK (ROWS 1-4) ---
  # Case-insensitive check to be slightly user-friendly, but enforces the word
  if (!grepl("^sample$", first_col_values[1], ignore.case = TRUE)) {
    stop(paste0("Row 1 must be 'Sample'. Found: '", first_col_values[1], "'"), call. = FALSE)
  }
  if (!grepl("^group$", first_col_values[2], ignore.case = TRUE)) {
    stop(paste0("Row 2 must be 'Group'. Found: '", first_col_values[2], "'"), call. = FALSE)
  }
  if (!grepl("^batch$", first_col_values[3], ignore.case = TRUE)) {
    stop(paste0("Row 3 must be 'Batch'. Found: '", first_col_values[3], "'"), call. = FALSE)
  }
  if (!grepl("^injection$", first_col_values[4], ignore.case = TRUE)) {
    stop(paste0("Row 4 must be 'Injection'. Found: '", first_col_values[4], "'"), call. = FALSE)
  }

  # --- FIND RESPONSE ROW ---
  # Must be after Injection (row 4)
  response_indices <- which(grepl("^response$", first_col_values, ignore.case = TRUE))

  if (length(response_indices) == 0) {
    stop("Could not find a 'Response' row in the first column.", call. = FALSE)
  }

  # We take the *first* occurrence of Response that is after row 4
  valid_response_row <- response_indices[response_indices > 4][1]

  if (is.na(valid_response_row)) {
    stop("'Response' row must appear AFTER the 'Injection' row (Row 4).", call. = FALSE)
  }

  # --- UNIQUE METADATA CHECK ---
  # Check all rows from 1 to Response for duplicates
  metadata_names <- first_col_values[1:valid_response_row]

  # Check if names are empty
  if (any(metadata_names == "" | is.na(metadata_names))) {
    stop("One or more metadata rows (between Sample and Response) have empty names.", call. = FALSE)
  }

  # Check for duplicates
  duplicates <- metadata_names[duplicated(metadata_names)]
  if (length(duplicates) > 0) {
    stop("Duplicate metadata names found between Row 1 and Response row: ",
         paste(unique(duplicates), collapse = ", "),
         ". Metadata names must be unique.", call. = FALSE)
  }

  # Calculate Optional Rows
  # Optional rows are anything between 4 (Injection) and the Response row
  optional_rows <- setdiff(5:(valid_response_row - 1), integer(0))
  optional_names <- if (length(optional_rows) > 0) first_col_values[optional_rows] else character(0)

  # Feature start
  feature_start_row <- valid_response_row + 1
  total_features <- nrow(raw_data) - valid_response_row

  if (total_features <= 0) {
    stop("No feature rows found after the 'Response' row.", call. = FALSE)
  }

  return(list(
    sample_row = 1,
    group_row = 2,
    batch_row = 3,
    injection_row = 4,
    response_row = valid_response_row,
    optional_rows = optional_rows,
    optional_names = optional_names,
    feature_start_row = feature_start_row,
    total_features = total_features
  ))
}

#' Transpose Data for Processing
#' @keywords internal
transpose_for_processing <- function(raw_data) {
  # Transpose
  transposed_data <- as.data.frame(t(raw_data), stringsAsFactors = FALSE)

  # Set column names from first row (which was raw_data[,1])
  new_colnames <- as.character(transposed_data[1, ])

  # Ensure column names are syntactically valid and unique
  new_colnames <- make.names(new_colnames, unique = TRUE)
  colnames(transposed_data) <- new_colnames

  # Remove the first row (which is now headers)
  transposed_data <- transposed_data[-1, , drop = FALSE]

  # Reset row names (which were sample names, not needed)
  rownames(transposed_data) <- NULL

  return(transposed_data)
}

#' Validate Uniqueness Constraints
#' @keywords internal
validate_uniqueness_constraints <- function(processed_data, merge_duplicate_samples_by,
                                            feature_cols_processed, verbose) {
  validation_results <- list()

  # 1. Check sample names
  if (!"Sample" %in% colnames(processed_data)) {
    stop("Internal error: 'Sample' column not found after transpose.", call. = FALSE)
  }
  sample_names <- as.character(processed_data$Sample)
  sample_duplicates <- find_duplicates(sample_names)

  if (length(sample_duplicates) > 0) {
    if (is.null(merge_duplicate_samples_by)) {
      stop("Duplicate sample names found: ", paste(sample_duplicates, collapse = ", "),
           ". Set merge_duplicate_samples_by to 'mean' or 'median' to merge.",
           call. = FALSE)
    } else {
      validation_results$sample_names <- paste0("OK (Duplicates found and will be merged: ",
                                                paste(sample_duplicates, collapse = ", "), ")")
    }
  } else {
    validation_results$sample_names <- "All sample names are unique"
  }

  # 2. Check injection sequence
  if ("Injection" %in% colnames(processed_data)) {
    injection_seq <- as.character(processed_data$Injection)

    non_blank_injections <- injection_seq[injection_seq != "" & !is.na(injection_seq)]
    suppressWarnings(numeric_injections <- as.numeric(non_blank_injections))

    if (any(is.na(numeric_injections))) {
      stop("Non-numeric values found in Injection sequence", call. = FALSE)
    }

    injection_duplicates <- find_duplicates(numeric_injections)
    if (length(injection_duplicates) > 0) {
      validation_results$injection_sequence <- paste0("Duplicate injection sequences found: ",
                                                      paste(injection_duplicates, collapse = ", "))
    } else {
      validation_results$injection_sequence <- "All injection sequences are unique"
    }
  }

  feature_duplicates <- find_duplicates(feature_cols_processed)
  if (length(feature_duplicates) > 0) {
    stop(paste("Internal error: Duplicate feature names detected after processing: ",
               paste(feature_duplicates, collapse = ", ")), call. = FALSE)
  }
  validation_results$feature_names <- "All feature/metabolite names are unique"

  return(validation_results)
}

#' Find Duplicates in Vector
#' @keywords internal
find_duplicates <- function(x) {
  unique(x[duplicated(x)])
}

#' Merge Duplicate Samples
#' @keywords internal
merge_duplicate_samples <- function(processed_data, method,
                                              metadata_cols_processed,
                                              feature_cols_processed,
                                              verbose) {

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for merging duplicates.", call. = FALSE)
  }

  present_metadata_cols <- intersect(metadata_cols_processed, colnames(processed_data))
  present_feature_cols <- intersect(feature_cols_processed, colnames(processed_data))

  if (length(present_feature_cols) == 0) {
    warning("No feature columns found during merge. Merging metadata only.")
  }

  agg_func <- if (method == "mean") mean else median

  suppressWarnings({
    processed_data <- processed_data %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(present_feature_cols), as.numeric))
  })

  metadata_for_first <- setdiff(present_metadata_cols, c("Sample", "Injection"))

  merged_data <- processed_data %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(metadata_for_first), dplyr::first),
      dplyr::across(dplyr::all_of(present_feature_cols), ~ agg_func(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  if ("Injection" %in% present_metadata_cols) {
    injection_summary <- processed_data %>%
      dplyr::group_by(Sample) %>%
      dplyr::summarise(
        Injection = min(as.numeric(as.character(Injection)), na.rm = TRUE),
        .groups = "drop"
      )
    merged_data <- dplyr::left_join(merged_data, injection_summary, by = "Sample")
  }

  if (verbose) {
    message("... Merged duplicate samples. New sample count: ", nrow(merged_data))
  }

  return(as.data.frame(merged_data))
}

#' Sort Data by Injection Sequence
#' @keywords internal
sort_by_injection_sequence <- function(processed_data) {
  if (!"Injection" %in% colnames(processed_data)) {
    return(processed_data)
  }
  tryCatch({
    processed_data$Injection <- as.numeric(as.character(processed_data$Injection))
    processed_data <- processed_data[order(processed_data$Injection), ]
  }, error = function(e) {
    warning("Could not sort by injection sequence: ", e$message)
  })
  rownames(processed_data) <- NULL
  return(processed_data)
}

#' Validate QC Rules
#' @keywords internal
validate_qc_rules <- function(processed_data) {
  if (!"Group" %in% colnames(processed_data)) {
    return(list(message = "No 'Group' column found, skipping QC validation", qc_samples = 0))
  }

  qc_indices <- which(grepl("QC", processed_data$Group, ignore.case = TRUE))
  if (length(qc_indices) == 0) {
    return(list(message = "No QC samples found", qc_samples = 0))
  }

  rows_to_check <- c("SubjectID", "Replicate", "Normalization", "Response")
  violations <- c()

  for (row_name in rows_to_check) {
    if (row_name %in% colnames(processed_data)) {
      qc_values <- processed_data[qc_indices, row_name]
      problematic <- !is.na(qc_values) & qc_values != "" & qc_values != "0"

      if (any(problematic)) {
        bad_samples <- processed_data$Sample[qc_indices[which(problematic)]]
        violations <- c(violations, paste0(row_name, " (samples: ", paste(unique(bad_samples), collapse = ", "), ")"))
      }
    }
  }

  if (length(violations) > 0) {
    stop("QC rule violations: QC samples must have empty/0 values in ",
         paste(rows_to_check, collapse = ", "),
         " columns. Violations found in: ", paste(violations, collapse = "; "),
         call. = FALSE)
  }

  return(list(
    message = "QC rules validated successfully",
    qc_samples_found = length(qc_indices)
  ))
}

#' Validate Numeric Features
#' @keywords internal
validate_numeric_features <- function(processed_data, feature_cols_processed) {
  if (length(feature_cols_processed) == 0) {
    stop("No feature/metabolite columns were identified for numeric validation.", call. = FALSE)
  }
  feature_data <- processed_data[, feature_cols_processed, drop = FALSE]
  if (nrow(feature_data) == 0) {
    stop("No sample rows found to validate.", call. = FALSE)
  }
  feature_matrix <- as.matrix(feature_data)
  original_na <- is.na(feature_matrix) | feature_matrix == "" | feature_matrix == "0"

  suppressWarnings({
    numeric_matrix <- matrix(
      as.numeric(feature_matrix),
      nrow = nrow(feature_matrix),
      dimnames = dimnames(feature_matrix)
    )
  })
  conversion_na <- is.na(numeric_matrix)
  problematic_cells <- conversion_na & !original_na

  if (any(problematic_cells)) {
    problem_indices <- which(problematic_cells, arr.ind = TRUE)
    problem_cols <- unique(problem_indices[, "col"])
    conversion_issues <- list()
    for (col_idx in problem_cols) {
      col_name <- colnames(feature_matrix)[col_idx]
      problem_rows <- problem_indices[problem_indices[, "col"] == col_idx, "row"]
      problem_values <- unique(feature_matrix[problem_rows, col_idx])
      conversion_issues[[col_name]] <- paste(problem_values, collapse = ", ")
    }
    stop_message <- paste(names(conversion_issues), conversion_issues, sep = ": ", collapse = "\n")
    stop("Non-numeric values found in feature columns:\n", stop_message, call. = FALSE)
  }

  return(list(
    message = "All feature data is numeric or convertible to numeric (NA/empty/0)",
    feature_count = ncol(feature_data),
    sample_count = nrow(feature_data)
  ))
}

#' Finalize Data Structure
#' @keywords internal
finalize_data_structure <- function(processed_data, feature_cols_processed) {
  if (length(feature_cols_processed) > 0) {
    processed_data[feature_cols_processed] <- lapply(
      processed_data[feature_cols_processed],
      function(x) as.numeric(as.character(x))
    )
    processed_data[feature_cols_processed] <- lapply(
      processed_data[feature_cols_processed],
      function(x) {
        x[is.na(x)] <- 0
        return(x)
      }
    )
  }
  if ("Injection" %in% colnames(processed_data)) {
    processed_data$Injection <- as.numeric(as.character(processed_data$Injection))
  }
  return(processed_data)
}

#' Generate Metadata Summary
#' @keywords internal
generate_metadata_summary <- function(processed_data, metadata_cols_processed,
                                      feature_cols_processed) {

  summary_stats <- list(
    total_samples = nrow(processed_data),
    total_features = length(feature_cols_processed),
    metadata_columns = metadata_cols_processed
  )

  # Group summary - CHANGED TO DATA FRAME
  if ("Group" %in% colnames(processed_data)) {
    group_tbl <- table(processed_data$Group)
    summary_stats$groups <- data.frame(
      Group = names(group_tbl),
      Count = as.numeric(group_tbl),
      stringsAsFactors = FALSE
    )
    summary_stats$qc_samples_found <- sum(grepl("QC", processed_data$Group, ignore.case = TRUE))
  }

  # Batch summary
  if ("Batch" %in% colnames(processed_data)) {
    summary_stats$batches <- length(unique(processed_data$Batch))
    summary_stats$unique_batches <- unique(processed_data$Batch)
  }

  # Injection summary
  if ("Injection" %in% colnames(processed_data)) {
    summary_stats$injection_range <- range(processed_data$Injection, na.rm = TRUE)
  }

  return(summary_stats)
}

# S3 Methods
#' @export
print.perform_DataQualityCheck <- function(x, ...) {
  cat("=== Metabolomics Data Quality Check ===\n")
  cat("File:      ", x$file_info$file_name, "\n")
  cat("Raw Dim:   ", x$file_info$dimensions_raw$rows, "rows x", x$file_info$dimensions_raw$columns, "cols\n")
  cat("Status:    ", ifelse(x$processing_log$success, "Pass", "Fail"), "\n")
  cat("Processed: ", format(x$processing_log$completion_time, "%Y-%m-%d %H:%M:%S"), "\n")

  if (!is.null(x$validation_report$numeric_data$message)) {
    cat("Validation:", x$validation_report$numeric_data$message, "\n")
  }
  invisible(x)
}

#' @export
summary.perform_DataQualityCheck <- function(object, ...) {
  ans <- list(
    file = object$file_info$file_name,
    samples = object$metadata_summary$total_samples,
    features = object$metadata_summary$total_features,
    groups = object$metadata_summary$groups,
    qc_count = object$metadata_summary$qc_samples_found,
    merge_info = object$validation_report$merge_info,
    duration = object$processing_log$TotalProcessingTime_in_seconds
  )
  class(ans) <- "summary.perform_DataQualityCheck"
  return(ans)
}

#' @export
print.summary.perform_DataQualityCheck <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("Data Quality Check Summary\n")
  cat("---------------------------------------\n")
  cat(sprintf("%-20s %s\n", "File Name:", x$file))
  cat(sprintf("%-20s %d\n", "Total Samples:", x$samples))
  cat(sprintf("%-20s %d\n", "Total Features:", x$features))
  cat(sprintf("%-20s %d\n", "QC Samples:", x$qc_count))
  cat(sprintf("%-20s %.2f seconds\n", "Processing Time:", x$duration))

  if (!is.null(x$merge_info$method)) {
    cat("\n-- Replicate Merging --\n")
    cat("Method:", x$merge_info$method, "\n")
    cat("Samples after merge:", x$merge_info$samples_after_merge, "\n")
  }

  if (!is.null(x$groups)) {
    cat("\n-- Group Distribution --\n")
    print(x$groups, row.names = FALSE)
  }
  invisible(x)
}
