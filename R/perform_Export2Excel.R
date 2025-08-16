#' Export Data Frames to Excel Workbook
#'
#' @description
#' Exports data frames from named lists or multiple lists to multi-sheet Excel workbooks.
#' The function automatically detects data frames and tibbles within the input, creates
#' appropriately named worksheets, and handles various edge cases including sheet name
#' length restrictions, duplicate names, and invalid characters. Each data frame becomes
#' a separate worksheet in the output Excel file. The function can accept either a single
#' named list or multiple lists, and will create separate Excel files for each list when
#' multiple file names are provided.
#'
#' @param results Either a named list containing data frames/tibbles, or multiple lists
#'   that can be passed as: results = list1, or results = c(list1, list2),
#'   or results = list(list1, list2, list3). Only data frames and tibbles will be exported,
#'   but empty data frames will still be included as worksheets.
#' @param folder_name Character string specifying the folder name where the Excel file(s)
#'   will be saved. The folder will be created if it doesn't exist. Default is "Results_Folder".
#' @param file_name Character vector specifying the base name(s) of the Excel file(s) (without
#'   extension). If a single name is provided, it will be used for all lists. If multiple names
#'   are provided, each list will get its own file with the corresponding name. A timestamp
#'   will be automatically appended. Default is "Results_Export".
#' @param include_timestamp Logical indicating whether to include a timestamp in the
#'   filename(s). Default is TRUE.
#' @param overwrite Logical indicating whether to overwrite existing files. If FALSE and file
#'   exists, an error will be thrown. If TRUE but file is currently open in Excel, the function
#'   will detect this and provide a specific error message. Default is TRUE.
#' @param row_names Logical indicating whether to include row names in the Excel output.
#'   Default is TRUE.
#' @param freeze_first_row Logical indicating whether to freeze the first row (headers)
#'   in each worksheet. Default is TRUE.
#' @param auto_width Logical indicating whether to auto-adjust column widths. Default is TRUE.
#' @param max_sheets Integer specifying maximum number of worksheets to create per file.
#'   Default is 255 (Excel limit).
#'
#' @return Invisibly returns a character vector of full paths to the created Excel file(s).
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Validates and processes input parameters (single list or multiple lists)
#'   \item Processes multiple lists and creates separate files when multiple file names provided
#'   \item Filters the list(s) to include only data frames and tibbles (including empty ones)
#'   \item Sanitizes sheet names to comply with Excel naming conventions
#'   \item Creates multi-sheet Excel workbook(s) with formatted headers
#'   \item Applies optional formatting (frozen headers, auto-width columns)
#'   \item Saves the workbook(s) with optional timestamps
#'   \item Detects if Excel files are currently open and provides specific error messages
#' }
#'
#' Multiple list and file handling:
#' \itemize{
#'   \item If single file_name provided: all lists merged into one Excel file
#'   \item If multiple file_names provided: each list becomes a separate Excel file
#'   \item If more lists than file names: remaining lists use the last file name with suffixes
#'   \item Empty data frames are included as worksheets (preserving structure information)
#' }
#'
#' Sheet names are automatically sanitized to:
#' \itemize{
#'   \item Remove or replace invalid characters
#'   \item Ensure uniqueness by appending numbers to duplicates
#'   \item Truncate names exceeding 31 characters (Excel limit)
#'   \item Preserve readability by intelligently shortening names
#' }
#'
#' @section Dependencies:
#' Requires the \code{openxlsx} package for Excel file operations.
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' list1 <- list(
#'   summary_stats = data.frame(mean = c(1, 2, 3), sd = c(0.1, 0.2, 0.3)),
#'   raw_data = data.frame(x = 1:10, y = rnorm(10)),
#'   empty_data = data.frame()  # This will still be included
#' )
#'
#' list2 <- list(
#'   processed_data = data.frame(a = 1:5, b = letters[1:5]),
#'   metadata = "This will be ignored"
#' )
#'
#' list3 <- list(
#'   final_results = data.frame(value = runif(8), category = LETTERS[1:8])
#' )
#'
#' # Method 1: Single list, single file
#' perform_Export2Excel(results = list1, file_name = "Analysis_Results")
#'
#' # Method 2: Multiple lists, single file (all merged)
#' perform_Export2Excel(
#'   results = c(list1, list2, list3),
#'   file_name = "Combined_Analysis"
#' )
#'
#' # Method 3: Multiple lists, multiple files (separate files)
#' perform_Export2Excel(
#'   results = c(list1, list2, list3),
#'   file_name = c("Analysis", "Processing", "Final_Results")
#' )
#'
#' # Method 4: Multiple lists, partial file names
#' perform_Export2Excel(
#'   results = c(list1, list2, list3),
#'   file_name = c("Primary_Analysis", "Secondary_Analysis")
#'   # list3 will use "Secondary_Analysis_2"
#' )
#'
#' # Advanced usage with custom options
#' perform_Export2Excel(
#'   results = c(list1, list2),
#'   folder_name = "Analysis_Results",
#'   file_name = c("Dataset_A", "Dataset_B"),
#'   include_timestamp = FALSE,
#'   row_names = TRUE,
#'   freeze_first_row = TRUE,
#'   auto_width = TRUE,
#'   max_sheets = 50
#' )
#' }
#'
#' @seealso
#' \code{\link[openxlsx]{createWorkbook}}, \code{\link[openxlsx]{addWorksheet}},
#' \code{\link[openxlsx]{writeData}}
#'
#' @export
perform_Export2Excel <- function(
    results,
    folder_name = "Results_Folder",
    file_name = "Results_Export",
    include_timestamp = TRUE,
    overwrite = TRUE,
    row_names = TRUE,
    freeze_first_row = TRUE,
    auto_width = TRUE,
    max_sheets = 255
) {

  # Start timing for performance monitoring
  start_time <- Sys.time()

  # Process and validate inputs
  processed_results <- .process_multiple_inputs(results)
  .validate_inputs_export2excel(processed_results, folder_name, file_name, include_timestamp,
                                overwrite, row_names, freeze_first_row, auto_width, max_sheets)

  # Check for openxlsx package
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required but not installed. Please install it using: install.packages('openxlsx')",
         call. = FALSE)
  }

  # Create output directory
  output_path <- .create_output_directory(folder_name)

  # Process file names and create corresponding files
  file_paths <- .create_multiple_excel_files(processed_results, output_path, file_name,
                                             include_timestamp, overwrite, row_names,
                                             freeze_first_row, auto_width, max_sheets)

  # Performance summary
  end_time <- Sys.time()
  processing_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

  message(sprintf("Processing completed in %s seconds", processing_time))
  message(sprintf("Total Excel files created: %d", length(file_paths)))

  return(invisible(file_paths))
}


#' Process Multiple Input Lists (Simplified)
#' @noRd
.process_multiple_inputs <- function(results) {
  if (missing(results) || is.null(results)) {
    stop("No input provided. Please provide at least one list containing data frames.",
         call. = FALSE)
  }

  if (!is.list(results)) {
    stop("'results' must be a list or multiple lists.", call. = FALSE)
  }

  # Check if this is a list of lists (multiple lists scenario)
  if (.is_list_of_lists(results)) {
    return(results)
  } else {
    # Single named list scenario - wrap in list to standardize processing
    if (is.null(names(results))) {
      # Unnamed list - assign generic names to data frame elements
      df_indices <- vapply(results, .is_dataframe_or_tibble, logical(1))
      if (any(df_indices)) {
        names(results)[df_indices] <- paste0("Dataset_", seq_len(sum(df_indices)))
      }
      if (any(!df_indices)) {
        names(results)[!df_indices] <- paste0("Other_", seq_len(sum(!df_indices)))
      }
    }
    return(list(results))  # Wrap single list
  }
}


#' Check if input is a list of lists
#' @noRd
.is_list_of_lists <- function(x) {
  if (!is.list(x) || length(x) == 0) return(FALSE)

  # Check if most elements are lists and at least one contains data frames
  list_elements <- vapply(x, is.list, logical(1))
  if (sum(list_elements) < length(x) * 0.5) return(FALSE)

  # Check if any of the list elements contain data frames
  contains_dfs <- any(vapply(x[list_elements], function(sublist) {
    any(vapply(sublist, .is_dataframe_or_tibble, logical(1)))
  }, logical(1)))

  return(contains_dfs)
}


#' Check if object is data frame or tibble
#' @noRd
.is_dataframe_or_tibble <- function(x) {
  is.data.frame(x) || (requireNamespace("tibble", quietly = TRUE) && tibble::is_tibble(x))
}


#' Create Multiple Excel Files
#' @noRd
.create_multiple_excel_files <- function(processed_results, output_path, file_names,
                                         include_timestamp, overwrite, row_names,
                                         freeze_first_row, auto_width, max_sheets) {

  n_lists <- length(processed_results)
  n_file_names <- length(file_names)

  # Generate file names for each list
  final_file_names <- character(n_lists)

  if (n_file_names == 1) {
    # Single file name - merge all lists into one file
    if (n_lists > 1) {
      # Merge all lists
      merged_list <- .merge_all_lists(processed_results)
      processed_results <- list(merged_list)
      final_file_names <- file_names[1]
    } else {
      final_file_names <- file_names[1]
    }
  } else {
    # Multiple file names
    for (i in seq_len(n_lists)) {
      if (i <= n_file_names) {
        final_file_names[i] <- file_names[i]
      } else {
        # Use last file name with suffix
        base_name <- file_names[n_file_names]
        suffix_num <- i - n_file_names + 1
        final_file_names[i] <- paste0(base_name, "_", suffix_num)
      }
    }
  }

  # Create Excel files
  created_files <- character(0)

  for (i in seq_along(processed_results)) {
    current_list <- processed_results[[i]]
    current_file_name <- final_file_names[i]

    # Filter and validate data frames (including empty ones)
    df_list <- .extract_dataframes_optimized(current_list, max_sheets, include_empty = TRUE)

    if (length(df_list) == 0) {
      message(sprintf("No data frames found in list %d. Skipping file: %s", i, current_file_name))
      next
    }

    # Generate file path
    full_path <- .generate_file_path(output_path, current_file_name, include_timestamp)

    # Check if file is open in Excel
    .check_file_availability(full_path, overwrite)

    # Create and populate workbook
    wb <- .create_workbook_optimized(df_list, row_names, freeze_first_row, auto_width, include_empty = TRUE)

    # Save workbook with error handling
    .save_workbook_safely(wb, full_path, overwrite)

    message(sprintf("Excel file created: %s (Worksheets: %d)", full_path, length(df_list)))
    created_files <- c(created_files, full_path)
  }

  return(created_files)
}


#' Merge All Lists into Single List
#' @noRd
.merge_all_lists <- function(list_collection) {
  merged_list <- list()

  for (i in seq_along(list_collection)) {
    current_list <- list_collection[[i]]
    list_name <- names(list_collection)[i]

    if (is.null(list_name) || list_name == "") {
      list_name <- paste0("List", i)
    }

    # Handle unnamed elements in current list
    current_names <- names(current_list)
    if (is.null(current_names)) {
      df_indices <- vapply(current_list, .is_dataframe_or_tibble, logical(1))
      current_names <- character(length(current_list))
      current_names[df_indices] <- paste0("Dataset_", seq_len(sum(df_indices)))
      current_names[!df_indices] <- paste0("Other_", seq_len(sum(!df_indices)))
      names(current_list) <- current_names
    } else {
      # Fix any unnamed elements
      unnamed_indices <- which(is.na(current_names) | current_names == "")
      if (length(unnamed_indices) > 0) {
        current_names[unnamed_indices] <- paste0("Item_", unnamed_indices)
        names(current_list) <- current_names
      }
    }

    # Add elements to merged list with conflict resolution
    for (name in names(current_list)) {
      new_name <- name

      # Handle name conflicts by prefixing with list name
      if (new_name %in% names(merged_list)) {
        new_name <- paste(list_name, name, sep = "_")
      }

      # Ensure uniqueness
      counter <- 1
      base_name <- new_name
      while (new_name %in% names(merged_list)) {
        new_name <- paste(base_name, counter, sep = "_")
        counter <- counter + 1
      }

      merged_list[[new_name]] <- current_list[[name]]
    }
  }

  return(merged_list)
}


#' Enhanced Input Validation
#' @noRd
.validate_inputs_export2excel <- function(results, folder_name, file_name, include_timestamp,
                                          overwrite, row_names, freeze_first_row, auto_width, max_sheets) {

  # Check results parameter
  if (is.null(results) || !is.list(results) || length(results) == 0) {
    stop("'results' must be a non-empty list.", call. = FALSE)
  }

  # Check folder_name
  if (!is.character(folder_name) || length(folder_name) != 1 || nchar(trimws(folder_name)) == 0) {
    stop("'folder_name' must be a non-empty character string.", call. = FALSE)
  }

  # Check file_name (can be vector)
  if (!is.character(file_name) || length(file_name) == 0 || any(nchar(trimws(file_name)) == 0)) {
    stop("'file_name' must be a non-empty character vector.", call. = FALSE)
  }

  # Check logical parameters
  logical_params <- list(
    include_timestamp = include_timestamp,
    overwrite = overwrite,
    row_names = row_names,
    freeze_first_row = freeze_first_row,
    auto_width = auto_width
  )

  for (param_name in names(logical_params)) {
    param_value <- logical_params[[param_name]]
    if (!is.logical(param_value) || length(param_value) != 1 || is.na(param_value)) {
      stop(sprintf("'%s' must be a single logical value (TRUE or FALSE).", param_name), call. = FALSE)
    }
  }

  # Check max_sheets
  if (!is.numeric(max_sheets) || length(max_sheets) != 1 || max_sheets < 1 || max_sheets > 1024) {
    stop("'max_sheets' must be a single integer between 1 and 1024.", call. = FALSE)
  }

  # Check for invalid characters in folder and file names
  invalid_chars <- "[<>:\"/\\|?*]"
  if (grepl(invalid_chars, folder_name)) {
    stop("'folder_name' contains invalid characters. Avoid: < > : \" / \\ | ? *", call. = FALSE)
  }
  if (any(grepl(invalid_chars, file_name))) {
    stop("'file_name' contains invalid characters. Avoid: < > : \" / \\ | ? *", call. = FALSE)
  }
}


#' Optimized Data Frame Extraction (Including Empty Ones)
#' @noRd
.extract_dataframes_optimized <- function(results, max_sheets, include_empty = TRUE) {
  # Pre-allocate logical vector for better performance
  n_items <- length(results)
  df_indices <- logical(n_items)

  # Vectorized check for data frames
  for (i in seq_len(n_items)) {
    df_indices[i] <- .is_dataframe_or_tibble(results[[i]])
  }

  df_list <- results[df_indices]

  # Apply sheet limit
  if (length(df_list) > max_sheets) {
    warning(sprintf("Found %d data frames but max_sheets is set to %d. Only the first %d will be exported.",
                    length(df_list), max_sheets, max_sheets), call. = FALSE)
    df_list <- df_list[1:max_sheets]
  }

  # Report excluded items
  if (length(df_list) < n_items) {
    excluded_count <- n_items - length(df_list)
    message(sprintf("Excluded %d non-data frame objects from export.", excluded_count))
  }

  return(df_list)
}


#' Create Output Directory
#' @noRd
.create_output_directory <- function(folder_name) {
  tryCatch({
    folder_name <- trimws(folder_name)
    if (!dir.exists(folder_name)) {
      dir.create(folder_name, recursive = TRUE)
      message(sprintf("Created directory: %s", folder_name))
    }
    return(normalizePath(folder_name, mustWork = FALSE))
  }, error = function(e) {
    stop(sprintf("Failed to create directory '%s': %s", folder_name, e$message), call. = FALSE)
  })
}


#' Generate File Path with Optional Timestamp
#' @noRd
.generate_file_path <- function(output_path, file_name, include_timestamp) {
  file_name <- trimws(file_name)

  if (include_timestamp) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    final_name <- paste0(file_name, "_", timestamp, ".xlsx")
  } else {
    final_name <- paste0(file_name, ".xlsx")
  }

  return(file.path(output_path, final_name))
}


#' Check if Excel File is Currently Open
#' @noRd
.check_file_availability <- function(file_path, overwrite) {
  if (!file.exists(file_path)) {
    return(TRUE)  # File doesn't exist, safe to create
  }

  if (!overwrite) {
    stop(sprintf("File already exists and overwrite is FALSE: %s", file_path), call. = FALSE)
  }

  # Test if file is locked (open in Excel)
  tryCatch({
    # Try to open file in append mode to test if it's locked
    con <- file(file_path, open = "r+b")
    close(con)
    return(TRUE)  # File is available
  }, error = function(e) {
    if (grepl("cannot open the connection|Permission denied|sharing violation", e$message, ignore.case = TRUE)) {
      stop(sprintf("Excel file is currently open in another application. Please close the file and try again: %s",
                   file_path), call. = FALSE)
    } else {
      stop(sprintf("Unable to access file: %s. Error: %s", file_path, e$message), call. = FALSE)
    }
  })
}


#' Optimized Sheet Name Sanitization
#' @noRd
.sanitize_sheet_name <- function(name, existing_names = character(0)) {
  # Remove invalid characters for Excel sheet names
  name <- gsub("[\\[\\]:*?/\\\\]", "_", name)

  # Remove leading/trailing whitespace and multiple underscores
  name <- trimws(name)
  name <- gsub("_{2,}", "_", name)
  name <- gsub("^_|_$", "", name)

  # Ensure non-empty name
  if (nchar(name) == 0) {
    name <- "Sheet"
  }

  # Truncate if too long (Excel limit is 31 characters)
  if (nchar(name) > 31) {
    name <- .intelligent_truncate(name, 31)
  }

  # Ensure uniqueness
  original_name <- name
  counter <- 1
  while (name %in% existing_names) {
    suffix <- paste0("_", counter)
    max_base_length <- 31 - nchar(suffix)
    if (nchar(original_name) > max_base_length) {
      base_name <- substr(original_name, 1, max_base_length)
    } else {
      base_name <- original_name
    }
    name <- paste0(base_name, suffix)
    counter <- counter + 1

    # Prevent infinite loops
    if (counter > 1000) {
      name <- paste0("Sheet_", sample.int(9999, 1))
      break
    }
  }

  return(name)
}


#' Intelligent Truncation of Long Names
#' @noRd
.intelligent_truncate <- function(name, max_length) {
  if (nchar(name) <= max_length) {
    return(name)
  }

  # Strategy 1: Keep first and last parts, remove middle
  if (max_length >= 10 && nchar(name) > max_length) {
    first_part_length <- floor(max_length / 2) - 1
    last_part_length <- max_length - first_part_length - 2
    first_part <- substr(name, 1, first_part_length)
    last_part <- substr(name, nchar(name) - last_part_length + 1, nchar(name))
    truncated <- paste0(first_part, "..", last_part)

    if (nchar(truncated) <= max_length) {
      return(truncated)
    }
  }

  # Strategy 2: Remove vowels (except first character)
  if (nchar(name) > 1) {
    first_char <- substr(name, 1, 1)
    rest_chars <- substr(name, 2, nchar(name))
    rest_no_vowels <- gsub("[aeiouAEIOU]", "", rest_chars)
    truncated <- paste0(first_char, rest_no_vowels)

    if (nchar(truncated) <= max_length) {
      return(truncated)
    }
  }

  # Strategy 3: Simple truncation
  return(substr(name, 1, max_length))
}


#' Optimized Workbook Creation (Including Empty Data Frames)
#' @noRd
.create_workbook_optimized <- function(df_list, row_names, freeze_first_row, auto_width, include_empty = TRUE) {
  wb <- openxlsx::createWorkbook()

  # Pre-create styles for better performance
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#D3D3D3",
    border = "TopBottomLeftRight"
  )

  # Track used sheet names to ensure uniqueness
  used_sheet_names <- character(0)
  successful_sheets <- 0

  for (i in seq_along(df_list)) {
    original_name <- names(df_list)[i]
    data <- df_list[[i]]

    # Sanitize sheet name
    sheet_name <- .sanitize_sheet_name(original_name, used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)

    # Notify if sheet name was changed
    if (sheet_name != original_name) {
      message(sprintf("Sheet name adjusted: '%s' -> '%s'", original_name, sheet_name))
    }

    tryCatch({
      # Add worksheet
      openxlsx::addWorksheet(wb, sheetName = sheet_name)

      # Write data (even if empty)
      openxlsx::writeData(wb, sheet = sheet_name, x = data, rowNames = row_names)

      # Apply formatting (even for empty data frames if they have column structure)
      if (nrow(data) > 0 || ncol(data) > 0) {
        if (freeze_first_row && ncol(data) > 0) {
          openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
        }

        if (auto_width && ncol(data) > 0) {
          col_range <- if (row_names) 1:(ncol(data) + 1) else 1:ncol(data)
          openxlsx::setColWidths(wb, sheet = sheet_name, cols = col_range, widths = "auto")
        }

        # Style header row (if columns exist)
        if (ncol(data) > 0) {
          col_range <- if (row_names) 1:(ncol(data) + 1) else 1:ncol(data)
          openxlsx::addStyle(wb, sheet = sheet_name, style = header_style,
                             rows = 1, cols = col_range, gridExpand = TRUE)
        }
      }

      successful_sheets <- successful_sheets + 1

    }, error = function(e) {
      warning(sprintf("Failed to create worksheet '%s': %s", sheet_name, e$message), call. = FALSE)
    })
  }

  if (successful_sheets == 0) {
    stop("Failed to create any worksheets. Please check your data.", call. = FALSE)
  }

  return(wb)
}


#' Safe Workbook Saving
#' @noRd
.save_workbook_safely <- function(wb, full_path, overwrite) {
  tryCatch({
    openxlsx::saveWorkbook(wb, file = full_path, overwrite = overwrite)
  }, error = function(e) {
    # Try to provide more specific error information
    if (grepl("permission|sharing violation|being used by another process", e$message, ignore.case = TRUE)) {
      stop(sprintf("Excel file is currently open in another application. Please close the file and try again: %s",
                   full_path), call. = FALSE)
    } else if (grepl("path", e$message, ignore.case = TRUE)) {
      stop(sprintf("Invalid file path: %s", full_path), call. = FALSE)
    } else {
      stop(sprintf("Failed to save Excel file: %s", e$message), call. = FALSE)
    }
  })
}
