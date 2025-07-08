#' Perform Export Data Frames to Excel
#'
#' @description
#' This function exports the data frames from a list to Excel. The data frames are detected automatically. All data frames
#' in a list are combined in one (1) Excel file.
#'
#' @param results List. A list of results from `performPreprocessingPeakData`, `performDimensionReduction`, or `performComparativeAnalysis`, etc. There must only be one at a time. The list may or may not contain a data frame.
#' @param folder_name String. Creates a folder in the current working directory where the Excel file will be exported.
#' @param file_name String. The file name of the results.
#'
#' @returns An Excel file.
#' @export
#'
#' @examples
#' \dontrun{
#' performExport2Excel(
#'   results     = results_from_performRegression_function,
#'   folder_name = "Regression Results",
#'   file_name   = "Regression Results")
#'}
#'
performExport2Excel <- function(
    results,
    folder_name = "Results_Folder",
    file_name = "Results Exports"
) {

  # NOTE: THIS WORKS BUT FOR RESULTS THAT ARE NOT IN A LIST

  library(openxlsx)

  # Generate a timestamp for the file name
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  # Create folder if it doesn't exist
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }

  # Initialize a workbook
  wb <- createWorkbook()

  # Loop through each element in the list
  for (original_name in names(results)) {
    data <- results[[original_name]]
    sheet_name <- original_name

    # Check if the name is too long
    if (nchar(sheet_name) > 31) {
      # Remove vowels (case-insensitive)
      sheet_name_short <- gsub("[aeiouAEIOU]", "", sheet_name)

      # If still too long, truncate to 31 characters
      if (nchar(sheet_name_short) > 31) {
        sheet_name_short <- substr(sheet_name_short, 1, 31)
      }

      message(sprintf(
        "Sheet name '%s' was too long. Adjusted to '%s'.",
        original_name, sheet_name_short
      ))

      sheet_name <- sheet_name_short
    }

    # Write to workbook if it's a data frame or tibble
    if (is.data.frame(data) || is_tibble(data)) {
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet = sheet_name, x = data, rowNames = TRUE)
    }
  }

  # Save the file with timestamp
  saveWorkbook(wb, file = file.path(folder_name, paste0(file_name, "_", timestamp, ".xlsx")), overwrite = TRUE)
}

# # Function to save selected results to Excel with row names
# performExport2Excel <- function(
    #     results,
#     file_name = "Results Exports"
#     ) {
#   library(openxlsx)
#
#   # Initialize a workbook
#   wb <- createWorkbook()
#
#   # Loop through each element in the list
#   for (name in names(results)) {
#     data <- results[[name]]
#
#     # Check if the element is a data frame
#     if (is.data.frame(data)) { #  || inherits(data, "confusionMatrix")
#       # Add a worksheet with the name
#       addWorksheet(wb, name)
#
#       # Write the data frame to the worksheet, including row names
#       writeData(wb, sheet = name, x = data, rowNames = TRUE)
#     }
#   }
#
#   # Save the workbook to the specified file
#   saveWorkbook(wb, file = paste0(file_name, ".xlsx"), overwrite = TRUE)
# }


# Usage
# performExport2Excel(myregression  , "Regression Results2")



# Other codes but may not be usable.

# # Function to save selected results to Excel with row names and create a uniquely named folder
# performExport2Excel <- function(
    #     results,
#     folder_name = "Results_Folder",  # Base folder name
#     file_name = "Results Exports"    # Excel file name
# ) {
#   library(openxlsx)
#
#   # Generate a unique folder name based on the current date and time
#   timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
#   unique_folder_name <- paste0(folder_name, "_", timestamp)
#
#   # Create the folder (including the timestamp) if it doesn't exist
#   dir.create(unique_folder_name)
#
#   # Initialize a workbook
#   wb <- createWorkbook()
#
#   # Loop through each element in the list
#   for (name in names(results)) {
#     data <- results[[name]]
#
#     # Check if the element is a data frame
#     if (is.data.frame(data)) {
#       # Add a worksheet with the name
#       addWorksheet(wb, name)
#
#       # Write the data frame to the worksheet, including row names
#       writeData(wb, sheet = name, x = data, rowNames = TRUE)
#     }
#   }
#
#   # Save the workbook to the new folder with the specified file name
#   saveWorkbook(wb, file = file.path(unique_folder_name, paste0(file_name, ".xlsx")), overwrite = TRUE)
# }


# # Function to save selected results to Excel with row names and create a uniquely named file
# performExport2Excel <- function(
    #     results,
#     folder_name = "Results_Folder",  # Base folder name
#     file_name = "Results Exports"    # Excel file name
# ) {
#   library(openxlsx)
#
#   # Generate a timestamp for the file name
#   timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
#
#   # Ensure the folder exists
#   if (!dir.exists(folder_name)) {
#     dir.create(folder_name)
#   }
#
#   # Initialize a workbook
#   wb <- createWorkbook()
#
#   # Loop through each element in the list
#   for (name in names(results)) {
#     data <- results[[name]]
#
#     # Check if the element is a data frame
#     if (is.data.frame(data)) {
#       # Add a worksheet with the name
#       addWorksheet(wb, name)
#
#       # Write the data frame to the worksheet, including row names
#       writeData(wb, sheet = name, x = data, rowNames = TRUE)
#     }
#   }
#
#   # Save the workbook with the timestamped file name
#   saveWorkbook(wb, file = file.path(folder_name, paste0(file_name, "_", timestamp, ".xlsx")), overwrite = TRUE)
# }


# performExport2Excel <- function(
#     results,
#     folder_name = "Results_Folder",
#     file_name   = "Results_Export"
# ) {
#   library(openxlsx)
#
#   # Helper function to sanitize sheet names
#   sanitize_sheet_name <- function(name) {
#     # Remove vowels first
#     name_short <- gsub("[aeiouAEIOU]", "", name)
#     # Truncate if still too long
#     if (nchar(name_short) > 31) {
#       name_short <- substr(name_short, 1, 31)
#     }
#     # Inform if name was changed
#     if (nchar(name) > 31) {
#       message(sprintf("Sheet name '%s' was too long. Adjusted to '%s'.", name, name_short))
#     }
#     return(name_short)
#   }
#
#   # Recursive function to flatten nested lists and return a list of data frames with unique names
#   flatten_results <- function(x, prefix = NULL) {
#     out <- list()
#     for (name in names(x)) {
#       current <- x[[name]]
#       full_name <- if (is.null(prefix)) name else paste(prefix, name, sep = "_")
#       if (is.data.frame(current)) {
#         out[[full_name]] <- current
#       } else if (is.list(current)) {
#         out <- c(out, flatten_results(current, prefix = full_name))
#       } else {
#         message(sprintf("Skipping non-data.frame element: %s", full_name))
#       }
#     }
#     return(out)
#   }
#
#   # Create the output folder if it doesn't exist
#   if (!dir.exists(folder_name)) {
#     dir.create(folder_name)
#   }
#
#   # Timestamp for filename
#   timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
#
#   # Flatten results and sanitize sheet names
#   flat_results <- flatten_results(results)
#
#   # Create workbook and add sheets
#   wb <- createWorkbook()
#   for (sheet in names(flat_results)) {
#     data <- flat_results[[sheet]]
#     sheet_name <- sanitize_sheet_name(sheet)
#     addWorksheet(wb, sheet_name)
#     writeData(wb, sheet = sheet_name, x = data, rowNames = TRUE)
#   }
#
#   # Save workbook
#   saveWorkbook(wb, file = file.path(folder_name, paste0(file_name, "_", timestamp, ".xlsx")), overwrite = TRUE)
# }

