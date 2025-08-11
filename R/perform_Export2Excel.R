#' Perform Export Data Frames to Excel
#'
#' @description
#' This function exports the data frames from a list to Excel. The data frames are detected automatically. All data frames
#' in a list are combined in one (1) Excel file.
#'
#' @param results List. A list of results from `perform_PreprocessingPeakData`, `perform_DimensionReduction`, or `perform_ComparativeAnalysis`, etc. There must only be one at a time. The list may or may not contain a data frame.
#' @param folder_name String. Creates a folder in the current working directory where the Excel file will be exported.
#' @param file_name String. The file name of the results.
#'
#' @returns An Excel file.
#' @export
#'
#' @examples
#' \dontrun{
#' perform_Export2Excel(
#'   results     = results_from_perform_Regression_function,
#'   folder_name = "Regression Results Folder",
#'   file_name   = "Regression_Results_Filename_no_extension")
#'}
#'
perform_Export2Excel <- function(
    results,
    folder_name = "Results_Folder",
    file_name = "Results Exports"
) {

  # # NOTE: THIS WORKS BUT FOR RESULTS THAT ARE NOT IN A LIST

  # Generate a timestamp for the file name
  timestamp <- base::format(base::Sys.time(), "%Y-%m-%d_%H-%M-%S")

  # Create folder if it doesn't exist
  if (!base::dir.exists(folder_name)) {
    base::dir.create(folder_name)
  }

  # Initialize a workbook
  wb <- openxlsx::createWorkbook()

  # Loop through each element in the list
  for (original_name in base::names(results)) {
    data       <- results[[original_name]]
    sheet_name <- original_name

    # Check if the name is too long
    if (base::nchar(sheet_name) > 31) {
      # Remove vowels (case-insensitive)
      sheet_name_short <- base::gsub("[aeiouAEIOU]", "", sheet_name)

      # If still too long, truncate to 31 characters
      if (base::nchar(sheet_name_short) > 31) {
        sheet_name_short <- base::substr(sheet_name_short, 1, 31)
      }

      message(base::sprintf(
        "Sheet name '%s' was too long. Adjusted to '%s'.",
        original_name, sheet_name_short
      ))

      sheet_name <- sheet_name_short
    }

    # Write to workbook if it's a data frame or tibble
    if (base::is.data.frame(data) || tibble::is_tibble(data)) {
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = data, rowNames = TRUE)
    }
  }

  # Save the file with timestamp
  openxlsx::saveWorkbook(wb, file = base::file.path(folder_name, base::paste0(file_name, "_", timestamp, ".xlsx")), overwrite = TRUE)
}
