#' Import an Excel (.xlsx or .csv, or .tsv) File or a Text (.txt) File
#'
#' @description
#' This function imports Excel or text files and performs quality controls prior to usage in the
#' succeeding functions. This is ideally used to check data requirements such as uniqueness of
#' inputs in the "Sample" row, uniqueness of Features/Metabolites, data types, and the
#' overall characteristics of the data.
#'
#' @param file_location String. A string of the location of the data in the device. Must have the following characteristics:
#'   \itemize{
#'     \item 1st row name is "Sample": The sample names, ideally shortened versions to make it easier to view in visualizations; without spaces.
#'     \item 2nd row name is "SubjectID": Must be numeric to be sorted correctly. Can contain non-unique identifiers of "Sample" names, usually used when "Replicate" is present. Also referred to as the individuals in the data.
#'     \item 3rd row name is "Replicate": The replicate IDs of Samples. Ideally non-unique identifiers.
#'     \item 4th row name is "Group": Contains QC and the groups (with/out spaces).
#'     \item 5th row name is "Batch": The batch numbers.
#'     \item 6th row name is "Injection": The injection sequence.
#'     \item 7th row name is "DilutionMarker": Concentration markers to be used during data normalization to correct for variations in e.g., urine, concentration between samples. Can be osmolality values, specific gravity, etc.
#'     \item 8th row name is "Response": Numeric. Ideally used when the the response variable (dependent variable) to be used is not "Group".
#'     \item Others requirement/s are:
#'     \itemize{
#'       \item The 1st 8 rows above must be present.
#'       \item After the 8th row (Response), are the features possibly in a mass-to-charge ratio (m/z) and retention time format. Example `199.101@0.111` assuming you rounded off the decimal places to 3 decimal places (suggested for a cleaner output in visualizations).
#'       \item If the DilutionMarker values, SubjectID, Replicate, and Response are not present, leave the rows blank except for the words in the 1st column (e.g., "DilutionMarker' ).
#'       \item Missing values must be left blank or encoded as 0. Blanks and 0s are identified as NAs or missing values.
#'       }
#'     }
#' @param sheet_name String. The name of the worksheet where the data is located. Ignored when data type is ".csv".
#' @param skip_rows Numeric. The number of rows in the data to skip reading.
#' @param separator String. The separator in the ".csv" or ".txt" file. If this fails, it's likely that the file is supposed to be delimited, which should be saved as ".csv", ".tsv", or ".txt".
#'   \itemize{
#'     \item ",": A comma. From a comma-separated file (csv).
#'     \item "/t": A tab From a tab-separated file (.txt or .tsv).
#'     \item Also accepts other delimiters, as long as it is consistent in your data.
#'     }
#'     Defaults to ",". Ignored when the file is ".xlsx".
#' @returns A data frame ready to be used in the `performPreprocessingPeakData` function in the `MetaboStatR` package.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using the default parameters
#' # The location typically looks like
#' # "C:/Users/Documents/Metabolomics Data Analysis/file_name_of_data.xlsx"
#' performDataQualityCheck(
#'  file_location = the_complete_location_of_your_data
#' )
#' }
#'
performDataQualityCheck <- function(
    file_location = NULL,
    sheet_name    = NULL,
    skip_rows     = 0,
    separator     = ","
) {

  importresults <- list()
  importresults$FunctionOrigin <- "performDataQualityCheck"
  importresults$FileLocation <- file_location
  importresults$SheetName <- sheet_name

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------- Load the data -----------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # If file_location is NULL, then choose a file.
  if (is.null(file_location)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      file_location <- rstudioapi::selectFile(
        caption = "Select a file"
        # ,filter = "CSV Files (*.csv)"
      )
    } else {
      # This part ideally should not be used, but just in case the above fails.
      file_location <- file.choose()
    }
  }

  # Check if file exists
  if (!file.exists(file_location)) {
    stop("File does not exist: ", file_location)
  }

  # Extract file extension
  file_ext <- tolower(tools::file_ext(file_location))

  # Read data based on file type
  raw_data <- base::suppressMessages(
    switch(file_ext,
           "xlsx" = {
             # Read Excel file
             readxl::read_excel(
               file_location,
               sheet = sheet_name,
               range = NULL,
               col_names = FALSE,
               skip = skip_rows
             )
           },
           "csv" = {
             # Read CSV file
             utils::read.csv(
               file_location,
               header = FALSE,
               sep = separator,
               skip = skip_rows
             )
           },
           "txt" = {
             # Read text file (assuming tab-delimited by default)
             utils::read.delim(
               file_location,
               header = FALSE,
               sep = separator,
               skip = skip_rows
             )
           },
           "tsv" = {
             # Read text file (assuming tab-delimited by default)
             utils::read.delim(
               file_location,
               header = FALSE,
               sep = separator,
               skip = skip_rows
             )
           },
           {
             # Unsupported file type
             stop("Unsupported file type: ", file_ext,
                  ". Supported types are: xlsx, csv, txt")
           }
    )
  )

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # ------------------------------------ Check if required rows are present -------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Transpose the data
  raw_data <- as.data.frame(t(raw_data)) %>%
    setNames(as.character(.[1, ])) %>% .[-1, ] %>% # Set column names from 1st row (features, etc.)
    # Arrange by Injection Sequence
    dplyr::mutate(Injection = as.numeric(Injection)) %>% # Convert to numeric so it will be arranged correctly in the line below
    dplyr::arrange(Injection) # Sort by injection, needed in drift-correction

  # Data checking steps. Stop if these columns do not exist in the data
  required_headers <- c("Sample", "SubjectID", "Replicate", "Group", "Batch", "Injection", "DilutionMarker", "Response")
  if (!all(required_headers %in% colnames(raw_data))) {
    stop("⚠️ Data preprocessing is halted. The raw data does not contain the expected column names (Sample, SubjectID, Replicate, Group, Batch, Injection, DilutionMarker, Response). Check for typo or the 1st 8 rows in your '.csv' file.")
  }

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------- Check if there are duplicates in some rows ---------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Function to check for duplicates
  check_duplicates <- function(data, data_name) {
    if (length(unique(data)) != length(data)) {
      duplicated_values <- unique(data[duplicated(data) | duplicated(data, fromLast = TRUE)])
      stop(paste0("There are duplicate ", data_name, " in your data. ", data_name, " must be unique to proceed. Here are the duplicated ", data_name, ": ",
                  paste(duplicated_values, collapse = ", ")))
    }
  }

  # Apply checks for sample names and injection orders
  check_duplicates(raw_data$Sample, "Sample names")
  check_duplicates(raw_data$Injection, "Injection orders")
  # Apply checks for features/metabolites
  all_column_names <- names(raw_data) # Get column names excluding required_headers and check for duplicates
  remaining_column_names <- all_column_names[!all_column_names %in% required_headers]
  check_duplicates(remaining_column_names, "Features/metabolites")

  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------- Check if there are not numbers in features/metabolites ---------------------------
  # -------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------

  # Metadata
  metadata <- raw_data %>%
    dplyr::select(c(required_headers))

  # Remove columns that are not part of the features/metabolites
  raw_data <- raw_data %>%
    dplyr::select(-c(required_headers))

  # Check for numeric-ness of data (updated to handle conversion)
  check_numeric <- function(data, data_name) {
    numeric_data <- data %>% as.data.frame()
    numeric_data[]        <- base::lapply(data, as.numeric) # Convert to numeric

    # Check if conversion introduced NAs (meaning non-convertible values)
    if (any(is.na(numeric_data) & !is.na(data))) {
      non_convertible <- unique(data[is.na(numeric_data) & !is.na(data)])
      stop(paste0("There are non-numeric values in your '", data_name, "' data that cannot be converted to numbers. Here are the non-convertible values: ",
                  paste(non_convertible, collapse = ", ")))
    }

    return(numeric_data)  # Return the converted data
  }

  # Combine the metadata and features/metabolites after all preprocessing
  raw_data <- cbind(metadata, raw_data) %>%
    t() %>%
    base::cbind(rownames(.), .) %>%
    `rownames<-`(NULL) %>%
    as.data.frame()

  importresults$raw_data <- raw_data

  return(importresults)
}
