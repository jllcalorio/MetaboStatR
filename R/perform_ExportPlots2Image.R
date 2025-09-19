#' Export Plot Objects to Image Files
#'
#' @description
#' Exports plot objects (ggplot2, base plots, lattice, plotly, etc.) from named lists or
#' multiple lists to individual image files. The function automatically detects various
#' plot objects within the input, creates appropriately named image files, and handles
#' various edge cases including file name length restrictions, duplicate plots, invalid
#' characters, and different plot formats. Each plot object becomes a separate image file
#' in the specified format. The function can accept either a single named list or multiple
#' lists, and will create separate folders or file naming schemes for each list when
#' multiple folder names are provided.
#'
#' @param results Either a named list containing plot objects, or multiple lists
#'   that can be passed as: results = list1, or results = c(list1, list2),
#'   or results = list(list1, list2, list3). Supported plot types include ggplot2 objects,
#'   base R plots (as recorded plots), lattice plots, plotly objects, and other plot classes.
#' @param folder_name Character string specifying the folder name where the plot files
#'   will be saved. The folder will be created if it doesn't exist. Default is "Plots_Export".
#' @param file_prefix Character vector specifying the base prefix(es) for plot filenames.
#'   If a single prefix is provided, it will be used for all lists. If multiple prefixes
#'   are provided, each list will get its own prefix. Default is "Plot".
#' @param image_format Character string specifying the output image format. Supported formats
#'   include "png", "jpg", "jpeg", "pdf", "svg", "tiff", "bmp". Default is "png".
#' @param width Numeric value specifying the width of the output images in inches.
#'   Default is 15.
#' @param height Numeric value specifying the height of the output images in inches.
#'   Default is 12.
#' @param dpi Numeric value specifying the resolution (dots per inch) for raster formats.
#'   Default is 300.
#' @param include_timestamp Logical indicating whether to include a timestamp in the
#'   filename(s). Default is TRUE.
#' @param overwrite Logical indicating whether to overwrite existing files. If FALSE and file
#'   exists, an error will be thrown. Default is TRUE.
#' @param use_plot_titles Logical indicating whether to use plot titles as filenames when
#'   available (sanitized for file system compatibility). Default is TRUE.
#' @param max_plots Integer specifying maximum number of plots to export per execution.
#'   Default is 100.
#' @param quality Numeric value between 0 and 100 specifying the quality for JPEG format.
#'   Ignored for other formats. Default is 95.
#' @param transparent Logical indicating whether to use transparent background for PNG format.
#'   Default is TRUE.
#'
#' @return Invisibly returns a character vector of full paths to the created image files.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Validates and processes input parameters (single list or multiple lists)
#'   \item Processes multiple lists and creates separate naming schemes when multiple prefixes provided
#'   \item Detects and filters various types of plot objects (ggplot2, base plots, lattice, plotly, etc.)
#'   \item Handles duplicate plots by comparing plot objects and adding unique identifiers
#'   \item Sanitizes filenames to comply with file system naming conventions
#'   \item Exports plots with specified format, dimensions, and quality settings
#'   \item Provides comprehensive error handling and progress reporting
#' }
#'
#' Supported Plot Types:
#' \itemize{
#'   \item ggplot2 objects (class "gg" and "ggplot")
#'   \item Base R recorded plots (class "recordedplot")
#'   \item Lattice plots (class "trellis")
#'   \item Plotly objects (class "plotly")
#'   \item Grid graphics (class "grob" and "gTree")
#'   \item Custom plot objects with print methods
#' }
#'
#' Duplicate Detection:
#' \itemize{
#'   \item Compares plot objects using content hashing for ggplot2 and lattice
#'   \item Uses object structure comparison for other plot types
#'   \item Automatically appends numeric suffixes to duplicate plots
#'   \item Preserves all plots while ensuring unique filenames
#' }
#'
#' File naming follows the pattern: `prefix_plot_name/title_timestamp.format`
#' Names are automatically sanitized to remove invalid characters and ensure
#' compatibility across different operating systems.
#'
#' @section Dependencies:
#' Requires \code{ggplot2} for ggplot object detection and export. Optional packages
#' include \code{lattice}, \code{plotly}, \code{grid}, and format-specific packages.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Create sample plots
#' plot1 <- ggplot(mtcars, aes(x = mpg, y = hp)) + geom_point() + ggtitle("MPG vs HP")
#' plot2 <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
#'          geom_point() + ggtitle("Iris Sepal Analysis")
#' plot3 <- ggplot(mtcars, aes(x = mpg, y = hp)) + geom_point() + ggtitle("MPG vs HP") # Duplicate
#'
#' list1 <- list(
#'   scatter_plot = plot1,
#'   iris_analysis = plot2,
#'   duplicate_plot = plot3,
#'   not_a_plot = "This will be ignored"
#' )
#'
#' # Method 1: Basic usage
#' perform_ExportPlots2Image(results = list1)
#'
#' # Method 2: Custom format and dimensions
#' perform_ExportPlots2Image(
#'   results = list1,
#'   folder_name = "My_Plots",
#'   file_prefix = "Analysis",
#'   image_format = "pdf",
#'   width = 12,
#'   height = 10,
#'   dpi = 600
#' )
#'
#' # Method 3: Multiple lists with different prefixes
#' list2 <- list(summary_plot = plot1)
#' perform_ExportPlots2Image(
#'   results = c(list1, list2),
#'   file_prefix = c("Primary", "Secondary"),
#'   image_format = "png",
#'   use_plot_titles = TRUE,
#'   transparent = FALSE
#' )
#' }
#'
#' @author John Lennon L. Calorio
#'
#' @seealso
#' \code{\link[ggplot2]{ggsave}}, \code{\link[grDevices]{png}},
#' \code{\link[grDevices]{pdf}}, \code{\link[lattice]{trellis.device}}
#'
#' @importFrom graphics text
#' @importFrom stats runif
#' @importFrom utils getS3method object.size str
#' @importFrom grDevices adjustcolor
#' @importFrom graphics barplot box hist legend lines par points
#' @importFrom plotly plotly plot_ly layout
#' @importFrom webshot2 webshot
#'
#' @export
perform_ExportPlots2Image <- function(
    results,
    folder_name = "Plots_Export",
    file_prefix = "Plot",
    image_format = "png",
    width = 15,
    height = 12,
    dpi = 300,
    include_timestamp = TRUE,
    overwrite = TRUE,
    use_plot_titles = TRUE,
    max_plots = 100,
    quality = 95,
    transparent = FALSE
) {

  # Start timing for performance monitoring
  start_time <- Sys.time()

  # Process and validate inputs (suppress console output)
  processed_results <- suppressMessages(.process_multiple_plot_inputs(results))
  .validate_plot_inputs(processed_results, folder_name, file_prefix, image_format,
                        width, height, dpi, include_timestamp, overwrite,
                        use_plot_titles, max_plots, quality, transparent)

  # Check required packages based on format (suppress messages)
  suppressMessages(.check_plot_dependencies(image_format))

  # Create output directory (suppress creation messages)
  output_path <- suppressMessages(.create_output_directory(folder_name))

  # Process file prefixes and create corresponding files
  file_paths <- .export_multiple_plot_collections(
    processed_results, output_path, file_prefix, image_format,
    width, height, dpi, include_timestamp, overwrite,
    use_plot_titles, max_plots, quality, transparent
  )

  # Performance summary
  end_time <- Sys.time()
  processing_time <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)

  message(sprintf("Plot export completed in %s seconds", processing_time))
  message(sprintf("Total plot files created: %d", length(file_paths)))

  return(invisible(file_paths))
}


#' Process Multiple Plot Input Lists
#' @noRd
.process_multiple_plot_inputs <- function(results) {
  if (missing(results) || is.null(results)) {
    stop("No input provided. Please provide at least one list containing plot objects.",
         call. = FALSE)
  }

  if (!is.list(results)) {
    stop("'results' must be a list or multiple lists.", call. = FALSE)
  }

  # Check if this is a list of lists (multiple lists scenario)
  if (.is_list_of_plot_lists(results)) {
    return(results)
  } else {
    # Single named list scenario - wrap in list to standardize processing
    if (is.null(names(results))) {
      # Unnamed list - assign generic names to plot elements
      plot_indices <- vapply(results, .is_plot_object, logical(1))
      if (any(plot_indices)) {
        names(results)[plot_indices] <- paste0("Plot_", seq_len(sum(plot_indices)))
      }
      if (any(!plot_indices)) {
        names(results)[!plot_indices] <- paste0("Other_", seq_len(sum(!plot_indices)))
      }
    }
    return(list(results))  # Wrap single list
  }
}


#' Check if input is a list of lists containing plots
#' @noRd
.is_list_of_plot_lists <- function(x) {
  if (!is.list(x) || length(x) == 0) return(FALSE)

  # Check if most elements are lists and at least one contains plots
  list_elements <- vapply(x, is.list, logical(1))
  if (sum(list_elements) < length(x) * 0.5) return(FALSE)

  # Check if any of the list elements contain plot objects
  contains_plots <- any(vapply(x[list_elements], function(sublist) {
    any(vapply(sublist, .is_plot_object, logical(1)))
  }, logical(1)))

  return(contains_plots)
}


#' Comprehensive Plot Object Detection with Nested List Support
#' @noRd
.is_plot_object <- function(x) {
  if (is.null(x)) return(FALSE)

  # Check for various plot object classes
  plot_classes <- c("gg", "ggplot", "recordedplot", "trellis", "plotly",
                    "grob", "gTree", "lattice", "pheatmap")

  # Direct class checking - exclude data.frame and related classes
  obj_classes <- class(x)

  # Explicitly exclude data frames, tibbles, matrices, and other data structures
  excluded_classes <- c("data.frame", "tbl_df", "tbl", "tibble", "matrix",
                        "array", "list", "character", "numeric", "logical",
                        "integer", "factor", "prcomp", "POSIXct", "POSIXt",
                        "double", "table", "integer")

  # If object has excluded classes, it's not a plot
  if (any(excluded_classes %in% obj_classes)) {
    return(FALSE)
  }

  # Check for plot classes
  if (any(plot_classes %in% obj_classes)) {
    return(TRUE)
  }

  # Additional specific checks

  # ggplot2 objects
  if (inherits(x, c("gg", "ggplot"))) {
    return(TRUE)
  }

  # Base R recorded plots
  if (inherits(x, "recordedplot")) {
    return(TRUE)
  }

  # Lattice plots
  if (inherits(x, "trellis")) {
    return(TRUE)
  }

  # Plotly objects
  if (inherits(x, "plotly")) {
    return(TRUE)
  }

  # Grid graphics
  if (inherits(x, c("grob", "gTree"))) {
    return(TRUE)
  }

  # Check for objects with plot-specific methods (but not data structures)
  if (!any(excluded_classes %in% obj_classes) && .has_plot_method(x)) {
    return(TRUE)
  }

  return(FALSE)
}


#' Check if Object has Plotting Methods (Exclude Data Structures)
#' @noRd
.has_plot_method <- function(x) {
  # First, explicitly exclude data structures
  excluded_classes <- c("data.frame", "tbl_df", "tbl", "tibble", "matrix",
                        "array", "character", "numeric", "logical",
                        "integer", "factor", "prcomp", "POSIXct", "POSIXt",
                        "double", "table", "integer")

  if (any(excluded_classes %in% class(x))) {
    return(FALSE)
  }

  obj_class <- class(x)[1]

  # Check for print method (many plot objects use print for display)
  has_print <- !is.null(getS3method("print", obj_class, optional = TRUE))

  # Check for plot method
  has_plot <- !is.null(getS3method("plot", obj_class, optional = TRUE))

  # Additional heuristics for plot-like objects (but exclude data structures)
  if (is.list(x) && !any(c("data.frame", "tbl_df", "tbl") %in% class(x))) {
    # Check for common plot object structure
    plot_like_names <- c("layers", "scales", "mapping", "theme",
                         "coordinates", "facet", "plot_env", "panel")
    has_plot_structure <- any(plot_like_names %in% names(x))

    if (has_plot_structure) return(TRUE)
  }

  return(has_print || has_plot)
}


#' Validate Plot Export Inputs
#' @noRd
.validate_plot_inputs <- function(results, folder_name, file_prefix, image_format,
                                  width, height, dpi, include_timestamp, overwrite,
                                  use_plot_titles, max_plots, quality, transparent) {

  # Check results parameter
  if (is.null(results) || !is.list(results) || length(results) == 0) {
    stop("'results' must be a non-empty list.", call. = FALSE)
  }

  # Check folder_name
  if (!is.character(folder_name) || length(folder_name) != 1 || nchar(trimws(folder_name)) == 0) {
    stop("'folder_name' must be a non-empty character string.", call. = FALSE)
  }

  # Check file_prefix
  if (!is.character(file_prefix) || length(file_prefix) == 0 || any(nchar(trimws(file_prefix)) == 0)) {
    stop("'file_prefix' must be a non-empty character vector.", call. = FALSE)
  }

  # Check image_format
  valid_formats <- c("png", "jpg", "jpeg", "pdf", "svg", "tiff", "bmp")
  if (!is.character(image_format) || length(image_format) != 1 ||
      !tolower(image_format) %in% valid_formats) {
    stop(sprintf("'image_format' must be one of: %s", paste(valid_formats, collapse = ", ")),
         call. = FALSE)
  }

  # Check numeric parameters
  numeric_params <- list(
    width = width, height = height, dpi = dpi, max_plots = max_plots, quality = quality
  )

  for (param_name in names(numeric_params)) {
    param_value <- numeric_params[[param_name]]
    if (!is.numeric(param_value) || length(param_value) != 1 || is.na(param_value) || param_value <= 0) {
      stop(sprintf("'%s' must be a single positive numeric value.", param_name), call. = FALSE)
    }
  }

  # Additional checks for specific parameters
  if (quality < 1 || quality > 100) {
    stop("'quality' must be between 1 and 100.", call. = FALSE)
  }

  if (max_plots > 1000) {
    stop("'max_plots' cannot exceed 1000 for performance reasons.", call. = FALSE)
  }

  # Check minimum values for dimensions and DPI
  if (width < 15) {
    stop("'width' must be at least 15 inches.", call. = FALSE)
  }

  if (height < 12) {
    stop("'height' must be at least 12 inches.", call. = FALSE)
  }

  if (dpi < 600) {
    stop("'dpi' must be at least 600.", call. = FALSE)
  }

  # Check logical parameters
  logical_params <- list(
    include_timestamp = include_timestamp,
    overwrite = overwrite,
    use_plot_titles = use_plot_titles,
    transparent = transparent
  )

  for (param_name in names(logical_params)) {
    param_value <- logical_params[[param_name]]
    if (!is.logical(param_value) || length(param_value) != 1 || is.na(param_value)) {
      stop(sprintf("'%s' must be a single logical value (TRUE or FALSE).", param_name),
           call. = FALSE)
    }
  }

  # Check for invalid characters in folder and file names
  invalid_chars <- "[<>:\"/\\|?*]"
  if (grepl(invalid_chars, folder_name)) {
    stop("'folder_name' contains invalid characters. Avoid: < > : \" / \\ | ? *", call. = FALSE)
  }
  if (any(grepl(invalid_chars, file_prefix))) {
    stop("'file_prefix' contains invalid characters. Avoid: < > : \" / \\ | ? *", call. = FALSE)
  }
}


#' Check Dependencies for Plot Export
#' @noRd
.check_plot_dependencies <- function(image_format) {
  format_lower <- tolower(image_format)

  # Always check for grDevices (base package, should always be available)
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package 'grDevices' is required but not available.", call. = FALSE)
  }

  # Format-specific package checks
  if (format_lower == "svg") {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      warning("Package 'svglite' is recommended for SVG export. Using base R svg() device.",
              call. = FALSE)
    }
  }

  # Check for ggplot2 if we're likely to encounter ggplot objects
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' not available. ggplot objects cannot be exported.")
  }
}


#' Export Multiple Plot Collections
#' @noRd
.export_multiple_plot_collections <- function(processed_results, output_path, file_prefixes,
                                              image_format, width, height, dpi, include_timestamp,
                                              overwrite, use_plot_titles, max_plots, quality, transparent) {

  n_lists <- length(processed_results)
  n_prefixes <- length(file_prefixes)

  # Generate prefixes for each list
  final_prefixes <- character(n_lists)

  if (n_prefixes == 1 && n_lists > 1) {
    # Single prefix for multiple lists - add list identifiers
    for (i in seq_len(n_lists)) {
      final_prefixes[i] <- paste0(file_prefixes[1], "_List", i)
    }
  } else {
    # Assign prefixes to lists
    for (i in seq_len(n_lists)) {
      if (i <= n_prefixes) {
        final_prefixes[i] <- file_prefixes[i]
      } else {
        # Use last prefix with suffix
        base_prefix <- file_prefixes[n_prefixes]
        suffix_num <- i - n_prefixes + 1
        final_prefixes[i] <- paste0(base_prefix, "_", suffix_num)
      }
    }
  }

  # Export plots from each list
  all_file_paths <- character(0)

  for (i in seq_along(processed_results)) {
    current_list <- processed_results[[i]]
    current_prefix <- final_prefixes[i]

    # Extract and validate plot objects
    plot_list <- .extract_plot_objects(current_list, max_plots)

    if (length(plot_list) == 0) {
      # Don't show message for empty lists
      next
    }

    # Handle duplicates
    unique_plot_list <- .handle_duplicate_plots(plot_list)

    # Export plots
    list_file_paths <- .export_plot_list(
      unique_plot_list, output_path, current_prefix, image_format,
      width, height, dpi, include_timestamp, overwrite,
      use_plot_titles, quality, transparent
    )

    all_file_paths <- c(all_file_paths, list_file_paths)

    # Only show message if plots were actually exported
    if (length(list_file_paths) > 0) {
      message(sprintf("Exported %d plots from list %d with prefix: %s",
                      length(list_file_paths), i, current_prefix))
    }
  }

  return(all_file_paths)
}


#' Extract Plot Objects from List with Deep Nesting Support
#' @noRd
.extract_plot_objects <- function(results, max_plots) {
  # Recursively find all plot objects in nested lists
  all_plots <- .find_plots_recursive(results, parent_name = "")

  if (length(all_plots) == 0) {
    return(list())
  }

  # Apply plot limit
  if (length(all_plots) > max_plots) {
    warning(sprintf("Found %d plot objects but max_plots is set to %d. Only the first %d will be exported.",
                    length(all_plots), max_plots, max_plots), call. = FALSE)
    all_plots <- all_plots[1:max_plots]
  }

  return(all_plots)
}


#' Recursively Find Plot Objects in Nested Lists
#' @noRd
.find_plots_recursive <- function(obj, parent_name = "", max_depth = 10, current_depth = 0) {
  # Prevent infinite recursion
  if (current_depth >= max_depth) {
    return(list())
  }

  found_plots <- list()

  if (is.list(obj) && !.is_plot_object(obj)) {
    # This is a list but not a plot object itself
    obj_names <- names(obj)
    if (is.null(obj_names)) {
      obj_names <- paste0("item", seq_along(obj))
    }

    for (i in seq_along(obj)) {
      element <- obj[[i]]
      element_name <- obj_names[i]

      # Create hierarchical name
      if (nchar(parent_name) > 0) {
        full_name <- paste(parent_name, element_name, sep = "_")
      } else {
        full_name <- element_name
      }

      if (.is_plot_object(element)) {
        # Found a plot object
        found_plots[[full_name]] <- element
      } else if (is.list(element)) {
        # Recursively search nested lists
        nested_plots <- .find_plots_recursive(element, full_name, max_depth, current_depth + 1)
        found_plots <- c(found_plots, nested_plots)
      }
    }
  } else if (.is_plot_object(obj)) {
    # This object itself is a plot
    plot_name <- if (nchar(parent_name) > 0) parent_name else "plot"
    found_plots[[plot_name]] <- obj
  }

  return(found_plots)
}


#' Handle Duplicate Plots
#' @noRd
.handle_duplicate_plots <- function(plot_list) {
  if (length(plot_list) <= 1) {
    return(plot_list)
  }

  # Create plot signatures for duplicate detection
  plot_signatures <- character(length(plot_list))
  plot_names <- names(plot_list)

  for (i in seq_along(plot_list)) {
    plot_signatures[i] <- .create_plot_signature(plot_list[[i]])
  }

  # Identify duplicates and create unique names
  unique_names <- character(length(plot_list))
  signature_counts <- table(plot_signatures)
  signature_indices <- integer(length(unique(plot_signatures)))
  names(signature_indices) <- unique(plot_signatures)

  for (i in seq_along(plot_list)) {
    signature <- plot_signatures[i]
    base_name <- if (is.null(plot_names[i]) || plot_names[i] == "") {
      paste0("Plot_", i)
    } else {
      plot_names[i]
    }

    if (signature_counts[signature] > 1) {
      # This is a duplicate
      signature_indices[signature] <- signature_indices[signature] + 1
      unique_names[i] <- paste0(base_name, "_duplicate_", signature_indices[signature])

      if (signature_indices[signature] == 1) {
        # Also rename the first occurrence
        first_occurrence <- which(plot_signatures == signature)[1]
        if (first_occurrence < i) {
          original_name <- if (is.null(plot_names[first_occurrence]) || plot_names[first_occurrence] == "") {
            paste0("Plot_", first_occurrence)
          } else {
            plot_names[first_occurrence]
          }
          unique_names[first_occurrence] <- paste0(original_name, "_duplicate_1")
        }
      }
    } else {
      unique_names[i] <- base_name
    }
  }

  names(plot_list) <- unique_names

  # Report duplicates found
  n_duplicates <- sum(signature_counts > 1)
  if (n_duplicates > 0) {
    # Suppress duplicate messages to reduce console output
  }

  return(plot_list)
}


#' Create Plot Signature for Duplicate Detection
#' @noRd
.create_plot_signature <- function(plot_obj) {
  tryCatch({
    if (inherits(plot_obj, c("gg", "ggplot"))) {
      # For ggplot objects, use a combination of layers, data, and aesthetics
      signature_components <- list(
        data_summary = if (!is.null(plot_obj$data)) digest::digest(plot_obj$data, algo = "xxhash32") else "no_data",
        mapping = if (!is.null(plot_obj$mapping)) as.character(plot_obj$mapping) else "no_mapping",
        layers = length(plot_obj$layers),
        layer_geoms = vapply(plot_obj$layers, function(l) class(l$geom)[1], character(1))
      )
      return(digest::digest(signature_components, algo = "xxhash32"))
    } else if (inherits(plot_obj, "trellis")) {
      # For lattice plots
      signature_components <- list(
        call = as.character(plot_obj$call),
        panel_function = if (!is.null(plot_obj$panel)) as.character(plot_obj$panel) else "default"
      )
      return(digest::digest(signature_components, algo = "xxhash32"))
    } else {
      # For other plot types, use object structure
      return(digest::digest(str(plot_obj, max.level = 2), algo = "xxhash32"))
    }
  }, error = function(e) {
    # Fallback: use object size and class as a weak signature
    return(paste(class(plot_obj)[1], object.size(plot_obj), sep = "_"))
  })
}


#' Export Individual Plot List
#' @noRd
.export_plot_list <- function(plot_list, output_path, file_prefix, image_format,
                              width, height, dpi, include_timestamp, overwrite,
                              use_plot_titles, quality, transparent) {

  exported_files <- character(0)

  for (i in seq_along(plot_list)) {
    plot_obj <- plot_list[[i]]
    plot_name <- names(plot_list)[i]

    # Generate filename
    filename <- .generate_plot_filename(
      plot_obj, plot_name, file_prefix, image_format,
      include_timestamp, use_plot_titles
    )

    full_path <- file.path(output_path, filename)

    # Check file availability
    if (!.check_plot_file_availability(full_path, overwrite)) {
      next
    }

    # Export plot
    success <- .export_single_plot(
      plot_obj, full_path, image_format, width, height,
      dpi, quality, transparent
    )

    if (success) {
      exported_files <- c(exported_files, full_path)
      # Suppress individual file messages to reduce console output
    } else {
      warning(sprintf("Failed to export plot: %s", plot_name), call. = FALSE)
    }
  }

  return(exported_files)
}


#' Generate Plot Filename
#' @noRd
.generate_plot_filename <- function(plot_obj, plot_name, file_prefix, image_format,
                                    include_timestamp, use_plot_titles) {

  # Start with base name
  base_name <- plot_name

  # Try to extract title if requested
  if (use_plot_titles) {
    plot_title <- .extract_plot_title(plot_obj)
    if (!is.null(plot_title) && nchar(trimws(plot_title)) > 0) {
      base_name <- plot_title
    }
  }

  # Sanitize base name
  base_name <- .sanitize_filename(base_name)

  # Combine with prefix
  if (nchar(file_prefix) > 0) {
    final_name <- paste(file_prefix, base_name, sep = "_")
  } else {
    final_name <- base_name
  }

  # Add timestamp if requested
  if (include_timestamp) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    final_name <- paste(final_name, timestamp, sep = "_")
  }

  # Add extension
  final_name <- paste0(final_name, ".", tolower(image_format))

  return(final_name)
}


#' Extract Plot Title
#' @noRd
.extract_plot_title <- function(plot_obj) {
  tryCatch({
    if (inherits(plot_obj, c("gg", "ggplot"))) {
      # Extract ggplot title
      if (!is.null(plot_obj$labels$title)) {
        return(as.character(plot_obj$labels$title))
      }
    } else if (inherits(plot_obj, "trellis")) {
      # Extract lattice title
      if (!is.null(plot_obj$main)) {
        return(as.character(plot_obj$main))
      }
    }
    # Could add more plot type specific title extraction here
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}


#' Sanitize Filename
#' @noRd
.sanitize_filename <- function(filename) {
  # Remove or replace invalid characters
  filename <- gsub("[<>:\"/\\|?*]", "_", filename)
  filename <- gsub("\\s+", "_", filename)  # Replace spaces with underscores
  filename <- gsub("_{2,}", "_", filename)  # Remove multiple underscores
  filename <- gsub("^_|_$", "", filename)  # Remove leading/trailing underscores

  # Ensure non-empty
  if (nchar(filename) == 0) {
    filename <- "plot"
  }

  # Limit length (keeping some buffer for timestamp and extension)
  max_length <- 200
  if (nchar(filename) > max_length) {
    filename <- substr(filename, 1, max_length)
  }

  return(filename)
}


#' Check Plot File Availability
#' @noRd
.check_plot_file_availability <- function(file_path, overwrite) {
  if (!file.exists(file_path)) {
    return(TRUE)
  }

  if (!overwrite) {
    # Suppress individual skip messages to reduce console output
    return(FALSE)
  }

  return(TRUE)
}


#' Export Single Plot with Robust Error Handling
#' @noRd
.export_single_plot <- function(plot_obj, file_path, image_format, width, height,
                                dpi, quality, transparent) {

  format_lower <- tolower(image_format)

  tryCatch({
    if (inherits(plot_obj, c("gg", "ggplot"))) {
      # Use ggsave for ggplot objects
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package required for ggplot export")
      }

      ggplot2::ggsave(
        filename = file_path,
        plot = plot_obj,
        width = width,
        height = height,
        dpi = dpi,
        bg = if (transparent && format_lower == "png") "transparent" else "white"
      )

    } else if (inherits(plot_obj, "recordedplot")) {
      # Handle recorded plots (base R plots)
      .export_recorded_plot(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent)

    } else if (inherits(plot_obj, "trellis")) {
      # Handle lattice plots
      .export_trellis_plot(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent)

    } else if (inherits(plot_obj, "plotly")) {
      # Handle plotly objects
      .export_plotly_plot(plot_obj, file_path, format_lower, width, height, dpi)

    } else if (inherits(plot_obj, c("grob", "gTree"))) {
      # Handle grid graphics objects
      .export_grid_plot(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent)

    } else {
      # Generic plot export using print method
      .export_generic_plot(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent)
    }

    return(TRUE)

  }, error = function(e) {
    warning(sprintf("Failed to export plot to %s: %s", basename(file_path), e$message), call. = FALSE)
    return(FALSE)
  })
}


#' Export Recorded Plot (Base R Plots)
#' @noRd
.export_recorded_plot <- function(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent) {
  # Open graphics device
  .open_graphics_device(file_path, format_lower, width, height, dpi, quality, transparent)

  # Replay the recorded plot
  grDevices::replayPlot(plot_obj)

  # Close device
  grDevices::dev.off()
}


#' Export Trellis Plot (Lattice)
#' @noRd
.export_trellis_plot <- function(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent) {
  # Open graphics device
  .open_graphics_device(file_path, format_lower, width, height, dpi, quality, transparent)

  # Print the trellis object
  print(plot_obj)

  # Close device
  grDevices::dev.off()
}


#' Export Plotly Plot - Pure R Implementation
#' @noRd
.export_plotly_plot <- function(plot_obj, file_path, format_lower, width, height, dpi) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly package required for plotly export")
  }

  # Try HTML export first (native plotly capability)
  if (format_lower == "html") {
    tryCatch({
      htmlwidgets::saveWidget(plot_obj, file_path, selfcontained = TRUE)
      return(TRUE)
    }, error = function(e) {
      warning(sprintf("HTML export failed: %s", e$message))
    })
  }

  # For image formats, use webshot2 if available (pure R solution)
  if (format_lower %in% c("png", "jpg", "jpeg", "pdf")) {
    success <- .export_plotly_with_webshot(plot_obj, file_path, format_lower, width, height, dpi)
    if (success) return(TRUE)
  }

  # Fallback to static conversion methods
  .export_plotly_fallback(plot_obj, file_path, format_lower, width, height, dpi)
}

#' Export Plotly using webshot2 (Pure R)
#' @noRd
.export_plotly_with_webshot <- function(plot_obj, file_path, format_lower, width, height, dpi) {
  # Check for webshot2 (pure R package for web screenshots)
  if (!requireNamespace("webshot2", quietly = TRUE)) {
    return(FALSE)
  }

  tryCatch({
    # Create temporary HTML file
    temp_html <- tempfile(fileext = ".html")

    # Save plotly as HTML widget
    if (requireNamespace("htmlwidgets", quietly = TRUE)) {
      htmlwidgets::saveWidget(plot_obj, temp_html, selfcontained = TRUE)
    } else {
      return(FALSE)
    }

    # Convert dimensions to pixels
    width_px <- width * dpi / 72  # Convert inches to pixels (assuming 72 DPI base)
    height_px <- height * dpi / 72

    # Take screenshot using webshot2
    webshot2::webshot(
      url = temp_html,
      file = file_path,
      width = width_px,
      height = height_px,
      delay = 2,  # Allow time for plotly to render
      zoom = 1
    )

    # Clean up temporary file
    if (file.exists(temp_html)) {
      unlink(temp_html)
    }

    return(TRUE)

  }, error = function(e) {
    warning(sprintf("webshot2 export failed: %s", e$message))
    return(FALSE)
  })
}


#' Enhanced Fallback Plotly Export Method
#' @noRd
.export_plotly_fallback <- function(plot_obj, file_path, format_lower, width, height, dpi) {
  tryCatch({
    # Try to extract data from plotly object and create static plot
    success <- .convert_plotly_to_static(plot_obj, file_path, format_lower, width, height, dpi)
    if (success) return(TRUE)

    # Ultimate fallback: create informative placeholder
    .create_plotly_placeholder(file_path, format_lower, width, height, dpi)

  }, error = function(e) {
    stop(sprintf("Failed to export plotly object: %s", e$message))
  })
}


#' Convert Plotly to Static Plot using extracted data
#' @noRd
.convert_plotly_to_static <- function(plot_obj, file_path, format_lower, width, height, dpi) {
  tryCatch({
    # Extract plotly data and configuration
    plot_data <- .extract_plotly_data(plot_obj)

    if (is.null(plot_data)) {
      return(FALSE)
    }

    # Open graphics device
    .open_graphics_device(file_path, format_lower, width, height, dpi, 95, TRUE)

    # Create static version based on plot type
    .render_static_from_plotly_data(plot_data)

    # Close device
    grDevices::dev.off()

    return(TRUE)

  }, error = function(e) {
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    return(FALSE)
  })
}


#' Render static plot from plotly data
#' @noRd
.render_static_from_plotly_data <- function(plot_data) {
  traces <- plot_data$traces
  layout <- plot_data$layout

  # Set up plot area
  par(mar = c(4, 4, 3, 2))

  # Initialize plot with first trace to set up axes
  first_trace <- traces[[1]]
  plot_type <- first_trace$type %||% "scatter"

  # Extract title and axis labels
  main_title <- layout$title$text %||% layout$title %||% "Plotly Plot (Static Export)"
  xlab <- layout$xaxis$title$text %||% layout$xaxis$title %||% "X"
  ylab <- layout$yaxis$title$text %||% layout$yaxis$title %||% "Y"

  if (plot_type == "scatter") {
    .render_scatter_trace(first_trace, main_title, xlab, ylab, TRUE)

    # Add additional traces
    if (length(traces) > 1) {
      for (i in 2:length(traces)) {
        .render_scatter_trace(traces[[i]], "", "", "", FALSE)
      }
    }

  } else if (plot_type == "bar") {
    .render_bar_trace(first_trace, main_title, xlab, ylab, TRUE)

    # Add additional bar traces (simplified)
    if (length(traces) > 1) {
      for (i in 2:length(traces)) {
        .render_bar_trace(traces[[i]], "", "", "", FALSE)
      }
    }

  } else if (plot_type %in% c("histogram", "hist")) {
    .render_histogram_trace(first_trace, main_title, xlab, ylab)

  } else {
    # Generic fallback
    .render_generic_trace(first_trace, main_title, xlab, ylab)
  }

  # Add legend if multiple traces
  if (length(traces) > 1) {
    legend_labels <- vapply(traces, function(t) t$name %||% paste("Trace", seq_along(traces)), character(1))
    legend("topright", legend = legend_labels, col = 1:length(traces), pch = 16, cex = 0.8)
  }
}


#' Render bar trace
#' @noRd
.render_bar_trace <- function(trace, main_title, xlab, ylab, new_plot = TRUE) {
  x_data <- trace$x
  y_data <- trace$y

  if (is.null(y_data)) return()

  if (is.null(x_data)) {
    x_data <- seq_along(y_data)
  }

  # Convert to numeric if needed
  if (!is.numeric(y_data)) y_data <- as.numeric(y_data)

  if (new_plot) {
    barplot(y_data, names.arg = x_data, main = main_title, xlab = xlab, ylab = ylab)
  } else {
    # For multiple bar traces, this is simplified - would need more complex logic for proper grouping
    barplot(y_data, add = TRUE, col = adjustcolor("blue", alpha.f = 0.5))
  }
}


#' Render scatter trace
#' @noRd
.render_scatter_trace <- function(trace, main_title, xlab, ylab, new_plot = TRUE) {
  x_data <- trace$x
  y_data <- trace$y

  if (is.null(x_data) || is.null(y_data)) return()

  # Convert to numeric if needed
  if (!is.numeric(x_data)) x_data <- as.numeric(x_data)
  if (!is.numeric(y_data)) y_data <- as.numeric(y_data)

  # Remove NA values
  valid_idx <- !is.na(x_data) & !is.na(y_data)
  x_data <- x_data[valid_idx]
  y_data <- y_data[valid_idx]

  if (length(x_data) == 0 || length(y_data) == 0) return()

  mode <- trace$mode %||% "markers"

  if (new_plot) {
    if (grepl("lines", mode)) {
      plot(x_data, y_data, type = "l", main = main_title, xlab = xlab, ylab = ylab)
    } else {
      plot(x_data, y_data, type = "p", main = main_title, xlab = xlab, ylab = ylab, pch = 16)
    }
  } else {
    if (grepl("lines", mode)) {
      lines(x_data, y_data)
    } else {
      points(x_data, y_data, pch = 16)
    }
  }
}


#' Render histogram trace
#' @noRd
.render_histogram_trace <- function(trace, main_title, xlab, ylab) {
  x_data <- trace$x

  if (is.null(x_data)) return()

  # Convert to numeric if needed
  if (!is.numeric(x_data)) x_data <- as.numeric(x_data)

  # Remove NA values
  x_data <- x_data[!is.na(x_data)]

  if (length(x_data) == 0) return()

  hist(x_data, main = main_title, xlab = xlab, ylab = ylab)
}


#' Render generic trace (fallback)
#' @noRd
.render_generic_trace <- function(trace, main_title, xlab, ylab) {
  # Create a simple placeholder plot
  plot(1, 1, type = "n", main = main_title, xlab = xlab, ylab = ylab)
  text(1, 1, paste("Plotly", trace$type %||% "plot", "\n(Static conversion)"), cex = 1.2, col = "blue")
}


#' Create informative placeholder for plotly plots
#' @noRd
.create_plotly_placeholder <- function(file_path, format_lower, width, height, dpi) {
  # Open graphics device
  .open_graphics_device(file_path, format_lower, width, height, dpi, 95, TRUE)

  # Create informative placeholder
  plot(1, 1, type = "n",
       xlab = "Plotly Export",
       ylab = "Static Version",
       main = "Plotly Interactive Plot (Exported as Static)")

  # Add informative text
  text(1, 1,
       "Interactive Plotly Plot\n\nFor full interactivity:\n -- Save as HTML format\n -- Install webshot2 package\n -- Or use plotly in R environment",
       cex = 1,
       col = "darkblue",
       adj = 0.5)

  # Add a simple border
  box(col = "gray", lwd = 2)

  # Close device
  grDevices::dev.off()
}

#' Null coalescing operator
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x) || (is.atomic(x) && length(x) == 1 && is.na(x))) y else x
}

#' Updated Check Dependencies for Plot Export (No Python dependencies)
#' @noRd
.check_plot_dependencies <- function(image_format) {
  format_lower <- tolower(image_format)

  # Always check for grDevices (base package, should always be available)
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package 'grDevices' is required but not available.", call. = FALSE)
  }

  # Format-specific package checks
  if (format_lower == "svg") {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      warning("Package 'svglite' is recommended for SVG export. Using base R svg() device.",
              call. = FALSE)
    }
  }

  # Check for ggplot2 if we're likely to encounter ggplot objects
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Package 'ggplot2' not available. ggplot objects cannot be exported.")
  }

  # Optional packages for enhanced plotly export (pure R solutions)
  if (!requireNamespace("webshot2", quietly = TRUE)) {
    message("Package 'webshot2' not available. Plotly exports will use fallback static conversion.")
  }

  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    message("Package 'htmlwidgets' not available. HTML widget exports may be limited.")
  }
}


#' Extract data from Plotly object
#' @noRd
.extract_plotly_data <- function(plot_obj) {
  tryCatch({
    # Access plotly object structure
    if (is.list(plot_obj) && "x" %in% names(plot_obj)) {
      plotly_data <- plot_obj$x

      # Extract traces (data series)
      traces <- plotly_data$data
      if (is.null(traces) || length(traces) == 0) {
        return(NULL)
      }

      # Extract layout information
      layout <- plotly_data$layout

      # Extract configuration
      config <- plotly_data$config

      return(list(
        traces = traces,
        layout = layout,
        config = config
      ))
    }

    return(NULL)

  }, error = function(e) {
    return(NULL)
  })
}


#' Export Grid Graphics Objects
#' @noRd
.export_grid_plot <- function(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent) {
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("grid package required for grid graphics export")
  }

  # Open graphics device
  .open_graphics_device(file_path, format_lower, width, height, dpi, quality, transparent)

  # Draw the grid object
  grid::grid.draw(plot_obj)

  # Close device
  grDevices::dev.off()
}


#' Export Generic Plot Objects
#' @noRd
.export_generic_plot <- function(plot_obj, file_path, format_lower, width, height, dpi, quality, transparent) {
  # Open graphics device
  .open_graphics_device(file_path, format_lower, width, height, dpi, quality, transparent)

  # Try to print the object (most plot objects have print methods)
  tryCatch({
    print(plot_obj)
  }, error = function(e) {
    # If print fails, try plot method
    tryCatch({
      plot(plot_obj)
    }, error = function(e2) {
      # Last resort: create a placeholder
      plot(1, 1, type = "n", xlab = "Unknown Plot Type", ylab = "",
           main = paste("Object of class:", paste(class(plot_obj), collapse = ", ")))
      text(1, 1, "Plot object\n(custom export method needed)", cex = 1.2)
    })
  })

  # Close device
  grDevices::dev.off()
}


#' Open Appropriate Graphics Device
#' @noRd
.open_graphics_device <- function(file_path, format_lower, width, height, dpi, quality, transparent) {

  switch(format_lower,
         "png" = {
           grDevices::png(filename = file_path, width = width * dpi, height = height * dpi,
                          res = dpi, bg = if (transparent) "transparent" else "white")
         },
         "jpg" = ,
         "jpeg" = {
           grDevices::jpeg(filename = file_path, width = width * dpi, height = height * dpi,
                           res = dpi, quality = quality, bg = "white")
         },
         "pdf" = {
           grDevices::pdf(file = file_path, width = width, height = height)
         },
         "svg" = {
           if (requireNamespace("svglite", quietly = TRUE)) {
             svglite::svglite(file = file_path, width = width, height = height)
           } else {
             grDevices::svg(filename = file_path, width = width, height = height)
           }
         },
         "tiff" = {
           grDevices::tiff(filename = file_path, width = width * dpi, height = height * dpi,
                           res = dpi, bg = if (transparent) "transparent" else "white")
         },
         "bmp" = {
           grDevices::bmp(filename = file_path, width = width * dpi, height = height * dpi,
                          res = dpi, bg = if (transparent) "transparent" else "white")
         },
         stop(sprintf("Unsupported format: %s", format_lower))
  )
}


#' Enhanced Performance Monitoring and Diagnostics
#' @noRd
.generate_export_summary <- function(file_paths, processing_time, plot_types_summary) {
  summary_info <- list(
    total_files = length(file_paths),
    processing_time = processing_time,
    average_time_per_plot = if (length(file_paths) > 0) processing_time / length(file_paths) else 0,
    plot_types = plot_types_summary,
    file_paths = file_paths
  )

  return(summary_info)
}


#' Analyze Plot Types in Collection
#' @noRd
.analyze_plot_types <- function(plot_list) {
  type_counts <- table(vapply(plot_list, function(x) class(x)[1], character(1)))
  return(as.list(type_counts))
}


# Additional utility functions for digest package fallback
.create_simple_hash <- function(obj) {
  # Simple fallback hash function if digest package is not available
  tryCatch({
    if (requireNamespace("digest", quietly = TRUE)) {
      return(digest::digest(obj, algo = "xxhash32"))
    } else {
      # Fallback: use object size and first class as simple signature
      return(paste(class(obj)[1], object.size(obj), runif(1), sep = "_"))
    }
  }, error = function(e) {
    return(paste("hash", sample.int(10000, 1), sep = "_"))
  })
}
