#' Perform Linear Mixed-Effects Model Analysis
#'
#' @description
#' This function automates the process of fitting Linear Mixed-Effects Models (LMEMs)
#' to a dataset, typically for analyzing paired or grouped metabolic data. It iterates
#' through a list of features/metabolites, prepares the data for each, fits an LMEM,
#' and extracts key results, including fixed effects, random effects, and model diagnostics.
#' The function is designed to handle common data structures from clinical studies,
#' where multiple measurements (e.g., 'Start' and 'End') are taken from the same subject.
#' It also includes data validation, error handling, and p-value adjustment.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param featMet_names A character vector of column names representing the
#'   features or metabolites to be analyzed.
#' @param model_formula A character string representing the LMEM formula (e.g.,
#'   "Response ~ Start + End + (1|SubjectID)"). This is a required parameter.
#' @param group_cat1 A character vector of the two expected categories in
#'   \code{group_cat1_col} that represent the paired measurements (e.g.,
#'   \code{c("Start", "End")}).
#' @param group_cat2 A character vector of categories for an optional second
#'   grouping variable, such as different lung sections (e.g., \code{c("Left", "Right")}).
#'   If \code{NULL}, this grouping is ignored.
#' @param group_cat1_col A character string specifying the column name containing
#'   the paired categories (default: "Group").
#' @param group_cat2_col A character string specifying the column name for the
#'   optional second grouping variable (default: "Group2").
#' @param subject_col A character string specifying the column name for subject
#'   identifiers (default: "SubjectID").
#' @param response_col A character string specifying the column name for the
#'   response variable (e.g., "IschemiaTime"; default: "Response"). This column
#'   must be numeric.
#' @param adj_p A character string specifying the p-value adjustment method to
#'   be passed to \code{\link[stats]{p.adjust}}. The default is "BH" (Benjamini-Hochberg).
#'   Use "none" for no adjustment.
#' @param stars A logical value. If \code{TRUE}, adds significance stars to the
#'   results based on the adjusted p-value (default: \code{TRUE}).
#' @param exclude_qc A logical value. If \code{TRUE}, removes rows where the
#'   value in \code{group_cat1_col} is "QC" (default: \code{TRUE}).
#' @param verbose A logical value. If \code{FALSE} (default), does not output which feature is being analyzed.
#'
#' @return A list containing the results of the analysis. The list includes:
#'   \describe{
#'     \item{Individual results for each feature/metabolite:}{Each element is a list
#'       with the fitted model object, fixed and random effects, and diagnostics.}
#'     \item{Combined_Fixed_Effects_Summary:}{A data frame summarizing the fixed
#'       effects from all successful models, including adjusted p-values.}
#'     \item{Processing_Summary:}{A list with a summary of the analysis, including
#'       the number of successful and failed models.}
#'   }
#'
#' @importFrom dplyr select filter group_by summarise rename across all_of
#' @importFrom tidyr pivot_wider
#' @importFrom lme4 VarCorr
#' @importFrom lmerTest lmer
#' @importFrom tibble rownames_to_column
#' @importFrom stats p.adjust
#' @importFrom rlang sym
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # Example data setup
#' library(dplyr)
#' df_lm <- data.frame(
#'   SubjectID = rep(1:20, each = 2),
#'   Group = rep(c("Start", "End"), 20),
#'   IschemiaTime = runif(40, 10, 60),
#'   MetaboliteA = rnorm(40, mean = c(5, 7), sd = 1),
#'   MetaboliteB = rnorm(40, mean = c(12, 11), sd = 2)
#' )
#'
#' # Perform LMEM analysis for two metabolites
#' results <- perform_LMEM(
#'   df = df_lm,
#'   featMet_names = c("MetaboliteA", "MetaboliteB"),
#'   model_formula = "IschemiaTime ~ Start + End + (1|SubjectID)",
#'   group_cat1 = c("Start", "End"),
#'   response_col = "IschemiaTime"
#' )
#'
#' # View combined results table
#' print(results$Combined_Fixed_Effects_Summary)
#'
#' # View details for a single metabolite
#' print(summary(results$MetaboliteA$Model_Object))
#' }
#'
#' @export
#'
perform_LMEM <- function(
    df,
    featMet_names,
    model_formula = NULL,
    group_cat1 = NULL,
    group_cat2 = NULL,
    group_cat1_col = "Group",
    group_cat2_col = "Group2",
    subject_col = "SubjectID",
    response_col = "Response",
    adj_p = "BH",
    stars = TRUE,
    exclude_qc = TRUE,
    verbose = FALSE
) {

  # Load required packages silently
  suppressPackageStartupMessages({
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required")
    if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package 'tidyr' is required")
    if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required")
    if (!requireNamespace("lmerTest", quietly = TRUE)) stop("Package 'lmerTest' is required")
    if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required")
  })

  # Input validation
  if (is.null(model_formula)) {
    stop("model_formula must be provided")
  }

  if (is.null(group_cat1) || length(group_cat1) == 0) {
    stop("group_cat1 must be provided (e.g., c('Start', 'End'))")
  }

  # Ensure necessary columns exist for the base dataframe
  base_required_cols <- c(subject_col, group_cat1_col, response_col)
  if (!is.null(group_cat2) && length(group_cat2) > 0) {
    base_required_cols <- c(base_required_cols, group_cat2_col)
  }

  if (!all(base_required_cols %in% colnames(df))) {
    missing_base <- setdiff(base_required_cols, colnames(df))
    stop(paste("Input data frame must contain all of the following base columns:",
               paste(missing_base, collapse = ", ")))
  }

  # Ensure all feature/metabolite columns exist
  if (!all(featMet_names %in% colnames(df))) {
    missing_cols <- setdiff(featMet_names, colnames(df))
    stop(paste("Input data frame is missing the following feature columns:",
               paste(missing_cols, collapse = ", ")))
  }

  # Exclude QC if requested
  if (exclude_qc) {
    original_rows <- nrow(df)
    df <- df %>% dplyr::filter(!!sym(group_cat1_col) != "QC")
    excluded_rows <- original_rows - nrow(df)
    if (excluded_rows > 0) {
      message(paste("Excluded", excluded_rows, "QC rows from analysis"))
    }
  }

  # Check if Start/End categories exist in the Group column
  available_groups <- unique(df[[group_cat1_col]])
  missing_start_end <- setdiff(group_cat1, available_groups)
  if (length(missing_start_end) > 0) {
    stop(paste("Missing required categories in", group_cat1_col, "column:",
               paste(missing_start_end, collapse = ", "),
               "\nAvailable categories:", paste(available_groups, collapse = ", ")))
  }

  # Check lung categories if provided
  use_lung_filter <- FALSE
  available_lung_cats <- NULL
  if (!is.null(group_cat2) && length(group_cat2) > 0 && group_cat2_col %in% colnames(df)) {
    available_group2 <- unique(df[[group_cat2_col]])
    available_lung_cats <- intersect(group_cat2, available_group2)
    if (length(available_lung_cats) == 0) {
      warning(paste("No specified lung categories found in", group_cat2_col, "column.",
                    "Expected one of:", paste(group_cat2, collapse = ", "),
                    "\nAvailable categories:", paste(available_group2, collapse = ", "),
                    "\nProceeding without lung filtering."))
      use_lung_filter <- FALSE
    } else {
      use_lung_filter <- TRUE
      message(paste("Found lung categories in", group_cat2_col, ":",
                    paste(available_lung_cats, collapse = ", ")))
    }
  }

  # Validate Response column is numeric
  if (!is.numeric(df[[response_col]])) {
    stop(paste("Response column", response_col, "must be numeric"))
  }

  # Validate adj_p method
  valid_adj_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!adj_p %in% valid_adj_methods) {
    stop(paste("adj_p method must be one of:", paste(valid_adj_methods, collapse = ", ")))
  }

  # Initialize results containers
  all_results <- list()
  successful_models <- 0
  failed_models <- 0
  processing_details <- list()

  # Convert model formula to formula object
  tryCatch({
    model_formula_parsed <- as.formula(model_formula)
  }, error = function(e) {
    stop(paste("Invalid model formula:", e$message))
  })

  for (met_name in featMet_names) {
    if (verbose) message(paste("Processing metabolite:", met_name))

    # Initialize variables for this iteration
    df_transformed <- NULL
    current_model <- NULL

    # Flag to control loop iteration
    skip_iteration <- FALSE

    # 1. Create and prepare data for the current metabolite
    tryCatch({
      # Filter for lung categories if they exist and are specified
      df_filtered <- if (use_lung_filter) {
        df %>% dplyr::filter(!!sym(group_cat2_col) %in% available_lung_cats)
      } else {
        df
      }

      # Check if we have data after filtering
      if (nrow(df_filtered) == 0) {
        warning(paste("No data remaining after filtering for metabolite:", met_name))
        all_results[[met_name]] <- list(
          Status = "Skipped: No data after filtering",
          Error_Details = "No rows remaining after lung category filtering"
        )
        failed_models <- failed_models + 1
        skip_iteration <- TRUE
      }

      if (!skip_iteration) {
        # Create column mapping based on what's available
        select_cols <- list(
          SubjectID = sym(subject_col),
          Group = sym(group_cat1_col),
          Response = sym(response_col),
          MetaboliteValue = sym(met_name)
        )

        # Add lung column if using lung filter
        if (use_lung_filter) {
          select_cols$Lung <- sym(group_cat2_col)
        }

        # Select and prepare data
        df_prep <- df_filtered %>%
          dplyr::select(!!!select_cols) %>%
          dplyr::filter(Group %in% group_cat1) %>%  # Filter for specified groups (Start/End)
          dplyr::filter(!is.na(MetaboliteValue), !is.na(Response)) %>%
          dplyr::filter(if (use_lung_filter) Lung %in% available_lung_cats else TRUE)

        # Check if we have both required groups after filtering
        available_groups_filtered <- unique(df_prep$Group)
        if (!all(group_cat1 %in% available_groups_filtered)) {
          missing_groups <- setdiff(group_cat1, available_groups_filtered)
          warning(paste("Missing groups for metabolite", met_name, ":",
                        paste(missing_groups, collapse = ", ")))
          all_results[[met_name]] <- list(
            Status = "Skipped: Missing required group",
            Error_Details = paste("Missing groups:", paste(missing_groups, collapse = ", "))
          )
          failed_models <- failed_models + 1
          skip_iteration <- TRUE
        }
      }

      if (!skip_iteration) {
        # Check for subjects with measurements in required groups
        subject_group_counts <- df_prep %>%
          dplyr::group_by(SubjectID) %>%
          dplyr::summarise(
            groups_present = list(unique(Group)),
            n_groups = length(unique(Group)),
            .groups = "drop"
          ) %>%
          dplyr::filter(n_groups == length(group_cat1))  # Subjects with all required groups

        paired_subjects <- subject_group_counts$SubjectID

        if (length(paired_subjects) == 0) {
          warning(paste("No subjects with all required group measurements for metabolite:", met_name))
          all_results[[met_name]] <- list(
            Status = "Skipped: No paired measurements",
            Error_Details = paste("No subjects have measurements for all groups:", paste(group_cat1, collapse = ", "))
          )
          failed_models <- failed_models + 1
          skip_iteration <- TRUE
        }
      }

      if (!skip_iteration) {
        # Filter to only paired subjects
        df_prep <- df_prep %>%
          dplyr::filter(SubjectID %in% paired_subjects)

        # For subjects with multiple measurements per group, take the mean
        group_cols <- c("SubjectID", "Group")
        if (use_lung_filter) group_cols <- c(group_cols, "Lung")
        group_cols <- c(group_cols, "Response")

        df_prep <- df_prep %>%
          dplyr::group_by(across(all_of(group_cols))) %>%
          dplyr::summarise(MetaboliteValue = mean(MetaboliteValue, na.rm = TRUE), .groups = "drop")

        # Transform to wide format for Start/End
        id_cols <- c("SubjectID", "Response")
        if (use_lung_filter) id_cols <- c(id_cols, "Lung")

        df_transformed <- df_prep %>%
          tidyr::pivot_wider(
            id_cols = all_of(id_cols),
            names_from = Group,
            values_from = MetaboliteValue
          )

        # Check that we have the required columns after pivoting
        if (!all(group_cat1 %in% colnames(df_transformed))) {
          missing_cols_after_pivot <- setdiff(group_cat1, colnames(df_transformed))
          warning(paste("Missing columns after pivot for metabolite", met_name, ":",
                        paste(missing_cols_after_pivot, collapse = ", ")))
          all_results[[met_name]] <- list(
            Status = "Skipped: Missing columns after pivot",
            Error_Details = paste("Missing after pivot:", paste(missing_cols_after_pivot, collapse = ", "))
          )
          failed_models <- failed_models + 1
          skip_iteration <- TRUE
        }
      }

      if (!skip_iteration) {
        # Remove rows with missing values in required columns
        required_cols <- c(group_cat1, "Response")
        df_transformed <- df_transformed %>%
          dplyr::filter(if_all(all_of(required_cols), ~ !is.na(.)))

        if (nrow(df_transformed) == 0) {
          warning(paste("No complete cases for metabolite:", met_name))
          all_results[[met_name]] <- list(
            Status = "Skipped: No complete cases",
            Error_Details = "No rows with all required values"
          )
          failed_models <- failed_models + 1
          skip_iteration <- TRUE
        }
      }

      if (!skip_iteration) {
        processing_details[[met_name]]$paired_subjects <- length(paired_subjects)
        processing_details[[met_name]]$complete_observations <- nrow(df_transformed)

        # Rename Response column to match what's expected in the model formula
        if ("Response" %in% colnames(df_transformed) && response_col != "Response") {
          # Check if the model formula uses the original response column name
          formula_vars <- all.vars(model_formula_parsed)
          if (response_col %in% formula_vars) {
            df_transformed <- df_transformed %>%
              dplyr::rename(!!response_col := Response)
          } else if ("IschemiaTime" %in% formula_vars) {
            df_transformed <- df_transformed %>%
              dplyr::rename(IschemiaTime = Response)
          }
        }
      }
    }, error = function(e) {
      warning(paste("Error in data preparation for", met_name, ":", e$message))
      all_results[[met_name]] <- list(
        Status = "Skipped: Data preparation failed",
        Error_Details = e$message
      )
      failed_models <<- failed_models + 1
      skip_iteration <<- TRUE
    })

    # Skip to next metabolite if data preparation failed
    if (skip_iteration || is.null(df_transformed) || nrow(df_transformed) == 0) {
      next
    }

    # 2. Fit Linear Mixed-Effects Models
    tryCatch({
      # Check if we have enough subjects for mixed effects
      n_subjects <- length(unique(df_transformed$SubjectID))
      if (n_subjects < 3) {
        warning(paste("Too few subjects (", n_subjects, ") for mixed effects model for metabolite:", met_name))
        all_results[[met_name]] <- list(
          Status = "Skipped: Too few subjects",
          Error_Details = paste("Only", n_subjects, "subjects available, need at least 3")
        )
        failed_models <- failed_models + 1
        skip_iteration <- TRUE
      }

      if (!skip_iteration) {
        # Ensure factors are properly coded
        if ("SubjectID" %in% colnames(df_transformed)) {
          df_transformed$SubjectID <- as.factor(df_transformed$SubjectID)
        }
        if (use_lung_filter && "Lung" %in% colnames(df_transformed)) {
          df_transformed$Lung <- as.factor(df_transformed$Lung)
        }

        # Fit the model using lmerTest for p-values
        if (verbose) {
          current_model <- lmerTest::lmer(model_formula_parsed, data = df_transformed, REML = TRUE)
        } else {
          current_model <- suppressMessages(suppressWarnings(
            lmerTest::lmer(model_formula_parsed, data = df_transformed, REML = TRUE)
          ))
        }

        # Check for convergence warnings
        if (length(current_model@optinfo$conv$lme4$messages) > 0) {
          warning(paste("Convergence warnings for", met_name, ":",
                        paste(current_model@optinfo$conv$lme4$messages, collapse = "; ")))
        }

        current_summary <- summary(current_model)
        successful_models <- successful_models + 1
      }
    }, error = function(e) {
      warning(paste("Error fitting lmer model for", met_name, ":", e$message))
      all_results[[met_name]] <- list(
        Status = "Skipped: Model fitting failed",
        Error_Details = e$message,
        Data_Summary = list(
          n_subjects = if (exists("df_transformed")) length(unique(df_transformed$SubjectID)) else NA,
          n_observations = if (exists("df_transformed")) nrow(df_transformed) else NA
        )
      )
      failed_models <<- failed_models + 1
      skip_iteration <<- TRUE
    })

    # Skip to next metabolite if model fitting failed
    if (skip_iteration || is.null(current_model)) {
      next
    }

    # 3. Extract results safely
    tryCatch({
      met_results <- list()
      met_results$Status <- "Success"
      met_results$Model_Object <- current_model
      met_results$Data_Summary <- list(
        n_subjects = length(unique(df_transformed$SubjectID)),
        n_observations = nrow(df_transformed),
        lung_types = if (use_lung_filter && "Lung" %in% colnames(df_transformed)) table(df_transformed$Lung) else "Not filtered"
      )

      # REML criterion
      met_results$REML <- if (!is.null(current_summary$devcomp$cmp["REML"])) {
        current_summary$devcomp$cmp["REML"]
      } else {
        deviance(current_model)
      }

      # Scaled Residuals
      model_residuals <- resid(current_model)
      met_results$Scaled_Residuals <- data.frame(
        Min = min(model_residuals, na.rm = TRUE),
        Q25 = quantile(model_residuals, prob = 0.25, na.rm = TRUE),
        Median = median(model_residuals, na.rm = TRUE),
        Q75 = quantile(model_residuals, prob = 0.75, na.rm = TRUE),
        Max = max(model_residuals, na.rm = TRUE)
      )

      # Random Effects
      met_results$Random_Effects <- as.data.frame(lme4::VarCorr(current_model))

      # Fixed Effects
      fixed_effects_df <- as.data.frame(current_summary$coefficients)

      # Standardize column names
      expected_cols <- c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")
      actual_cols <- colnames(fixed_effects_df)

      # Create mapping for common variations
      col_mapping <- list(
        "Estimate" = c("Estimate"),
        "Std..Error" = c("Std. Error", "Std.Error"),
        "df" = c("df"),
        "t.value" = c("t value", "t.value"),
        "Pr...t.." = c("Pr(>|t|)", "Pr(>|t|)", "p.value")
      )

      # Rename columns systematically
      for (expected in names(col_mapping)) {
        for (variation in col_mapping[[expected]]) {
          if (variation %in% actual_cols) {
            colnames(fixed_effects_df)[colnames(fixed_effects_df) == variation] <- expected
            break
          }
        }
      }

      # Ensure we have the p-value column
      if ("Pr...t.." %in% colnames(fixed_effects_df)) {
        fixed_effects_df$p_value <- fixed_effects_df$`Pr...t..`
        fixed_effects_df$`Pr...t..` <- NULL # Remove after renaming to avoid duplication
      } else if ("Pr(>|t|)" %in% colnames(fixed_effects_df)) {
        fixed_effects_df$p_value <- fixed_effects_df$`Pr(>|t|)`
        fixed_effects_df$`Pr(>|t|)` <- NULL # Remove after renaming to avoid duplication
      }

      fixed_effects_df <- fixed_effects_df %>%
        tibble::rownames_to_column(var = "Effect")

      # Confidence Intervals for Fixed Effects
      conf_int_df <- tryCatch({
        ci_result <- suppressMessages(confint(current_model, level = 0.95, method = "Wald"))
        if (is.matrix(ci_result)) {
          data.frame(ci_result, check.names = FALSE) %>%
            tibble::rownames_to_column(var = "Effect") %>%
            dplyr::rename_with(~ c("Effect", "CI95_Lower", "CI95_Upper"), everything())
        } else {
          data.frame(Effect = fixed_effects_df$Effect, CI95_Lower = NA, CI95_Upper = NA)
        }
      }, error = function(e) {
        warning(paste("Could not compute confidence intervals for", met_name, ":", e$message))
        data.frame(Effect = fixed_effects_df$Effect, CI95_Lower = NA, CI95_Upper = NA)
      })

      met_results$Fixed_Effects <- fixed_effects_df %>%
        dplyr::left_join(conf_int_df, by = "Effect") %>%
        tibble::column_to_rownames(var = "Effect")

      # Correlation of Fixed Effects (safely)
      tryCatch({
        vcov_summary <- vcov(current_model)
        if (!is.null(vcov_summary) && nrow(vcov_summary) > 1) {
          corr_matrix <- cov2cor(vcov_summary)
          met_results$Corr_Fixed_Effects <- as.data.frame(as.matrix(corr_matrix))
        } else {
          met_results$Corr_Fixed_Effects <- "Single parameter or not available"
        }
      }, error = function(e) {
        met_results$Corr_Fixed_Effects <- "Not available due to error"
        warning(paste("Could not extract correlation matrix for", met_name, ":", e$message))
      })

      all_results[[met_name]] <- met_results

    }, error = function(e) {
      warning(paste("Error extracting results for", met_name, ":", e$message))
      all_results[[met_name]] <- list(
        Status = "Partial failure: Results extraction failed",
        Error_Details = e$message,
        Model_Object = if (!is.null(current_model)) current_model else NULL
      )
      failed_models <- failed_models + 1
    })
  }

  # Create combined fixed effects summary
  tryCatch({
    successful_results <- all_results[sapply(all_results, function(x) {
      is.list(x) && "Fixed_Effects" %in% names(x) && "Status" %in% names(x) && x$Status == "Success"
    })]

    if (length(successful_results) > 0) {
      combined_fixed_effects <- do.call(rbind, lapply(names(successful_results), function(met_name) {
        fe_data <- successful_results[[met_name]]$Fixed_Effects
        if (!is.null(fe_data) && nrow(fe_data) > 0) {
          data.frame(
            Metabolite_Name = met_name,
            Effect = rownames(fe_data),
            fe_data,
            row.names = NULL,
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      }))

      # Add adjusted p-values and stars
      if (!is.null(combined_fixed_effects) && "p_value" %in% colnames(combined_fixed_effects)) {

        # Initialize the adj_p column with the original p-values
        combined_fixed_effects$adj_p <- combined_fixed_effects$p_value

        # Adjust p-values for each unique factor in the 'Effect' column
        if (!is.na(adj_p) && adj_p != "none") {
          unique_effects <- unique(combined_fixed_effects$Effect)
          for (effect in unique_effects) {
            indices <- combined_fixed_effects$Effect == effect
            combined_fixed_effects$adj_p[indices] <- p.adjust(combined_fixed_effects$p_value[indices], method = adj_p)
          }
        }

        if (stars) {
          combined_fixed_effects$Stars <- ifelse(is.na(combined_fixed_effects$adj_p), "",
                                                 ifelse(combined_fixed_effects$adj_p < 0.001, "***",
                                                        ifelse(combined_fixed_effects$adj_p < 0.01, "**",
                                                               ifelse(combined_fixed_effects$adj_p < 0.05, "*", ""))))
        }
      }

      all_results$Combined_Fixed_Effects_Summary <- combined_fixed_effects
    } else {
      all_results$Combined_Fixed_Effects_Summary <- data.frame()
      warning("No successful models to combine")
    }
  }, error = function(e) {
    warning(paste("Error creating combined summary:", e$message))
    all_results$Combined_Fixed_Effects_Summary <- data.frame()
  })

  # Add processing summary
  all_results$Processing_Summary <- list(
    Total_FeaturesMetabolites = length(featMet_names),
    Successful_Models = successful_models,
    Failed_Models = failed_models,
    Success_Rate = round(successful_models / length(featMet_names) * 100, 2),
    Processing_Details = processing_details
  )

  message(paste("\nProcessing complete!"))
  message(paste("Successful models:", successful_models, "out of", length(featMet_names)))
  message(paste("Success rate:", round(successful_models / length(featMet_names) * 100, 2), "%"))

  class(all_results) <- c("perform_LMEM", "list")
  return(all_results)
}

# -----------------------------------------------------------------------------

#' Plot Residuals and QQ-Plot for an LMER Model
#'
#' This function generates two diagnostic plots for a fitted linear mixed-effects
#' model object from the `lmerTest` or `lme4` package: a **residuals vs. fitted**
#' plot and a **Normal Q-Q plot**. These plots are essential for checking the
#' assumptions of linearity, homoscedasticity, and normality of residuals.
#'
#' @param lmer_model An object of class \code{lmerMod} or \code{lmerModLmerTest}
#'   from a fitted model.
#'
#' @return A list containing the two plot objects:
#'   \describe{
#'     \item{residuals_plot:}{A \code{ggplot} object for the residuals vs. fitted plot.}
#'     \item{qq_plot:}{A \code{ggplot} object for the Normal Q-Q plot.}
#'   }
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_hline labs
#'   ggtitle theme_minimal theme element_text
#' @importFrom stats fitted resid qt qqnorm qqline
#' @examples
#' \dontrun{
#' # Assuming `lmer_model` is a fitted model object
#' # For example:
#' # library(lmerTest)
#' # model <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
#' # plots <- plot_lmer_diagnostics(model)
#' # print(plots$residuals_plot)
#' # print(plots$qq_plot)
#' }
plot_lmer_diagnostics <- function(lmer_model) {
  # Check for correct input class
  if (!inherits(lmer_model, c("lmerMod", "lmerModLmerTest"))) {
    stop("Input must be an object of class 'lmerMod' or 'lmerModLmerTest'")
  }

  # Extract residuals and fitted values
  model_data <- data.frame(
    fitted = fitted(lmer_model),
    residuals = resid(lmer_model)
  )

  # 1. Residuals vs. Fitted plot
  residuals_plot <- ggplot(model_data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    labs(
      x = "Fitted Values",
      y = "Residuals"
    ) +
    ggtitle("Residuals vs. Fitted Values") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # 2. Normal Q-Q plot
  # Get standardized residuals (for better QQ plot interpretation)
  std_res <- resid(lmer_model, scaled = TRUE)

  # Create a data frame for the Q-Q plot
  qq_data <- data.frame(
    theoretical = qt(
      ppoints(length(std_res)),
      df = df.residual(lmer_model)
    ),
    sample = sort(std_res)
  )

  qq_plot <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_point() +
    labs(
      x = "Theoretical Quantiles",
      y = "Standardized Residuals"
    ) +
    ggtitle("Normal Q-Q Plot of Standardized Residuals") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  return(list(
    residuals_plot = residuals_plot,
    qq_plot = qq_plot
  ))
}

# S3 Methods
#' @export
print.perform_LMEM <- function(x, ...) {
  cat("=== Linear Mixed-Effects Models ===\n")
  cat("Total Features:  ", x$Processing_Summary$Total_FeaturesMetabolites, "\n")
  cat("Successful:      ", x$Processing_Summary$Successful_Models, "\n")
  cat("Success Rate:    ", x$Processing_Summary$Success_Rate, "%\n")
  invisible(x)
}

#' @export
summary.perform_LMEM <- function(object, ...) {
  # Get top significant effects
  top_sig <- head(object$Combined_Fixed_Effects_Summary[
    order(object$Combined_Fixed_Effects_Summary$adj_p),
    c("Metabolite_Name", "Effect", "Estimate", "adj_p", "Stars")], 10)

  ans <- list(
    summary = object$Processing_Summary,
    top_effects = top_sig
  )
  class(ans) <- "summary.perform_LMEM"
  return(ans)
}

#' @export
print.summary.perform_LMEM <- function(x, ...) {
  cat("---------------------------------------\n")
  cat("LMEM Analysis Summary\n")
  cat("---------------------------------------\n")
  cat("Models Converged:", x$summary$Successful_Models, "\n")
  cat("Models Failed:   ", x$summary$Failed_Models, "\n")

  cat("\n-- Top Significant Fixed Effects (Adjusted P) --\n")
  if (nrow(x$top_effects) > 0) {
    print(x$top_effects, row.names=FALSE)
  } else {
    cat("No significant effects found.\n")
  }
  invisible(x)
}
