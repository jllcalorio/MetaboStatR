# THIS IS NOT YET WORKING

# Visualize before and after normalization

# Combined Plot Function
plotBeforeAfter <- function(
    data,
    scaled = "OPLS-DA", # c("PCA", "OPLS-DA")
    # title_prefix,
    group_by = "Sample" # c("Sample", "Feature")
) {

  library(patchwork)

  # Parameter checks
  if (data$Metadata$FunctionOrigin != "performPreprocessingPeakData") {
    stop("Data did not come from 'performProcessingPeakData' function.")
  }

  if (!is.null(scaled)) {
    if (scaled != "PCA" && scaled != "OPLS-DA") {
      stop("scaled: Must be either PCA or OPLS-DA")
    }
  }


  beforeAfterResults <- list() # Empty list to store results
  beforeAfterResults$Class <- "plotBeforeAfter" # Update list


  # Prepare Data for Visualization

  non_qc_indices <- data$Metadata$Groups != "QC" # Non-QC indices

  original <- data$data_no_NA[non_qc_indices, ] # Select this as original data

  if (scaled == "PCA") {
    transformed <- data$data_scaledPCA[non_qc_indices, ]
  } else if (scaled == "OPLS-DA") {
    transformed <- data$data_scaledOPLSDA[non_qc_indices, ]
  } else {
    stop("scaled: Must be either PCA or OPLS-DA")
  }

  beforeAfterResults$data_before <- original # Update list
  beforeAfterResults$data_after <- transformed # Update list



  # Convert data to long format for ggplot2
  long_format <- function(data) {
    data %>%
      as.data.frame() %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(-Sample, names_to = "Feature", values_to = "Value")
  }

  # Convert data to long format
  orig_long  <- long_format(original)
  trans_long <- long_format(transformed)

  # Density plots
  p1 <- ggplot(orig_long, aes(x = Value)) +
    geom_density(alpha = 0.5) +
    # ggtitle(paste0(title_prefix, " - Before Scaling")) +
    theme_minimal() +
    theme(legend.position = "none")

  p2 <- ggplot(trans_long, aes(x = Value)) +
    geom_density(alpha = 0.5) +
    # ggtitle(paste0(title_prefix, " - After Scaling")) +
    theme_minimal() +
    theme(legend.position = "none")

  # Boxplots
  if (group_by == "Sample") {
    p3 <- ggplot(orig_long, aes(x = Sample, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot Before Scaling")) +
      theme_minimal() +
      coord_flip()

    p4 <- ggplot(trans_long, aes(x = Sample, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot After Scaling")) +
      theme_minimal() +
      coord_flip()
  } else if (group_by == "Feature") {

    random_features <- sample(colnames(transformed), 30)  # Select random 30 features

    p3 <- ggplot(orig_long %>% filter(Feature %in% random_features),
                 aes(x = Feature, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot Before Scaling")) +
      theme_minimal() +
      coord_flip()

    p4 <- ggplot(trans_long %>% filter(Feature %in% random_features),
                 aes(x = Feature, y = Value)) +
      geom_boxplot() +
      # ggtitle(paste0(title_prefix, " - Boxplot After Scaling")) +
      theme_minimal() +
      coord_flip()
  } else {
    stop("The 'group_by' parameter must be 'Sample' or 'Feature'")
  }


  beforeAfterResults$plot_density_before <- p1 # Update list
  beforeAfterResults$plot_density_after  <- p2 # Update list
  beforeAfterResults$plot_box_before     <- p3 # Update list
  beforeAfterResults$plot_box_after      <- p4 # Update list


  # Combine plots
  all_plots <- (p1 | p2) / (p3 | p4)


  beforeAfterResults$plot_4in1           <- all_plots # Update list


  return(beforeAfterResults)
}

# Usage
# myplotbeforeafter <- plotBeforeAfter(
#   data = mydata2,
#   scaled = "OPLS-DA",
#   # title_prefix = "XDD",
#   group_by = "Sample"
# )

