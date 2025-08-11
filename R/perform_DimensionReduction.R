#' Function to Perform Dimension Reduction Techniques on the Processed Data
#'
#' @description
#' This function performs dimension reduction techniques on the preprocessed data. These
#' techniques include Principal Component Analysis (PCA). and
#' Orthogonal Partial Least Squares-Discriminant Analysis (OPLS-DA). Techniques to be added
#' are Partial Least Squares-Discriminant Analysis (PLS-DA) and
#' sparse Partial Least Squares-Discriminant Analysis (sPLS-DA). Parameters will be ignored
#' if they are not usable in the 'type' indicated.
#'
#' @param data List. This list must be a result from the `performPreprocessingPeakData` function.
#' @param type String. The type of dimension reduction technique to perform. Choices are below. Future choices will include "PLS-DA" and "sPLS-DA".
#'   \itemize{
#'     \item "PCA": Perform Principal Component Analysis (PCA).
#'     \item "OPLS-DA": Perform Orthogonal Partial Least Squares-Discriminant Analysis (OPLS-DA).
#'     }
#'     Defaults to "PCA".
#' @param reduceWhat String. Controls what you want to perform the selected `type` on.
#'   \itemize{
#'     \item "QC": Perform the selected `type` on QC = Quality Control Samples.
#'     \item "BS": Perform the selected `type` on BS = Biological Samples.
#'     }
#'     Defaults to "QC".
#' @param arrangeLevels Vector. Determines how the groups will be arranged. The format could be "c('group1', 'group2', ...)". Defaults to `NULL` which sorts the groups in alphabetical order.
#' @param screeWhat String. Determines what will be plotted in the scree plot (also called elbow plot).
#'   \itemize{
#'     \item "variance": Plot the variance explained.
#'     \item "eigenvalue": Plot the eigen values.
#'     }
#'     Defaults to "variance".
#' @param screeType String. Controls how the scree plot is plotted.
#'   \itemize{
#'     \item "bar": Uses bar plots.
#'     \item "line": Uses line plots.
#'     \item "both": Uses both bar and line plots.
#'     }
#'     Defaults to "both".
#' @param screeLabs Boolean. If `TRUE` (default), adds labels.
#' @param screeTitle String. The scree plot title. If `NULL` or empty (""), defaults to "Scree plot".
#' @param screeXLab String. The scree plot x-axis label. Defaults to "Principal Components".
#' @param screeYLab String. The scree plot y-axis label. Defaults to "% Variance Explained".
#' @param scoresEllipse Boolean. If `TRUE` (default), adds an ellipse in the scores plot.
#' @param scoresTitle String. The scores plot title.
#' @param scoresLegend String. The title in the legend section of the scores plot. Defaults to `NULL` which means no legend title.
#'
#' @returns Returns a range of list of results from plots, data frames, etc.
#' @export
#'
#' @examples
#' \dontrun{
#' performDimensionReduction(
#'   data = data_from_performPreprocessingPeakData_function
#' )
#' }

performDimensionReduction <- function(
    data,
    type          = "PCA",
    reduceWhat    = "QC",
    arrangeLevels = NULL,
    screeWhat     = "variance",
    screeType     = "both",
    screeLabs     = TRUE,
    screeTitle    = "",
    screeXLab     = "Principal Components",
    screeYLab     = "% Variance Explained",
    scoresEllipse = TRUE,
    scoresTitle   = NULL,
    scoresLegend  = NULL

) {

  # Parameter Validation Function
  check_parameters <- function(
    arrangeLevels,
    screeWhat,
    screeType,
    screeLabs,
    screeTitle,
    screeXLab,
    screeYLab,
    scoresEllipse,
    scoresTitle,
    scoresLegend
  ) {
    errors <- character()

    # Check choices
    allowed_scales <- c("QC", "BS")
    if (!(reduceWhat %in% allowed_scales)) {
      errors <- c(errors, "reduceWhat: Must be either QC or BS. QC = Quality Control Samples; BS = Biological Samples.")
    }

    allowed_scales <- c("variance", "eigenvalue")
    if (!(screeWhat %in% allowed_scales)) {
      errors <- c(errors, "screeWhat: Must be either variance or eigenvalue.")
    }

    allowed_scales <- c("both", "bar", "line")
    if (!(screeType %in% allowed_scales)) {
      errors <- c(errors, "screeType: Must be either both, bar, or line.")
    }

    # Check boolean
    if (!is.logical(screeLabs)) {
      errors <- c(errors, "screeLabs: Not logical. Must be either TRUE or FALSE.")
    }

    if (!is.logical(scoresEllipse)) {
      errors <- c(errors, "scoresEllipse: Not logical. Must be either TRUE or FALSE.")
    }

    # Check custom names
    if (!is.null(scalePCA)) {
      if (!(scalePCA %in% allowed_scales)) {
        errors <- c(errors, "scalePCA: Must be either NULL, mean, meanSD, or meanSD2.")
      }
    }

  }

  dimensionReductionResults       <- list() # Empty list to store all results
  dimensionReductionResults$Class <- "performDimensionReduction"  # Update list

  qc_indices     <- data$Metadata$Groups == "QC"
  non_qc_indices <- !qc_indices

  # If PCA for QC samples
  if (type == "PCA" && reduceWhat == "QC") {

    if (is.null(data$data_scaledPCA)) {
      data_pca <- data$data_preprocessed[qc_indices, ]

    } else if (length(data$data_scaledPCA) == 0) {
      data_pca <- data$data_preprocessed[qc_indices, ]

    } else {
      data_pca <- data$data_scaledPCA[qc_indices, ]
    }

    dimensionReductionResults$PCAData <- data_pca # Update list

    # Perform PCA
    pca_res <- stats::prcomp(data_pca, center = FALSE, scale. = FALSE)

    dimensionReductionResults$All_PCAResults    <- pca_res # Update list
    dimensionReductionResults$SD_of_PrinComps   <- pca_res$sdev %>% as.data.frame() %>% `colnames<-`(NULL) # Update list
    dimensionReductionResults$Variable_Loadings <- pca_res$rotation %>% as.data.frame() # Update list

    if (screeType == "both") { # Set value of both
      screeType <- c("bar", "line")
    }

    # Scree plot
    scree_plot <- factoextra::fviz_eig(pca_res, choice = screeWhat, geom = screeType,
                                       barfill = "steelblue", barcolor = "steelblue", linecolor = "black", ncp = 10,
                                       addlabels = screeLabs,
                                       hjust = 0,
                                       main = screeTitle, xlab = screeXLab, ylab = screeYLab,
                                       ggtheme = theme_minimal())

    dimensionReductionResults$ScreePlot <- scree_plot # Update list

    # Scores plot using ggplot2
    # Prepare data
    df_scores_qc <- data.frame(PC1 = pca_res$x[, 1],
                               PC2 = pca_res$x[, 2],
                               Batch = data$Metadata$Batch[qc_indices],
                               Label = rownames(data_pca) # Add labels
    )
    variance_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
    x_label            <- paste0("PC1 (", variance_explained[1], "%)")
    y_label            <- paste0("PC2 (", variance_explained[2], "%)")

    # Plot
    scores_plot <- ggplot(df_scores_qc,
                          aes(x = PC1, y = PC2, color = factor(Batch), shape = factor(Batch))) +
      geom_point() +
      # Add labels to each data points
      geom_text(aes(label = Label), vjust = -0.5, hjust = 0.5, size = 3) +
      # Add circle at 95% confidence interval at 'level'
      stat_ellipse(aes(group = factor(Batch), fill = factor(Batch)),
                   # orig = .95 try 0 if results to no circles
                   level = 0.95, geom = "polygon", alpha = 0.5) +
      # Add x-, y-labels, title, etc.
      labs( x = x_label, y = y_label, title = scoresTitle, color = "Batch", shape = "Batch", fill = "Batch") +
      theme_minimal() +
      theme(legend.position = "bottom")

    dimensionReductionResults$ScoresPlot <- scores_plot # Update list

    # If PCA for Biological Samples
  } else if (type == "PCA" && reduceWhat == "BS") {

    if (is.null(data$data_scaledPCA)) {
      data_pca <- data$data_preprocessed[non_qc_indices, ]
    } else if (length(data$data_scaledPCA) == 0) {
      data_pca <- data$data_preprocessed[non_qc_indices, ]
    } else {
      data_pca <- data$data_scaledPCA[non_qc_indices, ]
    }

    dimensionReductionResults$PCAData <- data_pca # Update list

    # Perform PCA
    pca_res <- stats::prcomp(data_pca, center = FALSE, scale. = FALSE)

    dimensionReductionResults$All_PCAResults <- pca_res # Update list
    dimensionReductionResults$SD_of_PrinComps <- pca_res$sdev %>% as.data.frame() %>% `colnames<-`(NULL)# Update list
    dimensionReductionResults$Variable_Loadings <- pca_res$rotation %>% as.data.frame() # Update list

    if (screeType == "both") { # Set value of both
      screeType <- c("bar", "line")
    }

    # Scree plot
    scree_plot <- factoextra::fviz_eig(pca_res, choice = screeWhat, geom = screeType,
                                       barfill = "steelblue", barcolor = "steelblue", linecolor = "black", ncp = 10,
                                       addlabels = screeLabs,
                                       hjust = 0,
                                       main = screeTitle, xlab = screeXLab, ylab = screeYLab,
                                       ggtheme = theme_minimal())

    dimensionReductionResults$ScreePlot <- scree_plot # Update list

    # Scores plot using ggplot2
    # Prepare data

    df_scores_qc <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                               Group = data$Metadata$Groups[non_qc_indices],
                               # Group = if (!is.null(arrangeLevels)) {
                               #   factor(data$Metadata$Groups[non_qc_indices], levels = arrangeLevels)
                               # } else {
                               #   as.factor(data$Metadata$Groups[non_qc_indices])
                               # },
                               Batch = data$Metadata$Batches[non_qc_indices],
                               Label = rownames(data_pca) # Add labels
    )
    variance_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
    x_label            <- paste0("PC1 (", variance_explained[1], "%)")
    y_label            <- paste0("PC2 (", variance_explained[2], "%)")

    # Plot
    scores_plot <- ggplot(df_scores_qc,
                          aes(x = PC1, y = PC2, color = Group, shape = Group)) +
      geom_point() +
      # Add labels to each data points
      geom_text(aes(label = Label), vjust = -0.5, hjust = 0.5, size = 3) +
      # Add circle at 95% confidence interval at 'level'
      stat_ellipse(aes(group = Group, fill = Group),
                   # orig = .95 try 0 if results to no circles
                   level = 0.95, geom = "polygon", alpha = 0.5) +
      # Add x-, y-labels, title, etc.
      labs(x = x_label, y = y_label, title = scoresTitle, color = scoresLegend, shape = scoresLegend, fill = scoresLegend) +
      theme_minimal() +
      theme(legend.position = "bottom")

    dimensionReductionResults$ScoresPlot <- scores_plot # Update list

    # If OPLS-DA
  } else if (type == "OPLS-DA") {

    if (is.null(data$data_scaledOPLSDA)) {
      data_oplsda <- data$data_preprocessed
    } else if (length(data$data_scaledOPLSDA) == 0) {
      data_oplsda <- data$data_preprocessed
    } else {
      data_oplsda <- data$data_scaledOPLSDA
    }

    dimensionReductionResults$OPLSDAData <- data_oplsda[non_qc_indices, ] %>% as.data.frame() # Update list

    # # Perform OPLS-DA
    # oplsda_results <- ropls::opls(
    #   x         = data_oplsda[non_qc_indices, ],
    #   y         = data$Metadata$Groups[non_qc_indices],
    #   predI     = 1, # 1 for OPLS; NA for PLS
    #   orthoI    = NA, # NA for OPLS; 0 for PLS
    #   # nipals for OPLS; svd for PLS
    #   algoC     = "nipals", # c("default", "nipals", "svd")[1] - this is to replace NAs if any, but this is already done in preprocessing
    #   crossvalI = 10, # 10-fold cross validation
    #   log10L    = FALSE, # no need, already did transformation
    #   permI     = 20, # defaults to 20 if no train/test, 0 if otherwise
    #   # Already did data scaling
    #   scaleC    = "none", # c("none", "center", "pareto", "standard")[4],
    #   subset    = NULL, # integer vector of training data, NULL for no partition, 'odd' for equal size of train and test dataset
    #   plotSubC  = NA
    # )



    'NEW ANALYSIS'
    'NEW CODE'


    # Parameter-check for arrangeLevels
    if (is.null(arrangeLevels)) {

      # Get all unique group combinations for pairwise comparison
      group_combinations <- utils::combn(
        x        = unique(data$Metadata$Groups[non_qc_indices]),
        # x        = arrangeLevels,
        m        = 2, # number of elements to choose
        simplify = FALSE # Output as list
      )

    } else if (!is.null(arrangeLevels)) {
      if (length(setdiff(arrangeLevels, data$Metadata$Groups[non_qc_indices])) > 0) {
        stop(paste0("Check 'arrangeLevels' values. There might be (1) typos (2) Group not in the data or (3) Group is missing. Here are the groups provided: ", paste(arrangeLevels, collapse = ", ")))

      } else {

        # Get all unique group combinations for pairwise comparison
        group_combinations <- utils::combn(
          # x        = unique(data$Metadata$Groups[non_qc_indices] %>% factor(levels = arrangeLevels)),
          x        = arrangeLevels,
          m        = 2, # number of elements to choose
          simplify = FALSE # Output as list
        )
      }
    }

    # Loop through each pairwise group combination
    for (pair in group_combinations) {

      group1 <- pair[1]
      group2 <- pair[2]

      # Subset data for the current pair
      current_groups   <- data$Metadata$Groups %in% c(group1, group2)
      data_oplsda_pair <- data_oplsda[current_groups, ]
      y_pair           <- data$Metadata$Groups[current_groups]

      # Save results uniquely for each pair
      comparison_label <- paste0(group1, "vs", group2, sep = "_")

      # Perform OPLS-DA
      oplsda_results <- ropls::opls(
        x         = data_oplsda_pair,
        y         = y_pair,
        predI     = 1,
        orthoI    = NA,
        algoC     = "nipals",
        crossvalI = 10,
        log10L    = FALSE,
        permI     = 20,
        scaleC    = "none",
        subset    = NULL,
        plotSubC  = paste0(group1, " vs ", group2)
      )

      # Store model results
      dimensionReductionResults[[paste0("results_OPLSDA_", group1, " vs. ", group2)]] <- oplsda_results

      # Proceed only if a valid model was built
      if (inherits(oplsda_results, "opls") && length(oplsda_results@summaryDF) > 0) {

        # Extract top 20 VIP scores

        all_vip <- ropls::getVipVn(oplsda_results, orthoL = FALSE) %>%
          enframe(name = "Feature", value = "VIP") %>%
          arrange(desc(VIP)) %>%
          data.frame()

        dimensionReductionResults[[paste0("data_VIPScores_", group1, " vs. ", group2)]] <- all_vip%>% as.data.frame()

        top20_vip <- ropls::getVipVn(oplsda_results, orthoL = FALSE) %>%
          enframe(name = "Feature", value = "VIP") %>%
          arrange(desc(VIP)) %>%
          slice_head(n = 20) %>%
          data.frame()

        # Prepare abundance data
        abundance_data <- data_oplsda_pair %>%
          as.data.frame() %>%
          select(all_of(top20_vip$Feature)) %>%
          mutate(Group = factor(y_pair)) %>%
          pivot_longer(-Group, names_to = "Feature", values_to = "Abundance") %>%
          mutate(Feature = factor(Feature, levels = rev(top20_vip$Feature)))

        dimensionReductionResults[[paste0("data_Abundance_", group1, " vs. ", group2)]] <- abundance_data%>% as.data.frame()

        # VIP Score Plot
        vip_plot <- ggplot(top20_vip, aes(x = VIP, y = reorder(Feature, VIP))) +
          geom_point(size = 4, color = "steelblue") +
          labs(title = paste0("Top 20 Predictive Features by VIP Score", " (", group1, " vs. ", group2, ")"),
               x = "VIP Score", y = NULL) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8))

        # Abundance Heatmap
        abundance_plot <- ggplot(abundance_data, aes(x = Group, y = Feature, fill = Abundance)) +
          geom_tile() +
          scale_fill_gradient(low = "darkgreen", high = "red",
                              breaks = range(abundance_data$Abundance),
                              labels = round(range(abundance_data$Abundance), 1)) +
          labs(title = NULL, x = NULL, y = NULL) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text.y = element_blank(),
                aspect.ratio = 15/3)

        # Combine VIP and Abundance plots
        vip_abundance_plot <- vip_plot + abundance_plot

        print(vip_abundance_plot)

        dimensionReductionResults[[paste0("plot_VIPAbundance_", group1, " vs. ", group2)]] <- vip_abundance_plot

        # S-Plot Data
        splot_data <- tibble(
          Variable = colnames(data_oplsda_pair),
          Covariance = map_dbl(as.data.frame(data_oplsda_pair), ~ cov(.x, oplsda_results@scoreMN[, 1])),
          Correlation = map_dbl(as.data.frame(data_oplsda_pair), ~ cor(.x, oplsda_results@scoreMN[, 1]))
        )

        dimensionReductionResults[[paste0("data_SPlot_", group1, " vs. ", group2)]] <- splot_data%>% as.data.frame()

        # Create S-plot
        s_plot <- ggplot(splot_data,
                         aes(x = Covariance, y = Correlation)) +
          geom_point(aes(color = Correlation), size = 3) +
          scale_color_gradient2(low = "blue",
                                mid = "grey",
                                high = "red",
                                midpoint = 0) +
          theme_minimal() +
          labs(
            title = paste0("S-Plot of OPLS-DA ", "(", group1, " vs. ", group2, ")"),
            x = "Covariance [p1]",
            y = "Correlation [p(corr)1]",
            color = "Correlation"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
          )

        dimensionReductionResults[[paste0("plot_SPlot_", group1, " vs. ", group2)]] <- s_plot

        # # Store plots in results list
        # dimensionReductionResults[[paste0(group1, "vs", group2, "plots")]] <- list(
        #   VIP_Abundance = vip_abundance_plot,
        #   S_Plot = s_plot
        # )
      }
    }

  }

  return(dimensionReductionResults)

}
