# # Function for statistical analysis
performComparativeAnalysis <- function(
    data,
    sort_p         = TRUE,
    paired         = FALSE,
    plot_iden_met = NULL
) {

  comparativeAnalysisResults       <- list()
  comparativeAnalysisResults$Class <- "performComparativeAnalysis"

  non_qc_indices <- data$Metadata$Groups != "QC"
  df             <- data$data_scaledOPLSDA[non_qc_indices, ]
  groups         <- data$Metadata$Groups[non_qc_indices]

  num_groups <- length(unique(groups))

  # Two-group comparison
  if (num_groups == 2) {
    unique_groups     <- unique(groups)
    idx_group1        <- which(groups == unique_groups[1])
    idx_group2        <- which(groups == unique_groups[2])

    results           <- apply(df, 2, function(x) {
      residuals       <- residuals(lm(x ~ groups))
      # normality       <- shapiro.test(x)$p.value > 0.05
      normality       <- shapiro.test(residuals)$p.value > 0.05
      variance        <- var.test(x[idx_group1], x[idx_group2])$p.value > 0.05

      if (paired) {
        if (normality) {
          test_used   <- "Paired t-test"
          p_value     <- t.test(x[idx_group1], x[idx_group2], paired = TRUE)$p.value
        } else {
          test_used   <- "Wilcoxon Signed-Rank test"
          p_value     <- wilcox.test(x[idx_group1], x[idx_group2], paired = TRUE)$p.value
        }
      } else {
        if (normality) {
          if (variance) {
            test_used <- "Independent Samples t-test"
            p_value   <- t.test(x[idx_group1], x[idx_group2])$p.value
          } else {
            test_used <- "Welch's t-test"
            p_value   <- t.test(x[idx_group1], x[idx_group2], var.equal = FALSE)$p.value
          }
        } else {
          test_used   <- "Mann-Whitney U test"
          p_value     <- wilcox.test(x[idx_group1], x[idx_group2])$p.value
        }
      }
      return(c(test_used, p_value))
    })

  } else { # Multi-group comparison

    results <- apply(df, 2, function(x) {
      normality     <- all(by(x, groups, function(g) shapiro.test(g)$p.value > 0.05))
      variance      <- bartlett.test(x ~ groups)$p.value > 0.05

      if (paired) {
        test_used   <- "Repeated Measures ANOVA (MANOVA)"
        p_value     <- summary(aov(x ~ groups + Error(groups)))$p.value
      } else {
        if (normality && variance) {
          test_used <- "ANOVA"
          p_value   <- summary(aov(x ~ groups))[[1]][["Pr(>F)"]][1]
        } else {
          test_used <- "Kruskal-Wallis"
          p_value   <- kruskal.test(x ~ groups)$p.value
        }
      }
      return(c(test_used, p_value))
    })
  }

  results                <- as.data.frame(t(results), stringsAsFactors = FALSE)
  colnames(results)      <- c("Test used", "p-value")
  results$`p-value`      <- as.numeric(results$`p-value`)
  results$`adj. p-value` <- p.adjust(results$`p-value`, method = "BH")

  if(sort_p) {
    results <- results %>% dplyr::arrange(`adj. p-value`) # Sort p-values
  }

  results$`p-value2`      <- ifelse(results$`p-value`      < 0.001, "<.001", round(results$`p-value`,      3))
  results$`adj. p-value2` <- ifelse(results$`adj. p-value` < 0.001, "<.001", round(results$`adj. p-value`, 3))

  results$Significance <- case_when(
    as.numeric(results$`adj. p-value`) < 0.001 ~ "***",
    as.numeric(results$`adj. p-value`) < 0.01  ~ "**",
    as.numeric(results$`adj. p-value`) < 0.05  ~ "*",
    TRUE                                       ~ ""
  )

  comparativeAnalysisResults$results <- results

  # # Plotting identified metabolites using ggstatsplot
  # if (!is.null(plot_iden_met)) {
  #   for (metabolite in plot_iden_met) {
  #     if (metabolite %in% colnames(df)) {
  #
  #       # Set the test type dynamically
  #       test_type <- ifelse(results$`Test used` == "Mann-Whitney U test" | results$`Test used` == "Kruskal-Wallis",
  #                           "nonparametric", "parametric")
  #
  #       # Create the plot with dynamic test type
  #       library(ggstatsplot)
  #       p <- ggstatsplot::ggbetweenstats(
  #         data = data.frame(x = groups, y = df[[metabolite]]),
  #         x = "x", y = "y",
  #         type = test_type,
  #         pairwise = TRUE,
  #         p.adjust.method = "BH",  # Apply Benjamini-Hochberg adjustment
  #         title = paste("Comparative Analysis of", metabolite)
  #       )
  #     }
  #   }
  # }



  # Plotting identified metabolites using ggstatsplot
  if (!is.null(plot_iden_met)) {
    for (metabolite in plot_iden_met) {
      if (metabolite %in% rownames(results)) {  # Check if metabolite is in results

        # Get the test used for the metabolite
        test_used <- results[metabolite, "Test used"]

        # Set the test type dynamically based on the test used
        test_type <- ifelse(test_used %in% c("Mann-Whitney U test", "Kruskal-Wallis"),
                            "nonparametric", "parametric")

        # Create the plot with dynamic test type
        # library(ggstatsplot)
        plot_the_met <- NULL # Control, because if 'plot_the_met' will be returned as 'not found' otherwise
        plot_the_met <- ggstatsplot::ggbetweenstats(
          data = data.frame(x = groups, y = df[[metabolite]]),
          x    = "x",
          y    = "y",
          type = test_type,  # Use the dynamically determined test type
          # pairwise = TRUE,
          p.adjust.method = "BH",  # Apply Benjamini-Hochberg adjustment
          title = paste0("Comparative Analysis of ", metabolite)
        )


        if(!is.null(plot_the_met)) {

          comparativeAnalysisResults$plots[[metabolite]] <- plot_the_met

          print(plot_the_met)

        } else {


          next
        }

      } else {
        print(paste0("Metabolite '", metabolite, "' is not found in the preprocessed data (it has been removed in one of the data preprocessing steps). Boxplot is not generated."))

      }



    }
  }



  return(comparativeAnalysisResults)
}




# Usage
# myComparative <- performComparativeAnalysis(mydata)

# Add a "sort = TRUE" which sorts the results in descending order, i.e., lowest adjusted p-value first
