plotCorrelationHeatmap <- function(
    data,
    method                   = c("pearson", "spearman", "kendall")[1],
    plot_top_n               = 1000,
    plot_what                = c("Samples", "Features")[2],
    plot_type                = c("correlation", "hierarchical")[1],
    # non_significant_as_white = TRUE,
    show_rownames            = TRUE,
    show_colnames            = TRUE,
    clustering_distance_rows = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")[1],
    clustering_distance_cols = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")[1],
    clustering_method        = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")[1]





) {

  pacman::p_load(Hmisc, pheatmap)

  # Contains all results
  listResults <- list()

  # Set data to be used
  if (plot_what == "Samples") {
    df <- t(data$data_scaledOPLSDA)
  } else if (plot_what == "Features") {
    df <- data$data_scaledOPLSDA
  } else {
    stop("plot_what must be 'Samples' or 'Features'")
  }



  # Find interquartile ranges
  iqr_values <- apply(df, 2, IQR, na.rm = TRUE) %>% sort(decreasing = TRUE) %>% as.data.frame() %>% `colnames<-`("IQR")

  # Top N features
  top_features <- rownames(iqr_values)[1:min(plot_top_n, length(t(iqr_values)))]

  # Subset the data to these top features
  filtered_df <- df[, top_features]

  ## ============ Correlation Heatmap ============
  if (plot_type == "correlation") {

    # Compute correlation matrix
    # cor_matrix <- cor(filtered_df, use = "pairwise.complete.obs", method = method) %>% as.data.frame()

    # Compute correlation matrix with p-values
    cor_results <- Hmisc::rcorr(x = as.matrix(filtered_df), type = method)
    cor_results.matrix <- cor_results$r
    cor_results.pvalues <- cor_results$P

    # Masking correlation matrix
    cor_results.matrix_masked <- cor_results.matrix
    cor_results.matrix_masked[cor_results.pvalues >= 0.05] <- 0  # Zero correlation => neutral color in diverging palette


    # Step 3: Define color palette
    # Colors: blue (negative), white (zero; non-significant, i.e., p >= .05), red (positive)
    my_palette <- colorRampPalette(c("blue", "white", "red"))(50)


    # Plot heatmap
    # heatmap.plot <- pheatmap::pheatmap(cor_matrix,
    #                                    show_rownames            = show_rownames,
    #                                    show_colnames            = show_colnames,
    #                                    clustering_distance_rows = clustering_distance_rows,
    #                                    clustering_distance_cols = clustering_distance_cols,
    #                                    clustering_method        = clustering_method,
    #                                    main                     = paste0("Correlation Heatmap (Top ", length(top_features)," by IQR)"),
    #                                    silent                   = TRUE)

    corr_heatmap.plot <- pheatmap::pheatmap(cor_results.matrix_masked,
                                       color = my_palette,
                                       breaks = seq(-1, 1, length.out = 51),  # ensure zero stays white
                                       display_numbers = FALSE, #round(cor_results.matrix, 2),  # show original correlation (even if non-significant)
                                       #number_color = "grey30",
                                       show_rownames = ifelse(plot_what == "Samples", TRUE, show_rownames),
                                       show_colnames = ifelse(plot_what == "Samples", TRUE, show_rownames),
                                       fontsize = 10,
                                       fontsize_row = 5,
                                       fontsize_col = 5,
                                       clustering_distance_rows = clustering_distance_rows,
                                       clustering_distance_cols = clustering_distance_cols,
                                       clustering_method = clustering_method,
                                       main = paste0("Correlation Heatmap of ", plot_what," (Non-sig p > 0.05 shown as white)")
                                       # ,silent = TRUE
    )

    listResults$CorrelationMatrix         <- cor_results.matrix
    listResults$CorrelationMatrixPValues  <- cor_results.pvalues
    listResults$CorrelationMatrixMasked   <- cor_results.matrix_masked
    listResults$CorrelationHeatmap        <- corr_heatmap.plot
  }


  ## ============ Hierarchical Clustering Heatmap ============
  if (plot_type == "hierarchical") {

    # Define palette for raw values (blue to white to red)
    my_palette2 <- colorRampPalette(c("blue", "white", "red"))(50)

    hier_heatmap.plot <- pheatmap::pheatmap(filtered_df,
                                            color = my_palette2,
                                            show_rownames = ifelse(plot_what == "Samples", TRUE, show_rownames),
                                            show_colnames = ifelse(plot_what == "Samples", TRUE, show_rownames),
                                            fontsize = 10,
                                            fontsize_row = 5,
                                            fontsize_col = 5,
                                            clustering_distance_rows = clustering_distance_rows,
                                            clustering_distance_cols = clustering_distance_cols,
                                            clustering_method = clustering_method,
                                            main = paste0("Hierarchical Clustering Heatmap of ", plot_what)
                                            # , silent = TRUE
    )

    listResults$HierarchicalHeatmap <- hier_heatmap.plot

    listResults$HierarchicalHeatmap       <- hier_heatmap.plot
  }



  # Save results in 1 list
  listResults$FunctionOrigin            <- "plotCorrelationHeatmap"
  listResults$CorrelationMethod         <- method
  listResults$TopN                      <- plot_top_n
  listResults$TopFeatures               <- top_features
  listResults$FilteredData              <- filtered_df %>% as.data.frame()
  listResults$IQRValues                 <- iqr_values

  listResults$ClusteringMethodInRows    <- clustering_distance_rows
  listResults$ClusteringMethodInColumns <- clustering_distance_cols
  listResults$ClusteringMethod          <- clustering_method




  return(listResults)




}


# Usage
# myheatmap <- plotCorrelationHeatmap(data = mypreprocess)
