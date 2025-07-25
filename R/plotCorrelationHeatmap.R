#' Plot Correlational Heatmaps
#'
#' @description
#' This function plots correlational heatmaps based on several correlation and clustering methods.
#'
#' @param data List. This list must be a result from the `performPreprocessingPeakData` function. It can also be from a data scaling technique.
#' @param method String. Specifies the correlation method to be used.
#'   \itemize{
#'     \item "pearson": Uses Pearson r.
#'     \item "spearman": Uses Spearman rho.
#'     \item "kendall": Uses Kendall tau.
#'     }
#'     Defaults to "pearson".
#' @param plot_top_n Numeric. Specifies the top n features/samples to be plotted.
#' @param plot_what String. Specifies if Sample or Features are to be plotted.
#'   \itemize{
#'     \item "Samples": Plots the samples.
#'     \item "Features": Plots the features.
#'     }
#'     Defaults to "Features".
#' @param plot_type String. Specifies the type of plot.
#'   \itemize{
#'     \item "correlation": Correlational plot.
#'     \item "hierarchical": Hierarchical plot.
#'     }
#'     Defaults to "correlation".
#' @param show_rownames Boolean. If `TRUE` (default), shows the row names.
#' @param show_colnames Boolean. If `TRUE` (default), shows the column names.
#' @param clustering_distance_rows String. Specifies the clustering distance in the rows. See `?dist` for more info.
#'   \itemize{
#'     \item "euclidean": Usual distance between the two vectors.
#'     \item "maximum": Maximum distance between two components of `x` and `y` (supremum norm)
#'     \item "manhattan": Absolute distance between the two vectors (1 norm aka L_1). (L_1 = L subscript 1)
#'     \item "canberra": Terms with zero numerator and denominator are omitted from the sum and treated as if the values were missing. This is intended for non-negative values (e.g., counts), in which case the denominator can be written in various equivalent ways.
#'     \item "binary": (aka asymmetric binary): The vectors are regarded as binary bits, so non-zero elements are ‘on’ and zero elements are ‘off’. The distance is the proportion of bits in which only one is on amongst those in which at least one is on. This also called “Jaccard” distance in some contexts. Here, two all-zero observations have distance 0, whereas in traditional Jaccard definitions, the distance would be undefined for that case and give NaN numerically.
#'     \item "minkowski": The `p` norm, the `p-th` root of the sum of the `p-th` powers of the differences of the components.
#'     }
#'     Defaults to "euclidean".
#' @param clustering_distance_cols String. Specifies the clustering distance in the columns. Same as in `clustering_distance_rows`. Defaults to "euclidean".
#' @param clustering_method Specifies the clustering method. See `?hclust` for more info.
#'   \itemize{
#'     \item "ward.D": Ward's minimum variance method aims at finding compact, spherical clusters. Does not implement Ward's (1963) clustering criterion.
#'     \item "ward.D2" Implements that criterion (Murtagh and Legendre 2014). The dissimilarities are squared before cluster updating.
#'     \item "single": Closely related to the minimal spanning tree. Adopts a ‘friends of friends’ clustering strategy.
#'     \item "complete": Finds similar clusters.
#'     \item "average": Same as UPGMA.
#'     \item "mcquitty": Same as WPGMA.
#'     \item "median": Same as WPGMC. Does not lead to a monotone distance measure, or equivalently the resulting dendrograms can have so called inversions or reversals which are hard to interpret, but note the trichotomies in Legendre and Legendre (2012).
#'     \item "centroid": Same as UPGMC. Does not lead ... Legendre and Legendre (2012) (same as `median`).
#'     }
#'     Defaults to "ward.D".
#'
#' @returns A heatmap plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plotCorrelationHeatmap(data = results_from_performPreprocessingPeakData_function)
#' }
#'
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
