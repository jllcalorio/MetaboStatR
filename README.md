# Metabolomics-Data-Analysis

This repository contains R scripts which will analyze metabolomics data, where Samples are in the rows and Metabolites are in the columns.

The code does the ff.:

1.  Install and load R packages
2.  Load the data
3.  Pre-process the data
     -  Set groups and batches vectors
     -  Missing value filtering
     -  Missing value imputation
     -  Data normalization
     -  Data transformation
     -  Data scaling
5.  Batch correction
6.  Principal component analysis (PCA)
     -  Scree plots
     -  Scores plots
     -  Biplot
     -  Metabolites plot
7.  Partial least squares discriminant analysis (PLS-DA)
     -  Variable importance plot (VIP)
8.  Fold change analysis
9.  Statistics: Test for significant difference
     -  p-value correction
     -  p-value plots
     -  Volcano plot

# Format of data (required)

The data should be in a specified format for the R script to perform as intended.

