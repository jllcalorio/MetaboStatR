# Metabolomics-Data-Analysis

This repository contains a set of R scripts designed to analyze metabolomics data, where each row represents a sample and each column represents a metabolite. The repository aims to provide a comprehensive workflow for the preprocessing, analysis, and visualization of metabolomics data, with a focus on handling multiple groups and batches.

A readable documentation can be found [here](https://jllcalorio.github.io/MetaboStatR/reference/index.html) that was built using [pkgdown](https://pkgdown.r-lib.org/).

# Objective

The primary goal of this repository is to facilitate the analysis of metabolomics data, ensuring reproducibility and transparency. It covers key steps from raw data processing to statistical analysis and visualization, making it suitable for users new to metabolomics or experienced data analysts.

# Workflow

1.  Install and load required R packages
     -  Automatically check for and install any missing packages, then load all necessary libraries.
2.  Load the metabolomics data
     -  Load your data from a .csv, .xlsx, or other formats and ensure it is in the correct format for analysis (see below for formatting requirements).
3. Data preprocessing
     -  Missing value filtering: Identify and remove features or samples with a certain proportion of missing values in each group.
     -  Missing value imputation: Use imputation methods to fill in missing values. In this case, 1/5 of smallest absolute value in the row.
4. Set group and batch vectors
     -  Define group labels (e.g., control vs. treatment).
     -  Assign batch identifiers for batch correction steps.
5. Data normalization
     -  Apply a normalization technique to account for differences in sample concentration or technical variability.
6. Data transformation
     -  Apply log 10 transformation to reduce skewness and make the data more suitable for downstream statistical analysis.
7. Data scaling
     -  Scale the data using techniques such as mean-centering, autoscaling (unit variance scaling), or Pareto scaling to standardize the range of metabolite abundances.
8. Batch correction
     -  Use methods such as ComBat to remove batch effects and harmonize the data across different experimental runs.
9. Principal Component Analysis (PCA)
     -  Perform PCA for dimensionality reduction and exploratory data analysis.
     -  Generate scree plots, score plots, and biplots for data visualization.
10. Partial Least Squares Discriminant Analysis (PLS-DA)
     -  Perform PLS-DA to identify metabolites that discriminate between predefined groups.
     -  Generate variable importance in projection (VIP) scores to rank metabolites based on their contribution to group separation.
11. Fold change analysis
     -  Calculate fold changes between different conditions or groups to assess the relative changes in metabolite levels.
12. Statistical analysis
     -  Perform statistical tests (e.g., t-tests, ANOVA) to identify significantly different metabolites between groups.
     -  Correct for multiple testing using methods such as Bonferroni or Benjamini-Hochberg false discovery rate (FDR).
13. Visualization
     -  p-value plots: Display the distribution of p-values from statistical tests.
     -  Volcano plot: Visualize the fold change vs. significance level for metabolites, highlighting those that are statistically significant.

# Data Format

To use this R script, the data must be formatted as follows:
-  Samples are arranged in rows, with each row representing one sample.
-  Metabolites are arranged in columns, with each column representing a single metabolite.
-  The first column should contain the sample identifiers.
-  Group and batch information must be included as separate columns in the data file towards the end.

**Example format:**

Sample | Metabolite 1 | Metabolite 2 | ... | Metabolite N | Group | Batch

# 📄 License

MetaboStatR is available under a dual license:

### Open Source License (MIT)
- ✅ Free for open source projects
- ✅ Free for research and educational use
- ✅ Community contributions welcome
- ✅ Modify and distribute freely

### Commercial License
- 💼 For commercial applications and products
- 💼 For proprietary software integration
- 💼 Includes priority support
- 💼 Custom licensing terms available

**Need a commercial license?** Contact me at [jllcalorio@gmail.com]

### Contributing
By contributing to this project, you agree that your contributions will be licensed under the same dual license terms.
