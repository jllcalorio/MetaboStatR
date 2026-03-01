# Metabolomics-Data-Analysis

This repository contains a set of R scripts designed to analyze metabolomics data. It aims to provide a comprehensive workflow from data pre-preprocessing, statistical data analysis, and visualization, with a focus on handling multiple groups and batches, and injection sequences.

A readable documentation can be found [here](https://jllcalorio.github.io/MetaboStatR/reference/index.html) that was built using [pkgdown](https://pkgdown.r-lib.org/).

# Objective

The primary goal of this repository is to facilitate the analysis of metabolomics data, ensuring reproducibility and transparency. It covers key steps from raw data pre-processing to statistical analysis and visualization, making it suitable for users new to metabolomics or experienced data analysts.

# Workflow

1. Install and load required R packages
     -  Automatically check for and install any missing packages, then load all necessary libraries.
2. Load the metabolomics data and perform data quality check
     -  Load your data (a file location) of a .csv, .xlsx, or other formats and ensure it is in the correct format for analysis (see below for formatting requirements).
3. Data preprocessing
     -  Missing value filtering: Identify and remove features with a certain proportion of missing values across all groups (optional to include on QC groups).
     -  Missing value imputation: Use imputation methods to fill in missing values. The default and available for now is 1/n of smallest positive value per feature, where n is any integer (defaults to 5, which means missing values per feature are replaced with 1/5 of the smallest value in that feature).
4. Signal drift and batch correction
     -  Performs Quality Control-Robust Spline Correction (QC-RSC). 
5. Data normalization
     -  Apply a normalization technique to account for differences in sample concentration or technical variability.
6. Data transformation
     -  Apply a mathematical transformation to reduce skewness and make the data more suitable for downstream statistical analysis.
7. Data scaling
     -  Scale the data using techniques such as mean-centering, auto-scaling (unit variance scaling), or Pareto-scaling to standardize the range of metabolite abundances.
8. Feature filtering
     - Remove features with RSD (Relative Standard Deviation) >= 30%.
     - Remove features belonging to 10th percentile of lowest variability.
     - Optionally remove features that are not present in the data for NON-PLS (Partial Least Squares)-type analysis (usually auto-scaled data) and in PLS-type analysis (must be Pareto-scaled).
10. Principal Component Analysis (PCA)
     -  Perform PCA for dimensionality reduction and exploratory data analysis.
     -  Generate scree plots and score plots.
11. Partial Least Squares-Discriminant Analysis (PLS-DA)
     -  Perform Orthogonal PLS-DA to identify metabolites that discriminate between predefined groups.
     -  Can also perform PLS, PLS-DA, sPLS-DA.
     -  Generate variable importance in projection (VIP) scores to rank metabolites based on their contribution to group separation.
12. Fold change analysis
     -  Calculate fold changes between different conditions or groups to assess the relative changes in metabolite levels.
13. Comparative analysis
     -  Perform statistical tests (e.g., t-tests, ANOVA, or their non-parametric counterarts dynamically via assumptions testing) to identify significantly different metabolites between groups.
     -  Correct for multiple testing using methods such as Bonferroni or the default Benjamini-Hochberg a.k.a. false discovery rate (FDR).
14. Perform Area Under the Receiver Operating Characteristic Curve (AUROC).
    - This is to evaluate the ability of a binary classification model to distinguish between two groups.
14. Visualization
     -  Volcano plot: Visualize the fold change vs. comparative analysis p-values for features, highlighting those that are both clinically/biologicall and statistically significant features.

# Data Format

To use this R script, the CSV file must be formatted as follows:
-  1st row are the sample names, must be called "Sample". Sample names are suggested to be short yet informative to improve visualizations.
-  2nd row are the groupings of the samples, such as 'Has Disease,' 'Has no Disease,' 'EQC,' 'SQC', etc. These are later used as dependent variable in say Elastic Net Regression. This must be called "Group".
-  3rd row are the numeric batch numbers. This must be called "Batch".
-  4th row are the unique and numeric injection sequences. Can be of not in ascending order. This must be called "Injection".
-  5th row are the Subject IDs, which are only required when there are technical replicates. This must be named "SubjectID".
-  6th row are the technical replicates. Usually, these are non-unique identifiers to determine if 2 or more samples are just technical replicates. This must be named "Replicate".
-  7th row are the 2nd groupings. These are another groupings, like in 'Group' row, of the samples. This must be named "Group2".
-  8th row are the values to be used in data normalization step, typically specific gravities. This must be named "Normalization".
-  9th row are the values in case the dependent variable, in 2nd row 'Group', is not preferred, but a numeric dependent variable. This must be named "Response".
-  10th row and below are the metabolites/features. All must be in numeric format. Missing values are allowed. Metabolite/Feature names are also suggested to be short for visualization purposes.

**Example format:**

Sample | Sample_1 | Sample_2 | ... | QC_1 | QC_2 |
Group
Batch
Injection
Replicate
Group2
Normalization
Response
Feature_1
Feature_2
...
Feature_n

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
