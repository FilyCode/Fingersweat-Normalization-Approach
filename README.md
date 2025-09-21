# Fingersweat-Normalization-Approach

This repository provides an R-based framework for the normalization of omics datasets, specifically focusing on fingersweat metabolomics data. The approach integrates Bayesian optimization and combinatorial normalization techniques to improve data quality and enable more robust downstream analyses.

## Core Algorithm

The primary script for executing the normalization process is:

*   `normalization-algorithm.R`: This file orchestrates the Bayesian optimization and combinatorial normalization calculations. It also includes various options for visualizing the normalization results and assessing their impact.

## Supporting Functions

Essential functions called by the main algorithm are modularized for clarity and reusability:

*   **`general-files/feature-visualization-script.R`**: Contains functions dedicated to generating informative visualizations of features and normalization outcomes.
*   **`general-files/MS-data-analysis_functions.R`**: Houses general utility functions pertinent to MS data analysis.

## Data Preparation

Before applying the normalization algorithm, raw omics data files must be pre-processed into a standardized feature matrix and accompanied by relevant metadata.

### Required Input Files:

1.  **Feature Matrix:** A `.csv` file where samples are columns and features are rows. This file should ideally be supplemented with a SIRIUS summary file for feature annotation.
2.  **Metadata File:** A `.csv` file containing experimental metadata, with each row corresponding to a sample and columns providing relevant experimental details.

### Pre-processing Workflow (Raw Files to Feature Matrix):

The transformation from raw `.raw` files to the required feature matrix involves several key steps:

1.  **Raw to mzML Conversion:** Convert proprietary `.raw` files into the open-standard `.mzML` format.
2.  **Feature Detection & Quantification:** Process `.mzML` files using mzMine in Batch mode for robust feature detection and quantification.
3.  **SIRIUS Annotation:** Perform annotation of detected features using SIRIUS to identify molecular structures.

This entire pre-processing sequence can be automated using the script located at `raw-files-to-data-matrix-pipeline/get_untargeted_annotated_features_from_raw_files.R`.

Metadata files can be created or adapted using the `create_metadata-file.R` script (found in `various-other-files/`). This script will require modifications to align with the specific experimental design and metadata structure of your dataset.


## Authorship
This pipeline and its associated scripts were solely developed by Philipp Trollmann as part of his Master Thesis in Computational Science in the group of Dr. Jürgen Zanghellini and Dr. Samuel Meier-Menches at University of Vienna.
