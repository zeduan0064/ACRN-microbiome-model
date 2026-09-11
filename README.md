# ACRN Microbiome Model

Analysis code for the development and validation of a gut microbiome-based model for the early identification of advanced colorectal neoplasia (ACRN), including advanced adenoma (AA) and colorectal cancer (CRC).

## Study overview

This repository contains the main analysis and visualization code used in the study.

The microbiome analysis was based on fecal 16S rRNA gene sequencing. A microbiome panel was selected using iterative LASSO regression and subsequently evaluated using multiple machine-learning algorithms. Additional analyses included genus-pathway correlations, microbial co-occurrence networks, PRS/ERS model comparisons, decision curve analysis, reclassification metrics, and external validation.

## Repository structure

```text
ACRN_code/
├── 01_16S_processing/
│   ├── 01_QIIME2_pipeline.sh
│   ├── README.md
│   ├── manifest_example.tsv
│   └── metadata_example.tsv
├── 02_Feature_selection/
│   ├── 02_iterative_LASSO.R
│   └── README.md
├── 03_Correlation_network/
│   ├── 03_correlation_network_analysis.R
│   └── README.md
├── 04_Machine_learning/
│   ├── 04_six_machine_learning_models.R
│   └── README.md
├── 05_Model_comparison/
│   ├── 05_model_comparison.R
│   └── README.md
└── 06_External_validation/
    ├── 06_external_validation.R
    └── README.md
```

## Analysis modules

### 01. 16S processing

Processing of paired-end 16S rRNA V3-V4 sequencing data using QIIME 2.

Main software/resources:

- QIIME 2 v2024.2
- Greengenes2 database v2022.10
- DADA2
- 341F/805R primers

Only anonymized example manifest and metadata files are included. Raw sequencing data and participant-level metadata are not included in this repository.

### 02. Feature selection

Iterative LASSO feature selection.

Main procedure:

- CLR transformation of genus abundance
- 50 repeated LASSO runs
- 70% random subsampling in each iteration
- 10-fold cross-validation
- `lambda.min`
- ranking by selection frequency
- retention of the top 20% candidate features

### 03. Correlation and network analysis

Includes:

- genus-pathway Spearman correlation analysis
- FDR correction
- group-specific HC and ACRN correlation matrices
- microbial co-occurrence network analysis
- prevalence filtering
- Spearman correlation
- BH-FDR correction
- 1,000 bootstrap resamples
- bootstrap 95% confidence intervals

### 04. Machine-learning comparison

Comparison of six machine-learning algorithms:

1. Logistic regression
2. Neural network
3. Extreme gradient boosting (XGBoost)
4. Gradient boosting machine (GBM)
5. Support vector machine (SVM)
6. AdaBoost

The script includes:

- incremental 1-to-47 genus model evaluation
- repeated 10-fold cross-validation
- ROC/AUC comparison
- final 47-genus model comparison
- confusion matrices
- sensitivity, specificity, precision and F1 score
- DeLong tests comparing logistic regression with the other models

### 05. Model comparison

Comparison of microbiome, PRS/ERS, and combined models.

Includes:

- PRS and ERS distributions
- ROC/AUC comparison
- decision curve analysis
- continuous NRI and IDI
- PPV and NPV across assumed ACRN prevalence

### 06. External validation

External validation in two independent cohorts.

Includes:

- identification of available genera from the predefined microbiome panel
- median imputation for missing values
- removal of zero-variance features
- Z-score standardization
- multivariable logistic regression
- ROC/AUC and 95% CI
- Youden threshold
- sensitivity and specificity with 1,000 bootstrap resamples
- confusion matrices
- bootstrap feature importance based on median absolute standardized beta coefficients

## Software environment

The primary statistical analyses were conducted using:

- R v4.3.1
- QIIME 2 v2024.2
- Greengenes2 v2022.10
- PLINK v1.9

Major R packages include:

- glmnet
- compositions
- psych
- igraph
- ggraph
- caret
- pROC
- gbm
- kernlab
- xgboost
- adabag
- nnet
- rmda
- tidyverse
- ggplot2
- dplyr
- patchwork
- readxl

## Data availability and privacy

This repository contains analysis code only.

The following materials are intentionally excluded:

- personally identifiable information
- participant-level clinical data
- real sample identifiers
- raw FASTQ files
- genotype data
- private metadata
- local computer/server paths
- intermediate analysis outputs

Example identifiers, where provided, are anonymized.

16S rRNA sequencing data and genotyping data are available from the lead contact subject to the conditions described in the manuscript.

## Reproducibility

The scripts use project-relative paths such as:

```text
data/
results/
```

Users should prepare input files according to the README file provided in each analysis directory.

Because participant-level research data are not distributed with this repository, the scripts are intended to document and reproduce the analytical workflow when the corresponding authorized input data are available.

## Code availability

All original code used for the analyses and visualizations in this study is available in this GitHub repository.

Repository URL:

```text
https://github.com/USERNAME/ACRN-microbiome-model
```

Replace `USERNAME` with the GitHub account name after the repository is published.

## Citation

If this code is used, please cite the associated manuscript.

Manuscript citation information will be added after publication.