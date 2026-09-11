# External validation

This directory contains the external validation analysis for two independent cohorts.

## Core workflow

For each cohort:

1. Read the external validation dataset
2. Recode the outcome as HC or ACRN
3. Detect the available genera from the predefined microbiome panel
4. Convert genus abundance to numeric values
5. Impute missing values using the cohort-specific median
6. Remove zero-variance genera
7. Z-score standardize the available genera
8. Fit a multivariable logistic regression model
9. Calculate ROC, AUC and 95% CI
10. Select the classification threshold using the Youden index
11. Calculate sensitivity and specificity with 1000 bootstrap resamples
12. Generate the confusion matrix
13. Estimate feature importance using 1000 bootstrap logistic regressions

## Bootstrap feature importance

For every bootstrap resample:

- samples are drawn with replacement
- genus features are standardized again within the bootstrap sample
- the full multivariable logistic regression is refitted
- standardized beta coefficients are stored

Feature importance is defined as:

`median(abs(standardized beta))`

The public version reports:

- median absolute standardized beta
- bootstrap 95% interval
- median signed standardized beta

## Input

Place the two anonymized external validation files in:

`data/`

with the names:

- `external_cohort_1.xlsx`
- `external_cohort_2.xlsx`

Each file should contain:

- `Result`
- the available genus features from the predefined panel, with names beginning with `g__`

Outcome coding supported by the script:

- `0` / `HC` = healthy control
- `1` / `AA` = advanced adenoma
- `2` / `CRC` = colorectal cancer
- AA and CRC are combined as ACRN

Real participant identifiers should not be included in the public files.

## Output

Results are written to:

`results/external_validation/`

Main outputs:

- `External_validation_performance.csv`
- `Cohort1_bootstrap_importance.csv`
- `Cohort2_bootstrap_importance.csv`
- `Cohort1_confusion_matrix.csv`
- `Cohort2_confusion_matrix.csv`
- `External_validation_summary.pdf`

## Privacy

The public code does not contain real sample identifiers, participant information,
local computer paths, server paths, or other identifiable information.