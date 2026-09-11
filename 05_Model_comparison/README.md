# Model comparison

本目录包含与模型扩展和临床性能比较相关的分析。

## Models

比较以下三个模型：

1. Microbiome
2. PRS-ERS
3. PRS-ERS & Microbiome

## Analyses

### PRS / ERS distributions

- ACRN vs HC
- scatter plot
- mean ± SD
- two-sample t-test

### ROC / AUC

- AUC
- 95% CI
- three-model ROC comparison

### Decision curve analysis

- Microbiome
- PRS-ERS
- PRS-ERS & Microbiome
- 50 bootstrap resamples

### NRI / IDI

Reference model:

`Microbiome`

Comparisons:

- PRS-ERS vs. Microbiome
- PRS-ERS & Microbiome vs. Microbiome

Continuous NRI and IDI are estimated with 100 bootstrap resamples.

### PPV / NPV

The original analysis used:

- Sensitivity = 0.783
- Specificity = 0.847
- assumed ACRN prevalence = 3%–15%

Selected prevalence scenarios:

- 3%
- 6%
- 9%
- 12%
- 15%

## Input

`data/model_comparison_input.csv`

Required columns:

- optional anonymized `SampleID`
- `Result`: 0/1 or HC/ACRN
- `PRS`
- `ERS`
- 47 microbiome genera

The 47 genera are identified as all remaining feature columns.

## Output

Results are saved to:

`results/model_comparison/`

Main outputs:

- `PRS_ERS_summary.csv`
- `PRS_ERS_t_test.csv`
- `PRS_ERS_distribution.pdf`
- `AUC_3Models.csv`
- `ROC_AUC_3Models.pdf`
- `DCA_3Models.pdf`
- `NRI_IDI_results.csv`
- `PPV_NPV_selected_prevalence.csv`
- `PPV_NPV_curve.pdf`

## Privacy

The public version does not contain real participant identifiers, real sample IDs,
local computer paths, server paths, or other identifiable information.