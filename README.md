# Diabetes Transfer Learning Across Diverse Populations

**Research code for the manuscript:**

*Improving Prediction of Patient-Reported Health Outcomes Among Older Adults With Diabetes Across Diverse Populations Through Transfer Learning*

## Overview

This repository contains R code for a research project investigating whether transfer learning can improve the prediction of patient-reported health outcomes among older adults with diabetes across diverse populations.

The study examines two complementary settings:

- **Cross-country transfer learning:** Leveraging information across eight international aging cohorts.
- **Within-country transfer learning:** Leveraging information across racial and ethnic groups within the U.S. Health and Retirement Study (HRS).

The analyses compare transfer learning with target-only and pooled regression approaches, evaluate predictive performance across repeated train/test splits, and examine model coefficients using bootstrap resampling.

## Study Populations and Outcomes

### Cross-Country Analysis

The cross-country analysis includes eight international aging cohorts:

| Cohort | Study |
|---|---|
| HRS | Health and Retirement Study |
| MHAS | Mexican Health and Aging Study |
| LASI | Longitudinal Aging Study in India |
| CHARLS | China Health and Retirement Longitudinal Study |
| ELSA | English Longitudinal Study of Ageing |
| KLOSA | Korean Longitudinal Study of Aging |
| HAALSI | Health and Aging in Africa: A Longitudinal Study of an INDEPTH Community in South Africa |
| SHARE | Survey of Health, Ageing and Retirement in Europe |

The outcome is **binary diabetes complication burden**, defined as the presence of at least one reported pain or vision complication.

Each cohort serves as a target population, with the remaining seven cohorts considered as candidate sources.

### Within-Country Analysis

The within-country analysis uses data from the 2018 Health and Retirement Study and considers four racial and ethnic groups:

- Hispanic
- Non-Hispanic White
- Non-Hispanic Black
- Non-Hispanic Other

The outcome is **self-reported eyesight**, dichotomized as poor/fair versus good/very good/excellent.

Each group serves as a target population, with the remaining groups considered as candidate sources.

## Statistical Methods

Three predictive approaches are compared:

| Method | Description |
|---|---|
| Baseline | Target-only ridge logistic regression |
| Pooled | Ridge logistic regression combining target and source observations with population indicators |
| Transfer | Transfer learning using the `glmtrans` R package |

### Model Evaluation

Predictive performance is evaluated using:

- Stratified 70/30 training and testing splits
- 300 repeated evaluations for each target–source combination
- Accuracy, area under the ROC curve (AUC), and Brier score
- Paired performance comparisons against baseline and pooled models

The best source combination for each target is selected based on the largest mean AUC improvement relative to the target-only baseline.

### Coefficient Inference

For the selected source combinations, transfer learning models are fitted using the full analysis datasets.

Coefficient uncertainty is assessed using 1,000 bootstrap resamples of the target population, with source data held fixed.

Results are reported on both the log-odds and odds-ratio scales, together with bootstrap standard errors and percentile confidence intervals.

## Repository Structure

```text
Diabetes-Transfer-Learning/
├── README.md
├── Cross_Country/
│   ├── 01_setup_and_detection.r
│   ├── 02_evaluation_hpc.r
│   └── 03_coefficient_inference.r
├── Within-Country_RAND-HRS/
│   ├── 01_setup_and_detection.r
│   ├── 02_evaluation_hpc.r
│   └── 03_coefficient_inference.r
└── Figures/
    └── Figure3.r
```

Both analysis directories follow the same general workflow:

| Script | Purpose |
|---|---|
| `01_setup_and_detection.r` | Data setup, source detection, and enumeration of candidate source combinations |
| `02_evaluation_hpc.r` | Repeated evaluation of baseline, pooled, and transfer learning models |
| `03_coefficient_inference.r` | Bootstrap coefficient inference and export of regression results |

### Figure Generation

`Figures/Figure3.r` generates Figure 3, comparing the relative AUC of pooled and transfer learning models against the target-only baseline across both study settings.

The plotted intervals represent the 2.5th and 97.5th percentiles of relative AUC across 300 repeated train/test splits.

## Software and Usage

The analyses were developed in **R 4.5.1**, using packages including `glmtrans`, `glmnet`, `dplyr`, `tidyr`, `pROC`, `openxlsx`, `readxl`, `kableExtra`, `ggplot2`, `purrr`, `tibble`, and `patchwork`.

Within each analysis directory, scripts should be executed in numerical order:

1. `01_setup_and_detection.r`
2. `02_evaluation_hpc.r`
3. `03_coefficient_inference.r`

The analysis scripts use relative paths and should be run from their respective analysis directories. The Figure 3 script should be run from the repository root after the evaluation results have been generated.

The analyses were configured for parallel processing with 16 CPU cores. The cross-country evaluation is computationally intensive and may benefit from a high-performance computing environment. Results involving internal cross-validation may vary slightly across computing environments and parallel settings.

## Data Availability

Individual-level research data are not included in this repository.

The analyses use data from international aging studies and the RAND Health and Retirement Study. Access to the underlying datasets is subject to the respective studies' data access requirements and usage conditions.

The repository provides the core statistical analysis and figure-generation scripts. Data cleaning and preparation scripts are not included, and the analyses require appropriately prepared input datasets.

## Citation

This repository accompanies the following manuscript:

Wang X, Chen Y, Shao W, Beasley JM, Lee TF. *Improving Prediction of Patient-Reported Health Outcomes Among Older Adults With Diabetes Across Diverse Populations Through Transfer Learning.* Manuscript in preparation.
