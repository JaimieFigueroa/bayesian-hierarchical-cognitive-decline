# Bayesian Modeling of Cognitive Aging in Centenarian Offspring

## Overview

This project applies Bayesian statistical modeling to longitudinal cognitive data from a study of offspring of centenarians and a control group.

The primary outcome was the Telephone Interview for Cognitive Status (TICS), measured repeatedly over time. The analysis investigated whether individuals with a familial history of exceptional longevity demonstrated different patterns of cognitive aging compared with controls.

This was a collaborative academic project completed as part of graduate-level biostatistical training. The project applied Bayesian regression, hierarchical modeling, latent trajectory modeling, and missing-data analysis using R and JAGS.

## Research Questions

The project addressed five main questions:

1. Do centenarian offspring and controls differ in baseline cognitive function?
2. Does baseline cognitive function differ between groups after adjustment for potential confounders?
3. Do centenarian offspring and controls differ in their rate of cognitive change over time?
4. Are there distinct subgroups of participants with different cognitive trajectories?
5. How sensitive are the longitudinal findings to missing-data mechanisms and imputation?

## Data

The analysis used longitudinal TICS data from a study of centenarian offspring and controls.

The dataset included:

- TICS scores at five timepoints (TICS01–TICS05)
- Age
- BMI
- Sex
- Years of education
- Smoking history
- Hypertension
- Diabetes
- Stroke
- Coronary artery disease
- Other clinical characteristics

Because the underlying dataset contains participant-level research data, the raw data are not included in this repository.

## Statistical Methods

### Baseline Analysis

Baseline characteristics were compared between centenarian offspring and controls using:

- Descriptive statistics
- t-tests
- Wilcoxon rank-sum tests
- Chi-square tests

### Bayesian Baseline Model

A Bayesian normal regression model was used to evaluate differences in baseline TICS scores after adjustment for potential confounders.

### Bayesian Hierarchical Longitudinal Model

A Bayesian hierarchical linear model was used to evaluate differences in the rate of TICS change between groups.

The model incorporated:

- Repeated TICS measurements
- Subject-specific random intercepts
- Time
- Group
- Group-by-time interaction
- Baseline TICS
- Age
- BMI
- Sex
- Smoking
- Hypertension
- Diabetes
- Stroke
- Coronary artery disease

### Latent Trajectory Model

A Bayesian two-class latent trajectory mixture model was used to identify subgroups of participants with distinct longitudinal patterns of cognitive change.

### Missing Data Analysis

Missingness patterns were examined using descriptive summaries and logistic regression.

Mean imputation was then used as a sensitivity analysis to evaluate how conclusions from the longitudinal model changed under a different missing-data approach.

## Key Findings

The analysis suggested that familial longevity was more strongly associated with **the trajectory of cognitive aging than with baseline cognitive status**.

After adjustment for potential confounders, baseline TICS scores did not show strong evidence of a difference between centenarian offspring and controls.

In contrast, the hierarchical longitudinal model suggested that centenarian offspring experienced a **slower rate of cognitive decline** than controls.

The latent trajectory analysis identified two distinct longitudinal patterns, including a relatively stable group and a group demonstrating greater cognitive decline.

The missing-data analysis demonstrated that longitudinal conclusions were sensitive to the method used to handle missing observations, highlighting the importance of principled missing-data methods in longitudinal cognitive-aging research.

## Repository Structure

```text
bayesian-cognitive-aging/
│
├── README.md
│
├── report/
│   └── final_report.html
│
├── code/
│   ├── q1_baseline_analysis.R
│   ├── q2_baseline_bayesian.R
│   ├── q3_hierarchical_model.R
│   ├── q4_latent_trajectory.R
│   └── q5_missing_data.R
│
└── figures/
    ├── Q1_Table1.png
    ├── Q2_Summary.png
    ├── Q2_Diagnosis.png
    ├── Q2_gelman-rubin.png
    ├── q3_summary.png
    ├── q3_gelman-rubin.png
    ├── q3_trace_density.png
    ├── q3_acf.png
    ├── Q4_Summary.png
    ├── Q4_ClassSeperation.png
    ├── Q4_Diagnosis.png
    ├── Q5_na_summary.png
    ├── Q5_model_mar_summary.png
    ├── Q5_test_hier_imp_summary.png
    └── Q5_gelman_diag.png
