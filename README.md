# survivalContour2 <img src="https://img.shields.io/badge/R-package-blue" style="height:20px"/>

**survivalContour2** is an advanced visualization toolkit for **time-dependent effect size**, **AUC**, and **decision-curve analysis** across follow-up time in survival models.  
It extends the original *survivalContour* framework with support for machine-learning survival models and provides smooth and interpretable contour surfaces.

## Key Features

### 1. Effect-Size Contours (HR / logHR)
- Time-dependent **hazard ratio** visualization  
- Full support for:
  - Cox proportional hazards (CoxPH)
  - Random Survival Forests (RSF)
  - XGBoost AFT models
  - PyCox deep survival models
  - Fine & Gray competing-risk models (CSC)

### 2. Time-Dependent AUC Contours
- Time-Dependent AUC over:
  - Time (x-axis)
  - Continuous covariate threshold (y-axis)
- Supports Cox / RSF / XGB

### 3. ADE (Average Derivative Effect) Contours
- Shows marginal change in survival probability  
- `ADEContour()` supports Cox, RSF, XGB

### 4. Decision Curve Analysis (DCA) for Continuous Covariates
- Contour version of net benefit
- Supports model vs All / None / Pairwise comparison

## Installation

### From GitHub (recommended)

```r
# install.packages("devtools")
devtools::install_github("HanzhiChen02/survivalContour2")
