# PASC

This repository contains the R code used in the paper **"Label Efficient Phenotyping for Long COVID using Electronic Health Records"**, which explores a semi-supervised learning framework for phenotyping Long COVID using electronic health record (EHR) data. The code simulates a hybrid learning setting, where a surrogate model is trained using XGBoost on proxy labels (namely U09.9 status) and then combined with gold-standard labels via adaptive LASSO with leave-one-out cross-validation (LOOCV).

---

## Overview

The code pipeline consists of the following components:

1. **Helper Functions**: Various helper functions used throughout the rest of the code are defined. 
2. **Data Simulation**: Synthetic EHR features are generated with both surrogate labels (`u099.flag`) and true outcome labels (`Y`), including simulated missingness in the latter. This data simulates the preprocessed MVP data described in the Methods. 
3. **Surrogate Model Training**: An XGBoost model is trained under multiple cohort-specific settings (e.g., inpatient vs outpatient, time period variations), as described in step 1 of the "Semi-Supervised Phenotyping" section of the Methods. 
4. **Feature Integration via Adaptive LASSO**: Predictions from the surrogate model are aligned and combined with other EHR features. LOOCV is used with adaptive LASSO for supervised learning on the labeled subset. This corresponds to steps 2 and 3 of the "Semi-Supervised Phenotyping" section of the Methods.
5. **Evaluation**: Final predictions are evaluated using ROC curve metrics, including AUC, TPR, FPR, PPV, and NPV, as detailed in the "Methods for Comparison and Evaluation Metrics" section of the Methods. 

---

## Getting Started

### Prerequisites

This code requires R and the following packages:

```r
install.packages(c("xgboost", "pROC", "dplyr", "glmnet", "glmpath"))
```

### Running the Script

To run the entire pipeline:
```r
source("PASC_phenotype.R")
```

The script will:
- Simulate a dataset of 1000 patients and 50 features.
- Train multiple XGBoost models on surrogate labels (`u099.flag`) using different cohort settings.
- Align surrogate predictions with labeled data.
- Apply LOOCV with adaptive LASSO on the labeled subset (`Y`).
- Output ROC performance metrics including AUC.

### Output
The final output is a named list containing:

- `auc`: Area Under the ROC Curve
- `cut`: Classification cutoff value (e.g., 0.9)
- `p.pos`: Proportion of labeled patients classified as positive
- `fpr`: False Positive Rate
- `tpr`: True Positive Rate
- `ppv`: Positive Predictive Value
- `npv`: Negative Predictive Value
