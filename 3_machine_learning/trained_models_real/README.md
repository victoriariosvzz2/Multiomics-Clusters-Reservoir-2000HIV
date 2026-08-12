# Real Trained Models

These are the actual fitted XGBoost models trained on the real 2000HIV cohort
data, used to generate the pairwise classification results and SHAP figures
reported in the manuscript. They are provided so reviewers can inspect and
directly verify the published results without needing to retrain anything.

## Files

Cluster label mapping: **0 = Mixed, 1 = All Low, 2 = All High**

| File | Comparison |
|---|---|
| `XGBoost_1_vs_0_2025-12-12_..._Fixed.pkl` | All Low vs Mixed |
| `XGBoost_1_vs_2_2025-12-12_..._Fixed.pkl` | All Low vs All High |
| `XGBoost_0_vs_2_2025-12-12_..._Fixed.pkl` | Mixed vs All High |

Each filename documents the exact preprocessing pipeline applied before
training: no data leakage between train/test, quantile-transformed features,
KNN-imputed missing values, viral reservoir features included.

## What's in each file

Each `.pkl` contains a single fitted `xgboost.sklearn.XGBClassifier` object
(saved with Python's `pickle`), including its learned hyperparameters
(`max_depth`, `learning_rate`, `n_estimators`, etc., selected via
`RandomizedSearchCV` as described in the manuscript Methods) and its
`feature_names_in_` attribute, listing the exact features it expects as
input.

**These files contain only the fitted model (learned parameters) — no
patient-level training data is included or recoverable from them.**

## How to load

```python
import pickle

with open("XGBoost_1_vs_0_2025-12-12_..._Fixed.pkl", mode="rb") as f:
    model_1_vs_0 = pickle.load(f)

print(model_1_vs_0.feature_names_in_)  # features this model expects
```

## How to use with the reviewer package

The main pipeline script (`../train_and_evaluate_model.py`) can load these
models directly and run SHAP analysis + evaluation against them, without
retraining:

```bash
python train_and_evaluate_model.py \
    --data /path/to/real_2000hiv_data.csv \
    --skip-training \
    --models-dir trained_models_real
```

Note: this requires the real 2000HIV dataset (see the main README's
"Reproduction instructions" for access via the Radboud Data Repository) —
the models alone cannot be run without matching input data. Running this
against the synthetic demo dataset will fail or produce meaningless output,
since the demo data's feature values don't correspond to what these models
were actually trained on.
