# Pairwise Binary Classification of Immune-Reservoir Endotypes — Reviewer Package

This package accompanies the manuscript "Multi-Omics Clustering Differentiates
the Total and Intact HIV Reservoirs and Related Host Immune Mechanisms" (2000HIV
cohort, NCT03994835). It contains a self-contained, runnable version of the
machine learning pipeline used to train pairwise binary classifiers
distinguishing the three immune-reservoir endotypes identified by multi-omics
clustering, along with a synthetic demo dataset so the pipeline can be run and
inspected without access to patient-level data.

This package intentionally covers the **model training and evaluation** step
of the analysis only. Upstream steps (multi-omics clustering, differential
expression, GSEA, flow cytometry / cytokine regression models) are documented
separately in the full code repository (see "Code availability" below).

---

## 1. System requirements

**Software dependencies:**
- Python 3.12 (tested on 3.12.2, conda build)
- Package versions are pinned in `requirements.txt`, matching the exact
  environment used to produce the manuscript's results:
  - scikit-learn 1.5.1
  - xgboost 3.0.1
  - pandas 2.2.2
  - numpy 1.26.4
  - statsmodels 0.14.4
  - scipy 1.13.1
  - joblib 1.4.2
  - category_encoders 2.6.3
  - imbalanced-learn 0.13.0
  - shap 0.47.2
  - matplotlib 3.10.0


**If `pip install -r requirements.txt` reports a version conflict:** this
usually means your Python version is older than what a pinned package
requires (Python ≥3.10 is expected to work, but has not been independently
verified against every pinned version above). Run `python3 --version` to
check, and if you have multiple Python versions installed (e.g. via
Homebrew or Anaconda), create the virtual environment explicitly with a
matching version, e.g. `python3.12 -m venv venv` instead of
`python3 -m venv venv`.

**Operating systems tested:**
- Linux (Ubuntu 24.04, x86_64) — pip/venv installation path
- macOS (Python 3.12.2 via conda) —
  conda installation path recommended based on this experience;
  independently re-verified end-to-end in a clean environment by us

The code uses only standard Python/scikit-learn APIs and is not expected to be
platform-specific; it should also run on macOS and Windows with the same
package versions, though this has not been independently verified.

**Non-standard hardware:**
- None required. All models are small (tabular data, ≤1,230 samples, ≤hundreds
  of features) and train on CPU in seconds to minutes. No GPU is required.

---

## 2. Installation guide

**Instructions (conda — recommended):**

Since this pipeline depends on `shap`, which itself depends on `numba`/
`llvmlite` (packages that need to be compiled from source unless installed
via conda's pre-built binaries), installing via **conda** is strongly
recommended over plain `pip`/`venv`. A pip-only install can fail on some
systems with build errors for `numba`/`llvmlite` (missing `cmake`, missing
LLVM OpenMP runtime, etc.) — conda avoids this by installing pre-compiled
binaries directly.

```bash
# 1. Unzip this package and move into it
cd machine_learning

# 2. Create a conda environment with Python 3.12 (matching the environment
#    this pipeline was developed and tested on)
conda create -n reviewer_env python=3.12 -y
conda activate reviewer_env

# 3. Install all dependencies from requirements.txt via conda.
#    Conda resolves and installs numba/llvmlite (shap's own dependencies)
#    as pre-built binaries automatically -- this is what avoids the
#    compilation error.
conda install --file requirements.txt
```

If you don't have conda/Miniconda installed, get it from
[conda.io/miniconda](https://www.conda.io/miniconda.html) first.

**If any package fails to resolve via `--file`** (e.g. due to a name
mismatch between PyPI and conda): install that specific package by
name individually, e.g. `conda install shap=0.47.2 -y`, then
re-run the `--file` command for the rest.

**Instructions (pip/venv — alternative):** only use this if conda is not
available to you. This may fail with a build error for `numba`/`llvmlite`
on some systems (see above); if it does, either switch to the conda
instructions, or first run `xcode-select --install` (macOS) and
`brew install cmake libomp` before retrying.

```bash
cd machine_learning
python3.12 -m venv venv
source venv/bin/activate        # on Windows: venv\Scripts\activate
pip install -r requirements.txt
```

**If `python3.12` is not found:** check which Python versions are available
on your system with `ls /usr/local/bin/python3*` (macOS/Linux) or
`py --list` (Windows). If Python 3.12 is not available at all, install it
via [python.org](https://www.python.org/downloads/) or, on macOS,
`brew install python@3.12` — or switch to the conda instructions above,
which handle this automatically via `python=3.12`.

**Typical install time on a normal desktop computer:** approximately 2–5
minutes, depending on internet connection speed (dominated by downloading
xgboost and its dependencies).

---

## 3. Demo

The package includes a synthetic dataset (`demo/data/synthetic_demo_data.csv`,
300 simulated samples) that mimics the structure of the real dataset — the
same clinical, transcriptomic, proteomic, flow cytometry, methylation, ex
vivo cytokine, and viral reservoir feature names and approximate value
ranges — but contains **no real patient data**. A small artificial signal is
injected so the demo classifiers have something genuine to learn; the
resulting AUC values are not meaningful and should not be compared to the
manuscript's reported performance.

To regenerate this file (optional, not required to run the demo):
```bash
python demo/generate_synthetic_data.py
```

**Instructions to run the demo:**
```bash
python train_and_evaluate_model.py --data demo/data/synthetic_demo_data.csv --demo
```

This runs the full pipeline (preprocessing → VIF-based feature reduction →
hyperparameter search across 5 classifier types → best-model retraining per
class pair → evaluation with bootstrapped confidence intervals → SHAP
analysis) using a reduced hyperparameter search (5 iterations, 3-fold CV)
for a fast end-to-end test.

**Expected output:**
- Console log showing preprocessing steps, per-classifier AUC/F1 for each of
  the 3 pairwise comparisons, and the selected best-performing classifier
  type.
- Three trained model files saved to `trained_models/` (one per class pair:
  0 vs 2, 1 vs 0, 1 vs 2).
- `reports/model_performance_report_<date>.csv` — full hyperparameter search
  results for all classifiers and pairs.
- `reports/final_evaluation_auc.csv` — final train/test AUC with bootstrapped
  95% confidence intervals per pair.
- `shap_outputs/<pair>_beeswarm.png` — SHAP beeswarm plot per class pair, for
  visual comparison against the manuscript's published SHAP figures.
- `shap_outputs/<pair>_shap_values.csv` and `<pair>_feature_values.csv` — the
  numerical SHAP values and underlying feature values behind each beeswarm
  plot, so the plot's content can be checked numerically as well as visually.

**Expected run time for demo on a normal desktop computer:** under 2 minutes
(measured at ~20–30 seconds on a standard cloud CPU instance, including the
SHAP step).

**Checking your output against ours:** the folder `demo/expected_output/`
contains the actual files produced by a real run of the command above
(`console_log.txt`, `trained_models/`, `reports/`, `shap_outputs/`). After
running the demo yourself, compare your output against this folder:
- `reports/final_evaluation_auc.csv` — your AUC values should be reasonably
  close (exact digit-for-digit reproduction is not expected; see note in
  "Reproduction instructions" below on inherent non-determinism).
- `shap_outputs/*_beeswarm.png` — open both and confirm the same top-ranked
  features appear in a similar order.
- If your output looks substantially different (e.g. errors, empty files,
  wildly different AUCs), that's a sign something differs in your
  environment from what's pinned in `requirements.txt`.

---

## 4. Instructions for use

**To run on your own data:**
```bash
python train_and_evaluate_model.py --data /path/to/your_data.csv --demo
```

Your input CSV must:
- Have an ID column as the first column (used as the DataFrame index).
- Include a `cluster` column with integer class labels (e.g., 0, 1, 2).
- Include the remaining columns as features (numeric and/or categorical);
  column names do not need to match the synthetic dataset's, but the class 
  pairs defined in main() `([(1, 2), (1, 0), (0, 2)])` assume 3 classes 
  labeled 0, 1, 2, corresponding to the three immune-reservoir endotypes: 
  0 = Mixed, 1 = All Low, 2 = All High, adjust this list if your data has 
  different labels or a different number of classes.

Use `--full` instead of `--demo` to run with the same hyperparameter search
settings used in the manuscript (50 iterations, 10-fold CV per classifier per
pair) — this is substantially slower (expect tens of minutes depending on
dataset size and hardware) but matches the reported methodology exactly.

**To verify SHAP results against real, separately-provided trained models**
(without retraining anything): this script can load already-fitted model
files instead of training new ones, and run SHAP analysis directly on them.
This is the recommended way to check that the SHAP beeswarm plots and values
match the manuscript's published figures, using the actual models behind
those figures:

```bash
python train_and_evaluate_model.py \
    --data /path/to/data.csv \
    --skip-training \
    --models-dir /path/to/real_trained_models
```

`--models-dir` should point to a folder containing `.pkl` files named like
`<ModelType>_<a>_vs_<b>_<date>.pkl` (e.g. `XGBoost_1_vs_0_2025-12-12.pkl`) —
the same naming convention this script itself uses when saving models, and
the convention used for the manuscript's own trained models. `--data` should
be preprocessed-compatible data (the same column structure this script's
`preprocess_data()` produces); if run against the real cohort data with the
real trained models, the resulting `shap_outputs/*_beeswarm.png` plots and
`*_shap_values.csv` / `*_feature_values.csv` files should match the
manuscript's published SHAP figures.

### (Optional) Reproduction instructions

To reproduce the manuscript's reported pairwise classification results
exactly:

1. Obtain the real 2000HIV multi-omics dataset via the Radboud Data
   Repository (RDR) as described in the manuscript's Data Availability
   statement (DOIs: https://doi.org/10.34973/p96d-kz55,
   https://doi.org/10.34973/qk29-f305, https://doi.org/10.34973/2bf5-yc03),
   plus the flow cytometry collection (ID:
   ru.rumc.2000hiv_t0000321a_dsc_002), merged with clinical/demographic data
   into the same column structure as `data/synthetic_demo_data.csv`.
2. Run:
   ```bash
   python train_and_evaluate_model.py --data /path/to/real_2000hiv_data.csv --full
   ```
3. Random seed is fixed throughout (`RANDOM_SEED = 123`) for the train/test
   split and all hyperparameter searches, matching the original analysis.
4. Note: due to inherent non-determinism in some scikit-learn/XGBoost
   operations across different package versions/hardware (e.g., parallel
   execution order with `n_jobs=-1`), exact numerical reproduction of
   reported AUC values to the last decimal is not guaranteed, but results
   should be closely consistent with those reported in the `expected_output` folder.

---

## Code availability

The full analysis pipeline, including all upstream steps not covered by this
package (multi-omics clustering, differential expression, GSEA, flow
cytometry and cytokine regression models, SHAP interpretability analysis),
is available at:
https://github.com/victoriariosvzz2/Multiomics-Clusters-Reservoir-2000HIV

## Contact

For questions regarding this code package, please contact the first author: 
Victoria  Rios-Vazquez (Victoria.RiosVazquez@radboudumc.nl), or the 
corresponding authors: André van der Ven (Radboudumc) and Linos Vandekerckhove 
(Ghent University).
