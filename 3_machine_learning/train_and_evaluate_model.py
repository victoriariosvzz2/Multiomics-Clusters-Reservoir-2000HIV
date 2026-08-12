"""
train_and_evaluate_model.py

Self-contained, runnable version of the pairwise binary classification
pipeline used to distinguish the three immune-reservoir endotypes (cluster
0, 1, 2) in the manuscript "Multi-Omics Clustering Differentiates the Total
and Intact HIV Reservoirs and Related Host Immune Mechanisms" (2000HIV
cohort, NCT03994835).

This script reproduces the core modeling steps end to end:
    1. Load data (clinical + omics features + cluster label)
    2. Stratified train/test split
    3. Preprocessing: KNN imputation, log-transform of selected clinical
       variables, quantile normalization of numerical features, mode
       imputation + one-hot encoding of categorical features
    4. Iterative VIF-based multicollinearity reduction
    5. Pairwise binary classifier training with hyperparameter search
       (5 classifier types evaluated; best-performing type retrained per
       pair) -- exactly as in the original analysis notebook
    6. Evaluation: train/test AUC with bootstrapped 95% confidence intervals
    7. SHAP analysis: beeswarm plots + underlying values per pair, computed
       with shap.Explainer(model, X) -- the same call used to generate the
       manuscript's published SHAP figures

Two run modes are provided via command-line flags (see "Instructions for
use" in the README):
    --demo          Fast settings (small search space, few CV folds) -- for
                     quickly confirming the pipeline runs end to end on the
                     provided synthetic dataset. This is the default.
    --full          Settings matching the original manuscript analysis
                     (n_iter=50, cv=10) -- slower, intended for use with the
                     real dataset (see README "Reproduction instructions").

Two model sources are supported:
    (a) Train from scratch (default) -- trains new models on --data and
        immediately runs SHAP analysis on them.
    (b) Load pre-trained models with --models-dir -- points at a folder of
        already-fitted .pkl model files (e.g. the REAL models trained on
        the real cohort, provided separately alongside this package) and
        skips training entirely, going straight to evaluation + SHAP
        analysis using those real fitted models. This lets a reviewer (or
        you) confirm the real, published models' SHAP outputs without
        needing the real training data itself.

Usage:
    # Train from scratch on the synthetic demo data, then run SHAP:
    python train_and_evaluate_model.py --data data/synthetic_demo_data.csv --demo

    # Full manuscript settings on your own real data:
    python train_and_evaluate_model.py --data /path/to/your_data.csv --full

    # Use REAL pre-trained models (no training), running SHAP + evaluation
    # against whatever --data you point it at (e.g. the real dataset, run
    # locally by you -- not distributed with this reviewer package):
    python train_and_evaluate_model.py --data /path/to/real_data.csv \\
        --models-dir /path/to/real_trained_models --skip-training
"""

import argparse
import glob
import os
import warnings
from collections import defaultdict
from datetime import date

import joblib
import matplotlib
matplotlib.use("Agg")  # headless-safe backend for saving plots without a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from category_encoders import OneHotEncoder
from scipy.stats import randint, uniform
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (confusion_matrix, f1_score, precision_score,
                              recall_score, roc_auc_score)
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import QuantileTransformer
from sklearn.svm import SVC
from sklearn.utils import resample
from statsmodels.stats.outliers_influence import variance_inflation_factor
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")
RANDOM_SEED = 123
np.random.seed(RANDOM_SEED)


# =============================================================================
# 1. Preprocessing
# =============================================================================
def preprocess_data(X_train, X_test, log_scale_keywords=None, vif_threshold=10,
                     protected_features=None, corr_threshold=0.8, verbose=True):
    """
    Reproduces the preprocessing pipeline from the original analysis notebook:
      - split into numeric / categorical columns
      - KNN-impute numeric columns
      - log1p-transform selected clinical variables
      - quantile-normalize all numeric columns
      - iterative VIF-based multicollinearity reduction (protecting a
        pre-specified set of biologically important features)
      - drop features from remaining high-correlation pairs (|r| > threshold)
      - mode-impute categorical columns, one-hot encode
    """
    if log_scale_keywords is None:
        log_scale_keywords = ["CD4_NADIR", "VL_ZENITH", "CD4_LATEST",
                               "CD4CD8_LATEST", "CMV_IgG_IU/mL", "exposure"]
    if protected_features is None:
        protected_features = []

    num_cols = X_train.select_dtypes(include=[np.number]).columns.tolist()
    cat_cols = X_train.select_dtypes(exclude=[np.number]).columns.tolist()

    if verbose:
        print(f"Numeric columns: {len(num_cols)} | Categorical columns: {len(cat_cols)}")

    # --- Numeric: KNN impute ---
    if num_cols:
        imputer = KNNImputer(n_neighbors=5)
        X_train_num = pd.DataFrame(imputer.fit_transform(X_train[num_cols]),
                                    columns=num_cols, index=X_train.index)
        X_test_num = pd.DataFrame(imputer.transform(X_test[num_cols]),
                                   columns=num_cols, index=X_test.index)

        # --- Log-transform selected columns ---
        log_cols = [c for c in num_cols if any(k in c for k in log_scale_keywords)]
        if verbose and log_cols:
            print(f"Log-scaling: {log_cols}")
        for c in log_cols:
            X_train_num[c] = np.log1p(X_train_num[c].clip(lower=0))
            X_test_num[c] = np.log1p(X_test_num[c].clip(lower=0))

        # --- Quantile normalization ---
        qt = QuantileTransformer(output_distribution="normal", random_state=0)
        X_train_num = pd.DataFrame(qt.fit_transform(X_train_num),
                                    columns=num_cols, index=X_train.index)
        X_test_num = pd.DataFrame(qt.transform(X_test_num),
                                   columns=num_cols, index=X_test.index)
    else:
        X_train_num = pd.DataFrame(index=X_train.index)
        X_test_num = pd.DataFrame(index=X_test.index)

    # --- Iterative VIF reduction ---
    if len(X_train_num.columns) > 1:
        X_vif = X_train_num.copy()
        dropped = []
        while True:
            vif_vals = pd.Series(
                [variance_inflation_factor(X_vif.values, i) for i in range(X_vif.shape[1])],
                index=X_vif.columns,
            )
            if vif_vals.max() <= vif_threshold or len(X_vif.columns) <= 1:
                break
            sorted_vif = vif_vals.sort_values(ascending=False)
            candidate = next((f for f in sorted_vif.index if f not in protected_features), None)
            if candidate is None:
                break
            if verbose:
                print(f"Dropping '{candidate}' (VIF={vif_vals[candidate]:.2f})")
            dropped.append(candidate)
            X_vif = X_vif.drop(columns=[candidate])

        # --- High-correlation pair pruning on remaining features ---
        corr_matrix = X_vif.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        high_corr_drop = [col for col in upper.columns
                           if any(upper[col] > corr_threshold) and col not in protected_features]

        features_to_drop = dropped + high_corr_drop
        X_train_num = X_train_num.drop(columns=[c for c in features_to_drop if c in X_train_num.columns])
        X_test_num = X_test_num.drop(columns=[c for c in features_to_drop if c in X_test_num.columns])
        if verbose:
            print(f"Dropped {len(features_to_drop)} collinear feature(s) total.")

    # --- Categorical: mode-impute + one-hot encode ---
    if cat_cols:
        cat_imputer = SimpleImputer(strategy="most_frequent")
        X_train_cat = pd.DataFrame(cat_imputer.fit_transform(X_train[cat_cols]),
                                    columns=cat_cols, index=X_train.index)
        X_test_cat = pd.DataFrame(cat_imputer.transform(X_test[cat_cols]),
                                   columns=cat_cols, index=X_test.index)
        encoder = OneHotEncoder(cols=cat_cols, use_cat_names=True)
        X_train_cat = encoder.fit_transform(X_train_cat)
        X_test_cat = encoder.transform(X_test_cat)
    else:
        X_train_cat = pd.DataFrame(index=X_train.index)
        X_test_cat = pd.DataFrame(index=X_test.index)

    X_train_processed = pd.concat([X_train_num, X_train_cat], axis=1)
    X_test_processed = pd.concat([X_test_num, X_test_cat], axis=1)

    return X_train_processed, X_test_processed


# =============================================================================
# 2. Pairwise binary classifier training
# =============================================================================
def get_classifiers_and_grids(fast=True):
    classifiers = {
        "Random Forest": RandomForestClassifier(random_state=42, n_jobs=1),
        "XGBoost": XGBClassifier(eval_metric="logloss", random_state=42),
        "Logistic Regression": LogisticRegression(random_state=42, max_iter=1000, solver="liblinear"),
        "KNN": KNeighborsClassifier(),
        "SVM": SVC(random_state=42, probability=True),
    }
    param_distributions = {
        "Random Forest": {
            "max_depth": randint(5, 20) if fast else randint(10, 20),
            "n_estimators": randint(50, 200) if fast else randint(100, 500),
            "max_features": ["sqrt", "log2"],
            "min_samples_split": randint(2, 11),
            "min_samples_leaf": randint(1, 5),
            "bootstrap": [True, False],
        },
        "XGBoost": {
            "max_depth": randint(3, 10),
            "learning_rate": uniform(0.01, 0.3),
            "n_estimators": randint(50, 200) if fast else randint(100, 500),
            "subsample": uniform(0.5, 0.5),
            "colsample_bytree": uniform(0.5, 0.5),
            "gamma": uniform(0, 5),
        },
        "Logistic Regression": {"C": uniform(0.1, 10)},
        "KNN": {"n_neighbors": randint(3, 15), "weights": ["uniform", "distance"]},
        "SVM": {"C": uniform(0.1, 10), "kernel": ["rbf"], "gamma": ["scale"]},
    }
    return classifiers, param_distributions


def train_pairwise_models(X_train, X_test, y_train, y_test, pairs,
                           n_iter=50, cv=10, output_dir="trained_models",
                           reports_dir="reports", fast=True, extra_identifier=""):
    """
    Trains binary classifiers for each class pair, selects the best-performing
    classifier type by average AUC across pairs, retrains that type per pair,
    and saves the fitted models. Mirrors the original notebook's
    select_and_train_binary_classification_model() logic.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    classifiers, param_distributions = get_classifiers_and_grids(fast=fast)
    model_auc = defaultdict(list)
    model_f1 = defaultdict(list)
    all_results = []

    for pair in pairs:
        print(f"\n--- Evaluating classifiers for {pair[0]} vs {pair[1]} ---")
        train_mask = y_train.isin(pair)
        test_mask = y_test.isin(pair)
        X_tr = X_train[train_mask]
        y_tr = y_train[train_mask].apply(lambda v: 0 if v == pair[0] else 1)
        X_te = X_test[test_mask]
        y_te = y_test[test_mask].apply(lambda v: 0 if v == pair[0] else 1)

        if len(X_tr) < 10 or len(X_te) < 2:
            print(f"Skipping {pair}: insufficient data ({len(X_tr)} train, {len(X_te)} test)")
            continue

        for name, clf in classifiers.items():
            try:
                search = RandomizedSearchCV(
                    clf, param_distributions[name], n_iter=n_iter, scoring="roc_auc",
                    cv=min(cv, y_tr.value_counts().min()), n_jobs=1, random_state=RANDOM_SEED,
                )
                search.fit(X_tr, y_tr)
                y_prob = search.best_estimator_.predict_proba(X_te)[:, 1]
                y_pred = search.best_estimator_.predict(X_te)
                auc = roc_auc_score(y_te, y_prob)
                f1 = f1_score(y_te, y_pred)
                model_auc[name].append(auc)
                model_f1[name].append(f1)
                all_results.append({
                    "pair": f"{pair[0]} vs {pair[1]}", "model": name, "auc": auc,
                    "F1-Score": f1, "Precision": precision_score(y_te, y_pred),
                    "Recall": recall_score(y_te, y_pred),
                    "Confusion Matrix": str(confusion_matrix(y_te, y_pred)),
                    "best_params": search.best_params_,
                })
                print(f"  {name}: AUC={auc:.3f}, F1={f1:.3f}")
            except Exception as e:
                print(f"  {name} failed: {e}")

    avg_auc = {k: float(np.mean(v)) for k, v in model_auc.items() if v}
    if not avg_auc:
        raise RuntimeError("No classifiers trained successfully -- check input data.")
    best_type = max(avg_auc, key=avg_auc.get)
    print(f"\nBest-performing classifier type (avg AUC across pairs): {best_type} ({avg_auc[best_type]:.3f})")

    trained_models = {}
    for pair in pairs:
        train_mask = y_train.isin(pair)
        test_mask = y_test.isin(pair)
        X_tr = X_train[train_mask]
        y_tr = y_train[train_mask].apply(lambda v: 0 if v == pair[0] else 1)
        if len(X_tr) < 10:
            continue

        search = RandomizedSearchCV(
            classifiers[best_type], param_distributions[best_type], n_iter=n_iter,
            scoring="roc_auc", cv=min(cv, y_tr.value_counts().min()), n_jobs=1,
            random_state=RANDOM_SEED,
        )
        search.fit(X_tr, y_tr)
        model_path = os.path.join(
            output_dir, f"{best_type.replace(' ', '_')}_{pair[0]}_vs_{pair[1]}_{date.today()}{extra_identifier}.pkl"
        )
        joblib.dump(search.best_estimator_, model_path)
        trained_models[f"{pair[0]}_vs_{pair[1]}"] = {
            "model": search.best_estimator_, "hyperparameters": search.best_params_, "file": model_path,
        }
        print(f"Saved model for {pair[0]} vs {pair[1]} -> {model_path}")

    results_df = pd.DataFrame(all_results)
    results_path = os.path.join(reports_dir, f"model_performance_report_{date.today()}.csv")
    results_df.to_csv(results_path, index=False)
    print(f"\nPerformance report saved to {results_path}")

    return trained_models, best_type, avg_auc, results_df


# =============================================================================
# 3. Evaluation with bootstrapped CIs
# =============================================================================
def evaluate_models(trained_models, X_train, X_test, y_train, y_test, pairs, n_boot=1000):
    print("\n--- Final evaluation (train vs test AUC, bootstrapped 95% CI) ---")
    rows = []
    for pair in pairs:
        key = f"{pair[0]}_vs_{pair[1]}"
        if key not in trained_models:
            continue
        model = trained_models[key]["model"]
        features = model.feature_names_in_

        train_mask = y_train.isin(pair).values
        test_mask = y_test.isin(pair).values
        Xtr = X_train[features].iloc[train_mask].reset_index(drop=True)
        ytr = y_train.iloc[train_mask].apply(lambda v: 0 if v == pair[0] else 1).reset_index(drop=True)
        Xte = X_test[features].iloc[test_mask].reset_index(drop=True)
        yte = y_test.iloc[test_mask].apply(lambda v: 0 if v == pair[0] else 1).reset_index(drop=True)

        train_auc = roc_auc_score(ytr, model.predict_proba(Xtr)[:, 1])
        y_prob_te = model.predict_proba(Xte)[:, 1]
        test_auc = roc_auc_score(yte, y_prob_te)

        boot_aucs = []
        for b in range(n_boot):
            idx = resample(range(len(yte)), random_state=b)
            if len(set(yte.iloc[idx])) < 2:
                continue
            boot_aucs.append(roc_auc_score(yte.iloc[idx], y_prob_te[idx]))
        ci_lo, ci_hi = np.percentile(boot_aucs, [2.5, 97.5]) if boot_aucs else (np.nan, np.nan)

        print(f"{pair[0]} vs {pair[1]}: Train AUC={train_auc:.3f} | Test AUC={test_auc:.3f} "
              f"(95% CI: {ci_lo:.3f}-{ci_hi:.3f})")
        rows.append({"pair": f"{pair[0]} vs {pair[1]}", "train_auc": train_auc,
                      "test_auc": test_auc, "ci_lo": ci_lo, "ci_hi": ci_hi})
    return pd.DataFrame(rows)


# =============================================================================
# 4. Load pre-trained (e.g. REAL) models instead of training
# =============================================================================
def load_pretrained_models(models_dir, pairs):
    """
    Loads already-fitted model files from models_dir, matching filenames of
    the form "<AnyModelType>_<a>_vs_<b>_*.pkl" (the naming convention used by
    train_pairwise_models() / the original notebook's saved models). If
    multiple files match a given pair, the most recently modified one is
    used.

    This is how you point the script at your REAL trained models (e.g. the
    ones actually used to generate the manuscript's figures) instead of
    training new ones on synthetic data.
    """
    trained_models = {}
    for pair in pairs:
        pattern = os.path.join(models_dir, f"*_{pair[0]}_vs_{pair[1]}_*.pkl")
        matches = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)
        if not matches:
            print(f"  [warning] No model file found for pair {pair} matching '{pattern}' -- skipping.")
            continue
        model_path = matches[0]
        try:
            model = joblib.load(model_path)
        except Exception:
            import pickle
            with open(model_path, "rb") as f:
                model = pickle.load(f)
        trained_models[f"{pair[0]}_vs_{pair[1]}"] = {"model": model, "file": model_path}
        print(f"Loaded pre-trained model for {pair[0]} vs {pair[1]} <- {model_path}")

    if not trained_models:
        raise RuntimeError(f"No pre-trained model files found in {models_dir}. "
                            f"Expected filenames like 'XGBoost_1_vs_0_<date>.pkl'.")
    return trained_models


# =============================================================================
# 5. SHAP analysis -- same explainer call used to generate the manuscript's
#    published beeswarm figures: shap.Explainer(model, X) then explainer(X).
# =============================================================================
def run_shap_analysis(trained_models, X, pairs, output_dir="shap_outputs"):
    """
    For each pair's trained model, computes SHAP values on X (restricted to
    that model's own feature set), saves:
      - a beeswarm plot (PNG) for visual comparison against the manuscript's
        published SHAP figures
      - the underlying SHAP values and feature values (CSV) so the beeswarm
        plot's content can be inspected/verified numerically, not just
        visually

    IMPORTANT: uses shap.Explainer(model, X) + explainer(X) (interventional
    perturbation, background data = X) -- NOT shap.TreeExplainer(model)
    without background data, which is a different algorithm and can produce
    different feature rankings. This matches what generated the published
    figures.
    """
    os.makedirs(output_dir, exist_ok=True)

    for pair in pairs:
        key = f"{pair[0]}_vs_{pair[1]}"
        if key not in trained_models:
            continue
        model = trained_models[key]["model"]

        if not hasattr(model, "feature_names_in_"):
            print(f"  [warning] Model for {key} has no feature_names_in_ -- skipping SHAP.")
            continue

        features = [f for f in model.feature_names_in_ if f in X.columns]
        missing = set(model.feature_names_in_) - set(features)
        if missing:
            print(f"  [warning] {key}: {len(missing)} feature(s) expected by the model are "
                  f"not present in the provided data and will be skipped: {sorted(missing)[:5]}...")
        X_subset = X[features]

        print(f"\nRunning SHAP analysis for {pair[0]} vs {pair[1]} ({len(features)} features, "
              f"{X_subset.shape[0]} samples)...")

        explainer = shap.Explainer(model, X_subset)
        explanation = explainer(X_subset)
        values = explanation.values

        if values.ndim == 3:
            # Multiclass-shaped output from a binary model -- select the
            # positive-class slice (index 1) to match a standard binary
            # beeswarm plot.
            values = values[:, :, 1]

        # --- Beeswarm plot ---
        plt.figure()
        shap.summary_plot(values, X_subset, show=False, plot_size=(8, max(4, 0.3 * len(features))))
        plot_path = os.path.join(output_dir, f"{key}_beeswarm.png")
        plt.savefig(plot_path, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  Saved beeswarm plot: {plot_path}")

        # --- Underlying values (for numerical verification, not just visual) ---
        shap_df = pd.DataFrame(values, columns=features)
        shap_df.insert(0, "Sample_ID", X_subset.index)
        shap_path = os.path.join(output_dir, f"{key}_shap_values.csv")
        shap_df.to_csv(shap_path, index=False)

        feature_df = X_subset.reset_index().rename(columns={X_subset.index.name or "index": "Sample_ID"})
        feature_path = os.path.join(output_dir, f"{key}_feature_values.csv")
        feature_df.to_csv(feature_path, index=False)

        print(f"  Saved SHAP values: {shap_path}")
        print(f"  Saved feature values: {feature_path}")


# =============================================================================
# 6. Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="Train/evaluate pairwise binary classification models and run SHAP analysis.")
    parser.add_argument("--data", type=str, default="demo/data/synthetic_demo_data.csv",
                         help="Path to input CSV (must include a 'cluster' column and an ID index).")
    parser.add_argument("--demo", action="store_true", default=True,
                         help="Fast settings for quick end-to-end testing (default).")
    parser.add_argument("--full", action="store_true",
                         help="Full settings matching the manuscript's analysis (n_iter=50, cv=10). "
                              "Overrides --demo.")
    parser.add_argument("--output-dir", type=str, default="trained_models")
    parser.add_argument("--reports-dir", type=str, default="reports")
    parser.add_argument("--models-dir", type=str, default=None,
                         help="Path to a folder of already-fitted .pkl model files (e.g. the REAL "
                              "trained models). If given, use --skip-training to load these instead "
                              "of training new models.")
    parser.add_argument("--skip-training", action="store_true",
                         help="Skip training entirely and load pre-trained models from --models-dir. "
                              "Requires --models-dir.")
    parser.add_argument("--no-shap", action="store_true",
                         help="Skip the SHAP analysis step (it runs by default).")
    parser.add_argument("--shap-output-dir", type=str, default="shap_outputs")
    parser.add_argument("--shap-data", choices=["train", "test", "all"], default="train",
                         help="Which split to compute SHAP values on (default: train, matching the "
                              "manuscript's SHAP figures).")
    args = parser.parse_args()

    if args.skip_training and not args.models_dir:
        parser.error("--skip-training requires --models-dir to be set.")

    fast = not args.full
    n_iter = 5 if fast else 50
    cv = 3 if fast else 10

    print(f"Run mode: {'DEMO (fast)' if fast else 'FULL (manuscript settings)'}")
    if not args.skip_training:
        print(f"n_iter={n_iter}, cv={cv}\n")

    # --- 1. Load data ---
    df = pd.read_csv(args.data, index_col=0)
    if "cluster" not in df.columns:
        raise ValueError("Input data must contain a 'cluster' column with integer class labels.")
    X = df.drop(columns=["cluster"])
    y = df["cluster"]
    print(f"Loaded data: {X.shape[0]} samples, {X.shape[1]} features")
    print(f"Class distribution:\n{y.value_counts()}\n")

    # --- 2. Train/test split ---
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, stratify=y, random_state=RANDOM_SEED
    )
    print(f"Train: {X_train.shape}, Test: {X_test.shape}\n")

    # --- 3. Preprocessing ---
    protected_features = ["AGE", "SEX_BIRTH", "cg13452062", "CD4_NADIR",
                           "Total.million.avg", "intactDT.DSI"]
    X_train_processed, X_test_processed = preprocess_data(
        X_train, X_test, protected_features=protected_features
    )
    print(f"\nProcessed feature count: {X_train_processed.shape[1]}\n")

    class_pairs = [(1, 2), (1, 0), (0, 2)]  # 0=Mixed, 1=All Low, 2=All High

    # --- 4. Get models: either train new ones, or load real pre-trained ones ---
    if args.skip_training:
        print(f"\nSkipping training -- loading pre-trained models from {args.models_dir}")
        trained_models = load_pretrained_models(args.models_dir, class_pairs)
    else:
        trained_models, best_type, avg_auc, results_df = train_pairwise_models(
            X_train_processed, X_test_processed, y_train, y_test, pairs=class_pairs,
            n_iter=n_iter, cv=cv, output_dir=args.output_dir, reports_dir=args.reports_dir,
            fast=fast,
        )

    # --- 5. Evaluation ---
    os.makedirs(args.reports_dir, exist_ok=True)
    eval_df = evaluate_models(trained_models, X_train_processed, X_test_processed,
                               y_train, y_test, pairs=class_pairs)
    eval_path = os.path.join(args.reports_dir, "final_evaluation_auc.csv")
    eval_df.to_csv(eval_path, index=False)
    print(f"\nFinal evaluation table saved to {eval_path}")

    # --- 6. SHAP analysis ---
    if not args.no_shap:
        if args.shap_data == "train":
            X_shap = X_train_processed
        elif args.shap_data == "test":
            X_shap = X_test_processed
        else:
            X_shap = pd.concat([X_train_processed, X_test_processed], axis=0)

        run_shap_analysis(trained_models, X_shap, pairs=class_pairs, output_dir=args.shap_output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
