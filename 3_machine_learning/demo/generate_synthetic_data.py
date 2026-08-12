"""
generate_synthetic_data.py

Generates a fully synthetic dataset that mimics the structure (column names,
value ranges, and approximate class balance) of the real 2000HIV multi-omics
dataset used to train the pairwise binary classification models in the
manuscript. No real patient data is included or derivable from this file --
all values are randomly sampled.

Run this once to (re)create data/synthetic_demo_data.csv, which is the input
expected by train_and_evaluate_model.py.
"""

import numpy as np
import pandas as pd

RANDOM_SEED = 123
N_SAMPLES = 300  # smaller than the real cohort (n=1,230) for a fast demo


def generate_synthetic_data(n_samples=N_SAMPLES, random_state=RANDOM_SEED):
    rng = np.random.default_rng(random_state)

    # -------------------------------------------------------------------
    # Cluster / endotype label (3 classes, roughly matching real class
    # proportions: cluster 0 largest, cluster 2 smallest)
    # -------------------------------------------------------------------
    cluster = rng.choice([0, 1, 2], size=n_samples, p=[0.45, 0.35, 0.20])

    # Small class-dependent signal is injected into a handful of features
    # below, so the demo models have *something* real to learn -- this is
    # NOT meant to reproduce the manuscript's actual effect sizes.
    cluster_shift = {0: 0.0, 1: 0.4, 2: -0.4}
    shift = np.array([cluster_shift[c] for c in cluster])

    data = {"cluster": cluster}

    # -------------------------------------------------------------------
    # Clinical / demographic features
    # -------------------------------------------------------------------
    data["AGE"] = rng.normal(45, 12, n_samples).clip(18, 85).round(1)
    data["SEX_BIRTH"] = rng.choice(["Male", "Female"], size=n_samples, p=[0.75, 0.25])
    data["ETHNICITY"] = rng.choice(
        ["White", "Black", "Asian", "Other"], size=n_samples, p=[0.6, 0.2, 0.1, 0.1]
    )
    data["HIV_DURATION"] = rng.exponential(8, n_samples).clip(0, 35).round(1)
    data["CD4_NADIR"] = (rng.normal(300, 150, n_samples) + shift * 40).clip(1, 900).round(0)
    data["VL_ZENITH"] = (10 ** rng.normal(4.8, 0.8, n_samples)).clip(50, 1e7).round(0)
    data["CD4_LATEST"] = (rng.normal(650, 220, n_samples) - shift * 60).clip(50, 1800).round(0)
    data["EARLY_CART"] = rng.choice(["Yes", "No"], size=n_samples, p=[0.3, 0.7])
    data["CD4CD8_LATEST"] = rng.normal(1.0, 0.4, n_samples).clip(0.1, 3.0).round(2)
    data["CART_DURATION"] = rng.exponential(6, n_samples).clip(0, 25).round(1)
    data["CART_INTERRUPTED"] = rng.choice(["Yes", "No"], size=n_samples, p=[0.15, 0.85])
    data["CMV_IgG_Serology"] = rng.choice(["Yes", "No"], size=n_samples, p=[0.7, 0.3])
    data["CMV_IgG_IU/mL"] = rng.exponential(150, n_samples).clip(0, 2000).round(1)
    data["RISK_BEHAV"] = rng.choice(["MSM", "Heterosexual", "PWID", "Other"], size=n_samples)
    data["CONTROLLER_1B"] = rng.choice(["Non-EC", "EC"], size=n_samples, p=[0.95, 0.05])

    # Drug exposure (cumulative duration, years)
    for drug in ["-no ART-", "3TC", "ABC", "DTG", "EFV", "FTC", "RTV", "TAF", "TDF"]:
        col = f"{drug}_exposure"
        data[col] = rng.exponential(3, n_samples).clip(0, 20).round(1)

    # -------------------------------------------------------------------
    # Viral reservoir features (main total/intact reservoir markers)
    # -------------------------------------------------------------------
    data["Total.million.avg"] = (rng.normal(0, 1, n_samples) + shift * 1.5).round(3)
    data["intactDT.DSI"] = (rng.normal(0, 1, n_samples) - shift * 1.2).round(3)

    # -------------------------------------------------------------------
    # Bulk transcriptomics (synthetic ENSG-style feature names)
    # -------------------------------------------------------------------
    transcript_genes = [
        "ENSG00000117643:MAN1C1", "ENSG00000271447:MMP28", "ENSG00000081059:TCF7",
        "ENSG00000154027:AK5", "ENSG00000182983:ZNF662", "ENSG00000113319:RASGRF2",
        "ENSG00000092096:SLC22A17", "ENSG00000112394:SLC16A10", "ENSG00000135960:EDAR",
    ]
    for i, gene in enumerate(transcript_genes):
        s = shift if i % 3 == 0 else 0  # only some genes carry class signal
        data[gene] = (rng.normal(0, 1, n_samples) + 0.5 * s).round(3)

    # -------------------------------------------------------------------
    # DNA methylation (synthetic cg-site M-values)
    # -------------------------------------------------------------------
    data["cg13452062"] = rng.normal(0, 1, n_samples).round(3)

    # -------------------------------------------------------------------
    # Flow cytometry (Panel2 immune subsets, absolute counts, log-scale-ish)
    # -------------------------------------------------------------------
    flow_markers = [
        "Panel2.2408.CD4negCD8pos.CCR5pos", "Panel2.2338.CD8pos.Tc2.HLAnegDRpos",
        "Panel2.2162.CD8pos.Tem.PD1pos", "Panel2.2100.CD8posCXCR4posCCR5pos",
        "Panel2.2074.CD8pos.Tc17",
    ]
    for marker in flow_markers:
        data[marker] = rng.normal(0, 1, n_samples).round(3)

    data["CD4_LATEST"]  # already defined above; kept for clarity
    data["SLAMF7.Inflammation"] = rng.normal(0, 1, n_samples).round(3)

    # -------------------------------------------------------------------
    # Ex vivo stimulated cytokine production (naming pattern matches the
    # real script's feature-type detection logic: contains IFNy/IL17/IL22/
    # IL10/IL5/MIP1a/TNF/IL1b)
    # -------------------------------------------------------------------
    cytokine_stimuli = [
        "IL1b.LPS", "IL8.IL1a", "TNF.Spneu", "TNF.CMV", "IL1b.HIVENV",
        "MCP1.CMV", "S.pneu.IFNy", "IL1b.CMV", "E.coli.IL22", "MCP1.PolyIC",
        "PHA.IL5", "C.alb.con.IFNy",
    ]
    for cyt in cytokine_stimuli:
        data[cyt] = rng.lognormal(mean=1.5, sigma=0.8, size=n_samples).round(2)

    df = pd.DataFrame(data)

    # Inject a small amount of missingness (~5%) into a subset of numeric
    # columns, matching the real pipeline's need for imputation.
    numeric_cols = df.select_dtypes(include=[np.number]).columns.drop("cluster")
    missing_cols = rng.choice(numeric_cols, size=max(1, len(numeric_cols) // 4), replace=False)
    for col in missing_cols:
        mask = rng.random(n_samples) < 0.05
        df.loc[mask, col] = np.nan

    df.index.name = "ID"
    df.index = [f"SIM{idx:04d}" for idx in range(1, n_samples + 1)]

    return df


if __name__ == "__main__":
    df = generate_synthetic_data()
    out_path = "demo/data/synthetic_demo_data.csv"
    df.to_csv(out_path)
    print(f"Saved synthetic demo dataset: {out_path}")
    print(f"Shape: {df.shape}")
    print(f"Cluster distribution:\n{df['cluster'].value_counts()}")
