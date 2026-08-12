# Multiomics-Clusters-Reservoir-2000HIV

This repository hosts the analysis scripts used in the manuscript:

**"Multi-Omics Clustering Differentiates the Total and Intact HIV Reservoirs and Related Host Immune Mechanisms"**

submitted to *Nature Communications* (manuscript reference: NCOMMS-26-060652-T).

<!-- TODO: replace with the actual filename once uploaded, e.g. main_figures/graphical_abstract.png -->
![Graphical abstract](main_figures/graphical_abstract.png)

## Preprint

A preprint of this manuscript is available on bioRxiv:

> Rios-Vazquez V, Delporte M, Otten T, Matzaraki V, Vos WAJW, Blaauw MJT, van Eekeren LE, Groenendijk AL, Knoll R, Aschenbrenner AC, Arts RJW, dos Santos JC, Navas A, Gerlo S, van Lunzen J, Joosten LAB, Netea MG, van der Ven AJAM, Vandekerckhove L. Multi-Omics Clustering Differentiates the Total and Intact HIV Reservoirs and Related Host Immune Mechanisms. *bioRxiv*. 2026. doi: [10.64898/2026.06.24.734029](https://doi.org/10.64898/2026.06.24.734029)

## Overview

The analyses in this study were performed using custom R and Python scripts applying previously published and well-established statistical and machine-learning packages (DESeq2, MOVICS, CIMLR, MoCluster, limma, Rfit, clusterProfiler/ReactomePA, scikit-learn, XGBoost, SHAP). No new analytical methods or software tools were developed.

All analytical steps, parameters, and packages are described in detail in the Main Text and Methods sections of the manuscript to ensure reproducibility.

## Repository structure

```
1_qc/                       Quality control and preprocessing for each omics layer
2_clustering/                Multi-omics integrative clustering (MOVICS: SNF, PINSPlus, NEMO,
                              COCA, ConsensusClustering, CIMLR, MoCluster) and algorithm selection
3_machine_learning/           XGBoost/Random Forest/SVM/Logistic Regression/KNN classification
                              of the reservoir-immune endotypes, with SHAP feature importance
    demo/                        Runnable demo on synthetic data (no patient data)
    trained_models_real/        Models trained on the real study cohort
4_differential_expression/   Bulk RNA-seq differential expression (DESeq2) between endotypes,
                              and downstream GSEA (Reactome, KEGG, GO, MSigDB Hallmark)
5_flow_cytometry/            Immunophenotyping analysis (absolute counts and percentages)
                              per pairwise endotype comparison
6_ex_vivo_cytokines/         Ex-vivo cytokine production analysis (24h and 7-day stimulation)
                              and discovery/validation overlap analysis
main_figures/
    figure1/ ... figure6/         Final figure image, source data workbook, and source-data
                                  README for each main manuscript figure — see
                                  main_figures/README.md
LICENSE
```

Each numbered folder contains its own `README.md` with the scripts it contains, the order they need to be run in, required inputs, and a description of what each script outputs. `main_figures/README.md` links each figure to the analysis folder that generated its source data. Folders `3` through `6` all depend on the clustering result from `2_clustering/` but are otherwise independent of one another.

## Requirements

- R (v4.3.2) and Python (v3.9.12), as reported in the manuscript's Methods.
- Required packages are loaded at the top of each script; Python dependencies for `3_machine_learning/` are additionally listed in `3_machine_learning/requirements.txt`.
- File paths at the top of each script (`project_dir`, `data_dir`, `output_dir`) are placeholders — update these to match your local directory structure before running.

## Reproducibility

- `3_machine_learning/demo/` runs end-to-end on a synthetic dataset (no patient data) and can be used to verify the modeling pipeline runs correctly in a new environment.
- `3_machine_learning/trained_models_real/` contains the classifiers as trained on the real study cohort, for inspection without needing to re-run training.
- Scripts elsewhere in the repository require the underlying 2000HIV cohort data (see Data availability below) to run.

## Data availability

This repository contains **scripts only** — no patient-level data files are shared here. The datasets used in this study are available via the Radboud Data Repository (RDR):

| Data | Access |
|---|---|
| Bulk RNA-seq, DNA methylation, ex-vivo cytokine production | https://doi.org/10.34973/p96d-kz55 |
| Plasma proteomics | https://doi.org/10.34973/qk29-f305 |
| HIV-1 DNA reservoir quantification | https://doi.org/10.34973/2bf5-yc03 |
| Flow cytometry (controlled access) | RDR collection ID: `ru.rumc.2000hiv_t0000321a_dsc_002` |

Data access is subject to the 2000HIV consortium's data governance policy and approval by the data access committee. The 2000HIV cohort is registered at ClinicalTrials.gov ([NCT03994835](https://clinicaltrials.gov/study/NCT03994835)).

## License

See [LICENSE](./LICENSE).

## Contact

Author of the repository: Victoria Rios-Vazquez ([Victoria.RiosVazquez@radboudumc.nl](mailto:Victoria.RiosVazquez@radboudumc.nl))

For questions or requests during peer review, please contact:
- Corresponding author: linos.vandekerckhove@ugent.be
- Main author and analysis inquiries: Victoria.RiosVazquez@radboudumc.nl