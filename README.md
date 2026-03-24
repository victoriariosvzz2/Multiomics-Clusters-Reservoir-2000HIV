# Multiomics-Clusters-Reservoir-2000HIV

This repository hosts the analysis scripts used in the manuscript:

**"Multi-Omics Clustering Differentiates the Total and Intact HIV Reservoirs and Related Host Immune Mechanisms"**
submitted to *Nature Immunology* (manuscript reference: NI-A42808-T).

## Overview

The analyses in this study were performed using custom R and Python scripts applying previously published and well-established statistical and machine-learning packages. No new analytical methods or software tools were developed.

All analytical steps, parameters, and packages are described in detail in the Main Text and Methods sections of the manuscript to ensure reproducibility.

## Current status

This repository is actively being updated. Scripts are being added progressively, prioritizing those central to the main findings:

- **QC scripts** — currently available
- **Main analytical scripts** — being added during peer review
- **All scripts** — will be fully available upon acceptance of the manuscript

## Repository structure

Scripts will be organized by analysis type as they are added:

- `/QC` — quality control and preprocessing scripts for each omics layer
- `/clustering` — multi-omics clustering pipeline (MOVICS/MoCluster, outlier exclusion, cluster validation)
- `/machine_learning` — XGBoost training, RandomizedSearchCV, SHAP feature importance analyses

## Data availability

This repository contains **scripts only**. No data files are shared here. All datasets used in this study are available via the Radboud Data Repository (RDR) as described in the manuscript's Data Availability section.

## Reproducibility

The code relies on publicly available R (v4.3.2) and Python (v3.9.12) packages cited in the manuscript. Scripts are provided in a structured and documented form to allow reproduction of the main analyses.

## Contact

For questions or requests during peer review, please contact:
- Corresponding author: linos.vandekerckhove@ugent.be
- Main author and Analysis inquiries: Victoria.RiosVazquez@radboudumc.nl
