# Differential Expression & GSEA

Bulk RNA-seq differential expression analysis between the three MoCluster
immune-reservoir endotypes (AllHigh, AllLow, Mixed), followed by Gene Set
Enrichment Analysis (GSEA) on the resulting differential expression results.

## Contents

| File | Description |
|---|---|
| `deg_analysis.R` | DESeq2 differential expression pipeline: data loading, EDA/PCA diagnostics, LIMMA batch-effect correction, DE testing, and result visualization (volcano, MA, upset, Venn, hierarchical clustering plots). |
| `gsea_analysis.R` | GSEA on the AllHigh vs. Mixed comparison, across Reactome, KEGG, GO Biological Process, and MSigDB Hallmark gene sets. |
| `DESeq_functions_KD_RK_VR.R` | Custom helper functions sourced by `deg_analysis.R` (DESeq2 wrappers, plotting, batch correction). Not run directly. |

## Run order

1. **`clustering/`** must be run first — both scripts here depend on its output
   (`clustering/consensus_clusters_discovery.rds` and `..._validation.rds`).
2. **`deg_analysis.R`** — run **twice**, once per cohort. Open the script and
   set the `cohort` variable near the top:
   ```r
   cohort <- "discovery"   # then re-run with "validation"
   ```
   This is a manual edit, not a command-line argument — there is no way to
   pass the cohort in from outside the script.
3. **`gsea_analysis.R`** — run **once**, after `deg_analysis.R` has been run
   with `cohort <- "discovery"` (GSEA was performed on the discovery cohort
   only, as it is the largest). It loads the DE results that
   `deg_analysis.R` saves for the discovery cohort.

## Configuration

Before running either script, update the path placeholders near the top:

```r
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"
```

`DESeq_functions_KD_RK_VR.R` must be placed at
`project_dir/scripts/DESeq_functions_KD_RK_VR.R` — `deg_analysis.R` sources
it from that exact relative path.

## Outputs

Both scripts write to `output_dir/differential_expression/` (plots, tables,
RData/RDS objects, Excel workbooks); GSEA outputs additionally go to
`output_dir/differential_expression/gsea/`. Filenames are tagged with the
`cohort` value used for that run, so discovery and validation outputs do not
overwrite each other.
