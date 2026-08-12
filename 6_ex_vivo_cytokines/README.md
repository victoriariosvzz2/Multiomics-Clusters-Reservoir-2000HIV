# Ex-Vivo Cytokines

Ex-vivo cytokine production analysis (24h and 7-day stimulation panels)
comparing MoCluster endotypes, followed by a discovery-vs-validation
overlap/concordance analysis across all three pairwise comparisons.

## Contents

| File | Description |
|---|---|
| `ex_vivo_analysis.R` | Rank-based (Rfit) regression of each cytokine-stimulus pair against cluster membership, for both stimulation durations and both cohorts. Run once per pairwise comparison. |
| `ex_vivo_overlap_and_enrichment.R` | Loads all 3 comparisons' results, produces 4-quadrant discovery-vs-validation consistency plots, a validated-markers summary table, and the reviewer source-data export. Run once, after all 3 comparisons have been analyzed. |


## Run order

1. **`clustering/`** must be run first: `ex_vivo_analysis.R` depends on its
   output (`clustering/consensus_clusters_discovery.rds` and
   `..._validation.rds`).
2. Run `ex_vivo_analysis.R` **three times**, once per pairwise comparison:
   ```r
   group_a <- "All High"   # the group of interest
   group_b <- "All Low"    # the reference group — estimates are "group_a vs. group_b"
   ```
   - `group_a = "All Low"`, `group_b = "Mixed"`
   - `group_a = "Mixed"`, `group_b = "All High"`

   Each run analyzes both the 24h and 7-day panels and both cohorts in one
   pass — no separate timepoint switch needed.
3. Run `ex_vivo_overlap_and_enrichment.R` **once**, after all three
   comparisons above have completed — it reads all 6 comparison/timepoint
   result files written by step 2.

## Configuration

Update the path placeholders near the top of both scripts:

```r
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"
```

## Outputs

`ex_vivo_analysis.R` writes to `output_dir/ex_vivo_cytokines/` (regression
tables, effect/FDR matrices, and heatmaps, both raw and mean-normalized,
per comparison, per timepoint, per cohort).

`ex_vivo_overlap_and_enrichment.R` writes 4-quadrant and concordance plots
to the same folder, the combined validated-markers summary
(`validated_markers_summary_discovery_and_validation.xlsx`), and the shared
reviewer source-data workbook at
`output_dir/source_data/SourceData_ExVivo_4QuadrantPlots.xlsx` (one sheet
per comparison/timepoint panel, 6 total).
