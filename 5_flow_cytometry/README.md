# Flow Cytometry (Immunophenotyping)

Linear regression of flow cytometry cell populations against MoCluster
endotype membership, run separately for absolute cell counts (ABS) and
percentage of parent population (PER), for each pairwise endotype
comparison.

## Contents

| File | Description |
|---|---|
| `flow_abs_analysis.R` | Regression, heatmaps, 4-quadrant consistency plots, and confirmatory violin plots on absolute cell counts. |
| `flow_per_analysis.R` | Same analysis on percentage-of-parent data (no violin plots / broad-subset 4-quadrant panel — not part of the original PER analysis). |
| `export_4quadrant_source_data.R` | Shared helper sourced by both scripts above (kept in `scripts/source_data/`, not in this folder — see below). Appends each comparison's source data to one shared workbook per measurement type. |

## Run order

Each script must be run **once per pairwise endotype comparison** (3
comparisons), so **6 runs total** across both scripts:

1. **`clustering/`** must be run first — both scripts depend on its output
   (`clustering/consensus_clusters_discovery.rds` and `..._validation.rds`).
2. Open `flow_abs_analysis.R`, set the comparison near the top, and run:
   ```r
   group_a <- "All High"   # the group of interest
   group_b <- "All Low"    # the reference group — estimates are "group_a vs. group_b"
   ```
   Repeat for the other two comparisons:
   - `group_a = "All Low"`, `group_b = "Mixed"`
   - `group_a = "Mixed"`, `group_b = "All High"`
3. Repeat all three runs for `flow_per_analysis.R`.

Each run of `flow_abs_analysis.R` and `flow_per_analysis.R` sources `export_4quadrant_source_data.R` 
and appends a sheet to a shared workbook, run order between the three comparisons
doesn't matter, but all three should be run before treating that workbook as
complete.

## Configuration

Update the path placeholders near the top of both scripts:

```r
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"
```

`export_4quadrant_source_data.R` must be placed at
`project_dir/scripts/source_data/export_4quadrant_source_data.R` — both
analysis scripts source it from that exact relative path.


## Outputs

Both scripts write to `output_dir/flow_cytometry/` (per-comparison
regression tables, full and "highlights" heatmaps, 4-quadrant plots).
The shared 4-quadrant source data goes to
`output_dir/source_data/SourceData_4QuadrantPlots_FCM_ABS.xlsx` and
`..._FCM_PER.xlsx`, with one sheet per comparison added across the three
runs of each script. Filenames are tagged with the comparison
(`group_a` vs. `group_b`) so the three comparisons' outputs don't
overwrite each other.
