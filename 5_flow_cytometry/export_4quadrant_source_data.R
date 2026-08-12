#!/usr/bin/env Rscript
# =============================================================================
# Script: export_4quadrant_source_data.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Shared helper, sourced by each flow_abs_analysis.R /
#              flow_per_analysis.R run, that builds the FDR-based
#              discovery-vs-validation significance classification
#              (replicating the exact logic used for the 4-quadrant plots)
#              and appends it as one sheet to a single shared source-data
#              Excel workbook.
#
# WHY THIS APPROACH (instead of exporting `merged_data` directly): each
# analysis script builds the same-named `merged_data` object twice — once
# for the FDR-based panel, once for the nominal p-value panel — so capturing
# it from the environment risks grabbing whichever one ran last. This
# helper instead rebuilds the FDR-based classification itself from
# results_discovery/results_validation, so it's independent of what has run
# before it.
#
# HOW TO USE: call export_4qp_panel(...) from each of the six analysis
# scripts (3 pairwise comparisons x 2 measurement types) any time after
# results_discovery/results_validation exist. Order between scripts doesn't
# matter — each call appends its own sheet to the same shared workbook,
# creating the file on the first call.
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
# =============================================================================

library(openxlsx)
library(dplyr)

# Replicates the join + Groups classification used for the FDR-based
# 4-quadrant plot: "Validated" (significant + same direction in both
# cohorts), "Significant in discovery only", "Significant in validation
# only", or "Other".
build_fdr_classified_data <- function(results_discovery, results_validation, fdr_threshold = 0.05) {
  merged <- results_discovery %>%
    inner_join(results_validation, by = "variable", suffix = c("_disc", "_val")) %>%
    mutate(
      Groups = case_when(
        fdr_disc < fdr_threshold & p_value_val < fdr_threshold &
          ((estimate_disc < 0 & estimate_val < 0) | (estimate_disc > 0 & estimate_val > 0)) ~
          "Validated direction and significance",
        fdr_disc > fdr_threshold & p_value_val < fdr_threshold ~
          "Significant in validation (P-value<0.05)",
        fdr_disc < fdr_threshold & p_value_val > fdr_threshold ~
          "Significant in discovery (FDR<0.05)",
        TRUE ~ "Other"
      ),
      # Shortened display name (drops the "PanelN_NNNN_" prefix)
      Cell_Population = sub("^[^_]+_[^_]+_", "", variable)
    )

  merged$Groups <- factor(
    merged$Groups,
    levels = c("Validated direction and significance",
               "Significant in discovery (FDR<0.05)",
               "Significant in validation (P-value<0.05)",
               "Other")
  )
  merged
}

# Builds one panel's source data (optionally restricted to a curated marker
# subset) and appends it as a new sheet to the shared workbook at
# `output_path`, creating the file if it doesn't exist yet.
export_4qp_panel <- function(results_discovery, results_validation, sheet_name, output_path,
                              selected_populations = NULL, fdr_threshold = 0.05) {

  classified <- build_fdr_classified_data(results_discovery, results_validation, fdr_threshold = fdr_threshold)

  if (!is.null(selected_populations)) {
    classified <- classified %>% filter(variable %in% selected_populations)
  }

  cols_to_keep <- c("Cell_Population", "estimate_disc", "p_value_disc", "fdr_disc",
                     "estimate_val", "p_value_val")
  if ("fdr_val" %in% colnames(classified)) cols_to_keep <- c(cols_to_keep, "fdr_val")
  cols_to_keep <- c(cols_to_keep, "Groups")

  missing_cols <- setdiff(cols_to_keep, colnames(classified))
  if (length(missing_cols) > 0) {
    stop(paste0("export_4qp_panel: classified data is missing expected column(s): ",
                paste(missing_cols, collapse = ", "),
                ". Check that results_discovery/results_validation have columns: variable, estimate, p_value, fdr."))
  }

  export_df <- classified[, cols_to_keep]
  colnames(export_df)[colnames(export_df) == "estimate_disc"] <- "Beta_Discovery"
  colnames(export_df)[colnames(export_df) == "p_value_disc"]  <- "Pvalue_Discovery"
  colnames(export_df)[colnames(export_df) == "fdr_disc"]      <- "FDR_Discovery"
  colnames(export_df)[colnames(export_df) == "estimate_val"]  <- "Beta_Validation"
  colnames(export_df)[colnames(export_df) == "p_value_val"]   <- "Pvalue_Validation"
  if ("fdr_val" %in% colnames(export_df)) colnames(export_df)[colnames(export_df) == "fdr_val"] <- "FDR_Validation"
  colnames(export_df)[colnames(export_df) == "Groups"] <- "Significance_Category"

  if (nrow(export_df) == 0) {
    warning(paste0("export_4qp_panel: '", sheet_name, "' has 0 rows -- check selected_populations filtering."))
  }

  wb <- if (file.exists(output_path)) loadWorkbook(output_path) else createWorkbook()

  if (sheet_name %in% names(wb)) {
    warning(paste0("Sheet '", sheet_name, "' already exists in ", output_path, " -- overwriting its contents."))
    removeWorksheet(wb, sheet_name)
  }

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, export_df)
  saveWorkbook(wb, output_path, overwrite = TRUE)

  message(paste0("Appended sheet '", sheet_name, "' (", nrow(export_df), " rows) to ", output_path))
}
