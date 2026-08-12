#!/usr/bin/env Rscript
# =============================================================================
# Script: ex_vivo_overlap_and_enrichment.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Discovery-vs-validation overlap/concordance analysis of the
#              ex-vivo cytokine regression results (all 3 pairwise endotype
#              comparisons x 2 stimulation durations = 6 panels): 4-quadrant
#              consistency plots (FDR- and nominal-significance bases),
#              concordant-direction plots, and a summary table of markers
#              validated in both cohorts.
#
# NOTE: despite the "Enrichment" in this script's original filename, no
# pathway/GO enrichment analysis is actually performed here — this is an
# overlap/concordance analysis of the regression results themselves.
#
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Requires: the 6 discovery/validation CSV pairs written by ex_vivo_analysis.R
# (run once for each of the 3 comparisons, before running this script).
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(writexl)
library(rlang)
library(purrr)
library(openxlsx)

# Set paths — update these to match your local directory structure
# (must match the paths used in ex_vivo_analysis.R)
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"


# -----------------------------------------------------------------------------
# 1. Load and Merge Discovery/Validation Results (All 3 Comparisons)
# -----------------------------------------------------------------------------
# Each comparison's normalized-estimate results (from ex_vivo_analysis.R)
# are loaded for both cohorts and both timepoints, then merged by variable
# so discovery (.x) and validation (.y) columns sit side by side.

comparisons <- c("AllHigh_vs_AllLow", "AllLow_vs_Mixed", "Mixed_vs_AllHigh")
timepoints <- c("24h", "7d")

merged_results <- list()
for (comparison in comparisons) {
  for (timepoint in timepoints) {
    discovery <- read_csv(file.path(output_dir, paste0(
      "ex_vivo_cytokines/", timepoint, "_ex_vivo_results_normalized_estimates_", comparison,
      "_Season_Age_Timetolab_Sex_PC1_Covidvacc_discovery.csv"
    )))
    validation <- read_csv(file.path(output_dir, paste0(
      "ex_vivo_cytokines/", timepoint, "_ex_vivo_results_normalized_estimates_", comparison,
      "_Season_Age_Timetolab_Sex_PC1_Covidvacc_validation.csv"
    )))
    merged_results[[paste0(comparison, "_", timepoint)]] <- left_join(discovery, validation, by = "variable")
  }
}


# -----------------------------------------------------------------------------
# 2. 4-Quadrant Consistency Plots
# -----------------------------------------------------------------------------
# Compares discovery and validation normalized effect estimates for every
# cytokine-stimulus pair. Two significance bases are plotted: FDR-corrected
# discovery vs. nominal validation, and nominal p-value in both cohorts.
# Also produces a simpler "concordant direction" plot (same sign in both
# cohorts, regardless of significance) for each panel.

# Categorizes each point by discovery/validation significance and effect
# direction, using either the FDR-corrected discovery threshold or the
# nominal p-value threshold in both cohorts.
plot_quadrant <- function(df, comparison, time_point, sig_basis = c("fdr", "pvalue"),
                           x_col = "estimate_normalized.x", y_col = "estimate_normalized.y") {
  sig_basis <- match.arg(sig_basis)

  if (sig_basis == "fdr") {
    df <- df %>% mutate(Groups = case_when(
      fdr.x < 0.05 & p_value.y < 0.05 & ((!!sym(x_col) < 0 & !!sym(y_col) < 0) | (!!sym(x_col) > 0 & !!sym(y_col) > 0)) ~
        "Validated direction and significance",
      fdr.x > 0.05 & p_value.y < 0.05 ~ "Significant in validation (P<0.05)",
      fdr.x < 0.05 & p_value.y > 0.05 ~ "Significant in discovery (FDR<0.05)",
      TRUE ~ "Other"
    ))
    group_levels <- c("Validated direction and significance", "Significant in discovery (FDR<0.05)",
                       "Significant in validation (P<0.05)", "Other")
    subtitle <- "Highlighted: FDR < 0.05 in discovery and P-value < 0.05 in validation."
    file_suffix <- "FDRdisc_Pvalval"
  } else {
    df <- df %>% mutate(Groups = case_when(
      p_value.x < 0.05 & ((!!sym(x_col) < 0 & !!sym(y_col) < 0) | (!!sym(x_col) > 0 & !!sym(y_col) > 0)) ~
        "Validated direction and\n significance in discovery (P<0.05)",
      p_value.x > 0.05 & p_value.y < 0.05 ~ "Significant in validation (P<0.05)",
      p_value.x < 0.05 & p_value.y > 0.05 ~ "Significant in discovery (P<0.05)",
      TRUE ~ "Other"
    ))
    group_levels <- c("Validated direction and\n significance in discovery (P<0.05)", "Significant in discovery (P<0.05)",
                       "Significant in validation (P<0.05)", "Other")
    subtitle <- "Highlighted: P-value < 0.05 in discovery and concordant effect size with validation."
    file_suffix <- "Pvaldisc_Pvalval"
  }

  observed_levels <- unique(df$Groups)
  color_map <- setNames(c("#A93226", "#6aa84f", "#5DADE2", "#BEBEBE33"), group_levels)
  df$Groups <- factor(df$Groups, levels = group_levels)
  df_labeled <- df[df$Groups == group_levels[1], ]

  p <- ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = Groups, alpha = 0.5), size = 3) +
    scale_color_manual(values = color_map[observed_levels], breaks = observed_levels) +
    theme_classic(base_size = 10) +
    theme(legend.position = "top") +
    xlab("Normalized Estimate (Discovery)") +
    ylab("Normalized Estimate (Validation)") +
    geom_text_repel(data = df_labeled, aes(label = variable, hjust = 0.5), size = 4, box.padding = unit(0.4, "lines"), max.overlaps = 15) +
    geom_vline(xintercept = 0, linetype = "longdash", colour = "black", size = 0.4) +
    geom_hline(yintercept = 0, linetype = "longdash", colour = "black", size = 0.4) +
    guides(size = FALSE, alpha = FALSE, fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = paste0("Effect Size (Normalized Estimate) in Discovery vs. Validation - ", comparison, " (", time_point, ")"),
         subtitle = subtitle,
         caption = "Corrected for: Age, Sex, Season, Time to lab, genetic PC1, and Covid Vacc.")

  list(plot = p, file_suffix = file_suffix)
}

# Simpler plot: highlights points where discovery and validation agree in
# direction (sign), regardless of statistical significance.
plot_concordant <- function(df, comparison, time_point, x_col = "estimate_normalized.x", y_col = "estimate_normalized.y") {
  df <- df %>% mutate(
    Groups = ifelse((!!sym(x_col) < 0 & !!sym(y_col) < 0) | (!!sym(x_col) > 0 & !!sym(y_col) > 0), "Validated direction", "Other")
  )
  group_levels <- c("Validated direction", "Other")
  observed_levels <- unique(df$Groups)
  color_map <- setNames(c("#A93226", "#BEBEBE33"), group_levels)
  df$Groups <- factor(df$Groups, levels = group_levels)
  df_labeled <- df[df$Groups == "Validated direction", ]

  ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_point(aes(color = Groups, alpha = 0.5), size = 3) +
    scale_color_manual(values = color_map[observed_levels], breaks = observed_levels) +
    theme_classic(base_size = 10) +
    theme(legend.position = "top") +
    xlab("Normalized Estimate (Discovery)") +
    ylab("Normalized Estimate (Validation)") +
    geom_text_repel(data = df_labeled, aes(label = variable, hjust = 0.5), size = 4, box.padding = unit(0.4, "lines"), max.overlaps = 15) +
    geom_vline(xintercept = 0, linetype = "longdash", colour = "black", size = 0.4) +
    geom_hline(yintercept = 0, linetype = "longdash", colour = "black", size = 0.4) +
    guides(size = FALSE, alpha = FALSE, fill = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(title = paste0("Effect Size (Normalized Estimate) in Discovery vs. Validation - ", comparison, " (", time_point, ")"),
         subtitle = "Highlighted: Concordant effect size in discovery and validation.",
         caption = "Corrected for: Age, Sex, Season, Time to lab, genetic PC1, and Covid Vacc.")
}

# --- Generate and save every panel, for both significance bases ---
plot_labels <- list(
  "AllHigh_vs_AllLow_7d"  = "AllHigh vs AllLow", "AllHigh_vs_AllLow_24h" = "AllHigh vs AllLow",
  "AllLow_vs_Mixed_7d"    = "AllLow vs Mixed",   "AllLow_vs_Mixed_24h"   = "AllLow vs Mixed",
  "Mixed_vs_AllHigh_7d"   = "Mixed vs AllHigh",  "Mixed_vs_AllHigh_24h"  = "Mixed vs AllHigh"
)

for (panel_key in names(merged_results)) {
  df <- merged_results[[panel_key]]
  timepoint <- ifelse(grepl("24h", panel_key), "24h", "7d")
  label <- plot_labels[[panel_key]]

  for (basis in c("fdr", "pvalue")) {
    result <- plot_quadrant(df, label, timepoint, sig_basis = basis)
    print(result$plot)
    ggsave(file.path(output_dir, paste0("ex_vivo_cytokines/4QP_", panel_key, "_", result$file_suffix, ".pdf")), plot = result$plot)
  }

  p_concordant <- plot_concordant(df, label, timepoint)
  print(p_concordant)
  ggsave(file.path(output_dir, paste0("ex_vivo_cytokines/4QP_", panel_key, "_Concordant.pdf")), plot = p_concordant)
}


# -----------------------------------------------------------------------------
# 3. Validated Markers Summary
# -----------------------------------------------------------------------------
# One combined table of every cytokine-stimulus pair that was nominally
# significant in discovery (P<0.05) with a concordant effect direction in
# validation, across all 6 panels.

extract_validated <- function(df, comparison_label, timepoint) {
  if (nrow(df) == 0) {
    return(tibble(variable = character(), p_value_discovery = numeric(), p_value_validation = numeric(),
                  Discovery_Estimate = numeric(), Validation_Estimate = numeric(),
                  Comparison = character(), Stimulation_Time = character()))
  }
  df %>%
    mutate(Groups = case_when(
      p_value.x < 0.05 & ((estimate_normalized.x < 0 & estimate_normalized.y < 0) | (estimate_normalized.x > 0 & estimate_normalized.y > 0)) ~
        "Validated direction and significance",
      TRUE ~ "Other"
    )) %>%
    filter(Groups == "Validated direction and significance") %>%
    transmute(variable = variable, p_value_discovery = p_value.x, p_value_validation = p_value.y,
              Discovery_Estimate = estimate_normalized.x, Validation_Estimate = estimate_normalized.y,
              Comparison = comparison_label, Stimulation_Time = timepoint)
}

all_validated <- names(merged_results) %>%
  map(function(panel_key) {
    timepoint <- ifelse(grepl("24h", panel_key), "24h", "7d")
    extract_validated(merged_results[[panel_key]], plot_labels[[panel_key]], timepoint)
  }) %>%
  purrr::discard(~ nrow(.) == 0) %>%
  bind_rows()

write_xlsx(all_validated, path = file.path(output_dir, "ex_vivo_cytokines/validated_markers_summary_discovery_and_validation.xlsx"))


# -----------------------------------------------------------------------------
# 4. Export Source Data for Reviewers
# -----------------------------------------------------------------------------
# Builds source data for all 6 panels in one pass (unlike the flow cytometry
# scripts, all 6 comparison/timepoint combinations already coexist here as
# distinct merged_results entries, so there's no risk of a shared-variable
# collision and no need for a separately-sourced helper).

build_exvivo_export_df <- function(df) {
  classified <- df %>% mutate(Groups = case_when(
    fdr.x < 0.05 & p_value.y < 0.05 & ((estimate_normalized.x < 0 & estimate_normalized.y < 0) | (estimate_normalized.x > 0 & estimate_normalized.y > 0)) ~
      "Validated direction and significance",
    fdr.x > 0.05 & p_value.y < 0.05 ~ "Significant in validation (P<0.05)",
    fdr.x < 0.05 & p_value.y > 0.05 ~ "Significant in discovery (FDR<0.05)",
    TRUE ~ "Other"
  ))
  classified$Groups <- factor(classified$Groups, levels = c(
    "Validated direction and significance", "Significant in discovery (FDR<0.05)",
    "Significant in validation (P<0.05)", "Other"
  ))

  export_df <- classified[, c("variable", "estimate_normalized.x", "p_value.x", "fdr.x",
                               "estimate_normalized.y", "p_value.y", "Groups")]
  colnames(export_df) <- c("Cytokine_Stimulus", "Beta_Discovery", "Pvalue_Discovery", "FDR_Discovery",
                            "Beta_Validation", "Pvalue_Validation", "Significance_Category")
  export_df
}

export_exvivo_source_data <- function(output_path) {
  wb <- createWorkbook()
  for (panel_key in names(merged_results)) {
    export_df <- build_exvivo_export_df(merged_results[[panel_key]])
    addWorksheet(wb, panel_key)
    writeData(wb, panel_key, export_df)
    message(paste0("Added sheet '", panel_key, "' (", nrow(export_df), " rows)."))
  }
  saveWorkbook(wb, output_path, overwrite = TRUE)
  message(paste0("Saved: ", output_path))
}

export_exvivo_source_data(file.path(output_dir, "source_data/SourceData_ExVivo_4QuadrantPlots.xlsx"))
