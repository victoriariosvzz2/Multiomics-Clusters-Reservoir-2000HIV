#!/usr/bin/env Rscript
# =============================================================================
# Script: flow_abs_analysis.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Flow cytometry immunophenotyping analysis (absolute cell
#              counts) comparing two MoCluster endotypes at a time. For each
#              of the complete cohort, discovery cohort, and validation
#              cohort separately: linear regression of each cell population
#              against cluster membership (adjusted for season, age,
#              time-to-lab, COVID vaccination, sex, genetic PC1), summary
#              heatmaps, discovery-vs-validation 4-quadrant consistency
#              plots, and confirmatory violin plots of top markers.
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Run once per pairwise endotype comparison (3 total): set `group_a` and
# `group_b` below. `group_b` is always the reference level, so a positive
# estimate means "higher in group_a relative to group_b" — keep this
# consistent with the filenames/labels below when changing comparisons.
#   - AllHigh vs AllLow  (group_a = "All High", group_b = "All Low")
#   - AllLow  vs Mixed   (group_a = "All Low",  group_b = "Mixed")
#   - Mixed   vs AllHigh (group_a = "Mixed",    group_b = "All High")
#
# Continues in: flow_per_analysis.R (same comparisons, percentage data)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(lubridate)
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
library(writexl)
library(ggh4x)
library(ggrepel)
library(tidyr)
library(ggpubr)

# Set paths — update these to match your local directory structure
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

# Which two endotypes to compare — group_b is the reference level
group_a <- "All High"
group_b <- "All Low"
comparison_label <- paste0(gsub(" ", "", group_a), "_vs_", gsub(" ", "", group_b))


# -----------------------------------------------------------------------------
# 1. Load Data
# -----------------------------------------------------------------------------

sample_table <- read_excel(file.path(data_dir, "clinical/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.xlsx")) %>%
  as.data.frame()

sample_table$NIJMEGEN_DATE_COLLECTION_to_JAN01 <- as.Date(sample_table$DATE_VISIT) - as.Date("2019-01-01")
sample_table$NIJMEGEN_SEASON_SCORE_SINE <- sin((2 * pi * as.numeric(sample_table$NIJMEGEN_DATE_COLLECTION_to_JAN01)) / 365.2425)
sample_table$NIJMEGEN_SEASON_SCORE_COSINE <- cos((2 * pi * as.numeric(sample_table$NIJMEGEN_DATE_COLLECTION_to_JAN01)) / 365.2425)

# --- MoCluster cluster assignments (both cohorts, then restricted to the
# two endotypes being compared here) ---
clustering_results_discovery <- readRDS(file.path(output_dir, "clustering/consensus_clusters_discovery.rds"))
clustering_results_validation <- readRDS(file.path(output_dir, "clustering/consensus_clusters_validation.rds"))

clustering_results_MoCluster <- rbind(clustering_results_discovery$clust.res, clustering_results_validation$clust.res)
clustering_results_MoCluster <- clustering_results_MoCluster %>%
  mutate(ID = as.character(samID), MoCluster_Cluster = factor(clust))

clustering_results_MoCluster <- clustering_results_MoCluster %>%
  mutate(MoCluster_Cluster = recode(MoCluster_Cluster,
                            "1" = "Mixed",
                            "2" = "All Low",
                            "3" = "All High"))

# Restrict to the two endotypes being compared; group_b is the reference
# level, so downstream estimates are "group_a relative to group_b"
clustering_results_MoCluster <- clustering_results_MoCluster[clustering_results_MoCluster$MoCluster_Cluster %in% c(group_a, group_b), ]
clustering_results_MoCluster$MoCluster_Cluster <- relevel(factor(clustering_results_MoCluster$MoCluster_Cluster), ref = group_b)
clustering_results_MoCluster$MoCluster_Cluster <- factor(clustering_results_MoCluster$MoCluster_Cluster, levels = c(group_b, group_a))
table(clustering_results_MoCluster$MoCluster_Cluster)

sample_table <- right_join(sample_table, clustering_results_MoCluster, by = "ID")
sample_table$condition <- sample_table$MoCluster_Cluster

sample_table$CENTER_RUMC <- ifelse(sample_table$CENTER == "RUMC", 1, 0)
sample_table$Black_Ethnicity <- ifelse(sample_table$ETHNICITY == "Black", 1, 0)
sample_table$AGE <- as.numeric(sample_table$AGE)

# --- Genetic PCs (loaded per cohort, then recombined) ---
pcs_discovery <- read.table(file.path(data_dir, "genetics/2000HIV_genetic_pcs_discovery.pca.eigenvec"))
pcs_discovery <- pcs_discovery[, 2:7]
colnames(pcs_discovery) <- c("ID", "PC1_gen", "PC2_gen", "PC3_gen", "PC4_gen", "PC5_gen")

pcs_validation <- read.table(file.path(data_dir, "genetics/2000HIV_genetic_pcs_validation.pca.eigenvec"))
pcs_validation <- pcs_validation[, 2:7]
colnames(pcs_validation) <- c("ID", "PC1_gen", "PC2_gen", "PC3_gen", "PC4_gen", "PC5_gen")

sample_table_discovery <- subset(sample_table, COHORT == "DISCOVERY")
sample_table_validation <- subset(sample_table, COHORT == "VALIDATION")
sample_table_discovery <- merge(sample_table_discovery, pcs_discovery, by = "ID", all.x = TRUE)
sample_table_validation <- merge(sample_table_validation, pcs_validation, by = "ID", all.x = TRUE)
sample_table <- rbind(sample_table_discovery, sample_table_validation)

master_metadata <- sample_table

# --- Flow cytometry data (absolute counts) ---
df <- read_excel(file.path(data_dir, "omics/2000HIV_Flow_Cytometry/2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_356vars_Nov062023.xlsx"))
names(df)[names(df) == "SampleID"] <- "ID"

# --- Fix acquisition dates for 13 participants whose FCM run date differed
# from their eCRF visit date (see 230710_FCM_AcquisitionDates_CalculationFileREADME.xlsx) ---
acq_date <- read.csv2(file.path(data_dir, "clinical/230710_FCM_AcquisitionDates.csv"))
acq_date$RunDATE_FCM[acq_date$RunDATE_FCM == "#N/B"] <- NA
acq_date <- na.omit(acq_date)
acq_date$RunDATE_FCM <- dmy(acq_date$RunDATE_FCM)
master_metadata <- merge(master_metadata, acq_date, by = "ID")

# TIMETOLAB adjustment for these 13 participants: 6 set to <=24h (per known
# collection circumstances), 7 left unchanged after review
master_metadata <- master_metadata %>%
  mutate(
    RunDATE_FCM = ymd(RunDATE_FCM),
    TIMETOLAB_cat = case_when(TIMETOLAB <= 24 ~ "\u226424 hours", TIMETOLAB > 24 ~ ">24 hours"),
    TIMETOLAB_cat = ifelse(
      ID %in% c("RAD005", "RAD021", "RAD027", "RAD032", "RAD039", "EMC096"),
      "\u226424 hours",
      TIMETOLAB_cat
    )
  )
master_metadata$TIMETOLAB <- ifelse(
  master_metadata$ID %in% c("RAD005", "RAD021", "RAD027", "RAD032", "RAD039", "EMC096"),
  20.75,
  master_metadata$TIMETOLAB
)

master_metadata <- master_metadata[, c(
  "ID", "condition", "RunDATE_FCM", "AGE", "SEX_BIRTH", "COHORT", "BMI_BASELINE",
  "SMOKING_CURRENT", "PC1_gen", "Black_Ethnicity", "PANDEMIC_BEFOREAFTER",
  "CENTER_RUMC", "COVID19", "COVID_VACC", "TIMETOLAB", "TIMETOLAB_cat", "CART_DURATION"
)]

# Seasonality covariate (sine/cosine of FCM acquisition date)
master_metadata <- master_metadata %>%
  mutate(
    RunDATE_FCM = ymd(RunDATE_FCM),
    Days_since_Jan2020 = as.numeric(difftime(RunDATE_FCM, as.Date("2020-01-01"), units = "days")),
    SEASON_SIN = sin((2 * pi * Days_since_Jan2020) / 365.2425),
    SEASON_COS = cos((2 * pi * Days_since_Jan2020) / 365.2425)
  )

fcm_meta <- merge(master_metadata, df, by = "ID")


# -----------------------------------------------------------------------------
# 2. Data Transformation
# -----------------------------------------------------------------------------
# Inverse rank-based (normal quantile) transformation of every flow
# cytometry variable, to obtain approximately normally distributed inputs
# for linear regression.

table(fcm_meta$condition, useNA = "always")
fcm_meta$condition <- as.numeric(fcm_meta$condition) - 1
table(fcm_meta$condition, useNA = "always")

start_index <- grep("Panel1_1012_Eosinophils", colnames(fcm_meta))
end_index <- grep("Panel3_3116_Naive.B.cells_CD307d+", colnames(fcm_meta))

fcm_meta_transformed <- fcm_meta
for (i in start_index:end_index) {
  fcm_meta_transformed[, i] <- qnorm((rank(fcm_meta[, i], na.last = "keep") - 0.5) / sum(!is.na(fcm_meta[, i])))
}


# -----------------------------------------------------------------------------
# 3. Linear Regression and Heatmaps by Cohort Scope
# -----------------------------------------------------------------------------
# Run identically on three sample sets: the complete cohort, discovery only,
# and validation only. Each pass fits one linear model per flow marker
# (marker ~ condition + covariates), FDR-corrects across markers, then
# produces a full heatmap of all markers (faceted by broad cell-type
# category) and a "highlights" heatmap restricted to a curated marker list.

# Cell-type grouping used to facet the heatmaps below
classify_cell_type <- function(variable) {
  case_when(
    str_detect(variable, "Neutrophils|Eosinophils|Basophils") ~ "Granulocytes",
    str_detect(variable, "Mono|CD14") ~ "Monocytes",
    str_detect(variable, "NK|CD56\\+CD3\\+") ~ "NK cells",
    str_detect(variable, "DC") ~ "Dendritic cells",
    str_detect(variable, "Th|CD4\\+|Treg|Tfh") ~ "CD4+ T cells",
    str_detect(variable, "CD8\\+") ~ "CD8+ T cells",
    str_detect(variable, "TCR") ~ "TCR cells",
    str_detect(variable, c("Bc|B.cells|SMBC|AM|IM|RM|TLM|CD81\\+|CD19\\+|NBC|Plasma.Cells")) ~ "B cells",
    str_detect(variable, "Tc_CD86+|Tc_PDL1+|Tcells|Tconv") ~ "T cells"
  )
}
cell_type_levels <- c("Granulocytes", "Monocytes", "NK cells", "Dendritic cells",
                       "T cells", "CD4+ T cells", "CD8+ T cells", "TCR cells", "B cells")
heatmap_facet_design <- matrix(c(1, 2, 2, 3, 3, 4, 9, 9, 9,
                                  5, 8, 8, 7, 7, 7, 7, 7, 7,
                                  6, 6, 6, 6, 6, 6, 6, 6, 6),
                                9, 3, byrow = FALSE)

# Curated marker subset shown in the "highlights" heatmap and used to
# restrict the 4-quadrant "broad" panel later
selected_populations <- c(
  "Neutrophils", "Basophils", "Eosinophils", "Monocytes",
  "CM_CD14++CD16-", "NCM_CD14+CD16+", "IM_CD14++CD16+",
  "DC", "mDC", "mDC_CD16-CD1c+", "NKcells", "NK_CD56-CD16+",
  "NK_CD56++CD16-", "NK_CD56++CD16+", "NK_CD56+CD16-",
  "NK_CD56+CD16+", "CD56+CD3+", "CD19+", "Immature.B.cells",
  "Naive.B.cells", "Plasma.Cells", "Switched.Memory.B.cells",
  "Unswitched.Memory.B.cells", "Tcells", "Tconv", "TCRyd",
  "TCRvd1", "TCRvd2", "CD4-CD8+", "CD8+.Tc", "CD8+.Tc1",
  "CD8+.Tc2", "CD8+.Tc17", "CD8+.Tc1/17", "CD8+.Tcm",
  "CD8+.Tem", "CD8+.Temra", "CD8+.Tnaive", "CD4+CD8-",
  "CD4+CD8+", "Treg", "mTreg", "nTreg", "Tfh", "Tfh_1",
  "Tfh_2", "Tfh_1/17", "Tfh_17", "CD4+.Tcm", "CD4+.Tem",
  "CD4+.Temra", "CD4+.Tnaive", "CD4+.Th1", "CD4+.Th2",
  "CD4+.Th1/17", "CD4+.Th17"
)

# Fits marker ~ condition + covariates for every flow cytometry column in
# `data`, returning one row per marker with the condition effect estimate,
# p-value, and BH-adjusted FDR.
run_marker_regressions <- function(data) {
  column_names <- colnames(data)
  start_idx <- grep("Panel1_1012_Eosinophils", column_names)
  end_idx <- grep("Panel3_3116_Naive.B.cells_CD307d+", column_names)

  results <- data.frame(variable = character(), estimate = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
  for (i in start_idx:end_idx) {
    model <- lm(data[, i] ~ condition + SEASON_COS + SEASON_SIN + AGE + TIMETOLAB + COVID_VACC + SEX_BIRTH + PC1_gen, data = data)
    results <- rbind(results, data.frame(
      variable = column_names[i],
      estimate = summary(model)$coefficients[[2, "Estimate"]],
      p_value  = summary(model)$coefficients[[2, "Pr(>|t|)"]],
      stringsAsFactors = FALSE
    ))
  }
  results$fdr <- p.adjust(results$p_value, method = "BH")
  results
}

# Builds the full (all-markers) and highlights (curated subset) heatmaps
# for one set of regression results, and saves both to disk.
plot_regression_heatmaps <- function(results, scope_label) {
  sig <- results %>%
    mutate(
      stars = case_when(
        fdr < 0.05 & estimate > 0 ~ "*0.05 (pos)",
        fdr < 0.05 & estimate < 0 ~ "*0.05 (neg)",
        TRUE ~ NA_character_
      ),
      stars = case_when(
        fdr < 0.001 & estimate > 0 ~ "**0.001 (pos)",
        fdr < 0.001 & estimate < 0 ~ "**0.001 (neg)",
        TRUE ~ stars
      ),
      Cell_type = classify_cell_type(variable),
      variable_edited = sub("^(?:[^_]*_){2}(.*)", "\\1", variable)
    ) %>%
    mutate(
      variable = factor(variable, levels = unique(variable)),
      Cell_type = factor(Cell_type, levels = cell_type_levels)
    )

  # --- Full heatmap (all markers) ---
  p_full <- ggplot(sig, aes(y = variable_edited, x = paste0(group_a, " (compared to ", group_b, ")"), fill = estimate)) +
    geom_tile(colour = "black", size = 0.25) +
    scale_fill_gradient2(high = "#B05000", mid = "white", low = "#00598A", midpoint = 0, limits = c(-2, 2)) +
    geom_text(aes(label = stars), color = "black", size = 2) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw(base_size = 7) +
    theme(
      plot.caption = element_text(hjust = 0, face = "italic"),
      legend.key.size = unit(.4, 'cm'),
      panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    facet_manual(vars(Cell_type), scales = "free_y", design = heatmap_facet_design)
  print(p_full)

  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(output_dir, paste0("flow_cytometry/FCM_ABS_", comparison_label, "_", scope_label, "_AgeSexSeasonTimetolabPC1Covidvacc.", ext)),
           plot = p_full, height = 35, width = 20, units = "cm", dpi = 300)
  }

  # --- Highlights heatmap (curated marker subset) ---
  sig_highlights <- sig %>% filter(variable_edited %in% selected_populations)
  write.xlsx(sig_highlights, file.path(output_dir, paste0("flow_cytometry/FCM_ABS_Results_", comparison_label, "_", scope_label, "_AgeSexSeasonTimetolabPC1Covidvacc.xlsx")))

  p_highlights <- ggplot(sig_highlights, aes(y = variable_edited, x = group_a, fill = estimate)) +
    geom_tile(colour = "black", size = 0.25) +
    scale_fill_gradient2(high = "darkred", mid = "white", low = "turquoise4", midpoint = 0, limits = c(-2, 2)) +
    geom_text(aes(label = stars), color = "black", size = 4) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw(base_size = 7) +
    labs(x = " ", y = " ", caption = " ") +
    theme(
      plot.caption = element_text(hjust = 0, face = "italic"),
      legend.key.size = unit(.4, 'cm'),
      panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7)
    ) +
    facet_manual(vars(Cell_type), scales = "free_y", design = heatmap_facet_design)
  print(p_highlights)

  for (ext in c("tiff", "pdf")) {
    ggsave(file.path(output_dir, paste0("flow_cytometry/FCM_ABS_Highlights_", comparison_label, "_", scope_label, "_AgeSexSeasonTimetolabPC1Covidvacc.", ext)),
           plot = p_highlights, height = 12, width = 20, units = "cm", dpi = 300)
  }

  sig
}

# --- Run for all three cohort scopes ---
complete_cohort <- fcm_meta_transformed
results_complete <- run_marker_regressions(complete_cohort)
write_xlsx(results_complete, file.path(output_dir, paste0("flow_cytometry/FCM_ABS_", comparison_label, "_complete_AgeSexSeasonTimetolabPC1Covidvacc.xlsx")))
results_complete <- plot_regression_heatmaps(results_complete, "complete")

discovery_cohort <- subset(fcm_meta_transformed, COHORT == "DISCOVERY")
results_discovery <- run_marker_regressions(discovery_cohort)
write_xlsx(results_discovery, file.path(output_dir, paste0("flow_cytometry/FCM_ABS_", comparison_label, "_discovery_AgeSexSeasonTimetolabPC1Covidvacc.xlsx")))
results_discovery <- plot_regression_heatmaps(results_discovery, "discovery")

validation_cohort <- subset(fcm_meta_transformed, COHORT == "VALIDATION")
results_validation <- run_marker_regressions(validation_cohort)
write_xlsx(results_validation, file.path(output_dir, paste0("flow_cytometry/FCM_ABS_", comparison_label, "_validation_AgeSexSeasonTimetolabPC1Covidvacc.xlsx")))
results_validation <- plot_regression_heatmaps(results_validation, "validation")


# -----------------------------------------------------------------------------
# 4. 4-Quadrant Plots (Discovery vs. Validation Consistency)
# -----------------------------------------------------------------------------
# Compares discovery and validation effect estimates for every marker,
# classifying each as validated (significant + same direction in both),
# significant in one cohort only, or not significant. Run for both
# significance bases (FDR-corrected and nominal p-value) and both marker
# sets (all markers, and the curated "broad" subset from Section 3).

plot_4qp <- function(sig_basis = c("fdr", "pvalue"), marker_set = c("all", "broad")) {
  sig_basis <- match.arg(sig_basis)
  marker_set <- match.arg(marker_set)
  threshold <- 0.05

  merged_data <- results_discovery %>%
    inner_join(results_validation, by = "variable", suffix = c("_disc", "_val")) %>%
    mutate(variable_edited = sub("^[^_]+_[^_]+_", "", variable))

  if (sig_basis == "fdr") {
    merged_data <- merged_data %>%
      mutate(Groups = case_when(
        fdr_disc < threshold & p_value_val < threshold &
          ((estimate_disc < 0 & estimate_val < 0) | (estimate_disc > 0 & estimate_val > 0)) ~ "Validated direction and significance",
        fdr_disc > threshold & p_value_val < threshold ~ "Significant in validation (P-value<0.05)",
        fdr_disc < threshold & p_value_val > threshold ~ "Significant in discovery (FDR<0.05)",
        TRUE ~ "Other"
      ))
    group_levels <- c("Validated direction and significance", "Significant in discovery (FDR<0.05)",
                       "Significant in validation (P-value<0.05)", "Other")
    title_suffix <- ""
    file_suffix <- "FDRdisc_Pvalval"
  } else {
    merged_data <- merged_data %>%
      mutate(Groups = case_when(
        p_value_disc < threshold & p_value_val < threshold &
          ((estimate_disc < 0 & estimate_val < 0) | (estimate_disc > 0 & estimate_val > 0)) ~ "Validated direction and significance (nominal)",
        p_value_disc > threshold & p_value_val < threshold ~ "Significant in validation (P-value<0.05)",
        p_value_disc < threshold & p_value_val > threshold ~ "Significant in discovery (P-value<0.05)",
        TRUE ~ "Other"
      ))
    group_levels <- c("Validated direction and significance (nominal)", "Significant in discovery (P-value<0.05)",
                       "Significant in validation (P-value<0.05)", "Other")
    title_suffix <- " (Nominal Significance)"
    file_suffix <- "Pvaldisc_Pvalval"
  }

  if (marker_set == "broad") {
    merged_data <- merged_data %>% filter(variable_edited %in% selected_populations)
    file_suffix <- paste0(file_suffix, "_Broad")
  } else {
    file_suffix <- paste0(file_suffix, "_All")
  }

  merged_data$Groups <- factor(merged_data$Groups, levels = group_levels)
  color_palette <- setNames(c("#A93226", "#6aa84f", "#5DADE2", "#BEBEBE33"), group_levels)
  options(ggrepel.max.overlaps = 40)

  p <- ggplot(merged_data, aes(x = estimate_disc, y = estimate_val, color = Groups)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_color_manual(values = color_palette) +
    theme_classic(base_size = 10) +
    xlab("Effect Estimate (Discovery)") +
    ylab("Effect Estimate (Validation)") +
    geom_text_repel(
      data = merged_data %>% filter(Groups %in% group_levels[1:2]),
      aes(label = variable_edited), size = 3, box.padding = unit(0.4, "lines"), max.overlaps = 18
    ) +
    geom_vline(xintercept = 0, linetype = "longdash", color = "black", size = 0.4) +
    geom_hline(yintercept = 0, linetype = "longdash", color = "black", size = 0.4) +
    guides(size = FALSE, alpha = FALSE, fill = guide_legend(nrow = 4, byrow = TRUE)) +
    labs(title = paste0("Effect Estimates in Discovery and Validation", title_suffix))
  print(p)

  filename_base <- paste0("flow_cytometry/FCM_ABS_4QP_", comparison_label, "_discovery_and_validation_AgeSexSeasonTimetolabPC1Covidvacc_", file_suffix)
  for (ext in c("tiff", "pdf")) {
    ggsave(file.path(output_dir, paste0(filename_base, ".", ext)), plot = p, height = 7, width = 8, dpi = 300)
  }

  merged_data
}

plot_4qp(sig_basis = "fdr",    marker_set = "all")
plot_4qp(sig_basis = "pvalue", marker_set = "all")
plot_4qp(sig_basis = "fdr",    marker_set = "broad")
plot_4qp(sig_basis = "pvalue", marker_set = "broad")


# -----------------------------------------------------------------------------
# 5. Violin Plots (Confirmatory)
# -----------------------------------------------------------------------------
# Raw and covariate-adjusted (residualized) values for a hand-picked set of
# markers, to visually confirm the regression results above are not
# artifacts of the transformation or model fit.

violin_markers <- c(
  "Panel2_2147_CD4+.Th1_PD1+", "Panel2_2314_CD4+.Tnaive_HLA-DR+",
  "Panel2_2321_CD4+.Th1_HLA-DR+", "Panel3_3040_Naive.B.cells",
  "Panel3_3098_Naive.B.cells_CD81+", "Panel3_3068_NBC_CD21+",
  "Panel2_2082_CD8+PD1+", "Panel2_2078_CD4+PD1+",
  "Panel2_2524_CD4-CD8+_CCR7+", "Panel2_2523_CD4+CD8-_CCR7+",
  "Panel2_2525_CD4+CD8+_CCR7+"
)

complete_cohort <- complete_cohort %>%
  mutate(condition_label = recode(condition, "0" = group_b, "1" = group_a))

violin_data <- lapply(violin_markers, function(marker) {
  model <- lm(
    complete_cohort[[marker]] ~ AGE + SEX_BIRTH + COVID_VACC + PC1_gen + TIMETOLAB + SEASON_SIN + SEASON_COS,
    data = complete_cohort, na.action = na.exclude
  )
  bind_rows(
    data.frame(ID = complete_cohort$ID, condition = complete_cohort$condition, marker = marker, value = complete_cohort[[marker]], type = "Raw"),
    data.frame(ID = complete_cohort$ID, condition = complete_cohort$condition, marker = marker, value = resid(model), type = "Adjusted")
  )
}) %>% bind_rows() %>%
  mutate(condition_label = recode(condition, "0" = group_b, "1" = group_a))

p_violin <- ggplot(violin_data, aes(x = factor(condition_label), y = value, fill = factor(condition_label))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  stat_compare_means() +
  facet_grid(type ~ marker, scales = "free_y") +
  labs(x = "Condition", y = "Marker value") +
  theme_minimal()
print(p_violin)


# -----------------------------------------------------------------------------
# 6. Export Source Data for Reviewers
# -----------------------------------------------------------------------------
# Appends this comparison's discovery-vs-validation 4-quadrant source data
# to the shared cross-modality workbook (see source_data/export_4quadrant_source_data.R).

source(file.path(project_dir, "scripts/source_data/export_4quadrant_source_data.R"))

shared_output_path <- file.path(output_dir, "source_data/SourceData_4QuadrantPlots_FCM_ABS.xlsx")

# Raw variable names (not the shortened display names) corresponding to the
# curated "highlights" marker subset — same set regardless of cohort scope,
# since the mapping from raw name to display name doesn't depend on cohort.
highlighted_variable_names <- results_discovery %>%
  mutate(variable_edited = sub("^(?:[^_]*_){2}(.*)", "\\1", variable)) %>%
  filter(variable_edited %in% selected_populations) %>%
  pull(variable) %>%
  unique()

export_4qp_panel(results_discovery, results_validation,
                  paste0(comparison_label, "_ABS"), shared_output_path,
                  selected_populations = highlighted_variable_names)
