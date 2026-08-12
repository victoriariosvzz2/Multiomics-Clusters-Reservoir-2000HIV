#!/usr/bin/env Rscript
# =============================================================================
# Script: ex_vivo_analysis.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Ex-vivo cytokine production analysis comparing two MoCluster
#              endotypes at a time, for both stimulation durations (24h and
#              7-day). Rank-based regression (Rfit) of each cytokine-
#              stimulus pair against cluster membership (adjusted for
#              season, age, time-to-lab, sex, genetic PC1, COVID
#              vaccination), run separately for the discovery and
#              validation cohorts, followed by heatmap visualization (raw
#              and mean-normalized effect estimates).
# Original analysis by: Elise Meeder (elise.meeder@radboudumc.nl)
# Adapted by: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Run once per pairwise endotype comparison (3 total): set `group_a` and
# `group_b` below. `group_b` is always the reference level, so a positive
# estimate means "higher in group_a relative to group_b".
#   - AllHigh vs AllLow  (group_a = "All High", group_b = "All Low")
#   - AllLow  vs Mixed   (group_a = "All Low",  group_b = "Mixed")
#   - Mixed   vs AllHigh (group_a = "Mixed",    group_b = "All High")
#
# Unlike the flow cytometry scripts, this uses rank-based (Rfit) regression
# rather than OLS, since cytokine production data is highly right-skewed;
# "Bent1" scores are used per Rfit's recommendation for skewed data.
#
# For the discovery cohort, effect significance uses BH-adjusted FDR; for
# validation, the nominal p-value is used directly (stored in the same
# `fdr` column for downstream consistency, matching the original analysis).
#
# Continues in: ex_vivo_overlap_and_enrichment.R (all 3 comparisons combined)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)
library(limma)
library(Rfit)
library(reshape2)
library(RNOmni)
library(magrittr)
library(purrr)
library(ggplot2)

# Set paths — update these to match your local directory structure
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

# Which two endotypes to compare — group_b is the reference level
group_a <- "Mixed"
group_b <- "All Low"
comparison_label <- paste0(gsub(" ", "", group_a), "_vs_", gsub(" ", "", group_b))


# -----------------------------------------------------------------------------
# 1. Timepoint Configuration
# -----------------------------------------------------------------------------
# Cytokine/stimulus panels, column boundaries, and display orderings differ
# between the 24h and 7-day ex-vivo stimulation panels.

timepoint_config <- list(
  "24h" = list(
    start_marker = "IL1b_CMV",
    end_marker   = "MIP1a_IL1a",
    cytokines_extract_order = c("IL1b", "IL1RA", "TNF", "IL6", "IL8", "IL10", "MCP1", "MIP1a"),
    cytokines_matrix_order  = c("IL1RA", "IL10", "IL6", "IL8", "IL1b", "TNF", "MCP1", "MIP1a"),
    stimulus_patterns = c("*CMV*" = "CMV", "*Spneu*" = "Spneu", "*LPS*" = "LPS", "*IMQ*" = "IMQ",
                           "*HIVENV*" = "HIVENV", "*PolyIC*" = "PolyIC", "*IL1a*" = "IL1a"),
    heatmap_size = list(width = 7, height = 5),
    missing_cytstims = c('IL10_PolyIC', 'IL1b_PolyIC', 'IL6_PolyIC', 'IL8_PolyIC', 'TNF_PolyIC',
                          'IL10_IMQ', 'TNF_IMQ', 'IL10_IL1a', 'IL1b_IL1a', 'TNF_IL1a',
                          'IL10_HIVENV', 'IL10_CMV')
  ),
  "7d" = list(
    start_marker = "S.aureus_IFNy",
    end_marker   = "C.alb.hy_IL22",
    cytokines_extract_order = c("IL22", "IL5", "IL10", "IFNy", "IL17"),
    cytokines_matrix_order  = c("IL22", "IL5", "IL10", "IFNy", "IL17"),
    stimulus_patterns = c("*PHA*" = "PHA", "*MTB*" = "MTB", "*E.coli*" = "E.coli", "*S.pneu*" = "S.pneu",
                           "*S.aureus*" = "S.aureus", "*C.alb.con*" = "C.alb.con", "*C.alb.hy*" = "C.alb.hy"),
    heatmap_size = list(width = 5, height = 5),
    missing_cytstims = c('S.aureus_IL5', 'E.coli_IL5')
  )
)


# -----------------------------------------------------------------------------
# 2. Helper Functions
# -----------------------------------------------------------------------------

# Fits cytokine ~ condition + covariates (rank-based regression) for every
# cytokine-stimulus column between start_marker and end_marker.
run_cytokine_regressions <- function(data, config) {
  column_names <- colnames(data)
  start_idx <- grep(config$start_marker, column_names)
  end_idx   <- grep(config$end_marker, column_names)

  results <- data.frame(variable = character(), estimate = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
  for (i in start_idx:end_idx) {
    fit <- rfit(data[, i] ~ CU + season_cos + season_sin + AGE + TIMETOLAB + SEX_BIRTH + PC1_gen + COVID_VACC,
                scores = bentscores1, data = data)
    results <- rbind(results, data.frame(
      variable = column_names[i],
      estimate = summary(fit)$coefficients[[2, "Estimate"]],
      p_value  = summary(fit)$coefficients[[2, "p.value"]],
      stringsAsFactors = FALSE
    ))
  }
  results$sd <- apply(data[, start_idx:end_idx], 2, sd, na.rm = TRUE)
  results$cohen_d <- results$estimate / results$sd
  results
}

# Reshapes long-format results (one row per cytokine-stimulus variable) into
# a wide cytokine x stimulus matrix, cleaning stimulus names via the
# timepoint's regex pattern list.
build_wide_matrix <- function(long_df, value_col, config) {
  per_cytokine <- lapply(config$cytokines_matrix_order, function(cyt) {
    sub_df <- long_df[str_detect(long_df$variable, cyt), c("variable", value_col)]
    for (pattern in names(config$stimulus_patterns)) {
      sub_df$variable[grepl(pattern, sub_df$variable, ignore.case = TRUE)] <- config$stimulus_patterns[[pattern]]
    }
    colnames(sub_df)[colnames(sub_df) == value_col] <- cyt
    sub_df
  })
  purrr::reduce(per_cytokine, full_join, by = "variable")
}

# Raster heatmap of effect estimates with FDR/p-value significance asterisks
plot_fdr_heatmap <- function(estimate_matrix, fdr_matrix, out_path, width, height) {
  long_e <- melt(estimate_matrix); colnames(long_e) <- c("Stimulus", "Cytokine", "effect")
  long_p <- melt(fdr_matrix); colnames(long_p) <- c("Stimulus", "Cytokine", "Pvalue")
  long_p <- long_p %>% mutate(Asterisks = ifelse(Pvalue <= 0.0005, "***", ifelse(Pvalue <= 0.005, "**", ifelse(Pvalue <= 0.05, "*", NA))))

  tiff(file = out_path, width = width, height = height, units = "in", res = 600)
  print(
    ggplot() +
      geom_raster(data = long_e, aes(x = Cytokine, y = Stimulus, fill = effect)) +
      scale_y_discrete(limits = rev(levels(long_e$Stimulus))) +
      geom_text(data = long_p, aes(x = Cytokine, y = Stimulus, label = Asterisks), size = 10) +
      labs(title = " ", x = " ", y = " ", caption = " ") +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                            limits = c(-0.3, 0.3), na.value = "lightgray",
                            breaks = c(-0.3, 0, 0.3), guide = guide_colorbar(ticks = FALSE)) +
      theme(axis.text = element_text(size = 13, color = "black"),
            panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  )
  dev.off()
}

# Formats a long-results table with cytokine/stimulus columns and a
# significance-label column, filling in NA rows for cytokine-stimulus pairs
# that weren't measured at this timepoint (see timepoint_config).
step1_format_sumstat <- function(sumstat_table, config, use_fdr_label = TRUE) {
  missing_df <- data.frame(variable = config$missing_cytstims, estimate = NA, p_value = NA, sd = NA, cohen_d = NA, fdr = NA)
  sumstat_table %<>% rbind(missing_df)
  # 24h variables are "cytokine_stimulus"; 7d variables are "stimulus_cytokine"
  if (identical(config$start_marker, "IL1b_CMV")) {
    sumstat_table <- sumstat_table %>% mutate(cytokine = str_split_i(variable, "_", 1), stimulus = str_split_i(variable, "_", 2))
  } else {
    sumstat_table <- sumstat_table %>% mutate(cytokine = str_split_i(variable, "_", 2), stimulus = str_split_i(variable, "_", 1))
  }
  sig_col <- if (use_fdr_label) "fdr" else "p_value"
  sumstat_table %>% mutate(label = ifelse(.data[[sig_col]] < 0.05, "*", NA))
}

# Normalizes effect estimates by each variable's mean expression in the
# given sample subset (so estimates are comparable in magnitude across
# cytokines with very different absolute expression levels).
step2_normalize_estimate <- function(sumstat_table, cytdat, filter_list) {
  cytdat_filtered <- cytdat %>% filter(ID %in% filter_list) %>% tibble::column_to_rownames("ID")
  mean_cyt <- colMeans(cytdat_filtered, na.rm = TRUE) %>%
    as.data.frame() %>% cbind(colnames(cytdat_filtered)) %>% set_colnames(c("mean_expression", "variable"))
  sumstat_table %<>% left_join(mean_cyt) %>% mutate(estimate_normalized = estimate / mean_expression)
  sumstat_table
}

# Cytokine x stimulus heatmap of normalized effect estimates
step3_plot_heatmap <- function(sumstat, title, caption, timepoint) {
  if (timepoint == "24h") {
    sumstat %<>% mutate(cytokine = factor(cytokine, levels = c("TNF", "IL6", "IL1b", "IL1RA", "IL8", "MCP1", "MIP1a", "IL10"))) %>%
      mutate(stimulus = factor(stimulus, levels = c("PolyIC", "IL1a", "CMV", "HIVENV", "IMQ", "LPS", "Spneu")))
  } else if (timepoint == "7d") {
    sumstat %<>% mutate(cytokine = factor(cytokine, levels = c("IFNy", "IL5", "IL17", "IL22", "IL10"))) %>%
      mutate(stimulus = factor(stimulus, levels = c("S.pneu", "S.aureus", "MTB", "E.coli", "C.alb.hy", "C.alb.con", "PHA")))
  }
  plot <- ggplot(data = sumstat[!is.na(sumstat$estimate_normalized), ],
                 aes(x = cytokine, y = stimulus, fill = estimate_normalized)) +
    geom_tile() +
    geom_text(aes(label = label), size = 6) +
    scale_fill_gradient2(low = "#00598A", mid = "white", high = "#B05000", na.value = "grey", limits = c(-0.6, 0.6)) +
    theme_classic() + coord_fixed() +
    labs(title = title, caption = caption, fill = "Normalized\nestimate") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 11.5, vjust = 0.5), axis.text.y = element_text(size = 11.5))

  if (timepoint == "24h") {
    plot + scale_x_discrete(labels = c(expression(paste("TNF-", alpha)), "IL-6", expression(paste("IL-1", beta)),
                                        "IL-1RA", "IL-8", "MCP-1", expression(paste("MIP-1", alpha)), "IL-10")) +
      scale_y_discrete(labels = c("Poly(I:C)", expression(paste("IL-1", alpha)), "CMV", "HIVENV", "IMQ", "LPS", expression(italic("S. pneumoniae"))))
  } else {
    plot + scale_x_discrete(labels = c(expression(paste("IFN-", gamma)), "IL-5", "IL-17", "IL-22", "IL-10")) +
      scale_y_discrete(labels = c(expression(italic("S. pneumoniae")), expression(italic("S. aureus")), expression(italic("M. tuberculosis")),
                                   expression(italic("E. coli")), expression(italic("C. albicans hyphae")), expression(italic("C. albicans conidia")), "PHA"))
  }
}


# -----------------------------------------------------------------------------
# 3. Load Data
# -----------------------------------------------------------------------------

sample_table <- readRDS(file.path(data_dir, "clinical/2000HIV_bulk_transcriptomics_sample_table.RDS"))

clustering_results_discovery <- readRDS(file.path(output_dir, "clustering/consensus_clusters_discovery.rds"))
clustering_results_validation <- readRDS(file.path(output_dir, "clustering/consensus_clusters_validation.rds"))
clustering_results_MoCluster <- rbind(clustering_results_discovery$clust.res, clustering_results_validation$clust.res) %>%
  mutate(DONOR_ID = as.character(samID), MoCluster_Cluster = factor(clust)) %>%
  mutate(MoCluster_Cluster = recode(MoCluster_Cluster, "1" = "Mixed", "2" = "All Low", "3" = "All High"))

sample_table <- right_join(sample_table, clustering_results_MoCluster, by = "DONOR_ID")
sample_table$condition <- sample_table$MoCluster_Cluster

# Restrict to the two endotypes being compared; encode as numeric 0/1 with
# group_b = 0 (reference), group_a = 1
sample_table$condition <- factor(sample_table$condition, levels = c(group_b, group_a))
sample_table <- sample_table[sample_table$condition %in% c(group_a, group_b), ]
table(sample_table$condition)
sample_table$condition <- as.numeric(sample_table$condition) - 1
table(sample_table$condition)

sample_table$Black_Ethnicity <- ifelse(sample_table$ETHNICITY == "Black", 1, 0)
sample_table$CENTER_RUMC <- ifelse(sample_table$CENTER == "RUMC", 1, 0)
sample_table$AGE <- as.numeric(sample_table$AGE)
sample_table$ID <- sample_table$DONOR_ID

# --- Genetic PCs (loaded per cohort, then recombined) ---
pcs_discovery <- read.table(file.path(data_dir, "genetics/2000HIV_genetic_pcs_discovery.pca.eigenvec"))
pcs_discovery <- pcs_discovery[, 2:7]
colnames(pcs_discovery) <- c("ID", "PC1_gen", "PC2_gen", "PC3_gen", "PC4_gen", "PC5_gen")

pcs_validation <- read.table(file.path(data_dir, "genetics/2000HIV_genetic_pcs_validation.pca.eigenvec"))
pcs_validation <- pcs_validation[, 2:7]
colnames(pcs_validation) <- c("ID", "PC1_gen", "PC2_gen", "PC3_gen", "PC4_gen", "PC5_gen")

CA_discovery <- merge(subset(sample_table, COHORT == "DISCOVERY"), pcs_discovery, by = "ID", all.x = TRUE)
CA_validation <- merge(subset(sample_table, COHORT == "VALIDATION"), pcs_validation, by = "ID", all.x = TRUE)
CA <- rbind(CA_discovery, CA_validation)

# Seasonality covariate (sine/cosine of visit date, days since Jan 2019)
CA$dayssincejan2019 <- as.numeric(as.POSIXct(CA$DATE_VISIT) - as.POSIXct("2019-01-01"))
CA$season_sin <- sin((2 * pi * CA$dayssincejan2019) / 365.2425)
CA$season_cos <- cos((2 * pi * CA$dayssincejan2019) / 365.2425)

CA <- CA[, colnames(CA) %in% c("ID", "condition", "AGE", "SEX_BIRTH", "season_sin", "season_cos",
                                "Black_Ethnicity", "COVID_VACC", "COHORT", "CENTER_RUMC", "TIMETOLAB", "PC1_gen")]
CA <- CA %>% mutate(CU = case_when(condition == 1 ~ 1, is.na(condition) ~ NA_real_, TRUE ~ 0))

# --- Ex-vivo cytokine data (both stimulation panels) ---
cyt24 <- read_excel(file.path(data_dir, "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_RAW.xlsx"))
cyt7  <- read_excel(file.path(data_dir, "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_7DAY_38meas_1754samples_afterQC_RAW.xlsx"))


# -----------------------------------------------------------------------------
# 4. Run Regression, Build Matrices, and Plot Heatmaps (per Timepoint x Cohort)
# -----------------------------------------------------------------------------

for (timepoint in c("24h", "7d")) {
  config <- timepoint_config[[timepoint]]
  cyt_data <- if (timepoint == "24h") cyt24 else cyt7

  clincyt <- merge(CA, cyt_data, by = "ID")
  clincyt <- clincyt[!is.na(clincyt$CU), ]
  table(clincyt$CU, useNA = "always")

  clincyt_discovery <- subset(clincyt, COHORT == "DISCOVERY")
  clincyt_validation <- subset(clincyt, COHORT == "VALIDATION")

  for (cohort_name in c("discovery", "validation")) {
    cohort_data <- if (cohort_name == "discovery") clincyt_discovery else clincyt_validation
    results <- run_cytokine_regressions(cohort_data, config)

    # Discovery uses BH-adjusted FDR; validation uses the nominal p-value
    # directly (stored in the same `fdr` column for downstream consistency)
    results$fdr <- if (cohort_name == "discovery") p.adjust(results$p_value, method = "fdr") else results$p_value

    estimate_matrix <- build_wide_matrix(results[, c("variable", "cohen_d")], "cohen_d", config)
    fdr_matrix <- build_wide_matrix(results[, c("variable", "fdr")], "fdr", config)

    write.csv(fdr_matrix, file.path(output_dir, paste0("ex_vivo_cytokines/", timepoint, "_ex_vivo_FDRpvalues_", comparison_label, "_Season_Age_Timetolab_Sex_PC1_Covidvacc_", cohort_name, ".csv")))
    write.csv(estimate_matrix, file.path(output_dir, paste0("ex_vivo_cytokines/", timepoint, "_ex_vivo_estimates_", comparison_label, "_Season_Age_Timetolab_Sex_PC1_Covidvacc_", cohort_name, ".csv")))

    plot_fdr_heatmap(estimate_matrix, fdr_matrix,
                      file.path(output_dir, paste0("ex_vivo_cytokines/", timepoint, "_heatmap_", comparison_label, "_Season_Age_Timetolab_Sex_PC1_Covidvacc_", cohort_name, ".tiff")),
                      width = config$heatmap_size$width, height = config$heatmap_size$height)

    # --- Normalized-estimate heatmap ---
    sumstat <- step1_format_sumstat(results, config, use_fdr_label = (cohort_name == "discovery"))
    group_no_na <- cohort_data %>% filter(!is.na(condition)) %>% pull(ID)
    sumstat <- step2_normalize_estimate(sumstat, cyt_data, group_no_na)
    write.csv(sumstat, file.path(output_dir, paste0("ex_vivo_cytokines/", timepoint, "_ex_vivo_results_normalized_estimates_", comparison_label, "_Season_Age_Timetolab_Sex_PC1_Covidvacc_", cohort_name, ".csv")))

    step3_plot_heatmap(sumstat,
                        title = paste0(timepoint, " cytokines (", tools::toTitleCase(cohort_name), " Cohort) - ", group_a, " (compared to ", group_b, ")"),
                        caption = "Corrected for: Seasonality, Age, Time to lab, Sex, Genetic PC1, Covid Vacc.",
                        timepoint = timepoint)
    ggsave(file.path(output_dir, paste0("ex_vivo_cytokines/", timepoint, "_results_", comparison_label, "_Season_Age_Timetolab_Sex_PC1_Covidvacc_", cohort_name, ".pdf")),
           width = config$heatmap_size$width + 2, height = config$heatmap_size$height + 2)

    # Keep discovery/validation results in the environment for downstream use
    assign(paste0("results_", cohort_name), results)
  }
}
