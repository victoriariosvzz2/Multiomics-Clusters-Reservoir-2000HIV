#!/usr/bin/env Rscript
# =============================================================================
# Script: QC_and_preprocessing.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Load the six omics layers used in the final clustering model
#              (bulk transcriptome, proteome, DNA methylation, flow cytometry
#              immunophenotyping, ex-vivo cytokines, viral reservoir), apply
#              the cohort exclusion list, check sample overlap across omics
#              layers, and perform the per-layer normalization, transformation,
#              and feature selection needed before multi-omics clustering.
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Continues in: multiomics_clustering.R (MOVICS: CIMLR, MoCluster, +9 algos)
#
# Data types among the six layers used:
#   - Count data:       Bulk transcriptome, Flow cytometry (immunophenotyping)
#   - Continuous data:  Proteome, Ex-vivo cytokines, Viral reservoir,
#                        DNA Methylation M-values
#
# Recommended transformations (see COPS paper
# https://doi.org/10.1371/journal.pcbi.1012275 and COPS multi_omics.R
# https://github.com/UEFBiomedicalInformaticsLab/COPS):
#   - Methylation beta values: logit-transform toward normality
#   - Other continuous data:   log-transform
#   - Count/sequencing data:   log2-transform (size factor + VST)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(devtools)
library(MOVICS)  # installed from GitHub (devtools::install_github)
library(dplyr)
library(readr)
library(readxl)
library(openxlsx)
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(biomaRt)

# Suppress warnings (kept from original notebook)
options(warn = -1)

# NOTE on CIMLR installation (macOS gfortran issue):
# CIMLR (github.com/danro9685/CIMLR) can fail to build on macOS unless
# gfortran's LDFLAGS/CPPFLAGS/FLIBS are set and a ~/.R/Makevars file is
# written pointing at libgfortran before running:
#   devtools::install_github("danro9685/CIMLR", ref = "R")

# Set paths — update these to match your local directory structure
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

source(file.path(project_dir, "scripts/Functions.R"))

# -----------------------------------------------------------------------------
# 1. Load Clinical and Reservoir Data
# -----------------------------------------------------------------------------

# Main clinical sample table with cohort assignment (DISCOVERY/VALIDATION),
# demographics, and HIV-related clinical parameters
sample_table <- read_excel(file.path(data_dir, "clinical/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.xlsx"))
cat("Clinical sample table:", dim(sample_table), "\n")

# HIV-1 DNA reservoir quantification data (total + intact, DSI-corrected)
viral_reservoir <- read_csv(file.path(data_dir, "reservoir/Reservoir_Rscript_processed_2025.csv")) %>%
  as.data.frame()
cat("Viral reservoir (raw):", dim(viral_reservoir), "\n")

# -----------------------------------------------------------------------------
# 2. Apply Exclusion Criteria
# -----------------------------------------------------------------------------

# Load exclusion lists
# - "Exclusion": samples failing QC or inclusion criteria
# - "Separate_analysis": samples analyzed separately (e.g. on immunomodulatory drugs)
exclusion_IDs <- read_excel(
  file.path(data_dir, "clinical/20220322_Exclusion_2000HIV.xlsx"),
  sheet = "Exclusion", skip = 1
)
separate_analysis_IDs <- read_excel(
  file.path(data_dir, "clinical/20220322_Exclusion_2000HIV.xlsx"),
  sheet = "Separate_analysis", skip = 1
)

viral_reservoir_clean <- viral_reservoir[
  !viral_reservoir$SampleID %in% c(exclusion_IDs$`Study ID`, separate_analysis_IDs$`Study ID`),
]
cat("Reservoir samples before exclusion:", nrow(viral_reservoir), "\n")
cat("Reservoir samples after exclusion:", nrow(viral_reservoir_clean), "\n")

# Prepare the viral reservoir dataframe (keeps the full set of reservoir
# variables, not just the two used directly in clustering — the remaining
# columns are used later for the clinical characterization of the clusters).
# Column groups:
#   - DNA_conc, Input_volume:  assay input QC metrics
#   - Total_million_avg:       total HIV-1 DNA copies/million CD4+ T cells
#                               (used in clustering)
#   - IPDA_million_DSI, intactDT_DSI: intact HIV-1 DNA copies/million cells
#                               from the Intact Proviral DNA Assay (intactDT_DSI
#                               is used in clustering)
#   - TRIPGAG/TRIPPOL/D4PCR_million_DSI: per-probe PCR readouts (Gag/Pol
#                               regions) that feed into the intactness call
#   - RAINBOW_million:          total reservoir quantification from the
#                               Rainbow digital PCR assay (alternative total
#                               DNA measure)
#   - DECISION TREE:            intactness classification outcome per sample
#   - UD/UR/OK:                 sample quality flag (Undetermined/
#                               Unreportable/OK)
#   - DSI:                      input-adequacy/correction flag applied to
#                               the assay readouts
#   - Total_cells_tested:       number of cells assayed per sample
viral_reservoir_df <- viral_reservoir_clean[!is.na(viral_reservoir_clean$SampleID), ]
rownames(viral_reservoir_df) <- viral_reservoir_df$SampleID
viral_reservoir_df <- viral_reservoir_df[, colnames(viral_reservoir_df) %in% c(
  "DNA_conc",
  "Input_volume",
  "Total_million_avg",
  "IPDA_million_DSI",
  "intactDT_DSI",
  "TRIPGAG_million_DSI",
  "TRIPPOL_million_DSI",
  "D4PCR_million_DSI",
  "RAINBOW_million",
  "DECISION TREE",
  "UD/UR/OK",
  "DSI",
  "Total_cells_tested"
)]
viral_reservoir_df <- viral_reservoir_df %>% mutate(across(everything(), as.numeric))
cat("Viral reservoir matrix dimensions:", dim(viral_reservoir_df), "\n")

# Define discovery and validation cohort IDs
# Sub-cohort assignment was made prior to any analytical preprocessing
discovery_ids <- sample_table[sample_table$COHORT == "DISCOVERY", ]$ID
validation_ids <- sample_table[sample_table$COHORT == "VALIDATION", ]$ID

# -----------------------------------------------------------------------------
# 3. Load Omics Datasets
# -----------------------------------------------------------------------------

# --- Bulk RNA-seq (raw counts) ---
# The raw counts matrix is indexed by LAB_ID (sequencing lab identifier),
# not by the cohort's DONOR_ID/SampleID — so we first load a lookup table
# that maps LAB_ID to DONOR_ID, then use it to re-index the count matrix.
sample_table_bulkrnaseq <- readRDS(file.path(data_dir, "clinical/2000HIV_bulk_transcriptomics_sample_table.RDS"))
sample_table_bulkrnaseq$ID <- sample_table_bulkrnaseq$DONOR_ID

# Load raw counts (genes x samples) and transpose to samples x genes;
# drop the gene-annotation columns (DESCRIPTION, GENETYPE, etc.) first so
# the transpose leaves only numeric count data.
bulk_transcriptomics_df <- readRDS(
  file.path(
    data_dir,
    "omics/2000HIV_bulk_transcriptomics/2000HIV_bulk_transcriptomics_raw_counts.RDS"
  )
)
bulk_transcriptomics_df <- t(bulk_transcriptomics_df[, !names(bulk_transcriptomics_df) %in% c("DESCRIPTION", "GENETYPE", "SYMBOL", "GENEID", "CHR")])

# Re-index from LAB_ID to DONOR_ID via the lookup table, then drop the
# now-redundant ID columns so only gene counts remain.
sample_table_bulkrnaseq$LAB_ID <- as.character(sample_table_bulkrnaseq$LAB_ID)
bulk_transcriptomics_df <- as.data.frame(bulk_transcriptomics_df)
bulk_transcriptomics_df$LAB_ID <- as.character(rownames(bulk_transcriptomics_df))
bulk_transcriptomics_df <- merge(
  bulk_transcriptomics_df,
  sample_table_bulkrnaseq[c("LAB_ID", "DONOR_ID")],
  by = "LAB_ID",
  all.x = TRUE
)
rownames(bulk_transcriptomics_df) <- bulk_transcriptomics_df$DONOR_ID
bulk_transcriptomics_df <- bulk_transcriptomics_df[, !names(bulk_transcriptomics_df) %in% c("LAB_ID", "DONOR_ID")]

# Restrict to samples that passed the reservoir exclusion criteria above
bulk_transcriptomics_df <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% viral_reservoir_clean$SampleID, ]
cat("Bulk RNA-seq:", dim(bulk_transcriptomics_df), "\n")

# --- Proteomics (Olink NPX values, already log2-transformed) ---
# Proteins with value <LOD in >25% of individuals were already filtered
# out during QC. First column is the sample ID, used as rownames.
proteomics_df <- read_excel(file.path(data_dir, "omics/2000HIV_Proteomics/Olink.xlsx")) %>% as.data.frame()
rownames(proteomics_df) <- proteomics_df$SampleID
proteomics_df <- proteomics_df[, -1]
proteomics_df <- proteomics_df[rownames(proteomics_df) %in% viral_reservoir_clean$SampleID, ]
cat("Proteomics:", dim(proteomics_df), "\n")

# --- Flow Cytometry (absolute counts / immunophenotyping, QC-processed Nov 2023) ---
# Final dataset after 3 rounds of QC filtering — 355 immune cell subset
# variables. First column is the sample ID, used as rownames.
flow_abs_df <- read_excel(
  file.path(
    data_dir,
    "omics/2000HIV_Flow_Cytometry/2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_356vars_Nov062023.xlsx"
  )
) %>%
  as.data.frame()
rownames(flow_abs_df) <- flow_abs_df$SampleID
flow_abs_df <- flow_abs_df[, -1]
flow_abs_df <- flow_abs_df[rownames(flow_abs_df) %in% viral_reservoir_clean$SampleID, ]
cat("Flow cytometry (immunophenotyping):", dim(flow_abs_df), "\n")

# --- Ex Vivo Cytokine Data (24h and 7-day stimulations merged) ---
# Two separate assay files (different stimulation durations, different
# measured cytokines) are loaded, each restricted to samples that passed
# exclusion, then merged into one table keyed by sample ID.
ex_vivo_24h_df <- read_excel(
  file.path(
    data_dir,
    "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_RAW.xlsx"
  )
) %>% as.data.frame()
rownames(ex_vivo_24h_df) <- ex_vivo_24h_df$ID
ex_vivo_24h_df <- ex_vivo_24h_df[, -1]
ex_vivo_24h_df <- ex_vivo_24h_df[rownames(ex_vivo_24h_df) %in% viral_reservoir_clean$SampleID, ]

ex_vivo_7d_df <- read_excel(
  file.path(
    data_dir,
    "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_7DAY_38meas_1754samples_afterQC_RAW.xlsx"
  )
) %>% as.data.frame()
rownames(ex_vivo_7d_df) <- ex_vivo_7d_df$ID
ex_vivo_7d_df <- ex_vivo_7d_df[, -1]
ex_vivo_7d_df <- ex_vivo_7d_df[rownames(ex_vivo_7d_df) %in% viral_reservoir_clean$SampleID, ]

# Merge 24h and 7d stimulation data (all = TRUE to preserve all samples)
ex_vivo_df <- merge(ex_vivo_24h_df, ex_vivo_7d_df, by = "row.names", all = TRUE)
rownames(ex_vivo_df) <- ex_vivo_df$Row.names
ex_vivo_df <- ex_vivo_df[, -1]
cat("Ex vivo cytokine (24h + 7d merged):", dim(ex_vivo_df), "\n")

# --- DNA Methylation (M-values from EPIC array) ---
methylation_df <- readRDS(file.path(data_dir, "omics/2000HIV.Mvalue.rds")) %>% t()
methylation_df <- methylation_df[rownames(methylation_df) %in% viral_reservoir_clean$SampleID, ]
cat("DNA Methylation:", dim(methylation_df), "\n")

gc()

# -----------------------------------------------------------------------------
# 4. Sample Overlap Across Omics Layers
# -----------------------------------------------------------------------------
# Before proceeding with preprocessing, check the overlapping samples across
# the six layers to confirm none of them is disproportionately limiting the
# analysis (i.e. contributing far fewer samples than the rest).

samples_list <- list(
  Bulk_transcriptomics = rownames(bulk_transcriptomics_df),
  Proteomics           = rownames(proteomics_df),
  Flow_cytometry       = rownames(flow_abs_df),
  Ex_vivo              = rownames(ex_vivo_df),
  Methylation          = rownames(methylation_df),
  Viral_reservoir      = rownames(viral_reservoir_df[!is.na(viral_reservoir_df$Total_million_avg), ])
)

all_samples <- unique(unlist(samples_list))
# Build a binary sample x layer presence matrix (1 = sample has data in
# that layer, 0 = missing) so overlap can be visualized as a heatmap.
presence_matrix <- sapply(samples_list, function(s) as.integer(all_samples %in% s))
rownames(presence_matrix) <- all_samples

overlapping_samples <- Reduce(intersect, samples_list)
cat("\nSamples per layer:\n")
print(sapply(samples_list, length))
cat("Samples present in ALL six layers:", length(overlapping_samples), "\n")

pheatmap(
  t(presence_matrix),
  color = c("gray90", "steelblue"),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  legend = FALSE,
  main = "Sample Overlap Across Omics Layers"
)

# -----------------------------------------------------------------------------
# 5. Split into Discovery and Validation Cohorts
# -----------------------------------------------------------------------------
# Each layer is split into the two pre-assigned sub-cohorts so that every
# downstream normalization/transformation step (Section 6) can be fit
# separately per cohort, avoiding information leakage from validation
# into discovery.

bulk_disc   <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% discovery_ids, ]
bulk_val    <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% validation_ids, ]
prot_disc   <- proteomics_df[rownames(proteomics_df) %in% discovery_ids, ]
prot_val    <- proteomics_df[rownames(proteomics_df) %in% validation_ids, ]
flow_disc   <- flow_abs_df[rownames(flow_abs_df) %in% discovery_ids, ]
flow_val    <- flow_abs_df[rownames(flow_abs_df) %in% validation_ids, ]
exvivo_disc <- ex_vivo_df[rownames(ex_vivo_df) %in% discovery_ids, ]
exvivo_val  <- ex_vivo_df[rownames(ex_vivo_df) %in% validation_ids, ]
methyl_disc <- methylation_df[rownames(methylation_df) %in% discovery_ids, ]
methyl_val  <- methylation_df[rownames(methylation_df) %in% validation_ids, ]
res_disc    <- viral_reservoir_df[rownames(viral_reservoir_df) %in% discovery_ids, ]
res_val     <- viral_reservoir_df[rownames(viral_reservoir_df) %in% validation_ids, ]

# -----------------------------------------------------------------------------
# 6. Normalization and Transformation
# -----------------------------------------------------------------------------

# --- Bulk RNA-seq: DESeq2 Variance Stabilizing Transformation (VST) ---
# Applied independently per sub-cohort to prevent data leakage.

# Filter to low-expression genes: keep genes with at least 5 counts in at
# least half of all samples (both cohorts combined, before the DESeq2 fit)
sample_table_temp <- sample_table[sample_table$ID %in% c(discovery_ids, validation_ids), ]
bulk_t <- t(bulk_transcriptomics_df)
bulk_t <- bulk_t[, colnames(bulk_t) %in% sample_table$ID]
keep <- rowSums(bulk_t >= 5) >= (ncol(bulk_t) / 2)
bulk_t_filtered <- bulk_t[keep, ]

# Removes ENSEMBL version suffixes (e.g. ENSG00000001.5 -> ENSG00000001)
# so gene IDs match cleanly against the BioMart annotation in Section 7
strip_version <- function(mat) { colnames(mat) <- gsub("\\.\\d+$", "", colnames(mat)); mat }

# Fit DESeq2's size-factor normalization + variance-stabilizing transform
# for each cohort separately (design = ~1: no experimental groups, this is
# purely a normalization step, not a differential expression test)
dds_disc <- DESeqDataSetFromMatrix(
  countData = bulk_t_filtered[, colnames(bulk_t_filtered) %in% discovery_ids],
  colData   = sample_table_temp[sample_table_temp$ID %in% discovery_ids, ],
  design    = ~1
)
dds_disc <- DESeq(dds_disc)
bulk_vst_disc <- assay(varianceStabilizingTransformation(estimateSizeFactors(dds_disc), blind = TRUE)) %>% 
  t() %>% 
  strip_version()

dds_val <- DESeqDataSetFromMatrix(
  countData = bulk_t_filtered[, colnames(bulk_t_filtered) %in% validation_ids],
  colData   = sample_table_temp[sample_table_temp$ID %in% validation_ids, ],
  design    = ~1
)
dds_val <- DESeq(dds_val)
bulk_vst_val <- assay(varianceStabilizingTransformation(estimateSizeFactors(dds_val), blind = TRUE)) %>% 
  t() %>% 
  strip_version()

cat("Bulk RNA-seq (VST) — discovery:", dim(bulk_vst_disc), " validation:", dim(bulk_vst_val), "\n")

# --- Proteomics: no additional transformation ---
# NPX values are already log2-transformed by Olink.
prot_log_disc <- prot_disc
prot_log_val  <- prot_val

# --- Flow cytometry (immunophenotyping): log2(x + 1) ---
# Right-skewed absolute cell counts; log2 stabilizes variance.
flow_log_disc <- log2(flow_disc + 1)
flow_log_val  <- log2(flow_val + 1)

# --- Ex vivo cytokines: log2(x + 1) ---
exvivo_log_disc <- log2(exvivo_disc + 1)
exvivo_log_val  <- log2(exvivo_val+ 1)

# --- DNA Methylation: shift to non-negative ---
# M-values can be negative; add absolute minimum to shift without changing shape.
methyl_disc_nn <- methyl_disc + abs(min(methyl_disc))
methyl_val_nn  <- methyl_val + abs(min(methyl_val))

# --- Viral reservoir: log2(x + 1) ---
# HIV DNA copy numbers span orders of magnitude; log2 reduces skewness,
# +1 pseudocount avoids log(0).
res_log_disc <- log2(res_disc + 1)
res_log_val  <- log2(res_val + 1)

gc()

# -----------------------------------------------------------------------------
# 7. BioMart Annotation
# -----------------------------------------------------------------------------
# Map ENSEMBL IDs to HGNC gene symbols; used to exclude Y chromosome genes
# and biotype artifacts before feature selection, and to create
# ENSEMBLID:GENE_SYMBOL feature labels for the bulk RNA-seq layer.

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations_bm <- getBM(
  filters    = "ensembl_gene_id",
  values     = colnames(bulk_vst_disc),
  mart       = mart,
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "chromosome_name")
) %>% as.data.frame()

feature_to_label <- setNames(
  paste(annotations_bm$ensembl_gene_id, annotations_bm$hgnc_symbol, sep = ":"),
  annotations_bm$ensembl_gene_id
)

# -----------------------------------------------------------------------------
# 8. Feature Selection (top MAD features per omics layer)
# -----------------------------------------------------------------------------
# Feature selection is performed on the discovery cohort only; the same
# features are then applied to the validation cohort.
#
# Selected features per layer:
#   - Bulk RNA-seq:      top 5,000 genes by MAD (Y chr and biotype artifacts excluded)
#   - DNA Methylation:   top 1% variable CpG sites
#   - All other layers:  all features retained

select_top_mad <- function(data, top_n, annotation_bm = NULL) {
  if (!is.null(annotation_bm)) {
    # Remove Y chromosome genes and gene biotype artifacts
    keep_genes <- annotation_bm$ensembl_gene_id[
      annotation_bm$gene_biotype != "artifact" &
        annotation_bm$chromosome_name != "Y"
    ]
    data <- data[, colnames(data) %in% keep_genes, drop = FALSE]
  }
  mad_vals     <- apply(data, 2, mad)
  top_features <- names(sort(mad_vals, decreasing = TRUE)[1:min(top_n, ncol(data))])
  data[, top_features]
}

# Bulk RNA-seq: top 5,000 genes by MAD
bulk_vst_disc_top <- select_top_mad(bulk_vst_disc, top_n = 5000, annotation_bm = annotations_bm)
bulk_vst_val_top  <- bulk_vst_val[, colnames(bulk_vst_disc_top), drop = FALSE]

# Rename bulk RNA-seq features to ENSEMBLID:GENE_SYMBOL format
colnames(bulk_vst_disc_top) <- feature_to_label[colnames(bulk_vst_disc_top)]
colnames(bulk_vst_val_top)  <- feature_to_label[colnames(bulk_vst_val_top)]

# DNA Methylation: top 1% variable CpG sites
top_n_methyl    <- round(ncol(methyl_disc_nn) * 0.01)
methyl_disc_top <- select_top_mad(methyl_disc_nn, top_n = top_n_methyl)
methyl_val_top  <- methyl_val_nn[, colnames(methyl_disc_top), drop = FALSE]

cat("\nFeatures selected per layer (discovery cohort):\n")
cat("  Bulk RNA-seq:      ", ncol(bulk_vst_disc_top), "\n")
cat("  Proteomics:        ", ncol(prot_log_disc), "\n")
cat("  Flow cytometry:    ", ncol(flow_log_disc), "\n")
cat("  Ex vivo cytokines: ", ncol(exvivo_log_disc), "\n")
cat("  Methylation:       ", ncol(methyl_disc_top), "\n")
cat("  Viral reservoir:   ", ncol(res_log_disc), "\n")

# -----------------------------------------------------------------------------
# 9. Save Preprocessed Data
# -----------------------------------------------------------------------------
# Objects below feed directly into multiomics_clustering.R:
#   bulk_vst_disc_top / bulk_vst_val_top
#   prot_log_disc / prot_log_val
#   flow_log_disc / flow_log_val
#   exvivo_log_disc / exvivo_log_val
#   methyl_disc_top / methyl_val_top
#   res_log_disc / res_log_val

saveRDS(
  list(
    Bulk_transcriptomics = bulk_vst_disc_top,
    Proteomics            = prot_log_disc,
    Flowcytometry_abs     = flow_log_disc,
    Ex_Vivo               = exvivo_log_disc,
    Methylation           = methyl_disc_top,
    Viral_Reservoir       = res_log_disc
  ),
  file.path(output_dir, "preprocessing/data_omics_discovery_preprocessed.rds")
)

saveRDS(
  list(
    Bulk_transcriptomics = bulk_vst_val_top,
    Proteomics            = prot_log_val,
    Flowcytometry_abs     = flow_log_val,
    Ex_Vivo               = exvivo_log_val,
    Methylation           = methyl_val_top,
    Viral_Reservoir       = res_log_val
  ),
  file.path(output_dir, "preprocessing/data_omics_validation_preprocessed.rds")
)

# Full (un-subsetted) reservoir table with clinical variables, exclusion-
# filtered but not yet restricted to the two clustering variables — used
# downstream to compare clustering solutions against clinical/reservoir
# measurements (see multiomics_clustering.R, Section 4).
saveRDS(viral_reservoir_clean, file.path(output_dir, "preprocessing/viral_reservoir_clean.rds"))

cat("\nPreprocessed data saved successfully.\n")