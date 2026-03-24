# =============================================================================
# Script: 01_data_loading_preprocessing.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Load all omics datasets, apply exclusion criteria, perform
#              data preprocessing, normalization, feature selection, imputation,
#              sample overlap filtering, and prepare MOVICS input lists.
# Author: Victoria Rios (Victoria.RiosVazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
options(warn = -1)

library(dplyr)       # Data wrangling
library(readr)       # Fast CSV reading
library(readxl)      # Excel file reading
library(openxlsx)    # Excel file writing
library(MOVICS)      # Multi-omics integrative clustering pipeline
library(DESeq2)      # RNA-seq normalization (VST transformation)
library(impute)      # KNN imputation for missing values
library(biomaRt)     # Gene annotation (ENSEMBL to HGNC symbol mapping)

# Set paths — update these to match your local directory structure
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

# -----------------------------------------------------------------------------
# 1. Load Clinical and Reservoir Data
# -----------------------------------------------------------------------------

# Main clinical sample table with cohort assignment (DISCOVERY/VALIDATION),
# demographics, and HIV-related clinical parameters
sample_table <- read_excel(file.path(data_dir, "clinical/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.xlsx"))

# HIV-1 DNA reservoir quantification data
# Generated using the Rainbow proviral HIV-1 DNA digital PCR assay (Delporte et al. 2025)
viral_reservoir <- read_csv(file.path(data_dir, "reservoir/Reservoir_Rscript_processed_2025.csv")) %>%
  as.data.frame()

# Load exclusion lists
# - "Exclusion": samples failing QC or inclusion criteria
# - "Separate_analysis": samples on immunomodulatory drugs
exclusion_IDs <- read_excel(
  file.path(data_dir, "clinical/20220322_Exclusion_2000HIV.xlsx"),
  sheet = "Exclusion", skip = 1
)
separate_analysis_IDs <- read_excel(
  file.path(data_dir, "clinical/20220322_Exclusion_2000HIV.xlsx"),
  sheet = "Separate_analysis", skip = 1
)

# Apply exclusion criteria
viral_reservoir_clean <- viral_reservoir[
  !viral_reservoir$SampleID %in% c(exclusion_IDs$`Study ID`, separate_analysis_IDs$`Study ID`),
]
cat("Reservoir samples before exclusion:", nrow(viral_reservoir), "\n")
cat("Reservoir samples after exclusion:", nrow(viral_reservoir_clean), "\n")

# Prepare viral reservoir dataframe
# Select only the key quantitative variables used in clustering:
#   - Total_million_avg: total HIV-1 DNA copies per million CD4+ T cells
#   - intactDT_DSI: intact HIV-1 DNA copies per million CD4+ T cells (DSI-corrected)
# All columns coerced to numeric to avoid type issues downstream
viral_reservoir_df <- viral_reservoir_clean[!is.na(viral_reservoir_clean$SampleID), ]
rownames(viral_reservoir_df) <- viral_reservoir_df$SampleID
viral_reservoir_df <- viral_reservoir_df[, colnames(viral_reservoir_df) %in% c(
  "Total_million_avg",
  "intactDT_DSI"
)]
viral_reservoir_df <- viral_reservoir_df %>% mutate(across(everything(), as.numeric))
cat("Viral reservoir matrix dimensions:", dim(viral_reservoir_df), "\n")

# Define discovery and validation cohort IDs
# Sub-cohort assignment was made prior to any analytical preprocessing
discovery_ids <- sample_table[sample_table$COHORT == "DISCOVERY", ]$ID
validation_ids <- sample_table[sample_table$COHORT == "VALIDATION", ]$ID

# -----------------------------------------------------------------------------
# 2. Load Omics Datasets
# -----------------------------------------------------------------------------

# --- Bulk RNA-seq (raw counts) ---
sample_table_bulkrnaseq <- readRDS(file.path(data_dir, "clinical/2000HIV_bulk_transcriptomics_sample_table.RDS"))
sample_table_bulkrnaseq$ID <- sample_table_bulkrnaseq$DONOR_ID
bulk_transcriptomics_df <- readRDS(file.path(data_dir, "omics/2000HIV_bulk_transcriptomics/2000HIV_bulk_transcriptomics_raw_counts.RDS"))
bulk_transcriptomics_df <- t(bulk_transcriptomics_df[, !names(bulk_transcriptomics_df) %in% c("DESCRIPTION", "GENETYPE", "SYMBOL", "GENEID", "CHR")])
sample_table_bulkrnaseq$LAB_ID <- as.character(sample_table_bulkrnaseq$LAB_ID)
bulk_transcriptomics_df <- as.data.frame(bulk_transcriptomics_df)
bulk_transcriptomics_df$LAB_ID <- as.character(rownames(bulk_transcriptomics_df))
bulk_transcriptomics_df <- merge(bulk_transcriptomics_df, sample_table_bulkrnaseq[c("LAB_ID", "DONOR_ID")], by = "LAB_ID", all.x = TRUE)
rownames(bulk_transcriptomics_df) <- bulk_transcriptomics_df$DONOR_ID
bulk_transcriptomics_df <- bulk_transcriptomics_df[, !names(bulk_transcriptomics_df) %in% c("LAB_ID", "DONOR_ID")]
bulk_transcriptomics_df <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% viral_reservoir_clean$SampleID, ]
cat("Bulk RNA-seq:", dim(bulk_transcriptomics_df), "\n")

# --- Flow Cytometry (absolute counts, QC-processed Nov 2023) ---
# Final dataset after 3 rounds of QC filtering — 355 immune cell subset variables
flow_abs_df <- read_excel(file.path(data_dir, "omics/2000HIV_Flow_Cytometry/2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_356vars_Nov062023.xlsx")) %>%
  as.data.frame()
rownames(flow_abs_df) <- flow_abs_df$SampleID
flow_abs_df <- flow_abs_df[, -1]
flow_abs_df <- flow_abs_df[rownames(flow_abs_df) %in% viral_reservoir_clean$SampleID, ]
cat("Flow cytometry (absolute):", dim(flow_abs_df), "\n")

# --- Ex Vivo Cytokine Data (24h and 7-day stimulations merged) ---
ex_vivo_24h_df <- read_excel(file.path(data_dir, "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_RAW.xlsx")) %>% as.data.frame()
rownames(ex_vivo_24h_df) <- ex_vivo_24h_df$ID
ex_vivo_24h_df <- ex_vivo_24h_df[, -1]
ex_vivo_24h_df <- ex_vivo_24h_df[rownames(ex_vivo_24h_df) %in% viral_reservoir_clean$SampleID, ]

ex_vivo_7d_df <- read_excel(file.path(data_dir, "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_7DAY_38meas_1754samples_afterQC_RAW.xlsx")) %>% as.data.frame()
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
# 3. Split into Discovery and Validation Cohorts
# -----------------------------------------------------------------------------

bulk_disc    <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% discovery_ids, ]
bulk_val     <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% validation_ids, ]
flow_disc    <- flow_abs_df[rownames(flow_abs_df) %in% discovery_ids, ]
flow_val     <- flow_abs_df[rownames(flow_abs_df) %in% validation_ids, ]
exvivo_disc  <- ex_vivo_df[rownames(ex_vivo_df) %in% discovery_ids, ]
exvivo_val   <- ex_vivo_df[rownames(ex_vivo_df) %in% validation_ids, ]
methyl_disc  <- methylation_df[rownames(methylation_df) %in% discovery_ids, ]
methyl_val   <- methylation_df[rownames(methylation_df) %in% validation_ids, ]
res_disc     <- viral_reservoir_df[rownames(viral_reservoir_df) %in% discovery_ids, ]
res_val      <- viral_reservoir_df[rownames(viral_reservoir_df) %in% validation_ids, ]

# -----------------------------------------------------------------------------
# 4. Normalization and Transformation
# -----------------------------------------------------------------------------

# --- Bulk RNA-seq: DESeq2 Variance Stabilizing Transformation (VST) ---
# Applied independently per sub-cohort to prevent data leakage.
# Removes version suffixes from ENSEMBL IDs (e.g. ENSG00000001.5 -> ENSG00000001)
sample_table_temp <- sample_table[sample_table$ID %in% c(discovery_ids, validation_ids), ]
bulk_t <- t(bulk_transcriptomics_df)
bulk_t <- bulk_t[, colnames(bulk_t) %in% sample_table$ID]
keep <- rowSums(bulk_t >= 5) >= (ncol(bulk_t) / 2)
bulk_t_filtered <- bulk_t[keep, ]

strip_version <- function(mat) { colnames(mat) <- gsub("\\.\\d+$", "", colnames(mat)); mat }

dds_disc <- DESeqDataSetFromMatrix(
  countData = bulk_t_filtered[, colnames(bulk_t_filtered) %in% discovery_ids],
  colData   = sample_table_temp[sample_table_temp$ID %in% discovery_ids, ],
  design    = ~1
)
dds_disc <- DESeq(dds_disc)
bulk_vst_disc <- assay(varianceStabilizingTransformation(estimateSizeFactors(dds_disc), blind = TRUE)) %>% t() %>% strip_version()

dds_val <- DESeqDataSetFromMatrix(
  countData = bulk_t_filtered[, colnames(bulk_t_filtered) %in% validation_ids],
  colData   = sample_table_temp[sample_table_temp$ID %in% validation_ids, ],
  design    = ~1
)
dds_val <- DESeq(dds_val)
bulk_vst_val <- assay(varianceStabilizingTransformation(estimateSizeFactors(dds_val), blind = TRUE)) %>% t() %>% strip_version()

# --- BiomaRt annotation ---
# Map ENSEMBL IDs to HGNC gene symbols; used to exclude Y chromosome genes
# and artifacts, and to create ENSEMBLID:GENE_SYMBOL feature labels
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations_bm <- getBM(
  filters    = "ensembl_gene_id",
  values     = colnames(bulk_vst_disc),
  mart       = mart,
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "chromosome_name")
) %>% as.data.frame()

# --- Flow Cytometry: log2(x + 1) ---
# Right-skewed counts; log2 stabilizes variance
fillna <- function(df, value = 0) { df[is.na(df)] <- value; df }
flow_log_disc <- log2(fillna(flow_disc, 0) + 1)
flow_log_val  <- log2(fillna(flow_val, 0) + 1)

# --- Ex Vivo Cytokine: log2(x + 1) ---
exvivo_log_disc <- log2(fillna(exvivo_disc, 0) + 1)
exvivo_log_val  <- log2(fillna(exvivo_val, 0) + 1)

# --- DNA Methylation: shift to non-negative ---
# M-values can be negative; add absolute minimum to shift without changing shape
methyl_disc_nn <- methyl_disc + abs(min(methyl_disc))
methyl_val_nn  <- methyl_val + abs(min(methyl_val))

# --- Viral Reservoir: log2(x + 1) ---
# HIV DNA copy numbers span orders of magnitude; log2 reduces skewness
# +1 pseudocount avoids log(0)
res_log_disc <- log2(res_disc + 1)
res_log_val  <- log2(res_val + 1)

gc()

# -----------------------------------------------------------------------------
# 5. Feature Selection (top MAD features per omics layer)
# -----------------------------------------------------------------------------
# Feature selection is performed on discovery cohort only;
# the same features are then applied to the validation cohort.
#
# Selected features per layer (Methods):
#   - Bulk RNA-seq:      top 5,000 genes by MAD (Y chr and artifacts excluded)
#   - DNA Methylation:   top 1% variable CpG sites (~7,938 sites)
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

# Bulk RNA-seq
bulk_vst_disc_top <- select_top_mad(bulk_vst_disc, top_n = 5000, annotation_bm = annotations_bm)
bulk_vst_val_top  <- bulk_vst_val[, colnames(bulk_vst_disc_top), drop = FALSE]

# Create ENSEMBLID:GENE_SYMBOL label mapping for feature renaming
feature_to_label <- setNames(
  paste(annotations_bm$ensembl_gene_id, annotations_bm$hgnc_symbol, sep = ":"),
  annotations_bm$ensembl_gene_id
)

# DNA Methylation: top 1% variable CpG sites
top_n_methyl    <- round(ncol(methyl_disc_nn) * 0.01)
methyl_disc_top <- select_top_mad(methyl_disc_nn, top_n = top_n_methyl)
methyl_val_top  <- methyl_val_nn[, colnames(methyl_disc_top), drop = FALSE]

cat("Bulk RNA-seq features selected:", ncol(bulk_vst_disc_top), "\n")
cat("Methylation features selected:", ncol(methyl_disc_top), "\n")

# -----------------------------------------------------------------------------
# 6. Build MOVICS Input Matrices (features x samples)
# -----------------------------------------------------------------------------

# Helper: create feature x sample matrix, filling missing samples with NA
convert_to_matrix <- function(data, all_samples) {
  res_mat <- matrix(NA, nrow = ncol(data), ncol = length(all_samples),
                    dimnames = list(colnames(data), all_samples))
  common  <- intersect(rownames(data), all_samples)
  res_mat[, common] <- t(data[common, ])
  res_mat
}

all_samples_disc <- unique(c(rownames(bulk_vst_disc_top), rownames(flow_log_disc),
                             rownames(exvivo_log_disc), rownames(methyl_disc_top),
                             rownames(res_log_disc)))
all_samples_val  <- unique(c(rownames(bulk_vst_val_top), rownames(flow_log_val),
                             rownames(exvivo_log_val), rownames(methyl_val_top),
                             rownames(res_log_val)))

# Build matrices
bulk_mat_disc  <- convert_to_matrix(bulk_vst_disc_top, all_samples_disc)
flow_mat_disc  <- convert_to_matrix(flow_log_disc,      all_samples_disc)
exvivo_mat_disc <- convert_to_matrix(exvivo_log_disc,   all_samples_disc)
methyl_mat_disc <- convert_to_matrix(methyl_disc_top,   all_samples_disc)
res_mat_disc    <- convert_to_matrix(res_log_disc,       all_samples_disc)

bulk_mat_val   <- convert_to_matrix(bulk_vst_val_top,  all_samples_val)
flow_mat_val   <- convert_to_matrix(flow_log_val,       all_samples_val)
exvivo_mat_val <- convert_to_matrix(exvivo_log_val,     all_samples_val)
methyl_mat_val <- convert_to_matrix(methyl_val_top,     all_samples_val)
res_mat_val    <- convert_to_matrix(res_log_val,        all_samples_val)

# Rename bulk RNA-seq rows to ENSEMBLID:GENE_SYMBOL format
rownames(bulk_mat_disc) <- feature_to_label[rownames(bulk_mat_disc)]
rownames(bulk_mat_val)  <- feature_to_label[rownames(bulk_mat_val)]

# Keep only the two reservoir variables used in clustering
res_mat_disc <- res_mat_disc[c("Total_million_avg", "intactDT_DSI"), ]
res_mat_val  <- res_mat_val[c("Total_million_avg", "intactDT_DSI"), ]

# Assemble MOVICS lists (pre-imputation)
data_omics_disc_raw <- list(
  Bulk_transcriptomics = bulk_mat_disc,
  Flowcytometry_abs    = flow_mat_disc,
  Methylation          = methyl_mat_disc,
  Ex_Vivo              = exvivo_mat_disc,
  Viral_Reservoir      = res_mat_disc
)
data_omics_val_raw <- list(
  Bulk_transcriptomics = bulk_mat_val,
  Flowcytometry_abs    = flow_mat_val,
  Methylation          = methyl_mat_val,
  Ex_Vivo              = exvivo_mat_val,
  Viral_Reservoir      = res_mat_val
)

# -----------------------------------------------------------------------------
# 7. KNN Imputation of Missing Values
# -----------------------------------------------------------------------------
# Features with >= 80% missing values are excluded before imputation.
# KNN imputation (k=5) applied separately per layer and per sub-cohort
# to prevent any cross-cohort information leakage.

impute_matrix <- function(mat, missing_threshold = 0.8) {
  missing_prop <- rowMeans(is.na(mat))
  mat_filtered <- mat[missing_prop < missing_threshold, , drop = FALSE]
  cat("  Features removed (>= 80% missing):", sum(missing_prop >= missing_threshold), "\n")
  impute.knn(as.matrix(mat_filtered), k = 5)$data
}

cat("\nImputing discovery cohort...\n")
data_omics_disc_imputed <- lapply(data_omics_disc_raw, impute_matrix)

cat("\nImputing validation cohort...\n")
data_omics_val_imputed <- lapply(data_omics_val_raw, impute_matrix)

# -----------------------------------------------------------------------------
# 8. Filter to Overlapping Samples Across All Omics Layers
# -----------------------------------------------------------------------------
# Retain only samples present in ALL omics layers simultaneously.
# This ensures a complete data matrix for every sample entering clustering.

overlapping_disc <- Reduce(intersect, lapply(data_omics_disc_imputed, colnames))
overlapping_val  <- Reduce(intersect, lapply(data_omics_val_imputed, colnames))

cat("\nOverlapping samples — Discovery:", length(overlapping_disc), "\n")
cat("Overlapping samples — Validation:", length(overlapping_val), "\n")

data_omics_disc_filtered <- lapply(data_omics_disc_imputed, function(mat) mat[, overlapping_disc, drop = FALSE])
data_omics_val_filtered  <- lapply(data_omics_val_imputed,  function(mat) mat[, overlapping_val,  drop = FALSE])

# -----------------------------------------------------------------------------
# 9. Feature Name Cleaning for CIMLR Compatibility
# -----------------------------------------------------------------------------
# CIMLR does not handle special characters (_, +, -) in feature names.
# Replace: underscores -> dots, + -> "pos", - -> "neg"

clean_names <- function(x) {
  x <- gsub("_", ".", x)
  x <- gsub("\\+", "pos", x)
  x <- gsub("-", "neg", x)
  x
}

for (layer in c("Flowcytometry_abs", "Ex_Vivo", "Viral_Reservoir")) {
  rownames(data_omics_disc_filtered[[layer]]) <- clean_names(rownames(data_omics_disc_filtered[[layer]]))
  rownames(data_omics_val_filtered[[layer]])  <- clean_names(rownames(data_omics_val_filtered[[layer]]))
}

# Clean list names too
names(data_omics_disc_filtered) <- gsub("_", ".", names(data_omics_disc_filtered))
names(data_omics_val_filtered)  <- gsub("_", ".", names(data_omics_val_filtered))

cat("\nFinal list names:", names(data_omics_disc_filtered), "\n")

# -----------------------------------------------------------------------------
# 10. Save Final Preprocessed Data
# -----------------------------------------------------------------------------

saveRDS(data_omics_disc_filtered, file.path(output_dir, "data_omics_discovery_preprocessed.rds"))
saveRDS(data_omics_val_filtered,  file.path(output_dir, "data_omics_validation_preprocessed.rds"))

cat("\nPreprocessed data saved successfully.\n")
cat("Discovery dimensions per layer:\n")
lapply(data_omics_disc_filtered, dim)
cat("Validation dimensions per layer:\n")
lapply(data_omics_val_filtered, dim)
