# =============================================================================
# Script: 01_data_loading_preprocessing.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Load all omics datasets, apply exclusion criteria, and perform
#              data preprocessing and normalization per omics layer.
# Author: Victoria Rios (Victoria.RiosVazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#              Nature Immunology (NI-A42808-T)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------
options(warn = -1)

# Core packages for data manipulation and I/O
library(dplyr)       # Data wrangling
library(readr)       # Fast CSV reading
library(readxl)      # Excel file reading
library(openxlsx)    # Excel file writing

# Omics-specific packages
library(MOVICS)      # Multi-omics integrative clustering pipeline
library(DESeq2)      # RNA-seq normalization (VST transformation)
library(impute)      # KNN imputation for missing values

# Set paths — update these to match your local directory structure
# project_dir: root of the project
# data_dir: directory containing all raw/processed omics data
# output_dir: directory where intermediate and final results will be saved
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

# -----------------------------------------------------------------------------
# 1. Load Clinical and Reservoir Data
# -----------------------------------------------------------------------------

# Main clinical sample table containing cohort assignment (DISCOVERY/VALIDATION),
# demographic variables, and HIV-related clinical parameters
sample_table <- read_excel(file.path(data_dir, "clinical/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.xlsx"))

# HIV-1 DNA reservoir quantification data (total and intact HIV DNA per million CD4+ T cells)
# Generated using the Rainbow proviral HIV-1 DNA digital PCR assay (Delporte et al. 2025)
viral_reservoir <- read_csv(file.path(data_dir, "reservoir/Reservoir_Rscript_processed_2025.csv")) %>%
  as.data.frame()

# Load exclusion lists:
# - "Exclusion" sheet: samples failing QC or inclusion criteria
# - "Separate_analysis" sheet: samples on immunomodulatory drugs (analyzed separately)
exclusion_IDs <- read_excel(
  file.path(data_dir, "clinical/20220322_Exclusion_2000HIV.xlsx"),
  sheet = "Exclusion", skip = 1
)
separate_analysis_IDs <- read_excel(
  file.path(data_dir, "clinical/20220322_Exclusion_2000HIV.xlsx"),
  sheet = "Separate_analysis", skip = 1
)

# Apply exclusion criteria — remove samples flagged in either exclusion list
viral_reservoir_clean <- viral_reservoir[
  !viral_reservoir$SampleID %in% c(exclusion_IDs$`Study ID`, separate_analysis_IDs$`Study ID`),
]
cat("Reservoir samples before exclusion:", nrow(viral_reservoir), "\n")
cat("Reservoir samples after exclusion:", nrow(viral_reservoir_clean), "\n")

# Prepare viral reservoir matrix with only the two key variables used in clustering:
# - Total_million_avg: total HIV-1 DNA copies per million CD4+ T cells (RU5 region)
# - intactDT_DSI: intact HIV-1 DNA copies per million CD4+ T cells (DSI-corrected)
viral_reservoir_df <- viral_reservoir_clean %>%
  dplyr::select(SampleID, Total_million_avg, intactDT_DSI) %>%
  tibble::column_to_rownames("SampleID")

# Define discovery (n=1,027) and validation (n=203) cohort IDs
# Note: sub-cohort assignment was made prior to any analytical preprocessing
discovery_ids <- sample_table[sample_table$COHORT == "DISCOVERY", ]$ID
validation_ids <- sample_table[sample_table$COHORT == "VALIDATION", ]$ID

# -----------------------------------------------------------------------------
# 2. Load Omics Datasets
# -----------------------------------------------------------------------------

# --- Bulk RNA-seq ---
sample_table_bulkrnaseq <- readRDS(file.path(data_dir, "clinical/2000HIV_bulk_transcriptomics_sample_table.RDS"))
sample_table_bulkrnaseq$ID <- sample_table_bulkrnaseq$DONOR_ID

bulk_transcriptomics_df <- readRDS(file.path(data_dir, "omics/2000HIV_bulk_transcriptomics/2000HIV_bulk_transcriptomics_raw_counts.RDS"))
bulk_transcriptomics_df <- t(bulk_transcriptomics_df[, !names(bulk_transcriptomics_df) %in% c("DESCRIPTION", "GENETYPE", "SYMBOL", "GENEID", "CHR")])

sample_table_bulkrnaseq$LAB_ID <- as.character(sample_table_bulkrnaseq$LAB_ID)
bulk_transcriptomics_df <- as.data.frame(bulk_transcriptomics_df)
bulk_transcriptomics_df$LAB_ID <- as.character(rownames(bulk_transcriptomics_df))

bulk_transcriptomics_df <- merge(
  bulk_transcriptomics_df,
  sample_table_bulkrnaseq[c("LAB_ID", "DONOR_ID")],
  by = "LAB_ID", all.x = TRUE
)
rownames(bulk_transcriptomics_df) <- bulk_transcriptomics_df$DONOR_ID
bulk_transcriptomics_df <- bulk_transcriptomics_df[, !names(bulk_transcriptomics_df) %in% c("LAB_ID", "DONOR_ID")]

# Filter to reservoir samples only
bulk_transcriptomics_df <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% viral_reservoir_clean$SampleID, ]
cat("Bulk RNA-seq:", dim(bulk_transcriptomics_df), "\n")

# --- Flow Cytometry (Absolute counts) ---
flow_abs_df <- read_excel(file.path(
  data_dir,
  "omics/2000HIV_Flow_Cytometry/2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_356vars_Nov062023.xlsx"
)) %>% as.data.frame()
rownames(flow_abs_df) <- flow_abs_df$SampleID
flow_abs_df <- flow_abs_df[, -1]
flow_abs_df <- flow_abs_df[rownames(flow_abs_df) %in% viral_reservoir_clean$SampleID, ]
cat("Flow cytometry (absolute):", dim(flow_abs_df), "\n")

# --- Ex Vivo Cytokine Data ---
ex_vivo_24h_df <- read_excel(file.path(data_dir, "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_24H_52meas_1760samples_afterQC_RAW.xlsx")) %>%
  as.data.frame()
rownames(ex_vivo_24h_df) <- ex_vivo_24h_df$ID
ex_vivo_24h_df <- ex_vivo_24h_df[, -1]

ex_vivo_7d_df <- read_excel(file.path(data_dir, "omics/2000HIV_Ex_Vivo_Cytokine/2000HIV_EX_VIVO_7DAY_38meas_1754samples_afterQC_RAW.xlsx")) %>%
  as.data.frame()
rownames(ex_vivo_7d_df) <- ex_vivo_7d_df$ID
ex_vivo_7d_df <- ex_vivo_7d_df[, -1]

# Merge 24h and 7d cytokine data
common_ids <- intersect(rownames(ex_vivo_24h_df), rownames(ex_vivo_7d_df))
ex_vivo_df <- cbind(
  ex_vivo_24h_df[common_ids, ],
  ex_vivo_7d_df[common_ids, ]
)
ex_vivo_df <- ex_vivo_df[rownames(ex_vivo_df) %in% viral_reservoir_clean$SampleID, ]
cat("Ex vivo cytokine:", dim(ex_vivo_df), "\n")

# --- DNA Methylation (M-values) ---
methylation_df <- readRDS(file.path(data_dir, "omics/2000HIV.Mvalue.rds")) %>% t()
methylation_df <- methylation_df[rownames(methylation_df) %in% viral_reservoir_clean$SampleID, ]
cat("DNA Methylation:", dim(methylation_df), "\n")

gc()

# -----------------------------------------------------------------------------
# 3. Split into Discovery and Validation Cohorts
# -----------------------------------------------------------------------------

bulk_transcriptomics_discovery <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% discovery_ids, ]
bulk_transcriptomics_validation <- bulk_transcriptomics_df[rownames(bulk_transcriptomics_df) %in% validation_ids, ]

flow_abs_discovery <- flow_abs_df[rownames(flow_abs_df) %in% discovery_ids, ]
flow_abs_validation <- flow_abs_df[rownames(flow_abs_df) %in% validation_ids, ]

ex_vivo_discovery <- ex_vivo_df[rownames(ex_vivo_df) %in% discovery_ids, ]
ex_vivo_validation <- ex_vivo_df[rownames(ex_vivo_df) %in% validation_ids, ]

methylation_discovery <- methylation_df[rownames(methylation_df) %in% discovery_ids, ]
methylation_validation <- methylation_df[rownames(methylation_df) %in% validation_ids, ]

viral_reservoir_discovery <- viral_reservoir_df[rownames(viral_reservoir_df) %in% discovery_ids, ]
viral_reservoir_validation <- viral_reservoir_df[rownames(viral_reservoir_df) %in% validation_ids, ]

# -----------------------------------------------------------------------------
# 4. Normalization and Transformation
# -----------------------------------------------------------------------------
# Each omics layer requires a different transformation strategy to stabilize
# variance and meet the assumptions of the MOVICS clustering algorithms.
# All data types are treated as gaussian (continuous) in the clustering step.

# --- Bulk RNA-seq: DESeq2 Variance Stabilizing Transformation (VST) ---
# VST is applied separately per sub-cohort to avoid data leakage between
# discovery and validation sets. This removes the dependency of variance
# on the mean, producing approximately homoskedastic data suitable for
# downstream analyses. Normalization is performed blind to any conditions.
sample_table_temp <- sample_table[sample_table$ID %in% c(discovery_ids, validation_ids), ]

bulk_t <- t(bulk_transcriptomics_df)
bulk_t <- bulk_t[, colnames(bulk_t) %in% sample_table$ID]

# Filter low-count genes: keep genes with >= 5 counts in >= 50% of samples
# This removes noise from lowly expressed genes prior to normalization
keep <- rowSums(bulk_t >= 5) >= (ncol(bulk_t) / 2)
bulk_t_filtered <- bulk_t[keep, ]

# Discovery VST
dds_discovery <- DESeqDataSetFromMatrix(
  countData = bulk_t_filtered[, colnames(bulk_t_filtered) %in% discovery_ids],
  colData   = sample_table_temp[sample_table_temp$ID %in% discovery_ids, ],
  design    = ~1  # No experimental design — unsupervised normalization
)
dds_discovery <- DESeq(dds_discovery)
bulk_vst_discovery <- assay(varianceStabilizingTransformation(estimateSizeFactors(dds_discovery), blind = TRUE)) %>% t()

# Validation VST — applied independently to avoid data leakage
dds_validation <- DESeqDataSetFromMatrix(
  countData = bulk_t_filtered[, colnames(bulk_t_filtered) %in% validation_ids],
  colData   = sample_table_temp[sample_table_temp$ID %in% validation_ids, ],
  design    = ~1
)
dds_validation <- DESeq(dds_validation)
bulk_vst_validation <- assay(varianceStabilizingTransformation(estimateSizeFactors(dds_validation), blind = TRUE)) %>% t()

cat("Bulk RNA-seq VST discovery:", dim(bulk_vst_discovery), "\n")
cat("Bulk RNA-seq VST validation:", dim(bulk_vst_validation), "\n")

# --- Flow Cytometry: log2(x + 1) transformation ---
# Absolute immune cell counts are right-skewed; log2 transformation
# compresses the dynamic range and reduces the influence of extreme values.
# Adding 1 before log-transformation avoids log(0) issues (pseudocount).
fillna <- function(df, value = 0) { df[is.na(df)] <- value; df }

flow_log_discovery  <- log2(fillna(flow_abs_discovery, 0) + 1)
flow_log_validation <- log2(fillna(flow_abs_validation, 0) + 1)

# --- Ex Vivo Cytokine: log2(x + 1) transformation ---
# Cytokine concentrations are highly skewed due to stimulation-induced
# outliers; log2 transformation brings distributions closer to gaussian.
exvivo_log_discovery  <- log2(fillna(ex_vivo_discovery, 0) + 1)
exvivo_log_validation <- log2(fillna(ex_vivo_validation, 0) + 1)

# --- DNA Methylation: shift to non-negative ---
# M-values (logit-transformed beta values) can be negative.
# For MOVICS, we shift the entire matrix by its minimum value to ensure
# all values are non-negative, without altering the distribution shape.
methylation_discovery_nn  <- methylation_discovery + abs(min(methylation_discovery))
methylation_validation_nn <- methylation_validation + abs(min(methylation_validation))

# --- Viral Reservoir: log2(x + 1) transformation ---
# HIV DNA copy numbers span several orders of magnitude; log2 transformation
# stabilizes variance and reduces the influence of extreme reservoir sizes.
# Adding 1 avoids log(0) for samples with undetectable intact HIV DNA.
reservoir_log_discovery  <- log2(viral_reservoir_discovery + 1)
reservoir_log_validation <- log2(viral_reservoir_validation + 1)

gc()

# -----------------------------------------------------------------------------
# 5. Feature Selection (top MAD features per omics layer)
# -----------------------------------------------------------------------------
# To reduce noise from low-variability features and address dimensionality
# imbalance between omics layers, we select the most variable features per
# layer using Mean Absolute Deviation (MAD). MAD is preferred over variance
# as it is more robust to outliers.
#
# Feature selection is performed on the discovery cohort only.
# The same features are then applied to the validation cohort to ensure
# that no information from validation samples influences feature selection.
#
# Selected features per layer (as described in the Methods):
#   - Bulk RNA-seq: top 5,000 genes (by MAD, excluding Y chromosome genes)
#   - DNA Methylation: top 1% variable CpG sites (7,938 sites)
#   - Flow cytometry, ex vivo cytokines, viral reservoir: all features retained

select_top_mad <- function(data, top_n) {
  # Calculate MAD for each feature (column)
  mad_vals <- apply(data, 2, mad)
  # Select the top_n features with the highest MAD values
  top_features <- names(sort(mad_vals, decreasing = TRUE)[1:top_n])
  data[, top_features]
}

# Select top features from discovery cohort
bulk_vst_discovery_top  <- select_top_mad(bulk_vst_discovery, top_n = 5000)
methylation_discovery_top  <- select_top_mad(methylation_discovery_nn, top_n = 7938)

# Apply the same feature set to validation cohort (no re-selection)
bulk_vst_validation_top <- bulk_vst_validation[, colnames(bulk_vst_discovery_top), drop = FALSE]
methylation_validation_top <- methylation_validation_nn[, colnames(methylation_discovery_top), drop = FALSE]

cat("Feature selection complete.\n")
cat("Bulk RNA-seq features selected:", ncol(bulk_vst_discovery_top), "\n")
cat("Methylation features selected:", ncol(methylation_discovery_top), "\n")

# -----------------------------------------------------------------------------
# 6. Build MOVICS Input Data Lists
# -----------------------------------------------------------------------------

# Helper: align samples across layers and transpose to features x samples
build_omics_list <- function(bulk, flow, exvivo, methylation, reservoir) {
  # Find common samples across all layers
  common_samples <- Reduce(intersect, list(
    rownames(bulk),
    rownames(flow),
    rownames(exvivo),
    rownames(methylation),
    rownames(reservoir)
  ))
  cat("Common samples:", length(common_samples), "\n")

  list(
    Bulk.transcriptomics = t(bulk[common_samples, ]),
    Flowcytometry.abs    = t(flow[common_samples, ]),
    Ex.Vivo              = t(exvivo[common_samples, ]),
    Methylation          = t(methylation[common_samples, ]),
    Viral.Reservoir      = t(reservoir[common_samples, ])
  )
}

data_omics_discovery <- build_omics_list(
  bulk_vst_discovery_top,
  flow_log_discovery,
  exvivo_log_discovery,
  methylation_discovery_top,
  reservoir_log_discovery
)

data_omics_validation <- build_omics_list(
  bulk_vst_validation_top,
  flow_log_validation,
  exvivo_log_validation,
  methylation_validation_top,
  reservoir_log_validation
)

# -----------------------------------------------------------------------------
# 7. KNN Imputation of Missing Values
# -----------------------------------------------------------------------------

data_omics_discovery <- lapply(data_omics_discovery, function(mat) {
  impute.knn(as.matrix(mat), k = 5)$data
})

data_omics_validation <- lapply(data_omics_validation, function(mat) {
  impute.knn(as.matrix(mat), k = 5)$data
})

cat("Imputation complete.\n")
lapply(data_omics_discovery, dim)
lapply(data_omics_validation, dim)

# -----------------------------------------------------------------------------
# 8. Save Preprocessed Data
# -----------------------------------------------------------------------------

saveRDS(data_omics_discovery,
        file.path(output_dir, "data_omics_discovery_preprocessed.rds"))
saveRDS(data_omics_validation,
        file.path(output_dir, "data_omics_validation_preprocessed.rds"))

cat("Preprocessed data saved successfully.\n")
