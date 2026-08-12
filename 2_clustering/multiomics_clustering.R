#!/usr/bin/env Rscript
# =============================================================================
# Script: multiomics_clustering.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Build the MOVICS-formatted multi-omics input object from the
#              preprocessed omics layers, impute missing values, run the
#              final panel of MOVICS integrative clustering algorithms (SNF,
#              PINSPlus, NEMO, COCA, ConsensusClustering, CIMLR, MoCluster),
#              and visualize each solution with comprehensive heatmaps for
#              model comparison/selection.
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Requires: outputs from QC_and_preprocessing.R (normalized, log-transformed
#           omics layers; discovery/validation cohort splits; sample tables)
#
# MoCluster was selected as the reported clustering solution based on the
# downstream clinical-feature comparisons (see differential_expression/ and
# figures/ for those steps).
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(impute)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)

# Set paths: update these to match your local directory structure
# (must match the paths used in QC_and_preprocessing.R)
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"


# -----------------------------------------------------------------------------
# 1. Build MOVICS Input Matrices
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1.0 Load Preprocessed Data
# -----------------------------------------------------------------------------
# Loads the per-layer normalized, transformed, and feature-selected data
# produced by QC_and_preprocessing.R and unpacks it into the per-layer
# objects used throughout this section. Bulk RNA-seq features are already
# labeled as ENSEMBLID:GENE_SYMBOL at this point.

preproc_discovery  <- readRDS(file.path(output_dir, "preprocessing/data_omics_discovery_preprocessed.rds"))
preproc_validation <- readRDS(file.path(output_dir, "preprocessing/data_omics_validation_preprocessed.rds"))

# Clinical sample table (cohort assignment): used to align sample IDs below
sample_table <- read_excel(file.path(data_dir, "clinical/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.xlsx"))

# Full reservoir table with clinical variables (exclusion-filtered), used in
# Section 4 to compare clustering solutions against reservoir measurements
viral_reservoir_clean <- readRDS(file.path(output_dir, "preprocessing/viral_reservoir_clean.rds"))

bulk_transcriptomics_df_normalized_vst_high_variance_discovery  <- preproc_discovery$Bulk_transcriptomics
proteomics_df_log_discovery                                     <- preproc_discovery$Proteomics
flow_abs_df_log_discovery                                       <- preproc_discovery$Flowcytometry_abs
ex_vivo_df_log_discovery                                        <- preproc_discovery$Ex_Vivo
methylation_df_high_variance_discovery                          <- preproc_discovery$Methylation
viral_reservoir_df_log_discovery                                <- preproc_discovery$Viral_Reservoir

bulk_transcriptomics_df_normalized_vst_high_variance_validation <- preproc_validation$Bulk_transcriptomics
proteomics_df_log_validation                                    <- preproc_validation$Proteomics
flow_abs_df_log_validation                                      <- preproc_validation$Flowcytometry_abs
ex_vivo_df_log_validation                                       <- preproc_validation$Ex_Vivo
methylation_df_high_variance_validation                         <- preproc_validation$Methylation
viral_reservoir_df_log_validation                               <- preproc_validation$Viral_Reservoir


gc()

# MOVICS requires each omics layer as a matrix with features in rows and
# samples in columns, all layers sharing the same sample (column) order,
# and missing values represented as NA (rather than dropping samples with
# incomplete data). The block below builds that structure for one cohort
# at a time. Two of the reservoir table's columns are used for clustering
# (see column list explained in QC_and_preprocessing.R): the total and the
# intact HIV-1 DNA measures.
viral_reservoir_selected_vars <- c("Total_million_avg", "intactDT_DSI")

# Step 1: Union of sample IDs across all six layers (a sample missing from
# one layer will simply get NA columns for that layer once matrices are
# built), then restrict to IDs that exist in the clinical sample table
all_samples_discovery <- unique(c(
  rownames(bulk_transcriptomics_df_normalized_vst_high_variance_discovery),
  rownames(proteomics_df_log_discovery),
  rownames(ex_vivo_df_log_discovery),
  rownames(flow_abs_df_log_discovery),
  rownames(methylation_df_high_variance_discovery),
  rownames(viral_reservoir_df_log_discovery[,c(viral_reservoir_selected_vars)])
))

all_samples_discovery <- sample_table[sample_table$ID %in% all_samples_discovery,]$ID

# Step 2: Reindex a samples x features dataframe to the full sample list
# (introducing NA rows for missing samples), coerce to numeric, then
# transpose to the features x samples orientation MOVICS expects
convert_to_matrix <- function(df, all_samples) {
  df <- df[match(all_samples, rownames(df)), ]
  df[is.na(df)] <- NA  # Replace missing values with actual NA
  df <- apply(df, 2, as.numeric)  # Convert each column to numeric
  colnames(df) <- colnames(df)  # Set unique column names 
  rownames(df) <- all_samples  # Set row names to sample IDs
  res_mat <- t(as.matrix(df))
  rownames(res_mat) <- make.unique(rownames(res_mat), sep = "_")  # Make row names unique
  return(res_mat)
}

# Step 3: Apply the conversion to every layer
bulk_transcriptomics_mat_discovery <- convert_to_matrix(bulk_transcriptomics_df_normalized_vst_high_variance_discovery,
                                                        all_samples_discovery)
proteomics_mat_discovery <- convert_to_matrix(proteomics_df_log_discovery, all_samples_discovery)
flow_abs_mat_discovery <- convert_to_matrix(flow_abs_df_log_discovery, all_samples_discovery)
methylation_mat_discovery <- convert_to_matrix(methylation_df_high_variance_discovery, all_samples_discovery)
ex_vivo_mat_discovery <- convert_to_matrix(ex_vivo_df_log_discovery, all_samples_discovery)
viral_reservoir_mat_discovery <- convert_to_matrix(viral_reservoir_df_log_discovery, all_samples_discovery)

# Bulk transcriptomics features are already labeled as ENSEMBLID:GENE_SYMBOL
# by QC_and_preprocessing.R, so no additional renaming is needed here.

# Filter to keep only the selected reservoir variables
viral_reservoir_mat_discovery <- viral_reservoir_mat_discovery[c("Total_million_avg", "intactDT_DSI"), ]

# Step 4: Create a nested list of matrices
data_omics_discovery <- list(
  Bulk_transcriptomics = bulk_transcriptomics_mat_discovery,
  Proteomics = proteomics_mat_discovery,
  Flowcytometry_abs = flow_abs_mat_discovery,
  Methylation = methylation_mat_discovery,
  Ex_Vivo = ex_vivo_mat_discovery,
  Viral_Reservoir = viral_reservoir_mat_discovery
)

summary(data_omics_discovery)
saveRDS(data_omics_discovery, file.path(output_dir, "clustering/data_omics_discovery_preimputed.rds"))
gc()

# Same construction as above (Step 1-4), repeated for the validation cohort
viral_reservoir_selected_vars <- c("Total_million_avg", "intactDT_DSI")

# Step 1: Union of sample IDs across all six layers, restricted to IDs
# present in the clinical sample table
all_samples_validation <- unique(c(
  rownames(bulk_transcriptomics_df_normalized_vst_high_variance_validation),
  rownames(proteomics_df_log_validation),
  rownames(ex_vivo_df_log_validation),
  rownames(flow_abs_df_log_validation),
  rownames(methylation_df_high_variance_validation),
  rownames(viral_reservoir_df_log_validation[,c(viral_reservoir_selected_vars)])
))

all_samples_validation <- sample_table[sample_table$ID %in% all_samples_validation,]$ID

# Step 2: Same helper function as for the discovery cohort above
convert_to_matrix <- function(df, all_samples) {
  df <- df[match(all_samples, rownames(df)), ]
  df[is.na(df)] <- NA  # Replace missing values with actual NA
  df <- apply(df, 2, as.numeric)  # Convert each column to numeric
  colnames(df) <- colnames(df)  # Set unique column names 
  rownames(df) <- all_samples  # Set row names to sample IDs
  res_mat <- t(as.matrix(df))
  rownames(res_mat) <- make.unique(rownames(res_mat), sep = "_")  # Make row names unique
  return(res_mat)
}

# Step 3: Apply the conversion to every layer
bulk_transcriptomics_mat_validation <- convert_to_matrix(bulk_transcriptomics_df_normalized_vst_high_variance_validation,
                                                         all_samples_validation)
proteomics_mat_validation <- convert_to_matrix(proteomics_df_log_validation, all_samples_validation)
flow_abs_mat_validation <- convert_to_matrix(flow_abs_df_log_validation, all_samples_validation)
methylation_mat_validation <- convert_to_matrix(methylation_df_high_variance_validation, all_samples_validation)
ex_vivo_mat_validation <- convert_to_matrix(ex_vivo_df_log_validation, all_samples_validation)
viral_reservoir_mat_validation <- convert_to_matrix(viral_reservoir_df_log_validation, all_samples_validation)

# Bulk transcriptomics features are already labeled as ENSEMBLID:GENE_SYMBOL
# by QC_and_preprocessing.R, so no additional renaming is needed here.

# Filter to keep only the selected reservoir variables
viral_reservoir_mat_validation <- viral_reservoir_mat_validation[c("Total_million_avg", "intactDT_DSI"), ]

# Step 4: Create a nested list of matrices
data_omics_validation <- list(
  Bulk_transcriptomics = bulk_transcriptomics_mat_validation,
  Proteomics = proteomics_mat_validation,
  Flowcytometry_abs = flow_abs_mat_validation,
  Methylation = methylation_mat_validation,
  Ex_Vivo = ex_vivo_mat_validation,
  Viral_Reservoir = viral_reservoir_mat_validation
)

summary(data_omics_validation)
saveRDS(data_omics_validation, file.path(output_dir, "clustering/data_omics_validation_preimputed.rds"))
gc()

# Impute remaining missing values (from the union-of-samples step above)
# via KNN. Columns: i.e. samples: with more than 80% missingness in a
# given layer are dropped first, since imputing an almost-entirely-missing
# sample would mostly just reproduce the neighbor average.
impute_matrix <- function(mat, missing_threshold = 0.8) {
  # Calculate the proportion of missing values in each column
  missing_proportion <- colMeans(is.na(mat))
  
  # Identify columns to keep (those with less than the threshold of missing values)
  cols_to_keep <- missing_proportion < missing_threshold
  
  # Filter the matrix to keep only the columns with acceptable missing values
  mat_filtered <- mat[, cols_to_keep, drop = FALSE]
  
  # Apply KNN imputation to the filtered matrix
  imputed_data <- impute.knn(as.matrix(mat_filtered))$data
  
  # Return the imputed matrix
  return(imputed_data)
}

# Apply KNN imputation to each dataset in the data_omics list
data_omics_imputed_discovery <- lapply(data_omics_discovery, impute_matrix)
data_omics_imputed_validation <- lapply(data_omics_validation, impute_matrix)


# -----------------------------------------------------------------------------
# Restrict to samples present in every layer
# -----------------------------------------------------------------------------
# Some clustering algorithms (e.g. CIMLR) require a complete data matrix
# with no missing samples across layers, so: after imputing within each
# layer above: samples still missing an entire layer are dropped here.

# Step 1: Identify overlapping samples (column names) across all datasets
overlapping_samples_discovery <- Reduce(intersect, lapply(data_omics_imputed_discovery, colnames))
overlapping_samples_validation <- Reduce(intersect, lapply(data_omics_imputed_validation, colnames))

# Step 2: Filter each dataset to keep only the overlapping samples
data_omics_filtered_discovery <- lapply(data_omics_imputed_discovery, function(mat_discovery) {
  mat_discovery[, overlapping_samples_discovery, drop = FALSE]
})
data_omics_filtered_validation <- lapply(data_omics_imputed_validation, function(mat_validation) {
  mat_validation[, overlapping_samples_validation, drop = FALSE]
})

# Verify the filtering
lapply(data_omics_filtered_discovery, function(mat_discovery) head(mat_discovery[, 1:5]))
lapply(data_omics_filtered_validation, function(mat_validation) head(mat_validation[, 1:5]))

saveRDS(data_omics_filtered_discovery, 
        file.path(output_dir, "clustering/data_omics_discovery_filtered.rds"))
saveRDS(data_omics_filtered_validation, file.path(output_dir, "clustering/data_omics_validation_filtered.rds"))


# -----------------------------------------------------------------------------
# Clean feature and layer names for CIMLR compatibility
# -----------------------------------------------------------------------------
# CIMLR does not handle special characters in feature names, so underscores,
# "+" and "-" (common in flow cytometry marker names, e.g. "CD4+_CD8-") are
# replaced. Layer names in the list are also switched from underscores to
# dots to avoid conflicts once features are later addressed by dataset name.

# Get the current column names
current_names_flow_discovery <- rownames(data_omics_filtered_discovery$Flowcytometry_abs)
current_names_proteomics_discovery <- rownames(data_omics_filtered_discovery$Proteomics)
current_names_exvivo_discovery <- rownames(data_omics_filtered_discovery$Ex_Vivo)
current_names_reservoir_discovery <- rownames(data_omics_filtered_discovery$Viral_Reservoir)
current_names_methylation_discovery <- rownames(data_omics_filtered_discovery$Methylation)

# Function to replace all underscores with commas
replace_underscores <- function(name) {
  # Replace all underscores with commas
  new_name <- gsub("_", ".", name)
  new_name <- gsub("\\+", "pos", new_name)
  new_name <- gsub("-", "neg", new_name)
  return(new_name)
}

# Apply the function to all names
new_names_flow <- sapply(current_names_flow_discovery, replace_underscores)
new_names_proteomics <- sapply(current_names_proteomics_discovery, replace_underscores)
new_names_exvivo <- sapply(current_names_exvivo_discovery, replace_underscores)
new_names_reservoir <- sapply(current_names_reservoir_discovery, replace_underscores)
new_names_methylation <- sapply(current_names_methylation_discovery, replace_underscores)

# Replace the column names
rownames(data_omics_filtered_discovery$Flowcytometry_abs) <- new_names_flow
rownames(data_omics_filtered_discovery$Proteomics) <- new_names_proteomics
rownames(data_omics_filtered_discovery$Ex_Vivo) <- new_names_exvivo
rownames(data_omics_filtered_discovery$Viral_Reservoir) <- new_names_reservoir
rownames(data_omics_filtered_discovery$Methylation) <- new_names_methylation

# Print the new names to verify
data_omics_filtered_discovery$Flowcytometry_abs[1:5, 1:5]
data_omics_filtered_discovery$Proteomics[1:5,1:5]
data_omics_filtered_discovery$Ex_Vivo[1:5,1:5]
data_omics_filtered_discovery$Viral_Reservoir[1:2,1:5]
data_omics_filtered_discovery$Methylation[1:2,1:5]

# Replace underscores in dataset names to avoid downstream naming conflicts
# Remove underscores completely
names(data_omics_filtered_discovery) <- gsub("_", ".", names(data_omics_filtered_discovery))

# Print the new names to verify
print(names(data_omics_filtered_discovery))

saveRDS(data_omics_filtered_discovery, file.path(output_dir, "clustering/data_omics_discovery_filtered_renamed.rds"))


# Same name cleaning as above, repeated for the validation cohort
# Get the current column names
current_names_flow_validation <- rownames(data_omics_filtered_validation$Flowcytometry_abs)
current_names_proteomics_validation <- rownames(data_omics_filtered_validation$Proteomics)
current_names_exvivo_validation <- rownames(data_omics_filtered_validation$Ex_Vivo)
current_names_reservoir_validation <- rownames(data_omics_filtered_validation$Viral_Reservoir)
current_names_methylation_validation <- rownames(data_omics_filtered_validation$Methylation)

# Function to replace all underscores with commas
replace_underscores <- function(name) {
  # Replace all underscores with commas
  new_name <- gsub("_", ".", name)
  new_name <- gsub("\\+", "pos", new_name)
  new_name <- gsub("-", "neg", new_name)
  return(new_name)
}

# Apply the function to all names
new_names_flow <- sapply(current_names_flow_validation, replace_underscores)
new_names_proteomics <- sapply(current_names_proteomics_validation, replace_underscores)
new_names_exvivo <- sapply(current_names_exvivo_validation, replace_underscores)
new_names_reservoir <- sapply(current_names_reservoir_validation, replace_underscores)
new_names_methylation <- sapply(current_names_methylation_validation, replace_underscores)

# Replace the column names
rownames(data_omics_filtered_validation$Flowcytometry_abs) <- new_names_flow
rownames(data_omics_filtered_validation$Proteomics) <- new_names_proteomics
rownames(data_omics_filtered_validation$Ex_Vivo) <- new_names_exvivo
rownames(data_omics_filtered_validation$Viral_Reservoir) <- new_names_reservoir
rownames(data_omics_filtered_validation$Methylation) <- new_names_methylation

# Print the new names to verify
data_omics_filtered_validation$Flowcytometry_abs[1:5, 1:5]
data_omics_filtered_validation$Proteomics[1:5,1:5]
data_omics_filtered_validation$Ex_Vivo[1:5,1:5]
data_omics_filtered_validation$Viral_Reservoir[1:2,1:5]
data_omics_filtered_validation$Methylation[1:2,1:5]

# Replace underscores in dataset names to avoid downstream naming conflicts
# Remove underscores completely
names(data_omics_filtered_validation) <- gsub("_", ".", names(data_omics_filtered_validation))

# Print the new names to verify
print(names(data_omics_filtered_validation))

saveRDS(data_omics_filtered_validation, file.path(output_dir, "clustering/data_omics_validation_filtered_renamed.rds"))


# -----------------------------------------------------------------------------
# 2. Run Clustering Algorithms (GET Module)
# -----------------------------------------------------------------------------

# Identify optimal clustering number (may take a while)
optk.brca <- getClustNum(
  data        = data_omics_filtered_discovery,
  is.binary   = c(F, F, F, F, F, F),
  try.N.clust = 2:8,
  fig.path = file.path(output_dir, "clustering/"),
  # try cluster number from 2 to 8
  fig.name    = "cluster_number_assessment_discovery.pdf" 
)


# -----------------------------------------------------------------------------
# 2.1 Get Results from Single Algorithms
# -----------------------------------------------------------------------------

clustering_methods <- list("SNF", "PINSPlus", "NEMO", "COCA", "ConsensusClustering", "CIMLR", "MoCluster")

omics_type <- rep("gaussian", 6)

moic.res.list_discovery <- getMOIC(
  data        = data_omics_filtered_discovery,
  methodslist = clustering_methods,
  N.clust     = 3,
  type        = omics_type
)
saveRDS(moic.res.list_discovery, file = file.path(output_dir, "clustering/moic_results_discovery.rds"))

moic.res.list_validation <- getMOIC(
  data        = data_omics_filtered_validation,
  methodslist = clustering_methods,
  N.clust     = 3,
  type        = omics_type
)
saveRDS(moic.res.list_validation, file = file.path(output_dir, "clustering/moic_results_validation.rds"))


# -----------------------------------------------------------------------------
# 2.2 Get Consensus from Different Algorithms
# -----------------------------------------------------------------------------

cmoic.2000HIV_all_discovery <- getConsensusMOIC(
  moic.res.list = moic.res.list_discovery,
  fig.name      = "Consensus Heatmap (Discovery)",
  distance      = "euclidean",
  linkage       = "average"
)
saveRDS(cmoic.2000HIV_all_discovery, file = file.path(output_dir, "clustering/consensus_clusters_discovery.rds"))

cmoic.2000HIV_all_validation <- getConsensusMOIC(
  moic.res.list = moic.res.list_validation,
  fig.name      = "Consensus Heatmap (Validation)",
  distance      = "euclidean",
  linkage       = "average"
)
saveRDS(cmoic.2000HIV_all_validation, file = file.path(output_dir, "clustering/consensus_clusters_validation.rds"))


# -----------------------------------------------------------------------------
# 2.3 Assess Clustering Quality (Silhouette)
# -----------------------------------------------------------------------------

getSilhouette(sil      = cmoic.2000HIV_all_discovery$sil,
              fig.path = file.path(output_dir, "clustering/"),
              fig.name = "Silhouette - Consensus Clustering (Discovery)",
              height   = 5.5,
              width    = 5)

getSilhouette(sil      = cmoic.2000HIV_all_validation$sil,
              fig.path = file.path(output_dir, "clustering/"),
              fig.name = "Silhouette - Consensus Clustering (Validation)",
              height   = 5.5,
              width    = 5)


# -----------------------------------------------------------------------------
# 3. Visualize Clustering Solutions (Comprehensive Heatmaps)
# -----------------------------------------------------------------------------

# Data normalization for heatmap display
plotdata_discovery <- getStdiz(data       = data_omics_filtered_discovery,
                               halfwidth  = c(2,2,2,2,2,2), # no truncation for mutation
                               centerFlag = c(T,T,T,T,T,T), # no center for mutation
                               scaleFlag  = c(T,T,T,T,T,T)) # no scale for mutation

plotdata_validation <- getStdiz(data       = data_omics_filtered_validation,
                                halfwidth  = c(2,2,2,2,2,2), # no truncation for mutation
                                centerFlag = c(T,T,T,T,T,T), # no center for mutation
                                scaleFlag  = c(T,T,T,T,T,T)) # no scale for mutation

# Shared display settings for all six omics layers (used by every heatmap below)
row_titles <- c("Bulk_transcriptomics", "Proteomics", "Flowcytometry_abs",
                "Methylation", "Ex_Vivo", "Viral Reservoir")
legend_names <- c("bulkRNA-seq (log2)", "NPX (log2)",
                  "Flow Cytometry (Absolute, log2)", "DNA Methylation (M-values)",
                  "Ex-vivo Cytokine (log2)", "Viral Reservoir (log2)")
col_list <- list(
  c("#00FF00", "#008000", "#000000", "#800000", "#FF0000"),  # Bulk transcriptomics
  c("#1E90FF", "#87CEFA", "#FF4500"),                          # Proteomics
  c("#2C7BB6", "#FFFFBF", "#D7191C"),                          # Flow cytometry
  c("#0000FF", "#FFFFFF", "#FF0000"),                          # Methylation
  c("#0000FF", "#FFFFFF", "#FF0000"),                          # Ex vivo
  c("#0000FF", "#FFFFFF", "#FF0000")                           # Viral reservoir
)

# Reusable wrapper around getMoHeatmap(): all the plain per-algorithm
# comparisons below share identical row/legend/color settings and differ
# only in cluster assignment, feature annotation, and figure name.
plot_comprehensive_heatmap <- function(clust_res, data, annRow = NULL, fig_name) {
  getMoHeatmap(
    data          = data,
    row.title     = row_titles,
    is.binary     = rep(FALSE, 6),
    legend.name   = legend_names,
    clust.res     = clust_res,
    clust.dend    = NULL,
    show.rownames = rep(FALSE, 6),
    show.colnames = FALSE,
    annRow        = annRow,
    color         = col_list,
    annCol        = NULL,
    annColors     = NULL,
    width         = 10,
    height        = 5,
    fig.name      = fig_name
  )
}

# Helper: top-N features by absolute loading/importance per omics dataset,
# used to annotate rows for CIMLR and MoCluster heatmaps.
get_key_markers <- function(feat_res, dataset_names, n = 10) {
  lapply(dataset_names, function(ds) feat_res[feat_res$dataset == ds, "feature"][1:n])
}

# -----------------------------------------------------------------------------
# 3.1 SNF, PINSPlus, NEMO, COCA, ConsensusClustering, CIMLR
# -----------------------------------------------------------------------------
# Plain comprehensive heatmap per algorithm, discovery and validation cohorts.

annRow_CIMLR_discovery  <- get_key_markers(moic.res.list_discovery$CIMLR$feat.res,
                                           c("Bulk_transcriptomics", "Proteomics", "Flowcytometry_abs",
                                             "Methylation", "Ex_Vivo", "Viral_Reservoir"))
annRow_CIMLR_validation <- get_key_markers(moic.res.list_validation$CIMLR$feat.res,
                                           c("Bulk_transcriptomics", "Proteomics", "Flowcytometry_abs",
                                             "Methylation", "Ex_Vivo", "Viral_Reservoir"))

simple_algorithms <- c("SNF", "PINSPlus", "NEMO", "COCA", "ConsensusClustering", "CIMLR")

for (algo in simple_algorithms) {
  annRow_disc <- if (algo == "CIMLR") annRow_CIMLR_discovery  else NULL
  annRow_val  <- if (algo == "CIMLR") annRow_CIMLR_validation else NULL
  
  plot_comprehensive_heatmap(
    clust_res = moic.res.list_discovery[[algo]]$clust.res,
    data      = plotdata_discovery,
    annRow    = annRow_disc,
    fig_name  = paste0("COMPREHENSIVE HEATMAP OF ", algo, " (Discovery)")
  )
  plot_comprehensive_heatmap(
    clust_res = moic.res.list_validation[[algo]]$clust.res,
    data      = plotdata_validation,
    annRow    = annRow_val,
    fig_name  = paste0("COMPREHENSIVE HEATMAP OF ", algo, " (Validation)")
  )
}

# -----------------------------------------------------------------------------
# 3.2 MoCluster: Driver Marker Comparisons (Discovery vs. Validation)
# -----------------------------------------------------------------------------
# MoCluster was the final selected algorithm, so its solution is examined in
# more depth than the others: comparing feature loadings between cohorts and
# restricting to markers validated in both. All six layers are shown in
# every heatmap below (row.title/legend.name/color come from the shared
# 6-layer constants used by plot_comprehensive_heatmap()); what differs
# between the four heatmaps is only which features are annotated as driver
# markers (annRow), not which layers are displayed.

# Top 15 driver markers per layer from MoCluster's own feature loadings.
# Dataset name order matches row_titles/legend_names/col_list exactly.
mocluster_dataset_names <- c("Bulk.transcriptomics", "Proteomics", "Flowcytometry.abs",
                             "Methylation", "Ex.Vivo", "Viral.Reservoir")

annRow_MoCluster_discovery <- get_key_markers(MoCluster.res_discovery$feat.res, mocluster_dataset_names, n = 15)
annRow_MoCluster_validation <- get_key_markers(MoCluster.res_validation$feat.res, mocluster_dataset_names, n = 15)

# --- Discovery, MoCluster's own top driver markers per layer ---
plot_comprehensive_heatmap(
  clust_res = MoCluster.res_discovery$clust.res,
  data      = plotdata_discovery,
  annRow    = annRow_MoCluster_discovery,
  fig_name  = "COMPREHENSIVE HEATMAP OF MOCluster - Discovery (Driver Markers)"
)

# --- Discovery, restricted to markers validated in both cohorts ---
# For bulk RNA-seq, proteomics, ex-vivo cytokines, and flow cytometry, only features selected
# as drivers in BOTH the discovery and validation MoCluster solutions are
# annotated. 
markers_validated <- intersect(MoCluster.res_discovery$feat.res$feature, MoCluster.res_validation$feat.res$feature)

plotdata_discovery_temp <- plotdata_discovery
plotdata_discovery_temp$Bulk.transcriptomics <- plotdata_discovery_temp$Bulk.transcriptomics[rownames(plotdata_discovery_temp$Bulk.transcriptomics) %in% markers_validated, ]
plotdata_discovery_temp$Proteomics           <- plotdata_discovery_temp$Proteomics[rownames(plotdata_discovery_temp$Proteomics) %in% markers_validated, ]
plotdata_discovery_temp$Flowcytometry.abs    <- plotdata_discovery_temp$Flowcytometry.abs[rownames(plotdata_discovery_temp$Flowcytometry.abs) %in% markers_validated, ]
plotdata_discovery_temp$Ex.Vivo    <- plotdata_discovery_temp$Ex.Vivo[rownames(plotdata_discovery_temp$Ex.Vivo) %in% markers_validated, ]
plotdata_discovery_temp$Methylation    <- plotdata_discovery_temp$Methylation[rownames(plotdata_discovery_temp$Methylation) %in% markers_validated, ]

annRow_MoCluster_discovery_temp <- annRow_MoCluster_discovery
annRow_MoCluster_discovery_temp[[1]] <- rownames(plotdata_discovery_temp$Bulk.transcriptomics)
annRow_MoCluster_discovery_temp[[2]] <- rownames(plotdata_discovery_temp$Proteomics)
annRow_MoCluster_discovery_temp[[3]] <- rownames(plotdata_discovery_temp$Flowcytometry.abs)
annRow_MoCluster_discovery_temp[[4]] <- rownames(plotdata_discovery_temp$Ex.Vivo)
annRow_MoCluster_discovery_temp[[5]] <- rownames(plotdata_discovery_temp$Methylation)

plot_comprehensive_heatmap(
  clust_res = MoCluster.res_discovery$clust.res,
  data      = plotdata_discovery_temp,
  annRow    = annRow_MoCluster_discovery_temp,
  fig_name  = "COMPREHENSIVE HEATMAP OF MOCluster - Discovery (Markers Validated in Both Cohorts)"
)

# --- Validation, restricted to markers validated in both cohorts ---
# Cluster labels 1 and 3 are swapped in the validation solution so that
# cluster numbering is consistent with the discovery cohort.
MoCluster.res_validation_temp <- MoCluster.res_validation$clust.res
MoCluster.res_validation_temp$clust <- ifelse(
  MoCluster.res_validation_temp$clust == 1, 3,
  ifelse(MoCluster.res_validation_temp$clust == 3, 1, MoCluster.res_validation_temp$clust)
)

plotdata_validation_temp <- plotdata_validation
plotdata_validation_temp$Bulk.transcriptomics <- plotdata_validation_temp$Bulk.transcriptomics[rownames(plotdata_validation_temp$Bulk.transcriptomics) %in% markers_validated, ]
plotdata_validation_temp$Proteomics           <- plotdata_validation_temp$Proteomics[rownames(plotdata_validation_temp$Proteomics) %in% markers_validated, ]
plotdata_validation_temp$Flowcytometry.abs    <- plotdata_validation_temp$Flowcytometry.abs[rownames(plotdata_validation_temp$Flowcytometry.abs) %in% markers_validated, ]
plotdata_validation_temp$Ex.Vivo    <- plotdata_validation_temp$Ex.Vivo[rownames(plotdata_validation_temp$Ex.Vivo) %in% markers_validated, ]
plotdata_validation_temp$Methylation    <- plotdata_validation_temp$Methylation[rownames(plotdata_validation_temp$Methylation) %in% markers_validated, ]

annRow_MoCluster_validation_temp <- annRow_MoCluster_validation
annRow_MoCluster_validation_temp[[1]] <- rownames(plotdata_validation_temp$Bulk.transcriptomics)
annRow_MoCluster_validation_temp[[2]] <- rownames(plotdata_validation_temp$Proteomics)
annRow_MoCluster_validation_temp[[3]] <- rownames(plotdata_validation_temp$Flowcytometry.abs)
annRow_MoCluster_validation_temp[[4]] <- rownames(plotdata_validation_temp$Ex.Vivo)
annRow_MoCluster_validation_temp[[5]] <- rownames(plotdata_validation_temp$Methylation)

plot_comprehensive_heatmap(
  clust_res = MoCluster.res_validation_temp,
  data      = plotdata_validation_temp,
  annRow    = annRow_MoCluster_validation_temp,
  fig_name  = "COMPREHENSIVE HEATMAP OF MOCluster - Validation (Markers Validated in Both Cohorts)"
)

# --- Validation, restricted to the validation cohort's own driver markers ---
plotdata_validation_temp2 <- plotdata_validation
plotdata_validation_temp2$Bulk.transcriptomics <- plotdata_validation_temp2$Bulk.transcriptomics[rownames(plotdata_validation_temp2$Bulk.transcriptomics) %in% MoCluster.res_validation$feat.res$feature, ]
plotdata_validation_temp2$Proteomics    <- plotdata_validation_temp2$Proteomics[rownames(plotdata_validation_temp2$Proteomics) %in% MoCluster.res_validation$feat.res$feature, ]
plotdata_validation_temp2$Flowcytometry.abs    <- plotdata_validation_temp2$Flowcytometry.abs[rownames(plotdata_validation_temp2$Flowcytometry.abs) %in% MoCluster.res_validation$feat.res$feature, ]
plotdata_validation_temp2$Ex.Vivo    <- plotdata_validation_temp2$Ex.Vivo[rownames(plotdata_validation_temp2$Ex.Vivo) %in% MoCluster.res_validation$feat.res$feature, ]
plotdata_validation_temp2$Methylation    <- plotdata_validation_temp2$Methylation[rownames(plotdata_validation_temp2$Methylation) %in% MoCluster.res_validation$feat.res$feature, ]

plot_comprehensive_heatmap(
  clust_res = MoCluster.res_validation_temp,
  data      = plotdata_validation_temp2,
  annRow    = annRow_MoCluster_validation,
  fig_name  = "COMPREHENSIVE HEATMAP OF MOCluster - Validation (Own Driver Markers)"
)


# -----------------------------------------------------------------------------
# 4. Compare Clustering Algorithms Against Reservoir/Clinical Variables
# -----------------------------------------------------------------------------
# Evaluates how well each algorithm's clusters discriminate the viral
# reservoir variables (total and intact HIV-1 DNA), run separately on the
# discovery and validation cohorts so the algorithms' consistency across
# cohorts can be compared. This comparison is what determined the final
# algorithm choice reported in the manuscript.

viral_reservoir_clean$ID <- viral_reservoir_clean$SampleID
sample_table_viral_reservoir <- merge(sample_table, viral_reservoir_clean, by = "ID")
rownames(sample_table_viral_reservoir) <- sample_table_viral_reservoir$ID

nonnormal_vars <- c("Total_million (Average)", "IPDA_million_DSI", "intactDT_DSI")
factor_vars    <- c("Total_million (Average) Quantiles", "IPDA_million_DSI_quantile", "intactDT_DSI_quantile")

# Statistical test per variable: Kruskal-Wallis for continuous, chi-square
# for categorical, run against a given clustering assignment
test_variable_vs_cluster <- function(data, var, cluster_column, test_type) {
  result <- tryCatch({
    if (test_type == "continuous") {
      formula <- as.formula(paste0("`", var, "` ~ `", cluster_column, "`"))
      test <- kruskal_test(formula, data = data)
      list(test = "Kruskal-Wallis", statistic = test$statistic, p.value = test$p)
    } else {
      tbl <- table(data[[var]], data[[cluster_column]])
      test <- chisq.test(tbl)
      list(test = "Chi-square", statistic = test$statistic, p.value = test$p.value)
    }
  }, error = function(e) list(test = ifelse(test_type == "continuous", "Kruskal-Wallis", "Chi-square"), statistic = NA, p.value = NA))
  
  data.frame(variable = var, test = result$test, statistic = result$statistic, p.value = result$p.value)
}

# Runs the full comparison (clinical table build, formal MOVICS comparison,
# statistical tests, and ranking) for one cohort's clustering results
evaluate_clustering_vs_reservoir <- function(moic_res_list, cohort_label) {
  
  clustering_results <- list(
    "Cluster_SNF"       = moic_res_list$SNF,
    "Cluster_PINSPlus"  = moic_res_list$PINSPlus,
    "Cluster_NEMO"      = moic_res_list$NEMO,
    "Cluster_COCA"      = moic_res_list$COCA,
    "Cluster_Consensus" = moic_res_list$ConsensusClustering,
    "Cluster_CIMLR"     = moic_res_list$CIMLR,
    "Cluster_MoCluster" = moic_res_list$MoCluster
  )
  
  # Restrict to samples with omics data, in the algorithms' sample order
  sample_ids <- moic_res_list$ConsensusClustering$clust.res$samID
  clin_table <- sample_table_viral_reservoir[sample_ids, ]
  
  for (clust_name in names(clustering_results)) {
    clin_table[[clust_name]] <- as.factor(clustering_results[[clust_name]]$clust.res$clust)
  }
  
  clin_table <- clin_table[, c("Total_million_avg", "IPDA_million_DSI", "intactDT_DSI", names(clustering_results))]
  
  clin_table <- clin_table %>%
    mutate(
      `Total_million (Average)` = as.numeric(Total_million_avg),
      `Total_million (Average) Quantiles` = as.factor(ntile(as.numeric(Total_million_avg), 4)),
      IPDA_million_DSI = as.numeric(IPDA_million_DSI),
      IPDA_million_DSI_quantile = as.factor(ntile(as.numeric(IPDA_million_DSI), 4)),
      intactDT_DSI = as.numeric(intactDT_DSI),
      intactDT_DSI_quantile = as.factor(ntile(as.numeric(intactDT_DSI), 4)),
      across(all_of(names(clustering_results)), factor)
    )
  
  # Formal MOVICS clinical-variable comparison table per algorithm (Word output)
  for (clust_name in names(clustering_results)) {
    compClinvar(
      moic.res      = clustering_results[[clust_name]],
      var2comp      = clin_table,
      strata        = clustering_results[[clust_name]]$clust.res,
      factorVars    = factor_vars,
      nonnormalVars = nonnormal_vars,
      exactVars     = NULL,
      doWord        = TRUE,
      tab.name      = paste0("Clinical_Features_Comparison_", clust_name, "_", cohort_label)
    )
  }
  
  # Kruskal-Wallis / chi-square tests per algorithm x variable
  all_results <- bind_rows(lapply(names(clustering_results), function(clust_name) {
    continuous_tests <- lapply(nonnormal_vars, test_variable_vs_cluster, data = clin_table, cluster_column = clust_name, test_type = "continuous")
    categorical_tests <- lapply(factor_vars, test_variable_vs_cluster, data = clin_table, cluster_column = clust_name, test_type = "categorical")
    results <- bind_rows(continuous_tests, categorical_tests)
    results$clustering_method <- clust_name
    results
  })) %>%
    select(clustering_method, variable, test, statistic, p.value)
  
  # Heatmap of -log10(p-values) across algorithms and variables
  heatmap_data <- all_results %>%
    select(clustering_method, variable, p.value) %>%
    pivot_wider(names_from = clustering_method, values_from = p.value) %>%
    pivot_longer(cols = -variable, names_to = "clustering_method", values_to = "p_value")
  
  max_log_p <- ceiling(max(-log10(heatmap_data$p_value), na.rm = TRUE))
  
  p_heatmap <- ggplot(heatmap_data, aes(x = variable, y = clustering_method, fill = -log10(p_value))) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = -log10(0.05), limit = c(0, max_log_p),
                         space = "Lab", name = "-Log10 p-value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Variables", y = "Clustering Methods",
         title = paste0("Heatmap of -Log10 p-values for Viral Reservoir Variables (", cohort_label, ")"))
  ggsave(file.path(output_dir, paste0("clustering/algorithm_comparison_heatmap_", cohort_label, ".pdf")), plot = p_heatmap, width = 10, height = 6)
  
  # Bar plot of significant-variable counts per algorithm
  significant_vars <- all_results %>%
    group_by(clustering_method) %>%
    summarize(significant_count = sum(p.value < 0.05, na.rm = TRUE))
  
  p_barplot <- ggplot(significant_vars, aes(x = reorder(clustering_method, -significant_count), y = significant_count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Clustering Method", y = "Number of Significant Variables",
         title = paste0("Significant Viral Reservoir Variables by Clustering Method (", cohort_label, ")"))
  ggsave(file.path(output_dir, paste0("clustering/algorithm_comparison_barplot_", cohort_label, ".pdf")), plot = p_barplot, width = 10, height = 6)
  
  # Rank algorithms by overall discrimination of the reservoir variables
  method_ranking <- all_results %>%
    group_by(clustering_method) %>%
    summarize(
      mean_log10_p = mean(-log10(p.value), na.rm = TRUE),
      significant_count = sum(p.value < 0.05, na.rm = TRUE),
      total_count = n()
    ) %>%
    mutate(
      significant_proportion = significant_count / total_count,
      overall_score = (mean_log10_p + significant_proportion) / 2
    ) %>%
    arrange(desc(overall_score))
  
  list(results = all_results, ranking = method_ranking)
}

evaluation_discovery  <- evaluate_clustering_vs_reservoir(moic.res.list_discovery, "Discovery")
evaluation_validation <- evaluate_clustering_vs_reservoir(moic.res.list_validation, "Validation")

cat("\nAlgorithm ranking: Discovery:\n")
print(evaluation_discovery$ranking)

cat("\nAlgorithm ranking: Validation:\n")
print(evaluation_validation$ranking)

# Consistency check: does the same algorithm rank well in both cohorts?
ranking_comparison <- merge(
  evaluation_discovery$ranking[, c("clustering_method", "overall_score")],
  evaluation_validation$ranking[, c("clustering_method", "overall_score")],
  by = "clustering_method",
  suffixes = c("_discovery", "_validation")
) %>% arrange(desc(overall_score_discovery + overall_score_validation))

cat("\nAlgorithm ranking consistency (Discovery vs. Validation):\n")
print(ranking_comparison)
