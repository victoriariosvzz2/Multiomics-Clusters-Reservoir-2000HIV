#!/usr/bin/env Rscript
# =============================================================================
# Script: deg_analysis.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Differential expression analysis (bulk RNA-seq) between the
#              three MoCluster immune-reservoir endotypes (Mixed, AllLow, AllHigh),
#              including exploratory PCA/UMAP diagnostics, LIMMA-based batch
#              effect correction (season, plate, center, time-to-lab,
#              genetic PC1, age, sex, COVID vaccination status), DESeq2
#              differential testing, and summary visualizations (volcano,
#              MA, upset, Venn, hierarchical clustering).
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Run independently once per cohort: set `cohort` below to "discovery" or
# "validation" and re-run. Continues in: gsea_analysis.R (discovery only).
# =============================================================================

# ------------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------------

# Set paths — update these to match your local directory structure
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

# Which cohort to run this script on — "discovery" or "validation"
cohort <- "discovery"

knitr::opts_chunk$set(dev = "png",
                      dpi = 300,
                      echo = FALSE,
                      cache = TRUE)


# ------------------------------------------------------------------------------------
# 1. Install and Load Packages
# ------------------------------------------------------------------------------------

# CRAN packages
list.of.packages <- c(
  "apeglm",
  "factoextra",
  "ggbeeswarm",
  "ggplot2",
  "ggrepel",
  "gplots",
  "hexbin",
  "Hmisc",
  "openxlsx",
  "patchwork",
  "plotly",
  "reshape2",
  "scales",
  "tidyr",
  "VennDiagram",
  "UpSetR",
  "rstatix",
  "dplyr",
  "ggExtra",
  "readr"
)

new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0)
  install.packages(new.packages)

# Bioconductor packages (GSVA/GSEABase/DOSE/clusterProfiler/IHW are not
# needed here — they're only used by the GSEA/GSVA-specific functions in
# Functions.R, none of which this script calls; see gsea_analysis.R)
list.of.bioc.packages <- c(
  "biomaRt",
  "ComplexHeatmap",
  "DESeq2",
  "genefilter",
  "ggpubr",
  "limma",
  "org.Hs.eg.db",
  "pheatmap",
  "RColorBrewer",
  "rhdf5",
  "sva",
  "tximport",
  "vsn",
  "ggnewscale",
  "grid" ,
  "ggExtra",
  "reshape2"
)

new.packages.bioc <-
  list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[, "Package"])]
if (length(new.packages.bioc) > 0)
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(new.packages.bioc, update = FALSE)

# chunk options: load packages, results='hide',message=FALSE,warning=FALSE
lapply(c(list.of.packages, list.of.bioc.packages),
       require,
       character.only = TRUE)

rm(list.of.packages,
   list.of.bioc.packages,
   new.packages,
   new.packages.bioc)


# -----------------------------------------------------------------------------
# 1.2 Load Functions
# -----------------------------------------------------------------------------

source(file.path(project_dir, "scripts/DESeq_functions_KD_RK_VR.R"))


# -----------------------------------------------------------------------------
# 2. Project Information
# -----------------------------------------------------------------------------

organism = "human" 

# chunk options: , warning=F
dir <-
  file.path(project_dir, "scripts/")

# save for later
dir_out <-
  file.path(output_dir, "differential_expression/")


# -----------------------------------------------------------------------------
# 3.1 Load Gene Annotation
# -----------------------------------------------------------------------------

# chunk options: gene annotation import
# Specify the filename of your gene annotation file here: 
annotation_filename <- "ID2SYMBOL_gencode_v27_transcript.txt"
annotation_dir <- file.path(data_dir, "reference/GMTfiles")

tx_annotation <-
  read.delim(
    file.path(
      paste0(annotation_dir),
      annotation_filename
    ),
    header = F ,
    stringsAsFactors = F,
    col.names = c("GENEID", "TXNAME", "SYMBOL", "GENETYPE")
  )


# -----------------------------------------------------------------------------
# 3.2 Load Sample Table and Cluster Assignments
# -----------------------------------------------------------------------------

# chunk options: sample table import
library(readxl)
library(dplyr)


# Transcriptomics sample table
sample_table <-
  readRDS(file = file.path(data_dir, "clinical/2000HIV_bulk_transcriptomics_sample_table.RDS"))

# Left join the data frames based on the "ID" column
sample_table$NIJMEGEN_DATE_COLLECTION_to_JAN01 <-
  as.Date(sample_table$NIJMEGEN_DATE_COLLECTION) - as.Date("2019-01-01")
sample_table$NIJMEGEN_SEASON_SCORE_SINE <-
  sin((2 * pi * (
    as.numeric(sample_table$NIJMEGEN_DATE_COLLECTION_to_JAN01)
  )) / 365.2425)
sample_table$NIJMEGEN_SEASON_SCORE_COSINE <-
  cos((2 * pi * (
    as.numeric(sample_table$NIJMEGEN_DATE_COLLECTION_to_JAN01)
  )) / 365.2425)


####
# Import the multiomics clusters
clustering_results_discovery <- readRDS(file.path(output_dir, paste0("clustering/consensus_clusters_", cohort, ".rds")))
clustering_results_validation <- readRDS(file.path(output_dir, "clustering/consensus_clusters_validation.rds"))


# Check the new group distribution
table(clustering_results_discovery$clust.res$clust)
table(clustering_results_validation$clust.res$clust)

clustering_results_MoCluster <- rbind(clustering_results_discovery$clust.res, clustering_results_validation$clust.res)

clustering_results_MoCluster <- clustering_results_MoCluster %>% mutate(ID = as.character(samID),
                                                                        MoCluster_Cluster = factor(clust))
dim(clustering_results_MoCluster)

clustering_results_MoCluster <- clustering_results_MoCluster %>%
  mutate(MoCluster_Cluster = recode(MoCluster_Cluster, 
                            "1" = "Mixed",
                            "2" = "AllLow",
                            "3" = "AllHigh"))

clustering_results_MoCluster$MoCluster_Cluster <- relevel(clustering_results_MoCluster$MoCluster_Cluster, ref = "Mixed")
clustering_results_MoCluster$MoCluster_Cluster  <- factor(clustering_results_MoCluster$MoCluster_Cluster,
                                                          levels = c("Mixed", "AllLow", "AllHigh"))

clustering_results_MoCluster$MoCluster_Cluster <- relevel(clustering_results_MoCluster$MoCluster_Cluster, ref = "Mixed")

table(clustering_results_MoCluster$MoCluster_Cluster)

# Let's just keep the samples that were part of the clustering step and merge the clinical dataset
sample_table <- right_join(sample_table, clustering_results_MoCluster, by = c("DONOR_ID" = "ID"))
dim(sample_table)


# Let's keep only the cohort of interest
sample_table <- sample_table[sample_table$COHORT == toupper(cohort),]

#### INSERT here your condition you want to inspect and the colors of your condition, default is clinical group
sample_table$condition <-
  sample_table$MoCluster_Cluster  #insert your condition of interest here

####
col_condition <-
  c("#44015480",
    "#21908c80",
    "#fde72580")
names(col_condition) <-
  c("Mixed",
    "AllLow",
    "AllHigh")

table(sample_table$condition, useNA = "ifany") #we first check how many have NA in the column of interest

# subset the sample_table so only the one with present information are there (of course any other filtering can also be done):
# this code is not executed else we would mess up the data if we want to include all
sample_table <- sample_table[!is.na(sample_table$condition),]

table(sample_table$condition, useNA = "ifany") 
table(sample_table$ETHNICITY, useNA = "ifany") 


## Add columns with factors for comparisons in model
# sample_table$condition <- factor(sample_table$condition,
#                                  levels = unique(sample_table$condition))
sample_table$AGE <- as.numeric(sample_table$AGE)
sample_table$CENTER_RUMC <- as.factor(sample_table$CENTER_RUMC)
sample_table$ETHNICITY <- as.factor(sample_table$ETHNICITY)

# define factor for order of samples in plotting
plot_order <- c("condition")

col_CENTER <-
  c(c(
    "darkblue",
    "orange",
    "indianred",
    "pink",
    "lightblue",
    "#CF8EF8"
  ))
names(col_CENTER) <-
  c("OLV", "EMC", "RUMC", "ETZ", "NA", "batch_control")

col_sex <- c("deepskyblue2", "coral")
names(col_sex) <- c("M", "F")

col_iso2 <- c("cornflowerblue", "firebrick")
names(col_iso2) <- c("before", "after")

col_cov <- c("blue", "red")
names(col_cov) <- c("negative", "past_positive")

col_GEN <- c("grey", "deepskyblue2", "coral", "red", "blue")
names(col_GEN) <- c("Batch_Control", "M", "F", "Trans-F", "Trans-M")

# combine color code into list
ann_colors <-
  list(
    condition = col_condition,
    SEX_BIRTH = col_sex,
    CENTER = col_CENTER,
    iso2 = col_iso2,
    cov = col_cov
  )

## ALL ===================
plot_obj <-
  ggplot(sample_table, aes(x = condition, fill = condition)) +
  geom_bar(stat = "count", color = "black") +
  stat_count(
    geom = "text",
    colour = "black",
    size = 3.5,
    aes(label = ..count..),
    position = position_stack(vjust = 0.5)
  ) +
  theme_bw()+
  scale_fill_manual(values = col_condition)

# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", paste0("Barplot_MoClusters_k3_", cohort, ".png")),
  plot = plot_obj,
  device = "png"
)

print(plot_obj)


# -----------------------------------------------------------------------------
# 3.3 Import Genetic PCs (PLINK)
# -----------------------------------------------------------------------------

# Write out the sample IDs used for this cohort
sample_ids <- c(sample_table$DONOR_ID)

# Function to process genetic pcs data
process_genetic_pcs <- function(file_path, sample_table, pcs_to_select = 5) {
  # Read the data
  genetic_pcs <- read.table(file_path, row.names = 1)
  
  # Match the IDs and row names
  genetic_pcs <- genetic_pcs[match(sample_table$DONOR_ID, as.character(genetic_pcs$V2)), ]
  genetic_pcs$V2 <- sample_table$DONOR_ID
  row.names(genetic_pcs) <- row.names(sample_table)
  
  # Replace the NAs with 0
  genetic_pcs[is.na(genetic_pcs)] <- 0
  
  # Choose the first 5 PCs
  pcs_to_select <- pcs_to_select + 2
  genetic_pcs_covariates <- genetic_pcs[, 2:pcs_to_select]
  
  # Verify the order
  head(cbind(genetic_pcs$V2, sample_table$DONOR_ID))
  head(genetic_pcs)
  head(genetic_pcs_covariates)
  
  return(genetic_pcs_covariates)
}

genetic_pcs <- process_genetic_pcs(
  file.path(data_dir, paste0("genetics/2000HIV_genetic_pcs_", cohort, ".pca.eigenvec")),
  sample_table
)


# Function to add principal components to sample table
add_pcs_to_sample_table <- function(sample_table, genetic_pcs_covariates) {
  sample_table$PC1 <- genetic_pcs_covariates[, 1]
  sample_table$PC2 <- genetic_pcs_covariates[, 2]
  sample_table$PC3 <- genetic_pcs_covariates[, 3]
  sample_table$PC4 <- genetic_pcs_covariates[, 4]
  sample_table$PC5 <- genetic_pcs_covariates[, 5]
  
  return(sample_table)
}

# Add PCs to All samples sample table
sample_table <- add_pcs_to_sample_table(sample_table, genetic_pcs)


# -----------------------------------------------------------------------------
# 4. Import Count Data
# -----------------------------------------------------------------------------

# Function to read RDS file and align column names
read_and_align <- function(file_path, sample_table) {
  data <- readRDS(file = file_path)
  
  if (!identical(colnames(data), as.character(sample_table$ID))) {
    print(
      "colnames of data and sample_table are not identical. Two options: either wrongly ordered or you subsetted the sample_table for inspection of a certain comparison"
    )
  }
  
  if (length(colnames(data)) == length(sample_table$ID)) {
    sample_table <- sample_table[, match(colnames(data), as.character(sample_table$ID))]
  } else {
    data <- data[, colnames(data) %in% sample_table$ID]
    sample_table <- sample_table[match(colnames(data), as.character(sample_table$ID)), ]
  }
  
  identical(colnames(data), as.character(sample_table$ID))
  
  return(list(data = data, sample_table = sample_table))
}

path_star_count <- file.path(data_dir, "omics/2000HIV_bulk_transcriptomics_raw_counts.RDS")

# Read and align star.count
result_star <- read_and_align(
  file_path = path_star_count,
  sample_table = sample_table
)

star.count <- result_star$data
sample_table <- result_star$sample_table


# -----------------------------------------------------------------------------
# 5.1 Generate DESeq2 Object
# -----------------------------------------------------------------------------

# chunk options: DESeqDataSetFromMatrix, cache = TRUE
# Function to create ETHNICITY_Black variable and DESeqDataSet
create_ethnicity_and_deseq <- function(sample_table, star_count) {
  # Create ETHNICITY_Black variable
  sample_table$ETHNICITY_Black <- ifelse(sample_table$ETHNICITY == "Black", 1, 0)
  sample_table$ETHNICITY_Black[is.na(sample_table$ETHNICITY_Black)] <- 0
  
  # Also let's make the CENTER_RUMC numeric
  sample_table$CENTER_RUMC <- ifelse(sample_table$CENTER == "RUMC", 1, 0)
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = star_count,
                                colData = sample_table,
                                design = ~ condition)
  
  return(dds)
}

# Create ETHNICITY_Black and DESeqDataSet for sample_table
dds_txi <- create_ethnicity_and_deseq(sample_table, star.count)


# Create dummy variables for the ethnicities
library(fastDummies)

sample_table <- fastDummies::dummy_cols(sample_table,
                                        select_columns = "ETHNICITY",
                                        remove_selected_columns = FALSE)


# -----------------------------------------------------------------------------
# 5.2 Pre-Filtering
# -----------------------------------------------------------------------------

# chunk options: pre-filtering
# Function to filter genes based on counts and conditions
filter_genes <- function(dds, sample_table) {
  # Filtering based on condition
  genes_to_keep <- rowSums(counts(dds) >= 5) >= table(sample_table$condition) %>% min()
  
  # Display tables
  print(table(genes_to_keep))
  
  # Filter dds
  dds_filtered <- dds[genes_to_keep,]
  
  return(list(dds_filtered = dds_filtered, genes_to_keep = genes_to_keep))
}

# Filter genes for dds_txi
dds <- filter_genes(dds_txi, sample_table)[["dds_filtered"]]

# Save the names of the genes we kept in the analysis to use the same ones in the validation cohort
dds_genes <- dds %>% rownames()

write_lines(dds_genes, file = file.path(dir_out, paste0("Tables/dds_genes_to_keep_", cohort, ".txt")))


# -----------------------------------------------------------------------------
# 5.3 DESeq2 Calculations
# -----------------------------------------------------------------------------
dds <- DESeq(dds, parallel = F) # All


# -----------------------------------------------------------------------------
# 5.4 Normalized Counts
# -----------------------------------------------------------------------------
library(biomaRt)

# Function to generate normalization annotation and extract required columns
generate_and_extract_norm_anno <- function(dds_object) {
  # Generate normalization annotation
  norm_anno_list <- generate_norm_anno(dds_object = dds_object)
  
  # Extract required columns from norm_anno
  norm_anno <- norm_anno_list$norm_anno
  
  # Extract gene_annotation
  gene_annotation <- norm_anno_list$gene_annotation
  
  return(list(norm_anno = norm_anno, gene_annotation = gene_annotation))
}

# Generate and extract normalization annotation for dds_txi
result_txi <- generate_and_extract_norm_anno(dds)
norm_anno <- result_txi$norm_anno
gene_annotation <- result_txi$gene_annotation


# Function to generate and save boxplot
generate_and_save_boxplot <- function(norm_anno, sample_df, dir_out, title) {
  # Generate boxplot
  plot_obj <- boxplot_norm(norm_anno = norm_anno, sample_table = sample_df)
  
  # Save plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", title),
    plot = plot_obj,
    device = "pdf"
  )
  
  # Print plot
  print(plot_obj)
}

# Generate and save boxplot for dds_txi
generate_and_save_boxplot(norm_anno, sample_table, dir_out, "Boxplot_Normalized_Samples_MoClusters_k3.pdf")

norm_mean <- mean_function(input = norm_anno,
                           anno = sample_table,
                           condition = "condition")


# -----------------------------------------------------------------------------
# 5.5 Variance Stabilizing Transformation (VST)
# -----------------------------------------------------------------------------
dds_vst <- vst(dds, blind = TRUE) # All

# Function to generate and save meanSdPlot
generate_meanSdPlot <- function(dds_vst, dir_out, filename) {
  plot_obj <- meanSdPlot(as.matrix(assay(dds_vst)), ranks = FALSE) + theme_bw()
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", filename),
    plot = plot_obj,
    device = "pdf"
  )
  
  print(plot_obj)
}

# Apply the function
generate_meanSdPlot(dds_vst, dir_out, paste0("Row_STD_vs_Row_Means_MoClusters_k3_", cohort, ".pdf")) # All

# Function to process dds_vst and norm_anno
process_vst_data <- function(dds_vst, norm_anno) {
  # Data frame of log values
  vst_anno_log <- as.data.frame(assay(dds_vst))
  vst_anno_log$GENEID <- rownames(vst_anno_log)
  vst_anno_log <- merge(
    vst_anno_log,
    norm_anno[, c("GENEID", "SYMBOL", "GENETYPE", "DESCRIPTION", "CHR", "start", "end")],
    by = "GENEID"
  )
  rownames(vst_anno_log) <- vst_anno_log$GENEID
  
  # Data frame of unlog values
  vst_anno <- as.data.frame(assay(dds_vst))
  vst_anno <- 2 ^ vst_anno
  vst_anno$GENEID <- rownames(vst_anno)
  vst_anno <- merge(
    vst_anno,
    norm_anno[, c("GENEID", "SYMBOL", "GENETYPE", "DESCRIPTION", "CHR", "start", "end")],
    by = "GENEID"
  )
  rownames(vst_anno) <- vst_anno$GENEID
  
  # Print selected rows from both data frames
  print(vst_anno_log[1:3, c(1:2, (ncol(vst_anno_log) - 5):ncol(vst_anno_log))])
  print(vst_anno[1:3, c(1:2, (ncol(vst_anno) - 5):ncol(vst_anno))])
  
  return(list(vst_anno_log = vst_anno_log, vst_anno = vst_anno))
}

# Apply the function
# All
result <- process_vst_data(dds_vst, norm_anno)
vst_anno_log <- result$vst_anno_log
vst_anno <- result$vst_anno

# All
vst_mean <- mean_function(input = vst_anno,
                          anno = sample_table,
                          condition = "condition")


# -----------------------------------------------------------------------------
# 6. Exploratory Data Analysis
# -----------------------------------------------------------------------------

# Function to generate plot annotations
generate_plot_annotations <- function(sample_table) {
  # Choose columns from the sample table for the heatmap annotation
  plot_annotation <- sample_table[, c("condition"), drop = FALSE]
  rownames(plot_annotation) <- sample_table$ID
  
  # Choose unique condition values for the heatmap annotation for mean expression
  plot_annotation_mean <- sample_table[, c("condition"), drop = FALSE]
  rownames(plot_annotation_mean) <- NULL
  plot_annotation_mean <- unique(plot_annotation_mean)
  rownames(plot_annotation_mean) <- plot_annotation_mean$condition
  
  return(list(plot_annotation = plot_annotation, plot_annotation_mean = plot_annotation_mean))
}

# Apply the function
result_annotations <- generate_plot_annotations(sample_table)
plot_annotation <- result_annotations$plot_annotation
plot_annotation_mean <- result_annotations$plot_annotation_mean


# -----------------------------------------------------------------------------
# 6.1 Frequencies of Gene Types
# -----------------------------------------------------------------------------

# Function to generate and save frequency plot for GENETYPE
generate_type_frequency_plot <- function(norm_anno, plot_title) {
  TypeCounts <- as.data.frame(table(norm_anno$GENETYPE))
  colnames(TypeCounts) <- c("Type", "Frequency")
  TypeCounts <- subset(TypeCounts, Frequency > 0)
  
  plot_obj <-
    ggplot(TypeCounts, aes(x = Type, y = Frequency, label = Frequency)) +
    geom_bar(stat = "identity",
             fill = "grey",
             colour = "grey") +
    theme_bw() +
    geom_text(size = 3, position = position_stack(vjust = 1)) +
    guides(fill = FALSE) +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          )) +
    xlab("")
  
  # Save Plot as an image file
  ggsave(paste0(dir_out, "/Plots/", plot_title),
         plot = plot_obj,
         device = "pdf")
  
  print(plot_obj)
}

# Apply the function
generate_type_frequency_plot(norm_anno, paste0("Frequencies_Gene_Types_MoClusters_k3_", cohort, ".pdf"))


# -----------------------------------------------------------------------------
# 6.2 Histogram of Means per Gene
# -----------------------------------------------------------------------------

# Function to generate and save histogram plot for rowMeans
generate_rowMeans_histogram <- function(dds, plot_title) {
  rMeans <- as.data.frame(log(rowMeans(counts(dds, normalized = TRUE)), 10))
  colnames(rMeans) <- "rowMeans"
  
  plot_obj <- ggplot(rMeans, aes(x = rowMeans)) +
    geom_histogram(bins = 100) +
    xlab("log10(rowMeans)") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) +
    theme_bw()
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", plot_title),
    plot = plot_obj,
    device = "pdf"
  )
  
  print(plot_obj)
}

# Apply the function
generate_rowMeans_histogram(dds, paste0("Histogram_Means_Per_Gene_MoClusters_k3_", cohort, ".pdf"))


# -----------------------------------------------------------------------------
# 6.3 Boxplots of Highest Expressed Genes
# -----------------------------------------------------------------------------

# Function to generate and save plot for top 50 highly expressed genes
generate_top50_genes_plot <- function(norm_anno, sample_table, plot_title) {
  plot_obj <- highestGenes(numGenes = 50,
                           data = norm_anno,
                           sample_table = sample_table) + ylim(c(0, 0.5 * 10 ^ 6))
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", plot_title),
    plot = plot_obj,
    device = "pdf",
    width = 9,
    height = 9
  )
  
  print(plot_obj)
}

# Apply the function
generate_top50_genes_plot(norm_anno, sample_table, paste0("Top50_Highly_Expressed_Genes_MoClusters_k3_", cohort, ".pdf"))


# -----------------------------------------------------------------------------
# 6.4 PCA: Percentage of Explained Variance
# -----------------------------------------------------------------------------

# Function to generate and save PCA eigenvalues/variances plots
generate_pca_eigenvalues_plot <- function(dds_vst, plot_titles) {
  # Extract the eigenvalues/variances of the principal dimensions
  eigenvalue <- get_eig(prcomp(t(assay(dds_vst))))
  eigenvalue$dim <- as.numeric(c(1:nrow(eigenvalue)))
  
  # Plot for all dimensions
  plot_obj_all <- ggplot(eigenvalue, aes(dim, variance.percent)) +
    geom_bar(stat = "identity") +
    geom_line(aes(dim, variance.percent)) +
    geom_point(aes(dim, variance.percent)) +
    geom_line(aes(dim, cumulative.variance.percent), colour = "grey") +
    geom_point(aes(dim, cumulative.variance.percent), colour = "grey") +
    scale_x_continuous(breaks = c(1:nrow(eigenvalue))) +
    xlab("Dimensions") +
    ylab("Percentage of explained variances") +
    theme_bw()
  
  # Save Plot as an image file for all dimensions
  ggsave(
    paste0(dir_out, "/Plots/", plot_titles),
    plot = plot_obj_all,
    device = "pdf"
  )
  
  print(plot_obj_all)
}

# Apply the function
generate_pca_eigenvalues_plot(dds_vst, paste0("PCA_Eigenvalues_MoClusters_k3_", cohort, ".pdf"))


# -----------------------------------------------------------------------------
# 6.5 PCA: Sample Plot
# -----------------------------------------------------------------------------

# Function to generate and save PCA plot
generate_pca_plot <- function(pca_input_df, sample_table, col_condition, plot_title) {
  p <- plotPCA(pca_input = pca_input_df,
    ntop = "all",
    xPC = 1,
    yPC = 2,
    color = "condition",
    anno_colour = col_condition,
    point_size = 3,
    label = NULL,
    title = paste0("PCA of Clusters (", tools::toTitleCase(cohort), " Cohort)"),
    pca_sample_table = sample_table
    #max.overlaps = 20
  )
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", plot_title),
    plot = p,
    device = "pdf"
  )
  
  ggsave(
    paste0(dir_out, "/Plots/", plot_title),
    plot = p,
    device = "tiff"
  )
  
  print(p)
}

# Apply the function
generate_pca_plot(dds_vst, sample_table, col_condition, paste0("PCA_condition_Genes_PC1_and_PC2_MoClusters_k3_", cohort, ".pdf"))


# -----------------------------------------------------------------------------
# 6.6 PCA: Linear Regression of PCs vs. Metadata (Batch Covariate Screening)
# -----------------------------------------------------------------------------

# The final covariate set below was arrived at after screening candidate
# batch/clinical covariates against PC1-PC5 (Section 6.6) and correlation
# with each other (Section 6.7); earlier, broader candidate lists are not
# kept here.
meta_variables <- c(
    "AGE",
    "PLATE",
    "SEX_BIRTH",
    "BMI_BASELINE",
    "CENTER_RUMC",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "NIJMEGEN_SEASON_SCORE_SINE",
    "NIJMEGEN_SEASON_SCORE_COSINE",
    "ETHNICITY_Asian",
    "ETHNICITY_Black",
    "ETHNICITY_Hispanic",
    "ETHNICITY_Mixed",
    "ETHNICITY_White",
    "PANDEMIC_BEFOREAFTER",
    "COVID19",
    "COVID_VACC",
    "TIMETOLAB",
    "MUTHIV",
    "EARLY_CART",
    "CART_DURATION",
    "CD4_NADIR",
    "VL_ZENITH",
    "HEPC",
    "HEPA_BASELINE",
    "HEPB_BASELINE",
    "HIV_STAGE_CDC",
    "HIV_DURATION"
  )

plot_pc_regression <- function(input = removedbatch_dds_vst,
                               ntop = "all",
                               meta_variables, 
                               nPCs = 10,
                               title = "PC contribution",
                               sample_table = sample_table, 
                               minValue, 
                               maxValue){
  
  #compute the PCA embedding outside the original function (optionally select the most variable features of the data)
  
  if(class(input) == "DESeqTransform"){
    if(ntop == "all") {
      select <- order(rowVars(as.matrix(assay(input))), decreasing=TRUE)[1:nrow(as.matrix(assay(input)))]
    } else {
      select <- order(rowVars(as.matrix(assay(input))), decreasing=TRUE)[1:ntop]
    }
    pca <- prcomp(t(as.matrix(assay(input))[select,]))
    
  } else if(class(input) == "data.frame") {
    if(ntop == "all") {
      select <- order(rowVars(input), decreasing=TRUE)[1:nrow(input)]
    } else {
      select <- order(rowVars(input), decreasing=TRUE)[1:ntop]
    }
    pca <- prcomp(t(as.matrix(input)[select,]))
    
  } else { print("unknown input format")}
  df_pca <- as.data.frame(pca[["x"]])
  
  
  # Store PCA scores
  df_pca <- as.data.frame(pca[["x"]])
  
  # Calculate variance explained %
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  var_explained <- round(var_explained, 1)
  
  # Rename PC columns with variance explained
  colnames(df_pca) <- paste0(colnames(df_pca), " (", var_explained, "%)")
  
  
  # Calculate variance attributed to metadata
  M <- sample_table[ , which(colnames(sample_table) %in% meta_variables)]
  
  for(i in colnames(M)){
    if(length(unique(M[ ,i])) <2 ){
      print(paste("exclude", i, sep = " "))
    }
  }
  
  pc_adj_r_squared <- matrix(NA, ncol = dim(df_pca)[2], nrow = dim(M)[2])
  for(i in 1:dim(df_pca)[2]){
    for(j in 1:dim(M)[2]){
      pc_adj_r_squared[j,i] <- summary(lm(df_pca[,i] ~ M[,j], na.action = na.exclude))$adj.r.squared
    }
  }
  
  pc_adj_r_squared <- as.data.frame(pc_adj_r_squared)
  colnames(pc_adj_r_squared) <- colnames(df_pca)
  rownames(pc_adj_r_squared) <- colnames(M)
  
  df <- pc_adj_r_squared[, 1:nPCs]
  
  paletteLength<-50
  my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  breakList <- c(seq(minValue, 0, length.out = ceiling(paletteLength/2) + 1), 
                 seq(maxValue/paletteLength, maxValue, length.out = floor(paletteLength/2)))
  
  hm <- pheatmap(as.matrix(df),
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = T,
                 show_colnames = T,
                 scale = "none",
                 main = title,
                 display_numbers = T,
                 color= my_palette,
                 breaks = breakList)
}


# Function to generate and save plot for principal component regression
generate_pc_regression_plot <- function(dds_vst, meta_variables, sample_table, plot_title) {
  plot_obj <- plot_pc_regression(
    input = dds_vst,
    ntop = "all",
    nPCs = 10,
    meta_variables = meta_variables,
    sample_table = sample_table, 
    minValue = -0.5, 
    maxValue = 0.5
  )
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", plot_title),
    plot = plot_obj,
    device = "pdf"
  )
  
  print(plot_obj)
}

# Apply the function
generate_pc_regression_plot(dds_vst, meta_variables, sample_table, paste0("PC_Confounders_Before_Correction_MoClusters_k3_", cohort, ".pdf"))


# -----------------------------------------------------------------------------
# 6.7 Metadata Correlation Matrix (Multicollinearity Check)
# -----------------------------------------------------------------------------

# Load required packages
library(ggplot2)
library(reshape2)
library(pheatmap)

# Call the function to generate the correlation matrix heatmap
plot_obj <-
  plot_meta_correlation_VR(sample_table,
                           meta_variables,
                           title = "Metadata Correlation",
                           threshold = 0.2)

# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", paste0("metadata_correlation_before_correction_MoClusters_k3_", cohort, ".pdf")),
  plot = plot_obj$plot,
  device = "pdf"
)

print(plot_obj)


# Function to generate and save PCA plots
generate_multiple_pca_plots <- function(pca_input, sample_table, dir_out, suffix = "") {
  # PCA based on CENTER_RUMC
  plot_obj_center <- plotPCA(
    pca_input, 
    ntop = "all",
    xPC = 1,
    yPC = 2,
    color = "CENTER",
    anno_colour = "NULL",
    shape = "NULL",
    point_size = 3,
    title = paste("PCA based on variance-stabilized counts", suffix),
    label = "ID",
    pca_sample_table = sample_table
  )
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", "PCA_center_before_correction_MoClusters_k3", suffix, ".pdf"),
    plot = plot_obj_center,
    device = "pdf"
  )
  
  print(plot_obj_center)
  
  # PCA based on ETHNICITY
  plot_obj_ethnicity <- plotPCA(
    pca_input, 
    ntop = "all",
    xPC = 1,
    yPC = 2,
    color = "ETHNICITY",
    anno_colour = "NULL",
    shape = "NULL",
    point_size = 3,
    title = paste("PCA based on variance-stabilized counts", suffix),
    label = "ID",
    pca_sample_table = sample_table
  )
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", "PCA_ethnicity_before_correction_MoClusters_k3", suffix, ".pdf"),
    plot = plot_obj_ethnicity,
    device = "pdf"
  )
  
  print(plot_obj_ethnicity)
  
  # PCA based on SEX_BIRTH
  plot_obj_sex <- plotPCA(
    pca_input, 
    ntop = "all",
    xPC = 1,
    yPC = 2,
    color = "SEX_BIRTH",
    anno_colour = "NULL",
    shape = "NULL",
    point_size = 3,
    title = paste("PCA based on variance-stabilized counts", suffix),
    label = "ID",
    pca_sample_table = sample_table
  )
  
  # Save Plot as an image file
  ggsave(
    paste0(dir_out, "/Plots/", "PCA_sex_before_correction_MoClusters_k3", suffix, ".pdf"),
    plot = plot_obj_sex,
    device = "pdf"
  )
  
  print(plot_obj_sex)
}

# Apply the function
generate_multiple_pca_plots(dds_vst, sample_table, dir_out, "_")


# -----------------------------------------------------------------------------
# 6.8 PCA Loadings
# -----------------------------------------------------------------------------
plotLoadings(PC="PC1", ntop="all")


# -----------------------------------------------------------------------------
# 6.9 LIMMA: Known Batch Effects (Diagnostic)
# -----------------------------------------------------------------------------

# the recommended correction is for SEX_BIRTH, CENTER (RUMC vs all other), PLATE (technical) and the season in form of the sine and cosine scores defined before
# Assuming sample_table$CENTER_RUMC is a factor/character variable
# Convert to binary numeric representation

# Function to remove batch effects
remove_batch_effects <- function(dds_vst, sample_table, suffix = "") {
  removedbatch_dds_vst <- as.data.frame(
    removeBatchEffect_RK(
      x = as.matrix(assay(dds_vst)),
      batch = sample_table[, colnames(sample_table) == "PLATE"],
      batch2 = sample_table[, colnames(sample_table) == "SEX_BIRTH"],
      batch3 = sample_table[, colnames(sample_table) == "CENTER_RUMC"],
      covariates = sample_table[, colnames(sample_table) %in%  c(
        "NIJMEGEN_SEASON_SCORE_SINE",
        "NIJMEGEN_SEASON_SCORE_COSINE",
        "TIMETOLAB",
        "AGE",
        "PC1",
        "COVID_VACC"
      )],
      design = model.matrix( ~ condition, data = sample_table)
    )
  )
  
  # Save the adjusted data frame
  write.table(
    removedbatch_dds_vst,
    file = paste0(dir_out, "/RemovedBatchEffect", suffix, ".txt"),
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE
  )
  
  print(head(removedbatch_dds_vst))
  
  return(removedbatch_dds_vst)
}

# Apply the function
removedbatch_dds_vst <- remove_batch_effects(dds_vst, sample_table, paste0("_", cohort))

# All
bc_anno_list <- generate_bc_anno(dataframe = removedbatch_dds_vst)
bc_anno <- bc_anno_list$norm_anno
gene_annotation_bc <- bc_anno_list$gene_annotation


# Plot the regression plot after correction
generate_pc_regression_plot(removedbatch_dds_vst, meta_variables, sample_table, paste0("PC_Confounders_After_Correction_MoClusters_k3_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, ".pdf"))

# Apply the function
generate_multiple_pca_plots(removedbatch_dds_vst, sample_table, dir_out, "after_correction_All_MoClusters_k3_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc")

# Apply the function
generate_pca_plot(removedbatch_dds_vst, sample_table, col_condition, "PCA_condition_All_Genes_PC1_and_PC2_MoClusters_k3_after_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc.pdf")


# -----------------------------------------------------------------------------
# 7.2 Include Batch Effect Variables in the DESeq2 Model
# -----------------------------------------------------------------------------

# Function to preprocess and run DESeq
run_deseq <- function(dds, sample_table, suffix = "", include_center = TRUE) {
  
  # Convert categorical covariates to factor
  dds$SEX_BIRTH <- as.numeric(sample_table$SEX_BIRTH)
  dds$PLATE <- as.numeric(sample_table$PLATE)
  
  if (include_center == TRUE) {
    dds$CENTER_RUMC <- as.numeric(sample_table$CENTER_RUMC)
  }
  
  # Convert numeric covariates to numeric
  dds$NIJMEGEN_SEASON_SCORE_SINE <- as.numeric(sample_table$NIJMEGEN_SEASON_SCORE_SINE)
  dds$NIJMEGEN_SEASON_SCORE_COSINE <- as.numeric(sample_table$NIJMEGEN_SEASON_SCORE_COSINE)
  dds$PC1 <- as.numeric(sample_table$PC1)
  dds$PC2 <- as.numeric(sample_table$PC2)
  dds$PC3 <- as.numeric(sample_table$PC3)
  dds$AGE <- as.numeric(sample_table$AGE)
  dds$TIMETOLAB <- as.numeric(sample_table$TIMETOLAB)
  dds$ETHNICITY_Black <- as.numeric(sample_table$ETHNICITY_Black)
  
  
  # Update design
  if (include_center == TRUE) {
    design(dds) <-
      ~ NIJMEGEN_SEASON_SCORE_SINE + NIJMEGEN_SEASON_SCORE_COSINE + PLATE + TIMETOLAB + CENTER_RUMC + COVID_VACC + AGE + SEX_BIRTH + PC1 + condition
  } else {
    design(dds) <-
      ~ NIJMEGEN_SEASON_SCORE_SINE + NIJMEGEN_SEASON_SCORE_COSINE + PLATE + TIMETOLAB + COVID_VACC + AGE + SEX_BIRTH + PC1 + condition

    print(design(dds))
  }
  
  # Run DESeq
  dds <- DESeq(dds, parallel = FALSE)
  
  # Save DESeq object
  saveRDS(
    dds,
    file = paste0(dir_out, "/DESeq_Object", suffix, ".RDS")
  )
  
  return(dds)
}

# Run DESeq2 with the full batch-effect design for this cohort
gc()
dds_final <- run_deseq(dds, sample_table, include_center = TRUE, paste0("_MoClusters_k3_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort))
gc()


# Save the results
## Export the object to an RDS file and Excel file
saveRDS(dds, file = paste0(dir_out, "DESeq2_Results/DDS_Initial_", Sys.Date(), ".rds"))


# -----------------------------------------------------------------------------
# 7.3 Batch-Corrected Counts
# -----------------------------------------------------------------------------

# Function to process data frames
process_df <- function(df, norm_anno) {
  # Add GENEID column
  df$GENEID <- rownames(df)
  
  # Merge with normalization annotation
  df <- merge(df, norm_anno[, c("GENEID", "SYMBOL", "GENETYPE", "DESCRIPTION", "CHR")], by = "GENEID")
  
  # Set row names
  rownames(df) <- df$GENEID
  
  # Return processed data frame
  return(df)
}

# All
# Process log values
removed_batch_anno_log <- process_df(removedbatch_dds_vst, norm_anno)
# Process unlog values
removed_batch_anno <- process_df(2 ^ removedbatch_dds_vst, norm_anno)


# -----------------------------------------------------------------------------
# 7.4 Seasonality Check (Before vs. After Correction)
# -----------------------------------------------------------------------------

color_GROUP_CLINICAL <- c(
  "EC_PERSISTENT" = "#6a3d9a",
  "EC_TRANSIENT" = "#e31a1c",
  "HIV" = "#ff7f00",
  "RP" = "blue"
)

# chunk options: , fig.width=20, fig.height=14, cache = TRUE
plot_obj <- season_feature_plot(
  data_table = norm_anno,
  # your data table with values to show (e.g. normalized gene expression table)
  meta_table = sample_table,
  # your meta data table to add respective metadata to your data_table
  ID_in_meta_for_merge = "LAB_ID",
  # name of the ID in your metadata which is identical to the column names of your data_table, enables merging
  meta_columns_to_include = c(
    "DONOR_ID",
    "GROUP_CLINICAL",
    "NIJMEGEN_DATE_COLLECTION",
    "AGE",
    "SEX_BIRTH"
  ),
  # which metadata to add? the collection date must be included!!
  collection_date_name = "NIJMEGEN_DATE_COLLECTION",
  # the collection date name in your meta table
  color_by = "GROUP_CLINICAL",
  # choose what to color by, e.g. the clinical group or sex (must be included previously in meta_columns_to_include, for sex convert to factors first)
  colors = color_GROUP_CLINICAL,
  # named color vector for coloring the metadata
  feature_to_check = c(
    "ARSI",
    "ARHGAP8",
    "ERMN",
    "DUSP5",
    "SNX9",
    "FOS",
    "FOSB",
    "ALPL",
    "CXCR1",
    "MME",
    "CXCR2",
    "DUSP4"
  ),
  # features you want to check, these are row.names in your data_table
  ncol = 3
) # number of columns for plotting

ggsave(
  paste0(dir_out, "/plots/", "seasonality_plot_before_correction_DTG_TDF_FTC_vs_DTG_TAF_FTC_AgeSexSeasonCenterRUMCPlate.pdf"),
  plot = plot_obj,
  device = "pdf"
)

print(plot_obj)


plot_obj <- season_feature_plot(
  data_table = removed_batch_anno_log,
  # your data table with values to show (e.g. normalized gene expression table)
  meta_table = sample_table,
  # your meta data table to add respective metadata to your data_table
  ID_in_meta_for_merge = "LAB_ID",
  # name of the ID in your metadata which is identical to the column names of your data_table, enables merging
  meta_columns_to_include = c(
    "DONOR_ID",
    "GROUP_CLINICAL",
    "NIJMEGEN_DATE_COLLECTION",
    "AGE",
    "SEX_BIRTH"
  ),
  # which metadata to add? the collection date must be included!!
  collection_date_name = "NIJMEGEN_DATE_COLLECTION",
  # the collection date name in your meta table
  color_by = "GROUP_CLINICAL",
  # choose what to color by, e.g. the clinical group or sex (must be included previously in meta_columns_to_include, for sex convert to factors first)
  colors = color_GROUP_CLINICAL,
  # named color vector for coloring the metadata
  feature_to_check = c(
    "ARSI",
    "ARHGAP8",
    "ERMN",
    "DUSP5",
    "SNX9",
    "FOS",
    "FOSB",
    "ALPL",
    "CXCR1",
    "MME",
    "CXCR2",
    "DUSP4"
  ),
  # features you want to check, these are row.names in your data_table
  ncol = 3
) # number of columns for plotting

ggsave(
  paste0(dir_out, "/plots/", "seasonality_plot_after_correction_DTG_TDF_FTC_vs_DTG_TAF_FTC_AgeSexSeasonCenterRUMCPlate.pdf"),
  plot = plot_obj,
  device = "pdf"
)

print(plot_obj)


# -----------------------------------------------------------------------------
# 8.1 Define Relevant Comparisons
# -----------------------------------------------------------------------------

# insert here your comparison, "comparison" vs "control" for each DE call, e.g. here its HIV vs FAM or EC_PERSISTENT vs HIV
comparison_table <-
  data.frame(
    comparison = c("Mixed", "AllHigh", "AllHigh"),
    control = c("AllLow", "AllLow", "Mixed")
  )

comparison_table


# -----------------------------------------------------------------------------
# 8.2 Differential Expression Testing
# -----------------------------------------------------------------------------

# chunk options: , cache = TRUE
# Function to check if a data frame is empty
is_empty <- function(df) {
  nrow(df) == 0 && ncol(df) == 0
}

# Function to perform DE analysis and save results
perform_DE_analysis <-
  function(dds,
           comparison_table,
           dir_out,
           comparison_name) {
    dds_dea <- dds
    
    DEresults_list <- list()
    
    for (i in unique(comparison_table$control)) {
      dds_dea$condition <- relevel(dds_dea$condition, i)
      dds_dea <- nbinomWaldTest(object = dds_dea)
      comparison_table_subset <-
        comparison_table[comparison_table$control == i, ]
      
      DEresults <- DEAnalysis(
        input = dds_dea,
        comparison_table = comparison_table_subset,
        condition = "condition",
        alpha = 0.05 ,
        lfcThreshold = 0,
        sigFC = 1.5,
        multiple_testing = "BH",
        shrinkage = TRUE,
        shrinkType = "apeglm"
      )
      
      DEresults_list <- c(DEresults_list, DEresults)
    }
    
    DEresults <- DEresults_list[unique(names(DEresults_list))]
    
    # Save DE results
    save(
      DEresults,
      file = paste0(
        dir_out,
        "DESeq2_Results/",
        comparison_name,
        ".RData"
      )
    )
    
    # Create Excel workbook
    wb <- createWorkbook()
    
    addWorksheetWithData <- function(workbook, sheet_name, data) {
      addWorksheet(workbook, sheet_name)
      writeData(workbook, sheet_name, data)
    }
    
    # Loop over the comparisons and add results to the workbook
    comparisons <- c("Mixed_vs_AllLow", "AllHigh_vs_AllLow", "AllHigh_vs_Mixed")
    
    for (comparison in comparisons) {
      results <- as.data.frame(DEresults[[comparison]]@results)
      number_DE_genes <-
        as.data.frame(DEresults[[comparison]]@Number_DE_genes)
      upregulated_genes <-
        as.data.frame(DEresults[[comparison]]@DE_genes$up_regulated_Genes)
      downregulated_genes <-
        as.data.frame(DEresults[[comparison]]@DE_genes$down_regulated_Genes)
      
      # Check if data frames are empty before adding them to the workbook
      if (!is_empty(results)) {
        addWorksheetWithData(wb, paste0(comparison, "_1"), results)
      }
      
      if (!is_empty(number_DE_genes)) {
        addWorksheetWithData(wb, paste0(comparison, "_2"), number_DE_genes)
      }
      
      if (!is_empty(upregulated_genes)) {
        addWorksheetWithData(wb, paste0(comparison, "_3"), upregulated_genes)
      }
      
      if (!is_empty(downregulated_genes)) {
        addWorksheetWithData(wb, paste0(comparison, "_4"), downregulated_genes)
      }
    }
    
    # Save Excel file
    excel_file <-
      paste0(dir_out,
             "DESeq2_Results/",
             comparison_name,
             ".xlsx")
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    
    return(DEresults)
  }

# All dataset
DEresults <- perform_DE_analysis(dds_final, comparison_table, dir_out, paste0("DESeq2_DEresults_MoClusters_k3_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort))
gc()


# -----------------------------------------------------------------------------
# 9.1 Summary of DE Genes
# -----------------------------------------------------------------------------

# DEresults is already in memory from Section 8.2 above; no need to reload
# from disk here.

# Function to extract DE counts for each comparison
extract_DE_counts <- function(DEresults, comparison) {
  DEcounts <- NULL
  for (i in c(comparison)) {
    print(i)
    tmp <- unlist(DEresults[[i]]@Number_DE_genes)
    DEcounts <- rbind(DEcounts, tmp)
  }
  rownames(DEcounts) <- names(DEresults)[-1]
  return(DEcounts)
}

# Extract DE counts for all datasets
DEcounts <- extract_DE_counts(DEresults, c("Mixed_vs_AllLow", "AllHigh_vs_Mixed", "AllHigh_vs_AllLow"))

# Melt the data for plotting
DEcounts_melt_all <- reshape2::melt(DEcounts)


# Plot for all datasets
ggplot(DEcounts_melt_all, aes(x = Var1, fill = Var2)) +
  geom_col(aes(y = value), position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) +
  scale_fill_manual(
    values = c(
      "up_regulated_Genes" = "firebrick",
      "down_regulated_Genes" = "cornflowerblue"
    )
  ) +
  labs(title = "Differential Expression Counts (MoCluster)")


# -----------------------------------------------------------------------------
# 9.2 Upset Plot
# -----------------------------------------------------------------------------

library(UpSetR)
upset_list <-
  list(
    Mixed_vs_AllLow_up = DEresults[["Mixed_vs_AllLow"]]@DE_genes[["up_regulated_Genes"]][["SYMBOL"]],
    Mixed_vs_AllLow_down = DEresults[["Mixed_vs_AllLow"]]@DE_genes[["down_regulated_Genes"]][["SYMBOL"]],
    AllHigh_vs_Mixed_up = DEresults[["AllHigh_vs_Mixed"]]@DE_genes[["up_regulated_Genes"]][["SYMBOL"]],
    AllHigh_vs_Mixed_down = DEresults[["AllHigh_vs_Mixed"]]@DE_genes[["down_regulated_Genes"]][["SYMBOL"]],
    AllHigh_vs_AllLow_up = DEresults[["AllHigh_vs_AllLow"]]@DE_genes[["up_regulated_Genes"]][["SYMBOL"]],
    AllHigh_vs_AllLow_down = DEresults[["AllHigh_vs_AllLow"]]@DE_genes[["down_regulated_Genes"]][["SYMBOL"]]
  )

upset(fromList(upset_list), nsets = 6)


# -----------------------------------------------------------------------------
# 9.3 Hierarchical Clustering of the Union of DE Genes
# -----------------------------------------------------------------------------
# the uDEG() function produces the union of the DE genes from the specified comparisons
allDEgenes <-
  uDEG(
    comparisons = c(
      "Mixed_vs_AllLow",
      "AllHigh_vs_Mixed", 
      "AllHigh_vs_AllLow"
    )
  )

plot <- plotHeatmap(
  geneset = allDEgenes, #allDEgenes
  title = "Heatmap of all differentially expressed genes",
  key = "Ensembl",
  show_rownames = T,
  cluster_cols = F
)

ggsave(
  paste0(dir_out, "/Plots/", paste0("heatmap_degs_MoClusters_k3_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, ".pdf")),
  plot = plot,
  device = "pdf"
)

print(plot)

# -----------------------------------------------------------------------------
# 9.4 Hierarchical Clustering of DE Transcription Factors
# -----------------------------------------------------------------------------

# chunk options: , fig.height=6, fig.width=10, cache = TRUE

if(organism=="mouse") {
  DE_TF <-
    allDEgenes[which(allDEgenes %in% norm_anno[norm_anno$SYMBOL %in% TFlist_mm, ]$GENEID)]
} else if (organism == "human") {
  DE_TF <-
    allDEgenes[which(allDEgenes %in% norm_anno[norm_anno$SYMBOL %in% TFlist_hs$HGNC.symbol, ]$GENEID)]
}

TF_anno <- norm_anno[norm_anno$GENEID %in% DE_TF, ]

plotHeatmap(
  geneset = DE_TF,
  title = "Heatmap of all differentially expressed transcription factors",
  show_rownames = T,
  cluster_cols = F
  
)


# -----------------------------------------------------------------------------
# 9.5 Venn Diagrams
# -----------------------------------------------------------------------------

# chunk options: , fig.height=5, fig.width=10
plotVenn(
  comparisons =  c(
    "Mixed_vs_AllLow", 
    "AllHigh_vs_Mixed", 
    "AllHigh_vs_AllLow"
  ),
  regulation = "up"
)

plotVenn(
  comparisons =  c(
    "Mixed_vs_AllLow", 
    "AllHigh_vs_Mixed", 
    "AllHigh_vs_AllLow"
  ),
  regulation = "down"
)

# NOTE: DEresults is keyed as "<comparison>_vs_<control>" (matching how
# perform_DE_analysis() names its results, e.g. "Mixed_vs_AllLow" — there is no
# "AllLow_vs_Mixed" entry). Genes "up in AllLow vs Mixed" are equivalently "down
# in Mixed vs AllLow", so the existing Mixed_vs_AllLow object is used throughout
# with the up/down direction flipped accordingly.
DEGs_unique_Mixed_up <- rbind(DEresults$Mixed_vs_AllLow@DE_genes$up_regulated_Genes,
                             DEresults$AllHigh_vs_Mixed@DE_genes$down_regulated_Genes)

DEGs_unique_AllLow_up <- rbind(DEresults$Mixed_vs_AllLow@DE_genes$down_regulated_Genes,
                             DEresults$AllHigh_vs_AllLow@DE_genes$down_regulated_Genes)

DEGs_unique_AllHigh_up <- rbind(DEresults$AllHigh_vs_Mixed@DE_genes$up_regulated_Genes,
                             DEresults$AllHigh_vs_AllLow@DE_genes$up_regulated_Genes)

# Filter out genes that are upregulated in AllLow or AllHigh
DEGs_unique_Mixed_up <- DEGs_unique_Mixed_up %>%
  filter(!(SYMBOL %in% DEGs_unique_AllLow_up$SYMBOL)) %>%
  filter(!(SYMBOL %in% DEGs_unique_AllHigh_up$SYMBOL))

DEGs_unique_Mixed_up

# Filter out genes that are upregulated in Mixed or AllHigh
DEGs_unique_AllLow_up <- DEGs_unique_AllLow_up %>%
  filter(!(SYMBOL %in% DEGs_unique_Mixed_up$SYMBOL)) %>%
  filter(!(SYMBOL %in% DEGs_unique_AllHigh_up$SYMBOL))

DEGs_unique_AllLow_up

DEGs_unique_AllHigh_up <- DEGs_unique_AllHigh_up %>%
  filter(!(SYMBOL %in% DEGs_unique_Mixed_up$SYMBOL)) %>%
  filter(!(SYMBOL %in% DEGs_unique_AllLow_up$SYMBOL))

DEGs_unique_AllHigh_up

# Get all DEGs for each comparison (union of up + down is direction-invariant)
DEGs_Mixed_vs_AllLow <- c(DEresults$Mixed_vs_AllLow@DE_genes$up_regulated_Genes$SYMBOL,
                       DEresults$Mixed_vs_AllLow@DE_genes$down_regulated_Genes$SYMBOL)

DEGs_Mixed_vs_AllHigh <- c(DEresults$AllHigh_vs_Mixed@DE_genes$up_regulated_Genes$SYMBOL,
                       DEresults$AllHigh_vs_Mixed@DE_genes$down_regulated_Genes$SYMBOL)

DEGs_AllLow_vs_AllHigh <- c(DEresults$AllHigh_vs_AllLow@DE_genes$up_regulated_Genes$SYMBOL,
                       DEresults$AllHigh_vs_AllLow@DE_genes$down_regulated_Genes$SYMBOL)

# Find shared DEGs across all comparisons
shared_DEGs <- Reduce(intersect, list(DEGs_Mixed_vs_AllLow, DEGs_Mixed_vs_AllHigh, DEGs_AllLow_vs_AllHigh))

# If you want to find DEGs shared between any two comparisons:
shared_Mixed_vs_AllLow_Mixed_vs_AllHigh <- intersect(DEGs_Mixed_vs_AllLow, DEGs_Mixed_vs_AllHigh)
shared_Mixed_vs_AllLow_AllLow_vs_AllHigh <- intersect(DEGs_Mixed_vs_AllLow, DEGs_AllLow_vs_AllHigh)
shared_Mixed_vs_AllHigh_AllLow_vs_AllHigh <- intersect(DEGs_Mixed_vs_AllHigh, DEGs_AllLow_vs_AllHigh)

shared_Mixed_vs_AllLow_Mixed_vs_AllHigh
shared_Mixed_vs_AllLow_AllLow_vs_AllHigh
shared_Mixed_vs_AllHigh_AllLow_vs_AllHigh
# You can then use these shared gene lists for further analysis or visualization


# -----------------------------------------------------------------------------
# 9.6 MA Plots
# -----------------------------------------------------------------------------

# chunk options: , cache = TRUE
plot_obj <- plotMA(comparison = "Mixed_vs_AllLow", ylim=c(-2,2))
plot_obj
# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", "MAplot_Mixed_vs_AllLow_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc.pdf"),
  plot = plot_obj,
  device = "pdf"
)

print(plot_obj)


# -----------------------------------------------------------------------------
# 9.7 P-Value Distributions
# -----------------------------------------------------------------------------

# chunk options: , cache = TRUE
plot_obj <- plotPvalues(comparison="Mixed_vs_AllLow")
plot_obj
# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", "pvalue_histogram_Mixed_vs_AllLow_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVaccCovidvacc.pdf"),
  plot = plot_obj,
  device = "pdf"
)
print(plot_obj)

# chunk options: , cache = TRUE
plot_obj <- plotDEHeatmap(comparison = "AllHigh_vs_AllLow",
              factor = "condition",
              conditions= "all", sample_annotation = sample_table)

# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", paste0("heatmap_degs_AllHigh_vs_AllLow_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, ".pdf")),
  plot = plot_obj,
  device = "pdf"
)
print(plot_obj)


# -----------------------------------------------------------------------------
# 9.8 Volcano Plots
# -----------------------------------------------------------------------------

################
# AllLow vs Mixed #
################

# Plot Volcano Plot
plot_obj_AllLow_vs_Mixed <-
  plotVolcano_xlsx2(
    comparison = "Mixed_vs_AllLow",
    DE_results = DEresults$Mixed_vs_AllLow@results,
    labelnum = 20
  )
plot_obj_AllLow_vs_Mixed


# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", paste0("Volcano_plot_MoClusters_k3_Mixed_vs_AllLow_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, ".pdf")),
  plot = plot_obj_AllLow_vs_Mixed,
  width = 7,
  height = 7
)
print(plot_obj_AllLow_vs_Mixed)


################
# AllHigh vs Mixed #
################

# Plot Volcano Plot
plot_obj_AllHigh_vs_Mixed <-
  plotVolcano_xlsx2(comparison = "AllHigh_vs_Mixed",
                   DE_results = DEresults$AllHigh_vs_Mixed@results,
                   labelnum = 20)


# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", paste0("Volcano_plot_MoClusters_k3_AllHigh_vs_Mixed_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, ".pdf")),
  plot = plot_obj_AllHigh_vs_Mixed,
  width = 7,
  height = 7
)
print(plot_obj_AllHigh_vs_Mixed)


################
# AllHigh vs AllLow #
################

# Plot Volcano Plot
plot_obj_AllHigh_vs_AllLow <-
  plotVolcano_xlsx2(
    comparison = "AllHigh_vs_AllLow",
    DE_results = DEresults$AllHigh_vs_AllLow@results,
    labelnum = 20
  )

# Save Plot as an image file
ggsave(
  paste0(dir_out, "/Plots/", paste0("Volcano_plot_MoClusters_k3_AllHigh_vs_AllLow_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, ".pdf")),
  plot = plot_obj_AllHigh_vs_AllLow,
  width = 7,
  height = 7
)
print(plot_obj_AllHigh_vs_AllLow)

# chunk options: fig.height=12, fig.width=12
library(EnhancedVolcano)

# Single parameterized version of the two near-identical volcano-plot
# functions from the original notebook. `fc_fold` sets the fold-change
# threshold used both for point coloring/labeling and for EnhancedVolcano's
# own FCcutoff line; `label_size`/`connector_width`/plot `width`/`height`
# match the two original variants (labSize 5 vs 6, widthConnectors 1.0 vs
# 0.5, 12x12in vs 9x13in) and are exposed as arguments below instead.
plotEnhancedVolcano <-
  function(data,
           comparison,
           fc_fold = 1.5,
           xlim = c(-1.25, 1.5),
           ylim = c(-0.1, 17),
           label_size = 5,
           connector_width = 1.0,
           plot_width = 10,
           plot_height = 8,
           filename = NULL) {

    upFC <- log2(fc_fold)
    downFC <- -log2(fc_fold)

    # Filter up-regulated genes
    up_degs <- data[order(data$padj),]
    up_degs_labels <- up_degs[up_degs$log2FoldChange > upFC,]$SYMBOL
    
    # Filter down-regulated genes
    down_degs <- data[order(data$padj),]
    down_degs_labels <-
      down_degs[down_degs$log2FoldChange < downFC,]$SYMBOL
    
    # Create custom color labels
    keyvals <- ifelse(
      data$log2FoldChange < downFC & data$padj < 0.05,
      '#00468BFF',
      ifelse(
        data$log2FoldChange > upFC & data$padj < 0.05,
        '#ED0000FF',
        'grey'
      )
    )
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == '#ED0000FF'] <-
      'up-regulated in AllHigh'
    names(keyvals)[keyvals == 'grey'] <- 'n.s.'
    names(keyvals)[keyvals == '#00468BFF'] <-
      'up-regulated in AllLow'
    
    p <- EnhancedVolcano(
      data,
      lab = data$SYMBOL,
      x = 'log2FoldChange',
      y = 'padj',
      selectLab = c(up_degs_labels, down_degs_labels),
      xlab = bquote( ~ Log[2] ~ 'fold change'),
      ylab = bquote( ~ -Log[10] ~ 'FDR'),
      pCutoff = 0.05,
      FCcutoff = downFC,
      pointSize = 4.0,
      labSize = label_size,
      labCol = 'black',
      labFace = 'bold',
      boxedLabels = FALSE,
      colAlpha = 3 / 5,
      colCustom = keyvals,
      legendPosition = 'top',
      legendLabSize = 15,
      legendIconSize = 4.0,
      drawConnectors = TRUE,
      widthConnectors = connector_width,
      colConnectors = 'black',
      xlim = xlim,
      ylim = ylim,
      title = paste("Differentially Expressed Genes (DEG) of:", comparison),
      caption = "Corrected for Age, Sex, Seasonality, Plate, Time to lab, Center, and genetic PC1",
      subtitle = "FDR < 0.05"
    )
    
    # Save the plot as PNG if filename is provided
    if (!is.null(filename)) {
      ggsave(
        paste0(dir_out, "/Plots/",filename),
        p,
        width = plot_width,
        height = plot_height,
        units = "in",
        dpi = 300
      )
    }
    
    return(p)
  }

# Volcano plots at two fold-change thresholds (1.5x and 1.25x), each run
# across all three pairwise comparisons.
volcano_comparisons <- list(
  list(key = "Mixed_vs_AllLow", label = "AllLow vs Mixed"),
  list(key = "AllHigh_vs_AllLow", label = "AllLow vs AllHigh"),
  list(key = "AllHigh_vs_Mixed", label = "Mixed vs AllHigh")
)

for (vc in volcano_comparisons) {
  plotEnhancedVolcano(
    data = DEresults[[vc$key]]@results,
    comparison = vc$label,
    fc_fold = 1.5,
    label_size = 5,
    connector_width = 1.0,
    plot_width = 10, plot_height = 8,
    filename = paste0("VolcanoPlot_FDR_", vc$key, "_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, "_LFC1pt5.pdf")
  )
}


# chunk options: fig.height=13, fig.width=9
for (vc in volcano_comparisons) {
  plotEnhancedVolcano(
    data = DEresults[[vc$key]]@results,
    comparison = vc$label,
    fc_fold = 1.25,
    label_size = 6,
    connector_width = 0.5,
    plot_width = 9, plot_height = 13,
    filename = paste0("VolcanoPlot_FDR_", vc$key, "_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_", cohort, "_LFC1pt25.pdf")
  )
}


# chunk options: , cache = TRUE
my_comparisons <-
  list(c("Mixed", "AllLow"), c("Mixed", "AllHigh"), c("AllLow", "AllHigh"))

degs <- c(DEresults$Mixed_vs_AllLow@DE_genes$up_regulated_Genes$SYMBOL,
          DEresults$Mixed_vs_AllLow@DE_genes$down_regulated_Genes$SYMBOL,
          DEresults$AllHigh_vs_Mixed@DE_genes$up_regulated_Genes$SYMBOL,
          DEresults$AllHigh_vs_Mixed@DE_genes$down_regulated_Genes$SYMBOL,
          DEresults$AllHigh_vs_AllLow@DE_genes$up_regulated_Genes$SYMBOL,
          DEresults$AllHigh_vs_AllLow@DE_genes$down_regulated_Genes$SYMBOL) %>% unique()

degs_lfc_1p25 <- c(
  DEresults$Mixed_vs_AllLow@results[abs(DEresults$Mixed_vs_AllLow@results$log2FoldChange) > log2(1.25) &
                                   DEresults$Mixed_vs_AllLow@results$padj < 0.05, ]$SYMBOL,
  DEresults$AllHigh_vs_Mixed@results[abs(DEresults$AllHigh_vs_Mixed@results$log2FoldChange) > log2(1.25) &
                                   DEresults$AllHigh_vs_Mixed@results$padj < 0.05, ]$SYMBOL,
  DEresults$AllHigh_vs_AllLow@results[abs(DEresults$AllHigh_vs_AllLow@results$log2FoldChange) > log2(1.25) &
                                   DEresults$AllHigh_vs_AllLow@results$padj < 0.05, ]$SYMBOL
) %>% unique()


plot_gene(
  data = removed_batch_anno_log,
  symbols = degs_lfc_1p25,
  condition = "condition",
  anno_colour = col_condition,
  wilcox.test = T,
  my_comparisons = my_comparisons,
  ncol = 2,
  p.adjustment = "none"
)


# -----------------------------------------------------------------------------
# 10.1 Export Count Tables (TSV)
# -----------------------------------------------------------------------------

write.table(
  norm_anno,
  paste(
    dir_out,
    "Tables/",
    "DESeq2_norm_anno_", cohort, "_",
    Sys.Date(),
    ".tsv",
    sep = ""
  ),
  sep = "\t",
  quote = F,
  row.names = F
)

write.table(
  vst_anno_log,
  paste(
    dir_out,
    "Tables/",
    "DESeq2_vst_anno_log_", cohort, "_",
    Sys.Date(),
    ".tsv",
    sep = ""
  ),
  sep = "\t",
  quote = F,
  row.names = F
)

write.table(
  removed_batch_anno_log,
  paste(
    dir_out,
    "Tables/",
    "DESeq2_removedbatch_anno_log_", cohort, "_",
    Sys.Date(),
    ".tsv",
    sep = ""
  ),
  sep = "\t",
  quote = F,
  row.names = F
)
write.table(
  gene_annotation,
  paste(
    dir_out,
    "Tables/",
    "DESeq2_gene_annotation_", cohort, "_",
    Sys.Date(),
    ".tsv",
    sep = ""
  ),
  sep = "\t",
  quote = F,
  row.names = F
)
saveRDS(dds,
        paste0(dir_out, "Tables/", "DESeq2_dds_", cohort, "_", Sys.Date(), ".rds"))
save(
  DEresults,
  file = paste0(
    dir_out,
    "Tables/",
    "DESeq2_DEresults_simple_run_all_ethnicities_", cohort, "_",
    Sys.Date(),
    ".RData"
  )
)


# -----------------------------------------------------------------------------
# 10.2 Export Combined Results (Excel)
# -----------------------------------------------------------------------------

# create Workbook
ExcelOutput<-createWorkbook()

# add sample table
sheet <- addWorksheet(ExcelOutput, sheetName = "Samples")
writeDataTable(ExcelOutput, sheet, sample_table, withFilter=FALSE)

# add normalized counts
tmp <- norm_anno
sheet <- addWorksheet(ExcelOutput, sheetName = "Normalized counts & Annotation")
writeDataTable(ExcelOutput, sheet, tmp, withFilter=FALSE)

# add variance-stabilized counts
tmp <- as.data.frame(assay(dds_vst))
tmp$GENEID <- norm_anno$GENEID
tmp$SYMBOL <- norm_anno$SYMBOL
sheet <- addWorksheet(ExcelOutput, sheetName = "Variance-stabilized counts")
writeDataTable(ExcelOutput, sheet, tmp, withFilter=FALSE)

# add DE test parameters
tmp <- stack(unlist(DEresults$parameters))
colnames(tmp)<-c("value","parameter")
tmp <- rbind(tmp, data.frame(value = as.character(design(dds))[2], parameter = "design"))
sheet <- addWorksheet(ExcelOutput, sheetName = "DE parameters")
writeDataTable(ExcelOutput, sheet, tmp, withFilter=FALSE)

# add combined DE results
cResults <- NULL
for (i in 2:length(names(DEresults))) {
  if(i == 2){
    cResults<-data.frame(DEresults[[i]]@results[,c("GENEID","SYMBOL","baseMean","log2FoldChange","padj","regulation")])
    colnames(cResults)[4:6]<- paste(names(DEresults)[i],colnames(cResults)[4:6],sep="_")
  }else{
    tmp<-data.frame(DEresults[[i]]@results[,c("log2FoldChange","padj","regulation")])
    colnames(tmp)<- paste(names(DEresults)[i],colnames(tmp),sep="_")    
    cResults<-cbind(cResults,tmp)
  }
}
sheet <- addWorksheet(ExcelOutput, sheetName = "combined DE Results")
writeDataTable(ExcelOutput, sheet, cResults, withFilter=FALSE)

# add DE results in single sheets
for(i in 2:length(names(DEresults))){
  gc()
  sheet <- addWorksheet(ExcelOutput, sheetName = substr(names(DEresults[i]),1 , 30))
  writeDataTable(ExcelOutput, sheet, DEresults[[i]]@results, withFilter=FALSE)
  sheet <- addWorksheet(ExcelOutput, sheetName = paste(substr(names(DEresults[i]),1 , 22),"_upDEGs"))
  writeDataTable(ExcelOutput, sheet, DEresults[[i]]@DE_genes$up_regulated_Genes, withFilter=FALSE)
  sheet <- addWorksheet(ExcelOutput, sheetName = paste(substr(names(DEresults[i]),1 , 20),"_downDEGs"))
  writeDataTable(ExcelOutput, sheet, DEresults[[i]]@DE_genes$down_regulated_Genes, withFilter=FALSE)
  }

# Save Workbook
filename <- paste(dir_out,"Tables/AnalysisOutput_", cohort, "_",gsub(":","-",as.character(Sys.time())),".xlsx",sep="")
saveWorkbook(ExcelOutput, file=filename)


# -----------------------------------------------------------------------------
# 11. Save Image and Session Info
# -----------------------------------------------------------------------------

save.image(paste(dir_out, cohort, "_", Sys.Date(), "_Image.RData", sep = ""))
sessionInfo()

