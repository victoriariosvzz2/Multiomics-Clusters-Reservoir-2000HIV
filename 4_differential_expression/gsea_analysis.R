#!/usr/bin/env Rscript
# =============================================================================
# Script: gsea_analysis.R
# Project: Multi-Omics Clustering - 2000HIV Study
# Description: Gene Set Enrichment Analysis (GSEA) on the AllHigh_vs_Mixed
#              differential expression result (discovery cohort — the
#              largest cohort, so used here in place of running GSEA
#              separately per comparison/cohort), across four gene set
#              collections: Reactome, KEGG, GO Biological Process, and
#              MSigDB Hallmark.
# Author: Victoria Rios (victoria.riosvazquez@radboudumc.nl)
# Manuscript: "Multi-Omics Clustering Differentiates the Total and Intact HIV
#              Reservoirs and Related Host Immune Mechanisms"
#
# Requires: DEresults object saved by deg_analysis.R (Section 8.2), run with
#           cohort <- "discovery"
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(ReactomePA)
library(data.table)
library(biomaRt)
library(msigdbr)
library(openxlsx)
library(org.Hs.eg.db)

# Set paths — update these to match your local directory structure
# (must match the paths used in deg_analysis.R)
project_dir <- "path/to/project/"
data_dir    <- "path/to/data/"
output_dir  <- "path/to/output/"

# GSEA was run on the discovery cohort only (the largest cohort)
cohort <- "discovery"


# -----------------------------------------------------------------------------
# 1. Load DE Results and Select Comparison
# -----------------------------------------------------------------------------

load(file.path(output_dir, paste0(
  "differential_expression/DESeq2_Results/DESeq2_DEresults_MoClusters_k3_AgeSexSeasonPlateTimetolabCenterRUMCPC1CovidVacc_",
  cohort, ".RData"
)))

# AllHigh vs. Mixed was chosen as the comparison to analyze
results_dea <- DEresults$AllHigh_vs_Mixed@results


# -----------------------------------------------------------------------------
# 2. Map ENSEMBL IDs to ENTREZ IDs and Gene Symbols
# -----------------------------------------------------------------------------

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract Ensembl IDs without version numbers
ensembl_ids <- sub("\\..*", "", rownames(results_dea))

# Query BioMart: Ensembl -> Entrez + HGNC symbol
mapping <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

results_dea <- results_dea %>%
  tibble::rownames_to_column("ENSEMBL") %>%
  mutate(ENSEMBL_clean = sub("\\..*", "", ENSEMBL)) %>%
  left_join(mapping, by = c("ENSEMBL_clean" = "ensembl_gene_id"))


# -----------------------------------------------------------------------------
# 3. Build the Ranked Gene List
# -----------------------------------------------------------------------------
# Genes are ranked by -log10(unadjusted p-value) signed by the direction of
# log2 fold change (unadjusted p gives fewer ties than adjusted p). Genes
# without an Entrez ID are dropped. Where multiple ENSEMBL IDs map to the
# same Entrez ID, the one with the most extreme ranking metric is kept.

GSEA_gene_df <- results_dea %>%
  mutate(GSEA_ranking = -log10(pvalue) * sign(log2FoldChange)) %>%
  filter(!is.na(entrezgene_id)) %>%
  group_by(entrezgene_id) %>%
  slice_max(abs(GSEA_ranking)) %>%
  ungroup() %>%
  arrange(desc(GSEA_ranking))

# Ranked list keyed by Entrez ID (used for Reactome, KEGG, GO)
GSEA_gene_list <- GSEA_gene_df$GSEA_ranking
names(GSEA_gene_list) <- GSEA_gene_df$entrezgene_id

# Ranked list keyed by HGNC symbol (used for MSigDB Hallmark)
GSEA_gene_list_Hallmark <- GSEA_gene_df$GSEA_ranking
names(GSEA_gene_list_Hallmark) <- GSEA_gene_df$hgnc_symbol


# -----------------------------------------------------------------------------
# 4. Plotting and Annotation Helper Functions
# -----------------------------------------------------------------------------

# Extracts the top n most up/down-regulated significant gene sets from a
# GSEA result object for the overview dot plot below.
overview_GSEA_prep <- function(GSEA_result, n_features, up_or_down, P_val_thresh) {
  df <- rbind(GSEA_result@result %>% arrange(desc(NES)) %>% .[1:n_features,],
              GSEA_result@result %>% arrange(NES) %>% .[1:n_features,])
  df %<>% filter(p.adjust < P_val_thresh)
  df$type <- 'upregulated'
  df$type[df$NES < 0] <- 'downregulated'
  if (up_or_down == 'up') {
    return(subset(df, type == 'upregulated'))
  } else if (up_or_down == 'down') {
    return(subset(df, type == 'downregulated'))
  } else if (up_or_down == 'updown') {
    return(df)
  }
}

# Dot plot of normalized enrichment score (NES) vs. gene set, sized by set
# size and colored by significance.
overview_GSEA_plot <- function(dot_df, up_or_down, title) {
  dot_df %<>% mutate(log10P = -log10(p.adjust))
  p <- ggplot(dot_df) +
    theme_bw(base_size = 12) +
    scale_colour_gradient(low = "grey", high = '#A73030') +
    ylab(NULL) +
    ggtitle(title)
  if (up_or_down == 'up') {
    return(p + geom_point(aes(x = NES,
                          y = forcats::fct_reorder(Description, NES),
                          size = setSize,
                          color = log10P)))
  }
  if (up_or_down == 'down') {
    return(p + geom_point(aes(x = NES,
                          y = forcats::fct_reorder(Description, NES, .desc = TRUE),
                          size = setSize,
                          color = log10P))) +
      scale_x_reverse()
  }
  if (up_or_down == 'updown') {
    return(p + geom_point(aes(x = NES,
                          y = forcats::fct_reorder(Description, NES),
                          size = setSize,
                          color = log10P)) +
      facet_grid(.~type))
  }
}

# Combines the two functions above into a single call.
overview_GSEA <- function(GSEA_result, n_features, up_or_down, title, P_val_thresh = 0.05) {
  dot_df <- overview_GSEA_prep(GSEA_result, n_features, up_or_down, P_val_thresh = P_val_thresh)
  plot <- overview_GSEA_plot(dot_df, up_or_down, title)
  return(plot)
}

# Per-gene-set running enrichment score plot (not called by default below;
# kept as a utility for inspecting individual gene sets interactively).
GSEA_plot <- function(GSEA_object, geneset.ID) {
  enrichplot::gseaplot2(GSEA_object,
                        geneSetID = match(geneset.ID, GSEA_object@result[['ID']]),
                        title = GSEA_object@result[['Description']][match(geneset.ID, GSEA_object@result[['ID']])],
                        pvalue_table = TRUE)
}

# Converts a GSEA result's Entrez-ID core_enrichment genes to HGNC symbols,
# collapsed into one comma-separated "genes" column per gene set.
GSEA_res_Entrez_to_SYMBOL <- function(Enrichment_result) {
  res_enrich <- Enrichment_result@result
  res_enrich_long <- res_enrich %>%
    separate_rows(core_enrichment, sep = '/')
  genes_to_symbol <- getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
                           filters = 'entrezgene_id',
                           values = res_enrich_long$core_enrichment,
                           mart = mart) %>%
    mutate(entrezgene_id = as.character(entrezgene_id))
  res_enrich_long %<>% left_join(genes_to_symbol, by = c('core_enrichment' = 'entrezgene_id')) %>%
    distinct(ID, core_enrichment, .keep_all = TRUE) %>%
    dplyr::select(-core_enrichment)

  res_enrich_wide <- res_enrich_long %>%
    group_by(across(-hgnc_symbol)) %>%
    summarise(genes = paste(hgnc_symbol, collapse = ", "), .groups = "drop")
}

# Comparison label used in output filenames below (spaces -> underscores)
comparison_label <- gsub(" ", "_", results_dea$comparison[1])


# -----------------------------------------------------------------------------
# 5. Reactome
# -----------------------------------------------------------------------------

Reactome_GSEA <- gsePathway(geneList = GSEA_gene_list,
                            organism = 'human',
                            pvalueCutoff = 0.50,
                            pAdjustMethod = 'BH',
                            minGSSize = 10,
                            maxGSSize = 500)

overview_GSEA(Reactome_GSEA, up_or_down = 'updown', n_features = 15, P_val_thresh = 0.05,
              title = paste0('GSEA Reactome ', results_dea$comparison[1]))
ggsave(file.path(output_dir, paste0("differential_expression/gsea/GSEA_Reactome_", comparison_label, ".pdf")),
       height = 7, width = 9)

Reactome_GSEA %<>% GSEA_res_Entrez_to_SYMBOL %>%
  arrange(p.adjust)

fwrite(Reactome_GSEA,
       file.path(output_dir, paste0("differential_expression/gsea/GSEA_Reactome_", comparison_label, ".csv")),
       sep = '\t')


# -----------------------------------------------------------------------------
# 6. KEGG
# -----------------------------------------------------------------------------

KEGG_GSEA <- gseKEGG(geneList = GSEA_gene_list,
                     organism = 'human',
                     pvalueCutoff = 0.50,
                     pAdjustMethod = 'BH',
                     minGSSize = 10,
                     maxGSSize = 500)

overview_GSEA(KEGG_GSEA, up_or_down = 'updown', n_features = 15, P_val_thresh = 0.05,
              title = paste0('GSEA KEGG ', results_dea$comparison[1]))
ggsave(file.path(output_dir, paste0("differential_expression/gsea/GSEA_KEGG_", comparison_label, ".pdf")),
       height = 7, width = 9)

KEGG_GSEA %<>% GSEA_res_Entrez_to_SYMBOL %>%
  arrange(p.adjust)

fwrite(KEGG_GSEA,
       file.path(output_dir, paste0("differential_expression/gsea/GSEA_KEGG_", comparison_label, ".csv")),
       sep = '\t')


# -----------------------------------------------------------------------------
# 7. GO (Biological Process)
# -----------------------------------------------------------------------------

GO_GSEA <- gseGO(geneList = GSEA_gene_list,
                 ont = 'BP',
                 OrgDb = org.Hs.eg.db,
                 pvalueCutoff = 0.50,
                 pAdjustMethod = 'BH',
                 minGSSize = 10,
                 maxGSSize = 500)

overview_GSEA(GO_GSEA, up_or_down = 'updown', n_features = 15, P_val_thresh = 0.05,
              title = paste0('GSEA GO BP ', results_dea$comparison[1]))
ggsave(file.path(output_dir, paste0("differential_expression/gsea/GSEA_GOBP_", comparison_label, ".pdf")),
       height = 7, width = 9)

GO_GSEA %<>% GSEA_res_Entrez_to_SYMBOL %>%
  arrange(p.adjust)

fwrite(GO_GSEA,
       file.path(output_dir, paste0("differential_expression/gsea/GSEA_GOBP_", comparison_label, ".csv")),
       sep = '\t')


# -----------------------------------------------------------------------------
# 8. MSigDB Hallmark
# -----------------------------------------------------------------------------

hallmark_gene_sets <- msigdbr(species = 'Homo sapiens', category = "H")
TERM2GENE_Hallmark <- hallmark_gene_sets %>% dplyr::select(gs_name, human_gene_symbol)
TERM2NAME_Hallmark <- hallmark_gene_sets %>% dplyr::select(gs_name, gs_description)

Hallmark_GSEA <- clusterProfiler::GSEA(geneList = GSEA_gene_list_Hallmark,
                                        TERM2GENE = TERM2GENE_Hallmark,
                                        TERM2NAME = TERM2NAME_Hallmark,
                                        pvalueCutoff = 0.50,
                                        pAdjustMethod = 'BH',
                                        minGSSize = 10,
                                        maxGSSize = 500)

overview_GSEA(Hallmark_GSEA, up_or_down = 'updown', n_features = 15, P_val_thresh = 0.05,
              title = paste0('GSEA Hallmark ', results_dea$comparison[1]))
ggsave(file.path(output_dir, paste0("differential_expression/gsea/GSEA_Hallmark_", comparison_label, ".pdf")),
       height = 7, width = 9)

Hallmark_GSEA %<>% GSEA_res_Entrez_to_SYMBOL %>%
  arrange(p.adjust)

fwrite(Hallmark_GSEA,
       file.path(output_dir, paste0("differential_expression/gsea/GSEA_Hallmark_", comparison_label, ".csv")),
       sep = '\t')


# -----------------------------------------------------------------------------
# 9. Export Combined Results (Excel)
# -----------------------------------------------------------------------------

wb <- createWorkbook()

# Adds a gene set enrichment table as a worksheet, with genes collapsed
# into one semicolon-separated cell per gene set.
add_GSEA_sheet <- function(wb, df, sheetname) {
  df$genes <- gsub(", ", "; ", df$genes)
  addWorksheet(wb, sheetname)
  writeData(wb, sheetname, df)
}

add_GSEA_sheet(wb, Reactome_GSEA, "Reactome")
add_GSEA_sheet(wb, KEGG_GSEA,     "KEGG")
add_GSEA_sheet(wb, GO_GSEA,       "GO_BP")
add_GSEA_sheet(wb, Hallmark_GSEA, "Hallmark")

saveWorkbook(wb,
             file = file.path(output_dir, paste0(
               "differential_expression/gsea/GSEA_results_", comparison_label, ".xlsx"
             )),
             overwrite = TRUE)
