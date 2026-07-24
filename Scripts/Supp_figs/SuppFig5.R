library(tidyverse)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)

gene_count <- read.csv("../../Data/GSE157103/merged_genes.csv", row.names=1)
norm_gene_count <- read.csv("../../Data/GSE157103/normalized_genes.csv", row.names=1)
erv_count <- read.csv("../../Data/GSE157103/merged_ervs.csv", row.names=1)
coldata <- read.csv("../../Data/GSE157103/metadata.csv", row.names=1)
norm_erv_count <- read.csv("../../Data/GSE157103/normalized_ervs.csv", row.names=1)

# COVID-19 vs nonCOVID-19 ICU patients 
dir.create("COVID_vs_nonCOVID_icu", showWarnings = FALSE)
icu <- coldata
icu = icu %>% filter(icu == 'yes')

erv_count <- erv_count[colnames(erv_count) %in% rownames(icu)]
norm_erv_count <- norm_erv_count[colnames(norm_erv_count) %in% rownames(icu)]
gene_count <- gene_count[colnames(gene_count) %in% rownames(icu)]
norm_gene_count <- norm_gene_count[colnames(norm_gene_count) %in% rownames(icu)]

# Differential expression analysis
output_name <- "./COVID_vs_nonCOVID_icu/DESeqResults_ERVs.txt"
source("../wrapper_functions/DESeq2Analysis.r", echo=TRUE, max.deparse.length=10e3)
run_DESeq2(countData = erv_count, 
           coldata = icu, 
           keep_limit = 0, 
           output_name = output_name,
           conditions = c("condition", "COVID-19", "non-COVID-19"))
DESeqResults_erv <- read.table("./COVID_vs_nonCOVID_icu/DESeqResults_ERVs.txt", header=TRUE)
DESeqResults_erv <- DESeqResults_erv[DESeqResults_erv$baseMean >= 10,]

write.table(DESeqResults_erv, file="./COVID_vs_nonCOVID_icu/DESeqResults_ERVs_filtered.txt")

# Isolate up/down-regulated genes
source("../wrapper_functions/IsolateGenes.R", echo=TRUE, max.deparse.length=10e3)
output_dir <- "./COVID_vs_nonCOVID_icu"
run_IsolateGenes(res = DESeqResults_erv, 
                 l2FC_limit = 1, 
                 padj_limit = 0.05, 
                 output_dir = output_dir,
                 type = "erv")

# Define signature ERVs
erv_overlap_feature_summary <- read.delim("../../Assets/erv_overlap_summary.tsv", stringsAsFactors = FALSE)
upreg_erv <- scan("./COVID_vs_nonCOVID_icu/upreg_erv.txt", what = "", sep = ",", quiet = TRUE)
downreg_erv <- scan("./COVID_vs_nonCOVID_icu/downreg_erv.txt", what = "", sep = ",", quiet = TRUE)

upreg_erv_overlap_feature_summary <- erv_overlap_feature_summary[erv_overlap_feature_summary$erv_name %in% upreg_erv, ]
downreg_erv_overlap_feature_summary <- erv_overlap_feature_summary[erv_overlap_feature_summary$erv_name %in% downreg_erv, ]
DE_erv_overlap_summary <- rbind(upreg_erv_overlap_feature_summary, downreg_erv_overlap_feature_summary)
DE_ERVs <- DE_erv_overlap_summary$erv_name

# =====  ERVs that have overlapping genomic features (intron retention?)  =====
overlapped_ervs <- DE_erv_overlap_summary %>% filter (strand_overlap %in% c("same", "both"))
## Expand rows if multiple gene_ids per ERV
erv_pairs <- overlapped_ervs[, c("erv_name", "overlapping_gene_ids")]
erv_pairs <- separate_rows(erv_pairs, overlapping_gene_ids, sep = ",")

## Perform Spearman correlation for each ERV-gene pair
### Initialize a result df
cor_results <- data.frame(
  erv = character(),
  gene = character(),
  spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
### Loop through each pair
for (i in seq_len(nrow(erv_pairs))) {
  erv <- erv_pairs$erv_name[i]
  gene <- erv_pairs$overlapping_gene_ids[i]
  
  ## Make sure both are present in the count matrices
  if (erv %in% rownames(norm_erv_count) && gene %in% rownames(norm_gene_count)) {
    erv_expr <- as.numeric(norm_erv_count[erv, ])
    gene_expr <- as.numeric(norm_gene_count[gene, ])
    
    cor_test <- cor.test(erv_expr, gene_expr, method = "spearman")
    
    cor_results <- rbind(cor_results, data.frame(
      erv = erv,
      gene = gene,
      spearman_rho = cor_test$estimate,
      p_value = cor_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
cor_results$padj <- p.adjust(cor_results$p_value, method = "BH")

## Summarize the correlations
write.csv(cor_results, "./COVID_vs_nonCOVID_icu/correlation_sig_erv_and_overlap_features.csv", row.names = FALSE)

## Filter for significant (padj <0.05 and |rho| >0.5) correlations (i.e. ERVs to be excluded)
cor_results_no_na <- na.omit(cor_results)
potential_IR <- cor_results_no_na[cor_results_no_na$padj < 0.05 & (cor_results_no_na$spearman_rho > 0.5 | cor_results_no_na$spearman_rho < -0.5), ]
write.csv(potential_IR, "./COVID_vs_nonCOVID_icu/correlation_de_erv_and_overlap_features_significant.csv", row.names = FALSE)

## Filter out ERVs that had significant correlation with its overlap
signature_ervs <- DE_erv_overlap_summary %>%
  filter(!erv_name %in% potential_IR$erv)

# =====  ERVs w/o overlapping features, but have upstream protein-coding genes (Transcriptional readthrough?) =====
## Create ERV–upstreamGene pairs
erv_RT_pairs <- DE_erv_overlap_summary %>%
  filter(closest_upstream_pc_gene != ".") %>% 
  dplyr::select(erv_name, closest_upstream_pc_gene)

erv_RT_pairs <- erv_RT_pairs %>%
  mutate(
    closest_upstream_pc_gene = recode(
      closest_upstream_pc_gene,
      "RNF213" = "ENSG00000173821",
      "CMPK2" = "ENSG00000134326"
    )
  )

## Compute Spearman Correlation for every ERV-gene pair
### Initialize results df
cor_results_ERV_upstreamGene <- data.frame(
  erv_name = character(),
  gene = character(),
  spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

### Loop through each ERV–gene pair
for (i in seq_len(nrow(erv_RT_pairs))) {
  erv_name <- erv_RT_pairs$erv_name[i]
  gene <- erv_RT_pairs$closest_upstream_pc_gene[i]  
  # Check to make sure both exist in count matrices
  erv_expr <- as.numeric(norm_erv_count[erv_name, ])
  gene_expr <- as.numeric(norm_gene_count[gene, ])
  
  cor_test <- cor.test(erv_expr, gene_expr, method = "spearman")
  
  cor_results_ERV_upstreamGene <- rbind(cor_results_ERV_upstreamGene, data.frame(
    erv_name = erv_name,
    gene = gene,
    spearman_rho = cor_test$estimate,
    p_value = cor_test$p.value,
    stringsAsFactors = FALSE
  ))
}
cor_results_ERV_upstreamGene$padj <- p.adjust(cor_results_ERV_upstreamGene$p_value, method = "BH")
## Summarize the correlations
write.csv(cor_results_ERV_upstreamGene, "./COVID_vs_nonCOVID_icu/ERVs_nearby_genes_spearman_cor.csv", row.names = FALSE)
## filter for p value <0.05 and |rho| >0.5 (potential RT)
potential_readthrough <- cor_results_ERV_upstreamGene[cor_results_ERV_upstreamGene$padj < 0.05 & (cor_results_ERV_upstreamGene$spearman_rho > 0.5 | cor_results_ERV_upstreamGene$spearman_rho < -0.5), ]
# exclude potential IR from signature ERVs
signature_ervs <- signature_ervs %>%
  filter(!erv_name %in% potential_readthrough$erv)

# Volcano Plot
output_name <- "./COVID_vs_nonCOVID_icu/figures/VolcanoPlot_ERV.pdf"
# Highlight ERV signatures on the volcano plot
select_labels <- signature_ervs$erv_name
source("../wrapper_functions/VolcanoPlot.R", echo=TRUE, max.deparse.length=10e3)
run_VolcanoPlot(res = DESeqResults_erv, 
                select_labels = select_labels,
                output_name = output_name) 
