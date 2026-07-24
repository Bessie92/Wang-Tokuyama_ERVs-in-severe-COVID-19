library(tidyr)

erv_overlap_feature_summary <- read.delim("../../Assets/erv_overlap_summary.tsv", stringsAsFactors = FALSE)
upreg_erv <- scan("../../Data/GSE171110/upreg_erv.txt", what = "", sep = ",", quiet = TRUE)
downreg_erv <- scan("../../Data/GSE171110/downreg_erv.txt", what = "", sep = ",", quiet = TRUE)

norm_gene_count <- read.csv("../../Data/GSE171110/normalized_genes.csv", row.names=1)
norm_gene_count_name <- read.csv("../../Data/GSE171110/normalized_genes_with_name.csv", row.names=1)
norm_erv_count <- read.csv("../../Data/GSE171110/normalized_ervs.csv", row.names=1)

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
write.csv(cor_results, "../../Data/GSE171110/define_signature_ervs/correlation_sig_erv_and_overlap_features.csv", row.names = FALSE)
#cor_results <- read.csv("correlation_sig_erv_and_overlap_features.csv")

## Filter for significant (padj <0.05 and |rho| >0.5) correlations (i.e. ERVs to be excluded)
cor_results_no_na <- na.omit(cor_results)
potential_IR <- cor_results_no_na[cor_results_no_na$padj < 0.05 & (cor_results_no_na$spearman_rho > 0.5 | cor_results_no_na$spearman_rho < -0.5), ]
write.csv(potential_IR, "../../Data/GSE171110/define_signature_ervs/correlation_de_erv_and_overlap_features_significant.csv", row.names = FALSE)

## Filter out ERVs that had significant correlation with its overlap
signature_ervs <- DE_erv_overlap_summary %>%
  filter(!erv_name %in% potential_IR$erv)

# =====  ERVs w/o overlapping features, but have upstream protein-coding genes (Transcriptional readthrough?) =====
## Create ERV–upstreamGene pairs
erv_RT_pairs <- DE_erv_overlap_summary %>%
  filter(closest_upstream_pc_gene != ".") %>% 
  dplyr::select(erv_name, closest_upstream_pc_gene)

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
    gene_expr <- as.numeric(norm_gene_count_name[gene, ])
    
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
write.csv(cor_results_ERV_upstreamGene, "../../Data/GSE171110/define_signature_ervs/ERVs_nearby_genes_spearman_cor.csv", row.names = FALSE)
## filter for p value <0.05 and |rho| >0.5 (potential RT)
potential_readthrough <- cor_results_ERV_upstreamGene[cor_results_ERV_upstreamGene$padj < 0.05 & (cor_results_ERV_upstreamGene$spearman_rho > 0.5 | cor_results_ERV_upstreamGene$spearman_rho < -0.5), ]
# exclude potential IR from signature ERVs
signature_ervs <- signature_ervs %>%
  filter(!erv_name %in% potential_readthrough$erv)

write.csv(signature_ervs, "../../Data/GSE171110/define_signature_ervs/Raw_Table1_Signature_erv_overlap_summary.csv")
