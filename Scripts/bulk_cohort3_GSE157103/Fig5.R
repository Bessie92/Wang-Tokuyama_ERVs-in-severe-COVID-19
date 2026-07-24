library(tidyverse)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggrepel)

gene_count <- read.csv("../../Data/GSE157103/merged_genes.csv", row.names=1)

# Normalization
library(DESeq2)
conds <- factor( c(1:ncol(gene_count)))
dds <- DESeqDataSetFromMatrix(countData = as.matrix(gene_count), colData = as.matrix(conds), design = ~ 1)
dds <- estimateSizeFactors(dds)
normalizedCounts <- counts(dds, normalized=T)
write.table(normalizedCounts,
            file = "../../Data/GSE157103/normalized_genes.csv",
            sep = ",",
            quote = FALSE,
            col.names = NA)
write.table(sizeFactors(dds), "../../Data/GSE157103/normalized_factors.txt", sep="\t", quote=FALSE, col.names=NA)
norm_gene_count <- read.csv("../../Data/GSE157103/normalized_genes.csv", row.names=1)
size_factor <- read.table("../../Data/GSE157103/normalized_factors.txt")
colnames(size_factor) <- c("sizeFactor")
erv_count <- read.csv("../../Data/GSE157103/merged_ervs.csv", row.names=1)
coldata <- read.csv("../../Data/GSE157103/metadata.csv", row.names=1)

norm_erv_count <- erv_count
for (sampleName in colnames(erv_count)) {
  if (sampleName %in% rownames(size_factor)) {
    sizeFactor <- size_factor[sampleName, ]
    norm_erv_count[, sampleName] <- erv_count[, sampleName] / sizeFactor
  } else {
    warning(paste("Size factor for", sampleName, "not found!"))
  }
}

write.table(norm_erv_count, "../../Data/GSE157103/normalized_ervs.csv", sep = ",", quote=F, col.names = NA)

dir.create("COVID19_ICU_vs_nonICU", showWarnings = FALSE)
covid <- coldata
covid = covid %>% filter(condition == 'COVID-19')

erv_count <- erv_count[colnames(erv_count) %in% rownames(covid)]
norm_erv_count <- norm_erv_count[colnames(norm_erv_count) %in% rownames(covid)]
gene_count <- gene_count[colnames(gene_count) %in% rownames(covid)]
norm_gene_count <- norm_gene_count[colnames(norm_gene_count) %in% rownames(covid)]

# Differential expression analysis
output_name <- "./COVID19_ICU_vs_nonICU/DESeqResults_ERVs.txt"
source("../wrapper_functions/DESeq2Analysis.r", echo=TRUE, max.deparse.length=10e3)
run_DESeq2(countData = erv_count, 
           coldata = covid, 
           keep_limit = 0, 
           output_name = output_name,
           conditions = c("icu", "yes", "no"))


DESeqResults_erv <- read.table("./COVID19_ICU_vs_nonICU/DESeqResults_ERVs.txt", header=TRUE)
DESeqResults_erv <- DESeqResults_erv[DESeqResults_erv$baseMean >= 10,]
write.table(DESeqResults_erv, file="./COVID19_ICU_vs_nonICU/DESeqResults_ERVs_filtered.txt")

# Isolate up/down-regulated ERVs
source("../wrapper_functions/IsolateGenes.R", echo=TRUE, max.deparse.length=10e3)
output_dir <- "../../Data/GSE157103"
run_IsolateGenes(res = DESeqResults_erv, 
                 l2FC_limit = 1, 
                 padj_limit = 0.05, 
                 output_dir = output_dir,
                 type = "erv")

# Define signature ERVs
dir.create("../../Data/GSE157103/define_signature_ervs")
erv_overlap_feature_summary <- read.delim("../../Assets/erv_overlap_summary.tsv", stringsAsFactors = FALSE)
upreg_erv <- scan("../../Data/GSE157103/upreg_erv.txt", what = "", sep = ",", quiet = TRUE)
downreg_erv <- scan("../../Data/GSE157103/downreg_erv.txt", what = "", sep = ",", quiet = TRUE)

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
write.csv(cor_results, "../../Data/GSE157103/define_signature_ervs/correlation_sig_erv_and_overlap_features.csv", row.names = FALSE)
#cor_results <- read.csv("correlation_sig_erv_and_overlap_features.csv")

## Filter for significant (padj <0.05 and |rho| >0.5) correlations (i.e. ERVs to be excluded)
cor_results_no_na <- na.omit(cor_results)
potential_IR <- cor_results_no_na[cor_results_no_na$padj < 0.05 & (cor_results_no_na$spearman_rho > 0.5 | cor_results_no_na$spearman_rho < -0.5), ]
write.csv(potential_IR, "../../Data/GSE157103/define_signature_ervs/correlation_de_erv_and_overlap_features_significant.csv", row.names = FALSE)

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
      "CMPK2" = "ENSG00000134326",
      "HERC6" = "ENSG00000138642"
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
write.csv(cor_results_ERV_upstreamGene, "../../Data/GSE157103/define_signature_ervs/ERVs_nearby_genes_spearman_cor.csv", row.names = FALSE)
## filter for p value <0.05 and |rho| >0.5 (potential RT)
potential_readthrough <- cor_results_ERV_upstreamGene[cor_results_ERV_upstreamGene$padj < 0.05 & (cor_results_ERV_upstreamGene$spearman_rho > 0.5 | cor_results_ERV_upstreamGene$spearman_rho < -0.5), ]
# exclude potential IR from signature ERVs
signature_ervs <- signature_ervs %>%
  filter(!erv_name %in% potential_readthrough$erv)

# Figure 5A: Volcano Plot
output_name <- "./COVID19_ICU_vs_nonICU/VolcanoPlot_ERV.pdf"
# Highlight ERV signatures on the volcano plot
select_labels <- signature_ervs$erv_name
source("../wrapper_functions/VolcanoPlot.R", echo=TRUE, max.deparse.length=10e3)
run_VolcanoPlot(res = DESeqResults_erv, 
                select_labels = select_labels,
                output_name = output_name) 
# Figure 5B
erv_all  <- c("2124", "3743", "1874", "5302", "4830", "2724")

# Helper function to make ranked based spearman cor visualization
plot_rank_spearman <- function(df, x_col, y_col, title, outfile) {
  # Rank
  df <- df %>%
    dplyr::mutate(
      RankX = rank(.data[[x_col]], ties.method = "average"),
      RankY = rank(.data[[y_col]], ties.method = "average")
    )
  
  # Spearman correlation (Pearson on ranks)
  spearman_test <- cor.test(
    df$RankX,
    df$RankY,
    method = "pearson"
  )
  
  rho <- spearman_test$estimate
  p   <- spearman_test$p.value
  
  p_label <- ifelse(p < 1e-3,
                    "p < 0.001",
                    paste0("p = ", signif(p, 3)))
  
  # Plot
  p_out <- ggplot(df, aes(x = RankX, y = RankY)) +
    geom_point(size = 4, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(
      title    = title,
      subtitle = paste0("Spearman ρ = ", round(rho, 3), ", ", p_label),
      x = "Rank (ERV expression)",
      y = "Rank (HFD.45)"
    ) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(outfile, p_out, width = 6, height = 6)
  
  return(invisible(list(rho = rho, p = p)))
}

cor_data <- covid %>%
  rownames_to_column("sample") %>%
  left_join(
    norm_erv_count %>% 
      t() %>%
      as.data.frame() %>%
      rownames_to_column("sample"),
    by = "sample"
  )
# Sum of six ERVs
cor_data <- cor_data %>%
  dplyr::mutate(
    Sum_ERV_6 = rowSums(across(all_of(erv_all)))
  )

plot_rank_spearman(
  df      = cor_data,
  x_col  = "Sum_ERV_6",
  y_col  = "HFD.45",
  title  = "Sum of 6 ERVs vs HFD.45",
  outfile = "./COVID19_ICU_vs_nonICU/Sum_6_ERVs_vs_HFD45_rank_spearman.pdf"
)

# Figure 5C - Individual ERVs
for (erv in erv_all) {
  
  plot_rank_spearman(
    df      = data_base,
    x_col  = erv,
    y_col  = "HFD.45",
    title  = paste0("ERV ", erv, " vs HFD.45"),
    outfile = paste0("./COVID19_ICU_vs_nonICU/ERV_", erv, "_vs_HFD45_rank_spearman.pdf")
  )
}
