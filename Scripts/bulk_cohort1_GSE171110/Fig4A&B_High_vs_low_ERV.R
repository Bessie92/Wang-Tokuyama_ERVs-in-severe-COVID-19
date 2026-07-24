library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(grid)

# Inputs 
gene_count            <- read.csv("../../Data/GSE171110/merged_genes.csv", row.names = 1)
norm_gene_count       <- read.csv("../../Data/GSE171110/normalized_genes.csv", row.names = 1)
erv_count             <- read.csv("../../Data/GSE171110/merged_ervs.csv", row.names = 1)
norm_erv_count        <- read.csv("../../Data/GSE171110/normalized_ervs.csv", row.names = 1)
norm_gene_count_name  <- read.csv("../../Data/GSE171110/normalized_genes_with_name.csv", row.names = 1)
coldata               <- read.csv("../../Data/GSE171110/metadata.csv", row.names = 1)
sig_ervs <- read.csv("../../Data/GSE171110/define_signature_ervs/Raw_Table1_Signature_erv_overlap_summary.csv")$erv_name
sig_ervs_upreg <- sig_ervs[!sig_ervs %in% c("6217", "5208", "4352")]
reactome_top5 <- read.csv(
  "./WGCNA/Reactome/reactome_top_5_per_ME.csv"
)

# Helper functions
subset_samples <- function(mat, samples) {
  mat[, colnames(mat) %in% samples, drop = FALSE]
}

subset_features <- function(mat, features) {
  as.matrix(mat[rownames(mat) %in% features, , drop = FALSE])
}

to_long_df <- function(mat, group_name, feature_col = "erv", value_col = "expr") {
  mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(feature_col) %>%
    pivot_longer(cols = -all_of(feature_col), names_to = "sample", values_to = value_col) %>%
    mutate(group = group_name)
}

outdir <- "./COVID_subsets_based_on_ervs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Define sample groups
coldata_covid   <- coldata %>% filter(condition == "COVID-19")
low_erv_covid <- c(
  "SRX10470994", "SRX10470987", "SRX10470991", "SRX10470998",
  "SRX10470981", "SRX10470997", "SRX10470990", "SRX10471000"
)
high_erv_covid <- setdiff(rownames(coldata_covid), low_erv_covid)

coldata_covid <- coldata_covid %>%
  mutate(
    ERV_subset = ifelse(rownames(.) %in% low_erv_covid, "Low_ERV", "High_ERV")
  )

# Subset matrices by group
norm_erv_count_high_erv <- subset_samples(norm_erv_count, high_erv_covid)
norm_erv_count_low_erv  <- subset_samples(norm_erv_count, low_erv_covid)

norm_gene_count_name_high_erv <- subset_samples(norm_gene_count_name, high_erv_covid)
norm_gene_count_name_low_erv  <- subset_samples(norm_gene_count_name, low_erv_covid)

# === Fig 4A: Compare average signature ERV expression === 
selected_high <- subset_features(norm_erv_count_high_erv, sig_ervs_upreg)
selected_low  <- subset_features(norm_erv_count_low_erv, sig_ervs_upreg)

df_all <- bind_rows(
  to_long_df(selected_low,  "Low ERV"),
  to_long_df(selected_high, "High ERV")
)

stats_out <- df_all %>%
  group_by(erv) %>%
  summarize(
    p = wilcox.test(expr ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

mean_low_per_sample  <- colMeans(selected_low,  na.rm = TRUE)
mean_high_per_sample <- colMeans(selected_high, na.rm = TRUE)

df_global_sample <- data.frame(
  sample = c(names(mean_low_per_sample), names(mean_high_per_sample)),
  group = c(
    rep("Low ERV",  length(mean_low_per_sample)),
    rep("High ERV", length(mean_high_per_sample))
  ),
  mean_expr = c(mean_low_per_sample, mean_high_per_sample)
) %>%
  mutate(group = factor(group, levels = c("Low ERV", "High ERV")))

p_box <- ggplot(df_global_sample, aes(x = group, y = mean_expr, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  scale_fill_manual(values = c("Low ERV" = "#f26838", "High ERV" = "#991d20")) +
  theme_bw(base_size = 14) +
  labs(
    x = "",
    y = "Mean ERV expression per sample"
  )

p_box
ggsave(file.path(outdir, "Mean_sig_ERV_comparison.pdf"), p_box)

# === Fig 4B: Volcano between erv high and low with reactome gene overlayed  ===
gene_count_covid <- subset_samples(gene_count, rownames(coldata_covid))
# DESeq2 on cellular genes
output_name <- file.path(outdir, "DESeqResults_cellular_genes_ERV_subset.txt")
source("../wrapper_functions/DESeq2Analysis.r", echo = TRUE, max.deparse.length = 10e3)

run_DESeq2(
  countData   = gene_count_covid,
  coldata     = coldata_covid,
  keep_limit  = 0,
  output_name = output_name,
  conditions  = c("ERV_subset", "High_ERV", "Low_ERV")
)

DESeqResults_gene_covid_subset <- read.table(output_name, header = TRUE) %>%
  filter(baseMean >= 10)

write.table(
  DESeqResults_gene_covid_subset,
  file = file.path(outdir, "DESeqResults_cellular_genes_ERV_subset_filtered.txt")
)

# Build Reactome gene-term map
reactome_gene_map <- reactome_top5 %>%
  arrange(p.adjust, Description) %>%
  dplyr::select(Description, geneID, p.adjust) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(gene = geneID) %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene, reactome_term = Description)

# Convert DESeq ENSEMBL IDs to gene symbols
de_df <- DESeqResults_gene_covid_subset %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  mutate(
    ensembl_gene_id = sub("\\..*$", "", ensembl_gene_id)
  )

de_df$gene <- mapIds(
  org.Hs.eg.db,
  keys = de_df$ensembl_gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

de_annot <- de_df %>%
  left_join(reactome_gene_map, by = "gene") %>%
  mutate(
    plot_label = ifelse(is.na(gene), ensembl_gene_id, gene),
    sig = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
    padj_for_plot  = ifelse(is.na(padj) | padj <= 0, NA, padj),
    neg_log10_padj = -log10(padj_for_plot)
  )

reactome_terms <- sort(unique(de_annot$reactome_term[!is.na(de_annot$reactome_term) & de_annot$sig]))

de_annot <- de_annot %>%
  mutate(
    category = case_when(
      !sig ~ "non_sig",
      sig & is.na(reactome_term) ~ "other_sig",
      sig & !is.na(reactome_term) ~ reactome_term
    ),
    alpha_group = case_when(
      category %in% reactome_terms ~ 1,
      TRUE ~ 0.6
    ),
    category = factor(category, levels = c("non_sig", "other_sig", reactome_terms)),
    is_reactome = category %in% reactome_terms
  ) %>%
  arrange(is_reactome)

n_terms <- length(reactome_terms)
term_cols <- if (n_terms > 0) {
  brewer.pal(min(max(n_terms, 3), 12), "Paired")[seq_len(n_terms)]
} else {
  character(0)
}
names(term_cols) <- reactome_terms

color_values <- c(
  non_sig   = "grey40",
  other_sig = "grey75",
  term_cols
)

p_volcano <- ggplot(
  de_annot,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = category,
    alpha = alpha_group
  )
) +
  geom_point(size = 1.5, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = de_annot %>% filter(category %in% reactome_terms),
    aes(label = plot_label),
    size = 4,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.15,
    segment.color = "grey40",
    segment.size = 0.3
  ) +
  scale_color_manual(values = color_values, name = "Reactome term") +
  scale_alpha_identity() +
  labs(
    x = "log2 fold change",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.7, "lines"),
    plot.title = element_text(face = "bold")
  )

p_volcano

ggsave(
  file.path(outdir, "Volcano_DESeq2_with_Reactome_terms_sig_vs_nonsig_cutoffs.pdf"),
  p_volcano,
  width = 15,
  height = 7
)

