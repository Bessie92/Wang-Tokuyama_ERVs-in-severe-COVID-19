library(reshape2)
library(ggpubr)
library(dplyr)
library(ggplot2)
# Figure 2A: Visualize dysregulated epigenetic regulators in COVID-19
## Read in gene count table with gene names 
norm_gene_count_name <- read.csv("../../Data/GSE171110/normalized_genes_with_name.csv", row.names=1)
coldata <- read.csv("../../Data/GSE171110/metadata.csv", row.names=1)

## Subset and melt data
epigenetic_reulators <- c("PPHLN1", "TASOR", "MPHOSPH8", "MORC2", 
                          "TET1", "TET2", "TET3", 
                          "TRIM28", "SETDB1", "DNMT1", "DNMT3A", "DNMT3B",
                          'PADI1', 'PADI2', 'PADI3', 'PADI4', 'PADI6')
epigenetic_regulators_count <- norm_gene_count_name[epigenetic_reulators, ]

## Remove genes with 0 expression across all samples (otherwise wilcox won't perform)
epigenetic_regulators_count <- epigenetic_regulators_count[rowSums(epigenetic_regulators_count == 0) < ncol(epigenetic_regulators_count), ]

epigenetic_regulators_count <- as.data.frame(t(epigenetic_regulators_count))
epigenetic_regulators_count$sample <- rownames(epigenetic_regulators_count)
epigenetic_regulators_count$condition <- coldata[epigenetic_regulators_count$sample, "condition"]
expr_long <- melt(epigenetic_regulators_count, id.vars = c("sample", "condition"), variable.name = "gene", value.name = "expression")
expr_long$condition <- factor(expr_long$condition, levels = c("Healthy", "COVID-19"))

## Box-and-Whisker Plot
ggplot(expr_long, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, width = 0.5) +
  scale_fill_manual(values = c("Healthy" = "#26547d", "COVID-19" = "#ec476f")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
  facet_wrap(~ gene, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Normalized Expression") +
  xlab("")
ggsave("./figures/Box-and-Whisker Plot - epigenetic regulator between condition.pdf", width = 10, height = 9)

# Figure 2B: Spearman correlation between ERVs and epigenetic regulators in severe COVID-19 patients
norm_erv_count <- read.csv("../../Data/GSE171110/normalized_ervs.csv",row.names = 1)

## Signature ERVs in plotting order
sig_ervs <- c(
  "909", "4745", "1179", "2724", "3052", "5302","3743","3659","4426", "942", 
  "814", "882", "1879", 
  "2095", "W-8", 
  "3687", "4830", 
  "2124", "2637", 
  "1379", "4174", 
  "606", "1874", 
  "4896", "565", "3104", "3508", "4657", "6080", "1043", 
  "6217", "4352","5208"
)

## Epigenetic regulators in plotting order
epigenetic_regulators_significant <- c(
  "PPHLN1", "MORC2",
  "TET2", "TET3",
  "SETDB1", "DNMT3B",
  "PADI2", "PADI4"
)


## Pull ERV and gene counts separately
ervs <- norm_erv_count[sig_ervs, , drop = FALSE]
genes <- norm_gene_count_name[epigenetic_regulators_significant, , drop = FALSE]

## Filter for COVID-19 patients only
covid_samples <- rownames(coldata)[coldata$condition == "COVID-19"]
ervs <- ervs[, covid_samples, drop = FALSE]
genes <- genes[, covid_samples, drop = FALSE]

## Confirm sample order is identical
stopifnot(identical(colnames(ervs), colnames(genes)))

## Reorder gene matrix samples to match the ERV matrix
genes <- genes[, colnames(ervs), drop = FALSE]

## Preserve the feature order after filtering
erv_order <- rownames(ervs)
gene_order <- rownames(genes)

## Run ERV-gene correlations
results <- list()

for (erv in erv_order) {
  for (gene in gene_order) {
    
    erv_expr <- as.numeric(ervs[erv, ])
    gene_expr <- as.numeric(genes[gene, ])
    
    ## Keep samples with finite values in both vectors
    valid <- is.finite(erv_expr) & is.finite(gene_expr)
    
    ## Skip correlations when there are too few samples
    ## or either feature has no variation
    if (
      sum(valid) < 3 ||
      length(unique(erv_expr[valid])) < 2 ||
      length(unique(gene_expr[valid])) < 2
    ) {
      next
    }
    
    rho_test <- cor.test(
      erv_expr[valid],
      gene_expr[valid],
      method = "spearman",
      exact = FALSE
    )
    
    results[[length(results) + 1]] <- data.frame(
      ERV = erv,
      Gene = gene,
      Rho = unname(rho_test$estimate),
      p_value = rho_test$p.value,
      stringsAsFactors = FALSE
    )
  }
}

## Combine results and adjust across all ERV-gene tests
cor_df <- bind_rows(results) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

## Retain significant correlations for plotting
significant_cor_df <- cor_df %>%
  filter(p_adj < 0.05)

## Create complete ERV-gene grid so ERVs without significant
## correlations remain on the axes
full_grid <- expand.grid(
  ERV = erv_order,
  Gene = gene_order,
  stringsAsFactors = FALSE
)

full_cor_df <- full_grid %>%
  left_join(
    significant_cor_df,
    by = c("ERV", "Gene")
  ) %>%
  mutate(
    Gene = factor(Gene, levels = gene_order),
    ERV = factor(ERV, levels = rev(erv_order))
  )

## Plot
p <- ggplot(full_cor_df, aes(x = Gene, y = ERV)) +
  geom_point(
    aes(
      color = Rho,
      size = -log10(p_adj)
    ),
    na.rm = TRUE
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  labs(
    title = "Significant ERV–Gene Spearman Correlations",
    x = "Epigenetic regulators",
    y = "ERV",
    color = "Spearman rho",
    size = "-log10 adjusted p-value"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p

ggsave(
  "./figures/Dotplot_ERV_epigenetic_regulator_significant_only.pdf",
  plot = p,
  width = 9,
  height = 7
)
