library(dplyr)
library(ggplot2)
library(ggpattern)

reactome_dir <- "./WGCNA/Reactome"

module_order_plot <- c(
  "ME5", "ME48", "ME2", "ME61", "ME9", "ME0"
)

# Read previously saved gene–ERV correlations
cor_df <- read.csv(
  file.path(
    reactome_dir,
    "Reactome_gene_ERV_spearman_correlations_mergedTerms_pruned.csv"
  ),
  stringsAsFactors = FALSE
)

# Check required columns
required_cols <- c(
  "module", "gene", "erv", "rho", "padj"
)

missing_cols <- setdiff(required_cols, colnames(cor_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Count significant ERV correlations per gene
gene_sig_counts <- cor_df %>%
  filter(
    module %in% module_order_plot,
    !is.na(padj),
    padj < 0.05,
    !is.na(rho),
    rho != 0
  ) %>%
  mutate(
    sign = if_else(rho > 0, "pos", "neg")
  ) %>%
  group_by(module, gene) %>%
  summarise(
    n_sig_erv = n_distinct(erv),
    n_pos = n_distinct(erv[sign == "pos"]),
    n_neg = n_distinct(erv[sign == "neg"]),
    .groups = "drop"
  ) %>%
  mutate(
    bin = case_when(
      n_sig_erv == 1 ~ "1 ERV",
      n_sig_erv <= 5 ~ "2-5 ERVs",
      TRUE ~ ">5 ERVs"
    ),
    bin = factor(
      bin,
      levels = c("1 ERV", "2-5 ERVs", ">5 ERVs")
    ),
    module = factor(
      module,
      levels = module_order_plot
    )
  )

write.csv(
  gene_sig_counts,
  file.path(reactome_dir, "gene_sig_counts.csv"),
  row.names = FALSE
)

# Calculate the proportion of genes in each bin within each module
bin_counts_prop <- gene_sig_counts %>%
  group_by(module, bin) %>%
  summarise(
    n_genes = n_distinct(gene),
    total_pos = sum(n_pos, na.rm = TRUE),
    total_neg = sum(n_neg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    module,
    bin,
    fill = list(
      n_genes = 0,
      total_pos = 0,
      total_neg = 0
    )
  ) %>%
  group_by(module) %>%
  mutate(
    prop = n_genes / sum(n_genes)
  ) %>%
  ungroup() %>%
  mutate(
    neg_frac = if_else(
      total_pos + total_neg == 0,
      0,
      total_neg / (total_pos + total_neg)
    ),
    prop_neg = prop * neg_frac
  )

# Plot
p <- ggplot(
  bin_counts_prop,
  aes(
    x = bin,
    y = prop,
    fill = module
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.75
  ) +
  geom_col_pattern(
    aes(
      y = prop_neg,
      group = module
    ),
    position = position_dodge(width = 0.8),
    width = 0.75,
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.35,
    pattern_spacing = 0.01,
    pattern_colour = "black",
    fill = "black",
    alpha = 0.25
  ) +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Proportion of genes in module",
    fill = "Module"
  )

p

ggsave(
  file.path(
    reactome_dir,
    "Reactome_gene_sigERV_bins_percentage_by_module_with_sign.pdf"
  ),
  plot = p,
  width = 9,
  height = 5
)