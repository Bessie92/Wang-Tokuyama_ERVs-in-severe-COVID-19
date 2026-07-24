library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)

# Make sure 07_Reactome_enrichment.R has been run!

# Inputs
reactome_dir <- "./WGCNA/Reactome"
outdir <- reactome_dir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

modules_of_interest <- c("ME0", "ME9", "ME61", "ME2", "ME48", "ME5")
module_order_plot <- c("ME5", "ME48", "ME2", "ME61", "ME9", "ME0")

coldata <- read.csv("../../Data/GSE171110/metadata.csv", row.names = 1)
norm_gene_count_name <- read.csv("../../Data/GSE171110/normalized_genes_with_name.csv", row.names = 1)
norm_erv_count <- read.csv("../../Data/GSE171110/normalized_ervs.csv", row.names = 1)
sig_ervs <- read.csv("../../Data/GSE171110/define_signature_ervs/Raw_Table1_Signature_erv_overlap_summary.csv")$erv_name
sig_up_ervs <- sig_ervs[!sig_ervs %in% c("6217", "5208", "4352")]

# Helper functions
read_reactome_module <- function(module, dir = reactome_dir) {
  f <- file.path(dir, paste0("Reactome_", module, ".csv"))
  stopifnot(file.exists(f))
  read.csv(f, stringsAsFactors = FALSE) %>%
    mutate(module = module)
}

subset_samples <- function(mat, samples) {
  mat[, colnames(mat) %in% samples, drop = FALSE]
}

subset_features <- function(mat, features) {
  mat[rownames(mat) %in% features, , drop = FALSE]
}

# === Figure 3B: Top 5 Reactome terms per selected module ===
reactome_all <- purrr::map_dfr(modules_of_interest, read_reactome_module)

reactome_top5 <- reactome_all %>%
  dplyr::filter(!is.na(p.adjust)) %>%
  dplyr::group_by(module) %>%
  dplyr::slice_min(
    order_by = p.adjust,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

write.csv(
  reactome_top5,
  file.path(outdir, "reactome_top_5_per_ME.csv"),
  row.names = FALSE
)

reactome_top5_plot <- reactome_top5 %>%
  group_by(module) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  mutate(
    Description = factor(Description, levels = rev(unique(Description)))
  ) %>%
  ungroup()

reactome_top5_plot$module <- factor(
  reactome_top5_plot$module,
  levels = c("ME5", "ME48", "ME2", "ME61", "ME9", "ME0")
)

dotplot_out <- file.path(outdir, "Reactome_top5_terms_dotplot_selected_modules.pdf")

p_dot <- ggplot(
  reactome_top5_plot,
  aes(
    x = module,
    y = Description,
    size = Count,
    color = p.adjust
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    name = "Gene count",
    range = c(2, 8)
  ) +
  scale_color_viridis_c(
    option = "D",
    direction = -1,
    name = "Adjusted p-value"
  ) +
  labs(
    x = "Module",
    y = "Reactome pathway (top 5 by adj. p)",
    title = "Top 5 Reactome pathways per selected module"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )

p_dot
ggsave(dotplot_out, p_dot, width = 15, height = 10)

# Figure 3C: Spearman correlation between reactome term genes and upreg signatrue ERVs === 
## Prepare merged Reactome terms
padj_col <- "p.adjust"

terms0 <- reactome_top5 %>%
  filter(module %in% modules_of_interest) %>%
  transmute(
    module,
    ID,
    Description,
    term_padj = .data[[padj_col]],
    gene_set = str_split(geneID, "/", simplify = FALSE)
  ) %>%
  mutate(
    gene_set = map(gene_set, ~ sort(unique(.x[.x != "" & !is.na(.x)])))
  )

## Each (module, gene) kept only in the term with the smallest adjusted p-value
gene2term <- terms0 %>%
  dplyr::select(module, ID, Description, term_padj, gene_set) %>%
  unnest(gene_set, names_repair = "minimal") %>%
  rename(gene = gene_set) %>%
  group_by(module, gene) %>%
  slice_min(order_by = term_padj, with_ties = FALSE) %>%
  ungroup()

terms_pruned <- gene2term %>%
  group_by(module, ID, Description, term_padj) %>%
  summarise(
    gene_set = list(sort(unique(gene))),
    n_genes = length(gene_set[[1]]),
    .groups = "drop"
  ) %>%
  filter(n_genes > 0)

## Merge identical gene sets
terms_unique <- terms_pruned %>%
  mutate(
    gene_set_str = map_chr(gene_set, ~ paste(.x, collapse = ";"))
  ) %>%
  group_by(module, gene_set_str) %>%
  summarise(
    IDs = paste(unique(ID), collapse = ","),
    Description = paste(unique(Description), collapse = " / "),
    term_padj = min(term_padj, na.rm = TRUE),
    gene_set = list(first(gene_set)),
    n_genes = length(first(gene_set)),
    .groups = "drop"
  ) %>%
  group_by(module) %>%
  mutate(term_id = row_number()) %>%
  ungroup()

## If a larger term fully contains a smaller term within the same module,
### keep the larger term and append the smaller term name to it
merge_terms_within_module <- function(dfm) {
  n_terms <- nrow(dfm)
  if (n_terms == 0) return(dfm[0, , drop = FALSE])

  final_parent <- seq_len(n_terms)
  final_name_vec <- dfm$Description

  for (i in seq_len(n_terms)) {
    genes_i <- dfm$gene_set[[i]]
    for (j in seq_len(n_terms)) {
      if (i == j) next
      genes_j <- dfm$gene_set[[j]]

      if (length(genes_i) < length(genes_j) && all(genes_i %in% genes_j)) {
        final_parent[i] <- final_parent[j]
        final_name_vec[j] <- paste0(
          final_name_vec[j],
          " (includes ",
          dfm$Description[i],
          ")"
        )
      }
    }
  }

  dfm %>%
    mutate(
      final_parent = final_parent,
      final_name = final_name_vec
    ) %>%
    group_by(module, final_parent) %>%
    summarise(
      term_name = first(final_name),
      IDs = paste(unique(IDs), collapse = ","),
      term_padj = min(term_padj, na.rm = TRUE),
      gene_set = list(first(gene_set)),
      n_genes = length(first(gene_set)),
      .groups = "drop"
    )
}

merged_terms <- terms_unique %>%
  group_split(module) %>%
  map_dfr(merge_terms_within_module) %>%
  mutate(
    term_uid = paste(module, term_name, sep = " | ")
  )

merged_terms_out <- merged_terms %>%
  mutate(
    genes = map_chr(gene_set, ~ paste(.x, collapse = "/"))
  ) %>%
  dplyr::select(module, term_name, IDs, term_padj, n_genes, genes)

write.csv(
  merged_terms_out,
  file.path(outdir, "Reactome_merged_terms_summary_pruned.csv"),
  row.names = FALSE
)

## Calculate spearman correlations
covid_samples <- rownames(coldata)[coldata$condition == "COVID-19"]

norm_gene_count_name_covid <- subset_samples(norm_gene_count_name, covid_samples)
norm_erv_count_covid <- subset_samples(norm_erv_count, covid_samples)
norm_erv_count_up_sig_covid <- subset_features(norm_erv_count_covid, sig_up_ervs)

common_samples <- intersect(
  colnames(norm_gene_count_name_covid),
  colnames(norm_erv_count_up_sig_covid)
)

gene_mat <- as.matrix(norm_gene_count_name_covid[, common_samples, drop = FALSE])
erv_mat  <- as.matrix(norm_erv_count_up_sig_covid[, common_samples, drop = FALSE])

cor_list <- vector("list", length = 0)

for (t in seq_len(nrow(merged_terms))) {
  term_uid  <- merged_terms$term_uid[t]
  term_name <- merged_terms$term_name[t]
  mod       <- merged_terms$module[t]
  genes     <- merged_terms$gene_set[[t]]

  genes_in_mat <- intersect(genes, rownames(gene_mat))
  if (length(genes_in_mat) == 0) next

  for (g in genes_in_mat) {
    x_all <- as.numeric(gene_mat[g, ])

    for (erv in rownames(erv_mat)) {
      y_all <- as.numeric(erv_mat[erv, ])
      ok <- complete.cases(x_all, y_all)
      if (sum(ok) < 3) next

      ct <- suppressWarnings(
        cor.test(x_all[ok], y_all[ok], method = "spearman", exact = FALSE)
      )

      cor_list[[length(cor_list) + 1]] <- data.frame(
        module = mod,
        term_uid = term_uid,
        term_name = term_name,
        gene = g,
        erv = erv,
        rho = unname(ct$estimate),
        pvalue = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

cor_df <- bind_rows(cor_list) %>%
  group_by(term_uid, gene) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  ungroup()

write.csv(
  cor_df,
  file.path(outdir, "Reactome_gene_ERV_spearman_correlations_mergedTerms_pruned.csv"),
  row.names = FALSE
)

## Heatmap visualization: top 20 genes per merged term
gene_ranking_term <- cor_df %>%
  group_by(module, term_uid, term_name, gene) %>%
  summarise(
    n_sig_erv = sum(padj < 0.05, na.rm = TRUE),
    max_abs_rho = max(abs(rho), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_sig_erv > 0)

top20_genes_term <- gene_ranking_term %>%
  mutate(module = factor(module, levels = module_order_plot)) %>%
  group_by(module, term_uid) %>%
  arrange(desc(n_sig_erv), desc(max_abs_rho), gene, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  top20_genes_term,
  file.path(outdir, "Reactome_top20_genes_per_merged_term.csv"),
  row.names = FALSE
)

cor_top <- cor_df %>%
  semi_join(top20_genes_term, by = c("module", "term_uid", "gene")) %>%
  filter(padj < 0.05) %>%
  mutate(module = factor(module, levels = module_order_plot))

y_levels <- top20_genes_term %>%
  mutate(module = factor(module, levels = module_order_plot)) %>%
  arrange(module, term_uid, desc(n_sig_erv), desc(max_abs_rho), gene) %>%
  mutate(term_gene = paste0(term_name, ": ", gene)) %>%
  pull(term_gene)

cor_top <- cor_top %>%
  mutate(
    term_gene = paste0(term_name, ": ", gene),
    term_gene = factor(term_gene, levels = rev(y_levels))
  )

erv_order <- cor_top %>%
  group_by(erv) %>%
  summarise(n_pts = n(), .groups = "drop") %>%
  arrange(desc(n_pts)) %>%
  pull(erv)

cor_top <- cor_top %>%
  mutate(erv = factor(erv, levels = erv_order))

heatmap_out <- file.path(
  outdir,
  "Reactome_gene_ERV_heatmap_top20_per_mergedTerm_byModule.pdf"
)

p_heat <- ggplot(cor_top, aes(x = erv, y = term_gene, fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Spearman rho"
  ) +
  labs(
    x = "Signature ERVs",
    y = "Merged Reactome gene (top 20 per term)",
    title = "Gene–ERV correlations across merged Reactome terms"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  )

p_heat
ggsave(heatmap_out, p_heat, width = 10, height = 17)

