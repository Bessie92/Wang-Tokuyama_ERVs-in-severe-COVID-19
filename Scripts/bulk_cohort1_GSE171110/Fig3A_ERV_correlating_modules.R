coldata <- read.csv("../../Data/GSE171110/metadata.csv", row.names=1)

MEs <- readRDS("./WGCNA/MEs.rds")
gene_module_df <- readRDS("./WGCNA/gene_module_df.rds")

# Identify COVID-19 associated modules 

## double check samples are still in order
all.equal(rownames(coldata), rownames(MEs))
## Create the design matrix from the `condition` variable
coldata$condition <- relevel(factor(coldata$condition), ref = "Healthy")
des_mat <- model.matrix(~ coldata$condition)

## Run linear model on each module
fit <- limma::lmFit(t(MEs), design = des_mat)
fit <- limma::eBayes(fit)

stats_df <- limma::topTable(fit, number = ncol(MEs)) %>%
  tibble::rownames_to_column("module")

module_covid <- stats_df %>% filter(adj.P.Val < 0.05)
write.csv(module_covid, "./WGCNA/COVID-19 associated modules.csv", row.names = FALSE)

# Figure 3A left panel: Visualize LogFC of all modules
df_plot <- stats_df %>%
  mutate(
    direction = case_when(
      adj.P.Val >= 0.05 ~ "Not significant",
      logFC > 0         ~ "Upregulated",
      logFC < 0         ~ "Downregulated"
    )
  )

df_plot$module <- factor(
  df_plot$module,
  levels = df_plot$module[order(df_plot$logFC)]
)

ggplot(df_plot, aes(x = module, y = logFC, fill = direction)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "Not significant" = "grey70",
      "Downregulated"   = "#377eb8",  # blue
      "Upregulated"     = "#e41a1c"   # red
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  labs(
    x = "Module",
    y = "log2 Fold Change",
    fill = "Module status",
    title = "Module differential expression (limma)"
  ) +
  theme_classic(base_size = 12)
ggsave("./WGCNA/Module_logFC_bar_plot.pdf", width = 10, height = 7)

# Correlate signature ERVs with COVID-associated MEs 
MEs_COVID <- MEs[ , colnames(MEs) %in% module_covid$module]
expr_sig_erv <- norm_erv_count[rownames(norm_erv_count) %in% sig_ervs, ]
## Match sample order to MEs (rows of MEs = samples)
all(colnames(expr_sig_erv) == rownames(MEs_COVID))
nSamples <- nrow(MEs_COVID)
# cor expects samples x variables, so transpose ERV matrix
erv_MM <- cor(t(expr_sig_erv), MEs_COVID, method = "spearman")
erv_MM_p <- corPvalueStudent(erv_MM, nSamples)
erv_MM_p_adj <- matrix(
  p.adjust(as.vector(erv_MM_p), method = "BH"),
  nrow = nrow(erv_MM_p),
  dimnames = dimnames(erv_MM_p)
)
hits <- which(erv_MM_p_adj < 0.05, arr.ind = TRUE)
erv_module_pairs <- data.frame(
  ERV    = rownames(erv_MM)[hits[, 1]],
  Module = colnames(erv_MM)[hits[, 2]],
  rho    = erv_MM[hits],
  padj   = erv_MM_p_adj[hits],
  row.names = NULL
)
erv_MM_sig <- erv_MM[ , colnames(erv_MM) %in% erv_module_pairs$Module]

# Figure 3A right panel: significant ERV-module correlations on a heatmap
library(pheatmap)
pheatmap(
  erv_MM_sig,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  border_color = NA,
  filename = "./WGCNA/ERV-module-correlation.pdf"
)

saveRDS(stats_df, "./WGCNA/module_limma_stats.rds")
saveRDS(module_covid, "./WGCNA/module_covid.rds")
saveRDS(erv_module_pairs, "./WGCNA/erv_module_pairs.rds")
saveRDS(erv_MM_sig, "./WGCNA/erv_MM_sig.rds")
