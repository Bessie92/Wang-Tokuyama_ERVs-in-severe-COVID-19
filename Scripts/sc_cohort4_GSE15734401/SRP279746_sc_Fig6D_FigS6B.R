library(edgeR)
library(MAST)

## Shown below is the code for one cell type
### Define fit_model
fit_model <- function(adata_){
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group = colData(adata_)$condition)
    # filter out genes with low counts
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    #keep <- filterByExpr(y) # Default EdgeR filtering,keeps genes with at least ~10 counts in a few   (e.g., ≥2–3) samples from the smallest group
    keep <- filterByExpr(y, min.count = 5) # try a less stringent filtering
    y <- y[keep, , keep.lib.sizes=FALSE]
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    # normalize
    y <- calcNormFactors(y)
    # create a vector that is a concatenation of condition and cell type that we will later use with contrasts
    group <- paste0(colData(adata_)$condition, ".", colData(adata_)$cell_type)
    replicate <- colData(adata_)$replicate
    # create a design matrix: here we have multiple donors so also consider that in the design matrix
    design <- model.matrix(~ 0 + group + replicate)
    # estimate dispersion
    y <- estimateDisp(y, design = design)
    # fit the model
    fit <- glmQLFit(y, design)
    return(list("fit"=fit, "design"=design, "y"=y))
}

library(zellkonverter)
sce <- readH5AD(
  "./Pseudobulk_H5ADs/pseudobulk_Erythroid.h5ad",
  reader   = "R",    # <- important
  use_hdf5 = FALSE   # or TRUE if the file is huge
)

outs <- fit_model(sce)
fit  <- outs$fit
y    <- outs$y
design <- outs$design

plotMDS(y, col=ifelse(y$samples$group == "severe COVID", "red", "blue"))
plotBCV(y)

colnames(y$design)
# Make them compatible with R
colnames(y$design) <- make.names(colnames(y$design))
colnames(y$design)

myContrast <- makeContrasts('groupsevere.COVID.Erythroid.like.and.erythroid.precursor.cells-grouphealthy.Erythroid.like.and.erythroid.precursor.cells', levels = y$design)
qlf <- glmQLFTest(fit, contrast=myContrast)
# get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
tt <- topTags(qlf, n = Inf)
tt <- tt$table
tr <- glmTreat(fit, contrast=myContrast, lfc=1.5)
print(head(topTags(tr)))

# ============== Fig. S6B ==============
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(tidyr)

# Extract significant up and downregulated genes
sig_up_genes <- tt %>%
  filter(logFC > 1, FDR < 0.05) %>%
  rownames()
sig_down_genes <- tt %>%
  filter(logFC < -1, FDR < 0.05) %>%
  rownames()

# Convert gene symbols to Entrez IDs
gene_entrez_up <- bitr(sig_up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_entrez_down <- bitr(sig_down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up <- enrichGO(gene         = gene_entrez_up$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable     = TRUE)
ego_down <- enrichGO(gene         = gene_entrez_down$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable     = TRUE)
write.csv(ego_up@result, file = "./GO_cellular/GO_enrichment_results_upreg_Ery.csv", row.names = FALSE)
write.csv(ego_down@result, file = "./GO_cellular/GO_enrichment_results_downreg_Ery.csv", row.names = FALSE)

# Dotplot - top 5 enriched pathways
dotplot(ego_up, showCategory = 5, title = "Upreg GO Biological Processes")
ggsave("./GO_cellular/dotplot_upreg_Ery.pdf", width = 7, height=5)
dotplot(ego_down, showCategory = 5, title = "Downreg GO Biological Processes")
ggsave("./GO_cellular/dotplot_downreg_Ery.pdf", width = 7, height=5)

# ============== Fig. 6D ==============
library(tidyverse)
indir <- "./DE_ERV"
files <- list.files(indir, pattern="\\.csv$", full.names=TRUE)

# read one file + tag cell type
read_edger <- function(f) {
  cell_type <- tools::file_path_sans_ext(basename(f))
  df <- readr::read_csv(f, show_col_types = FALSE)
  if (!"...1" %in% names(df)) stop("Expected first column to be ...1 in: ", f)
  df %>%
    rename(erv_raw = `...1`) %>%
    mutate(
      cell_type = cell_type,
      erv_raw = as.character(erv_raw),
      erv = sub("\\.\\d+$", "", erv_raw),       # 4720.1 -> 4720
      had_suffix = grepl("\\.\\d+$", erv_raw)
    ) %>%
    arrange(had_suffix) %>%
    group_by(cell_type, erv) %>%
    slice(1) %>%
    ungroup()
}

all_res <- purrr::map_dfr(files, read_edger)

plot_df <- all_res %>%
  mutate(sig = FDR < 0.05)

# Group sig ERVs together
erv_group <- c("1874", "2124", "2637", "3052", "4830", "4896", "W-8")

plot_df2 <- plot_df %>%
  mutate(
    erv = as.character(erv),
    erv = factor(
      erv,
      levels = c(erv_group, setdiff(sort(unique(erv)), erv_group))
    )
  )

ggplot(plot_df2, aes(x = cell_type, y = erv, fill = logFC)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "logFC") +
  geom_text(data = subset(plot_df2, sig), aes(label="*"), size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("Tile_plot_DE_ERV_across_cell_types.pdf", width = 7, height = 15)
