library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(patchwork)

# Read in ERVmap2 output files
gene_count <- read.csv("../../Data/GSE163317/merged_genes.csv", row.names=1)
erv_count <- read.csv("../../Data/GSE163317/merged_ervs.csv", row.names=1)

# Normalize raw cellular genes using DEseq2 and generate size factors per sample
library(DESeq2)
conds <- factor( c(1:ncol(gene_count)))
dds <- DESeqDataSetFromMatrix(countData = as.matrix(gene_count), colData = as.matrix(conds), design = ~ 1)
dds <- estimateSizeFactors(dds)
normalizedCounts <- counts(dds, normalized=T)
write.table(normalizedCounts,
            file = "../../Data/GSE163317/normalized_genes.csv",
            sep = ",",
            quote = FALSE,
            col.names = NA)
write.table(sizeFactors(dds), "../../Data/GSE163317/normalized_factors.txt", sep="\t", quote=FALSE, col.names=NA)

## Normalize raw ERV counts using cellular gene size factors
size_factor <- read.table("../../Data/GSE163317/normalized_factors.txt")
colnames(size_factor) <- c("sizeFactor")

norm_erv_count <- erv_count
# Loop through each column in erv_count
for (sampleName in colnames(erv_count)) {
  # Check if the sample name exists in the row names of size_factor
  if (sampleName %in% rownames(size_factor)) {
    # Retrieve the corresponding size factor
    sizeFactor <- size_factor[sampleName, ]
    
    # Apply the size factor to normalize ERV counts
    norm_erv_count[, sampleName] <- erv_count[, sampleName] / sizeFactor
  } else {
    warning(paste("Size factor for", sampleName, "not found!"))
  }
}
# Output the normalized ERV counts as a .csv file
write.table(norm_erv_count, "../../Data/GSE163317/normalized_ervs.csv", sep = ",", quote=F, col.names = NA)

# Differential expression analysis - paired
coldata <- read.csv("../../Data/GSE163317/coldata.csv", row.names=1)

erv_count <- erv_count[, rownames(coldata)]
gene_count <- gene_count[, rownames(coldata)]

coldata$Patient <- factor(coldata$Patient)
coldata$Day <- factor(coldata$Day, levels = c("0", "7"))

source("../wrapper_functions/run_DESeq2_paired.r", echo=TRUE, max.deparse.length=10e3)
DESeqResults_erv <- run_DESeq2_paired(
  countData = erv_count,
  coldata = coldata,
  paired_var = "Patient",
  condition_var = "Day",
  reference_level = "0",
  test_level = "7",
  output_prefix = "ERV_paired"
)

DESeqResults_gene <- run_DESeq2_paired(
  countData = gene_count,
  coldata = coldata,
  paired_var = "Patient",
  condition_var = "Day",
  reference_level = "0",
  test_level = "7",
  output_prefix = "Gene_paired"
)

# Isolate sig DE erv/genes
DESeqResults_erv <- DESeqResults_erv$filtered
DESeqResults_gene <- DESeqResults_gene$filtered

source("../wrapper_functions/IsolateGenes.R", echo=TRUE, max.deparse.length=10e3)
output_dir <- "../../Data/GSE163317/"

run_IsolateGenes(res = DESeqResults_erv, 
                 l2FC_limit = 1, 
                 padj_limit = 0.05, 
                 output_dir = output_dir,
                 type = "erv")

run_IsolateGenes(res = DESeqResults_gene, 
                 l2FC_limit = 1, 
                 padj_limit = 0.05, 
                 output_dir = output_dir,
                 type = "genes")

# Figure 4C: Volcano plot - DE ERVs
dir.create("figures")
output_name <- "./figures/VolcanoPlot_ERV.pdf"
# Highlight upregulated severe COVID-19 signature ERVs on the plot
ervs_to_plot <- c("2724", "2124", "4745", "942", "4896", "4830")
source("../wrapper_functions/VolcanoPlot.R", echo=TRUE, max.deparse.length=10e3)
run_VolcanoPlot(res = DESeqResults_erv, 
                select_labels = ervs_to_plot,
                output_name = output_name) 
# Figure 4D: Heatmap, individual patients, before and after treatment 
library(pheatmap)
# ERVs of interest
# Subset ERVs
heatmap_mat <- norm_erv_count[rownames(norm_erv_count) %in% ervs_to_plot, ]

# Make sure metadata matches heatmap columns
coldata_sub <- coldata[colnames(heatmap_mat), ]

# Make Day an ordered factor
coldata_sub$Day <- factor(coldata_sub$Day, levels = c("0", "7"))

# Order samples by Patient, then Day
sample_order <- rownames(coldata_sub)[order(coldata_sub$Patient, coldata_sub$Day)]

# Reorder heatmap matrix and metadata
heatmap_mat <- heatmap_mat[, sample_order]
coldata_sub <- coldata_sub[sample_order, ]

# Z-score by ERV / row
heatmap_mat_scaled <- t(scale(t(as.matrix(heatmap_mat))))

# Column annotation
annotation_col <- data.frame(
  Patient = coldata_sub$Patient,
  Day = coldata_sub$Day
)

rownames(annotation_col) <- rownames(coldata_sub)

# Plot heatmap
pheatmap(
  heatmap_mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,   # keeps Day 0 and Day 7 side by side
  annotation_col = annotation_col,
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  main = "Selected ERV expression: Day 0 vs Day 7 paired", 
  filename = "./figures/sig_erv_heatmap.pdf"
)

