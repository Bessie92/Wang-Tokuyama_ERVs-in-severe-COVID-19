library(pheatmap)
norm_erv_count <- read.csv("../../Data/GSE171110/normalized_ervs.csv", row.names=1)
sig_ervs <- read.csv("../../Data/GSE171110/define_signature_ervs/Raw_Table1_Signature_erv_overlap_summary.csv")$erv_name
coldata <- read.csv("../../Data/GSE171110/metadata.csv", row.names=1)

sig_norm_erv_count <- norm_erv_count[rownames(norm_erv_count) %in% sig_ervs, ]

# Convert to matrix 
mat <- as.matrix(sig_norm_erv_count)
# log2-transform counts
mat_log2 <- log2(mat + 1)

# Create annotation for the heatmap
cond_annotation <- data.frame(condition = coldata$condition)
rownames(cond_annotation) <- rownames(coldata)
# Assign color for conditions
cond_colors <- list(
  condition = c("Healthy" = "#26547d", "COVID-19" = "#ec476f"))

# heatmap on log2 counts + row z-scoring
heatmap_overview <- pheatmap(
  mat_log2,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = cond_annotation,
  annotation_colors = cond_colors,
  show_rownames = TRUE, 
  show_colnames = FALSE,
  border_color = NA,
  main = "Normalized ERV heatmap (log2 counts, row-scaled)",
  filename = "./figures/Normalized_ERV_heatmap.pdf", 
  width = 7, 
  height = 5
)
