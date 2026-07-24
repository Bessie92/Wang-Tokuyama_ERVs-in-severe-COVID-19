# Read in ERVmap2 output files
gene_count <- read.csv("../../Data/GSE171110/merged_genes.csv", row.names=1)
erv_count <- read.csv("../../Data/GSE171110/merged_ervs.csv", row.names=1)

# Normalize raw cellular genes using DEseq2 and generate size factors per sample
library(DESeq2)
conds <- factor( c(1:ncol(gene_count)))
dds <- DESeqDataSetFromMatrix(countData = as.matrix(gene_count), colData = as.matrix(conds), design = ~ 1)
dds <- estimateSizeFactors(dds)
normalizedCounts <- counts(dds, normalized=T)
write.table(normalizedCounts,
            file = "../../Data/GSE171110/normalized_genes.csv",
            sep = ",",
            quote = FALSE,
            col.names = NA)
write.table(sizeFactors(dds), "../../Data/GSE171110/normalized_factors.txt", sep="\t", quote=FALSE, col.names=NA)

## Normalize raw ERV counts using cellular gene size factors
size_factor <- read.table("../../Data/GSE171110/normalized_factors.txt")
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
write.table(norm_erv_count, "../../Data/GSE171110/normalized_ervs.csv", sep = ",", quote=F, col.names = NA)

