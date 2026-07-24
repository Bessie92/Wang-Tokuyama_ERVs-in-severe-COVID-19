# Raw counts 
gene_count <- read.csv("../../Data/GSE171110/merged_genes.csv", row.names=1)
erv_count <- read.csv("../../Data/GSE171110/merged_ervs.csv", row.names=1)
coldata <- read.csv("../../Data/GSE171110/metadata.csv", row.names=1)

# DEseq for ERV
output_name <- "../../Data/GSE171110/DESeqResults_ERVs.txt"
source("../wrapper_functions/DESeq2Analysis.r", echo=TRUE, max.deparse.length=10e3)
run_DESeq2(countData = erv_count, 
           coldata = coldata, 
           keep_limit = 0, 
           output_name = output_name,
           conditions = c("condition", "COVID-19", "Healthy"))
# DEseq for cellular genes
output_name <- "../../Data/GSE171110/DESeqResults_cellular_genes.txt"
source("../wrapper_functions/DESeq2Analysis.r", echo=TRUE, max.deparse.length=10e3)
run_DESeq2(countData = gene_count, 
           coldata = coldata, 
           keep_limit = 0, 
           output_name = output_name,
           conditions = c("condition", "COVID-19", "Healthy"))

DESeqResults_erv <- read.table("../../Data/GSE171110/DESeqResults_ERVs.txt", header=TRUE)
DESeqResults_gene <- read.table("../../Data/GSE171110/DESeqResults_cellular_genes.txt", header=TRUE)

# QC: remove DEseq results with base mean <10
DESeqResults_erv <- DESeqResults_erv[DESeqResults_erv$baseMean >= 10,]
DESeqResults_gene <- DESeqResults_gene[DESeqResults_gene$baseMean >= 10,]

write.table(DESeqResults_erv, file="../../Data/GSE171110/DESeqResults_ERVs_filtered.txt")
write.table(DESeqResults_gene, file="../../Data/GSE171110/DESeqResults_cellular_genes_filtered.txt")
