DESeqResults_erv <- read.table("../../Data/GSE171110/DESeqResults_ERVs_filtered.txt", header=TRUE)
DESeqResults_gene <- read.table("../../Data/GSE171110/DESeqResults_cellular_genes_filtered.txt", header=TRUE)

# Isolate up/down-regulated ERVs/Genes
source("../wrapper_functions/IsolateGenes.R", echo=TRUE, max.deparse.length=10e3)
output_dir <- "../../Data/GSE171110"
run_IsolateGenes(res = DESeqResults_erv, 
                 l2FC_limit = 1, 
                 padj_limit = 0.05, 
                 output_dir = output_dir,
                 type = "erv")
# Isolate Regulated Cellular Genes:
source("../wrapper_functions/IsolateGenes.R", echo=TRUE, max.deparse.length=10e3)
output_dir <- "../../Data/GSE171110"
run_IsolateGenes(res = DESeqResults_gene, 
                 l2FC_limit = 1, 
                 padj_limit = 0.05, 
                 output_dir = output_dir,
                 type = "genes")
