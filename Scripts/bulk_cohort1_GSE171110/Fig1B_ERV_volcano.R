sig_ervs <- read.csv("../../Data/GSE171110/define_signature_ervs/Raw_Table1_Signature_erv_overlap_summary.csv")$erv_name
DESeqResults_erv <- read.table("../../Data/GSE171110/DESeqResults_ERVs_filtered.txt", header=TRUE)
DESeqResults_gene <- read.table("../../Data/GSE171110/DESeqResults_cellular_genes_filtered.txt", header=TRUE)

# ERV Volcano Plot
output_name <- "./figures/VolcanoPlot_ERV.pdf"
# Highlight ERV signatures on the volcano plot
select_labels <- sig_ervs
source("../wrapper_functions/VolcanoPlot.R", echo=TRUE, max.deparse.length=10e3)
run_VolcanoPlot(res = DESeqResults_erv, 
                select_labels = select_labels,
                output_name = output_name) 
