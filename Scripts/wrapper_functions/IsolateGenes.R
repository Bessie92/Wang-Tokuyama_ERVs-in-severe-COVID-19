library(tidyverse)

# Creates text files for the upregulated and downregulated ervs/genes
run_IsolateGenes <- function(res, l2FC_limit, padj_limit, output_dir, type) {

  # Up-regulated Genes
  upreg <- filter(res, log2FoldChange > l2FC_limit, padj < padj_limit)
  upreg <- upreg[order(upreg$padj),]
  upreg_names <- as.vector(rownames(upreg))
  cat(paste(shQuote(upreg_names, type="cmd"), collapse=","), file=paste(output_dir, "/upreg_", type, ".txt", sep=""))

  # Down-regulated Genes
  downreg <- filter(res, log2FoldChange < -l2FC_limit, padj < padj_limit)
  downreg <- downreg[order(downreg$padj),]
  downreg_names <- as.vector(rownames(downreg))
  cat(paste(shQuote(downreg_names, type="cmd"), collapse=","), file=paste(output_dir, "/downreg_", type, ".txt", sep=""))
}
