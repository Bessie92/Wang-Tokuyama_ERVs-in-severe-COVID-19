# Runs DESeq on input data
run_DESeq2 <- function(countData, coldata, keep_limit, output_name, conditions) {
  
  # Libraries
  library(ggplot2)
  library(DESeq2)
  
  cts <- as.matrix(countData)

  # Check if data is columns/rows of data is aligned
  if (!all(rownames(coldata) == colnames(cts))) {
    print("Data un-aligned")
    break
  }
  
  # DESeqDataSetFromMatrix input
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~condition)

  # Keeps significant rows
  keep <- rowSums(counts(dds)) > keep_limit
  dds <- dds[keep,]

  # Differential Expression
  dds <- DESeq(dds) 
  res <- results(dds, contrast=conditions) # Change for disease
  res <- na.omit(res)
  write.table(res, file=output_name)
  res
}
