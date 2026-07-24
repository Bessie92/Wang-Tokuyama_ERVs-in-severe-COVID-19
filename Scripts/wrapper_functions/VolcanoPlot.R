run_VolcanoPlot <- function(res, select_labels, output_name) {
  
  # Libraries
  library(ggplot2)
  library(dplyr)
  library(EnhancedVolcano)
  library(ggdogs)
  options(ggrepel.max.overlaps=Inf)
  
  Log2FC <-res$log2FoldChange
  Padj <-res$padj
  
  ## make sure the folder exists
  dir.create(dirname(output_name), showWarnings = FALSE, recursive = TRUE)
  
  ## open PDF device
  pdf(file = output_name, width = 8, height = 7)
  
  upreg  <- filter(res, Log2FC > 1, Padj < 0.05)
  downreg <- filter(res, Log2FC < -1, Padj < 0.05)
  
  keyvals <- ifelse(
    row.names(res) %in% select_labels, '#ec476f',
      ifelse(row.names(res) %in% row.names(upreg) |
               row.names(res) %in% row.names(downreg), 'royalblue', 'grey30'))
  names(keyvals)[keyvals == '#ec476f'] <- 'Severe COVID-19 Signatures'
  names(keyvals)[keyvals == 'royalblue'] <- 'Significant'
  names(keyvals)[keyvals == 'grey30'] <- 'Not Signficant'
  
  # Volcano Plot
  P1 <- EnhancedVolcano(res, x='log2FoldChange', y='padj', lab=row.names(res),
                    border="full",
                    selectLab = select_labels,
                    ylim=c(0,14),        
                    #xlim=c(-4,4),
                    labSize=4,
                    pCutoff=0.05,
                    #legendPosition='right',
                    #legendLabSize=12,
                    #legendIconSize=3,
                    #legendLabels=c('Not sig.','Signature ERVs', 'Log (base 2) FC', 'p-value & Log (base 2) FC'),
                    legendPosition = 'none', # REMOVE TO SHOW LEGEND
                    gridlines.minor=FALSE,
                    titleLabSize=12,
                    subtitleLabSize=2,
                    axisLabSize=12,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5,
                    typeConnectors = "open",
                    #col=c('grey30', 'grey30', 'grey30', 'royalblue'),
                    colCustom = keyvals,
                    colAlpha=0.5)
  P1+ geom_point(aes(select_labels, colour="#ec476f"))
  P1+ theme(plot.margin=unit(c(0,1,0,1), "cm"))
  # Highlight selected lables with alpha =1
  P1 <- P1 + 
    geom_point(
      data = subset(res, row.names(res) %in% select_labels),
      aes(x = log2FoldChange, y = -log10(padj)),
      color = "#ec476f",
      alpha = 1,
      size = 2.5
    )
  
  print(P1)
  
  dev.off()
}
