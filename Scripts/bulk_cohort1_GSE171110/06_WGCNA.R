library(WGCNA)
library(flashClust)
library(curl)

norm_gene_count <- read.csv("../../Data/GSE171110/normalized_genes.csv", row.names=1)
norm_all_count <- norm_gene_count
norm_all_count_t <- as.data.frame(t(norm_all_count))

dir.create("WGCNA")

# Identify outlier genes
gsg <- goodSamplesGenes(norm_all_count_t)
summary(gsg)
gsg$allOK

## If gsg$allOK = FALSE, filter outlier genes out with:
if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0) 
printFlush(paste("Removing genes:", paste(names(norm_all_count_t)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(norm_all_count_t)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
norm_all_count_t <- norm_all_count_t[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

# Identify outlier samples
sampleTree <- hclust(dist(norm_all_count_t), method = "average") #Clustering samples based on distance 

## Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

## Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

spt <- pickSoftThreshold(norm_all_count_t) 

# Determining the Soft Power Threshold
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

softPower = 3

# Calling the Adjacency Function
datExpr <- norm_all_count_t
print(dim(datExpr))    # nSamples x nGenes

## Main WGCNA network construction
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  numericLabels = TRUE,
  saveTOMs = FALSE,         # don’t save TOM, too big
  maxBlockSize = 8000,
  randomSeed = 1234
)

moduleColors <- labels2colors(net$colors)
## Extract module eigengenes
MEs <- net$MEs
saveRDS(MEs, "./WGCNA/MEs.rds")

## Save a gene-to-module table from net
gene_module_df <- data.frame(
  Gene = names(net$colors),
  module_number = net$colors,
  module = paste0("ME", net$colors),
  module_color = labels2colors(net$colors),
  stringsAsFactors = FALSE
)
saveRDS(gene_module_df, "./WGCNA/gene_module_df.rds")

# Hierarchical clustering analysis
## Creating the dendogram 
geneTree <- net$dendrograms[[1]]
## Plot the dendogram 
plotDendroAndColors(
  geneTree,
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03
)

####### Next script in chronological order: Fig3A_ERV_correlating_modules.R #######
