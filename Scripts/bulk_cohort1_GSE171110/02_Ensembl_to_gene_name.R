norm_gene_count <- read.csv("../../Data/GSE171110/normalized_genes.csv", row.names=1)

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- rownames(norm_gene_count)

genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
               filters = 'ensembl_gene_id', 
               values = gene_ids, 
               mart = ensembl)
name_mapping <- setNames(genes$external_gene_name, genes$ensembl_gene_id)
new_row_names <- name_mapping[rownames(norm_gene_count)]

# Replace NA or empty gene names with original identifiers
new_row_names[is.na(new_row_names) | new_row_names == ""] <- rownames(norm_gene_count)[is.na(new_row_names) | new_row_names == ""]

# Function to append a suffix to duplicate names
make_unique <- function(names) {
  counts <- table(names)
  names <- as.character(names)
  for (name in names(counts[counts > 1])) {
    inds <- which(names == name)
    names[inds] <- paste0(name, "_", seq_along(inds))
  }
  names
}

# Apply the function to make row names unique
unique_row_names <- make_unique(new_row_names)

# Replace row names
rownames(norm_gene_count) <- unique_row_names

# Save this count table with names
write.csv(norm_gene_count, "../../Data/GSE171110/normalized_genes_with_name.csv")
