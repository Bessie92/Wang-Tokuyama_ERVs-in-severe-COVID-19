library("org.Hs.eg.db")
library("dplyr")

gene_module_df <- readRDS("./WGCNA/gene_module_df.rds")
erv_MM_sig <- readRDS("./WGCNA/erv_MM_sig.rds")

# Pull genes from ERV-correlating modules
modules_with_sig_cor_ERVs <- colnames(erv_MM_sig)
moduleLabels <- net$colors 

get_genes_in_ME <- function(ME_name) {
  module_number <- as.numeric(sub("^ME", "", ME_name))   # "ME4" -> 4
  names(moduleLabels)[moduleLabels == module_number]
}

# Build list: for each ME, which genes are inside?
module_gene_list <- list()
for (ME in modules_with_sig_cor_ERVs) {
  module_gene_list[[ME]] <- get_genes_in_ME(ME)
}

# Long data.frame with all genes in ERV-associated modules
module_gene_df <- bind_rows(
  lapply(names(module_gene_list), function(ME) {
    data.frame(
      Module = ME,
      Gene   = module_gene_list[[ME]],
      stringsAsFactors = FALSE
    )
  })
)

# Map to gene symbols
map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(module_gene_df$Gene),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

# Merge mapping back
module_gene_df2 <- module_gene_df %>%
  left_join(map_df, by = c("Gene" = "ENSEMBL")) %>% mutate(Symbol = ifelse(is.na(SYMBOL) | SYMBOL == "", Gene, SYMBOL))
module_gene_df2 <- module_gene_df2[ , -3]

module_gene_wide <- module_gene_df2 %>%
  group_by(Module) %>%
  summarize(
    N_genes = n_distinct(Gene),
    Symbols = paste(Symbol, collapse = ","),
    Ensembl = paste(Gene, collapse = ","),
    .groups = "drop"
  )

# Write one CSV listing all genes in all ERV-associated modules
write.csv(
  module_gene_wide,
  file = "./WGCNA/softpower3.csv",
  row.names = FALSE
)

# Map all unique symbols to Entrez
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = unique(module_gene_df2$Symbol),
  columns  = "ENTREZID",
  keytype  = "SYMBOL"
)

# Clean mapping (drop NAs)
symbol_map <- symbol_map %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# Join Entrez back to module_gene_df2
module_gene_annot <- module_gene_df2 %>%
  left_join(symbol_map, by = c("Symbol" = "SYMBOL"))

modules_to_annotate <- sort(unique(module_gene_annot$Module))

# Reactome analysis on ERV-correlating MEs
dir.create("./WGCNA/Reactome", showWarnings = FALSE)

library(ReactomePA)
for (ME in modules_to_annotate) {
  df_mod <- module_gene_annot %>%
    filter(Module == ME, !is.na(ENTREZID)) %>%
    distinct(ENTREZID)
  
  entrez_ids <- df_mod$ENTREZID
  
  ereact <- enrichPathway(
    gene          = entrez_ids,
    organism      = "human",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE   # converts Entrez back to symbols in result
  )
  
  ereact_df <- as.data.frame(ereact)
  
  if (nrow(ereact_df) == 0) {
    message("Reactome: ", ME, " has no significant pathways.")
    next
  }
  
  out_file <- file.path("./WGCNA/Reactome",
                        paste0("Reactome_", ME, ".csv"))
  write.csv(ereact_df, out_file, row.names = FALSE)
}
