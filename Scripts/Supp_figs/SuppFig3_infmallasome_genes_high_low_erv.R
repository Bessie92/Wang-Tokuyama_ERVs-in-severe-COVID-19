library(dplyr)
library(tidyr)
library(corrplot)
library(AnnotationDbi)
library(org.Hs.eg.db)

outdir <- "./COVID_subsets_based_on_ervs"
reactome_dir <- "./WGCNA/Reactome"

norm_erv_count <- read.csv("../../Data/GSE171110/normalized_ervs.csv",row.names = 1)
norm_gene_count <- read.csv("../../Data/GSE171110/normalized_genes.csv",row.names = 1,)
norm_gene_count_name <- read.csv("../../Data/GSE171110/normalized_genes_with_name.csv",row.names = 1)
coldata <- read.csv("../../Data/GSE171110/metadata.csv",row.names = 1)
sig_ervs <- read.csv(paste0("../../Data/GSE171110/define_signature_ervs/","Raw_Table1_Signature_erv_overlap_summary.csv"),stringsAsFactors = FALSE)$erv_name
reactome_top5 <- read.csv(file.path(reactome_dir, "reactome_top_5_per_ME.csv"))
de_results <- read.table(
  file.path(
    outdir,
    "DESeqResults_cellular_genes_ERV_subset_filtered.txt"
  ),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

low_erv_covid <- c(
  "SRX10470994", "SRX10470987", "SRX10470991", "SRX10470998",
  "SRX10470981", "SRX10470997", "SRX10470990", "SRX10471000"
)
covid_samples <- rownames(coldata)[
  coldata$condition == "COVID-19"
]
high_erv_covid <- setdiff(
  covid_samples,
  low_erv_covid
)
healthy_samples <- rownames(coldata)[
  coldata$condition == "Healthy"
]


# Select genes upregulated in ERV-high and present in Reactome terms

reactome_genes <- reactome_top5 %>%
  dplyr::select(geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::filter(!is.na(geneID), geneID != "") %>%
  dplyr::distinct(geneID) %>%
  dplyr::pull(geneID)

de_df <- de_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl") %>%
  dplyr::mutate(
    ensembl = sub("\\..*$", "", ensembl)
  )

de_df$symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = de_df$ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

selected_genes <- de_df %>%
  dplyr::filter(
    padj < 0.05,
    log2FoldChange > 1,
    symbol %in% reactome_genes
  ) %>%
  dplyr::pull(ensembl) %>%
  unique()

# Helper for building one combined gene–ERV matrix

make_combined_matrix <- function(gene_counts, genes, ervs, samples) {
  
  samples <- samples[
    samples %in% colnames(gene_counts) &
      samples %in% colnames(norm_erv_count)
  ]
  
  gene_mat <- as.matrix(
    gene_counts[
      rownames(gene_counts) %in% genes,
      samples,
      drop = FALSE
    ]
  )
  
  erv_mat <- as.matrix(
    norm_erv_count[
      rownames(norm_erv_count) %in% ervs,
      samples,
      drop = FALSE
    ]
  )
  
  rbind(gene_mat, erv_mat)
}

# Helper for correlation and plotting

make_corrplot <- function(
    combined_mat,
    feature_order = NULL,
    output_file,
    width = 15,
    height = 15
) {
  
  cor_mat <- cor(
    t(combined_mat),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  p_mat <- cor.mtest(
    t(combined_mat),
    method = "spearman"
  )$p
  
  if (is.null(feature_order)) {
    hc_order <- corrplot::corrMatOrder(
      cor_mat,
      order = "hclust"
    )
    
    feature_order <- rownames(cor_mat)[hc_order]
  }
  
  cor_mat <- cor_mat[
    feature_order,
    feature_order,
    drop = FALSE
  ]
  
  p_mat <- p_mat[
    feature_order,
    feature_order,
    drop = FALSE
  ]
  
  plot_cols <- colorRampPalette(
    c("#2166AC", "#67A9CF", "white", "#EF8A62", "#B2182B")
  )(20)
  
  pdf(output_file, width = width, height = height)
  
  corrplot(
    cor_mat,
    p.mat = p_mat,
    order = "original",
    method = "square",
    diag = FALSE,
    insig = "blank",
    addrect = 7,
    rect.col = "black",
    rect.lwd = 3,
    tl.pos = "full",
    col = plot_cols
  )
  
  dev.off()
  
  feature_order
}

# Panel A: ERV-high-upregulated Reactome genes

combined_high <- make_combined_matrix(
  norm_gene_count,
  selected_genes,
  sig_up_ervs,
  high_erv_covid
)

combined_low <- make_combined_matrix(
  norm_gene_count,
  selected_genes,
  sig_up_ervs,
  low_erv_covid
)

combined_healthy <- make_combined_matrix(
  norm_gene_count,
  selected_genes,
  sig_up_ervs,
  healthy_samples
)

# Convert ENSEMBL row names to symbols
gene_rows <- rownames(combined_high) %in% selected_genes

symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(combined_high)[gene_rows],
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

new_names <- rownames(combined_high)
new_names[gene_rows] <- ifelse(
  is.na(symbols),
  new_names[gene_rows],
  symbols
)

rownames(combined_high) <- make.unique(new_names)
rownames(combined_low) <- rownames(combined_high)
rownames(combined_healthy) <- rownames(combined_high)

feature_order <- make_corrplot(
  combined_high,
  output_file = file.path(
    outdir,
    "corrplot_highERV_Reactome_genes.pdf"
  )
)

make_corrplot(
  combined_low,
  feature_order = feature_order,
  output_file = file.path(
    outdir,
    "corrplot_lowERV_Reactome_genes.pdf"
  )
)

make_corrplot(
  combined_healthy,
  feature_order = feature_order,
  output_file = file.path(
    outdir,
    "corrplot_healthy_Reactome_genes.pdf"
  )
)

# Panel B: Inflammasome genes
inflammasome_genes <- c(
  "IL1B", "IL1R1", "IL1RAP", "TNF", "TNFRSF1A",
  "TLR4", "TLR2", "TLR8", "TICAM1", "MYD88",
  "IRAK1", "IRAK4", "FADD", "CASP8", "NOD1", "NOD2",
  "MAPK8", "PRKD1", "BRCC3", "NFKB1", "NLRP3",
  "CASR", "PLCG1", "P2RX7", "MAVS", "MFN2",
  "PYCARD", "CASP1", "GSDMD"
)

# Define this ERV vector before running Panel B
inflammasome_ervs <- c(
  "2124", "1874", "5302", "1179","3104", "3743"
)

combined_high_inf <- make_combined_matrix(
  norm_gene_count_name,
  inflammasome_genes,
  inflammasome_ervs,
  high_erv_covid
)

combined_low_inf <- make_combined_matrix(
  norm_gene_count_name,
  inflammasome_genes,
  inflammasome_ervs,
  low_erv_covid
)

combined_healthy_inf <- make_combined_matrix(
  norm_gene_count_name,
  inflammasome_genes,
  inflammasome_ervs,
  healthy_samples
)

feature_order_inf <- make_corrplot(
  combined_high_inf,
  output_file = file.path(
    outdir,
    "corrplot_highERV_inflammasome.pdf"
  ),
  width = 10,
  height = 10
)

make_corrplot(
  combined_low_inf,
  feature_order = feature_order_inf,
  output_file = file.path(
    outdir,
    "corrplot_lowERV_inflammasome.pdf"
  ),
  width = 10,
  height = 10
)

make_corrplot(
  combined_healthy_inf,
  feature_order = feature_order_inf,
  output_file = file.path(
    outdir,
    "corrplot_healthy_inflammasome.pdf"
  ),
  width = 10,
  height = 10
)