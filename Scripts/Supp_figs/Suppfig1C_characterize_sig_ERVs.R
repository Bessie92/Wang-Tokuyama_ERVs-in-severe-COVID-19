library(dplyr)
library(tibble)
library(ggplot2)

signature_ervs <- read.csv("../../Data/GSE171110/define_signature_ervs/Raw_Table1_Signature_erv_overlap_summary.csv", row.names = 1)
Blomberg_annot <- read.csv("../../Assets/Blomberg_hg38_lifted.csv")

B_sig_ervs <- Blomberg_annot[Blomberg_annot$rvnr %in% signature_ervs$erv_name, ]
B_sig_ervs <- B_sig_ervs[, c("rvnr", "subgenes", "noncanon", "canon", "supergroup", "class", "comment")]
B_sig_ervs <- B_sig_ervs %>% 
  mutate(rvnr = as.character(rvnr)) %>%
  add_row(rvnr = "W-8", canon = "HERVW") %>%
  dplyr::rename(erv_name = rvnr)

B_sig_ervs <- left_join(B_sig_ervs, signature_ervs)
# Combine canon and noncanon columns into ERV_group
B_sig_ervs <- B_sig_ervs %>%
  mutate(
    noncanon = na_if(noncanon, ""),
    canon    = na_if(canon, "")
  )%>%
  mutate(erv_group = coalesce(noncanon, canon))
write.csv(B_sig_ervs, "../../Data/GSE171110/signature_ERV_blomberg_feature_table.csv", row.names = F)

# 1. Histogram to visualize chromosome distribution 
## Order the chromosomes
chrom_order <- c(as.character(1:22), "X", "Y")
B_sig_ervs$chromosome <- factor(B_sig_ervs$chromosome, levels = chrom_order)
## Plot
ggplot(B_sig_ervs, aes(x = chromosome)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "ERV count per chromosome", x = "Chromosome", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./figures/ERV signature count per chromosome.pdf", width = 7, height = 5)

# 2. Histogram to visualize ERV locus subfamily distribution
ggplot(B_sig_ervs, aes(x = erv_group)) +
  geom_bar(fill = "darkorange") +
  theme_minimal() +
  labs(title = "ERV count per ERV group", x = "Group", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./figures/ERV signature group distribution.pdf", width = 7, height = 5)

# 3. histogram to visualize ERV locus length distribution
B_sig_ervs$length <- B_sig_ervs$end - B_sig_ervs$start
# Use bins to show distribution 
B_sig_ervs$length_group <- cut(
  B_sig_ervs$length,
  breaks = c(0, 5000, 7500, 10000, 12500, 15000, 20000),
  labels = c("0-5k", "5-7.5k", "7.5-10k", "10-12.5k", "12.5-15k", "15-20k"),
  right = FALSE  # interval notation: [a, b)
)

# Bar plot
ggplot(B_sig_ervs, aes(x = length_group)) +
  geom_bar(fill = "darkorange") +
  theme_minimal() +
  labs(title = "ERV Length Distribution",
       x = "Length Range (bp)",
       y = "Number of ERVs") +
  theme(axis.text.x = element_text(size = 10))
ggsave("./figures/ERV signature Length Distribution by Range.pdf", width = 7, height = 5)

# Calculate the average length
mean_length <- mean(B_sig_ervs$length)
mean_length
