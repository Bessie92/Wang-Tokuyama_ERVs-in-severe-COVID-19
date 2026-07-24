library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)

norm_gene_count <- read.csv("../../Data/GSE163317/normalized_genes.csv",row.names = 1)
norm_erv_count <- read.csv("../../Data/GSE163317/normalized_ervs.csv",row.names = 1)
coldata <- read.csv("../../Data/GSE163317/coldata.csv",row.names = 1)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

ervs <- c("2724", "2124", "4745", "942", "4896", "4830")
il1_genes <- c("IL1R1", "IL1R2", "IL1RAP", "IRAK3", "PELI1")

inflammasome_genes <- c(
  # Priming / signal 1
  "IL1B", "IL1R1", "IL1RAP", "TNF", "TNFRSF1A",
  "TLR4", "TLR2", "TLR8", "TICAM1", "MYD88",
  "IRAK1", "IRAK4", "FADD", "CASP8", "NOD1", "NOD2",
  "MAPK8", "PRKD1", "BRCC3", "NFKB1", "NLRP3",
  
  # Activation / signal 2
  "CASR", "PLCG1", "P2RX7", "MAVS", "MFN2",
  "PYCARD", "CASP1", "GSDMD"
)

# Helpers

# Map gene symbols to ENSEMBL IDs present in the count matrix
map_symbols_to_ensembl <- function(symbols, count_mat) {
  
  count_ids <- sub("\\..*$", "", rownames(count_mat))
  
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENSEMBL")
  ) %>%
    dplyr::filter(
      !is.na(ENSEMBL),
      ENSEMBL %in% count_ids
    ) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

# Convert an expression matrix to long format and attach metadata
expression_to_long <- function(
    count_mat,
    features,
    feature_col,
    metadata,
    feature_levels = features
) {
  
  features_present <- features[
    features %in% rownames(count_mat)
  ]
  
  count_mat[
    features_present,
    ,
    drop = FALSE
  ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(feature_col) %>%
    tidyr::pivot_longer(
      cols = -all_of(feature_col),
      names_to = "sample",
      values_to = "normalized_count"
    ) %>%
    dplyr::left_join(
      metadata %>%
        tibble::rownames_to_column("sample"),
      by = "sample"
    ) %>%
    dplyr::mutate(
      "{feature_col}" := factor(
        .data[[feature_col]],
        levels = feature_levels
      )
    )
}

# Convert genes to long format
gene_expression_to_long <- function(
    count_mat,
    symbols,
    metadata
) {
  
  mapping <- map_symbols_to_ensembl(
    symbols,
    count_mat
  )
  
  count_ids <- sub("\\..*$", "", rownames(count_mat))
  rownames(count_mat) <- count_ids
  
  count_mat[
    mapping$ENSEMBL,
    ,
    drop = FALSE
  ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ENSEMBL") %>%
    tidyr::pivot_longer(
      cols = -ENSEMBL,
      names_to = "sample",
      values_to = "normalized_count"
    ) %>%
    dplyr::left_join(
      mapping %>%
        dplyr::select(ENSEMBL, SYMBOL),
      by = "ENSEMBL"
    ) %>%
    dplyr::left_join(
      metadata %>%
        tibble::rownames_to_column("sample"),
      by = "sample"
    ) %>%
    dplyr::mutate(
      SYMBOL = factor(
        SYMBOL,
        levels = symbols
      )
    )
}

# Calculate the mean feature z-score per sample
calculate_score <- function(
    long_df,
    feature_col,
    score_name
) {
  
  long_df %>%
    dplyr::group_by(
      .data[[feature_col]]
    ) %>%
    dplyr::mutate(
      feature_z = as.numeric(
        scale(log2(normalized_count + 1))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
      sample,
      Patient,
      Day
    ) %>%
    dplyr::summarise(
      score = mean(feature_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(
      "{score_name}" := score
    )
}

# Create a patient spaghetti plot
make_spaghetti_plot <- function(
    data,
    value_col,
    title,
    y_label
) {
  
  ggplot(
    data,
    aes(
      x = Day,
      y = .data[[value_col]],
      group = Patient,
      color = Patient
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      size = 2.8,
      na.rm = TRUE
    ) +
    labs(
      title = title,
      x = "Treatment day",
      y = y_label,
      color = "Patient"
    ) +
    theme_classic(base_size = 12)
}


erv_long <- expression_to_long(
  count_mat = norm_erv_count,
  features = ervs,
  feature_col = "erv",
  metadata = coldata
)

il1_long <- gene_expression_to_long(
  count_mat = norm_gene_count,
  symbols = il1_genes,
  metadata = coldata
)

inflammasome_long <- gene_expression_to_long(
  count_mat = norm_gene_count,
  symbols = inflammasome_genes,
  metadata = coldata
)


erv_score_long <- calculate_score(
  long_df = erv_long,
  feature_col = "erv",
  score_name = "ERV_score"
)

il1_score_long <- calculate_score(
  long_df = il1_long,
  feature_col = "SYMBOL",
  score_name = "IL1_score"
)

inflammasome_score_long <- calculate_score(
  long_df = inflammasome_long,
  feature_col = "SYMBOL",
  score_name = "Inflammasome_score"
)

# Top panel: spaghetti plots
p_erv_score <- make_spaghetti_plot(
  data = erv_score_long,
  value_col = "ERV_score",
  title = "Combined ERV expression score",
  y_label = "Mean feature z-score"
)

p_il1_score <- make_spaghetti_plot(
  data = il1_score_long,
  value_col = "IL1_score",
  title = "Combined IL-1 signaling score",
  y_label = "Mean feature z-score"
)

p_inflammasome_score <- make_spaghetti_plot(
  data = inflammasome_score_long,
  value_col = "Inflammasome_score",
  title = "Combined inflammasome activation score",
  y_label = "Mean feature z-score"
)

crp_long <- coldata %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::select(
    sample,
    Patient,
    Day,
    CRP
  )

ldh_long <- coldata %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::select(
    sample,
    Patient,
    Day,
    LDH
  )

p_crp <- make_spaghetti_plot(
  data = crp_long,
  value_col = "CRP",
  title = "CRP from Day 0 to Day 7",
  y_label = "CRP"
)

p_ldh <- make_spaghetti_plot(
  data = ldh_long,
  value_col = "LDH",
  title = "LDH from Day 0 to Day 7",
  y_label = "LDH"
)

p_summary <- patchwork::wrap_plots(
  p_erv_score +
    coord_cartesian(ylim = c(-1.7, 1.3)),
  
  p_il1_score +
    coord_cartesian(ylim = c(-1.7, 1.3)),
  
  p_inflammasome_score +
    coord_cartesian(ylim = c(-1.7, 1.3)),
  
  p_crp,
  p_ldh,
  nrow = 1
)

p_summary

ggsave(
  filename = "figures/combined_score.pdf",
  plot = p_summary,
  height = 5,
  width = 12
)

# Bottom panel: Day 7 minus Day 0 changes

molecular_long <- dplyr::bind_rows(
  erv_score_long %>%
    dplyr::transmute(
      sample,
      Patient,
      Day = as.integer(as.character(Day)),
      measurement = "ERV score",
      value = ERV_score
    ),
  
  il1_score_long %>%
    dplyr::transmute(
      sample,
      Patient,
      Day = as.integer(as.character(Day)),
      measurement = "IL-1 score",
      value = IL1_score
    ),
  
  inflammasome_score_long %>%
    dplyr::transmute(
      sample,
      Patient,
      Day = as.integer(as.character(Day)),
      measurement = "Inflammasome score",
      value = Inflammasome_score
    )
)

clinical_long <- coldata %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::select(
    sample,
    Patient,
    Day,
    CRP,
    Ferritin,
    LDH
  ) %>%
  dplyr::mutate(
    Day = as.integer(as.character(Day))
  ) %>%
  tidyr::pivot_longer(
    cols = c(CRP, Ferritin, LDH),
    names_to = "measurement",
    values_to = "raw_value"
  ) %>%
  dplyr::group_by(measurement) %>%
  dplyr::mutate(
    value = as.numeric(
      scale(log2(raw_value + 1))
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    sample,
    Patient,
    Day,
    measurement,
    value
  )

measurement_long <- dplyr::bind_rows(
  molecular_long,
  clinical_long
)

measurement_order <- c(
  "ERV score",
  "IL-1 score",
  "Inflammasome score",
  "CRP",
  "Ferritin",
  "LDH"
)

measurement_delta <- measurement_long %>%
  dplyr::select(
    Patient,
    Day,
    measurement,
    value
  ) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from = Day,
    values_from = value,
    names_prefix = "Day"
  ) %>%
  dplyr::mutate(
    delta = Day7 - Day0,
    measurement = factor(
      measurement,
      levels = measurement_order
    )
  )

p_delta_patient <- ggplot(
  measurement_delta,
  aes(
    x = measurement,
    y = delta
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_col(
    width = 0.7,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ Patient,
    nrow = 1
  ) +
  labs(
    title = "Within-patient changes from Day 0 to Day 7",
    x = NULL,
    y = expression(
      Delta * " standardized value (Day 7 - Day 0)"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

p_delta_patient

ggsave(
  filename = "figures/combined_score_delta.pdf",
  plot = p_delta_patient,
  height = 5,
  width = 8
)