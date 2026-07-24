run_DESeq2_paired <- function(
    countData,
    coldata,
    paired_var,
    condition_var,
    reference_level,
    test_level,
    output_prefix,
    keep_limit = 0,
    baseMean_cutoff = 10,
    extra_design_vars = NULL
) {
  
  library(DESeq2)
  
  # Convert count table to matrix
  cts <- as.matrix(countData)
  
  # Check count matrix and metadata alignment
  if (!all(rownames(coldata) == colnames(cts))) {
    stop("Data unaligned: rownames(coldata) must match colnames(countData)")
  }
  
  # Convert paired/blocking variable to factor
  coldata[[paired_var]] <- factor(coldata[[paired_var]])
  
  # Convert condition variable to factor with specified reference level
  coldata[[condition_var]] <- factor(
    coldata[[condition_var]],
    levels = c(reference_level, test_level)
  )
  
  # Build design formula
  # Basic paired design: ~ paired_var + condition_var
  if (is.null(extra_design_vars)) {
    design_terms <- c(paired_var, condition_var)
  } else {
    design_terms <- c(paired_var, extra_design_vars, condition_var)
  }
  
  design_formula <- as.formula(
    paste("~", paste(design_terms, collapse = " + "))
  )
  
  message("Using DESeq2 design: ", deparse(design_formula))
  message("Testing contrast: ", condition_var, " ", test_level, " vs ", reference_level)
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design_formula
  )
  
  # Pre-filter low-count features before DESeq2
  keep <- rowSums(counts(dds)) > keep_limit
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results
  contrast_vector <- c(condition_var, test_level, reference_level)
  res <- results(dds, contrast = contrast_vector)
  
  # Convert to data frame
  res_df <- as.data.frame(res)
  res_df$feature <- rownames(res_df)
  
  # Save unfiltered results
  unfiltered_file <- paste0(
    output_prefix,
    "_DESeq2_",
    test_level,
    "_vs_",
    reference_level,
    "_unfiltered.txt"
  )
  
  write.table(
    res_df,
    file = unfiltered_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Filter by baseMean
  res_filtered <- res_df[
    !is.na(res_df$baseMean) & res_df$baseMean >= baseMean_cutoff,
  ]
  
  # Save filtered results
  filtered_file <- paste0(
    output_prefix,
    "_DESeq2_",
    test_level,
    "_vs_",
    reference_level,
    "_baseMean",
    baseMean_cutoff,
    "_filtered.txt"
  )
  
  write.table(
    res_filtered,
    file = filtered_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Saved unfiltered results to: ", unfiltered_file)
  message("Saved filtered results to: ", filtered_file)
  
  return(list(
    dds = dds,
    design = design_formula,
    contrast = contrast_vector,
    unfiltered = res_df,
    filtered = res_filtered
  ))
}