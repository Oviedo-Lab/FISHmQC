
# Function to load summary stats
summary_stats <- function(
    id,
    path
  ) {
    # id = m1_blanks, m1_genes, m2_blanks, etc ...
    pattern_var <- paste0(id, ".*\\.csv$")
    files <- list.files(
      path = path,
      pattern = pattern_var,
      full.names = TRUE
    )
    if (length(files) > 1) stop("Multiple files found for pattern: ", pattern_var)
    return(read.csv(files))
  }

# This function can be used to extract summary statistics from your own data
# ...  it will return a list holding two data.frames formatted like those imported below from the csv files
summary_stats_from_data <- function(
    obs_data,         # data frame, columns as barcodes (genes or blanks), rows as cells, elements as spot counts
    blank_mask,       # logical vector with length equal to number of columns in obs_data, TRUE for blank barcodes, FALSE for gene barcodes
    cell_mask = TRUE  # optional row mask, logical vector with length equal to the number of rows in obs_data, rows marked FALSE will be discarded
  ) {
    
    # Parse provided data
    obs_data_genes <- obs_data[cell_mask, !blank_mask]
    obs_data_blanks <- obs_data[cell_mask, blank_mask]
   
    # Extract summary statistics
    observed_gene_rates <- colMeans(obs_data_genes)
    observed_blank_rates <- colMeans(obs_data_blanks)
    observed_gene_variance <- rep(NA, length(observed_gene_rates))
    observed_blank_variance <- rep(NA, length(observed_blank_rates))
    for (i in 1:length(observed_gene_variance)) {observed_gene_variance[i] <- sd(obs_data_genes[,i])^2}
    for (i in 1:length(observed_blank_variance)) {observed_blank_variance[i] <- sd(obs_data_blanks[,i])^2}
    observed_gene_counts <- colSums(obs_data_genes)
    observed_blank_counts <- colSums(obs_data_blanks)
    
    # Collect summary stats
    gene_data <- data.frame(
      rates = observed_gene_rates,
      variance = observed_gene_variance,
      counts = observed_gene_counts
    )
    rownames(gene_data) <- colnames(obs_data_genes)
    blank_data <- data.frame(
      rates = observed_blank_rates,
      variance = observed_blank_variance,
      counts = observed_blank_counts
    )
    rownames(blank_data) <- colnames(obs_data_blanks)
    
    return(
      list(gene = gene_data, blank = blank_data)
    )
    
  }

# ... This was produced from the parse_hdf5 and make_count_data functions, using the two good P18 WT mice and 
#      the P18 Mecp2 mouse. The entire brain is included, with L1, and with all genes and all blanks.

make_csv_from_hdf5 <- function(
    hdf5_data_path
  ) {
    
    countdata <- make_count_data(
      data_path = hdf5_data_path, 
      remove_L1 = FALSE,   # Exclude transcripts in layer 1?
      ROIname = "Primary auditory area",
      raw = TRUE,             # Grab normalized transcript counts, or raw ones?
      verbose = TRUE,
      drop_blanks = FALSE
    )$count_data
    
    mice <- unique(countdata$mouse)
    
    for (m in mice) {
      
      data_m <- countdata[countdata$mouse == m,11:829]
      blank_mask_m <- c(rep(FALSE, 829-98-10), rep(TRUE, 98))
      
      summary_data_m <- summary_stats_from_data(
        obs_data = data_m,
        blank_mask = blank_mask_m,
        cell_mask = TRUE
      )
      gene_data_m <- summary_data_m$gene
      blank_data_m <- summary_data_m$blank
      
      write.csv(
        gene_data_m,
        file = paste0("summary_stats/summary_stats_m", m, "_genes.csv"),
        row.names = FALSE
      )
      write.csv(
        blank_data_m,
        file = paste0("summary_stats/summary_stats_m", m, "_blanks.csv"),
        row.names = FALSE
      )
      
    }
  }



