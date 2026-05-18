
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
    obs_data         # data frame, columns as barcodes (genes or blanks), rows as cells, elements as spot counts
  ) {
    
    # Extract summary statistics
    observed_rates <- colMeans(obs_data)
    observed_variance <- rep(NA, length(observed_rates))
    for (i in 1:length(observed_variance)) {observed_variance[i] <- sd(obs_data[,i])^2}
    observed_counts <- colSums(obs_data)
    
    # Collect summary stats
    data <- data.frame(
      rates = observed_rates,
      variance = observed_variance,
      counts = observed_counts
    )
    rownames(data) <- colnames(obs_data)
    
    return(data)
    
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



