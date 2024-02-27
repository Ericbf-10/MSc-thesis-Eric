########################
### FUNCTIONS SCRIPT ###
########################

### Functions needed for the project

# Set data.frame columns of type list to character
set_lists_to_chars <- function(x) {
  if(class(x) == 'list') {
    y <- paste(unlist(x[1]), sep='', collapse=', ')
  } else {
    y <- x 
  }
  return(y)
}

# Transform a list of data tables coming from the fora results into a big data.frame maintaining the id
merge_fora_data_tables_adding_id <- function(data_tables_list) {
  # Check if input is a list
  if (!is.list(data_tables_list)) {
    stop("The input must be a list of data tables or data frames.")
  }
  
  # Iterate over the list to modify each data table
  modified_list <- lapply(names(data_tables_list), function(id) {
    # Ensure the element is a data frame or data table
    if (!is.data.frame(data_tables_list[[id]])) {
      stop("All elements of the list must be data frames or data tables.")
    }
    
    # Copy the data table to avoid modifying the original list in place
    dt_copy <- data_tables_list[[id]]
    # Add a new column with the ID
    dt_copy$sample_id <- id
    return(dt_copy)
  })
  
  # Combine all modified data tables into one dataframe
  combined_df <- do.call(rbind, modified_list)
  combined_df <- combined_df %>% 
    mutate(overlapGenes = map_chr(overlapGenes, ~ paste(.x, collapse = ", ")))
  return(combined_df)
}

# Convert Hugo Symbols to Entrez IDs
SYMBOLtoENTREZ <- function(x){
  vec=mapIds(org.Hs.eg.db,
             keys=x,
             column="ENTREZID",
             keytype="SYMBOL",
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
    }
  return(vec)
}

# Convert Entrez IDs to Hugo Symbols
ENTREZtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,
             keys=x,
             column="SYMBOL",
             keytype="ENTREZID",
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
  }
  return(vec)
}

# Convert Entrez IDs to Ensembl Gene IDs
ENTREZtoENSEMBL <- function(x){
  vec=mapIds(org.Hs.eg.db,
             keys=x,
             column="ENSEMBL",
             keytype="ENTREZID",
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
    }
  return(vec)
}

# Plot top 5 Gene Sets across samples stratified by response
plot_top5 <- function(response_df, mut_list, plot_path, fig_title, fig_y_label) {
  
  # Prepare data for all samples in the list: Filter padjust < 0.5 & calculate Gene ratio per sample from a mut_list
  prepared_data_list <- lapply(names(mut_list), function(sample_id) {
    mut_list[[sample_id]] %>%
      dplyr::filter(padj < 0.05) %>% 
      dplyr::mutate(Gene_ratio = overlap / size,
                    Count = overlap,
                    sample_id = sample_id) %>%
      dplyr::slice_min(order_by = padj, n = 5) %>%
      dplyr::arrange(desc(Gene_ratio))
  })

  # Combine all prepared datasets into one
  combined_data <- bind_rows(prepared_data_list)

  # Merge response info
  combined_data <- combined_data %>%
    dplyr::left_join(response_df, by = "sample_id")

  # Step 1: Create an ordered factor for sample_id based on patient_response
  # First, create a mapping of sample_id to patient_response
  response_order <- combined_data %>%
    dplyr::select(sample_id, patient_response) %>%
    dplyr::distinct() %>%
    dplyr::arrange(patient_response, sample_id) # Arrange so that "R" comes before "NR"

  # Use this ordering to reorder sample_id in the combined data
  combined_data$sample_id <- factor(combined_data$sample_id,
                                    levels = response_order$sample_id)

  # Determine dynamic breaks based on quantiles of Gene Ratio
  gene_ratio_breaks <- quantile(combined_data$Gene_ratio, probs = seq(0, 1, by = 0.5), na.rm = TRUE)

  # Step 2: Plot with ordered sample_id and facets to distinguish "R" and "NR"
  top5_plot <- ggplot(combined_data, aes(x = sample_id,
                                              y = reorder(pathway, Gene_ratio),
                                              size = Gene_ratio,
                                              color = padj)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red", trans = "log10",
                         limits = c(1e-6, 1), oob = scales::oob_squish) +
    scale_size_continuous(name = "Gene Ratio",
                          range = c(1, 6),  # Adjust size range to match your preference
                          breaks = gene_ratio_breaks,  # Define breaks based on "Gene Ratio" distribution
                          labels = format(gene_ratio_breaks, digits = 1)) +  # Label breaks as needed
    facet_grid(. ~ patient_response, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", colour = "black"),
      panel.grid.major = element_line(color = "grey", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = fig_title,
         x = "Sample ID",
         y = fig_y_label,
         color = "p.adjust (log10)") +
    guides(size = guide_legend(title = "Gene Ratio"))

  # Save plot
  ggsave(plot_path, top5_plot, width = 12, height = 8, dpi = 300)
  
  # Return
  return(top5_plot)
}