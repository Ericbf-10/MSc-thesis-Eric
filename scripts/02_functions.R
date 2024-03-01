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

# Convert Hugo Symbols to Ensembl Gene IDs
SYMBOLtoENSEMBL <- function(x){
  vec=mapIds(org.Hs.eg.db,
             keys=x,
             column="ENSEMBL",
             keytype="SYMBOL",
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
  }
  return(vec)
}

# Plot that summarizes the mutation count for each patient
sum_plot <- function(all_mut_df, nonsyn_mut_df, lof_mut_df, plot_path) {
  # Prepare a combined_df
  combined_df <- bind_rows(
    all_mut_df %>% mutate(source = 'All'),
    nonsyn_mut_df %>% mutate(source = 'Nonsyn'),
    lof_mut_df %>% mutate(source = 'LoF')
  ) %>%
    dplyr::select(sample_id, source)
  
  # Adjust the factor levels for 'source' to control the order
  combined_df$source <- factor(combined_df$source, levels = c("All", "Nonsyn", "LoF"))
  
  # Recalculate the counts
  count_df <- combined_df %>%
    group_by(sample_id, source) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the maximum count to set the y-axis limit dynamically
  max_count <- max(count_df$count, na.rm = TRUE)
  buffer <- max_count * 0.1 # Add 10% buffer to the maximum count for the labels
  
  # Generate the plot with the adjusted source ordering
  bar_plot <- ggplot(count_df, aes(x = source, y = count, fill = source)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = count), vjust = -0.5, position = position_dodge(0.9)) +
    scale_fill_manual(values = c("All" = "darkblue", "Nonsyn" = "darkred", "LoF" = "darkgreen"),
                      name = "Mutation Type") +
    facet_wrap(~ sample_id, scales = "free_x") +
    labs(x = "Mutation Type", 
         y = "Count", 
         title = "Number of mutations per patient") +
    ylim(0, max_count + buffer) + # Adjust the y-axis to include the buffer
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) # Attempt to visually center the title, though limitations exist
  
  # Save
  ggsave(plot_path, bar_plot, width = 15, height = 15, dpi = 300)
  
  # Return
  return(bar_plot)
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

# Filter for ubiquitin genes in a mutation list
ubiquitin_filter <- function(mut_list, ubiquitin_df){
  filtered_list <- lapply(names(mut_list), function(sample_id) {
    mut_list[[sample_id]] %>%
      dplyr::inner_join(ubiquitin_df, by = "ensembl_id")
  })
  
  # Set the names of the list elements to be the sample_id
  names(filtered_list) <- names(mut_list)
  
  return(filtered_list)
}

# Perform FORA analysis for a specific Human Gene Set collection
perform_fora_analysis <- function(all_mut_list, nonsyn_mut_list, lof_mut_list, fora_results_all_csv_path, fora_results_nosnyn_csv_path, fora_results_lof_csv_path, category=NULL, subcategory=NULL, gene_sets_list=NULL) {
  # Load required libraries
  library(msigdbr)
  library(fgsea)
  library(EnsDb.Hsapiens.v75)
  
  if (!is.null(category)){
    # Get a list of gene sets for the specified category with their associated genes
    gene_sets_df <- msigdbr(species = "human", category = category, subcategory = subcategory)
    gene_sets_list <- split(x = gene_sets_df$ensembl_gene, f = gene_sets_df$gs_name)
  } else if (is.null(gene_sets_list)) {
    # Raise exception & stop execution
    stop("category and gene_sets_list are both NULL. Please provide a 'category' or 'gene_sets_list'.")
  }
  
  # Get all the Ensembl human Gene IDs
  genes_info <- genes(EnsDb.Hsapiens.v75)
  ensembl_hgene_ids <- unique(genes_info$gene_id)
  
  # Initialize lists to store FORA results
  fora_results_all_list <- list()
  fora_results_nonsyn_list <- list()
  fora_results_lof_list <- list()
  
  for(sample_id in sample_id_list) {
    # Perform FORA for all mutations
    mutated_genes_all_ids <- all_mut_list[[sample_id]]$ensembl_id
    fora_results_all_list[[sample_id]] <- fora(gene_sets_list, mutated_genes_all_ids, ensembl_hgene_ids)
    
    # Perform FORA for nonsynonymous mutations
    mutated_genes_nonsyn_ids <- nonsyn_mut_list[[sample_id]]$ensembl_id
    fora_results_nonsyn_list[[sample_id]] <- fora(gene_sets_list, mutated_genes_nonsyn_ids, ensembl_hgene_ids)
    
    # Perform FORA for LoF mutations
    mutated_genes_lof_ids <- lof_mut_list[[sample_id]]$ensembl_id
    fora_results_lof_list[[sample_id]] <- fora(gene_sets_list, mutated_genes_lof_ids, ensembl_hgene_ids)
  }
  
  # Save the results to CSV files using the provided variable names
  fora_results_all_df <- merge_fora_data_tables_adding_id(fora_results_all_list)
  write.csv(fora_results_all_df, file = fora_results_all_csv_path, row.names = TRUE)
  
  fora_results_nonsyn_df <- merge_fora_data_tables_adding_id(fora_results_nonsyn_list)
  write.csv(fora_results_nonsyn_df, file = fora_results_nosnyn_csv_path, row.names = TRUE)
  
  fora_results_lof_df <- merge_fora_data_tables_adding_id(fora_results_lof_list)
  write.csv(fora_results_lof_df, file = fora_results_lof_csv_path, row.names = TRUE)
}

# Plot top 5 Gene Sets across samples stratified by response
plot_top5 <- function(response_df, mut_df, plot_path, fig_title, fig_y_label) {
  
  # Prepare data for all samples in the df: Filter padjust < 0.5 & calculate Gene ratio per sample
  prepared_data_df <- mut_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::mutate(Gene_ratio = overlap / size,
                  Count = overlap) %>%
    dplyr::slice_min(order_by = padj, n = 5) %>%
    dplyr::arrange(sample_id, desc(Count)) %>%
    dplyr::ungroup()

  # Merge response info
  combined_data <- prepared_data_df %>%
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
  
  # Calculate dynamic breaks and labels for Count
  count_breaks <- quantile(combined_data$Count, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  count_breaks <- round(count_breaks) # Ensure count_breaks are integers by rounding
  count_breaks <- count_breaks[count_breaks > 0] # Remove 0 values from count_breaks
  count_breaks <- unique(count_breaks) # Make sure breaks are unique
  
  # Dynamically generate count_labels based on the available count_breaks
  count_labels <- vector("character", length(count_breaks))
  for (i in seq_along(count_breaks)) {
    count_labels[i] <- as.character(count_breaks[i])
  }
  
  # Step 2: Plot with ordered sample_id and facets to distinguish "R" and "NR"
  top5_plot <- ggplot(combined_data, aes(x = sample_id,
                                              y = reorder(pathway, Count),
                                              size = Count,
                                              color = padj)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red", trans = "log10", # Log10 transformation for a better visualization
                         limits = c(1e-6, 1), oob = scales::oob_squish) +
    scale_size_continuous(name = "Count",
                          range = c(1, 6),  # Adjust the visual size range as needed
                          breaks = count_breaks,
                          labels = count_labels) +
    facet_grid(. ~ patient_response, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = fig_title,
         x = "Sample ID",
         y = fig_y_label,
         color = "p.adjust (log10)",
         size = "Count")

  # Save plot
  ggsave(plot_path, top5_plot, width = 15, height = 20, dpi = 300)
  
  # Return
  return(top5_plot)
}

# Plot (a) specific Gene Set(s) across samples stratified by response
plot_geneset <- function(geneset, response_df, mut_df, plot_path, fig_title, fig_y_label) {
  # geneset is a vector i.e. c("GS1", "GS2")
  
  # Prepare data for all samples in the df: Filter Gene Set & calculate Gene ratio per sample
  prepared_data_df <- mut_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::filter(pathway %in% geneset) %>% 
    dplyr::mutate(Gene_ratio = overlap / size,
                  Count = overlap,
                  sample_id = sample_id) %>%
    dplyr::arrange(desc(Count)) %>%
    dplyr::ungroup()
  
  # Merge response info
  combined_data <- prepared_data_df %>%
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
  
  # Calculate dynamic breaks and labels for Count
  count_breaks <- quantile(combined_data$Count, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  count_breaks <- round(count_breaks) # Ensure count_breaks are integers by rounding
  count_breaks <- count_breaks[count_breaks > 0] # Remove 0 values from count_breaks
  count_breaks <- unique(count_breaks) # Make sure breaks are unique
  
  # Dynamically generate count_labels based on the available count_breaks
  count_labels <- vector("character", length(count_breaks))
  for (i in seq_along(count_breaks)) {
    count_labels[i] <- as.character(count_breaks[i])
  }
  
  # Step 2: Plot with ordered sample_id and facets to distinguish "R" and "NR"
  geneset_plot <- ggplot(combined_data, aes(x = sample_id, 
                                                    y = reorder(pathway, Count),
                                                    size = Count,
                                                    color = padj)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red", trans = "log10",
                         limits = c(1e-6, 1), oob = scales::oob_squish) +
    scale_size_continuous(name = "Count",
                          range = c(1, 6),  # Adjust the visual size range as needed
                          breaks = count_breaks,
                          labels = count_labels) +
    facet_grid(. ~ patient_response, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = fig_title,
         x = "Sample ID",
         y = fig_y_label,
         color = "p.adjust (log10)",
         size = "Count")
  
  # Save plot
  ggsave(plot_path, geneset_plot, width = 12, height = 8, dpi = 300)
  
  # Return
  return(geneset_plot)
}