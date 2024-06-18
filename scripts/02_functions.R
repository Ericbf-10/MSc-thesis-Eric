########################
### FUNCTIONS SCRIPT ###
########################

### Functions needed for the project

# Color blind friendly palette
cbPalette <- c("#56B4E9", "#999999", "#8A2BE2", "#E69F00", "#D55E00", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#FF69B4", "#FFD700", "#000000", "#00CED1", "#4B0082")

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
  vec=mapIds(EnsDb.Hsapiens.v75,
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
  vec=mapIds(EnsDb.Hsapiens.v75,
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
  vec=mapIds(EnsDb.Hsapiens.v75,
             keys=x,
             column="GENEID", # = ENSEMBL
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
  vec=mapIds(EnsDb.Hsapiens.v75,
             keys=x,
             column="GENEID", # = ENSEMBL
             keytype="SYMBOL",
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
  }
  return(vec)
}

# Convert Hugo Symbols to Ensembl Gene IDs
SYMBOLtoENSEMBL38 <- function(x){
  vec=mapIds(EnsDb.Hsapiens.v86,
             keys=x,
             column="GENEID", # = ENSEMBL
             keytype="SYMBOL",
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
  }
  return(vec)
}

# Convert Ensembl Gene IDs to Hugo Symbols
ENSEMBLtoSYMBOL <- function(x){
  vec=mapIds(EnsDb.Hsapiens.v75,
             keys=x,
             column="SYMBOL",
             keytype="GENEID", # = ENSEMBL
             multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){
    vec[y]=x[y]
  }
  return(vec)
}

# Define a function to correct date-like HUGO symbols
correct_date_like_symbols <- function(gene_ids) {
  sapply(gene_ids, function(gene_id) {
    if (str_detect(gene_id, "^\\w{3}-\\d{2}$")) {  # Detects patterns like Mar-01
      parts <- unlist(str_split(gene_id, "-"))
      # Reverse the parts and remove leading zeros from the day part
      corrected_id <- paste0(as.integer(parts[2]), "-", parts[1])
      return(corrected_id)
    } else {
      return(gene_id)  # Return the original ID if it doesn't match the pattern
    }
  })
}

# Function to make the patients characteristics plot (patchwork)
make_annotation_plot <- function(df, var_name, title, col_start, custom_labels, show_sample_id = FALSE) {
  plot <- ggplot(df, aes_string(x = reorder(df$sample_id, df$mut_load), y = "1", fill = var_name)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_fill_manual(values=cbPalette[col_start:length(cbPalette)], labels = custom_labels) +
    labs(x = NULL, y = NULL, fill = title) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right")
  
  # Conditionally add sample ID labels
  if (show_sample_id) {
    plot <- plot + 
      geom_text(aes(label = reorder(df$sample_id, df$mut_load)), y = 0, vjust = 0, size = 3) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }
  
  return(plot)
}

# Function to search GO terms by keyword
searchGOTerms <- function(keyword) {
  all_terms <- as.list(GOTERM)
  matches <- sapply(all_terms, function(term) grepl(keyword, Term(term), ignore.case = TRUE))
  matched_terms <- all_terms[matches]
  data.frame(
    go_id = names(matched_terms),
    go_term = sapply(matched_terms, Term),
    definition = sapply(matched_terms, Definition),
    go_category = sapply(matched_terms, Ontology),
    stringsAsFactors = FALSE
  )
}

# Bar Plot that summarizes the mutation count for each patient
sum_plot <- function(response_df, all_mut_df, nonsyn_mut_df, lof_mut_df, plot_path, select_r = FALSE) {
  # Ensure all sample IDs are included, even those with no mutations
  all_sample_ids <- unique(c(all_mut_df$sample_id, nonsyn_mut_df$sample_id, lof_mut_df$sample_id, response_df$sample_id))
  
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
    dplyr::group_by(sample_id, source) %>%
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    tidyr::complete(sample_id = all_sample_ids, source, fill = list(count = 0)) %>% # Ensure every sample_id and source combination is represented
    dplyr::ungroup()
  
  # Merge response info
  combined_data <- count_df %>%
    dplyr::left_join(response_df, by = "sample_id")
  
  if (select_r) {
    combined_data <- combined_data %>% 
      dplyr::filter(patient_response == "R")
    plot_title = "Number of mutations in Responders"
  } else {
    combined_data <- combined_data %>% 
      dplyr::filter(patient_response == "NR")
    plot_title = "Number of mutations in Non-responders"
  }
  
  # Define order R > NR
  response_order <- combined_data %>%
    dplyr::select(sample_id, patient_response) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(patient_response)) # Arrange so that "R" comes before "NR"
  
  # Step 1: Create an ordered factor in 'response_order'
  response_order <- response_order %>%
    dplyr::mutate(order = as.integer(factor(sample_id, levels = unique(sample_id))))
  
  # Step 2: Join this order to 'combined_data'
  # Step 3: Arrange 'combined_data' using the new order and then by 'source' to keep its internal order
  ordered_data <- combined_data %>%
    dplyr::left_join(response_order %>% 
                       dplyr::select(sample_id, order), by = "sample_id") %>% 
    dplyr::arrange(order, source) %>%
    dplyr::select(-order) %>% # Optionally remove the 'order' column if it's no longer needed
    dplyr::mutate(sample_id = factor(sample_id, levels = unique(sample_id)))
  
  # Calculate the maximum count and buffer again just in case
  max_count <- max(ordered_data$count, na.rm = TRUE)
  buffer <- max_count * 0.1
  
  # Generate the plot with the adjusted source ordering
  bar_plot <- ggplot(ordered_data, aes(x = source, y = count, fill = source)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = count), vjust = -0.5, position = position_dodge(0.9)) +
    scale_fill_manual(values = c("All" = "#0072B2", "Nonsyn" = "#CC79A7", "LoF" = "#FFD700"),
                      labels = c("All" = "All", "Nonsyn" = "Nonsynonymous", "LoF" = "Loss of Function"),
                      name = "Mutation type") +
    facet_wrap(~ sample_id, scales = "free_x") +
    labs(x = NULL, 
         y = "Mutation count", 
         title = plot_title) +
    ylim(0, max_count + buffer) + # Adjust the y-axis to include the buffer
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) # Attempt to visually center the title, though limitations exist
  
  # Save
  ggsave(plot_path, bar_plot, width = 10, height = 8, dpi = 300)
  
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
perform_fora_analysis <- function(sample_id_list, all_mut_list, nonsyn_mut_list, lof_mut_list, fora_results_all_csv_path, fora_results_nosnyn_csv_path, fora_results_lof_csv_path, category=NULL, subcategory=NULL, gene_sets_list=NULL) {
  
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

# Perform FORA analysis for a specific Human Gene Set collection using clusterProfiler
perform_cP_fora_analysis <- function(sample_id_list, all_mut_list, nonsyn_mut_list, lof_mut_list, fora_results_all_csv_path, fora_results_nosnyn_csv_path, fora_results_lof_csv_path, category=NULL, subcategory=NULL, gene_sets_df=NULL) {
  
  if (!is.null(category)){
    # Get a list of gene sets for the specified category with their associated genes
    gene_sets_df <- msigdbr(species = "human", category = category, subcategory = subcategory)
  } else if (is.null(gene_sets_df)) {
    # Raise exception & stop execution
    stop("category and gene_sets_list are both NULL. Please provide a 'category' or 'gene_sets_df'.")
  }
  
  # Get all the Ensembl human Gene IDs
  genes_info <- genes(EnsDb.Hsapiens.v75)
  ensembl_hgene_ids <- unique(genes_info$gene_id)
  
  # Initialize dataframes to store FORA results
  fora_results_all_df <- fora_results_nonsyn_df <- fora_results_lof_df <- data.frame() # WOW, didn't know this was possible in R!
  
  # Initialize lists to store FORA results
  fora_results_all_list <- fora_results_nonsyn_list <- fora_results_lof_list <- list()

  for(sample_id in sample_id_list) {
    # Perform FORA for all mutations
    mutated_genes_all_ids <- as.character(all_mut_list[[sample_id]]$ensembl_id)
    fora_results_all <- enricher(
      gene = mutated_genes_all_ids, # A character vector of your genes of interest
      pvalueCutoff = 1, # Can choose a FDR cutoff
      pAdjustMethod = "BH", # Method to be used for multiple testing correction
      universe = ensembl_hgene_ids, # A character vector containing your background set genes
      TERM2GENE = dplyr::select( # Pathway information: data frame with a term name or identifier and the gene identifiers
        gene_sets_df,
        gene_set_name,
        ensembl_id
      )
    )
    # Store results in a list
    fora_results_all_list[[sample_id]] <- fora_results_all
    
    # Add sample_id
    temp_fora_results_all_df <- as.data.frame(fora_results_all@result)
    temp_fora_results_all_df$sample_id <- sample_id
    fora_results_all_df <- rbind(fora_results_all_df, temp_fora_results_all_df)
    
    # Perform FORA for nonsynonymous mutations
    mutated_genes_nonsyn_ids <- as.character(nonsyn_mut_list[[sample_id]]$ensembl_id)
    fora_results_nonsyn <- enricher(
      gene = mutated_genes_nonsyn_ids, # A character vector of your genes of interest
      pvalueCutoff = 1, # Can choose a FDR cutoff
      pAdjustMethod = "BH", # Method to be used for multiple testing correction
      universe = ensembl_hgene_ids, # A character vector containing your background set genes
      TERM2GENE = dplyr::select( # Pathway information: data frame with a term name or identifier and the gene identifiers
        gene_sets_df,
        gene_set_name,
        ensembl_id
      )
    )
    # Store results in a list
    fora_results_nonsyn_list[[sample_id]] <- fora_results_nonsyn
    
    # Add sample_id
    temp_fora_results_nonsyn_df <- as.data.frame(fora_results_nonsyn@result)
    temp_fora_results_nonsyn_df$sample_id <- sample_id
    fora_results_nonsyn_df <- rbind(fora_results_nonsyn_df, temp_fora_results_nonsyn_df)
    
    # Perform FORA for LoF mutations
    mutated_genes_lof_ids <- as.character(lof_mut_list[[sample_id]]$ensembl_id)
    fora_results_lof <- enricher(
      gene = mutated_genes_lof_ids, # A character vector of your genes of interest
      pvalueCutoff = 1, # Can choose a FDR cutoff
      pAdjustMethod = "BH", # Method to be used for multiple testing correction
      universe = ensembl_hgene_ids, # A character vector containing your background set genes
      TERM2GENE = dplyr::select( # Pathway information: data frame with a term name or identifier and the gene identifiers
        gene_sets_df,
        gene_set_name,
        ensembl_id
      )
    )
    # Store results in a list
    fora_results_lof_list[[sample_id]] <- fora_results_lof
    
    # Add sample_id
    temp_fora_results_lof_df <- as.data.frame(fora_results_lof@result)
    temp_fora_results_lof_df$sample_id <- sample_id
    fora_results_lof_df <- rbind(fora_results_lof_df, temp_fora_results_lof_df)
  }
  
  # Save the results to CSV files using the provided variable names
  write.csv(fora_results_all_df, file = fora_results_all_csv_path, row.names = TRUE)
  write.csv(fora_results_nonsyn_df, file = fora_results_nosnyn_csv_path, row.names = TRUE)
  write.csv(fora_results_lof_df, file = fora_results_lof_csv_path, row.names = TRUE)
  
  # Return
  results_lists <- list(
    all = fora_results_all_list,
    nonsyn = fora_results_nonsyn_list,
    lof = fora_results_lof_list
  )
  return(results_lists)
}

# Bubble Plot top 5 Gene Sets across samples stratified by response
plot_top5 <- function(response_df, mut_df, plot_path, fig_title, fig_y_label) {

  # Prepare data for all samples in the df: Filter padjust < 0.05 & calculate Gene ratio per sample
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

# Bubble Plot all Gene Sets across samples with pvalue significance
bubble_plot_pval <- function(response_df, mut_df, plot_path, fig_title, fig_y_label) {
  
  # Prepare data for all samples in the df: Filter padjust < 0.05 & calculate Gene ratio per sample
  prepared_data_df <- mut_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::mutate(Gene_ratio = overlap / size,
                  Count = overlap) %>%
    dplyr::arrange(sample_id, desc(Count)) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(short_sample_id = sub(".*_(\\d+)$", "\\1", sample_id))  # Extract the numeric part from sample_id
  
  # Complete the data to ensure all pathways are present
  complete_data_df <- prepared_data_df %>%
    tidyr::complete(pathway = unique(mut_df$pathway), fill = list(Gene_ratio = NA, Count = NA, padj = NA))
  
  # Calculate dynamic breaks and labels for Count
  count_breaks <- quantile(prepared_data_df$Count, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  count_breaks <- round(count_breaks) # Ensure count_breaks are integers by rounding
  count_breaks <- count_breaks[count_breaks > 0] # Remove 0 values from count_breaks
  count_breaks <- unique(count_breaks) # Make sure breaks are unique
  
  # Dynamically generate count_labels based on the available count_breaks
  count_labels <- vector("character", length(count_breaks))
  for (i in seq_along(count_breaks)) {
    count_labels[i] <- as.character(count_breaks[i])
  }
  
  # Step 2: Plot
  bubble_plot <- ggplot(complete_data_df, aes(x = Gene_ratio,
                                              y = reorder(pathway, Count),
                                              size = Count,
                                              fill = padj)) +
    geom_point(shape = 21, stroke = 0.5, color = "black") +  # Add black outline to the dots
    geom_label_repel(aes(label = short_sample_id), size = 3, color = "black", 
                     box.padding = 0.8, point.padding = 0.3, segment.color = "black", 
                     max.overlaps = 50, max.time = 2, alpha = 0.5, fill = "white") +  # Add sample_id labels with repel and edge
    scale_fill_gradient(low = "#D55E00", high = "#F0E442") +  # Fill gradient
    scale_size_continuous(name = "Count",
                          range = c(1, 6),  # Adjust the visual size range as needed
                          breaks = count_breaks,
                          labels = count_labels) +
    scale_x_continuous(breaks = seq(0, 0.08, by = 0.01)) +
    xlim(0,0.08) +
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
         x = "Gene Ratio",
         y = fig_y_label,
         fill = "Adjusted p-value",
         size = "Gene count")
  
  # Save plot
  ggsave(plot_path, bubble_plot, width = 10, height = 4, dpi = 300)
  
  # Return
  return(bubble_plot)
}

# Bubble Plot all Gene Sets across samples stratified by response (no pval)
bubble_plot_all <- function(response_df, mut_df, plot_path, fig_title, fig_y_label) {
  
  # Prepare data for all samples in the df: calculate Gene ratio per sample, do not filter by p.adjust and keep all GS
  prepared_data_df <- mut_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(Gene_ratio = overlap / size,
                  Count = overlap) %>%
    dplyr::arrange(sample_id, desc(Count)) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(Count > 0)  # Exclude rows where Count is 0
  
  # Merge response info
  combined_data <- prepared_data_df %>%
    dplyr::left_join(response_df, by = "sample_id")
  
  # Step 1: Create an ordered factor for sample_id based on patient_response
  response_order <- combined_data %>%
    dplyr::select(sample_id, patient_response) %>%
    dplyr::distinct() %>%
    dplyr::arrange(patient_response, sample_id) # Arrange so that "R" comes before "NR"
  
  combined_data$sample_id <- factor(combined_data$sample_id,
                                    levels = response_order$sample_id)
  
  # Calculate dynamic breaks and labels for Count
  count_breaks <- quantile(combined_data$Count, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  count_breaks <- round(count_breaks)
  count_breaks <- count_breaks[count_breaks > 0]
  count_breaks <- unique(count_breaks)
  
  count_labels <- vector("character", length(count_breaks))
  for (i in seq_along(count_breaks)) {
    count_labels[i] <- as.character(count_breaks[i])
  }
  
  # Step 2: Plot with ordered sample_id and facets to distinguish "R" and "NR"
  bubble_plot <- ggplot(combined_data, aes(x = sample_id,
                                           y = reorder(pathway, Count),
                                           size = Count,
                                           color = patient_response)) +
    geom_point() +
    scale_color_manual(values = c("R" = "darkgreen", "NR" = "darkred"), 
                       name = "Response",
                       labels = c("R" = "Responders", "NR" = "Non-responders")) +
    scale_size_continuous(name = "Count",
                          range = c(1, 6),
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
         color = "Response",
         size = "Count")
  
  # Save plot
  ggsave(plot_path, bubble_plot, width = 15, height = 20, dpi = 300)
  
  # Return
  return(bubble_plot)
}

# Bubble Plot for (a) specific Gene Set(s) across samples stratified by response
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

# Bar Plot Combined Gene Sets for 1 sample
plot_combinedGS <- function(sample_id, response_df, mut_df, plot_path, fig_title, max_x, dynamic_width) {
  
  # Prepare data for all samples in the df: Filter padjust < 0.05 & calculate Gene ratio per sample
  prepared_data_df <- mut_df %>%
    dplyr::filter(sample_id == !!sample_id) %>%
    # dplyr::filter(padj < 0.05) %>%
    dplyr::mutate(Gene_ratio = overlap / size,
                  Count = overlap) %>%
    dplyr::filter(Count > 0) %>%
    dplyr::arrange(desc(Count))
  
  # Merge response info
  combined_data <- prepared_data_df %>%
    dplyr::left_join(response_df, by = "sample_id")
  
  # Find the maximum value in the 'Count' column
  max_value <- max(combined_data$Count, na.rm = TRUE)
  
  # Define base plot width and scale factor
  base_width <- 10 # Base width for plots
  scale_factor <- 0.2 # Width increase per unit in max_value beyond a threshold
  
  # Calculate dynamic width based on max_value
  # dynamic_width <- base_width + (max_value * scale_factor)
  
  combinedGS_plot <- ggplot(combined_data, aes(x = Count,
                                               y = reorder(pathway, Count),
                                               fill = padj)) +
    geom_col() +
    scale_fill_gradient(low = "blue", high = "red", trans = "log10", # Log10 transformation for a better visualization
                        limits = c(1e-6, 1), oob = scales::oob_squish) +
    # scale_x_continuous(breaks = seq(0, max_value, by = 1)) +
    scale_x_continuous(
      breaks = seq(0, max_x, by = 5),
      limits = c(NA, max_x) # Set upper limit, lower limit is NA for auto-calculation
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
    ) +
    labs(title = fig_title,
         x = "Gene Count",
         y = "Enriched Gene Sets",
         fill = "p.adjust (log10)")
  
  # Save plot
  ggsave(plot_path, combinedGS_plot, width = dynamic_width, height = 15, dpi = 300)
  
  # Return
  return(combinedGS_plot)
}

# Box Plot for (a) specific Gene Set(s) across samples stratified by response
boxplot_genesets <- function(geneset, response_df, mut_df, plot_path, fig_title, max_y) {
  # geneset is a vector i.e. c("GS1", "GS2")
  
  # Prepare data: Filter Gene Set & calculate Gene ratio per sample
  prepared_data_df <- mut_df %>%
    dplyr::filter(pathway %in% geneset) %>%
    dplyr::mutate(Gene_ratio = overlap / size, Count = overlap) %>%
    dplyr::ungroup()
  
  # Merge response info
  combined_data <- prepared_data_df %>%
    dplyr::left_join(response_df, by = "sample_id")
  
  # Plot
  geneset_boxplot <- ggplot(combined_data, aes(x = pathway, 
                                               y = Count, 
                                               fill = patient_response)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +  # Adjust dodge width as necessary
    scale_fill_manual(values = c("NR" = "#E41A1C", "R" = "#377EB8"), 
                      labels = c("NR" = "Non-Responders", "R" = "Responders")) +
    scale_y_continuous(
      breaks = seq(0, max_y, by = 10),
      limits = c(NA, max_y) # Set upper limit, lower limit is NA for auto-calculation
    ) +
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
         x = "Gene Set",
         y = "Count",
         fill = "Patient Response")
  
  # Save plot
  ggsave(plot_path, geneset_boxplot, width = 20, height = 12, dpi = 300)
  
  # Return the plot object
  return(geneset_boxplot)
}

# Heat Map of combined GS across samples using ComplexHeatmap
complex_heatmap_genesets <- function(response_df, mut_df, plot_path, fig_title, width, height) {
  
  # Join response info and prepare data
  prepared_data_df <- mut_df %>%
    dplyr::left_join(response_df, by = "sample_id") %>%
    dplyr::mutate(Gene_ratio = overlap / size, Count = overlap)
  
  # Convert the data frame to a wide format and then to a matrix
  wide_data <- reshape2::dcast(prepared_data_df, pathway ~ sample_id, value.var = "Count", fill = 0)
  rownames(wide_data) <- wide_data$pathway
  wide_data$pathway <- NULL
  heat_matrix <- data.matrix(wide_data)
  
  # Make sure response_df is ordered to match the sample_id order in heat_matrix
  sample_responses <- response_df %>%
    dplyr::filter(sample_id %in% colnames(heat_matrix)) %>%
    dplyr::select(sample_id, patient_response) %>%
    dplyr::distinct() %>%
    dplyr::arrange(factor(sample_id, levels = colnames(heat_matrix)))
  
  # Create annotations for samples based on response
  response_annotation <- HeatmapAnnotation(df = data.frame(Response = factor(sample_responses$patient_response, levels = c("R", "NR"))),
                                           col = list(Response = c("R" = "darkgreen", "NR" = "darkred")),
                                           which = "col", 
                                           annotation_name_side = "left")
  
  # Draw the heatmap
  pdf(file = plot_path, width = width, height = height)
  heat_map_plot <- Heatmap(heat_matrix, name = "Count", 
                           bottom_annotation = response_annotation,
                           column_dend_side = "bottom",
                           row_names_side = "left",
                           show_row_names = TRUE, 
                           show_column_names = TRUE, 
                           column_title = gt_render(
                             paste0(fig_title, "<br>",
                                    "Samples")),
                           column_title_side = "top",
                           # row_title = "Ubiquitin Gene Sets", # Can't avoid overlapping
                           row_title_side = "left",  # Place row title on the left side
                           col = colorRampPalette(c("white", "blue"))(100),
                           column_names_side = "bottom",
                           show_row_dend = FALSE,
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "black", lwd = 0.5))
                           })
  
  # Print the heatmap and save it
  draw(heat_map_plot, heatmap_legend_side = 'right', annotation_legend_side = 'right', padding = unit(c(10, 70, 10, 10), "mm"))
  
  # Decorate the annotation to add vertical lines
  decorate_annotation("Response", {
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1))  # Draw a border around the annotation
    for(i in seq_along(sample_responses$sample_id)[-1]) {
      grid.lines(x = unit(i+11.5, "native") - unit(0.5, "npc"), y = unit(c(0, 1), "npc"),
                 gp = gpar(col = "black", lwd = 1))  # Thicker lines
    }
  })
  
  dev.off()
}

# Heat Map of cellular location across samples using ComplexHeatmap
complex_heatmap_cell_loc <- function(response_df, mut_df, plot_path, fig_title, width, height) {
  
  # Split the cell_loc by '|', and unnest to long format
  count_cell_loc_df <- mut_df %>%
    dplyr::mutate(cell_loc = strsplit(cell_loc, "\\|")) %>%
    tidyr::unnest(cell_loc) %>%
    dplyr::group_by(sample_id, cell_loc) %>%
    dplyr::summarise(count = n(), .groups = 'drop')
  
  # Convert the data frame to a wide format and then to a matrix
  wide_data <- reshape2::dcast(count_cell_loc_df, cell_loc ~ sample_id, value.var = "count", fill = 0)
  rownames(wide_data) <- wide_data$cell_loc
  wide_data$cell_loc <- NULL
  heat_matrix <- data.matrix(wide_data)
  
  # Make sure response_df is ordered to match the sample_id order in heat_matrix
  sample_responses <- response_df %>%
    dplyr::filter(sample_id %in% colnames(heat_matrix),
                  patient_response %in% mut_df$patient_response) %>%
    dplyr::select(sample_id, patient_response) %>%
    dplyr::distinct() %>%
    dplyr::arrange(factor(sample_id, levels = colnames(heat_matrix)))
  
  # Create annotations for samples based on response
  response_annotation <- HeatmapAnnotation(df = data.frame(Response = factor(sample_responses$patient_response, levels = c("R", "NR"))),
                                           col = list(Response = c("R" = "darkgreen", "NR" = "darkred")),
                                           which = "col", 
                                           annotation_name_side = "left")
  
  # Draw the heatmap
  png(file = plot_path, width = width, height = height)
  heat_map_plot <- Heatmap(heat_matrix, 
                           name = "Count", 
                           bottom_annotation = response_annotation,
                           column_dend_side = "bottom",
                           row_names_side = "left",
                           show_row_names = TRUE, 
                           show_column_names = TRUE, 
                           column_title = gt_render(
                             paste0(fig_title, "<br>",
                                    "Sample ID")),
                           column_title_side = "top",
                           row_title = "Cellular Location",
                           row_title_side = "left",  # Place row title on the left side
                           col = colorRampPalette(c("white", "blue"))(100),
                           column_names_side = "bottom",
                           show_row_dend = FALSE,
                           cell_fun = function(j, i, x, y, width, height, fill) {
                             grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "black", lwd = 0.5))
                           })
  
  # Print the heatmap and save it
  draw(heat_map_plot, heatmap_legend_side = 'right', annotation_legend_side = 'right', padding = unit(c(10, 70, 10, 10), "mm"))
  
  # Decorate the annotation to add vertical lines
  decorate_annotation("Response", {
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 1))  # Draw a border around the annotation
    for(i in seq_along(sample_responses$sample_id)[-1]) {
      grid.lines(x = unit(i+10, "native") - unit(0.5, "npc"), y = unit(c(0, 1), "npc"),
                 gp = gpar(col = "black", lwd = 1))  # Thicker lines
    }
  })
  
  dev.off()
}

# Heat Map of combined GS across samples
heat_map_genesets <- function(response_df, mut_df, plot_path, fig_title) {
  
  # Assuming response_df has 'sample_id' and 'patient_response' where patient_response are "R" or "NR"
  
  # Ensure that patient_response is a factor with levels in the correct order for sorting
  response_df$patient_response <- factor(response_df$patient_response, levels = c("R", "NR"))
  
  # Prepare data for all samples
  prepared_data_df <- mut_df %>%
    dplyr::mutate(Gene_ratio = overlap / size, Count = overlap) %>% 
    dplyr::select(sample_id, pathway, Count) %>%
    dplyr::left_join(response_df, by = "sample_id") %>%
    dplyr::mutate(sample_response = paste(sample_id, patient_response, sep = " - ")) %>%
    # Create an ordering column based on patient_response to sort "R" first then "NR"
    dplyr::mutate(ordering = as.numeric(patient_response))
  
  # Plot with ordered sample_response
  heat_map_plot <- ggplot2::ggplot(prepared_data_df, ggplot2::aes(x = reorder(sample_response, ordering), 
                                                                  y = pathway, 
                                                                  fill = Count)) +
    ggplot2::geom_tile(color = "black", linewidth = 0.2) +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = fig_title, y = "Gene Set", x = "Sample ID", fill = "Count") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Save plot
  ggplot2::ggsave(plot_path, heat_map_plot, width = 15, height = 10, dpi = 300)
  
  # Return
  return(heat_map_plot)
}

# Box Plot for (a) specific Gene Set(s) across samples stratified by response
boxplot_genesets <- function(geneset, response_df, mut_df, plot_path, fig_title, max_y) {
  # geneset is a vector i.e. c("GS1", "GS2")
  
  # Prepare data: Filter Gene Set & calculate Gene ratio per sample
  prepared_data_df <- mut_df %>%
    dplyr::filter(pathway %in% geneset) %>%
    dplyr::mutate(Gene_ratio = overlap / size, Count = overlap) %>%
    dplyr::ungroup()
  
  # Merge response info
  combined_data <- prepared_data_df %>%
    dplyr::left_join(response_df, by = "sample_id")
  
  # Plot
  geneset_boxplot <- ggplot(combined_data, aes(x = pathway, 
                                               y = Count, 
                                               fill = patient_response)) +
    geom_boxplot(position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = c("NR" = "#E41A1C", "R" = "#377EB8"), 
                      labels = c("NR" = "Non-Responders", "R" = "Responders")) +
    scale_y_continuous(
      breaks = seq(0, max_y, by = 10),
      limits = c(NA, max_y) # Set upper limit, lower limit is NA for auto-calculation
    ) +
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
         x = "Gene Set",
         y = "Count",
         fill = "Patient Response")
  
  # Save plot
  ggsave(plot_path, geneset_boxplot, width = 20, height = 12, dpi = 300)
  
  # Return the plot object
  return(geneset_boxplot)
}

# Heat Map of top5 GS across samples stratified by response
top5_heat_map <- function(response_df, mut_df, plot_path, fig_title) {
  
  # Prepare data for all samples
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
  
  # Find the total number of pathways
  max_value <- length(unique(combined_data$pathway))
  
  # Define base plot width and scale factor
  base_width <- 10 # Base width for plots
  scale_factor <- 0.01 # Width increase per unit in max_value beyond a threshold
  
  # Calculate dynamic width based on max_value
  dynamic_width <- base_width + (max_value * scale_factor)
  
  # Step 2: Plot with ordered sample_id and facets to distinguish "R" and "NR"
  heat_map_plot <- ggplot(combined_data, aes(x = pathway, 
                                             y = sample_id, 
                                             fill = Count)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradient(low = "yellow", high = "red") + # Customize gradient colors as needed
    theme_minimal() +
    labs(fill = "Count", x = "Pathway", y = "Sample ID") +
    theme(
      panel.grid.major = element_line(color = "grey", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
    labs(title = fig_title,
         x = "Gene Set",
         y = "Sample ID",
         fill = "Count")

  # Save plot
  ggsave(plot_path, heat_map_plot, width = dynamic_width, height = 10, dpi = 300, limitsize = FALSE)
  
  # Return
  return(heat_map_plot)
}

# Function to extract PFAM info
extract_family_info <- function(file_path) {
  # Load and parse the JSON data from file
  json_data <- fromJSON(file_path, flatten = TRUE)
  
  # Initialize an empty data frame for collecting FAMILY type entries
  family_entries_df <- data.frame()
  
  # Check if results or matches exist and is not empty
  if (!is.null(json_data$results) && "matches" %in% names(json_data$results)) {
    # Get matches data
    matches <- json_data$results$matches
    
    for (i in seq_along(matches)) {
      p_id <- json_data$results$xref[[i]]$id
      
      # Error handling
      if (length(p_id) > 1) {
        p_id <- p_id[1]
      }
      
      # Add p_id
      matches[[i]] <- matches[[i]] %>% 
        dplyr::mutate(Protein_accession = p_id) %>% 
        dplyr::select(Protein_accession, everything())
      
      # Append this match to the main data frame
      family_entries_df <- bind_rows(family_entries_df, matches[[i]])
    }
  }
  
  return(family_entries_df)
}

# Tile plot of PFAM domains using ComplexHeatmap
tile_plot_pfam <- function(mut_df, plot_path, fig_title, width, height) {
  
  # Convert to factors
  mut_df$corrected_hugo_symbol <- as.factor(mut_df$corrected_hugo_symbol)
  mut_df$dom_name <- as.factor(mut_df$dom_name)
  mut_df$source <- as.factor(mut_df$source)
  
  # Define the color mapping for 'source'
  color_mapping <- c("R" = "#009E73", "NR" = "#D55E00", "Shared" = "#8A2BE2")
  
  # Create the plot
  tile_plot <- ggplot(mut_df, aes(x = dom_name, y = corrected_hugo_symbol, fill = source)) +
    geom_tile(color = "white") +  # Add white borders for clarity
    scale_fill_manual(values = color_mapping) +  # Use custom colors
    labs(x = "Domain Name", y = "HUGO Symbol", fill = "Source", title = fig_title) +
    theme_minimal() +  # Clean minimalistic theme
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x labels for better readability
          axis.title = element_text(size = 12, face = "bold"), # Bold axis titles
          legend.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
  
  # Save plot
  ggsave(plot_path, tile_plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  
  # Return
  return(tile_plot)
}

# Define the function to create a boxplot for a given biomarker
create_biomarker_boxplot <- function(df, biomarkers, response_col = "Response", fill_colors = c("R" = "#009E73", "NR" = "#D55E00")) {
  # Check if df is a data frame
  if (!is.data.frame(df)) {
    stop("The input data (df) should be a data frame.")
  }

  # Check if biomarkers is a character vector
  if (!is.character(biomarkers)) {
    stop("The biomarkers argument should be a character vector.")
  }

  # Check if all biomarkers exist in the data frame
  missing_biomarkers <- setdiff(biomarkers, colnames(df))
  if (length(missing_biomarkers) > 0) {
    stop("The following biomarkers are not found in the data frame: ", paste(missing_biomarkers, collapse = ", "))
  }

  # Normalize the biomarker values
  df[biomarkers] <- lapply(df[biomarkers], function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Pivot the data longer for all biomarkers
  long_data <- df %>%
    pivot_longer(cols = all_of(biomarkers), 
                 names_to = "Measurement", values_to = "Value")
  
  # Create the boxplot
  p <- ggplot(long_data, aes(x = Measurement, y = Value, fill = !!sym(response_col))) +
    geom_boxplot() +
    scale_fill_manual(values = fill_colors) +
    labs(x = "Measurement", y = "Normalized Value", title = "Boxplots for Biomarkers") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
  
  # Add p-values to the plot
  p + stat_compare_means(aes(group = !!sym(response_col)), method = "wilcox.test", label = "p.format")
}

# Define the function to create a boxplot for proteasome biomarkers
create_proteas_boxplot <- function(df, biomarkers, response_col = "Response", fill_colors = c("R" = "#009E73", "NR" = "#D55E00")) {
  # Check if df is a data frame
  if (!is.data.frame(df)) {
    stop("The input data (df) should be a data frame.")
  }
  
  # Check if biomarkers is a character vector
  if (!is.character(biomarkers)) {
    stop("The biomarkers argument should be a character vector.")
  }
  
  # Check if all biomarkers exist in the data frame
  missing_biomarkers <- setdiff(biomarkers, colnames(df))
  if (length(missing_biomarkers) > 0) {
    stop("The following biomarkers are not found in the data frame: ", paste(missing_biomarkers, collapse = ", "))
  }
  
  # Normalize the biomarker values
  df[biomarkers] <- lapply(df[biomarkers], function(x) (x - min(x)) / (max(x) - min(x)))
  
  # Pivot the data longer for all biomarkers
  long_data <- df %>%
    pivot_longer(cols = all_of(biomarkers), 
                 names_to = "Measurement", values_to = "Value")
  
  # Set the order of the Measurement factor
  long_data$Measurement <- factor(long_data$Measurement, levels = biomarkers)
  
  # Calculate medians for annotation
  medians <- long_data %>%
    group_by(Measurement, !!sym(response_col)) %>%
    summarise(median_value = median(Value), .groups = 'drop')
  
  # Create the boxplot
  p <- ggplot(long_data, aes(x = Measurement, y = Value, fill = !!sym(response_col))) +
    geom_boxplot() +
    geom_text(data = medians, aes(x = Measurement, y = median_value, label = round(median_value, 2)),
              position = position_dodge(width = 0.75), size = 4, vjust = -0.5, color = "white") +
    scale_fill_manual(values = c(R = "#009E73", NR = "#D55E00"),
                      labels = c("R" = "Responders", "NR" = "Non-responders")) +
    labs(x = "Biomarker", y = "Expression level (normalized)", title = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          plot.background = element_rect(fill = "white", colour = NA))
  
  # Add p-values to the plot
  p + stat_compare_means(aes(group = !!sym(response_col)), method = "wilcox.test", label = "p.format", vjust = -1)
}

# Define biomarkers clustering order
biomarkers_clustering_order <- function(df, biomarkers, measure = "median", response_col = "Response", sample_col = "SampleID") {
  # Ensure the response column is a factor
  df <- df %>%
    mutate(!!sym(response_col) := as.factor(!!sym(response_col)),
           !!sym(sample_col) := as.factor(!!sym(sample_col)))
  
  # Calculate the median and categorize each biomarker
  for (biomarker in biomarkers) {
    category_col <- paste0(biomarker, "_Category")
    if (measure == "median") {
      statistic_value <- median(df[[biomarker]], na.rm = TRUE)
    } else {
      statistic_value <- mean(df[[biomarker]], na.rm = TRUE)
    }
    df <- df %>%
      dplyr::mutate(!!sym(category_col) := ifelse(!!sym(biomarker) >= statistic_value, "high", "low"))
  }
  
  # Convert the categorical values to numerical for clustering
  df_for_clustering <- df %>%
    mutate(across(ends_with("_Category"), ~ as.numeric(. == "high")))
  
  # Combine the categorical numeric columns into a matrix for clustering
  binary_matrix <- as.matrix(df_for_clustering %>%
                               dplyr::select(ends_with("_Category")))
  
  # Perform hierarchical clustering on the combined binary matrix
  dist_matrix <- dist(binary_matrix)
  clustering <- hclust(dist_matrix)
  
  # Get the order of SampleID based on clustering
  sample_order <- df[[sample_col]][clustering$order]
  
return(list("df_biomarkers_category" = df, "sample_order" = sample_order, "clustering" = clustering))
}

# Function to make ubiquitination and deubiquitination categories plot (patchwork)
make_ubi_annotation_plot <- function(df, var_name, fill_label, custom_labels) {
  # Ensure the var_name column is a factor with levels 0 and 1
  df[[var_name]] <- factor(df[[var_name]], levels = c(0, 1), labels = custom_labels)
  
  plot <- ggplot(df, aes_string(x = "SampleID", y = "1", fill = var_name)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_fill_manual(values = c("Present" = "darkgrey", "Absent" = "white"), na.value = "white") +
    labs(x = NULL, y = NULL, fill = fill_label) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right",
          plot.margin = margin(0, 0, 0, 0, "pt"))
  
  return(plot)
}

# Function to make ecotype and tme categories plot (patchwork)
make_ecotyper_tme_annotation_plot <- function(df, var_name, fill_label, col_start, custom_labels) {
  plot <- ggplot(df, aes_string(x = "sample_id", y = "1", fill = var_name)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_fill_manual(values=cbPalette[col_start:length(cbPalette)], labels = custom_labels, na.value = "white") +
    labs(x = NULL, y = NULL, fill = fill_label) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right",
          plot.margin = margin(0, 0, 0, 0, "pt"))
  
  return(plot)
}

# Function to make the biomarkers plot (patchwork)
make_biomarker_annotation_plot <- function(df, var_name, fill_label, show_sample_id = FALSE) {
  plot <- ggplot(df, aes_string(x = "sample_id", y = "1", fill = var_name)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    scale_fill_manual(values = c("high" = "#009E73", "low" = "#D55E00"),
                      labels = c("high" = "High", "low" = "Low")) +
    labs(x = NULL, y = NULL, fill = "Biomarker score") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = ifelse(show_sample_id == TRUE, "right", "none"),
          plot.margin = margin(0, 0, 0, 0, "pt")) +
  
  # Conditionally add sample ID labels
  if (show_sample_id) {
    plot <- plot + 
      geom_text(aes(label = "sample_id"), y = 0, vjust = 0, size = 3) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.ticks.x = element_blank())
  }
  
  return(plot)
}

# Define function to generate combined plot for all biomarkers
make_combined_biomarker_plots <- function(df, biomarkers, ecotype_plot, tme_plot, ubi_plot, deubi_plot, dend_data) {
  
  # Create dendrogram
  dend_plot <- ggplot() +
    geom_segment(data = dend_data$segments, aes(x, y, xend = xend, yend = yend)) +
    scale_y_reverse() +
    coord_fixed(ratio = 1.45, expand = FALSE) + 
    labs(x = NULL, 
         y = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 12.2, 0, 0), "cm"))
  
  # Fix the width by using cowplot's plot_grid
  dend_plot <- plot_grid(dend_plot, ncol = 1, rel_widths = c(0.1))
  
  # Create the response plot
  response_plot <- ggplot(df, aes(x = sample_id, y = "1", fill = patient_response)) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("R" = "#0072B2", "NR" = "#F0E442"),
                      labels = c("R" = "Responders", "NR" = "Non-responders")) +
    labs(x = NULL, 
         y = NULL,
         fill = "Patient Response",
         title = "Cancer-Immunity Cycle biomarkers") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right",
          plot.margin = margin(0, 0, 0, 0, "pt"))
  
  # Create first label plot
  manual_label_plots <- list()
  manual_label_plots[[1]] <- ggdraw() + 
    draw_label("Patient response", x = 0, hjust = 0, angle = 0, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  manual_label_plots[[2]] <- ggdraw() + 
    draw_label("Ecotype category", x = 0, hjust = 0, angle = 0, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  manual_label_plots[[3]] <- ggdraw() + 
    draw_label("TME category", x = 0, hjust = 0, angle = 0, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  manual_label_plots[[4]] <- ggdraw() + 
    draw_label("Ubi mutation", x = 0, hjust = 0, angle = 0, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  manual_label_plots[[5]] <- ggdraw() + 
    draw_label("Deubi mutation", x = 0, hjust = 0, angle = 0, size = 12) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  
  # Create plots for each biomarker
  biomarker_plots <- list()
  label_plots <- list()
  for (i in seq_along(biomarkers)) {
    show_sample_id <- ifelse(i == length(biomarkers), TRUE, FALSE)
    biomarker_plots[[i]] <- make_biomarker_annotation_plot(
      df, paste0(biomarkers[i], "_Category"), biomarkers[i], show_sample_id)
    label_plots[[i]] <- ggdraw() + 
      draw_label(biomarkers[i], x = 0, hjust = 0, angle = 0, size = 12) +
      theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  }
  
  # Extract legends
  legends <- list(
    p1_legend = get_legend(response_plot),
    p2_legend = get_legend(ecotype_plot),
    p3_legend = get_legend(tme_plot),
    p4_legend = get_legend(ubi_plot),
    p5_legend = get_legend(deubi_plot),
    p6_legend = get_legend(biomarker_plots[[length(biomarkers)]])
  )
  
  # Combine all plots
  combined_plot <- response_plot + theme(legend.position = "none") + manual_label_plots[[1]]
  combined_plot <- combined_plot / (ecotype_plot + theme(legend.position = "none") + manual_label_plots[[2]])
  combined_plot <- combined_plot / (tme_plot + theme(legend.position = "none") + manual_label_plots[[3]])
  combined_plot <- combined_plot / (ubi_plot + theme(legend.position = "none") + manual_label_plots[[4]])
  combined_plot <- combined_plot / (deubi_plot + theme(legend.position = "none") + manual_label_plots[[5]])
  for (i in seq_along(biomarkers)) {
    combined_plot <- combined_plot / (biomarker_plots[[i]] + theme(legend.position = "none") + label_plots[[i]])
  }
  combined_plot <- combined_plot / (dend_plot + theme(legend.position = "none"))
  
  combined_legends <- (as_ggplot(legends$p1_legend) + theme(plot.margin = unit(c(50,50,0,20), "pt")) + coord_cartesian(clip = "off")) +
    (as_ggplot(legends$p3_legend) + theme(plot.margin = unit(c(50,10,0,50), "pt")) + coord_cartesian(clip = "off")) +
    (as_ggplot(legends$p2_legend) + theme(plot.margin = unit(c(70,20,0,80), "pt")) + coord_cartesian(clip = "off")) +
    (as_ggplot(legends$p4_legend) + theme(plot.margin = unit(c(50,10,10,80), "pt")) + coord_cartesian(clip = "off")) /
    (as_ggplot(legends$p5_legend) + theme(plot.margin = unit(c(50,10,0,80), "pt")) + coord_cartesian(clip = "off")) +
    (as_ggplot(legends$p6_legend) + theme(plot.margin = unit(c(0,50,0,100), "pt")) + coord_cartesian(clip = "off")) +
    plot_layout(nrow = 1, ncol = 6)
  
  final_plot <- combined_plot / (combined_legends + coord_cartesian(clip = "off")) +
    plot_layout(nrow = length(biomarkers) + 7, ncol = 1, heights = c(rep(0.3, length(biomarkers) + 5), 2))
  
  return(final_plot)
}

# Function to remove duplicate genes from the TPM matrix
removeDuplicateGenesDF <- function(geneMatrix){
  # From gene matrix with HGNC symbols
  # The next three lines are for when making tpm to dataframe, that also changes gene names dash to dots and additionally adds numbers to duplicated genes.
  genes <- rownames(geneMatrix)
  geneMatrix <- geneMatrix %>% as.data.frame() %>% rownames_to_column(var="Gene")
  if(!is.null(genes)){
    geneMatrix$Gene <- genes
    # Write gene expression matrix to tab-delimited file - TPM, AND annotation file
    # Gene names should be unique.
    geneMatrix.dedup <- IOBR::remove_duplicate_genes(eset = geneMatrix , column_of_symbol = "Gene", method = "mean")
    geneMatrix.dedup
  }else{
    stop("No gene ids in the matrix as rownames, please provide gene expression matrix in correct format")
  }
}

# Function to make a tile plot using pheatmap and patchwork for biomarkers, much cleaner and easier solution
pheatmap_biomarkers <- function(df, numerical_biomarkers, categorical_biomarkers, response_col, ubi_counts = FALSE, split_response = FALSE) {
  
  # Set row names
  rownames(df) <- df$SampleID
  df$SampleID <- NULL
  
  # Separate numerical and categorical data
  numeric_data <- df %>% dplyr::select(one_of(numerical_biomarkers))
  ann_df <- df %>% dplyr::select(!!sym(response_col), one_of(categorical_biomarkers))
  ann_df$SampleID <- rownames(ann_df)
  
  # Get ubi data
  if (ubi_counts) {
    ubi_data <- df %>% dplyr::select("ubi_counts", "deubi_counts")
  } else {
    ubi_data <- df %>% dplyr::select("ubi_lof", "deubi_lof")
  }
  
  # Merge
  ann_df <- bind_cols(ann_df, ubi_data)
  
  # Normalize numerical data
  normalized_data <- as.data.frame(lapply(numeric_data, function(x) {
    if (all(is.na(x))) {
      return(rep(NA, length(x)))
    } else if (min(x, na.rm = TRUE) == max(x, na.rm = TRUE)) {
      return(rep(0.5, length(x)))  # Arbitrary constant for no variance columns
    } else {
      return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
    }
  }))
  normalized_data$SampleID <- rownames(numeric_data)
  
  # Melt the normalized data for ggplot
  melted_data <- melt(normalized_data, id.vars = "SampleID")
  
  # Convert to factors ann_df
  ann_df <- ann_df %>% 
    mutate(TME = as.factor(TME),
           Ecotype = factor(Ecotype, levels = c(paste0("CE", 1:10), NA)),
           !!sym(response_col) := as.factor(!!sym(response_col)))
  
  if (response_col == "resp_3_cat") {
    ann_df <- ann_df %>% 
      dplyr::rename(Response = resp_3_cat) %>% 
      dplyr::mutate(Response = recode(Response, "R" = "Responders", "SD" = "Stable Disease", "NR" = "Non-responders"))
  } else {
    ann_df <- ann_df %>% 
      dplyr::mutate(Response = recode(Response, "R" = "Responders", "NR" = "Non-responders"))
  }
  
  if (!ubi_counts) {
    ann_df <- ann_df %>% 
      dplyr::rename(Ubiquitin_LoF = ubi_lof,
                    Deubiquitin_LoF = deubi_lof) %>% 
      dplyr::mutate(Ubiquitin_LoF = as.factor(ifelse(Ubiquitin_LoF == 1, "Present", "Absent")),
                    Deubiquitin_LoF = as.factor(ifelse(Deubiquitin_LoF == 1, "Present", "Absent")))
  } else {
    ann_df <- ann_df %>% 
      dplyr::rename(Ubiquitin_LoF_counts = ubi_counts,
                    Deubiquitin_LoF_counts = deubi_counts)
  }
  
  # Add annotation bars
  # For this example, we will use pheatmap to simplify annotation integration
  # Define the initial categorical colors list
  ann_colors <- list(
    TME = c(D = "#FFD700", F = "#000000", IE = "#56B4E9", `IE/F` = "#0072B2"),
    Ecotype = c(CE1 = "#E69F00", CE2 = "#CC79A7", CE3 = "#F0E442", CE4 = "#FF69B4", CE5 = "#009E73", CE6 = "#8A2BE2", CE7 = "#D55E00", CE8 = "#000000", CE9 = "#00CED1", CE10 = "#4B0082", "NA" = "lightgrey")
  )
  
  # Conditionally add the Response colors
  if (response_col == "Response") {
    ann_colors$Response <- c("Responders" = "#009E73", "Non-responders" = "#D55E00")
  } else if (response_col == "resp_3_cat") {
    ann_colors$Response <- c("Responders" = "#009E73", "Stable Disease" = "#8A2BE2", "Non-responders" = "#D55E00")
  }
  
  # Conditionally add the ubiquitin info
  gradient_ubi <- colorRampPalette(c("white", "#E69F00"))(100)
  gradient_deubi <- colorRampPalette(c("white", "#CC79A7"))(100)
  if (!ubi_counts) {
    ann_colors$Ubiquitin_LoF <- c("Present" = "#E69F00", "Absent" = "white")
    ann_colors$Deubiquitin_LoF <- c("Present" = "#CC79A7", "Absent" = "white")
  } else {
    ann_colors$Ubiquitin_LoF_counts <- gradient_ubi
    ann_colors$Deubiquitin_LoF_counts <- gradient_deubi
  }
  
  # Plot heatmap with annotations
  df.fin <- as.data.frame(as.matrix(normalized_data[, -ncol(normalized_data)]))
  rownames(df.fin) <- ann_df$SampleID
  ann_df$SampleID <- NULL
  
  # Transpose
  df.fin.t <- t(df.fin)
  
  if (!split_response){
    p <- pheatmap::pheatmap(df.fin.t, 
                     cluster_rows = FALSE, 
                     cluster_cols = TRUE, 
                     annotation_col = ann_df, 
                     annotation_colors = ann_colors,
                     show_rownames = TRUE, 
                     show_colnames = TRUE,
                     scale = "none", 
                     border_color = "black",
                     angle_col = "45",
                     legend_breaks = c(0, 1),
                     legend_labels = c("Low", "High"),
                     display_numbers = T,
                     number_color = "black",
                     fontsize_number = 12,
                     color = colorRampPalette(c("#F0E442", "#0072B2"))(50)) %>% as.ggplot() +
      coord_cartesian(clip = "off")
    return(p)
  } else if (response_col == "Response") {
    p.r <- pheatmap::pheatmap(df.fin.t %>% as.data.frame() %>% dplyr::select(all_of(ann_df %>% dplyr::filter(Response=="Responders") %>% rownames())) %>% as.matrix(), 
                  cluster_rows = FALSE, 
                  cluster_cols = TRUE, 
                  annotation_col = ann_df %>% dplyr::filter(Response=="Responders"), 
                  annotation_colors = ann_colors,
                  annotation_legend = FALSE,
                  show_rownames = TRUE, 
                  show_colnames = TRUE,
                  scale = "none", 
                  border_color = "black",
                  angle_col = "45",
                  legend_breaks = c(0, 1),
                  legend_labels = c("Low", "High"),
                  display_numbers = T,
                  number_color = "black",
                  fontsize_number = 12,
                  color = colorRampPalette(c("#F0E442", "#0072B2"))(50),
                  legend = FALSE) %>% as.ggplot() +
      coord_cartesian(clip = "off")
    
    p.nr <- pheatmap::pheatmap(df.fin.t %>% as.data.frame()%>% dplyr::select(all_of(ann_df %>% dplyr::filter(Response=="Non-responders") %>% rownames())) %>% as.matrix(), 
                     cluster_rows = FALSE, 
                     cluster_cols = TRUE, 
                     annotation_col = ann_df %>% dplyr::filter(Response=="Non-responders"), 
                     annotation_colors = ann_colors,
                     show_rownames = TRUE, 
                     show_colnames = TRUE,
                     scale = "none", 
                     border_color = "black",
                     angle_col = "45",
                     legend_breaks = c(0, 1),
                     legend_labels = c("Low", "High"),
                     display_numbers = T,
                     number_color = "black",
                     fontsize_number = 12,
                     color = colorRampPalette(c("#F0E442", "#0072B2"))(50)) %>% as.ggplot() +
      coord_cartesian(clip = "off")
    combined_plot <- p.r + plot_spacer() + p.nr + plot_layout(widths = c(0.9, 0, 1.2))  # Adjust the widths as needed
    return(combined_plot)
  } else {
    p.r <- pheatmap::pheatmap(df.fin.t %>% as.data.frame() %>% dplyr::select(all_of(ann_df %>% dplyr::filter(Response=="Responders") %>% rownames())) %>% as.matrix(), 
                  cluster_rows = FALSE, 
                  cluster_cols = TRUE, 
                  annotation_col = ann_df %>% dplyr::filter(Response=="Responders"), 
                  annotation_colors = ann_colors,
                  annotation_legend = FALSE,
                  show_rownames = TRUE, 
                  show_colnames = TRUE,
                  scale = "none", 
                  border_color = "black",
                  angle_col = "45",
                  legend_breaks = c(0, 1),
                  legend_labels = c("Low", "High"),
                  display_numbers = T,
                  number_color = "black",
                  fontsize_number = 12,
                  color = colorRampPalette(c("#F0E442", "#0072B2"))(50),
                  legend = FALSE) %>% as.ggplot() +
      coord_cartesian(clip = "off")
    
    p.sd <- pheatmap::pheatmap(df.fin.t %>% as.data.frame() %>% dplyr::select(all_of(ann_df %>% dplyr::filter(Response=="Stable Disease") %>% rownames())) %>% as.matrix(), 
                  cluster_rows = FALSE, 
                  cluster_cols = TRUE, 
                  annotation_col = ann_df %>% dplyr::filter(Response=="Stable Disease"), 
                  annotation_colors = ann_colors,
                  annotation_legend = FALSE,
                  show_rownames = TRUE, 
                  show_colnames = TRUE,
                  scale = "none", 
                  border_color = "black",
                  angle_col = "45",
                  legend_breaks = c(0, 1),
                  legend_labels = c("Low", "High"),
                  display_numbers = T,
                  number_color = "black",
                  fontsize_number = 12,
                  color = colorRampPalette(c("#F0E442", "#0072B2"))(50),
                  legend = FALSE) %>% as.ggplot() +
      coord_cartesian(clip = "off")
    
    p.nr <- pheatmap::pheatmap(df.fin.t %>% as.data.frame()%>% dplyr::select(all_of(ann_df %>% dplyr::filter(Response=="Non-responders") %>% rownames())) %>% as.matrix(), 
                     cluster_rows = FALSE, 
                     cluster_cols = TRUE, 
                     annotation_col = ann_df %>% dplyr::filter(Response=="Non-responders"), 
                     annotation_colors = ann_colors,
                     show_rownames = TRUE, 
                     show_colnames = TRUE,
                     scale = "none", 
                     border_color = "black",
                     angle_col = "45",
                     legend_breaks = c(0, 1),
                     legend_labels = c("Low", "High"),
                     display_numbers = T,
                     number_color = "black",
                     fontsize_number = 12,
                     color = colorRampPalette(c("#F0E442", "#0072B2"))(50)) %>% as.ggplot() +
      coord_cartesian(clip = "off")
    combined_plot <- p.r + plot_spacer() + p.sd + plot_spacer() + p.nr + plot_layout(widths = c(0.9, 0, 0.9, 0, 1.2))  # Adjust the widths as needed
    return(combined_plot)
  }
}


barplot_mean_sd <- function(df, var_term, x_label, y_label, width, height, vjust, size, plot_path) {
  
  # Count occurrences of each var_term per sample_id
  counts_df <- df %>%
    group_by(sample_id, !!sym(var_term)) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Get the number of samples for each patient_response category
  sample_counts <- df %>%
    dplyr::select(sample_id, patient_response) %>%
    distinct() %>%
    group_by(patient_response) %>%
    summarise(sample_count = n(), .groups = 'drop')
  
  # Summarize by patient_response to calculate the mean and standard deviation
  summary_df <- counts_df %>%
    left_join(df %>% dplyr::select(sample_id, patient_response) %>% distinct(), by = "sample_id") %>%
    group_by(patient_response, !!sym(var_term)) %>%
    summarise(
      total_count = sum(count),
      std_dev = sd(count),
      .groups = 'drop'
    ) %>%
    left_join(sample_counts, by = "patient_response") %>%
    mutate(mean_count = total_count / sample_count)
  
  # Calculate p-values using Wilcoxon test
  wilcox_results <- counts_df %>%
    left_join(df %>% dplyr::select(sample_id, patient_response) %>% distinct(), by = "sample_id") %>%
    group_by(!!sym(var_term)) %>%
    summarise(
      p_value = if (n_distinct(patient_response) == 2) 
        wilcox.test(count ~ patient_response, exact = FALSE)$p.value 
      else NA,
      .groups = 'drop'
    ) %>%
    mutate(p_value = ifelse(is.nan(p_value), NA, p_value))
  
  # Combine summary statistics and p-values
  final_stats <- summary_df %>%
    left_join(wilcox_results, by = var_term) %>% 
    mutate(p_value_label = case_when(
      is.na(p_value) ~ NA_character_,
      TRUE ~ paste0("p=", format(p_value, digits = 2))
    ))
  
  # Calculate the maximum mean count to set the y-axis limit dynamically
  max_mean_count <- max(final_stats$mean_count, na.rm = TRUE)
  buffer <- max_mean_count * 0.1 # Add 10% buffer to the maximum mean count for the labels
  
  # Plotting
  plot <- ggplot(final_stats, aes(x = !!sym(var_term), y = mean_count, fill = patient_response)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.8) +
    geom_errorbar(aes(ymin = pmax(0, mean_count - std_dev), ymax = mean_count + std_dev), 
                  position = position_dodge(0.8), width = 0.25) +
    scale_fill_manual(values = c("R" = "#009E73", "NR" = "#D55E00"),
                      labels = c("R" = "Responders", "NR" = "Non-responders")) +
    labs(title = NULL,
         x = x_label,
         y = y_label,
         fill = "Patient Response") +
    ylim(0, max_mean_count + buffer) + # Adjust the y-axis to include the buffer
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  # Rotate X-axis labels for better readability
    scale_y_continuous(breaks = seq(0, max_mean_count + 25, by = 1),  # Ticks every unit
                       minor_breaks = seq(0, max_mean_count + 25, by = 1))  # Minor ticks every 1 unit
  
  # Add the p-value text labels
  plot <- plot + geom_text(data = final_stats, aes(x = !!sym(var_term), y = max_mean_count + buffer / 2, label = p_value_label), 
                           vjust = vjust, hjust = 0.5, size = size, color = "black")
  
  ggsave(plot_path, plot, width = width, height = height, dpi = 300)
  return(plot)
}