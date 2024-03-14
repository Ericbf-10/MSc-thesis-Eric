########################
### FUNCTIONS SCRIPT ###
########################

### Functions needed for the project

## Load libraries
library(msigdbr)
library(fgsea)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)

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

# Bar Plot that summarizes the mutation count for each patient
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

# Heat Map of combined GS across samples
heat_map_genesets <- function(response_df, mut_df, plot_path, fig_title) {
  
  # Prepare data for all samples
  prepared_data_df <- mut_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(Gene_ratio = overlap / size,
                  Count = overlap) %>%
    dplyr::arrange(sample_id, desc(Count)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(sample_id, pathway, Count)
  
  # Merge response info
  combined_data <- prepared_data_df %>%
    dplyr::left_join(response_df, by = "sample_id")
  
  # Step 2: Plot with ordered sample_id and facets to distinguish "R" and "NR"
  heat_map_plot <- ggplot(combined_data, aes(x = pathway, 
                                             y = sample_id, 
                                             fill = Count)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradient(low = "yellow", high = "red") + # Customize gradient colors as needed
    theme_minimal() +
    labs(fill = "Count", x = "Pathway", y = "Sample ID") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Improve axis label readability
    labs(title = fig_title,
         x = "Gene Set",
         y = "Sample ID",
         fill = "Count")
  
  # Save plot
  ggsave(plot_path, heat_map_plot, width = 15, height = 10, dpi = 300)
  
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