####################
### PATHS SCRIPT ###
####################

## Set all paths

# Find the root of the project
root <- rprojroot::find_root(rprojroot::is_rstudio_project)

# Define the current directory as the root directory
current_dir <- root

### Define paths
data_path <- file.path(current_dir,
                       "./data/_raw")
results_pc_path <- file.path(current_dir,
                             "./results/patient_characteristics/output")
figures_pc_path <- file.path(current_dir,
                             "./results/patient_characteristics/figures")
results_fip_path <- file.path(current_dir,
                              "./results/somatic_mutation_analysis/functional_impact_prediction/output")
figures_fip_path <- file.path(current_dir, 
                              "./results/somatic_mutation_analysis/functional_impact_prediction/figures")
results_gor_path <- file.path(current_dir, 
                              "./results/somatic_mutation_analysis/ora_exploration/output")
figures_gor_path <- file.path(current_dir, 
                              "./results/somatic_mutation_analysis/ora_exploration/figures")
interesting_gs_path <- file.path(data_path,
                                 "../interesting_genesets.xlsx")
results_ubi_ora_exp_path <- file.path(current_dir, 
                                      "./results/somatic_mutation_analysis/ubi_ora_exploration/output")
figures_ubi_ora_exp_path <- file.path(current_dir, 
                                      "./results/somatic_mutation_analysis/ubi_ora_exploration/figures")
results_ubi_db_path <- file.path(current_dir, 
                                 "./results/somatic_mutation_analysis/ubi_db_creation/output")
figures_ubi_db_path <- file.path(current_dir, 
                                 "./results/somatic_mutation_analysis/ubi_db_creation/figures")
results_ubfora_path <- file.path(current_dir, 
                                 "./results/somatic_mutation_analysis/ubiquitin_ora_class_go/output")
figures_ubfora_path <- file.path(current_dir, 
                                 "./results/somatic_mutation_analysis/ubiquitin_ora_class_go/figures")
results_protein_seqs_path <- file.path(current_dir, 
                                       "./results/somatic_mutation_analysis/protein_seqs/output")
results_cellular_loc_path <- file.path(current_dir, 
                                       "./results/somatic_mutation_analysis/cellular_loc/output")
figures_cellular_loc_path <- file.path(current_dir, 
                                       "./results/somatic_mutation_analysis/cellular_loc/figures")
results_pfam_path <- file.path(current_dir, 
                               "./results/somatic_mutation_analysis/pfam_analysis/output")
figures_pfam_path <- file.path(current_dir, 
                               "./results/somatic_mutation_analysis/pfam_analysis/figures")
results_expression_path <- file.path(current_dir, 
                                     "./results/expression_analysis/output")
figures_expression_path <- file.path(current_dir, 
                                     "./results/expression_analysis/figures")

# List of paths
paths <- c(
  data_path, results_pc_path, figures_pc_path, results_fip_path, figures_fip_path, 
  results_gor_path, figures_gor_path, results_ubi_ora_exp_path, figures_ubi_ora_exp_path, 
  results_ubi_db_path, figures_ubi_db_path, results_ubfora_path, figures_ubfora_path, 
  results_protein_seqs_path, results_cellular_loc_path, figures_cellular_loc_path, 
  results_pfam_path, figures_pfam_path, results_expression_path, figures_expression_path
)

# Create directories if they don't exist
lapply(paths, function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
})