######################
### INSTALL SCRIPT ###
######################

## Install the required packages for the whole pipeline

# Check and install packages if not already installed
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# List of packages to check and install
packages <- c("rprojroot", "readxl", "tidyverse", "dplyr", "purrr", "knitr", "DT", "httr", "jsonlite", "stringr", 
              "msigdbr", "HGNChelper", "VennDiagram", "circlize", "gridtext", "patchwork", "ggpubr", 
              "ggdendro", "reshape2", "ggnewscale", "cowplot", "pheatmap", "ggplotify", "ggrepel", 
              "broom", "corrplot", "stats", "dendextend", "magrittr")

# Iterate over the list and install missing packages
lapply(packages, install_if_missing)

# Specific check for BiocManager and Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of Bioconductor packages
bioconductor_packages <- c("VariantAnnotation", "biomaRt", "maftools", "GSEABase", "EnsDb.Hsapiens.v75", 
                           "EnsDb.Hsapiens.v86", "fgsea", "limma", "clusterProfiler", "ggupset", "mygene", 
                           "GO.db", "AnnotationDbi", "ComplexHeatmap")

# Iterate over the list and install missing Bioconductor packages
lapply(bioconductor_packages, install_if_missing)

# Load all packages
all_packages <- c(packages, bioconductor_packages)

# Load all packages using library()
lapply(all_packages, function(pkg) {
  library(pkg, character.only = TRUE)
})
