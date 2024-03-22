######################
### INSTALL SCRIPT ###
######################

## Install the required packages for the whole pipeline

if (!require("readxl")) {
  install.packages("readxl")
}
if (!require("tidyverse")) {
  install.packages("tidyverse")
}
if (!require("dplyr")) {
  install.packages("dplyr")
}
if (!require("purrr")) {
  install.packages("purrr")
}
if (!require("knitr")) {
  install.packages("knitr")
}
if (!require("DT")) {
  install.packages("DT")
}
if (!require("httr")) {
  install.packages("httr")
}
if (!require("jsonlite")) {
  install.packages("jsonlite")
}
if (!require("stringr")) {
  install.packages("stringr")
}
if (!require("msigdbr")) {
  install.packages("msigdbr")
}
if (!require("HGNChelper")) {
  install.packages("HGNChelper")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!("VariantAnnotation" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("VariantAnnotation", update = FALSE)
}
if (!("biomaRt" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("biomaRt", update = FALSE)
}
if (!("maftools" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("maftools", update = FALSE)
}
if (!("GSEABase" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GSEABase", update = FALSE)
}
if (!("EnsDb.Hsapiens.v75" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnsDb.Hsapiens.v75", update = FALSE) # From Ensembl
}
if (!("fgsea" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("fgsea", update = FALSE)
}
if (!("limma" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("limma", update = FALSE)
}
if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}
if (!("ggupset" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ggupset", update = FALSE)
}
if (!("mygene" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("mygene", update = FALSE)
}
if (!("GO.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GO.db", update = FALSE)
}
if (!("AnnotationDbi" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("AnnotationDbi", update = FALSE)
}