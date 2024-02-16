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
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("VariantAnnotation")
BiocManager::install("biomaRt")
BiocManager::install("maftools")
BiocManager::install("GSEABase")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("fgsea")
