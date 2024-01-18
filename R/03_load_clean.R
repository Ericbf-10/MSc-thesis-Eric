###########################
### LOAD + CLEAN SCRIPT ###
###########################

## Load the data and clean/wrangle

# Imports
library(readxl)
library(tidyverse)
library(dplyr)

# Functions
source(file = "R/02_functions.R")

# Define paths
current_dir <- getwd()
sup1_data <- "data/_raw/41467_2017_1460_MOESM4_ESM_clinical.xlsx"
sup1_path <- file.path(current_dir, sup1_data)

mm909_data <- "data/_raw/MM909_data.RData"
mm909_path <- file.path(current_dir, rdata_data)

provean_zip_path <- "data/_raw/provean_output.zip"
data_path <- "data"
unzip(zipfile = provean_zip_path, exdir = data_path)
provean_path <- "data/provean_output/provean_output"

# Read data
sup1_table <- read_excel(sup1_path, 
                 skip=1, 
                 col_names=TRUE)
colnames(sup1_table)[1] <- "SampleID"

load(mm909_path)

# Subset sup1_table
sup1_subset <- sup1_table %>%
  select(Sample.ID, 
         Sample.Timepoint, 
         AJCC.Stage, 
         PFS.Time, 
         PFS.Event, 
         OS.Time, 
         OS.Event, 
         Type.of.Lesion, 
         Type.of.Primary, RECIST)

sup1_response <- sup1_subset %>%
  mutate(Patient.Response = case_when(
    RECIST %in% c("CR", "PR", "SD") & OS.Time >= 24 ~ "R",
    RECIST %in% c("SD", "PD") & OS.Time < 24 ~ "NR",
    TRUE ~ NA_character_  # This line is for any cases that don't match above conditions
  ))

sup1_final <- sup1_response %>%
  left_join(response %>% 
              select(patient, mut_load), 
            by = c("Sample.ID" = "patient"))

head(sup1_final)

