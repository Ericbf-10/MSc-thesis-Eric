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