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