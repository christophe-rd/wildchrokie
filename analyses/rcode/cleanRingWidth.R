## Started 3 April 2025 ##
# Christophe

# goal of this script is to import the cookie measurements, clean them and change unit from inches to cm

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## Load Libraries
library(ggplot2)

# Set Working Directory
# Christophe's directory
directory <- "/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/input"

# list files in directory
file_list <- list.files(path = directory, pattern = "\\.csv$", full.names = TRUE)

# read all files
data_list <- lapply(file_list, read.csv)

# combine tables
combined_data <- do.call(rbind, data_list)

