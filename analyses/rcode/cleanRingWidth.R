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
c <- do.call(rbind, data_list)

### === === === === === ###
# Clean labels and years #
### === === === === === ###
# cleaned name column
c$Name <- gsub("_guides(_[0-9]+)?\\.tif", "", sub(":.*", "", c$Label))
# cleaned year column
c$Year <- sub(".*:", "", c$Label)
# Create a new table with only year values in the Year column i.e. excluding comments
e <- c[grepl("^\\d{2,4}$", c$Year), ] 
# Convert year column to numeric
e$Year <- as.numeric(e$Year)
# Create a new column with the length in cm with a conversion factor of 2.54cm because 1 inch is 2.540005 cm
e$LengthCM <- e$Length*2.54
# Create a new df with only necessary columns
d <- e[, c("Name", "Year", "LengthCM")]

### Select random measurements that I will physically measure on the cookies to verify the accuracy of the measurements
# Select 20 random measurements
set.seed(123) # for reproducibility0
# Select 20 random rows from the data frame
random_rows <- sample(nrow(d), 10)
# Create a new data frame with the selected rows  
vec <- d$Name[random_rows]
# Create a new data frame with the selected rows
random_measurements <- d[vec, ]

newdf <- subset(d, Name %in% vec)
