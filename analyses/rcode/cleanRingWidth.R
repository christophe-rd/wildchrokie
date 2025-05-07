## Started 3 April 2025 ##
# Christophe

# goal of this script is to import the cookie measurements, clean them and change unit from inches to cm

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## Load Libraries
library(ggplot2)
library(xlsx)

# Set main directory
directory <- "/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/input"

# List all items and keep only directories
items <- list.files(path = directory, full.names = TRUE)
folders <- items[file.info(items)$isdir]

# Get folder names and exclude undesired ones
folder_names <- basename(folders)
valid_folders <- folder_names[!folder_names %in% c("_notcookies", "_readme.md")]

# Initialize list to store all data
all_data <- list()

# Loop through each folder
for (folder in valid_folders) {
  folder_path <- file.path(directory, folder)
  file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Read all CSVs and tag each row with the folder name
  data_list <- lapply(file_list, function(file) {
    df <- read.csv(file)
    df$sourceFolder <- folder
    return(df)
  })
  
  all_data <- c(all_data, data_list)
}

# Combine all into a single data frame
c <- do.call(rbind, all_data)

### === === === === === ###
##### Read initial table that compiled what tree we had cookies/cores ##### 
# it involves a couple of mistakes when I entered the data that I will fix bellow #
### === === === === === ###
# read table
setwd(directory)
list.files()
d <- read.xlsx("_notcookies/treecookies.xlsx", sheetName = "Sheet1")
# remove _ between provenance and number
d$id <- gsub("(_)([A-Z]+)_([0-9]+[A-Z]?)$", "\\1\\2\\3", d$id)
# paste id and plot
d$idfull <- paste(d$id, paste0("P", d$Plot), sep = "_")
# subset only for species I want
vec <- c("ALNINC", "BETALL", "BETPAP", "BETPOP", "QUERUB")
dsub <- subset(d, species %in% vec )

# List the cores for which we have #no cookies
# ALNINC_SH_6_P5
# BETALL_SH_4A_9
# BETPAP HF16 P2
# BETPOP_HF_5_6
# ALNINC HF8 P12

# List the cores for which we have cookies but were not entered as such!
# BETALL_GR_9_P3

# List the cores for which i dont have a core matching this name yet. 
# BETALL SH9 P9, I have a PNA TO CHECK
# BETALL_WM_8_9
# BETPOP_GR_5_P6

# Clean cookie and cores when I made mistakes in entering them
### Remove cookie for ALNINC HF9 P6
dsub$cookie.[which(dsub$idfull == "ALNINC_HF_9_6")] <- "0" 
### Remove cookie for ALNINC WM2B P1
dsub$cookie.[which(dsub$idfull == "ALNINC_WM_2B_1")] <- "0"
### Change number of cookies to 2 for BETALL GR13 P15
dsub$cookie.[which(dsub$idfull == "BETALL_GR_13_15")] <- "2"
### Change number of cookies to 2 for BETALL GR13 P1
dsub$cookie.[which(dsub$idfull == "BETALL_GR_13_1")] <- "0"
### Change number of cookie for BETALL_GR9_P3
dsub$cookie.[which(dsub$idfull == "BETALL_GR9_P3")] <- "1"
### Change number of cookie for ALNINC_GR8B_P5
dsub$core[which(dsub$idfull == "ALNINC_GR8B_P5")] <- "1"

cookiesOG <- subset(dsub, cookie. == "1")

coresOG <- subset(dsub, core %in% c("1","2"))

### === === === === === ###
##### Clean labels and years #####
### === === === === === ###
# cleaned name column
c$Name <- gsub("[-_]?guides([-_][0-9]+)?\\.tif$", "", sub(":.*", "", c$Label))
# cleaned year column
c$Year <- sub(".*:", "", c$Label)
# Create a new table with only year values in the Year column i.e. excluding comments
e <- c[grepl("^\\d{2,4}$", c$Year), ] 
# Convert year column to numeric
e$Year <- as.numeric(e$Year)
# Create a new column with the length in cm with a conversion factor of 2.54cm because 1 inch is 2.540005 cm
e$LengthCM <- e$Length*2.54
# Create a new df with only necessary columns
d <- e[, c("Name", "Year", "LengthCM", "sourceFolder")]

#### Extra cleaning steps ####
# clean names that I entered wrong when scanning
# BETALL_WM8B_P95: decision made in CoreXCookies.md
d$Name[which(d$Name == "BETALL_WM8B_P95")] <- "BETALL_WM8B_P5" 
# BETALL_WM8D_P5: decision made in CoreXCookies.md
d$Name[which(d$Name == "BETALL_WM8D_P5")] <- "BETALL_WM8B_P5" 
# BETALL_WM_9_P5
d$Name[which(d$Name == "BETALL_WM_9_P5")] <- "BETALL_WM9_P5" 
# BETALL_WM_8B_P12
d$Name[which(d$Name == "BETALL_WM_8B_P12")] <- "BETALL_WM8B_P12" 
# BETPOP_WM-7_P2
d$Name[which(d$Name == "BETPOP_WM-7_P2")] <- "BETPOP_WM7_P2" 
# BETPOP_GRB_P12: decision made in CoreXCookies.md
d$Name[which(d$Name == "BETPOP_GRB_P12")] <- "BETPOP_GR5B_P12" 
# BETPOP_GR_5B_P12
d$Name[which(d$Name == "BETPOP_GR_5B_P12")] <- "BETPOP_GR5B_P12" 
# ALNINC_WM_2A_P1
d$Name[which(d$Name == "ALNINC_WM_2A_P1")] <- "ALNINC_WM2A_P1" 
# ALNINC_WM_5B_P2
d$Name[which(d$Name == "ALNINC_WM_5B_P2")] <- "ALNINC_WM5B_P2" 
# BETPAP_HF_16A_P6
d$Name[which(d$Name == "BETPAP_HF_16A_P6")] <- "BETPAP_HF16A_P6" 
# BETPOP_GR_5A_P3
d$Name[which(d$Name == "BETPOP_GR_5A_P3")] <- "BETPOP_GR5A_P3" 
# BETPOP_GR_5B_P11
d$Name[which(d$Name == "BETPOP_GR_5B_P11")] <- "BETPOP_GR5B_P11" 
# BETPOP_GR_5C_P12: need to find what it is
d$Name[which(d$Name == "BETPOP_GR_5C_P12")] <- "BETPOP_GR5C_P12" 
# BETPOP_HF__P6
d$Name[which(d$Name == "BETPOP_HF__P6")] <- "BETPOP_HFNA_P6"
# BETAL_GR9_P3
d$Name[which(d$Name == "BETAL_GR9_P3")] <- "BETALL_GR9_P3"
# BETALL_SH9_PNA: decision made in CoreXCookies.md
d$Name[which(d$Name == "BETALL_SH9_PNA")] <- "BETALL_SH9_P9"
# BETALL_WM8_P12: 
d$Name[which(d$Name == "BETALL_WM8_P12")] <- "BETALL_WM8B_P12" 

### === === === === === ###
# Verification steps #
### === === === === === ###
# compare if I have the data from all cookies in the og dataset
vcook <- c("cookies", "cookiesUnconfident")
cookies <- subset(d, sourceFolder %in% vcook)
listCookieNames <- unique(cookies$Name)
cookieInOG <- listCookieNames[which(!listCookieNames%in%cookiesOG$idfull)]
### investigate the mistmatches
# ALNINC_HF3_P16: not in table: to check if the cookie really exist
# ALNINC_SH5_P2: not added in table : to check
# ALNINC_WM8_P1: not in table: to check if the cookie really exist
# BETPAP_GR5B_P2: not added in table : to check
# BETPOP_GR5_P6: not added in table : to check

# compare if I have scanned the cookies from the data from og dataset 
OGtoCookies <- cookiesOG$idfull[which(!cookiesOG$idfull%in%listCookieNames)]

### investigate the mistmatches
# ALNINC_GR11B_P3: scanned, but not sure and have to look under microscope or perhaps comparing with other alninc for which I am more certain of the calls
# ALNINC_GR8_P3: scanned, but not sure and have to look under microscope
# ALNINC_HF9_P6: look if cookie exists
# ALNINC_SH4_P2: scanned, but not sure and have to look under microscope or perhaps comparing with other alninc for which I am more certain of the calls
# ALNINC_SH5_P3: could be the same as SH5_P2, TO LOOK at the core and at the cookie
# ALNINC_WM1_P1: scanned, but not sure and have to look under microscope or perhaps comparing with other alninc for which I am more certain of the calls. Also, entered as WM_P1 not WM1
# ALNINC_WM2B_P1: need to add guidelines
# ALNINC_WM5_P1: look for cookie
# ALNINC_WM5A_P16: look for cookie
# BETALL_GR13_P1: core without cookie
# BETALL_SH9_P6: scanned, but not sure and have to look under microscope or perhaps comparing with other alninc for which I am more certain of the calls.
# BETPAP_GR8_P6: look for cookie
# BETPAP_HF16_P3: look for cookie
# BETPOP_GR3_P6: look for cookie

# compare if I have the data from all cores in the og dataset
vcores <- c("coresUnconfident", "coresWithCookies", "coresWithoutCookies")
cores <- subset(d, sourceFolder %in% vcores)
listCoreNames <- unique(cores$Name)
coresinOG <- listCoreNames[which(!listCoreNames%in%coresOG$idfull)]

# BETALL_SH5_P6: look if its not SH6 P5 instead
# ALNINC_HF9_P6: verify if it's the right tag.

# compare if I have all the cores in the og dataset have been scanned
Ogscanned <- coresOG$idfull[which(!coresOG$idfull%in%listCoreNames)]

# BETALL_SH9_P6: verifiy if core exists
# BETALL_WM8_P1: OK, core exists but poor quality so no csv
# BETPAP_HF16_P2: will be rescanned 
# BETPOP_GR3_P16: RESAND and find rings. It doesn't work!
# BETPOP_GR5_P6: RESAND and RESCAN
# BETPOP_GR5C_P12: look if it exists
# BETPOP_HF3_P1: look if it exists
# BETPOP_WM7_P1: look if it exists

### === === === === === ###
# Start playing with the data #
### === === === === === ###
# Take the mean of each year and each tree replicate
agg <- aggregate(d$LengthCM, by = list(d$Name, d$Year), FUN = mean)
# rename the columns
colnames(agg) <- c("Name", "Year", "LengthCM")
# select one replicate
ALNINC_WM_2A_P1 <- subset(agg, Name == "ALNINC_WM_2A_P1")
# calculate diameter by adding the mean ring width of each year
radius <- sum(ALNINC_WM_2A_P1$LengthCM)
radius*2 
### Select random measurements that I will physically measure on the cookies to verify the accuracy of the measurements
# Select 20 random measurements
set.seed(123) # for reproducibility0
# Select 20 random rows from the data frame
random_rows <- sample(nrow(d), 10)
# Create a new data frame with the selected rows  
vec <- d$Name[random_rows]
# Create a new data frame with the selected rows
random_measurements <- subset(agg, Name %in% vec)

newdf <- subset(d, Name %in% vec)

