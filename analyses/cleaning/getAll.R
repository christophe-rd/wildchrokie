# Started 16 May 2025
# ÇRD

# write out a raw data file for phenology and ring width measurements

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

### === === === === === === === === === === ###
# Get cleaned data from this repo #
### === === === === === === === === === === ###
# 1. Get the data from the cleaning ring width file
source("cleaning/source/cleanRingWidth.R") 

### === === === === === === === === === === ###
# Get cleaned data from main repo #
### === === === === === === === === === === ###
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
# 1. Get observation data from main repo
source("fromMainRepo/source/cleaning_obsdata.R")
# 3. Get GDD data from main wildhill repo
source("fromMainRepo/source/combineWeather.R") # for that csv: gddData.csv

# 2. Grab climate data
# source("fromMainRepo/climatedata.R") # commenting it out for now as there are only 4 years and might be an old script
# 3. Get male and female flowering time
# source("fromMainRepo/dichogamy.R")


# merge columns from cgclean and d by Name and drop all names that are not in d
merged_df <- merge(d, obsdata, by = c("Name", "Year"), all.x = TRUE)
# Add species column
merged_df$species<-substr(merged_df$Name, 0,6)

