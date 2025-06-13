# Started 16 May 2025
# CRD

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
# 2. Get observation data from main repo
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses") # for now... its bad but i need to back down one folder in order for this code to work
source("cleaning/source/cleaning_obsdata.R")
# 3. Get GDD data from main wildhill repo
# source("cleaning/source/combineWeather.R") # may be not needed 


# 2. Grab climate data
# source("fromMainRepo/climatedata.R") # commenting it out for now as there are only 4 years and might be an old script
# 3. Get male and female flowering time
# source("fromMainRepo/dichogamy.R")


# merge columns from cgclean and d by Name and drop all names that are not in d
merged_df <- merge(rw, obsdata, by = c("name", "year"), all.x = TRUE)
# Add species column
# merged_df$species<-substr(merged_df$Name, 0,6)

