# Started 16 May 2025
# ÇRD

# write out a raw data file for phenology and ring width measurements

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
# grab a vec of interested species
vec <- c("ALNINC", "BETALL", "BETPAP", "BETPOP")
### === === === === === === === === === === ###
# Get cleaned data from this repo #
### === === === === === === === === === === ###
# 1. Get the data from the cleaning ring width file
source("cleaning/source/cleanRingWidth.R") 
# 2. Get phenology observations from main wildhill repo

# 3. Get GDD data from main wildhill repo


### === === === === === === === === === === ###
# Get cleaned data from main repo #
### === === === === === === === === === === ###
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
# 1. Get observation data from main repo
source("fromMainRepo/cleaning_obsdata.R")
# get only 4 species
cgclean2 <- subset(cgclean, spp %in% vec)
# select columns
cgclean3 <- cgclean2[, c(1:2, 6:ncol(cgclean2))]
# create individual name to fit the one I am using 

# 2. Grab climate data
# source("fromMainRepo/climatedata.R") # commenting it out for now as there are only 4 years and might be an old script
# 3. Get male and female flowering time
source("fromMainRepo/dichogamy.R")


# merge columns from cgclean and d by Name and drop all names that are not in d
merged <- merge(cgclean3, d, by = "Name")

