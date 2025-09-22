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
# 1. Get the data from the cleaning ring width file: wildchrokie_rw
source("cleaning/source/cleanRingWidth.R") 
# 2. Get observation data from main repo: obsdata
source("cleaning/source/cleanObsData.R")

temp <- merge(wildchrokie_rw, obsdata, by = c("id", "spp", "year"))
# get only the years we have data for
temp2 <- subset(temp, year %in% c(2018, 2019, 2020))

# 3. Get GDD data from main wildhill repo
source("cleaning/source/combineWeather.R") # may be not needed

# 4. Calculate primary and full growing season GDD
source("cleaning/source/growingSeasonGDD.R")

temp3 <- merge(temp2,
               obsdataWithGDD[, c("id", "year", "pgsGDD", "fgsGDD")],
               by = c("id", "year")
               )

# 3. Get male and female flowering time
# source("fromMainRepo/dichogamy.R")


# merge columns from cgclean and d by Name and drop all names that are not in d
merged_df <- merge(rw, obsdata, by = c("name", "year"), all.x = TRUE)


# Add species column
# merged_df$species<-substr(merged_df$Name, 0,6)

