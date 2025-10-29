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

# 3. Clean climate datasets: weldhillcomp
source("cleaning/source/cleanClimate.r") 

# 4. Calculate primary and full growing season GDD
source("cleaning/source/calculateGrowingSeasonGDD.R")

# merge dfs!
temp <- merge(wildchrokie_rw, obsdata, by = c("treeid", "spp", "year"))

# get only the years we have data for
temp2 <- subset(temp, year %in% c(2018, 2019, 2020))

temp3 <- merge(temp2,
               obsdataWithGDD[, c("treeid", "year", "pgsGDD", "fgsGDD")],
               by = c("treeid", "year")
               )

# make some checks
unique(temp3$sampleType)
# remove 2 samples for a given individuals
temp4 <- subset(temp3, sampleType != "coresWithCookies") 
temp5 <- temp4[!is.na(temp4$pgsGDD),]
temp5$idyear <- paste(temp5$treeid, temp5$year, sep = "_")
temp6 <- temp5[!duplicated(temp5$idyear),] 

# write csv
nrow(temp6)
y2020 <- subset(temp4, year == "2020")
y2020$idyear <- paste(y2020$treeid, y2020$year, sep = "_")
y2020 <- y2020[which(!is.na(y2020$leafout) & !is.na(y2020$budset)),]
y2020 <- y2020[!duplicated(y2020$idyear),]
nrow(y2020)
unique(y2020$id)

# merged_df$species<-substr(merged_df$Name, 0,6)

