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

# merge dfs for ring width with phenology data 
nrow(wildchrokie_rw)
nrow(obsdata)
temp <- merge(wildchrokie_rw, obsdata, by = c("treeid", "spp", "year"))
rw <- subset(wildchrokie_rw, year %in% c(2018, 2019, 2020))
nrow(rw)

# get only the years we have data for
temp2 <- subset(temp, year %in% c(2018, 2019, 2020))
nrow(temp2)

temp3 <- merge(temp2,
               obsdataWithGDD[, c("treeid", "year", "pgsGDD", "pgsGDDAVG", "fgsGDD", "fullGDD")],
               by = c("treeid", "year")
               )
nrow(temp3)
temp3$idyear <- paste(temp3$treeid, temp3$year, sep = "_")
temp4 <- temp3[!duplicated(temp3$idyear),] 
nrow(temp4)

temp4$latbi <- NA
temp4$latbi[which(temp4$spp == "ALNINC")] <- "Alnus incana"
temp4$latbi[which(temp4$spp == "BETPAP")] <- "Betula papyrifera"
temp4$latbi[which(temp4$spp == "BETPOP")] <- "Betula populifolia"
temp4$latbi[which(temp4$spp == "BETALL")] <- "Betula alleghaniensis"

# write csv
write_csv(temp4, "output/empiricalDataMAIN.csv")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Now without tree rings!
obsdataWithGDD$siteplot <- sub("^[^_]*_(.*?)_.*$", "\\1", obsdataWithGDD$treeid)
obsdataWithGDD$site <- substr(obsdataWithGDD$siteplot, 0,2)
obsdataWithGDD$plot <- substr(obsdataWithGDD$siteplot, 3,4)
obsdataWithGDD$replicate <- sub(".*_", "", obsdataWithGDD$treeid)

obsdataWithGDD <- obsdataWithGDD[order(
  obsdataWithGDD$spp,
  obsdataWithGDD$site,
  obsdataWithGDD$replicate,
  obsdataWithGDD$year
), ]

write_csv(obsdataWithGDD, "output/empiricalDataNORing.csv")
