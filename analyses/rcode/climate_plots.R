# Wildchrokie (and eventually coringTreespotters) climate data with phenology
# CRD 12 March 2026

# housekeeping
# rm(list=ls())
# options(stringsAsFactors = FALSE)
# options(max.print = 150)
# options(digits = 3)

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(wesanderson)
library(patchwork) 
library(sharpshootR)

if (length(grep("christophe_rouleau-desrochers", getwd())) > 0) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if (length(grep("lizzie", getwd())) > 0) {
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
} else  {
  setwd("/home/crouleau/wildchrokie/analyses")
}
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
source('rcode/tools.R')

# flags
makeplots <- FALSE

empir <- read.csv("output/empiricalDataMAIN.csv")
climatesum <- read.csv("output/climateSummariesYear.csv")
climatesummonth <- read.csv("output/climateSummariesByMonth.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

# Day met data. From Wildhell garden repo. It's the climate data for each provenance
climprov <- read.csv("input/_notcookies/dayMet.csv")

# transform my groups to numeric values
empir$site_num <- match(empir$site, unique(empir$site))
empir$spp_num <- match(empir$spp, unique(empir$spp))
empir$year_num <- match(empir$year, unique(empir$year))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate and Phenology summaries ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
empir$leafout <- as.integer(empir$leafout)
empir$budset <- as.integer(empir$budset)

empir$yeardoybudburst <- paste(empir$year, empir$budburst, sep = "_")
empir$yeardoyleafout <- paste(empir$year, empir$leafout, sep = "_")
empir$yeardoybudset <- paste(empir$year, empir$budset, sep = "_")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
##### Temperature min and max #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# common objects across budset and leafout
clim_vars  <- c("TempMeanMax", "TempMeanMean","TempMeanMin")

emp_clim <- merge(empir, climatesum, by = "year", all.x = TRUE)

colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmax")] <- "TempMeanMax"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmean")] <- "TempMeanMean"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmin")] <- "TempMeanMin"

emp_climlo <- emp_clim[!is.na(emp_clim$leafout),]
emp_climbs <- emp_clim[!is.na(emp_clim$budset),]

emp_climlo$anomleafout <- emp_climlo$leafout - mean(emp_climlo$leafout)
emp_climbs$anombudset <- emp_climbs$budset - mean(emp_climbs$budset)

nrow(emp_climlo)
nrow(emp_climbs)
nrow(emp_clim)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
##### Phenology #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
leafoutbyyr <- aggregate(leafout ~ year + latbi, empir, FUN = mean)
budsetbyyr <- aggregate(budset ~ year + latbi, empir, FUN = mean)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
##### Droughts #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
yrwc <- c(2018:2020)
yrts <- c(2016:2024)
climatesummonth$monthname <- month.name[climatesummonth$month]

moderatedrought <- subset(climatesummonth, pdsi < -2 & pdsi > -3)
severedrought <- subset(climatesummonth, pdsi < -3)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Temperatures and PDSI that are used in MS #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
wcgsl <- aggregate(pgsGSL ~ year, empir, FUN = mean, na.rm =TRUE)
wcgsl[order(wcgsl$pgsGSL),]

wcsos <- aggregate(leafout ~ year, empir, FUN = mean, na.rm =TRUE)
wcsos[order(wcsos$leafout),]

wceos <- aggregate(budset ~ year, empir, FUN = mean, na.rm =TRUE)
wceos[order(wceos$budset),]

wcclim <- subset(climatesum, year %in% 2018:2020)
wcdjf <- subset(wcclim, period %in% "DJF")
wcdjf[order(wcdjf$tmeanmean),]
wcmam <- subset(wcclim, period %in% "MAM")
wcmam[order(wcmam$tmeanmean),]
wcjja <- subset(wcclim, period %in% "JJA")
wcjja[order(wcjja$tmeanmean),]
wcson <- subset(wcclim, period %in% "SON")
wcson[order(wcson$tmeanmean),]
wcmam[order(wcmam$pdsi),]
wcjja[order(wcjja$pdsi),]
wcson[order(wcson$pdsi),]

wcrw <- aggregate(lengthCM ~ year, empir, FUN = mean, na.rm =TRUE)
wcrw[order(wcrw$lengthCM),]

agggdd <- aggregate(pgsGDD5 ~ year, empir, FUN = mean)
agggdd[order(agggdd$pgsGDD5),]

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
##### Late spring frosts #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
combined <- merge(empir[, c("treeid", "year", "budburst")], 
                  weldhillclim[, c("year", "doy", "minTempC")], 
                  by = "year")
# Define frost threshold (usually 0 or -2.2 for hardy buds)
frost_threshold <- 0

# Subset to days after budburst where it froze
late_frosts <- combined[combined$doy > combined$budburst & 
                          combined$minTempC <= frost_threshold, ]
combined$year[which(combined$budburst %in% min(combined$budburst, na.rm = TRUE))]
# Remove rows with NA budburst (if any)
late_frosts <- late_frosts[!is.na(late_frosts$budburst), ]

# Create a summary data frame
frost_summary <- aggregate(minTempC ~ treeid + year, 
                           data = late_frosts, 
                           FUN = function(x) c(Count = length(x), MinTemp = min(x)))

# Flatten the aggregate result into a clean data frame
frost_summary <- do.call(data.frame, frost_summary)
colnames(frost_summary) <- c("treeid", "year", "frost_days_after_burst", "minimum_frost_temp")

# View the results
head(frost_summary)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate overlap in 2020 ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (makeplots){
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/input/_notcookies/")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Start with years 2012 to 2020
#from https://labs.arboretum.harvard.edu/weather/
# Uncleaned weldhill climate

weldhill_climate <- read.csv("weldhill.csv")

# pull the years until 2019 --- see issue #8 for selecting until 2019 only
weld_15to20 <- weldhill_climate[grepl("2020", weldhill_climate$Eastern.daylight.time), ]

# make date and time separate cols
weld_15to20$date <- substr(weld_15to20$Eastern.daylight.time,
                           1, nchar(weld_15to20$Eastern.daylight.time) - 9)
weld_15to20$time <- substr(weld_15to20$Eastern.daylight.time,
                           nchar(weld_15to20$Eastern.daylight.time) - 8,
                           nchar(weld_15to20$Eastern.daylight.time))
weld_15to20$date <- as.Date(weld_15to20$date, format = "%m/%d/%Y")

# convert far to celcius: https://www.metric-conversions.org/temperature/fahrenheit-to-celsius.htm
weld_15to20$tempCelcius <- (weld_15to20$Temp..F - 32) / 1.8
weld_15to20$pptMM <- weld_15to20$Rain.in * 25.4

# summarize with min mean and max
weld_15to20_temp <- aggregate(
  tempCelcius ~ date,
  data = weld_15to20,
  FUN = function(x) c(max = max(x, na.rm = TRUE),
                      mean = mean(x, na.rm = TRUE),
                      min = min(x, na.rm = TRUE))
)
weld_15to20_daily <- data.frame(
  date = weld_15to20_temp$date,
  maxTempC = weld_15to20_temp$tempCelcius[, "max"],
  meanTempC = weld_15to20_temp$tempCelcius[, "mean"],
  minTempC = weld_15to20_temp$tempCelcius[, "min"]
)

# summarize for precipitation
weld_15to20_ppt <- aggregate(
  pptMM ~ date,
  data = weld_15to20,
  FUN = sum
)
weld_15to20_daily$pptMM <- weld_15to20_ppt$pptMM[match(weld_15to20_daily$date, weld_15to20_ppt$date)]

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Then follow with years 2020 to 2024
# From https://newa.cornell.edu/all-weather-data-query/
# precipitation in inches
newa_files <- list.files(pattern = "\\.csv$")[2:6] # 2 and 10 is wrong
newa_list <- lapply(newa_files, read.csv, stringsAsFactors = FALSE)
newa_20to25 <- do.call(rbind, newa_list)

# select cols of interest
newa_20to25_sub <- newa_20to25[, c("date",
                                   "Avg.Air.Temp...F.",
                                   "Max.Air.Temp...F.",
                                   "Min.Air.Temp...F.",
                                   "Total.Precipitation")]
newa_20to25_daily <- data.frame(
  date      = newa_20to25_sub$date,
  maxTempC  = (newa_20to25_sub$Max.Air.Temp...F. - 32) * 5/9,
  meanTempC = (newa_20to25_sub$Avg.Air.Temp...F. - 32) * 5/9,
  minTempC  = (newa_20to25_sub$Min.Air.Temp...F. - 32) * 5/9,
  pptMM     = newa_20to25_sub$Total.Precipitation * 25.4
)
newa_20to25_daily$date <- as.Date(newa_20to25_daily$date, format = "%m/%d/%Y")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Bind both dfs together but first give them a unique id
weld_15to20_daily$source <- "weldhillsite"
newa_20to25_daily$source <- "newa"
weldhill_comp <- rbind(weld_15to20_daily, newa_20to25_daily)

# removing temp of -500
weldhill_comp <- subset(weldhill_comp, minTempC > -50)
weldhill_comp$year <- substr(weldhill_comp$date, 1,4)
weldhill_comp$doy  <- as.numeric(format(weldhill_comp$date, "%j"))

# make a quick and dirty plot to compare overlapping data for 2020
y2020 <- subset(weldhill_comp, year == 2020)
y2020 <- subset(y2020, doy < 240)

src_levels <- unique(y2020$source)
src_cols <- setNames(c("#FF0000", "#00A08A")[seq_along(src_levels)], src_levels)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")
pdf("figures/climate/climateComparison2020.pdf", width = 8, height = 5)
par(mar = c(2,4,4,4))
plot(y2020$doy, y2020$meanTempC,
     ylim = c(min(y2020$minTempC), max(y2020$maxTempC)),
     xlim = c(0, 250),
     type = "n",
     xlab = "Day of year",
     ylab = expression(paste("Temperature (", degree, "C)")),
     frame = TRUE, bty = "l")

for (s in src_levels) {
  y2020_s <- y2020[y2020$source == s, ]
  y2020_s <- y2020_s[order(y2020_s$doy), ]
  
  polygon(
    c(y2020_s$doy, rev(y2020_s$doy)),
    c(y2020_s$minTempC, rev(y2020_s$maxTempC)),
    col = adjustcolor(src_cols[s], alpha.f = 0.2),
    border = NA
  )
  
  lines(y2020_s$doy, y2020_s$meanTempC, col = src_cols[s], lwd = 2)
}

legend("topleft",
       legend = src_levels,
       col = src_cols,
       lwd = 2, bty = "n",
       title = "Dataset")

dev.off()
}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate at each provenance ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
climprov$tmean <- (climprov$tmin + climprov$tmax) / 2

# Mean annual temperature --- --- --- --- ------ --- --- --- ------ --- --- 
mat <- aggregate(tmean ~ site, climprov, FUN = mean)

# Frost free days (from kens code in wildhellgarden repo) --- --- --- --- ---
climprov$datetime <- as.Date(sprintf('%s %s', climprov$year, climprov$doy), format="%Y %j")
climprov$value <- climprov$tmin
climprov <- subset(climprov, year < 2022)
sites <- unique(climprov$site)

ffdData <- vector()
ffdMeta <- vector()

for(s in 1:length(sites)){
  
  siteData <- subset(climprov, site == sites[s])
  siteYears <- unique(siteData$year)
  # cat(s, sites[s], nrow(siteData), "\n")
  
  for(y in 1:length(siteYears)){
    
    dayData <- subset(siteData, year == siteYears[y])
    
    if(nrow(dayData) == 0) next
    
    ffd <- FFD(dayData, returnDailyPr = TRUE, frostTemp = 0)
    
    ffdData <- rbind(ffdData, ffd$summary)
    ffdMeta <- rbind(ffdMeta, data.frame(site = sites[s], year = siteYears[y]))
  }
}

f <- cbind(ffdMeta, ffdData)

# calculate average ffd
ffd <- aggregate(ffd.90 ~ site, f, FUN = mean)
