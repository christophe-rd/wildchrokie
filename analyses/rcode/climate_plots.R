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
gddyr <- read.csv("output/gddByYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

# Day met data. From Wildhell garden repo. It's the climate data for each provenance
climprov <- read.csv("input/_notcookies/dayMet.csv")

# transform my groups to numeric values
empir$site_num <- match(empir$site, unique(empir$site))
empir$spp_num <- match(empir$spp, unique(empir$spp))
empir$year_num <- match(empir$year, unique(empir$year))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate data #### 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
empir$leafout <- as.integer(empir$leafout)
empir$budset <- as.integer(empir$budset)

gddyr$yeardoy <- paste(gddyr$year, gddyr$doy, sep = "_")
empir$yeardoybudburst <- paste(empir$year, empir$budburst, sep = "_")
empir$yeardoyleafout <- paste(empir$year, empir$leafout, sep = "_")
empir$yeardoybudset <- paste(empir$year, empir$budset, sep = "_")

yearcolors <- c("#931e18", "#da7901", "#247d3f")
par(mfrow = c(1,1))
years <- sort(unique(empir$year))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
years      <- sort(unique(empir$year))
firststeps <- colorRampPalette(c("#9cc184", "#192813"))(length(years))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate summaries ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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


# define plot objects and stuff
species_order <- c(
  "Alnus incana", 
  "Betula alleghaniensis", 
  "Betula papyrifera", 
  "Betula populifolia")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Climate overlap in 2020 #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
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
# Phenology ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
leafoutbyyr <- aggregate(leafout ~ year + latbi, empir, FUN = mean)
budsetbyyr <- aggregate(budset ~ year + latbi, empir, FUN = mean)

leafoutbyyrts <- aggregate(leafout ~ year + latbi, empts, FUN = mean)
colleavesbyyrts <- aggregate(coloredLeaves ~ year + latbi, empts, FUN = mean)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Droughts ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
yrwc <- c(2018:2020)
yrts <- c(2016:2024)
climatesummonth$monthname <- month.name[climatesummonth$month]

moderatedrought <- subset(climatesummonth, pdsi < -2 & pdsi > -3)
severedrought <- subset(climatesummonth, pdsi < -3)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Tests with stan the figures with stan ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Wildchrokie #####
if (FALSE){

clim_vars <- c("TempMeanMax", "TempMeanMin", "TempMeanMean")
climvar <- clim_vars[1] 
period <- "MAM"

d <- emp_climlo[emp_climlo$period == period & !is.na(emp_climlo$TempMeanMax) & 
                !is.na(emp_climlo$anomleafout), ]

d$TempMeanMean <- (d$TempMeanMean - mean(d$TempMeanMean)) / sd(d$TempMeanMean)

ainc <- subset(d, spp %in% "ALNINC")
ball <- subset(d, spp %in% "BETALL")
bpap <- subset(d, spp %in% "BETPAP")
bpop <- subset(d, spp %in% "BETPOP")

# transform data in vectors for gsl
data <- list(
  y = d$anomleafout,
  N = nrow(d),
  year = as.numeric(as.character(d$year_num)),
  species = as.numeric(as.character(d$spp_num)),
  site = as.numeric(as.character(d$site_num)),
  Nspp = length(unique(d$spp_num)),
  Nsite = length(unique(d$site_num)),
  Nyear = length(unique(d$year_num)),
  climpredictor = d$TempMeanMean
)

datainc <- list(
  y = ainc$anomleafout,
  N = nrow(ainc),
  year = as.numeric(as.character(ainc$year_num)),
  species = as.numeric(as.character(ainc$spp_num)),
  site = as.numeric(as.character(ainc$site_num)),
  Nspp = length(unique(ainc$spp_num)),
  Nsite = length(unique(ainc$site_num)),
  Nyear = length(unique(ainc$year_num)),
  climpredictor = ainc$TempMeanMean
)

dataall <- list(
  y = ball$anomleafout,
  N = nrow(ball),
  year = as.numeric(as.character(ball$year_num)),
  species = as.integer(factor(ball$spp_num)),
  site = as.numeric(as.character(ball$site_num)),
  Nspp = length(unique(ball$spp_num)),
  Nsite = length(unique(ball$site_num)),
  Nyear = length(unique(ball$year_num)),
  climpredictor = ball$TempMeanMean
)

datapap <- list(
  y = bpap$anomleafout,
  N = nrow(bpap),
  year = as.numeric(as.character(bpap$year_num)),
  species = as.integer(factor(bpap$spp_num)),
  site = as.numeric(as.character(bpap$site_num)),
  Nspp = length(unique(bpap$spp_num)),
  Nsite = length(unique(bpap$site_num)),
  Nyear = length(unique(bpap$year_num)),
  climpredictor = bpap$TempMeanMean
)

datapop <- list(
  y = bpop$anomleafout,
  N = nrow(bpop),
  year = as.numeric(as.character(bpop$year_num)),
  species =  as.integer(factor(bpop$spp_num)),
  site = as.numeric(as.character(bpop$site_num)),
  Nspp = length(unique(bpop$spp_num)),
  Nsite = length(unique(bpop$site_num)),
  Nyear = length(unique(bpop$year_num)),
  climpredictor = bpop$TempMeanMean
)


# Fit models
climmodel <- stan_model("stan/climatePredictors.stan")
fit <- sampling(climmodel, data = data, iter = 2000, chains = 4, cores = 4)
fitainc <- sampling(climmodel, data = datainc, iter = 2000, chains=4, cores = 4)
fitball <- sampling(climmodel, data = dataall, iter = 2000, chains=4, cores = 4)
fitbpap <- sampling(climmodel, data = datapap, iter = 2000, chains=4, cores = 4)
fitbpop <- sampling(climmodel, data = datapop, iter = 2000, chains=4, cores = 4)

diagnostics <- util$extract_hmc_diagnostics(fitball) 
util$check_all_hmc_diagnostics(diagnostics)

df_fitainc <- as.data.frame(fitainc)
df_fitball <- as.data.frame(fitball)
df_fitbpap <- as.data.frame(fitbpap)
df_fitbpop <- as.data.frame(fitbpop)

###### Plot posterior vs priors for gdd fit ######
pdf(file = "figures/climate/climateModelPriorVSPosterior.pdf", 
    width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(2, 2))
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
ayear_df <- df_fit[, columns[grepl("ayear", columns)]]

# sigma_y
plot(density(df_fit[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fit[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fit[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-100, 100), ylim = c(0, 0.05))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fit[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fit[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-50, 50),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 0.5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

###### Plot posterior vs priors for gdd fit AINC######
pdf(file = "figures/climate/climateModelPriorVSPosteriorAinc.pdf", 
    width = 8, height = 10)

columns <- colnames(df_fitainc)[!grepl("prior", colnames(df_fitainc))]
sigma_df <- df_fitainc[, columns[grepl("sigma", columns)]]
aspp_df <- df_fitainc[, columns[grepl("aspp", columns)]]
bspp_df <- df_fitainc[, columns[grepl("bsp", columns)]]
site_df <- df_fitainc[, columns[grepl("asite", columns)]]
ayear_df <- df_fitainc[, columns[grepl("ayear", columns)]]

par(mfrow = c(2, 2))

# sigma_y
plot(density(df_fitainc[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitainc[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitainc[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-100, 100), ylim = c(0, 0.05))
lines(density(aspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fitainc[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitainc[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-50, 50),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 0.5))
lines(density(bspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

###### Plot posterior vs priors for gdd fit BALL ######
pdf(file = "figures/climate/climateModelPriorVSPosteriorBall.pdf", 
    width = 8, height = 10)

par(mfrow = c(2, 2))

columns <- colnames(df_fitball)[!grepl("prior", colnames(df_fitball))]
sigma_df <- df_fitball[, columns[grepl("sigma", columns)]]
aspp_df <- df_fitball[, columns[grepl("aspp", columns)]]
bspp_df <- df_fitball[, columns[grepl("bsp", columns)]]
site_df <- df_fitball[, columns[grepl("asite", columns)]]
ayear_df <- df_fitball[, columns[grepl("ayear", columns)]]

# sigma_y
plot(density(df_fitball[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitball[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitball[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-100, 100), ylim = c(0, 0.05)) 
lines(density(aspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fitball[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitball[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-50, 50),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 0.5))
lines(density(bspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()


###### Plot posterior vs priors for gdd fit BPAP######
pdf(file = "figures/climate/climateModelPriorVSPosteriorBpap.pdf", 
    width = 8, height = 10)

par(mfrow = c(2, 2))

columns <- colnames(df_fitbpap)[!grepl("prior", colnames(df_fitbpap))]
sigma_df <- df_fitbpap[, columns[grepl("sigma", columns)]]
aspp_df <- df_fitbpap[, columns[grepl("aspp", columns)]]
bspp_df <- df_fitbpap[, columns[grepl("bsp", columns)]]
site_df <- df_fitbpap[, columns[grepl("asite", columns)]]
ayear_df <- df_fitbpap[, columns[grepl("ayear", columns)]]

# sigma_y
plot(density(df_fitbpap[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitbpap[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitbpap[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-100, 100), ylim = c(0, 0.05))
lines(density(aspp_df), col = pal[2], lwd = 1) 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fitbpap[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitbpap[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-50, 50),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 0.5))
lines(density(bspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()


###### Plot posterior vs priors for gdd fit BPOP ######
pdf(file = "figures/climate/climateModelPriorVSPosteriorBpop.pdf", 
    width = 8, height = 10)

par(mfrow = c(2, 2))

columns <- colnames(df_fitbpop)[!grepl("prior", colnames(df_fitbpop))]
sigma_df <- df_fitbpop[, columns[grepl("sigma", columns)]]
aspp_df <- df_fitbpop[, columns[grepl("aspp", columns)]]
bspp_df <- df_fitbpop[, columns[grepl("bsp", columns)]]
site_df <- df_fitbpop[, columns[grepl("asite", columns)]]
ayear_df <- df_fitbpop[, columns[grepl("ayear", columns)]]

# sigma_y
plot(density(df_fitbpop[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitbpop[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitbpop[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-100, 100), ylim = c(0, 0.05))
lines(density(aspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fitbpop[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitbpop[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-50, 50),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 0.5))
lines(density(bspp_df), col = pal[2], lwd = 1)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### PPC #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
samples <- util$extract_expectand_vals(fit)

# whole distribution
par(mfrow = c(1,1))
util$plot_hist_quantiles(samples, "y_rep", 
                         -50, # lower x axis limit
                         50, # upper x axis limit
                         3, # binning
                         baseline_values = data$y,
                         xlab = "anom Leafout")

par(mfrow = c(1, data$Nspp))
for (s in 1:data$Nspp) {
  util$plot_expectand_pushforward(samples[[paste0("aspp[", s, "]")]],
                                  B = 50,
                                  main = paste("aspp species", s))
}
jpeg(
  filename = "figures/climate/retrodictiveDiskSpp.jpeg",
  width = 3600, height = 2000, res = 300          
)
# discs by species
par(mfrow = c(1,data$Nspp))
for (s in unique(data$species)) { # s = 1
  idxs <- which(data$species == s)
  util$plot_disc_pushforward_quantiles(samples,
                                       paste0("y_rep[", idxs, "]"),
                                       baseline_values = data$y[idxs],
                                       ylab = "Leafout",
                                       main = paste("Spp", s))
}
dev.off()
# discs by year
par(mfrow = c(1,data$Nyear))
for (y in unique(data$year)) { # s = 1
  idxs <- which(data$year == y)
  util$plot_disc_pushforward_quantiles(samples,
                                       paste0("y_rep[", idxs, "]"),
                                       baseline_values = data$y[idxs],
                                       ylab = "Leafout",
                                       main = paste("Year", y))
}
# Mu plots 
df_fit <- as.data.frame(fit)

# full posterior
columns <- colnames(df_fit)
columns <- columns[!grepl("prior", columns)]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
ayear_df <- df_fit[, grepl("year", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(ayear) <- 1:ncol(year_df)
colnames(aspp_df) <- 1:ncol(aspp_df)

# posterior summaries
sigma_df2_ts  <- extract_params(df_fit, "sigma", "mean", "sigma")
bspp_df2_ts   <- extract_params(df_fit, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
year_df2_ts <- extract_params(df_fit, "ayear", "fit_ayear", "year", "ayear\\[(\\d+)\\]")
aspp_df2_ts   <- extract_params(df_fit, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")

year_df2_ts$year_name <- emp_climtslo$year[match(year_df2_ts$year, emp_climtslo$year_num)]
aspp_df2_ts$spp_name <- emp_climtslo$latbi[match(aspp_df2_ts$spp, emp_climtslo$spp_num)]
bspp_df2_ts$spp_name <- emp_climtslo$latbi[match(bspp_df2_ts$spp, emp_climtslo$spp_num)]
# jpeg(file = "figures/growthModelsMain/muALLbspp.jpeg",
#      width = 1800, height = 2500, res = 300)
n_spp <- length(unique(emp_climtslo$latbi))
y_pos <- rev(1:n_spp)

# set margins throught

# Row 1: GDD
par(mfrow = c(1,2))
plot(bspp_df2_ts$fit_bspp, y_pos,
     xlim = c(-10, 10), ylim = c(0.5, n_spp + 0.5), 
     xlab = "", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = tscolslatbi, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_ts$fit_bspp_per5,  y_pos, bspp_df2_ts$fit_bspp_per95, y_pos,
         col = tscolslatbi, lwd = 1.5)
segments(bspp_df2_ts$fit_bspp_per25, y_pos, bspp_df2_ts$fit_bspp_per75, y_pos,
         col = tscolslatbi, lwd = 3)
mtext("", side = 3, adj = 0, font = 2, cex = 0.9)


# Slot 5: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("right",
       legend = sapply(unique(bspp_df2_ts$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(tscolslatbi),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)

##### CoringTreespotters #####
clim_vars <- c("TempMeanMax", "TempMeanMin", "TempMeanMean")
climvar <- clim_vars[1] 
period <- "DJF"

d <- emp_climtslo[emp_climtslo$period == period & !is.na(emp_climtslo$TempMeanMax) & 
                  !is.na(emp_climtslo$anomleafout), ]

data <- list(
  y = d$anomleafout,
  N = nrow(d),
  year = as.numeric(as.character(d$year_num)),
  species = as.numeric(as.character(d$spp_num)),
  Nspp = length(unique(d$spp_num)),
  Nyear = length(unique(d$year_num)),
  climpredictor = d[[climvar]]
)
data

# Fit models
climmodelts <- stan_model("stan/TSclimatePredictors.stan")
fit <- sampling(climmodelts, data = data, iter = 2000, chains=4, cores = 4)

# Plot posterior vs priors for gdd fit 
pdf(file = "figures/climate/climateModelPriorVSPosteriorTS.pdf", 
    width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
ayear_df <- df_fit[, columns[grepl("ayear", columns)]]

# # a
# plot(density(df_fit[, "a_prior"]), 
#      col = pal[1], lwd = 2, 
#      main = "priorVSposterior_a", 
#      xlab = "a", ylim = c(0,0.5))
# lines(density(df_fit[, "a"]), col = pal[2], lwd = 2)
# legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fit[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fit[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fit[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fit[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fit[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-10, 10),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Plot spp parameter output #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
df_fit <- as.data.frame(fit)

# full posterior
columns <- colnames(df_fit)
columns <- columns[!grepl("prior", columns)]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
ayear_df <- df_fit[, grepl("year", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(ayear) <- 1:ncol(year_df)
colnames(aspp_df) <- 1:ncol(aspp_df)

# posterior summaries
sigma_df2_ts  <- extract_params(df_fit, "sigma", "mean", "sigma")
bspp_df2_ts   <- extract_params(df_fit, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
year_df2_ts <- extract_params(df_fit, "ayear", "fit_ayear", "year", "ayear\\[(\\d+)\\]")
aspp_df2_ts   <- extract_params(df_fit, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")

year_df2_ts$year_name <- emp_climtslo$year[match(year_df2_ts$year, emp_climtslo$year_num)]
aspp_df2_ts$spp_name <- emp_climtslo$latbi[match(aspp_df2_ts$spp, emp_climtslo$spp_num)]
bspp_df2_ts$spp_name <- emp_climtslo$latbi[match(bspp_df2_ts$spp, emp_climtslo$spp_num)]
# jpeg(file = "figures/growthModelsMain/muALLbspp.jpeg",
#      width = 1800, height = 2500, res = 300)
n_spp <- length(unique(emp_climtslo$latbi))
y_pos <- rev(1:n_spp)

# set margins throught

# Row 1: GDD
par(mfrow = c(1,2))
plot(bspp_df2_ts$fit_bspp, y_pos,
     xlim = c(-10, 10), ylim = c(0.5, n_spp + 0.5), 
     xlab = "", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = tscolslatbi, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_ts$fit_bspp_per5,  y_pos, bspp_df2_ts$fit_bspp_per95, y_pos,
         col = tscolslatbi, lwd = 1.5)
segments(bspp_df2_ts$fit_bspp_per25, y_pos, bspp_df2_ts$fit_bspp_per75, y_pos,
         col = tscolslatbi, lwd = 3)
mtext("", side = 3, adj = 0, font = 2, cex = 0.9)


# Slot 5: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("right",
       legend = sapply(unique(bspp_df2_ts$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(tscolslatbi),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)

# dev.off()
}

}

# Checks for temperatures and PDSI ####
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

# Late spring frosts ####
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
