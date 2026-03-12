# Climate comparison
# CRD 24 September 2025
# Weldhill climate station changed, so I want to look at the overlap

# Load library 
library(ggplot2)
# library(wesanderson)
library(tidyverse)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/input/_notcookies/")

makeplot <- FALSE

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Start with years 2012 to 2020
#from https://labs.arboretum.harvard.edu/weather/
weldhillclimate <- read.csv("weldhill.csv") 

yr15TO19 <- filter(weldhillclimate, grepl("2015|2016|2017|2018|2019",weldhillclimate$Eastern.daylight.time)) ## pull the years until 2019 --- see issue #8 for selecting until 2019 only

## get max and min T
yr15TO19 <- yr15TO19 %>% 
  separate(Eastern.daylight.time, into = c('date', 'time'), sep = -9, convert = TRUE) ### make date and time separate cols

yr15TO19$date <- as.Date(yr15TO19$date, format = "%m/%d/%Y")

# convert far to celcius: https://www.metric-conversions.org/temperature/fahrenheit-to-celsius.htm
yr15TO19$tempCelcius <- (yr15TO19$Temp..F-32)/1.8
yr15TO19$pptMM <- yr15TO19$Rain.in * 25.4

# summarize with min mean and max
yr15TO19_2 <- aggregate(
  tempCelcius ~ date,
  data = yr15TO19,
  FUN = function(x) c(max = max(x, na.rm = TRUE),
                      mean = mean(x, na.rm = TRUE),
                      min = min(x, na.rm = TRUE))
)

yr15TO19_3 <- data.frame(
  date = yr15TO19_2$date,
  maxTempC = yr15TO19_2$tempCelcius[, "max"],
  meanTempC = yr15TO19_2$tempCelcius[, "mean"],
  minTempC = yr15TO19_2$tempCelcius[, "min"]
)

# summarize for precipitation
# summarize with min mean and max
yr15TO19_ppt <- aggregate(
  pptMM ~ date,
  data = yr15TO19,
  FUN = sum)

yr15TO19_3 <- data.frame(
  date = yr15TO19_2$date,
  maxTempC = yr15TO19_2$tempCelcius[, "max"],
  meanTempC = yr15TO19_2$tempCelcius[, "mean"],
  minTempC = yr15TO19_2$tempCelcius[, "min"]
)

yr15TO19_3$pptMM <- yr15TO19_ppt$pptMM[match(yr15TO19_3$date, yr15TO19_ppt$date)]

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Then follow with years 2020 to 2024
# From https://newa.cornell.edu/all-weather-data-query/
# precipitation in inches
files <- list.files(pattern = "\\.csv$")[2:6] # 2 and 10 is wrong

# Read them all in as a list of data.frames
df_list <- lapply(files, read.csv, stringsAsFactors = FALSE)

# Bind them together (row-wise)
yr20TO25 <- do.call(rbind, df_list)

# select cols of interest
yr20TO25_2 <- yr20TO25[,c("date", 
                          "Avg.Air.Temp...F.", 
                          "Max.Air.Temp...F.", 
                          "Min.Air.Temp...F.",
                          "Total.Precipitation")]

yr20TO25_3 <- data.frame(
  date      = yr20TO25_2$date,
  maxTempC  = (yr20TO25_2$Max.Air.Temp...F. - 32) * 5/9,
  meanTempC = (yr20TO25_2$Avg.Air.Temp...F. - 32) * 5/9,
  minTempC  = (yr20TO25_2$Min.Air.Temp...F. - 32) * 5/9,
  pptMM  = yr20TO25_2$Total.Precipitation * 25.4
)

yr20TO25_3$date <- as.Date(yr20TO25_3$date, format = "%m/%d/%Y")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

# Bind both dfs together but first give them a unique id
yr15TO19_3$source <- "weldhillsite"
yr20TO25_3$source <- "newa"

weldhillcomp <- rbind(yr15TO19_3, yr20TO25_3)

# removing temp of -500
weldhillcomp <- subset(weldhillcomp, minTempC > -50) 

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")

# add years and doy
weldhillcomp$year<-substr(weldhillcomp$date, 1, 4)
weldhillcomp$doy <- yday(weldhillcomp$date)

# write csv
write.csv(weldhillcomp, "output/weldhillClimateCleaned.csv", row.names = FALSE)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if(makeplot) {
# make a quick and dirty plot to compare overlapping data for 2020
y2020 <- subset(weldhillcomp, year == 2020)
y2020 <- subset(y2020, doy < 240)

ggplot(y2020, aes(x = doy, y = meanTempC, color = source, fill = source)) +
  geom_ribbon(aes(ymin = minTempC, ymax = maxTempC), alpha = 0.2, color = NA) +
  geom_line() +
  labs(
    title = "2020 overlapping climate data",
    x = "doy",
    y = "TempC",
    color = "Dataset",
    fill = "Dataset"
  ) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2)) +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 2)) +
  theme_minimal()
# save ggplot!
ggsave("figures/climate/climateComparison2020.jpeg", width = 8, height = 5, dpi = 300)

 # do one with 1:1 line
df_wide <- reshape(y2020,
                  timevar = "source",
                  idvar = "date",
                  direction = "wide")

ggplot(df_wide, aes(x = meanTempC.weldhillsite, y = meanTempC.newa)) +
 geom_point(color = "#046C9A", size = 2) +
 geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
 labs(x = "weldhillsite mean temp", y = "newa mean temp", title = "") +
 theme_minimal()
ggsave("figures/climate/climateComparison2020_abline.jpeg", width = 8, height = 5, dpi = 300)
 
weldhillcomp2 <- subset(weldhillcomp, year %in% c(2016:2019, 2021:2024))

ggplot(weldhillcomp2, aes(x = doy, y = meanTempC, color = source, fill = source)) +
  geom_ribbon(aes(ymin = minTempC, ymax = maxTempC), alpha = 0.2, color = NA) +
  geom_smooth() +
  labs(
    title = "2016 to 2019 and 2021 to 2024 climate data overlaps",
    x = "doy",
    y = "TempC",
    color = "Dataset",
    fill = "Dataset"
  ) +
  scale_color_manual(values = wes_palette("Darjeeling1", n = 2)) +
  scale_fill_manual(values = wes_palette("Darjeeling1", n = 2)) +
  theme_minimal()
# save ggplot!
ggsave("figures/climate/climateComparisonMultipleyears.jpeg", width = 8, height = 5, dpi = 300)
}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# PDSI ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
readLines("input/_notcookies/terraclimate_42.2960N_71.1206W.csv", n = 20)
terra <- read.csv("input/_notcookies/terraclimate_42.2960N_71.1206W.csv", skip = 14)
str(terra)

# change colnames 
colnames(terra) <- c(
  "year",
  "month",
  "tmax",
  "tmin",
  "ppt",
  "soilMoisture",
  "vpd",
  "pdsi",
  "downwardRadiation"
  )
terra

terra$tmean <- (terra$tmax + terra$tmin)/2

# create a new column that attributes december as the following season winter
terra$season_year <- ifelse(terra$month == 12, terra$year + 1, terra$year)

# dec, jan, feb
terra_djf <- subset(terra, season_year > 2015 & month %in% c(12, 1, 2))
# march, apr, may
terra_mam <- subset(terra, year > 2015 & month %in% c(3, 4, 5))
# june, july, aug
terra_jja <- subset(terra, year > 2015 & month %in% c(6, 7, 8))
# sept, oct, nov
terra_son <- subset(terra, year > 2015 & month %in% c(9, 10, 11))

# pdsi
pdsidjf <- aggregate(pdsi ~ year, terra_djf, mean)
pdsimam <- aggregate(pdsi ~ year, terra_mam, mean)
pdsijja <- aggregate(pdsi ~ year, terra_jja, mean)
pdsison <- aggregate(pdsi ~ year, terra_son, mean)
colnames(pdsidjf) <- c("year", "pdsi")
colnames(pdsimam) <- c("year", "pdsi")
colnames(pdsijja) <- c("year", "pdsi")
colnames(pdsison) <- c("year", "pdsi")
pdsidjf$period <- "DJD"
pdsimam$period <- "MAM"
pdsijja$period <- "JJA"
pdsison$period <- "SON"
meanpdsi <- rbind(pdsidjf,pdsimam,pdsijja,pdsison)

# mean max
meanmaxdjf <- aggregate(tmax ~ year, terra_djf, mean)
meanmaxmam <- aggregate(tmax ~ year, terra_mam, mean)
meanmaxjja <- aggregate(tmax ~ year, terra_jja, mean)
meanmaxson <- aggregate(tmax ~ year, terra_son, mean)
colnames(meanmaxdjf) <- c("year", "tmeanmax")
colnames(meanmaxmam) <- c("year", "tmeanmax")
colnames(meanmaxjja) <- c("year", "tmeanmax")
colnames(meanmaxson) <- c("year", "tmeanmax")
meanmaxdjf$period <- "DJD"
meanmaxmam$period <- "MAM"
meanmaxjja$period <- "JJA"
meanmaxson$period <- "SON"
meanmax <- rbind(meanmaxdjf,meanmaxmam,meanmaxjja,meanmaxson)

# mean mean 
meanmeandjf <- aggregate(tmean ~ year, terra_djf, mean)
meanmeanmam <- aggregate(tmean ~ year, terra_mam, mean)
meanmeanjja <- aggregate(tmean ~ year, terra_jja, mean)
meanmeanson <- aggregate(tmean ~ year, terra_son, mean)
colnames(meanmeandjf) <- c("year", "tmeanmean")
colnames(meanmeanmam) <- c("year", "tmeanmean")
colnames(meanmeanjja) <- c("year", "tmeanmean")
colnames(meanmeanson) <- c("year", "tmeanmean")
meanmeandjf$period <- "DJD"
meanmeanmam$period <- "MAM"
meanmeanjja$period <- "JJA"
meanmeanson$period <- "SON"
meanmean <- rbind(meanmeandjf,meanmeanmam,meanmeanjja,meanmeanson)

# mean min
meanmindjf <- aggregate(tmin ~ year, terra_djf, mean)
meanminmam <- aggregate(tmin ~ year, terra_mam, mean)
meanminjja <- aggregate(tmin ~ year, terra_jja, mean)
meanminson <- aggregate(tmin ~ year, terra_son, mean)
colnames(meanmindjf) <- c("year", "tmeanmin")
colnames(meanminmam) <- c("year", "tmeanmin")
colnames(meanminjja) <- c("year", "tmeanmin")
colnames(meanminson) <- c("year", "tmeanmin")
meanmindjf$period <- "DJD"
meanminmam$period <- "MAM"
meanminjja$period <- "JJA"
meanminson$period <- "SON"
meanmin <- rbind(meanmindjf,meanminmam,meanminjja,meanminson)

# precipitation
meanpptdjf <- aggregate(ppt ~ year, terra_djf, mean)
meanpptmam <- aggregate(ppt ~ year, terra_mam, mean)
meanpptjja <- aggregate(ppt ~ year, terra_jja, mean)
meanpptson <- aggregate(ppt ~ year, terra_son, mean)
colnames(meanpptdjf) <- c("year", "ppt")
colnames(meanpptmam) <- c("year", "ppt")
colnames(meanpptjja) <- c("year", "ppt")
colnames(meanpptson) <- c("year", "ppt")
meanpptdjf$period <- "DJD"
meanpptmam$period <- "MAM"
meanpptjja$period <- "JJA"
meanpptson$period <- "SON"
meanppt <- rbind(meanpptdjf,meanpptmam,meanpptjja,meanpptson)

# radiation
meanraddjf <- aggregate(downwardRadiation ~ year, terra_djf, mean)
meanradmam <- aggregate(downwardRadiation ~ year, terra_mam, mean)
meanradjja <- aggregate(downwardRadiation ~ year, terra_jja, mean)
meanradson <- aggregate(downwardRadiation ~ year, terra_son, mean)
colnames(meanraddjf) <- c("year", "rad")
colnames(meanradmam) <- c("year", "rad")
colnames(meanradjja) <- c("year", "rad")
colnames(meanradson) <- c("year", "rad")
meanraddjf$period <- "DJD"
meanradmam$period <- "MAM"
meanradjja$period <- "JJA"
meanradson$period <- "SON"
meanrad <- rbind(meanraddjf,meanradmam,meanradjja,meanradson)

merged <- Reduce(function(x, y) merge(x, y, by = c("year", "period")),
                 list(meanpdsi, meanmax, meanmean, meanmin,
                      meanppt, meanrad))
merged <- subset(merged, year != 2015)

write_csv(merged, "output/climateSummariesYear.csv")
