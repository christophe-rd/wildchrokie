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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Then follow with years 2020 to 2024
# From https://newa.cornell.edu/all-weather-data-query/
files <- list.files(pattern = "\\.csv$")[2:6] # 2 and 10 is wrong

# Read them all in as a list of data.frames
df_list <- lapply(files, read.csv, stringsAsFactors = FALSE)

# Bind them together (row-wise)
yr20TO25 <- do.call(rbind, df_list)

# select cols of interest
yr20TO25_2 <- yr20TO25[,c("date", 
                          "Avg.Air.Temp...F.", 
                          "Max.Air.Temp...F.", 
                          "Min.Air.Temp...F.")]

yr20TO25_3 <- data.frame(
  date      = yr20TO25_2$date,
  maxTempC  = (yr20TO25_2$Max.Air.Temp...F. - 32) * 5/9,
  meanTempC = (yr20TO25_2$Avg.Air.Temp...F. - 32) * 5/9,
  minTempC  = (yr20TO25_2$Min.Air.Temp...F. - 32) * 5/9
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
 