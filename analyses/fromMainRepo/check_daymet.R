## Started 17 November 2024
## By Ken
## annotated by Deirdre Nov 18 2024

## check climate data

rm(list = ls())
options(stringsAsFactors = FALSE)

#setwd
if(length(grep("Ken", getwd())) > 0){
  setwd("/Users/Ken Michiko Samson/Documents/Temporal Ecology Lab/misc")
}else {   setwd("~/Documents/github/wildhellgarden")}

library(daymetr)
library(sf)
library(terra)

tmin <- rast("2129_daymet_v4_daily_na_tmin_1980.nc4")

coords <- data.frame(lat = c(42.5315, 45.9310, 44.92466697, 43.99837498), long = c(-72.1900, -74.0248, -71.09543597, -71.37388296))
coords <- st_as_sf(coords, coords = c("long", "lat"), crs = 4326)
coords <- st_transform(coords, crs = crs(tmin))

data <- extract(tmin, vect(coords), xy = TRUE, ID = TRUE, time = TRUE, na.rm = FALSE)

daymet <- download_daymet("daymet", 42.5315, -72.1900, 1980, 1980)
sample <- subset(daymet$data, select = c("year", "yday", "tmin..deg.c."))

data[,1:6]
head(sample)
# data from daymetr and manual subsetting seem to agree

rm(list = ls())

coords <- data.frame(ID = seq(1,4), 
                     lat = c(42.5315, 45.9310, 44.92466697, 43.99837498),
                     long = c(-72.1900, -74.0248, -71.09543597, -71.37388296))
leap <- seq(1980, 2020, 4)

data <- list()
for(n in 1:nrow(coords)){
  daymet <- download_daymet("daymet", coords$lat[n], coords$long[n], 1980, 2023)
  temp <- subset(daymet$data, select = c("year", "yday", "tmax..deg.c.", "tmin..deg.c."))
  
  for(y in leap){
    tmin_366 <- mean(temp$tmax..deg.c.[temp$year == y & temp$yday == 365],
                     temp$tmax..deg.c.[temp$year == (y + 1) & temp$yday == 1])
    tmax_366 <- mean(temp$tmin..deg.c.[temp$year == y & temp$yday == 365],
                     temp$tmin..deg.c.[temp$year == (y + 1) & temp$yday == 1])
    
    yday_366 <- data.frame(year = y, yday = 366,
                           tmax..deg.c. = tmax_366, tmin..deg.c. = tmin_366)
    
    temp <- rbind(temp, yday_366)
  }
  
  colnames(temp) <- c("year", "doy", "tmax", "tmin")
  temp$siteID <- coords$ID[n]
  temp <- temp[,c("siteID", "year", "doy", "tmax", "tmin")]
  
  data[[n]] <- temp
}
data <- do.call(rbind, data)
data <- data[order(data$siteID, data$year, data$doy),]

write.csv(coords, "sites.csv", row.names = FALSE)
write.csv(data, "temps.csv", row.names = FALSE)

colnames(temps)[colnames(temps) == "siteID"] <- "ID"

dayMet <- merge(coords, temps, by = "ID")

write.csv(dayMet, "analyses/daymetData/dayMet.csv", row.names = F)
