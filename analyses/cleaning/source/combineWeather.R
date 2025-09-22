# started November 11th by Deirdre
# I guess kinda updated by CRD for wildchrokie on 23 April 2025
# this source script cleans and merges the relevant weather data:
library(pollen)
library(lubridate)
library(tidyverse)


weldhillclimate <- read.csv("input/_notcookies/weldhill.csv") #from https://labs.arboretum.harvard.edu/weather/

weldhillclimate <- filter(weldhillclimate,grepl("2018|2019|2020",weather$Eastern.daylight.time)) ## pull the years we need

##get max and min T
weldhillclimate2 <- weldhillclimate %>% 
  separate(Eastern.daylight.time, into = c('date', 'time'), sep = -8, convert = TRUE) ### make date and time separate cols

# convert far to celcius: https://www.metric-conversions.org/temperature/fahrenheit-to-celsius.htm
weldhillclimate2$tempCelcius <- (weldhillclimate2$Temp..F-32)/1.8

# summarize with min mean and max
weldhillclimate3 <- weldhillclimate2 %>% 
  group_by(date) %>% 
  summarise(maxT=max(tempCelcius ),
            meanT=mean(tempCelcius ), 
            minT=min(tempCelcius ))

weldhillclimate3$date <- as.Date(weldhillclimate3$date, format = "%d/%m/%Y") # make dates dates again

weldhillclimate3$year<-substr(weldhillclimate3$date, 1, 4)
weldhillclimate3$doy <- yday(weldhillclimate3$date)

###give each year a seperate data frame
y18<-filter(weldhillclimate3,year=="2018")
y19<-filter(weldhillclimate3,year=="2019")
y20<-filter(weldhillclimate3,year=="2020")

###Calculate gdd
y18$GDD_10 <- gdd(tmax = y18$maxT, tmin = y18$minT, tbase = 10, type = "B")
y19$GDD_10 <- gdd(tmax = y19$maxT, tmin = y19$minT, tbase = 10, type = "B")
y20$GDD_10 <- gdd(tmax = y20$maxT, tmin = y20$minT, tbase = 10, type = "B")
 
allyr <- rbind(y18, y19, y20)

# rename
gdd <- allyr
write.csv(allyr, "output/gddData.csv", row.names = F)

