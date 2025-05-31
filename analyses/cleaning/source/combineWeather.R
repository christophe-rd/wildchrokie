# started November 11th by Deirdre
# I guess kinda updated by CRD for wildchrokie on 23 April 2025
# this source script cleans and merges the relevant weather data:
library("pollen")
library(lubridate)
library(tidyverse)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")

weldhill <- FALSE
if(weldhill) {

y<-read.csv("input/_notcookies/weldhill.csv") #from https://labs.arboretum.harvard.edu/weather/

weather<-filter(weather,grepl("2015|2016|2017|2018|2019|2020|2021|2022|2023",weather$Eastern.daylight.time)) ## pull the years we need

##get max and min T
goo<-weather %>% separate(Eastern.daylight.time, into = c('date', 'time'), sep = -8, convert = TRUE) ### make date and time separate cols

goo<-goo%>% group_by(date) %>% summarise(maxTf=max(Temp..F ),minTf=min(Temp..F )) ###make daily mins and maxT

}

### whoops, The weld hill temperature data ends, use the closest weather station
bostonairport<-read.csv("input/_notcookies/bostonAirportClimate.csv")

###get ready to mush them together
str(bostonairport)
bostonairport$date<-as.Date(bostonairport$DATE) # make dates dates again

bostonairport$doy<-yday(bostonairport$date) ##convert to doy

# grab columns of interest
bostonairportsub <- bostonairport[, c("date", "doy", "TAVG", "TMAX", "TMIN")]
# rename to something less CRAZY
colnames(bostonairportsub) <- c("date", "doy", "aveT", "maxT", "minT")

#new column for year
bostonairportsub$year<-substr(bostonairportsub$date, 1, 4)
###give each year a seperate data frame
y15<-filter(bostonairportsub,year=="2015")
y16<-filter(bostonairportsub,year=="2016")
y17<-filter(bostonairportsub,year=="2017")
y18<-filter(bostonairportsub,year=="2018")
y19<-filter(bostonairportsub,year=="2019")
y20<-filter(bostonairportsub,year=="2020")
y21<-filter(bostonairportsub,year=="2021")
y22<-filter(bostonairportsub,year=="2022")
y23<-filter(bostonairportsub,year=="2023")
y24<-filter(bostonairportsub,year=="2024")
###Calculate gdd
y15$GDD_10<- pollen::gdd(tmax = y15$maxT,tmin = y15$minT,tbase = 10,type = "B")
y16$GDD_10<- pollen::gdd(tmax = y16$maxT,tmin = y16$minT,tbase = 10,type = "B")
y17$GDD_10<- pollen::gdd(tmax = y17$maxT,tmin = y17$minT,tbase = 10,type = "B")
y18$GDD_10<- pollen::gdd(tmax = y18$maxT,tmin = y18$minT,tbase = 10,type = "B")
y19$GDD_10<- pollen::gdd(tmax = y19$maxT,tmin = y19$minT,tbase = 10,type = "B")
y20$GDD_10<- pollen::gdd(tmax = y20$maxT,tmin = y20$minT,tbase = 10,type = "B")
y21$GDD_10<- pollen::gdd(tmax = y21$maxT,tmin = y21$minT,tbase = 10,type = "B")
y22$GDD_10<- pollen::gdd(tmax = y22$maxT,tmin = y22$minT,tbase = 10,type = "B")
y23$GDD_10<- pollen::gdd(tmax = y23$maxT,tmin = y23$minT,tbase = 10,type = "B")
y24$GDD_10<- pollen::gdd(tmax = y24$maxT,tmin = y24$minT,tbase = 10,type = "B")


allyr<-rbind(y15, y16, y17, y18, y19, y20, y21, y22, y23, y24)

write.csv(allyr, "output/gddData.csv", row.names = F)

