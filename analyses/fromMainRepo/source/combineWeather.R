# started November 11th by Deirdre
# this source script cleans and merges the relevant weather data:
library("pollen")
library(lubridate)
weather<-read.csv("input/weldhill.csv") #from https://labs.arboretum.harvard.edu/weather/

weather<-filter(weather,grepl("2018|2019|2020",weather$Eastern.daylight.time)) ## pull the years we need

##get max and min T
goo<-weather %>% separate(Eastern.daylight.time, into = c('date', 'time'), sep = -8, convert = TRUE) ### make date and time separate cols

goo<-goo%>% group_by(date) %>% summarise(maxTf=max(Temp..F ),minTf=min(Temp..F )) ###make daily mins and maxT

###whoops, The weld hill temperature data ends, use the closest weather station
more20<-read.csv("input/Daily Data for Boston (Weld Hill), MA - [id=ma_weld nwon, lat=42.2953, lon=-71.1337].csv")

###get ready to mush them together
goo2<-goo # make new copy df
goo2$date<-as.Date(goo2$date,"%m/%d/%Y") # make dates dates again
more20$date<-as.Date(more20$date,"%m/%d/%Y")

goo2$doy<-yday(goo2$date) ##convert to doy
more20$doy<-yday(more20$date) ##convert to doy
###put it in order
goo2<-goo2[order(goo2$date), ]
more20<-more20[order(more20$date), ]
colnames(more20)
more20<-dplyr::select(more20,date,Max.Air.Temp...F.,Min.Air.Temp...F.,doy)
colnames(goo2)
colnames(more20)<-colnames(goo2)
goo2<-rbind(goo2,more20)

goo2$maxTc<-(goo2$maxTf-32)*(5/9) ##convert to  C
goo2$minTc<-(goo2$minTf-32)*(5/9)

#new column for year
goo2$year<-substr(goo2$date, 1, 4)
###give each year a seperate data frame
goo18<-filter(goo2,year=="2018")
goo19<-filter(goo2,year=="2019")
goo20<-filter(goo2,year=="2020")

###Calculate gdd
goo18$GDD_10<- pollen::gdd(tmax = goo18$maxTc,tmin = goo18$minTc,tbase = 10,type = "B")
goo19$GDD_10<- pollen::gdd(tmax = goo19$maxTc,tmin = goo19$minTc,tbase = 10,type = "B")
goo20$GDD_10<- pollen::gdd(tmax = goo20$maxTc,tmin = goo20$minTc,tbase = 10,type = "B")
goober<-rbind(goo18,goo19,goo20)

#write.csv(goober, "output/gddData.csv", row.names = F)
