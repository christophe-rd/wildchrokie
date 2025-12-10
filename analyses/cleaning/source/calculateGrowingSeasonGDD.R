## From main wildhell repo
# CRD on 22 Sept 2025

# Goal is calculate primary and full growing season

# Clear workspace 
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")

# Load libraries
library(pollen)

# copy of obsdata2
obsdata2 <- obsdata

obsdata2$leafout <- round(obsdata2$leafout)

###calculate primary growth season and full growing season in days
obsdata2$pgs <- obsdata2$budset-obsdata2$leafout
obsdata2$fgs <- obsdata2$leafcolor-obsdata2$leafout

# calculate GDD
###give each year a seperate data frame
y18 <- subset(weldhillcomp, year=="2018")
y19 <- subset(weldhillcomp, year=="2019")
y20 <- subset(weldhillcomp, year=="2020")

###Calculate gdd
y18$GDD_10 <- gdd(tmax = y18$maxT, tmin = y18$minT, tbase = 10, type = "B")
y19$GDD_10 <- gdd(tmax = y19$maxT, tmin = y19$minT, tbase = 10, type = "B")
y20$GDD_10 <- gdd(tmax = y20$maxT, tmin = y20$minT, tbase = 10, type = "B")

gdd <- rbind(y18, y19, y20)
str(gdd)

write_csv(gdd, "output/gddByYear.csv")

# replace NaN by NA. Hopefully this is ok
str(obsdata2)
for (i in 2:10) {
  obsdata2[[i]][is.nan(obsdata2[[i]])] <- NA
}
str(obsdata2)
nrow(obsdata)
nrow(obsdata2[!is.na(obsdata2$leafout),]) # missing 9 rows of NaN

unique(obsdata2$budburst)
unique(obsdata2$budset)
unique(obsdata2$leafout)
unique(obsdata2$leafcolor)

# create new GDD empty columns
obsdata2$budburstGDD <- NA
obsdata2$budsetGDD <- NA
obsdata2$leafoutGDD <- NA
obsdata2$leafcolorGDD <- NA


gdd$year <- as.numeric(gdd$year)
# Loop over rows
for(i in 1:nrow(obsdata2)) {
  yr <- obsdata2$year[i]
  
  # Budburst
  if(!is.na(obsdata2$budburst[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata2$budburst[i])
    if(length(idx) == 1) obsdata2$budburstGDD[i] <- gdd$GDD_10[idx]
  }
  
  # Budset
  if(!is.na(obsdata2$budset[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata2$budset[i])
    if(length(idx) == 1) obsdata2$budsetGDD[i] <- gdd$GDD_10[idx]
  }
  
  # Leafout
  if(!is.na(obsdata2$leafout[i])) { # 49
    idx <- which(gdd$year == yr & gdd$doy == obsdata2$leafout[i])
    if(length(idx) == 1) obsdata2$leafoutGDD[i] <- gdd$GDD_10[idx]
  }
  
  # Leafcolor
  if(!is.na(obsdata2$leafcolor[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata2$leafcolor[i])
    if(length(idx) == 1) obsdata2$leafcolorGDD[i] <- gdd$GDD_10[idx]
  }
  
}
str(obsdata2)

# check that the loop didn't miss any rows
nrow(obsdata2[!is.na(obsdata2$budburst),]) 
nrow(obsdata2[!is.na(obsdata2$budburstGDD),])

nrow(obsdata2[!is.na(obsdata2$budset),]) 
nrow(obsdata2[!is.na(obsdata2$budsetGDD),])

nrow(obsdata2[!is.na(obsdata2$leafout),]) 
nrow(obsdata2[!is.na(obsdata2$leafoutGDD),])

nrow(obsdata2[!is.na(obsdata2$leafcolor),]) 
nrow(obsdata2[!is.na(obsdata2$leafcolorGDD),])

# add primary GS and full GS cols
obsdata2$pgsGDD <- obsdata2$budsetGDD - obsdata2$leafoutGDD
nrow(obsdata2[!is.na(obsdata2$pgsGDD),]) 
obsdata2$fgsGDD <- obsdata2$leafcolorGDD - obsdata2$leafoutGDD

# add max gdd per year
obsdata2$fullGDD <- NA
y18maxGDD <- max(y18$GDD_10)
y19maxGDD <- max(y19$GDD_10)
y20maxGDD <- max(y20$GDD_10)
obsdata2$fullGDD[which(obsdata2$year == "2018")] <- y18maxGDD
obsdata2$fullGDD[which(obsdata2$year == "2019")] <- y19maxGDD
obsdata2$fullGDD[which(obsdata2$year == "2020")] <- y20maxGDD

obsdataWithGDD <- obsdata2

pgslookup <- obsdataWithGDD[!is.na(obsdataWithGDD$pgsGDD),]
nrow(pgslookup)
