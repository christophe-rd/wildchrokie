## From main wildhell repo
# CRD on 22 Sept 2025

# Goal is calculate primary and full growing season

# Clear workspace 
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")

# Load libraries
library(pollen)

# copy of obsdata2
obsdata2 <- obsdata

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

write_csv(gdd, "output/gddByYear.csv")

# replace NaN by NA. Hopefully this is ok
obsdata2$budburst <- gsub("NaN", "NA", obsdata2$budburst) 
obsdata2$budset <- gsub("NaN", "NA", obsdata2$budset)
obsdata2$leafout <- gsub("NaN", "NA", obsdata2$leafout)
obsdata2$leafcolor <- gsub("NaN", "NA", obsdata2$leafcolor)

# convert to integer to avoid downstream problems with matching gdd
obsdata2$budburst2 <- round(obsdata2$budburst)
obsdata2$budset2 <- round(obsdata$budset)

unique(as.numeric(obsdata2$leafout))

obsdata2$leafout2 <- obsdata2$leafout
obsdata2$leafout2 <- gsub("135.5", 136, obsdata2$leafout2)
obsdata2$leafout2 <- gsub("130.5", 131, obsdata2$leafout2)
obsdata2$leafout2 <- gsub("132.5", 133, obsdata2$leafout2)
obsdata2$leafout2 <- gsub("134.5", 135, obsdata2$leafout2)
obsdata2$leafout2 <- gsub("131.5", 132, obsdata2$leafout2)

unique(obsdata2$leafout2)
obsdata2$leafout3 <- as.numeric((obsdata2$leafout2))

nrow(obsdata2[!is.na(obsdata2$leafout),])
nrow(obsdata2[!is.na(obsdata2$leafout2),])

obsdata2$leafcolor2 <- as.integer(obsdata$leafcolor)


unique(obsdata2$budburst)
unique(obsdata2$budset)
unique(obsdata2$leafout)
unique(obsdata2$leafcolor)

# create new GDD empty columns
obsdata2$budburstGDD <- NA
obsdata2$budsetGDD <- NA
obsdata2$leafoutGDD <- NA
obsdata2$leafcolorGDD <- NA

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

# add primary GS and full GS cols
obsdata2$pgsGDD <- obsdata2$budsetGDD - obsdata2$leafoutGDD
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
