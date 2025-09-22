## From main wildhell repo
# CRD on 22 Sept 2025

# Goal is calculate primary and full growing season

# Clear workspace 
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")
# Load libraries

###calculate primary growth season and full growing season in days
obsdata$pgs <- obsdata$budset-obsdata$leafout
obsdata$fgs <- obsdata$leafcolor-obsdata$leafout

# Initialize new GDD columns
obsdata$budburstGDD <- NA
obsdata$budsetGDD <- NA
obsdata$leafoutGDD <- NA
obsdata$leafcolorGDD <- NA

# Loop over rows
for(i in 1:nrow(obsdata)) {
  yr <- obsdata$year[i]
  
  # Budburst
  if(!is.na(obsdata$budburst[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata$budburst[i])
    if(length(idx) == 1) obsdata$budburstGDD[i] <- gdd$GDD_10[idx]
  }
  
  # Budset
  if(!is.na(obsdata$budset[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata$budset[i])
    if(length(idx) == 1) obsdata$budsetGDD[i] <- gdd$GDD_10[idx]
  }
  
  # Leafout
  if(!is.na(obsdata$leafout[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata$leafout[i])
    if(length(idx) == 1) obsdata$leafoutGDD[i] <- gdd$GDD_10[idx]
  }
  
  # Leafcolor
  if(!is.na(obsdata$leafcolor[i])) {
    idx <- which(gdd$year == yr & gdd$doy == obsdata$leafcolor[i])
    if(length(idx) == 1) obsdata$leafcolorGDD[i] <- gdd$GDD_10[idx]
  }
}

# add primary GS and full GS cols
obsdata$pgsGDD <- obsdata$budsetGDD - obsdata$leafoutGDD
obsdata$fgsGDD <- obsdata$leafcolorGDD - obsdata$leafoutGDD

obsdataWithGDD <- obsdata

pgslookup <- obsdataWithGDD[!is.na(obsdataWithGDD$pgsGDD),]
nrow(pgslookup)
