## Started 29 May 2022 ##
## Day 2 of transit back to Vancouer ## 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# setwd and libraries
setwd("~/Documents/git/projects/treegarden/wildhellgarden/analyses")
library(ggplot2)

# get heights data, where we mentioned dead and alive
d <- read.csv("input/heights2021.csv", header=TRUE)
unique(d$Notes)
d <- subset(d, ID!="")

# make a column of dead versus alive
d$status <- 1
d$status[grep("Dead", d$Notes)] <- 0
d$status[grep("dead", d$Notes)] <- 0 # this grabs 'almost dead' and 'probably dead'

# now extract species etc. 
d$species <- unlist(lapply(strsplit(d$ID, "_"), "[", 1))
d$site <- unlist(lapply(strsplit(d$ID, "_"), "[", 2))
d$num <- unlist(lapply(strsplit(d$ID, "_"), "[", 3))

sort(unique(d$species))
unique(d$site)

# Hacky way to make some groupings...
# first, by plot
dplotalive <- aggregate(d["status"], d[c("Plot")], FUN=sum)
dplotall <- aggregate(d["status"], d[c("Plot")], FUN=length)
dplot <- merge(dplotalive, dplotall, by=c("Plot"), suffixes=c("alive", "total"), all.y=TRUE)
dplot$percentalive <- dplot$statusalive/dplot$statustotal # 11-12 seem fine
dplot[order(dplot$percentalive),] # I suspect differences here are related to species in plots


# Okay, how about site?
dsitealive <- aggregate(d["status"], d[c("site")], FUN=sum)
dsiteall <- aggregate(d["status"], d[c("site")], FUN=length)
dsite <- merge(dsitealive, dsiteall, by=c("site"), suffixes=c("alive", "total"), all.y=TRUE)
dsite$percentalive <- dsite$statusalive/dsite$statustotal
dsite[order(dsite$percentalive),]
# sites are similar!

# by species and site
dalive <- aggregate(d["status"], d[c("species", "site")], FUN=sum)
names(dalive)[names(dalive)=="status"] <- "alive"
dall <- aggregate(d["status"], d[c("species", "site")], FUN=length)
names(dall)[names(dall)=="status"] <- "total"
dsum <- merge(dalive, dall, by=c("species", "site"), all.y=TRUE)

dsum$percentalive <- dsum$alive/dsum$total
dsumsp <- aggregate(dsum["percentalive"], dsum["species"], FUN=mean)
dtotsp <- aggregate(dsum["total"], dsum["species"], FUN=mean)
dsumsp[order(dsumsp$percentalive),]
# So definitely some species differences...

# So let's do species alone... 
dspeciesalive <- aggregate(d["status"], d[c("species")], FUN=sum)
dspeciesall <- aggregate(d["status"], d[c("species")], FUN=length)
dspecies <- merge(dspeciesalive, dspeciesall, by=c("species"), suffixes=c("alive", "total"), all.y=TRUE)
dspecies$percentalive <- dspecies$statusalive/dspecies$statustotal 

dspecies[order(dspecies$percentalive),]

# I never figured a good way to plot ... 
ggplot(d, aes(y=status)) +
    geom_bar() +
    facet_grid(species~site)
