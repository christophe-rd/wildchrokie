# code started November 18, 2024 by Deirdre
# aim: to use tmin data from daymet to calculate frost free days for our four common garden sites
# note: daymet data downloaded in check_daymet.R
require(sharpshootR)
require(ggplot2)

rm(list = ls())
options(stringsAsFactors = FALSE)

#setwd
if(length(grep("deirdre", getwd())) > 0){
  setwd("~/Documents/github/wildhellgarden")
}else {   setwd("~/Documents/github/wildhellgarden")}

d <- read.csv("analyses/daymetData/dayMet.csv")

d$datetime <- as.Date(sprintf('%s %s', d$year, d$doy), format="%Y %j")
d$value <- d$tmin

d <- subset(d, year < 2022)
coords <- c("Harvard forest", "White Mountains", "Grant", "St. Hippolyte")
timeSeries <- unique(d$year)
frostFree <- data.frame(site = rep(sites, each = length(timeSeries)), year = rep(timeSeries, times = 4))
  
ffdData <-  vector()

for(s in 1: length(coords)){

  temp <- subset(d, site == coords[s])
  
  for(y in 1: length(unique(d$year))){
    tempY <- subset(temp, year == timeSeries[y])
    
    ffd <- FFD(tempY, returnDailyPr = TRUE, frostTemp = 0)
    
    ffdData <- rbind(ffdData, ffd$summary)
    
  }
}

f <- cbind(frostFree, ffdData)

write.csv(f, "analyses/daymetData/frostFreeDays.csv", row.names = F)

siteOrder <- c("Harvard forest", "White Mountain","Grant", "St. Hippolyte" )

pdf("analyses/figures/ffdBySite.pdf", width =6, height = 5)
ggplot(f, aes(x = factor(site, levels = siteOrder), y = ffd.90, fill= factor(site, siteOrder))) +
        geom_boxplot() + geom_point() +
        labs ( x = "Site", y = "FFD", fill = "Site") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) 
dev.off()

spring <- ggplot(f, aes(x = factor(site, levels = siteOrder), y = spring.90)) +
  geom_boxplot() + geom_point() +
  labs ( x = "Site", y = "Spring FFD") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) 

fall <- ggplot(f, aes(x = factor(site, levels = siteOrder), y = fall.90)) +
  geom_boxplot() + geom_point() +
  labs ( x = "Site", y = "Fall FFD") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) 

pdf("analyses/figures/ffdOverTime.pdf", width =12, height = 5)
ggplot(f, aes(x = year, y = ffd.90, fill= factor(site, siteOrder))) + 
  geom_point() +
  geom_smooth(method='lm', col = "black") +
  labs ( x = "Year", y = "FFD", fill = "Site") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  facet_wrap(~ factor(site, siteOrder), ncol = 4)
dev.off()





