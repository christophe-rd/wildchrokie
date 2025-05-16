## Integrating climate data into CG data
# Cat - 11 Feb 2021

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(egg)
library(rstan)

# Set Working Directory
setwd("~/Documents/git/wildhellgarden/analyses")

cols <-viridis_pal(option="viridis")(3)


### Let's add in Climate data now
clim <- read.csv("fromMainRepo/input/clean_addinclimate.csv", header=TRUE)
clim <- clim[(clim$climatetype=="weldhill"),]
clim <- clim[(clim$year>=2017),]
clim <- clim[!duplicated(clim),]

spring <- clim[(clim$doy>=60 & clim$doy<=181),]
spring <- spring[(spring$year>=2018),]

winter <- clim[(clim$doy>=244 | clim$doy<=60),]
winter <- winter[!(winter$doy<=60 & winter$year==2017),]
winter$chillyr <- ifelse(winter$doy>=244, (winter$year + 1), winter$year)

climate <- ggplot(spring, aes(x=doy, y=tmean, col=as.factor(year))) + #geom_point(aes(col=as.factor(year)), alpha=0.1) +
  geom_smooth(aes(col=as.factor(year), fill=as.factor(year)), stat="smooth", method="loess", se=TRUE, span=0.9) + 
  scale_color_manual(name = "Year", values=cols, labels = c("2018", "2019", "2020")) +
  scale_fill_manual(name = "Year", values=cols, labels = c("2018", "2019", "2020")) +
  theme_classic() + xlab("Day of Year") + ylab("Mean \n Temperature (°C)") +
  coord_cartesian(ylim=c(-8, 18), expand=0) + scale_x_continuous(breaks = seq(min(60), max(140), by=15)) +
  scale_y_continuous(breaks=seq(min(-8), max(18), by=4)) + theme(panel.spacing = unit(c(0,0,5,5),"cm"),
                                                                 legend.text = element_text(size=7),
                                                                 legend.title = element_text(size=8),
                                                                 legend.key.size = unit(0.8,"line"))

#quartz()
#climate

winter$doy2 <- ifelse(winter$doy<=60, winter$doy+365, winter$doy)
winclim <- ggplot(winter, aes(x=doy2, y=tmin, col=as.factor(chillyr))) + #geom_point(aes(col=as.factor(year)), alpha=0.1) +
  geom_smooth(aes(col=as.factor(chillyr), fill=as.factor(chillyr)), stat="smooth", method="loess", se=TRUE, span=0.9) + 
  scale_color_manual(name = "Year", values=cols, labels = c("2018", "2019", "2020")) +
  scale_fill_manual(name = "Year", values=cols, labels = c("2018", "2019", "2020")) +
  theme_classic() + xlab("Day of Year") + ylab("Mean \n Temperature (°C)") +
  coord_cartesian(ylim=c(-8, 24), expand=0) + scale_x_continuous(breaks = seq(min(245), max(425), by=30),
                                                                 labels = c(245, 275, 305, 335, 365, 30, 60))+
  scale_y_continuous(breaks=seq(min(-8), max(24), by=4)) + theme(panel.spacing = unit(c(0,0,5,5),"cm"),
                                                                 legend.text = element_text(size=7),
                                                                 legend.title = element_text(size=8),
                                                                 legend.key.size = unit(0.8,"line"))

#quartz()
#winclim

clim2018 <- ggplot(spring[(spring$year==2018),], aes(x=doy, y=tmean), col=cols[1]) + #geom_point(aes(col=as.factor(year)), alpha=0.1) +
  geom_smooth(col=cols[1], fill=cols[1], stat="smooth", method="loess", se=TRUE, span=0.9) + 
  #scale_color_manual(name = "Year", values=cols, labels = c("2016", "2017", "2018", "2019", "2020")) +
  #scale_fill_manual(name = "Year", values=cols, labels = c("2016", "2017", "2018", "2019", "2020")) +
  theme_classic() + xlab("Day of Year") + ylab("Mean \n Temperature (°C)") +
  coord_cartesian(ylim=c(-8, 18), expand=0) + scale_x_continuous(breaks = seq(min(0), max(140), by=30)) +
  scale_y_continuous(breaks=seq(min(-8), max(18), by=4)) + theme(panel.spacing = unit(c(0,0,5,5),"cm"),
                                                                 legend.text = element_text(size=7),
                                                                 legend.title = element_text(size=8),
                                                                 legend.key.size = unit(0.8,"line"),
                                                                 legend.position = "none")

clim2019 <- ggplot(spring[(spring$year==2019),], aes(x=doy, y=tmean), col=cols[2]) + #geom_point(aes(col=as.factor(year)), alpha=0.1) +
  geom_smooth(col=cols[2], fill=cols[2], stat="smooth", method="loess", se=TRUE, span=0.9) + 
  #scale_color_manual(name = "Year", values=cols, labels = c("2016", "2017", "2018", "2019", "2020")) +
  #scale_fill_manual(name = "Year", values=cols, labels = c("2016", "2017", "2018", "2019", "2020")) +
  theme_classic() + xlab("Day of Year") + ylab("Mean \n Temperature (°C)") +
  coord_cartesian(ylim=c(-8, 18), expand=0) + scale_x_continuous(breaks = seq(min(0), max(140), by=30)) +
  scale_y_continuous(breaks=seq(min(-8), max(18), by=4)) + theme(panel.spacing = unit(c(0,0,5,5),"cm"),
                                                                 legend.text = element_text(size=7),
                                                                 legend.title = element_text(size=8),
                                                                 legend.key.size = unit(0.8,"line"),
                                                                 legend.position = "none")

clim2020 <- ggplot(spring[(spring$year==2020),], aes(x=doy, y=tmean), col=cols[3]) + #geom_point(aes(col=as.factor(year)), alpha=0.1) +
  geom_smooth(col=cols[3], fill=cols[3], stat="smooth", method="loess", se=TRUE, span=0.9) + 
  #scale_color_manual(name = "Year", values=cols, labels = c("2016", "2017", "2018", "2019", "2020")) +
  #scale_fill_manual(name = "Year", values=cols, labels = c("2016", "2017", "2018", "2019", "2020")) +
  theme_classic() + xlab("Day of Year") + ylab("Mean \n Temperature (°C)") +
  coord_cartesian(ylim=c(-8, 18), expand=0) + scale_x_continuous(breaks = seq(min(0), max(140), by=30)) +
  scale_y_continuous(breaks=seq(min(-8), max(18), by=4)) + theme(panel.spacing = unit(c(0,0,5,5),"cm"),
                                                                 legend.text = element_text(size=7),
                                                                 legend.title = element_text(size=8),
                                                                 legend.key.size = unit(0.8,"line"),
                                                                 legend.position = "none")

#quartz()
#climall <- ggarrange(clim2018, clim2019, clim2020,  ncol=3)

winter$mwt <- ave(winter$tmin, winter$chillyr)
spring$mst <- ave(spring$tmean, spring$year)

###### Now let's add in the clean bbch data
bbch <- read.csv("output/clean_obs_allyrs.csv")

bbch$dvr <- bbch$leafout - bbch$budburst
bbch$flowering <- bbch$flowers - bbch$flobudburst
bbch$fruiting <- bbch$ripefruit - bbch$fruit
bbch$allflofruit <- bbch$ripefruit - bbch$flobudburst

bbch$fls <- bbch$flowers - bbch$leafout

mst <- subset(spring, select=c("year", "mst"))
mst <- mst[!duplicated(mst),]

mwt <- subset(winter, select=c("year", "mwt"))
mwt <- mwt[!duplicated(mwt),]

bbchmst <- full_join(bbch, mst)
bbchmt <- full_join(bbchmst, mwt)

if(FALSE){
dvr.stan <- subset(bbchmt, select=c("spp", "year", "site", "ind", "plot", "dvr", "mst"))
dvr.stan <- dvr.stan[complete.cases(dvr.stan),]


getpop <- paste(dvr.stan$spp, dvr.stan$site)
dvr.stan$pophere <- as.numeric(as.factor(getpop))
dvr.stan$spp <- as.numeric(as.factor(dvr.stan$spp))

datalist.dvr.pop <- with(dvr.stan, 
                        list(y = dvr,  
                             mst = mst,
                             #mwt = mwt,
                             sp = spp,
                             pop = pophere,
                             N = nrow(dvr.stan),
                             n_sp = length(unique(dvr.stan$spp)),
                             n_pop = length(unique(dvr.stan$pophere))
                        )
)


m3l.ni = stan('stan/nointer_3levelwpop_mst_ncp.stan', data = datalist.dvr.pop,
              iter = 2000, warmup=1500, chains=2)

}



fls.stan <- subset(bbchmt, select=c("spp", "year", "site", "ind", "plot", "fls", "mst", "mwt"))
fls.stan <- fls.stan[complete.cases(fls.stan),]


getpop <- paste(fls.stan$spp, fls.stan$site)
fls.stan$pophere <- as.numeric(as.factor(getpop))
fls.stan$spp <- as.numeric(as.factor(fls.stan$spp))

datalist.fls.pop <- with(fls.stan, 
                         list(y = fls,  
                              mst = mst,
                              mwt = mwt,
                              sp = spp,
                              pop = pophere,
                              N = nrow(fls.stan),
                              n_sp = length(unique(fls.stan$spp)),
                              n_pop = length(unique(fls.stan$pophere))
                         )
)


m3l.ni = stan('stan/nointer_3levelwpop_mst&mwt_ncp.stan', data = datalist.fls.pop,
              iter = 2000, warmup=1500, chains=4, control=list(adapt_delta=0.999,max_treedepth = 15))




