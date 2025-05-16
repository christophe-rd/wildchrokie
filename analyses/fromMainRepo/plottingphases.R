## Looking at each major phase individually
## Let's make some simple plots to start and then we can think about stan
# 12 Feb 2021 by Cat

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Load libraries
library(gridExtra)
library(ggplot2)
library(viridis)


# Set Working Directory
setwd("~/Documents/git/wildhellgarden/analyses")

my.pal <-rep(viridis_pal(option="viridis")(9),2)
my.pch <- rep(16:17, each=9)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

###### Now let's add in the clean bbch data
cg <- read.csv("output/clean_obs_allyrs.csv")
cg$spp <- ifelse(cg$spp=="VACMR", "VACMYR", cg$spp)

cg$species<- NA
cg$genus<- NA

cg$genus<-ifelse(cg$spp=="ACEPEN", "Acer", cg$genus)
cg$species<-ifelse(cg$spp=="ACEPEN", "pensylvanicum", cg$species)
cg$genus<-ifelse(cg$spp=="ACESPI", "Acer", cg$genus)
cg$species<-ifelse(cg$spp=="ACESPI", "spicatum", cg$species)
cg$genus<-ifelse(cg$spp=="ALNINC", "Alnus", cg$genus)
cg$species<-ifelse(cg$spp=="ALNINC", "incana", cg$species)
cg$genus<-ifelse(cg$spp=="AMECAN", "Amelanchier", cg$genus)
cg$species<-ifelse(cg$spp=="AMECAN", "canadensis", cg$species)
cg$genus<-ifelse(cg$spp=="AROMEL", "Aronia", cg$genus)
cg$species<-ifelse(cg$spp=="AROMEL", "melonacarpa", cg$species)
cg$genus<-ifelse(cg$spp=="BETALL", "Betula", cg$genus)
cg$species<-ifelse(cg$spp=="BETALL", "alleghaniensis", cg$species)
cg$genus<-ifelse(cg$spp=="BETPAP", "Betula", cg$genus)
cg$species<-ifelse(cg$spp=="BETPAP", "papyrifera", cg$species)
cg$genus<-ifelse(cg$spp=="BETPOP", "Betula", cg$genus)
cg$species<-ifelse(cg$spp=="BETPOP", "populifolia", cg$species)
cg$genus<-ifelse(cg$spp=="BETPOPX", "Betula", cg$genus)
cg$species<-ifelse(cg$spp=="BETPOPX", "populifolia", cg$species)
cg$genus<-ifelse(cg$spp=="DIELON", "Diervilla", cg$genus)
cg$species<-ifelse(cg$spp=="DIELON", "lonicera", cg$species)
cg$genus<-ifelse(cg$spp=="MYRGAL", "Myrica", cg$genus)
cg$species<-ifelse(cg$spp=="MYRGAL", "gale", cg$species)
cg$genus<-ifelse(cg$spp=="QUERUB", "Quercus", cg$genus)
cg$species<-ifelse(cg$spp=="QUERUB", "rubra", cg$species)
cg$genus<-ifelse(cg$spp=="SAMRAC", "Sambucus", cg$genus)
cg$species<-ifelse(cg$spp=="SAMRAC", "racemosa", cg$species)
cg$genus<-ifelse(cg$spp=="SORAME", "Sorbus", cg$genus)
cg$species<-ifelse(cg$spp=="SORAME", "americana", cg$species)
cg$genus<-ifelse(cg$spp=="SPIALB", "Spiraea", cg$genus)
cg$species<-ifelse(cg$spp=="SPIALB", "alba", cg$species)
cg$genus<-ifelse(cg$spp=="SPITOM", "Spiraea", cg$genus)
cg$species<-ifelse(cg$spp=="SPITOM", "tomentosa", cg$species)
cg$genus<-ifelse(cg$spp=="VACMYR", "Vaccinium", cg$genus)
cg$species<-ifelse(cg$spp=="VACMYR", "myrtilloides", cg$species)
cg$genus<-ifelse(cg$spp=="VIBCAS", "Viburnum", cg$genus)
cg$species<-ifelse(cg$spp=="VIBCAS", "cassinoides", cg$species)
cg$genus<-ifelse(cg$spp=="QUEALB", "Quercus", cg$genus)
cg$species<-ifelse(cg$spp=="QUEALB", "alba", cg$species)

cg$species.name <- paste(cg$genus, cg$species, sep=" ")

bb <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "budburst", "provenance.lat", "species.name"))
bb <- bb[complete.cases(bb),]
lo <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "leafout", "provenance.lat", "species.name"))
lo <- lo[complete.cases(lo),]
flo <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "flowers", "provenance.lat", "species.name"))
flo <- flo[complete.cases(flo),]
fru <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "ripefruit", "provenance.lat", "species.name"))
fru <- fru[complete.cases(fru),]
bset <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "budset", "provenance.lat", "species.name"))
bset <- bset[complete.cases(bset),]


bbplot <- ggplot(cg[!is.na(cg$budburst),], aes(x=provenance.lat, y=budburst, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(cg$species.name[!is.na(cg$budburst)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(cg$species.name[!is.na(cg$budburst)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")


#quartz()
#bbplot

loplot <- ggplot(cg[!is.na(cg$leafout),], aes(x=provenance.lat, y=leafout, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(lo$species.name[!is.na(cg$budburst)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(lo$species.name[!is.na(cg$budburst)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")


#quartz()
#loplot

floplot <- ggplot(cg[!is.na(cg$flowers),], aes(x=provenance.lat, y=flowers, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(flo$species.name[!is.na(cg$flowers)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(flo$species.name[!is.na(cg$flowers)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position="none")


#quartz()
#floplot

fruplot <- ggplot(cg[!is.na(cg$ripefruit),], aes(x=provenance.lat, y=ripefruit, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(fru$species.name[!is.na(cg$ripefruit)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(fru$species.name[!is.na(cg$ripefruit)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")


#quartz()
#fruplot


bsetplot <- ggplot(cg[!is.na(cg$budset),], aes(x=provenance.lat, y=budset, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(bset$species.name[!is.na(cg$budset)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(bset$species.name[!is.na(cg$budset)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 9), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 8, face="italic"), legend.margin=margin(c(1,1,1,1)))


#quartz()
#bsetplot

mylegend<-g_legend(bsetplot)

bsetplot <- ggplot(cg[!is.na(cg$budset),], aes(x=provenance.lat, y=budset, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(bset$species.name[!is.na(cg$budset)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(bset$species.name[!is.na(cg$budset)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.position="none",
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)))

quartz()
grid.arrange(bbplot, loplot, floplot, fruplot, bsetplot, mylegend, ncol=3)

png("figures/phenophase_sppop_raw.png", 
    width=10,
    height=6, units="in", res = 350 )
grid.arrange(bbplot, loplot, floplot, fruplot, bsetplot, mylegend, ncol=3)
dev.off()


cg$dvr <- cg$leafout - cg$budburst
cg$flowering <- cg$flowers - cg$flobudburst
cg$fruiting <- cg$ripefruit - cg$fruit
cg$allflofruit <- cg$ripefruit - cg$flobudburst
cg$fls <- cg$flowers - cg$leafout
cg$gs <- cg$budset - cg$leafout


dvrplot <- ggplot(cg[!is.na(cg$dvr),], aes(x=provenance.lat, y=dvr, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(fru$species.name[!is.na(cg$dvr)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(fru$species.name[!is.na(cg$dvr)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")

floingplot <- ggplot(cg[!is.na(cg$flowering),], aes(x=provenance.lat, y=flowering, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(fru$species.name[!is.na(cg$flowering)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(fru$species.name[!is.na(cg$flowering)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")

fruingplot <- ggplot(cg[!is.na(cg$fruiting),], aes(x=provenance.lat, y=fruiting, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(fru$species.name[!is.na(cg$fruiting)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(fru$species.name[!is.na(cg$fruiting)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")

flsplot <- ggplot(cg[!is.na(cg$fls),], aes(x=provenance.lat, y=fls, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(fru$species.name[!is.na(cg$fls)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(fru$species.name[!is.na(cg$fls)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")

gsplot <- ggplot(cg[!is.na(cg$gs),], aes(x=provenance.lat, y=gs, col=as.factor(species.name), shape=as.factor(species.name))) + 
  geom_jitter(aes(col=as.factor(species.name)), height=0.25) +
  scale_color_manual(name = "Species", values=my.pal, labels = unique(fru$species.name[!is.na(cg$gs)])) +
  scale_shape_manual(name = "Species", values=my.pch, labels = unique(fru$species.name[!is.na(cg$gs)])) +
  theme_classic() + xlab("Provenance Latitude") + ylab("Day of Year") +
  theme(legend.title = element_text(size = 7), legend.key.size=unit(0.5, "line"),
        legend.text = element_text(size = 6, face="italic"), legend.margin=margin(c(1,1,1,1)),
        legend.position = "none")


png("figures/interphenophase_sppop_raw.png", 
    width=10,
    height=6, units="in", res = 350 )
grid.arrange(dvrplot, floingplot, fruingplot, flsplot, gsplot, mylegend, ncol=3)
dev.off()


write.csv(cg, file="output/clean_obsdata_phases.csv", row.names=FALSE)
