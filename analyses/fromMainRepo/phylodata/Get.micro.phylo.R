## Script to build a phylogeny containing the species in Ospree departing from 
## the vascular plant megatree by Zanne et al. (2014);Nature

## Started by Ignacio Morales-Castilla on November 2018

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(pez)
library(caper)
library(phangorn)

## set your wd here:
setwd("~/Documents/git/microclimates")


# get the data 
gdd.stan <- read.csv("analyses/output/clean_gdd_bbanddvr.csv", header=TRUE)

gdd.stan$spp <- paste(substr(gdd.stan$genus, 0, 3), substr(gdd.stan$species,0,3), sep="")
gdd.stan$site <- NA
gdd.stan$site <- ifelse(gdd.stan$type=="Treespotters", "arb", gdd.stan$site)
gdd.stan$site <- ifelse(gdd.stan$type=="Common Garden", "cg", gdd.stan$site)
gdd.stan$site <- ifelse(gdd.stan$type=="Harvard Forest", "hf", gdd.stan$site)

gdd.stan <- subset(gdd.stan, select=c("id", "provenance.lat", "spp", "site",
                                      "gdd_bb", "gdd_dvr", "fs.count", "genus", "species"))

gdd.stan <- gdd.stan[!is.na(gdd.stan$provenance.lat),]
gdd.stan <- gdd.stan[!is.na(gdd.stan$gdd_bb),]
gdd.stan <- gdd.stan[!is.na(gdd.stan$gdd_dvr),]
gdd.stan <- gdd.stan[!is.na(gdd.stan$spp),]

gdd.stan <- gdd.stan[(gdd.stan$gdd_bb<1000),]
gdd.stan <- gdd.stan[(gdd.stan$gdd_dvr<1000),]

gdd.stan$spp <- ifelse(gdd.stan$spp=="NANA", "Quealb", gdd.stan$spp)


# geting a list of all species in ospree
gdd.stan$latbi <- paste(gdd.stan$genus, gdd.stan$species, sep="_")
#sort(unique(bb$latbi))

sps.list=sort(unique(gdd.stan$latbi))
genus.list=sort(unique(gdd.stan$genus))


## load phylogeny

## load phylo (from Zanne et al. 2014)
#phy.plants<-read.tree("data/phylogeny/Vascular_Plants_rooted.dated.tre")
phy.plants<-read.tree("~/Documents/git/microclimates/phylodata/Vascular_Plants_rooted.dated.tre")


## getting a list of genera in Zanne's phylo
phy.genera<-unlist(
  lapply(strsplit(phy.plants$tip.label, "_"),function(x){return(x[1])})
)

phy.genera.uniq<-sort(unique(phy.genera))



## how many microclimate genera are in the phylogeny?
microgenus.inphylo<-genus.list[which(genus.list%in%phy.genera.uniq)]


## first prune the phylogeny to include only these genera
phy.genera.micro<-drop.tip(phy.plants,
                            which(!phy.genera%in%microgenus.inphylo))


## we can add species that may not be present according to their genera
names.to.add=sps.list[which(!sps.list%in%phy.genera.micro$tip.label)]
phy.micro.clean<-congeneric.merge(phy.genera.micro,names.to.add,split="_")


## prunning the generated phylogeny to include micro species only
phy.plants.micro<-drop.tip(phy.micro.clean,
                            which(!phy.micro.clean$tip.label%in%sps.list))
# only 32 species are in the phylogeny

plot(phy.plants.micro,cex=.5)

## save phylogeny
write.tree(phy.plants.micro,"~/Documents/git/microclimates/phylodata/micro.phylogeny.tre")


### Let's build a prettier phylogeny now..
library(ggplot2)
library(ggtree)
library(phytools)
library(cowplot)

phy.plants.micro <- read.tree("~/Documents/git/microclimates/phylodata/micro.phylogeny.tre")

newnames <- data.frame(tiplabel=phy.plants.micro$tip.label)
newnames$genus <- gsub("_.*", "", newnames$tiplabel)
newnames$species <- gsub(".*_", "", newnames$tiplabel)
newnames$labels <- paste(newnames$genus, newnames$species)

tt <- as.data.frame(table(gdd.stan$site, gdd.stan$latbi))
names(tt) <- c("site", "phylo", "freq")
tt$freq <- ifelse(tt$freq==0, NA, tt$freq)
tt$freq <- ifelse(tt$site=="cg" & tt$phylo=="Quercus_alba", 1, tt$freq)

tt <- tt[!is.na(tt$freq),]
tt$freq <- NULL

arb <- tt[(tt$site=="arb"),]
arbids <- unique(arb$phylo)

hf <- tt[(tt$site=="hf"),]
hfids <- unique(hf$phylo)

cg <- tt[(tt$site=="cg"),]
cgids <- unique(cg$phylo)

tiplabels <- data.frame(tiplabel=unique(tt$phylo))
tiplabels$arb <- ifelse(tiplabels$tiplabel%in%arbids, 1, 0)
tiplabels$hf <- ifelse(tiplabels$tiplabel%in%hfids, 1, 0)
tiplabels$cg <- ifelse(tiplabels$tiplabel%in%cgids, 1, 0)
#tiplabels <- tiplabels[(tiplabels$tiplabel!="Quercus_alba"),]

arb.df<-subset(tiplabels, select=c("arb", "tiplabel"))
arb.df$tiplabel<-as.character(arb.df$tiplabel)
arb.df$arbcol <- ifelse(arb.df$arb==0, "white", "red")
cg.df<-subset(tiplabels, select=c("cg", "tiplabel"))
cg.df$tiplabel<-as.character(cg.df$tiplabel)
cg.df$cgcol <- ifelse(cg.df$cg==0, "white", "darkgreen")
hf.df<-subset(tiplabels, select=c("hf", "tiplabel"))
hf.df$tiplabel<-as.character(hf.df$tiplabel)
hf.df$hfcol <- ifelse(hf.df$hf==0, "white", "blue4")

foo <- left_join(newnames, arb.df)
foo <- left_join(foo, cg.df)
foo <- left_join(foo, hf.df)

colslegend <- data.frame(collabels=c("red", "darkgreen", "blue4"), site=c("Arboretum", "Common Garden", "Harvard Forest"))

leg <- ggplot(colslegend, aes(col=collabels)) + geom_point(aes(x=site, y=collabels)) + 
  theme(legend.key = element_rect(fill="transparent"), legend.text = element_text(size=7), 
        legend.title = element_text(size=8),
        plot.margin = unit(c(0,0,0,0), "lines")) +
  scale_colour_manual(name="Site", labels=c("Arboretum", "Common Garden", "Harvard Forest"),
                      values=c("red", "darkgreen", "blue4"))

legend <- cowplot::get_legend(leg)

##plot.margin: (t, r, b, l) *think 'trouble' ## 
tree.plot <- ggtree(phy.plants.micro)  + geom_tiplab(size=3, label=foo$labels, fontface='italic') + 
  theme(plot.margin = unit(c(1,0,1,1), "lines"), legend.position = "right") + geom_tippoint(col=foo$arbcol, x=175) + geom_tippoint(col=foo$cgcol, x=180) +
  geom_tippoint(col=foo$hfcol, x=185) +
  xlim(c(0, 200)) 

#tree.plot <- tree.plot %<+% tt + geom_tippoint(aes(group=c(tt$phylo, tt$site)))
#tree.plot+theme(legend.position="right")

quartz()
grid.arrange(tree.plot, legend, ncol=2, widths=c(3, 0.75))


