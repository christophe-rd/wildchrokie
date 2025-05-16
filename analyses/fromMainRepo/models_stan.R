## Looking at each major phase individually running stan
## Look at inter vs intra specific variation, is there local adaptation??
# 12 Feb 2021 by Cat

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

# Load libraries
library(gridExtra)
library(ggplot2)
library(viridis)
library(rstan)
library(dplyr)

# Set Working Directory
setwd("~/Documents/git/wildhellgarden/analyses")

cg <- read.csv("output/clean_obsdata_phases.csv")


###### Budburst model

bb.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "budburst"))
bb.stan <- bb.stan[complete.cases(bb.stan),]

###DB tries remove species without 3 or more populations represented
bb.stan<-filter(bb.stan, !spp %in% c("ACEPEN","ACESPI","AROMEL","QUERUB","QURALB","SAMRAC","SORAME","VACMYR"))

bb.stan2018<-dplyr::filter(bb.stan,year=="2018")
### do yearifying here.
getpop <- paste(bb.stan2018$spp, bb.stan2018$site)
bb.stan2018$pophere <- as.numeric(as.factor(getpop))
bb.stan2018$latbi<-bb.stan2018$spp
bb.stan2018$spp <- as.numeric(as.factor(bb.stan2018$spp))


datalist.bb.pop <- with(bb.stan2018, 
                         list(y = budburst,  
                              sp = spp,
                              pop = pophere,
                              N = nrow(bb.stan2018),
                              n_sp = length(unique(bb.stan2018$spp)),
                              n_pop = length(unique(bb.stan2018$pophere))
                         )
)


m3l.ni.bb = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.bb.pop,
                  iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.bb
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 

tidybayes::get_variables(m3l.ni.bb)
mod.sum[grep("a_sppop\\[", rownames(mod.sum)),] 
esto<-mod.sum[grep("a_sppop0\\[", rownames(mod.sum)),] 
## Get pop conversion
getpopall <- subset(bb.stan2018, select=c("spp", "site", "pophere"))
getpop <- getpopall[!duplicated(getpopall),]

concordance<-dplyr::select(bb.stan2018,spp,latbi)
concordance <- concordance[!duplicated(concordance),]

outy1<-left_join(getpop,concordance)
  outy1<-cbind(outy1,esto)

  level_order <- c("HF", "GR", "WM","SH")
pdf("figures/bb2018.pdf")  
  ggplot(outy1,aes(x=factor(site,level=level_order),mean))+geom_point()+geom_errorbar(aes(ymin=`25%`,ymax=`75%`))+facet_wrap(~latbi,scales = "free_y")
  dev.off()
  ggplot(outy1,aes(x=factor(site,level=level_order),mean))+geom_point()

###2019
  bb.stan2019<-dplyr::filter(bb.stan,year=="2019")
  ### do yearifying here.
  getpop <- paste(bb.stan2019$spp, bb.stan2019$site)
  bb.stan2019$pophere <- as.numeric(as.factor(getpop))
  bb.stan2019$latbi<-bb.stan2019$spp
  bb.stan2019$spp <- as.numeric(as.factor(bb.stan2019$spp))
  
  
  datalist.bb.pop <- with(bb.stan2019, 
                          list(y = budburst,  
                               sp = spp,
                               pop = pophere,
                               N = nrow(bb.stan2019),
                               n_sp = length(unique(bb.stan2019$spp)),
                               n_pop = length(unique(bb.stan2019$pophere))
                          )
  )
  
  
  m3l.ni.bb = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.bb.pop,
                   iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))
  
  modelhere <- m3l.ni.bb
  mod.sum <- summary(modelhere)$summary
  modtosave <- mod.sum[c(1:4),]
  mod.sum[grep("sigma", rownames(mod.sum)),] 
  
  tidybayes::get_variables(m3l.ni.bb)
  mod.sum[grep("a_sppop\\[", rownames(mod.sum)),] 
  esto19<-mod.sum[grep("a_sppop0\\[", rownames(mod.sum)),] 
  ## Get pop conversion
  getpopall <- subset(bb.stan2019, select=c("spp", "site", "pophere"))
  getpop <- getpopall[!duplicated(getpopall),]
  
  concordance<-dplyr::select(bb.stan2019,spp,latbi)
  concordance <- concordance[!duplicated(concordance),]
  
  outy2<-left_join(getpop,concordance)
  outy2<-cbind(outy2,esto19)
  
  level_order <- c("HF", "GR", "WM","SH")
  pdf("figures/bb2019.pdf") 
   ggplot()+geom_point(data=outy2,aes(x=factor(site,level=level_order),mean),color="blue")+geom_errorbar(data=outy2,aes(x=factor(site,level=level_order),ymin=`25%`,ymax=`75%`),color="blue")+facet_wrap(~latbi,scales = "free_y")
  
  dev.off()  
  ggplot(outy2,aes(x=factor(site,level=level_order),mean))+geom_point()
  
  ###2020
  bb.stan2020<-dplyr::filter(bb.stan,year=="2020")
  ### do yearifying here.
  getpop <- paste(bb.stan2020$spp, bb.stan2020$site)
  bb.stan2020$pophere <- as.numeric(as.factor(getpop))
  bb.stan2020$latbi<-bb.stan2020$spp
  bb.stan2020$spp <- as.numeric(as.factor(bb.stan2020$spp))
  
  
  datalist.bb.pop <- with(bb.stan2020, 
                          list(y = budburst,  
                               sp = spp,
                               pop = pophere,
                               N = nrow(bb.stan2020),
                               n_sp = length(unique(bb.stan2020$spp)),
                               n_pop = length(unique(bb.stan2020$pophere))
                          )
  )
  
  
  m3l.ni.bb = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.bb.pop,
                   iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))
  
  modelhere <- m3l.ni.bb
  mod.sum <- summary(modelhere)$summary
  modtosave <- mod.sum[c(1:4),]
  mod.sum[grep("sigma", rownames(mod.sum)),] 
  
  tidybayes::get_variables(m3l.ni.bb)
  mod.sum[grep("a_sppop\\[", rownames(mod.sum)),] 
  esto20<-mod.sum[grep("a_sppop0\\[", rownames(mod.sum)),] 
  ## Get pop conversion
  getpopall <- subset(bb.stan2020, select=c("spp", "site", "pophere"))
  getpop <- getpopall[!duplicated(getpopall),]
  
  concordance<-dplyr::select(bb.stan2020,spp,latbi)
  concordance <- concordance[!duplicated(concordance),]
  
  outy3<-left_join(getpop,concordance)
  outy3<-cbind(outy3,esto20)
  
  level_order <- c("HF", "GR", "WM","SH")
  pdf("figures/bb2020.pdf") 
  ggplot()+geom_point(data=outy3,aes(x=factor(site,level=level_order),mean),color="red")+geom_errorbar(data=outy3,aes(x=factor(site,level=level_order),ymin=`25%`,ymax=`75%`),color="red")+facet_wrap(~latbi,scales = "free_y")
  
  dev.off()
  
  
####Dan's going to do the same with bud set  
  bset.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "budset"))
  bset.stan <- bset.stan[complete.cases(bset.stan),]
  
  bset.stan2018<-dplyr::filter(bset.stan,year=="2019")
  ### do yearifying here.
  getpop <- paste(bset.stan2018$spp, bset.stan2018$site)
  bset.stan2018$pophere <- as.numeric(as.factor(getpop))
  bset.stan2018$latbi<-bset.stan2018$spp
  bset.stan2018$spp <- as.numeric(as.factor(bset.stan2018$spp))
  
  
  datalist.bset.pop <- with(bset.stan2018, 
                            list(y = budset,  
                                 sp = spp,
                                 pop = pophere,
                                 N = nrow(bset.stan2018),
                                 n_sp = length(unique(bset.stan2018$spp)),
                                 n_pop = length(unique(bset.stan2018$pophere))
                            )
  )
  
  
  m3l.ni.bset = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.bset.pop,
                     iter = 4000, warmup=3500, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))
  
  modelhere <- m3l.ni.bset
  mod.sum <- summary(modelhere)$summary
  modtosave <- mod.sum[c(1:4),]
  mod.sum[grep("sigma", rownames(mod.sum)),] 
  
  tidybayes::get_variables(m3l.ni.bset)
  mod.sum[grep("a_sppop\\[", rownames(mod.sum)),] 
  esto18bset<-mod.sum[grep("a_sppop0\\[", rownames(mod.sum)),] 
  ## Get pop conversion
  getpopall <- subset(bset.stan2018, select=c("spp", "site", "pophere"))
  getpop <- getpopall[!duplicated(getpopall),]
  
  concordance<-dplyr::select(bset.stan2018,spp,latbi)
  concordance <- concordance[!duplicated(concordance),]
  
  outy4<-left_join(getpop,concordance)
  outy4<-cbind(outy4,esto18bset)
  pdf("figures/bset2019.pdf") 
  ggplot()+geom_point(data=outy4,aes(x=factor(site,level=level_order),mean),color="black")+geom_errorbar(data=outy4,aes(x=factor(site,level=level_order),ymin=`25%`,ymax=`75%`),color="black")+facet_wrap(~latbi,scales = "free_y")
  dev.off()
  
  
  
  ###### Leafout model
lo.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "leafout"))
lo.stan <- lo.stan[complete.cases(lo.stan),]

getpop <- paste(lo.stan$spp, lo.stan$site)
lo.stan$pophere <- as.numeric(as.factor(getpop))
lo.stan$spp <- as.numeric(as.factor(lo.stan$spp))

datalist.lo.pop <- with(lo.stan, 
                        list(y = leafout,  
                             sp = spp,
                             pop = pophere,
                             N = nrow(lo.stan),
                             n_sp = length(unique(lo.stan$spp)),
                             n_pop = length(unique(lo.stan$pophere))
                        )
)


m3l.ni.lo = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.lo.pop,
                 iter = 4000, warmup=3500, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.lo
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 

###### Flowers model
flo.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "flowers"))
flo.stan <- flo.stan[complete.cases(flo.stan),]

getpop <- paste(flo.stan$spp, flo.stan$site)
flo.stan$pophere <- as.numeric(as.factor(getpop))
flo.stan$spp <- as.numeric(as.factor(flo.stan$spp))

datalist.flo.pop <- with(flo.stan, 
                        list(y = flowers,  
                             sp = spp,
                             pop = pophere,
                             N = nrow(flo.stan),
                             n_sp = length(unique(flo.stan$spp)),
                             n_pop = length(unique(flo.stan$pophere))
                        )
)


m3l.ni.flo = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.flo.pop,
                 iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.flo
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 

###### Fruits model
fru.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "ripefruit"))
fru.stan <- fru.stan[complete.cases(fru.stan),]

getpop <- paste(fru.stan$spp, fru.stan$site)
fru.stan$pophere <- as.numeric(as.factor(getpop))
fru.stan$spp <- as.numeric(as.factor(fru.stan$spp))

datalist.fru.pop <- with(fru.stan, 
                         list(y = ripefruit,  
                              sp = spp,
                              pop = pophere,
                              N = nrow(fru.stan),
                              n_sp = length(unique(fru.stan$spp)),
                              n_pop = length(unique(fru.stan$pophere))
                         )
)


m3l.ni.fru = stan('stan/nointer_3levelwpop_pheno_ncp.stan', data = datalist.fru.pop,
                  iter = 5000, warmup=4500, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.fru
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 

###### Budset model


##### DVR model
dvr.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "dvr"))
dvr.stan <- dvr.stan[complete.cases(dvr.stan),]


getpop <- paste(dvr.stan$spp, dvr.stan$site)
dvr.stan$pophere <- as.numeric(as.factor(getpop))
dvr.stan$spp <- as.numeric(as.factor(dvr.stan$spp))

datalist.dvr.pop <- with(dvr.stan, 
                         list(y = dvr,  
                              sp = spp,
                              pop = pophere,
                              N = nrow(dvr.stan),
                              n_sp = length(unique(dvr.stan$spp)),
                              n_pop = length(unique(dvr.stan$pophere))
                         )
)


m3l.ni.dvr = stan('stan/nointer_3levelwpop_interpheno_ncp.stan', data = datalist.dvr.pop,
              iter = 4000, warmup=3500, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.dvr
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 


##### Flowering model
floing.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "flowering"))
floing.stan <- floing.stan[complete.cases(floing.stan),]


getpop <- paste(floing.stan$spp, floing.stan$site)
floing.stan$pophere <- as.numeric(as.factor(getpop))
floing.stan$spp <- as.numeric(as.factor(floing.stan$spp))

datalist.floing.pop <- with(floing.stan, 
                         list(y = flowering,  
                              sp = spp,
                              pop = pophere,
                              N = nrow(floing.stan),
                              n_sp = length(unique(floing.stan$spp)),
                              n_pop = length(unique(floing.stan$pophere))
                         )
)


m3l.ni.floing = stan('stan/nointer_3levelwpop_interpheno_ncp.stan', data = datalist.floing.pop,
                  iter = 4000, warmup=3500, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.floing
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 


##### Fruiting model
fruing.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "fls"))
fruing.stan <- fruing.stan[complete.cases(fruing.stan),]


getpop <- paste(fruing.stan$spp, fruing.stan$site)
fruing.stan$pophere <- as.numeric(as.factor(getpop))
fruing.stan$spp <- as.numeric(as.factor(fruing.stan$spp))

datalist.fruing.pop <- with(fruing.stan, 
                            list(y = fls,  
                                 sp = spp,
                                 pop = pophere,
                                 N = nrow(fruing.stan),
                                 n_sp = length(unique(fruing.stan$spp)),
                                 n_pop = length(unique(fruing.stan$pophere))
                            )
)


m3l.ni.fruing = stan('stan/nointer_3levelwpop_interpheno_ncp.stan', data = datalist.fruing.pop,
                     iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.fruing
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 


##### Flower-Leaf Sequences model
fls.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "fls"))
fls.stan <- fls.stan[complete.cases(fls.stan),]


getpop <- paste(fls.stan$spp, fls.stan$site)
fls.stan$pophere <- as.numeric(as.factor(getpop))
fls.stan$spp <- as.numeric(as.factor(fls.stan$spp))

datalist.fls.pop <- with(fls.stan, 
                            list(y = fls,  
                                 sp = spp,
                                 pop = pophere,
                                 N = nrow(fls.stan),
                                 n_sp = length(unique(fls.stan$spp)),
                                 n_pop = length(unique(fls.stan$pophere))
                            )
)


m3l.ni.fls = stan('stan/nointer_3levelwpop_interpheno_ncp.stan', data = datalist.fls.pop,
                     iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.fls
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 


##### Growing Season model
gs.stan <- subset(cg, select=c("spp", "year", "site", "ind", "plot", "gs"))
gs.stan <- gs.stan[complete.cases(gs.stan),]
gs.stan<-dplyr::filter(gs.stan,year=="2020")

getpop <- paste(gs.stan$spp, gs.stan$site)
gs.stan$pophere <- as.numeric(as.factor(getpop))
gs.stan$spp <- as.numeric(as.factor(gs.stan$spp))

datalist.gs.pop <- with(gs.stan, 
                         list(y = gs,  
                              sp = spp,
                              pop = pophere,
                              N = nrow(gs.stan),
                              n_sp = length(unique(gs.stan$spp)),
                              n_pop = length(unique(gs.stan$pophere))
                         )
)


m3l.ni.gs = stan('stan/nointer_3levelwpop_interpheno_ncp.stan', data = datalist.gs.pop,
                  iter = 5000, warmup=4000, chains=4, control=list(adapt_delta=0.999, max_treedepth = 15))

modelhere <- m3l.ni.gs
mod.sum <- summary(modelhere)$summary
modtosave <- mod.sum[c(1:4),]
mod.sum[grep("sigma", rownames(mod.sum)),] 

tidybayes::get_variables(m3l.ni.gs)
mod.sum[grep("a_sppop\\[", rownames(mod.sum)),] 
gs18<-mod.sum[grep("a_sppop0\\[", rownames(mod.sum)),] 
## Get pop conversion
getpopall <- subset(gs.stan, select=c("spp", "site", "pophere"))
getpop <- getpopall[!duplicated(getpopall),]

concordance<-dplyr::select(gs.stan,spp,latbi)
concordance <- concordance[!duplicated(concordance),]

outygs2018<-left_join(getpop,concordance)
outygs2018<-cbind(outygs2018,gs18)
pdf("figures/bset2019.pdf") 
ggplot(data=outygs2018,aes(x=factor(site,level=level_order),mean))+geom_point()+geom_errorbar(aes(x=factor(site,level=level_order),ymin=`25%`,ymax=`75%`),color="black")+facet_wrap(~latbi,scales = "free_y")
dev.off()



