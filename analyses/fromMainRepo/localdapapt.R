#### local adaptation

###Aspriational titles
## Multi-species common garden study of temperate trees and shrubs
# reveals little evidence of local adaptation in phenology.

### Multi-species common garden study of temperate trees and shrubs reveals contrasting pattern of local adaptaiton
#in vegetative and reproductive phenology

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
library(brms)
library(tidybayes)


# Set Working Directory
setwd("~/Documents/git/wildhellgarden/analyses")
cg <- read.csv("output/clean_obs_allyrs.csv")
cg2<-tidyr::gather(cg,"phenophase","doy",6:14)

sites1<-dplyr::select(cg2,site,provenance.lat,provenance.long )
colnames()
unique(cg2$site)
lat<-unique(cg2$provenance.lat)
lon<-unique(cg2$provenance.lon)
sites<-data.frame(lat=lat,lon=lon)

AA<-data.frame(lat=42.31, lon=-71.12)
AA<-rbind(AA,AA,AA,AA)
library(geodist)
geodist::geodist(x =AA,y=sites,paired = TRUE ,measure = "haversine")

colnames(sites)[1:2]<-c("provenance.lat","provenance.long" )



#cg2<-left_join(cg2,sites)

library("leaflet")
library(maps)
map("state" , interior = FALSE)
?map()
points(sites$provenance.long, sites$provenance.lat, cex=.7,pch=16,col="skyblue3")

m<-leaflet()
m<-addTiles()

unique(cg2$phenophase)
cg2<- filter(cg2, phenophase %in% c(  "budburst" ,   "flowers" ,    "leafout"  ,     "ripefruit"  ,
 "budset"  ,    "leafcolor" ))
table(cg2$spp)


## select only species at least at 3 site
all3<-c("VIBCAS","ALNINC","SPIALB","SPITOM","DIELON","AMECAN","BETALL","MYRGAL","BETPAP","BETPOP","ACESPI")
All4<-c("VIBCAS","ALNINC","SPIALB","SPITOM","DIELON")
#cg2<-filter(cg2,spp %in% all3)
cg2$year<-cg2$year-2018 ##rescale year to 0-2
#cg2$transfer<-cg2$transfer/1000

##SPECIES AT 3 SITES
library(lme4)
mod.T.er<-lmer(doy~phenophase+(phenophase|site)+(phenophase|spp)+(phenophase|year),data=cg2)
mod.Trans.er<-lmer(doy~year+phenophase*site+(phenophase*site|spp),data=cg2)

summary(mod.T.er)
summary(mod.Trans.er)
coef(mod.T.er)


mod.T<-brm(doy~phenophase+(phenophase|site)+(phenophase|spp)+(phenophase|year),data=cg2,
           warmup=2000,iter=3000,chains =4,control=list(adapt_delta=.99),backend = "cmdstanr")

All4.dat<-filter(cg2,spp %in% All4)
mod.Tr<-brm(doy~phenophase*site+(phenophase*site|spp)+(phenophase*site|year),data=All3.dat,
           warmup=2000,iter=3000,chains =4,control=list(adapt_delta=.99),backend = "cmdstanr")

mod.Tr2<-brm(doy~phenophase*site*spp+(1|year),data=All4.dat,
            warmup=2000,iter=3000,chains =4,control=list(adapt_delta=.99),backend = "cmdstanr")

fixef(mod.Tr2)
#mod1<-brm(doy~year+phenophase+(phenophase|spp)+(phenophase|site),data=cg2,
 #              warmup=2000,iter=3000,chains =4, control=list(adapt_delta=.99),backend = "cmdstanr")

save.image("localadapt.Rda")

load("localadapt.Rda")

coef(mod1)
coef(mod.T,probs = c(.25,.75))$site
get_variables(mod.T)
get_variables(mod.Tr)


new.dater<-data.frame(spp=rep(unique(All4.dat$spp),each=4),site=rep(unique(All4.dat$site),5),phenophase=rep(unique(All4.dat$phenophase),each=20))
spsps<-mod.Tr2 %>% 
  epred_draws(newdata =new.dater,ndraws = 1000,re_formula = NA )

spsps2<-filter(spsps,phenophase=="leafout")
ggplot(spsps2,aes(.epred,spp))+stat_pointinterval(aes(color=site),.width=.5,position=pozy)+facet_wrap(~phenophase,scales = "free_x")

pop2<-mod.T %>% spread_draws(r_site[site,phenophase],r_spp[spp,phenophase])
unique(cg2$spp)


windy<-c("ALNINC","BETALL","BETPAP" ,"BETPOP","MYRGAL", "QUEALB" ,"QUERUB")
pop2$pol<-ifelse(pop2$spp %in% c(windy),"abiotic","biotic")
pop2$est<-pop2$r_site+pop2$r_spp
fruit<-filter(pop2,phenophase=="phenophaseripefruit")
ggplot(fruit,aes(est,spp))+stat_pointinterval(aes(color=site),position=pozy,.width=c(0.5,.89))+
  scale_color_viridis_d()+geom_vline(xintercept=0)+ggthemes::theme_base()+facet_wrap(~phenophase)


ggplot(pop2,aes(r_site,site))+stat_interval(aes(),position=pozy,.width=c(0.5,.89))+
  scale_color_viridis_d()+geom_vline(xintercept=0)+ggthemes::theme_base()+facet_wrap(~phenophase)


pop<-mod.T %>% spread_draws(c(sd_site__Intercept,sd_site__phenophaseleafout,sd_site__phenophaseflowers,
                              sd_site__phenophaseripefruit,sd_site__phenophasebudset,sd_site__phenophaseleafcolor))
spp<-mod.T %>% spread_draws(c(sd_spp__Intercept,sd_spp__phenophaseleafout,sd_spp__phenophaseflowers,
                             sd_spp__phenophaseripefruit,sd_spp__phenophasebudset,sd_spp__phenophaseleafcolor))

yrs<-mod.T %>% spread_draws(c(sd_year__Intercept,sd_year__phenophaseleafout,sd_year__phenophaseflowers,
                              sd_year__phenophaseripefruit,sd_year__phenophasebudset,sd_year__phenophaseleafcolor))

colnames(pop)[4:9]<-c("budburst"  ,  "budset","flowers"  ,  "leafcolor",  "leafout", "ripefruit")
colnames(spp)[4:9]<-c("budburst"  ,  "budset","flowers"  ,  "leafcolor",  "leafout", "ripefruit")
colnames(yrs)[4:9]<-c("budburst"  ,  "budset","flowers"  ,  "leafcolor",  "leafout", "ripefruit")
pop<-tidyr::gather(pop,"phenophase","estimate",4:9)
spp<-tidyr::gather(spp,"phenophase","estimate",4:9)
yrs<-tidyr::gather(yrs,"phenophase","estimate",4:9)
pop$var<-"intra-specific"
spp$var<-"inter-specific"
yrs$var<-"inter-annual"
var.data<-rbind(pop,spp)

pozy<-position_dodge(width=0.6)
  ggplot(pop,aes(reorder(phenophase,estimate),estimate))+stat_interval(size=3,.width = c(.5,.89),position=pozy)+
  coord_cartesian(ylim =c(0, 20))+
  scale_fill_viridis_d(option = "magma",begin=.2,alpha=.6)+ggthemes::theme_base()

ggplot(var.data,aes(reorder(phenophase,-estimate),estimate))+
  stat_interval(aes(color=var),.width = c(.5,.8,.95),position=pozy,alpha=0.3,size=10)+
  stat_pointinterval(aes(,color=var),size=5,.width = c(0),position=pozy)+
  #coord_cartesian(ylim =c(0, 10))+
  scale_color_viridis_d(option = "magma",begin=.2,end=.7)+ggthemes::theme_base()+geom_hline(yintercept=0)

ggplot(pop,aes(reorder(phenophase,-estimate),estimate))+
  stat_interval(aes(),.width = c(.1,.3,.5,.7,.9),position=pozy,size=20)+
  stat_pointinterval(aes(),size=10,.width = c(0),position=pozy)+xlab("phenophase")+
  ylab("population varience (stdv.)")+
  #coord_cartesian(ylim =c(0, 10))+
 ggthemes::theme_base(base_size = 11)+geom_hline(yintercept=0)

pop2<-left_join(pop2,sites1)
pop2<-left_join(pop2,sites)

ggplot(pop2,aes(r_site,site))+stat_interval(position=pozy,.width=c(.1,.3,.5,.7,.9),size=12)+
  stat_pointinterval(aes(),size=5,.width = c(0),position=pozy)+
  scale_color_viridis_d()+geom_vline(xintercept=0)+ggthemes::theme_base()+facet_wrap(~phenophase)








save.image("localadapt.Rda")
