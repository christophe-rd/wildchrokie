## Loking at each major phase individually running stan
## Look at inter vs intra specific variation, is there local adaptation??
### Dan August 2022

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
cg <- read.csv("output/clean_obsdata_phases.csv")

surv<-read.csv("input/height_data.csv", header=TRUE)

#unique(surv$Notes)
#deathlist<-c("Dead","dead","Almost dead","...almost dead","photo attached. Probably dead","dead)","Dead)","almost dead","dead sampled?")
surv<-filter(surv,Year=="2022")
surv$alive<-NA
surv$alive<-ifelse(is.na(surv$Height)==TRUE,0,1)
table(surv$alive)
surv<-dplyr::select(surv,ID,Plot,alive)
surv<-filter(surv,!is.na(Plot))
cg$ID<-NA
cg$ID<-paste(cg$spp,cg$site,cg$ind,sep="_")

dat<-left_join(surv,cg)

dat<-dplyr::select(dat,ID,Plot,alive,spp,site,ind,provenance.lat)
dat<-distinct(dat)

surv.mod<-brm(alive~(1|spp)+(1|site),data=dat,family="bernoulli", warmup=3000, iter=4000,control=list(adapt_delta=.99))


get_variables(surv.mod)
survout<-surv.mod%>% spread_draws(r_site[site,Intercept])
survout<-left_join(survout,lat)
jpeg("figures/survplot.jpeg",width=8,height=8,units = "in",res=300)
ggplot(survout,aes(r_site,reorder(site,provenance.lat)))+stat_halfeye(.width = c(.5),size=2)+geom_vline(xintercept=0)+
  xlim(-3,3)+xlab("estimated effect")+ggthemes::theme_few()+ylab("")
dev.off()

survout2<-surv.mod %>% spread_draws(r_spp[spp,Intercept])
jpeg("figures/survplotspp.jpeg",width=8,height=8,units = "in",res=300)
ggplot(survout2,aes(r_spp,reorder(spp,r_spp)))+stat_halfeye(.width = c(.5),size=2)+geom_vline(xintercept=0)+
  xlim(-10,10)+xlab("estimated effect")+ggthemes::theme_few()+ylab("")
dev.off()



ranef(surv.mod,prob=c(.25,.75))$site
fixef(surv.mod)
newbers<-data.frame(spp=rep(unique(dat$spp),length(unique(dat$site))),site=rep(unique(dat$site),each=length(unique(dat$spp))))
newbers<-newbers[complete.cases(newbers),]
newbers<-data.frame(site=unique(dat$site))
newbers<-filter(newbers,!is.na(site))
pred<-fitted(surv.mod,newdata=newbers,probs = c(.05,.25,.75,.95),re_formula = NULL)

newbers<-cbind(newbers,pred)
pozz<-position_dodge(width=0.4)
ggplot(newbers,aes(reorder(site,Estimate),Estimate))+geom_point(aes(group=site),position=pozz)
  geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`),position=pozz,width=0)+
  geom_errorbar(aes(ymin=`Q5`,ymax=`Q95`),position=pozz,width=0,linetype="dotted")+
  ggthemes::theme_few(base_size = 9)



###phenology
checkers<-cg %>% group_by(spp,site) %>% count()

checkers<-checkers %>% group_by(spp,site) %>% count()

##SPECIES AT 4 SITES
All4<-c("VIBCAS","ALNINC","SPIALB","SPITOM","DIELON")

##SPECIES AT 3 SITES
GR.SH.WM<-c("ACESPI", "AMECAN","BETALL","MYRGAL")
GR.HF.SH<-c("BETPAP")
GR.HF.WM<-c("BETPOP")

#SPECIES AT 2 sites
HF.WM<-c("ACEPEN","AROMEL","QUERUB")
GR.HF<-c("SAMRAC")
GR.SH<-c("SORAME")
GR.WM<-c("VACMYR")

Big5<-filter(cg,spp %in% All4)

table(is.na(Big5$budburst)) #342 obs
table(is.na(Big5$leafout)) # 450
table(is.na(Big5$flowers)) #196
table(is.na(Big5$budset)) # 288
table(is.na(Big5$fruit)) # 228

table(is.na(Big5$flowering)) #87
table(is.na(Big5$fruiting)) #88
table(is.na(Big5$allflofruit)) #66

table(is.na(Big5$dvr)) # 337
table(is.na(Big5$fls)) # 191
table(is.na(Big5$gs)) #283

Med8<-filter(cg, spp %in% c(All4,GR.SH.WM))

Med10<-filter(cg, spp %in% c(All4,GR.SH.WM,GR.HF.SH,GR.HF.WM))


table(is.na(Med10$budburst)) #575 obs
table(is.na(Med10$leafout)) # 739
table(is.na(Med10$flowers)) #267
table(is.na(Med10$budset)) # 525
table(is.na(Med10$gs)) # 514


table(is.na(cg$budburst)) #674 obs
table(is.na(cg$leafout)) # 885


one<-ggplot(Med10,aes(spp,leafout))+geom_jitter(height=0,width=0.1,size=0.5)+facet_wrap(~year,nrow=3)
two<-ggplot(Med10,aes(spp,budset))+geom_jitter(height=0,width=0.1,size=0.5)+facet_wrap(~year,nrow=3)
three<-ggplot(Med10,aes(spp,flowers))+geom_jitter(height=0,width=0.1,size=0.5)+facet_wrap(~year,nrow=3)
four<-ggplot(Med10,aes(spp,fruit))+geom_jitter(height=0,width=0.1,size=0.5)+facet_wrap(~year,nrow=3)
five<-ggplot(Med10,aes(spp,gs))+geom_jitter(height=0,width=0.1,size=0.5)+facet_wrap(~year,nrow=3)
six<-ggplot(Med10,aes(spp,fls))+geom_jitter(height=0,width=0.1,size=0.5)+facet_wrap(~year,nrow=3)

jpeg("figures/data_patches.jpeg")
ggpubr::ggarrange(one,two,three,four)
dev.off()

###try a crossed model first
####do it with all species for now if we are interested in species level variation

Med10a<-filter(Med10,!is.na(gs))
Med10a<-filter(Med10a,spp!="ACESPI")
#mod.bb.10<-brm(budburst~year+(1|spp)+(1|site),data=Med10,
 #              warmup=3000,iter=4000, control=list(adapt_delta=.99))

#coef(mod.bb.10)$site
#mod.bb.10.fix<-brm(budburst~year+spp+(spp|site),data=Med10,
                #   warmup=4000,iter=5000, control=list(adapt_delta=.99))

use.data<-filter(cg,!is.na(gs))
use.data<-filter(use.data,!spp %in% c("ACEPEN","ACESPI","QUEALB","QUERUB","VACMYR")) ##remove species with less than 10 data points
table(use.data$spp)
#summary(mod.bb.10.fix,probs=c(.25,.75))
mod.lo.10<-brm(leafout~as.factor(year)+(1|spp)+(1|site),data=use.data,
              warmup=4000,iter=5000, control=list(adapt_delta=.99))

#mod.lo.10.fix<-brm(leafout~year+spp+(spp|site),data=Med10,
 #              warmup=3000,iter=4000, control=list(adapt_delta=.99))

#summary(mod.lo.10.fix,probs=c(.25,.75))
#mod.lo.cg<-brm(leafout~(1|spp)+(1|site)+(1|plot)+(1|year),data=cg,
          #     warmup=3000,iter=4000, control=list(adapt_delta=.99))

mod.bs.10<-brm(budset~as.factor(year)+(1|spp)+(1|site),data=use.data,
               warmup=4000,iter=5000, control=list(adapt_delta=.99))

#mod.bs.10.fix<-brm(budset~year+spp+(spp|site),data=Med10,
 #              warmu=4000,iter=5000, control=list(adapt_delta=.99))


mod.gs.10<-brm(gs~as.factor(year)+(1|spp)+(1|site),data=use.data,
               warmup=4000,iter=5000, control=list(adapt_delta=.99))



save.image("dansscratch.Rda")



#mod.lo.10<-brm(leafout~year+(1|spp)+(1|site),data=Med10,
 #              warmup=3000,iter=4000, control=list(adapt_delta=.99))

#mod.lo.10.fix<-brm(leafout~year+spp+(spp|site),data=Med10,
 #                  warmup=3000,iter=4000, control=list(adapt_delta=.99))

#summary(mod.lo.10.fix,probs=c(.25,.75))
#mod.lo.cg<-brm(leafout~(1|spp)+(1|site)+(1|plot)+(1|year),data=cg,
#     warmup=3000,iter=4000, control=list(adapt_delta=.99))

#mod.bs.10<-brm(budset~year+(1|spp)+(1|site),data=Med10,
   #            warmu=4000,iter=5000, control=list(adapt_delta=.99))

#mod.bs.10.fix<-brm(budset~year+spp+(spp|site),data=Med10,
  #                 warmu=4000,iter=5000, control=list(adapt_delta=.99))


#mod.gs.10<-brm(gs~year+(1|spp)+(1|site),data=Med10,
 #              warmup=4000,iter=5000, control=list(adapt_delta=.99))


#mod.fl.10<-brm(flowers~(1|spp)+(1|site)+(1|plot)+(1|year),data=Med10,
#               warmu=4000,iter=5000, control=list(adapt_delta=.995))


#mod.fr.10<-brm(fruit~(1|spp)+(1|site)+(1|plot)+(1|year),data=Med10,
 #              warmu=4000,iter=5000, control=list(adapt_delta=.995))

#launch_shinystan(mod.lo.10)
#get_variables(mod.lo.10)
lat<-dplyr::select(cg,provenance.lat,site)
lat<-distinct(lat)


unique(use.data$spp)
shrub<-c("ACEPEN" ,"ACESPI", "AMECAN" ,"AROMEL","DIELON","MYRGAL","SAMRAC" ,"SORAME", "SPIALB" ,"SPITOM" ,"VACMYR", "VIBCAS")
FLS<-c("ALNINC","AMECAN","BETALL","BETPAP","BETPOP","MYRGAL")


bb.sppout<-mod.lo.10 %>% spread_draws(r_spp[spp,Intercept])
bb.sppout$type<-ifelse(bb.sppout$spp %in% shrub,"shrub","tree")
bb.sppout$FLS<-ifelse(bb.sppout$spp %in% FLS,"flower-first","leaf-first")



p1s<-ggplot(bb.sppout,aes(r_spp,reorder(spp,r_spp)))+
  stat_halfeye(.width = c(.5),fill="lightgreen",alpha=0.6)+geom_vline(xintercept=0)+xlim(-15,15)+
  xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()

ggplot(bb.sppout,aes(r_spp,reorder(FLS,r_spp)))+
  stat_halfeye(.width = c(.5),fill="lightgreen",alpha=0.6)+geom_vline(xintercept=0)+xlim(-15,15)+
  xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()

ggplot(bb.sppout,aes(r_spp,reorder(type,r_spp)))+
  stat_halfeye(.width = c(.5),fill="lightgreen",alpha=0.6)+geom_vline(xintercept=0)+xlim(-15,15)+
  xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()


#lo.sppeout<-mod.lo.10 %>% spread_draws(r_spp[sspp,Intercept])
#p2<-ggplot(lo.siteout,aes(r_site,site))+stat_halfeye(.width = c(.5,.89))+geom_vline(xintercept=0)+xlim(-2.5,2.5)+xlab("Difference from grand intercept")



bs.sppout<-mod.bs.10 %>% spread_draws(r_spp[spp,Intercept])


bs.sppout$type<-ifelse(bs.sppout$spp %in% shrub,"shrub","tree")
bs.sppout$FLS<-ifelse(bs.sppout$spp %in% FLS,"flower-first","leaf-first")






p3s<-ggplot(bs.sppout,aes(r_spp,reorder(spp,r_spp)))+
  stat_halfeye(.width = c(.5),fill="darkorchid4",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-30,20)+xlab("estimated effect")+ylab("")+
  ggthemes::theme_few()

ggplot(bs.sppout,aes(r_spp,reorder(FLS,r_spp)))+
  stat_halfeye(.width = c(.5),fill="darkorchid4",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-30,20)+xlab("estimated effect")+ylab("")+
  ggthemes::theme_few()

ggplot(bs.sppout,aes(r_spp,reorder(type,r_spp)))+
  stat_halfeye(.width = c(.5),fill="darkorchid4",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-30,20)+xlab("estimated effect")+ylab("")+
  ggthemes::theme_few()


gs.sppout<-mod.gs.10 %>% spread_draws(r_spp[spp,Intercept])

gs.sppout$type<-ifelse(gs.sppout$spp %in% shrub,"shrub","tree")
gs.sppout$FLS<-ifelse(gs.sppout$spp %in% FLS,"flower-first","leaf-first")


p5a<-ggplot(gs.sppout,aes(r_spp,reorder(spp,r_spp)))+
  stat_halfeye(.width = c(.5),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-30,20)+xlab("estimated effect")+
ylab("")+
  ggthemes::theme_few()



pFLS1<-ggplot(bb.sppout,aes(r_spp,FLS))+
  stat_interval(.width = c(.5,.8,.975),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-35,35)+xlab("leafout")+
  ylab("")+
  ggthemes::theme_few()

pFLS2<-ggplot(gs.sppout,aes(r_spp,FLS))+
  stat_interval(.width = c(.5,.8,.975),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-35,35)+xlab("growing season")+
  ylab("")+
  ggthemes::theme_few()

ggpubr::ggarrange(pFLS1,pFLS2,common.legend=TRUE,ncol=1)



T1<-ggplot(bs.sppout,aes(r_spp,type))+
  stat_interval(.width = c(.5,.8,.975),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-35,35)+xlab("leafout")+
  ylab("")+
  ggthemes::theme_few()

T2<-ggplot(gs.sppout,aes(r_spp,type))+
  stat_interval(.width = c(.5,.8,.975),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-35,35)+xlab("Growing season")+
  ylab("")+
  ggthemes::theme_few()

jpeg("figures/FLS_n_type.jpeg",width=10,height=3,units = "in",res=300)
ggpubr::ggarrange(T1,pFLS1,T2,pFLS2,common.legend=TRUE,ncol=2,nrow=2,labels=c("a)","b)","",""))
dev.off()


ggplot(gs.sppout,aes(r_spp,reorder(type,r_spp)))+
  stat_interval(.width = c(.5,.80,.975),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-30,20)+xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()

ggplot(gs.sppout,aes(r_spp,reorder(FLS,r_spp)))+
  stat_halfeye(.width = c(.5),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-30,20)+xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()


p4s<-ggpubr::ggarrange(p1s,p3s,ncol=2,labels = c("b)","c)"))

jpeg("figures/sppplots.jpeg",width=8,height=8,units = "in",res=300)
ggpubr::ggarrange(p5a,p4s,nrow=2,labels=c("a)"))
dev.off()
cor(use.data$budburst,use.data$budset,use = "complete.obs")


#####species effects
bb.siteout<-mod.lo.10 %>% spread_draws(r_site[site,Intercept])
bb.siteout<-left_join(bb.siteout,lat)
p1<-ggplot(bb.siteout,aes(r_site,reorder(site,provenance.lat)),)+
  stat_halfeye(.width = c(.5),fill="lightgreen",alpha=0.6)+geom_vline(xintercept=0)+xlim(-2.5,2.5)+
  xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()




lo.siteout<-mod.lo.10 %>% spread_draws(r_site[site,Intercept])
#p2<-ggplot(lo.siteout,aes(r_site,site))+stat_halfeye(.width = c(.5,.89))+geom_vline(xintercept=0)+xlim(-2.5,2.5)+xlab("Difference from grand intercept")

bs.siteout<-mod.bs.10 %>% spread_draws(r_site[site,Intercept])
bs.siteout<-left_join(bs.siteout,lat)

p3<-ggplot(bs.siteout,aes(r_site,reorder(site,provenance.lat)))+
  stat_halfeye(.width = c(.5),fill="darkorchid4",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-10,10)+xlab("estimated effect")+ylab("")+
  ggthemes::theme_few()

gs.siteout<-mod.gs.10 %>% spread_draws(r_site[site,Intercept])
gs.siteout<-left_join(gs.siteout,lat)
p5<-ggplot(gs.siteout,aes(r_site,reorder(site,provenance.lat)))+
  stat_halfeye(.width = c(.5),fill="brown",alpha=0.6)+geom_vline(xintercept=0)+
  xlim(-10,10)+xlab("estimated effect")+
  ylab("")+
  ggthemes::theme_few()

bb.siteout %>% group_by(site) %>% summarise(meangs=mean(r_site))
0.722+.668 ##1.5

bs.siteout %>% group_by(site) %>% summarise(meangs=mean(r_site))
1.99+1.37 ##3-4 days diff
gs.siteout %>% group_by(site) %>% summarise(meangs=mean(r_site))
1.29+2.52 #### 3-4day difference


bb.sppout %>% group_by(spp) %>% summarise(meangs=mean(r_spp))
11.3+10.1 ##21. 4 days

bs.sppout %>% group_by(spp) %>% summarise(meangs=mean(r_spp))
10.7+16.3## 27 days

gs.sppout %>% group_by(spp) %>% summarise(meangs=mean(r_spp))
10.5+11 #(21 daty)

p4<-ggpubr::ggarrange(p1,p3,ncol=2,labels = c("b)","c)"))

jpeg("figures/siteplots.jpeg",width=8,height=8,units = "in",res=300)
ggpubr::ggarrange(p5,p4,nrow=2,labels=c("a)"))
dev.off()





bb.pop<-mod.lo.10 %>% spread_draws(sd_site__Intercept)
bb.sp<-mod.lo.10 %>% spread_draws(sd_spp__Intercept)
#bb.yr<-mod.bb.10 %>% spread_draws(sd_year__Intercept)
#bb.plot<-mod.bb.10 %>% spread_draws(sd_plot__Intercept)
colnames(bb.pop)[4]<-"sd"
bb.pop$var<-"pop"
colnames(bb.sp)[4]<-"sd"
bb.sp$var<-"spp"
#colnames(bb.yr)[4]<-"sd"
#bb.yr$var<-"yrs"
#colnames(bb.plot)[4]<-"sd"
#bb.plot$var<-"plot"

var.bb<-rbind(bb.sp,bb.pop)



lo.pop<-mod.gs.10 %>% spread_draws(sd_site__Intercept)
lo.sp<-mod.gs.10 %>% spread_draws(sd_spp__Intercept)
#lo.yr<-mod.lo.10 %>% spread_draws(sd_year__Intercept)
#lo.plot<-mod.lo.10 %>% spread_draws(sd_plot__Intercept)
colnames(lo.pop)[4]<-"sd"
lo.pop$var<-"pop"
colnames(lo.sp)[4]<-"sd"
lo.sp$var<-"spp"
#colnames(lo.yr)[4]<-"sd"
#lo.yr$var<-"yrs"
#colnames(lo.plot)[4]<-"sd"
#lo.plot$var<-"plot"

var<-rbind(lo.sp,lo.pop)#,lo.yr,lo.plot)

bs.pop<-mod.bs.10 %>% spread_draws(sd_site__Intercept)
bs.sp<-mod.bs.10 %>% spread_draws(sd_spp__Intercept)
#bs.yr<-mod.bs.10 %>% spread_draws(sd_year__Intercept)
#bs.plot<-mod.bs.10 %>% spread_draws(sd_plot__Intercept)

colnames(bs.pop)[4]<-"sd"
bs.pop$var<-"pop"
colnames(bs.sp)[4]<-"sd"
bs.sp$var<-"spp"
#colnames(bs.yr)[4]<-"sd"
#bs.yr$var<-"yrs"
#colnames(bs.plot)[4]<-"sd"
#bs.plot$var<-"plot"

var.bs<-rbind(bs.sp,bs.pop)#,bs.yr,bs.plot)
var$phenophases<-"season length"
var.bs$phenophases<-"budset"
var.bb$phenophases<-"leafout"
var.veg<-rbind(var.bs,var.bb)
ggplot(var.veg,aes(sd,var))+stat_halfeye(aes(fill=phenophases,color=phenophases),alpha=0.5)+xlim(0,30)+ggthemes::theme_few()
ggplot(var.bs,aes(sd,var))+stat_halfeye()


##reproductive


ggplot(var.all,aes(sd,var))+stat_halfeye(aes(fill=phenophases,color=phenophases),alpha=0.5)+xlim(0,40)+ggthemes::theme_few()

pozy<-position_dodge(width=0.6)
ggplot(var.veg,aes(phenophases,sd))+stat_halfeye(aes(color=var),alpha=0.9,position=pozy,.width = c(.5))+ggthemes::theme_few()+scale_color_viridis_d(end=0.8,option = "C")+
  scale_fill_viridis_d(end=0.8)+
  scale_x_discrete(limits=c("budburst","leafout","set"))+ylim(0,15)+ylab("standard deviation")

jpeg("figures/varpar_veg.jpeg",width=7,height=8,units = "in",res=300)
ggplot(var.veg,aes(var,sd))+stat_halfeye(aes(shape=phenophases,fill=phenophases),size=2,alpha=0.6,position=pozy,.width = c(.5))+ggthemes::theme_few()+scale_color_viridis_d(end = 0.9)+
  scale_fill_viridis_d(end=0.8,direction = -1)+ylim(0,15)+xlab("source of variance")+ylab("standard deviation")+geom_hline(yintercept=0)+scale_shape_manual(values=c(1,8))
dev.off()
#+
  scale_x_discrete(limits=c("budburst","leafout","budset"))
  
  
ggplot(var.rep,aes(phenophases,sd))+stat_pointinterval(aes(shape=var),position=pozy,.width = c(.5,.8))+ggthemes::theme_few()+scale_color_viridis_d()+
  scale_x_discrete(limits=c("flowering","fruiting"))


stop("not an error, just stop")
mod.lo.5<-brm(leafout~(1|spp)+(1|site)+(1|plot)+(1|year),data=Big5,
               warmup=3000,iter=4000, control=list(adapt_delta=.99))

mod.bs.5<-brm(budset~(1|spp)+(1|site)+(1|plot)+(1|year),data=Big5,
               warmu=4000,iter=5000, control=list(adapt_delta=.99))

mod.fl.5<-brm(flowers~(1|spp)+(1|site)+(1|plot)+(1|year),data=Big5,
               warmu=4000,iter=5000, control=list(adapt_delta=.995))


mod.fr.5<-brm(fruit~(1|spp)+(1|site)+(1|plot)+(1|year),data=Big5,
               warmu=4000,iter=5000, control=list(adapt_delta=.995))







     

     
     
     
     ####### TRY AGAIN WITH LATITUDE AS FIXEFF INSTEAD
mod.lo.10.lat<-brm(leafout~provenance.lat+(1|year)+(provenance.lat|spp)+(1|plot),data=Med10,
                    warmup=4000,iter=5000, control=list(adapt_delta=.99))
fixef(mod.lo.10.lat,probs = c(.25,.75))       
coef(mod.lo.10.lat,probs = c(.25,.75))    
     

     ranef(mod.lo.5,probs = c(.25,.75))$site
     ranef(mod.fl.5,probs = c(.25,.75))$site
     ranef(mod.fr.5,probs = c(.25,.75))$site
     ranef(mod.bs.5,probs = c(.25,.75))$site
   
     ranef(mod.lo.10,probs = c(.25,.75))$site
     ranef(mod.fl.10,probs = c(.25,.75))$site
     ranef(mod.fr.10,probs = c(.25,.75))$site
     ranef(mod.bs.10,probs = c(.25,.75))$site


newbies<-data.frame(spp=rep(unique(Med10$spp),length(unique(Med10$site))),site=rep(unique(Med10$site),each=length(unique(Med10$spp))))

newbies2<-rbind(newbies,newbies,newbies)
newbies2$year<-rep(unique(Med10$year),each=length(newbies$spp))
newbies2$plot<-as.factor(100)

pred<-predict(mod.lo.10,newdata = newbies2,probs = c(.25,.75),allow_new_levels=TRUE)
pred<-cbind(newbies2,pred)
pred2<-predict(mod.fl.10,newdata = newbies2,probs = c(.25,.75),allow_new_levels=TRUE)
pred2<-cbind(newbies2,pred2)
ggplot(pred,aes(reorder(site,-Estimate),Estimate))+geom_point(aes(color=spp))+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=spp),width=0)+facet_grid(spp~year)
ggplot(pred2,aes(reorder(site,-Estimate),Estimate))+geom_point(aes(color=spp))+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=spp),width=0)+facet_grid(spp~year)

pred3<-predict(mod.bs.10,newdata = newbies2,probs = c(.25,.75),allow_new_levels=TRUE)
pred3<-cbind(newbies2,pred3)


poz<-position_dodge(width=10)
ggplot(pred,aes(0,Estimate))+geom_point(aes(color=reorder(spp,Estimate),shape=site),position=poz,size=2)+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=reorder(spp,Estimate),linetype=site),width=0,position=poz,alpha=0.3)+facet_wrap(~year,nrow=3)+
  scale_shape_manual(values=c(15,16,17,18))+scale_linetype_manual(values=c(rep("solid",4)))+scale_color_viridis_d(option = "A")+scale_x_discrete(name="Leafout")+ggthemes::theme_few()

ggplot(pred2,aes(0,Estimate))+geom_point(aes(color=reorder(spp,Estimate),shape=site),position=poz,size=2)+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=reorder(spp,Estimate),linetype=site),width=0,position=poz,alpha=0.3)+facet_wrap(~year,nrow=3)+
  scale_shape_manual(values=c(15,16,17,18))+scale_linetype_manual(values=c(rep("solid",4)))+scale_color_viridis_d()+scale_x_discrete(name="Flowering")+ggthemes::theme_few()

ggplot(pred3,aes(0,Estimate))+geom_point(aes(color=reorder(spp,Estimate),shape=site),position=poz,size=2)+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=reorder(spp,Estimate),linetype=site),width=0,position=poz,alpha=0.3)+facet_wrap(~year,nrow=3)+
  scale_shape_manual(values=c(15,16,17,18))+scale_linetype_manual(values=c(rep("solid",4)))+scale_color_viridis_d(option="C")+scale_x_discrete(name="Budset")+ggthemes::theme_few()
fl.pop<-mod.fl.10 %>% spread_draws(sd_site__Intercept)
fl.sp<-mod.fl.10 %>% spread_draws(sd_spp__Intercept)
fl.yr<-mod.fl.10 %>% spread_draws(sd_year__Intercept)
fl.plot<-mod.fl.10 %>% spread_draws(sd_plot__Intercept)

colnames(fl.pop)[4]<-"sd"
fl.pop$var<-"pop"
colnames(fl.sp)[4]<-"sd"
fl.sp$var<-"spp"
colnames(fl.yr)[4]<-"sd"
fl.yr$var<-"yrs"
colnames(fl.plot)[4]<-"sd"
fl.plot$var<-"plot"
var.fl<-rbind(fl.sp,fl.pop,fl.yr,fl.plot)

fr.pop<-mod.fr.10 %>% spread_draws(sd_site__Intercept)
fr.sp<-mod.fr.10 %>% spread_draws(sd_spp__Intercept)
fr.yr<-mod.fr.10 %>% spread_draws(sd_year__Intercept)
fr.plot<-mod.fr.10 %>% spread_draws(sd_plot__Intercept)

colnames(fr.pop)[4]<-"sd"
fr.pop$var<-"pop"
colnames(fr.sp)[4]<-"sd"
fr.sp$var<-"spp"
colnames(fr.yr)[4]<-"sd"
fr.yr$var<-"yrs"
colnames(fr.plot)[4]<-"sd"
fr.plot$var<-"plot"
var.fr<-rbind(fr.sp,fr.pop,fr.yr,fr.plot)


var.fl$phenophases<-"flowering"
var.fr$phenophases<-"fruiting"

var.rep<-rbind(var.fl,var.fr)

var.all<-rbind(var.rep,var.veg)

#mod.lo.10 %>% spread_draws(sd_year__Intercept)
#mod.lo.10 %>% spread_draws(sd_plot__Intercept) %>% head(10)

save.image("dansscratch.Rda")
