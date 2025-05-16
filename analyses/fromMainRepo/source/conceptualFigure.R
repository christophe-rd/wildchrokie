#started Nov 11 by Deirdre
# Code to create the conceptual figure copied into a separate sourced file 

library(gridExtra)
library(ggplot2)
library(viridis)
library(ggpubr)
#### below makes the concept figure 
if(FALSE){
  delta<-diff(goo19$GDD_10)
  delta<-c(0,delta)
  goo19$ggdperday<-delta
  
  early<-filter(goo19,doy %in% c(115:240))
  early$GDD_10<-early$GDD_10-early$GDD_10[1]
  240-115
  142-115 #27 day diff
  253-240 #13
  27-13
  # day diff
  late<-filter(goo19,doy %in% c(142:253))
  late$GDD_10<-late$GDD_10-late$GDD_10[1]
  
  all<-filter(goo19,doy %in% c(115:253))
  
  early2<-filter(goo19,doy %in% c(127:240))
  early2$GDD_10<-early2$GDD_10-early2$GDD_10[1]
  
  late2<-filter(goo19,doy %in% c(147:245))
  late2$GDD_10<-late2$GDD_10-late2$GDD_10[1]
  
  ##start a month differnce
  ##end a week difference
  p1<-ggplot(data=all,aes(x=doy,y=ggdperday))+
    geom_rect(xmin=142,xmax=253, ymin=-2, ymax=0, alpha=0.02,fill="green3")+
    geom_rect(xmin=115,xmax=240, ymin=-3, ymax=-1, alpha=0.02,fill="darkgreen")+
    geom_rect(xmin=147,xmax=245, ymin=-6, ymax=-4, alpha=0.02,fill="skyblue1")+
    geom_rect(xmin=127,xmax=240, ymin=-7, ymax=-5, alpha=0.02,fill="skyblue3")+
    geom_line()+ggthemes::theme_few()+ylab("GDD's per day")+ylim(-8,20)+xlab("day of season")
  
  p2<-ggplot()+geom_line(data=early,aes(x=doy,y=GDD_10),color="darkgreen", linewidth=2)+
    geom_line(data=late,aes(x=doy,y=GDD_10),color="green3",size=2)+
    ggthemes::theme_few()+ylab("Cumulative GDD's")+xlab("day of season")
  p3<-ggplot()+geom_line(data=early2,aes(x=doy,y=GDD_10),color="skyblue1", linewidth=2)+
    geom_line(data=late2,aes(x=doy,y=GDD_10),color="skyblue3",size=2)+
    ggthemes::theme_few()+ylab("Cumulative GDD's")+xlab("day of season")
  
  
  concept2<-ggpubr::ggarrange(p2,p3,labels=c("b)","c)"))
  
  jpeg("figures/aronia_examp.jpeg", height=7,width=7,unit='in',res=200)
  ggpubr::ggarrange(p1,concept2,ncol=1,labels=c("a)",""))
  dev.off()
}

##do some basic plots
if(FALSE){
  ggplot(cg1,aes(pgsGDD,leafout))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~spp)
  ggplot(cg1,aes(pgs,leafout))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~spp)
  
  ggplot(cg1,aes(leafout,pgsGDD))+geom_point()+geom_smooth(method="lm")
  
  ggplot(cg1,aes(pgsGDD,leafout))+geom_point(aes(color=spp))+geom_smooth(method="lm",aes(color=spp))
  ggplot(cg1,aes(pgs,leafout))+geom_point(aes(color=spp))+geom_smooth(method="lm",aes(color=spp))
}

if(FALSE){#more plots
  jpeg("figures/gs_vs_leafout.jpeg",width=10,height=8,units = 'in',res=300)
  ggpubr::ggarrange(ggplot(cg1,aes(leafout,pgsGDD))+geom_point(aes(color=spp),size=.4)+geom_smooth(method="lm",aes(color=spp),se=FALSE)+geom_smooth(method="lm",se=FALSE,size=2),
    ggplot(cg1,aes(leafout,pgs))+geom_point(aes(color=spp),size=.4)+geom_smooth(method="lm",aes(color=spp),se=FALSE)+geom_smooth(method="lm",se=FALSE,size=2),
    ggplot(cg2,aes(leafout,fgsGDD))+geom_point(aes(color=spp),size=.4)+geom_smooth(method="lm",aes(color=spp),se=FALSE)+geom_smooth(method="lm",se=FALSE,size=2),
    ggplot(cg2,aes(leafout,fgs))+geom_point(aes(color=spp),size=.4)+geom_smooth(method="lm",aes(color=spp),se=FALSE)+geom_smooth(method="lm",se=FALSE,size=2),common.legend = TRUE)
  dev.off()
  
  ggpubr::ggarrange(ggplot(cg1,aes(pgsGDD,leafout))+geom_point(size=.4)+geom_smooth(method="lm",se=TRUE,size=1),
    ggplot(cg1,aes(pgs,leafout))+geom_point(size=.4)+geom_smooth(method="lm",se=TRUE,size=1),
    ggplot(cg2,aes(fgsGDD,leafout))+geom_point(size=.4)+geom_smooth(method="lm",se=TRUE,size=1),
    ggplot(cg2,aes(fgs,leafout))+geom_point(size=.4)+geom_smooth(method="lm",se=TRUE,size=1),common.legend = TRUE)
  dev.off()
  
  cg1 %>% group_by(spp,year) %>% summarise(meanpgs=mean(pgs),meanpgsGDD=mean(pgsGDD))
  
  ggplot(cg1,aes(pgs,pgsGDD))+geom_point()+geom_smooth(method="lm",aes(color=spp),se=FALSE)
  ggplot(cg2,aes(fgs,fgsGDD))+geom_point()+geom_smooth(method="lm",aes(color=spp),se=FALSE)
  
  jpeg("figures/gs_vs_leafout_site.jpeg",width=10,height=8,units = 'in',res=300)
  ggpubr::ggarrange(ggplot(cg1,aes(leafout,pgsGDD))+geom_point(aes(color=site),size=.4)+geom_smooth(method="lm",aes(color=site),se=FALSE),
    ggplot(cg1,aes(leafout,pgs))+geom_point(aes(color=site),size=.4)+geom_smooth(method="lm",aes(color=site),se=FALSE),
    ggplot(cg2,aes(leafout,fgsGDD))+geom_point(aes(color=site),size=.4)+geom_smooth(method="lm",aes(color=site),se=FALSE),
    ggplot(cg2,aes(leafout,fgs))+geom_point(aes(color=site),size=.4)+geom_smooth(method="lm",aes(color=site),se=FALSE),common.legend = TRUE)
  dev.off()
}
