##DB dichogamy

### Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

### First start with 2018 data...
# Set Working Directory
cg18 <-read.csv("fromMainRepo/2018_data/2018_CG_datasheet.csv", header=TRUE)
## Clean data
cg18<-gather(cg18, "date","bbch", -Ind, -Plot)
cg18<-na.omit(cg18)
cg18$date<-substr(cg18$date, 2,8)
cg18$date<-as.character(as.Date(cg18$date,"%m.%d.%y"))
cg18$doy<-yday(cg18$date)
cg18$species<-substr(cg18$Ind, 0,6)
cg18<-dplyr::select(cg18, -date)
cg18$species<-ifelse(cg18$species=="betpap", "BETPAP", cg18$species)



cg18m<-cg18 %>% dplyr::filter(grepl("-m", cg18$bbch))###only budsed

cg18f<-cg18 %>% dplyr::filter(grepl("-f", cg18$bbch))
cg18f<-cg18f %>% dplyr::filter(grepl("6", cg18f$bbch))

cg18M<-cg18 %>% dplyr::filter(grepl("-M", cg18$bbch))
cg18M<-cg18M %>% dplyr::filter(grepl("6", cg18M$bbch))

cg18F<-cg18 %>% dplyr::filter(grepl("-F", cg18$bbch))
cg18F<-cg18F %>% dplyr::filter(grepl("6", cg18F$bbch))

cg18<-rbind(cg18F,cg18M)

cg18<-cg18[!(cg18$bbch==""),]
###fornow remove 69


dx<-separate(cg18, bbch, into = c("first", "second"), sep = "\\,")
dx<-separate(dx, first, into = c("first", "third"), sep = "\\,")

dx$first<-ifelse(grepl("69", dx$first)==TRUE,NA,dx$first)
dx$first<-ifelse(grepl("67", dx$first)==TRUE,NA,dx$first)

dx$second<-ifelse(grepl("69", dx$second)==TRUE,NA,dx$second)
dx$second<-ifelse(grepl("67", dx$second)==TRUE,NA,dx$second)

dx$stamen<-NA
dx$stamen<-ifelse(grepl("M", dx$first)==TRUE,dx$doy, dx$stamen)
dx$stamen<-ifelse(grepl("M", dx$second)==TRUE, dx$doy, dx$stamen)

dx$pistil<-NA
dx$pistil<-ifelse(grepl("F", dx$first)==TRUE, dx$doy, dx$pistil)
dx$pistil<-ifelse(grepl("F", dx$second)==TRUE, dx$doy, dx$pistil)


female<-dx %>%group_by(Ind,Plot,species) %>% slice_min(pistil)
female<-dplyr::select(female,Ind,Plot,species,pistil)
male<-dx %>%group_by(Ind,Plot,species) %>% slice_min(stamen)                                                      
male<-dplyr::select(male,Ind,Plot,species,stamen)
goober<-full_join(male, female)

goober<-tidyr::gather(goober,"sexphase","day",4:5)

ggplot(goober,aes(species,day))+geom_point(aes(color=sexphase))

rm(cg20)
cg20 <-read.csv("fromMainRepo/2020_data/2020_CG_datasheet.csv", header=TRUE) 

## Now let's clean 2019 data
#cg20$id <- paste(cg20$ID, cg20$Plot, sep="_")
#cg20$Ind<-NULL
#cg20$Plot<-NULL
cg20 <- tidyr::gather(cg20, "date", "bbch",4:45)
cg20 <- na.omit(cg20)
cg20 <- cg20[!(cg20$bbch==""),]

cg20$date <- gsub("X", "", cg20$date)
cg20$date <- as.Date(cg20$date, format="%m.%d.%Y")
cg20$doy <- yday(cg20$date)
cg20flowers <- cg20[(cg20$Phase%in%c("Flowers")),]

cg20flowers<-cg20flowers %>% dplyr::filter(grepl("6", cg20flowers$bbch))

cg20flowers <- cg20flowers %>% separate(bbch,sep = ";", c('bbchM', 'bbchF'))

cg20flowers$bbchM<-ifelse(grepl("69", cg20flowers$bbchM)==TRUE,NA,cg20flowers$bbchM)
cg20flowers$bbchM<-ifelse(grepl("67", cg20flowers$bbchM)==TRUE,NA,cg20flowers$bbchM)

cg20flowers$bbchF<-ifelse(grepl("69", cg20flowers$bbchF)==TRUE,NA,cg20flowers$bbchF)
cg20flowers$bbchF<-ifelse(grepl("67", cg20flowers$bbchF)==TRUE,NA,cg20flowers$bbchF)

cg20flowers$stamen<-NA
cg20flowers$stamen<-ifelse(grepl("6", cg20flowers$bbchM)==TRUE, cg20flowers$doy, NA)

cg20flowers$pistil<-NA
cg20flowers$pistil<-ifelse(grepl("6", cg20flowers$bbchF)==TRUE, cg20flowers$doy, NA)

cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("SPI", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("DIE", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("VIB", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("ARO", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("SAMR", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("SOR", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("ACE", cg20flowers$ID))
cg20flowers<-cg20flowers %>% dplyr::filter(!grepl("AME", cg20flowers$ID))


female20<-cg20flowers %>%group_by(ID,Plot) %>% slice_min(pistil)

male20<-cg20flowers %>%group_by(ID,Plot) %>% slice_min(stamen)                                                      

female20<-dplyr::select(female20,ID,Plot,pistil)
male20<-dplyr::select(male20,ID,Plot,stamen)
goober2<-full_join(female20,male20)
goober2<-filter(goober2,!is.na(pistil))

goober2<-tidyr::gather(goober2,"sexphase","day",3:4)
goober2<-separate(data = goober2, col = ID, into = c("species", "site","ind"), sep = "\\_")

Alninc<-filter(goober2,species=="ALNINC")
ggplot(goober2,aes(site,day))+stat_summary(aes(color=sexphase),alpha=.7)+facet_wrap(~species)

