## 16 July 2019 - Cat
# Aim is to investigate the relationship between growth traits and the duration of vegetative risk
# Below is code for 2018 and then initial code for 2019

### Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
options(max.print = 200) 

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")

### First start with 2018 data...
# Set Working Directory
cg18 <-read.csv("fromMainRepo/2018_data/2018_CG_datasheet.csv", header=TRUE)

## Clean data
cg18<-gather(cg18, "date","bbch", -Ind, -Plot) # bbch: international convention of dev stages of plants
cg18<-na.omit(cg18)
cg18$date<-substr(cg18$date, 2,8) # removes characters from position 2 and 8
cg18$date<-as.character(as.Date(cg18$date,"%m.%d.%y"))
cg18$doy<-yday(cg18$date)
cg18$species<-substr(cg18$Ind, 0,6)
cg18<-dplyr::select(cg18, -date)
cg18$species<-ifelse(cg18$species=="betpap", "BETPAP", cg18$species)
cg18$bbch<-gsub("-m", " ", cg18$bbch, fixed=TRUE)
cg18$bbch<-gsub("-M", " ", cg18$bbch, fixed=TRUE)
cg18$bbch<-gsub("-f", " ", cg18$bbch, fixed=TRUE)
cg18$bbch<-gsub("-F", " ", cg18$bbch, fixed=TRUE)


cg18<-cg18[!(cg18$bbch==""),]
dx<-separate(cg18, bbch, into = c("first", "second"), sep = "\\,") # why is 2 phenophases?
dx<-separate(dx, first, into = c("first", "third"), sep = "\\,")

#dx$first <- substr(dx$first, 0, 2)
#dx$second <- substr(dx$second, 0, 2)
#dx$third <- substr(dx$third, 0, 2)

dx$bb<-NA
dx$bb<-ifelse(dx$first=="9" | dx$first=="9-" | dx$first=="11" | dx$second=="9" | dx$second=="9-" |
                dx$second=="11" | dx$third=="9" | dx$third=="9-" | dx$third=="11", dx$doy, dx$bb) # it checks whether the first,2nd or 3d column match the condition, and if it does, it assigns a date. 9 and 11 likely refers to early budburst
dx$lo<-NA
dx$lo<-ifelse(dx$first=="19" | dx$second=="19" | dx$third=="19", dx$doy, dx$lo)

dx$first <- as.numeric(dx$first)
dx$second <- as.numeric(dx$second)
dx$third <- as.numeric(dx$third)
dx$flobuds<-NA
dx$flobuds<-ifelse(dx$first%in%c(51:59) | dx$second%in%c(51:59) | dx$third%in%c(51:59), dx$doy, dx$flobuds)

dx$flobudburst<-NA
dx$flobudburst<-ifelse(dx$first%in%c(60:67) | dx$second%in%c(60:67) | dx$third%in%c(60:67), dx$doy, dx$flobudburst)

dx$flowers<-NA
dx$flowers<-ifelse(dx$first==69 | dx$second==69 | dx$third==69, dx$doy, dx$flowers)

dx$fruit<-NA
dx$fruit<-ifelse(dx$first%in%c(70:73) | dx$second%in%c(70:73) | dx$third%in%c(70:73), dx$doy, dx$fruit)

dx$ripefruit<-NA
dx$ripefruit<-ifelse(dx$first==79 | dx$second==79 | dx$third==79, dx$doy, dx$ripefruit)

dx$budset<-NA
dx$budset<-ifelse(dx$first=="102" | dx$second=="102" | dx$third=="102", dx$doy, dx$budset)

if(FALSE){
### Some checks added by Cat 2 May 2023 - all about the same amount of data so adding both
### Check to see leaf coloration (92) vs leaf drop (93, 95, 97)
nrow(dx[dx$first==92,]) #532
nrow(dx[dx$second==92,]) # 9800

## Leaf drop
nrow(dx[dx$first==93,]); nrow(dx[dx$first==95,]); nrow(dx[dx$first==97,]) #426, 481, 353
nrow(dx[dx$second==93,]); nrow(dx[dx$second==95,]); nrow(dx[dx$second==97,]) # 9748, 9729, 9698
}


dx$leafcolor<-NA
dx$leafcolor<-ifelse(dx$first==92 | dx$second==92 | dx$third==92, dx$doy, dx$leafcolor)

dx$leafdrop<-NA
dx$leafdrop<-ifelse(dx$first%in%c(93:97) | dx$second%in%c(93:97) | dx$third%in%c(93:97), dx$doy, dx$leafdrop)


drisk<-dx%>%dplyr::select(Ind, Plot, bb, lo, flobuds, flobudburst, flowers, fruit, ripefruit, budset, leafcolor, leafdrop, species)
#drisk<-drisk[!(is.na(drisk$bb) & is.na(drisk$lo)),]

bb<-drisk[!is.na(drisk$bb),]
bb$budburst<-ave(bb$bb, bb$Ind, bb$Plot, FUN=min) # averaging the column bb for grouping variable ind and plot
bb<-subset(bb, select=c("Ind", "Plot", "budburst"))
bb<-bb[!duplicated(bb),]
lo<-drisk[!is.na(drisk$lo),]
lo$leafout<- ave(lo$lo, lo$Ind, lo$Plot, FUN=min) 
lo<-subset(lo, select=c("Ind", "Plot", "leafout"))
lo<-lo[!duplicated(lo),]
fbud<-drisk[!is.na(drisk$flobuds),]
fbud$flobuds<- ave(fbud$flobuds, fbud$Ind, fbud$Plot, FUN=min) 
fbud<-subset(fbud, select=c("Ind", "Plot", "flobuds"))
fbud<-fbud[!duplicated(fbud),]
fbb<-drisk[!is.na(drisk$flobudburst),]
fbb$flobudburst<- ave(fbb$flobudburst, fbb$Ind, fbb$Plot, FUN=min) 
fbb<-subset(fbb, select=c("Ind", "Plot", "flobudburst"))
fbb<-fbb[!duplicated(fbb),]
flos<-drisk[!is.na(drisk$flowers),]
flos$flowers<- ave(flos$flowers, flos$Ind, flos$Plot, FUN=min) 
flos<-subset(flos, select=c("Ind", "Plot", "flowers"))
flos<-flos[!duplicated(flos),]
fru<-drisk[!is.na(drisk$fruit),]
fru$fruit<- ave(fru$fruit, fru$Ind, fru$Plot, FUN=min) 
fru<-subset(fru, select=c("Ind", "Plot", "fruit"))
fru<-fru[!duplicated(fru),]
ripe<-drisk[!is.na(drisk$ripefruit),]
ripe$ripefruit<- ave(ripe$ripefruit, ripe$Ind, ripe$Plot, FUN=min) 
ripe<-subset(ripe, select=c("Ind", "Plot", "ripefruit"))
ripe<-ripe[!duplicated(ripe),]
bset<-drisk[!is.na(drisk$budset),]
bset$budset<- ave(bset$budset, bset$Ind, bset$Plot, FUN=min) 
bset<-subset(bset, select=c("Ind", "Plot", "budset"))
bset<-bset[!duplicated(bset),]
lcolor<-drisk[!is.na(drisk$leafcolor),]
lcolor$leafcolor<- ave(lcolor$leafcolor, lcolor$Ind, lcolor$Plot, FUN=min) 
lcolor<-subset(lcolor, select=c("Ind", "Plot", "leafcolor"))
lcolor<-lcolor[!duplicated(lcolor),]
ldrop<-drisk[!is.na(drisk$leafdrop),]
ldrop$leafdrop<- ave(ldrop$leafdrop, ldrop$Ind, ldrop$Plot, FUN=min) 
ldrop<-subset(ldrop, select=c("Ind", "Plot", "leafdrop"))
ldrop<-ldrop[!duplicated(ldrop),]

cg18clean<-full_join(bb, lo)
cg18clean <- full_join(cg18clean, fbud)
cg18clean <- full_join(cg18clean, fbb)
cg18clean <- full_join(cg18clean, flos)
cg18clean <- full_join(cg18clean, fru)
cg18clean <- full_join(cg18clean, ripe)
cg18clean<-full_join(cg18clean, bset)
cg18clean<-full_join(cg18clean, lcolor)
cg18clean<-full_join(cg18clean, ldrop)
cg18clean$risk<-cg18clean$leafout-cg18clean$budburst 
cg18clean$year <- 2018
#write.csv(cg18clean, file="output/clean_cg_2018.csv", row.names=FALSE)

### Now some starter code for 2019!
# Set Working Directory
cg19 <-read.csv("fromMainRepo/2019_data/2019_CG_dataupdated.csv", header=TRUE) 

## Now let's clean 2019 data
cg19$treeid <- paste(cg19$ID, cg19$Plot, sep="_")
cg19$Ind<-NULL
cg19$Plot<-NULL
cg19 <- gather(cg19, "date", "bbch", -treeid, -Phase)
cg19 <- cg19[!(cg19$date=="ID"),]
cg19 <- na.omit(cg19)
cg19 <- cg19[!(cg19$bbch==""),]

cg19$date <- gsub("X", "", cg19$date)
cg19$date <- as.Date(cg19$date, format="%m.%d.%Y")
cg19$doy <- yday(cg19$date)

cg19leaves <- cg19[(cg19$Phase=="Leaves"),]

cg19leaves <- subset(cg19leaves, select=c("treeid", "doy", "bbch"))
cg19leaves <- separate(data = cg19leaves, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg19leaves$ind <- ifelse(is.na(cg19leaves$ind), substr(cg19leaves$spp, 7,8), cg19leaves$ind)
cg19leaves$ind <- ifelse(cg19leaves$ind=="", "XX", cg19leaves$ind)
cg19leaves$spp <- substr(cg19leaves$spp, 0, 6)
cg19leaves$year <- 2019

cg19leaves<-cg19leaves[!duplicated(cg19leaves),]

cg19leaves$bb <- ifelse(cg19leaves$bbch%in%c(9:11), cg19leaves$doy, NA) # bb likely stands for bud burst
cg19leaves$lo <- ifelse(cg19leaves$bbch==19, cg19leaves$doy, NA) # lo stands for leaf out


cg19leaves$spindplot <- paste(cg19leaves$spp, cg19leaves$site, cg19leaves$ind, cg19leaves$plot)
cg19leaves$budburst <- NA
for(i in c(unique(cg19leaves$spindplot))){ 
  
  budburst <- cg19leaves$bb[i==cg19leaves$spindplot][1]
  cg19leaves$budburst[i==cg19leaves$spindplot] <- budburst
  
}
cg19leaves$leafout <- cg19leaves$lo

cg19leaves <- subset(cg19leaves, select=c("spp", "year", "site", "ind", "budburst", "leafout", "plot"))
cg19leaves$plot <- as.character(cg19leaves$plot)

cg19budset <- cg19[(cg19$Phase=="Budset"),]

cg19budset <- subset(cg19budset, select=c("treeid", "doy", "bbch"))
cg19budset <- separate(data = cg19budset, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg19budset$ind <- ifelse(is.na(cg19budset$ind), substr(cg19budset$spp, 7,8), cg19budset$ind)
cg19budset$ind <- ifelse(cg19budset$ind=="", "XX", cg19budset$ind)
cg19budset$spp <- substr(cg19budset$spp, 0, 6)
cg19budset$year <- 2019

cg19budset <-cg19budset%>% 
  group_by(spp, site, ind, plot, bbch, year) %>% 
  slice(which.min(doy))
cg19budset<-cg19budset[!duplicated(cg19budset),]

cg19budset$bbch <- gsub("\\.", "\\,", cg19budset$bbch)
cg19budset <- separate(data = cg19budset, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg19budset$budset <- ifelse(cg19budset$bbch==102, cg19budset$doy, NA)

cg19budset <- cg19budset %>% 
  group_by(spp, site, ind, plot, year) %>% 
  summarise_all(list(~first(na.omit(.))))

cg19budset <- subset(cg19budset, select=c("spp", "year", "site", "ind", "budset", "plot"))
cg19budset$plot <- as.character(cg19budset$plot)

cg19budset <- cg19budset[!duplicated(cg19budset),]

cg19clean <- full_join(cg19leaves, cg19budset)

cg19flowers <- cg19[(cg19$Phase=="Flowers"),]

cg19flowers <- subset(cg19flowers, select=c("treeid", "doy", "bbch"))
cg19flowers <- separate(data = cg19flowers, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg19flowers$ind <- ifelse(is.na(cg19flowers$ind), substr(cg19flowers$spp, 7,8), cg19flowers$ind)
cg19flowers$ind <- ifelse(cg19flowers$ind=="", "XX", cg19flowers$ind)
cg19flowers$spp <- substr(cg19flowers$spp, 0, 6)
cg19flowers$year <- 2019

cg19flowers<-cg19flowers[!duplicated(cg19flowers),]

cg19flowers$bbch <- gsub("\\.", "\\,", cg19flowers$bbch)
cg19flowers <- separate(data = cg19flowers, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg19flowers$fbud <- ifelse(cg19flowers$bbch%in%c(51:59), cg19flowers$doy, NA)
cg19flowers$fbb <- ifelse(cg19flowers$bbch%in%c(60:62), cg19flowers$doy, NA)
cg19flowers$flos <- ifelse(cg19flowers$bbch==69, cg19flowers$doy, NA)

cg19flowers$spindplot <- paste(cg19flowers$spp, cg19flowers$site, cg19flowers$ind, cg19flowers$plot)
cg19flowers$flobuds <- NA
for(i in c(unique(cg19flowers$spindplot))){ 
  
  flobuds <- cg19flowers$fbud[i==cg19flowers$spindplot][1]
  cg19flowers$flobuds[i==cg19flowers$spindplot] <- flobuds
  
}
cg19flowers$flobudburst <- NA
for(i in c(unique(cg19flowers$spindplot))){ 
  
  flobudburst <- cg19flowers$fbb[i==cg19flowers$spindplot][1]
  cg19flowers$flobudburst[i==cg19flowers$spindplot] <- flobudburst
  
}
cg19flowers$flowers <- cg19flowers$flos

cg19flowers <- subset(cg19flowers, select=c("spp", "year", "site", "ind", "flobuds", "flobudburst", "flowers", "plot"))
cg19flowers$plot <- as.character(cg19flowers$plot)

cg19flowers <- cg19flowers[!duplicated(cg19flowers),]

cg19clean <- full_join(cg19clean, cg19flowers)

cg19fruits <- cg19[(cg19$Phase=="Fruits"),]

cg19fruits <- subset(cg19fruits, select=c("treeid", "doy", "bbch"))
cg19fruits <- separate(data = cg19fruits, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg19fruits$ind <- ifelse(is.na(cg19fruits$ind), substr(cg19fruits$spp, 7,8), cg19fruits$ind)
cg19fruits$ind <- ifelse(cg19fruits$ind=="", "XX", cg19fruits$ind)
cg19fruits$spp <- substr(cg19fruits$spp, 0, 6)
cg19fruits$year <- 2019

cg19fruits<-cg19fruits[!duplicated(cg19fruits),]

cg19fruits$bbch <- gsub("\\.", "\\,", cg19fruits$bbch)
cg19fruits <- separate(data = cg19fruits, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg19fruits$fru <- ifelse(cg19fruits$bbch%in%c(70:74), cg19fruits$doy, NA)
cg19fruits$ripe <- ifelse(cg19fruits$bbch==79, cg19fruits$doy, NA)

cg19fruits$spindplot <- paste(cg19fruits$spp, cg19fruits$site, cg19fruits$ind, cg19fruits$plot)
cg19fruits$fruit <- NA
for(i in c(unique(cg19fruits$spindplot))){ 
  
  fruit <- cg19fruits$fru[i==cg19fruits$spindplot][1]
  cg19fruits$fruit[i==cg19fruits$spindplot] <- fruit
  
}
cg19fruits$ripefruit <- cg19fruits$ripe

cg19fruits <- subset(cg19fruits, select=c("spp", "year", "site", "ind", "fruit", "ripefruit", "plot"))
cg19fruits$plot <- as.character(cg19fruits$plot)

cg19fruits <- cg19fruits[!duplicated(cg19fruits),]
cg19clean <- full_join(cg19clean, cg19fruits)


cg19fall <- cg19[(cg19$Phase=="Fall Colors"),]

cg19fall <- subset(cg19fall, select=c("treeid", "doy", "bbch"))
cg19fall <- separate(data = cg19fall, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg19fall$ind <- ifelse(is.na(cg19fall$ind), substr(cg19fall$spp, 7,8), cg19fall$ind)
cg19fall$ind <- ifelse(cg19fall$ind=="", "XX", cg19fall$ind)
cg19fall$spp <- substr(cg19fall$spp, 0, 6)
cg19fall$year <- 2019

cg19fall<-cg19fall[!duplicated(cg19fall),]

cg19fall$bbch <- gsub("\\.", "\\,", cg19fall$bbch)
cg19fall <- separate(data = cg19fall, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg19fall$leafcolor <- ifelse(cg19fall$bbch%in%c(90:99), cg19fall$doy, NA)

cg19fall$spindplot <- paste(cg19fall$spp, cg19fall$site, cg19fall$ind, cg19fall$plot)
cg19fall$lcolor <- NA
for(i in c(unique(cg19fall$spindplot))){ 
  
  lcolor <- cg19fall$leafcolor[i==cg19fall$spindplot][1]
  cg19fall$lcolor[i==cg19fall$spindplot] <- lcolor
  
}
cg19fall$leafcolor <- cg19fall$lcolor

cg19fall <- subset(cg19fall, select=c("spp", "year", "site", "ind", "leafcolor", "plot"))
cg19fall$plot <- as.character(cg19fall$plot)

cg19fall <- cg19fall[!duplicated(cg19fall),]
cg19clean <- full_join(cg19clean, cg19fall)


cg18clean <- separate(data = cg18clean, col = Ind, into = c("spp", "site", "ind"), sep = "\\_")
cg18clean$plot <- as.character(cg18clean$Plot)
cg18clean$Plot <- NA

cg <- full_join(cg19clean, cg18clean)

### Now some starter code for 2020!
# Set Working Directory
cg20 <-read.csv("fromMainRepo/2020_data/2020_CG_datasheet.csv", header=TRUE) 

## Now let's clean 2019 data
cg20$treeid <- paste(cg20$ID, cg20$Plot, sep="_")
cg20$Ind<-NULL
cg20$Plot<-NULL
cg20 <- gather(cg20, "date", "bbch", -treeid, -Phase)
cg20 <- na.omit(cg20)
cg20 <- cg20[!(cg20$bbch==""),]

cg20$date <- gsub("X", "", cg20$date)
cg20$date <- as.Date(cg20$date, format="%m.%d.%Y")
cg20$doy <- yday(cg20$date)

cg20leaves <- cg20[(cg20$Phase=="Leaves"),]

cg20leaves <- subset(cg20leaves, select=c("treeid", "doy", "bbch"))
cg20leaves <- separate(data = cg20leaves, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg20leaves$ind <- ifelse(is.na(cg20leaves$ind), substr(cg20leaves$spp, 7,8), cg20leaves$ind)
cg20leaves$ind <- ifelse(cg20leaves$ind=="", "XX", cg20leaves$ind)
cg20leaves$spp <- substr(cg20leaves$spp, 0, 6)
cg20leaves$year <- 2020

cg20leaves<-cg20leaves[!duplicated(cg20leaves),]
cg20leaves$bbch <- as.numeric(cg20leaves$bbch)

cg20leaves$bb <- ifelse(cg20leaves$bbch%in%c(9:11), cg20leaves$doy, NA)
cg20leaves$lo <- ifelse(cg20leaves$bbch==19, cg20leaves$doy, NA)


cg20leaves$spindplot <- paste(cg20leaves$spp, cg20leaves$site, cg20leaves$ind, cg20leaves$plot)
cg20leaves$budburst <- NA
for(i in c(unique(cg20leaves$spindplot))){ #i="QUERUB HF 99 5"
  
  budburst <- cg20leaves$bb[!is.na(cg20leaves$bb) & i==cg20leaves$spindplot][1]
  cg20leaves$budburst[i==cg20leaves$spindplot] <- budburst
  
}
cg20leaves$leafout <- cg20leaves$lo

cg20leaves <- subset(cg20leaves, select=c("spp", "year", "site", "ind", "budburst", "leafout", "plot"))
cg20leaves$plot <- as.character(cg20leaves$plot)

cg20budset <- cg20[(cg20$Phase%in%c("Leaves", "Flowers")),]

cg20budset <- subset(cg20budset, select=c("treeid", "doy", "bbch"))
cg20budset <- separate(data = cg20budset, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg20budset$ind <- ifelse(is.na(cg20budset$ind), substr(cg20budset$spp, 7,8), cg20budset$ind)
cg20budset$ind <- ifelse(cg20budset$ind=="", "XX", cg20budset$ind)
cg20budset$spp <- substr(cg20budset$spp, 0, 6)
cg20budset$year <- 2020

cg20budset <-cg20budset%>% 
  group_by(spp, site, ind, plot, bbch, year) %>% 
  slice(which.min(doy))
cg20budset<-cg20budset[!duplicated(cg20budset),]

cg20budset$bbch <- gsub("\\.", "\\,", cg20budset$bbch)
cg20budset <- separate(data = cg20budset, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg20budset$budset <- ifelse(cg20budset$bbch==102, cg20budset$doy, NA)

cg20budset <- cg20budset %>% 
  group_by(spp, site, ind, plot, year) %>% 
  summarise_all(list(~first(na.omit(.))))

cg20budset <- subset(cg20budset, select=c("spp", "year", "site", "ind", "budset", "plot"))
cg20budset$plot <- as.character(cg20budset$plot)

cg20budset <- cg20budset[!duplicated(cg20budset),]

cg20clean <- full_join(cg20leaves, cg20budset)

cg20flowers <- cg20[(cg20$Phase%in%c("Leaves", "Flowers")),]

cg20flowers <- subset(cg20flowers, select=c("treeid", "doy", "bbch"))
cg20flowers <- separate(data = cg20flowers, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg20flowers$ind <- ifelse(is.na(cg20flowers$ind), substr(cg20flowers$spp, 7,8), cg20flowers$ind)
cg20flowers$ind <- ifelse(cg20flowers$ind=="", "XX", cg20flowers$ind)
cg20flowers$spp <- substr(cg20flowers$spp, 0, 6)
cg20flowers$year <- 2020

cg20flowers<-cg20flowers[!duplicated(cg20flowers),]

cg20flowers$bbch <- gsub("\\.", "\\,", cg20flowers$bbch)
cg20flowers <- separate(data = cg20flowers, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg20flowers$bbch <- as.numeric(cg20flowers$bbch)

cg20flowers$fbud <- ifelse(cg20flowers$bbch%in%c(51:59), cg20flowers$doy, NA)
cg20flowers$fbb <- ifelse(cg20flowers$bbch%in%c(60:62), cg20flowers$doy, NA)
cg20flowers$flos <- ifelse(cg20flowers$bbch==69, cg20flowers$doy, NA)

cg20flowers$spindplot <- paste(cg20flowers$spp, cg20flowers$site, cg20flowers$ind, cg20flowers$plot)
cg20flowers$flobuds <- NA
for(i in c(unique(cg20flowers$spindplot))){ 
  
  flobuds <- cg20flowers$fbud[i==cg20flowers$spindplot][1]
  cg20flowers$flobuds[i==cg20flowers$spindplot] <- flobuds
  
}
cg20flowers$flobudburst <- NA
for(i in c(unique(cg20flowers$spindplot))){ 
  
  flobudburst <- cg20flowers$fbb[i==cg20flowers$spindplot][1]
  cg20flowers$flobudburst[i==cg20flowers$spindplot] <- flobudburst
  
}
cg20flowers$flowers <- cg20flowers$flos

cg20flowers <- subset(cg20flowers, select=c("spp", "year", "site", "ind", "flobuds", "flobudburst", "flowers", "plot"))
cg20flowers$plot <- as.character(cg20flowers$plot)

cg20flowers <- cg20flowers[!duplicated(cg20flowers),]

cg20clean <- full_join(cg20clean, cg20flowers)

cg20fruits <- cg20[(cg20$Phase%in%c("Leaves", "Flowers")),]

cg20fruits <- subset(cg20fruits, select=c("treeid", "doy", "bbch"))
cg20fruits <- separate(data = cg20fruits, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg20fruits$ind <- ifelse(is.na(cg20fruits$ind), substr(cg20fruits$spp, 7,8), cg20fruits$ind)
cg20fruits$ind <- ifelse(cg20fruits$ind=="", "XX", cg20fruits$ind)
cg20fruits$spp <- substr(cg20fruits$spp, 0, 6)
cg20fruits$year <- 2020

cg20fruits<-cg20fruits[!duplicated(cg20fruits),]

cg20fruits$bbch <- gsub("\\.", "\\,", cg20fruits$bbch)
cg20fruits <- separate(data = cg20fruits, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg20fruits$bbch <- as.numeric(cg20fruits$bbch)

cg20fruits$fru <- ifelse(cg20fruits$bbch%in%c(70:74), cg20fruits$doy, NA)
cg20fruits$ripe <- ifelse(cg20fruits$bbch==79, cg20fruits$doy, NA)

cg20fruits$spindplot <- paste(cg20fruits$spp, cg20fruits$site, cg20fruits$ind, cg20fruits$plot)
cg20fruits$fruit <- NA
for(i in c(unique(cg20fruits$spindplot))){ 
  
  fruit <- cg20fruits$fru[i==cg20fruits$spindplot][1]
  cg20fruits$fruit[i==cg20fruits$spindplot] <- fruit
  
}
cg20fruits$ripefruit <- cg20fruits$ripe

cg20fruits <- subset(cg20fruits, select=c("spp", "year", "site", "ind", "fruit", "ripefruit", "plot"))
cg20fruits$plot <- as.character(cg20fruits$plot)

cg20fruits <- cg20fruits[!duplicated(cg20fruits),]
cg20clean <- full_join(cg20clean, cg20fruits)


cg20fall <- cg20[(cg20$Phase%in%c("Leaves", "Flowers")),]

cg20fall <- subset(cg20fall, select=c("treeid", "doy", "bbch"))
cg20fall <- separate(data = cg20fall, col = treeid, into = c("spp", "site", "ind", "plot"), sep = "\\_")
cg20fall$ind <- ifelse(is.na(cg20fall$ind), substr(cg20fall$spp, 7,8), cg20fall$ind)
cg20fall$ind <- ifelse(cg20fall$ind=="", "XX", cg20fall$ind)
cg20fall$spp <- substr(cg20fall$spp, 0, 6)
cg20fall$year <- 2020

cg20fall<-cg20fall[!duplicated(cg20fall),]

cg20fall$bbch <- gsub("\\.", "\\,", cg20fall$bbch)
cg20fall <- separate(data = cg20fall, col = bbch, into = c("bbch", "bbchmf"), sep = "\\,")

cg20fall$leafcolor <- ifelse(cg20fall$bbch%in%c(90:99), cg20fall$doy, NA)

cg20fall$spindplot <- paste(cg20fall$spp, cg20fall$site, cg20fall$ind, cg20fall$plot)
cg20fall$lcolor <- NA
cg20fall <- cg20fall[!is.na(cg20fall$leafcolor),]
for(i in c(unique(cg20fall$spindplot))){ 
  
  lcolor <- cg20fall$leafcolor[i==cg20fall$spindplot][1]
  cg20fall$lcolor[i==cg20fall$spindplot] <- lcolor
  
}
cg20fall$leafcolor <- cg20fall$lcolor

cg20fall <- subset(cg20fall, select=c("spp", "year", "site", "ind", "leafcolor", "plot"))
cg20fall$plot <- as.character(cg20fall$plot)

cg20fall <- cg20fall[!duplicated(cg20fall),]
cg20clean <- full_join(cg20clean, cg20fall)



cg <- full_join(cg, cg20clean)
cg$Plot <- NULL
cg$risk <- NULL

vec <- c("ALNINC", "BETALL", "BETPAP", "BETPOP")
# get only 4 species
cg <- subset(cg, spp %in% vec)

lookupbb      <- aggregate(cg[c("budburst")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookupflowers <- aggregate(cg[c("flowers")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookuplo      <- aggregate(cg[c("leafout")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookupfbud    <- aggregate(cg[c("flobuds")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookupfbb     <- aggregate(cg[c("flobudburst")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookupfru     <- aggregate(cg[c("fruit")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookupripe    <- aggregate(cg[c("ripefruit")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookupbset    <- aggregate(cg[c("budset")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)
lookuplcolor  <- aggregate(cg[c("leafcolor")], cg[c("spp", "year", "site", "ind", "plot")], FUN=mean, na.rm=TRUE)

lookupbb      <- lookupbb     [which(lookupbb     $budburst != "NaN"),]
lookupflowers <- lookupflowers[which(lookupflowers$flowers != "NaN"),]
lookuplo      <- lookuplo     [which(lookuplo     $leafout != "NaN"),]
lookupfbud    <- lookupfbud   [which(lookupfbud   $flobuds != "NaN"),]
lookupfbb     <- lookupfbb    [which(lookupfbb    $flobudburst != "NaN"),]
lookupfru     <- lookupfru    [which(lookupfru    $fruit != "NaN"),]
lookupripe    <- lookupripe   [which(lookupripe   $ripefruit != "NaN"),]
lookupbset    <- lookupbset   [which(lookupbset   $budset != "NaN"),]
lookuplcolor  <- lookuplcolor [which(lookuplcolor $leafcolor != "NaN"),]

cgclean <- merge(lookupbb, lookupflowers, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookuplo, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookupfbud, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookupfbb, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookupfru, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookupripe, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookupbset, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)
cgclean <- merge(cgclean, lookuplcolor, by=c("spp", "year", "site", "ind", "plot"), all.x=TRUE, all.y=TRUE)

nrow(cgclean[which(!is.na(cgclean$leafout) & !is.nan(cgclean$leafout)),])
nrow(cgclean[which(!is.na(cgclean$budset) & !is.nan(cgclean$budset)),])

cgclean$provenance.lat <- NA
cgclean$provenance.long <- NA

cgclean$provenance.lat <- ifelse(cgclean$site == "HF", 42.531705, cgclean$provenance.lat)
cgclean$provenance.long <- ifelse(cgclean$site == "HF", -72.189920, cgclean$provenance.long)
cgclean$provenance.lat <- ifelse(cgclean$site == "WM", 44.112337, cgclean$provenance.lat)
cgclean$provenance.long <- ifelse(cgclean$site == "WM", -71.230138, cgclean$provenance.long)
cgclean$provenance.lat <- ifelse(cgclean$site == "GR", 44.794942, cgclean$provenance.lat)
cgclean$provenance.long <- ifelse(cgclean$site == "GR", -71.146683, cgclean$provenance.long)
cgclean$provenance.lat <- ifelse(cgclean$site == "SH", 45.932675, cgclean$provenance.lat)
cgclean$provenance.long <- ifelse(cgclean$site == "SH", -74.025070, cgclean$provenance.long)


# Match the style of Names I am using:
cgclean$name <- paste0(cgclean$spp, "_", cgclean$site, cgclean$ind, "_P", cgclean$plot)

#write.csv(cgclean, file="~/Documents/git/wildhellgarden/analyses/output/clean_obs_allyrs.csv", row.names=FALSE)

if(FALSE){
cgclean$dvr <- cgclean$leafout - cgclean$budburst

foo <- subset(cgclean, select=c("spp", "year", "site", "ind", "plot", "dvr"))
foo <- foo[!duplicated(foo),]
foo <- foo[complete.cases(foo),]

#moddvr <- rstanarm::stan_glmer(dvr ~ as.factor(year) + (as.factor(year) | spp/site), data=foo)

}

#### Prepare to join to source file ####
dtemp <- cgclean 
str(dtemp)

dtemp3 <- dtemp[, c(1:2, 6:ncol(dtemp))]
# remove spp column because I will add it in the merge file
dtemp4 <- dtemp3[, names(dtemp3) != "spp"]

# rename df 
obsdata <- dtemp4
# add spp name
obsdata$spp <- sub("_.*", "", obsdata$name)
names(obsdata)[names(obsdata) == "name"] <- "treeid"

write.csv(obsdata, file="output/obsData.csv", row.names=FALSE)




