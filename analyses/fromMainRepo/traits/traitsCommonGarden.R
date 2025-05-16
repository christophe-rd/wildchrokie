# started October 4, 2022

# aim of this code is to test whether there are differences in traits across the populations

rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
options(mc.cores = parallel::detectCores())
#graphics.off()

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/wildhellgarden/analyses")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/projects/others/deirdre/synchrony")
} else{
  setwd("/home/deirdre/WildHellCommonGarden") # for midge
}

# Load libraries
library(gridExtra)
library(ggplot2)
library(viridis)
library(rstan)
library(dplyr)
require(stringr)

# read in height data
ht <- read.csv("input/height_data.csv")

temp <- str_split_fixed(ht$ID, "_", 4); head(temp)
ht$spp<- temp[,1]
ht$pop <- temp[,2]
ht$indiv <- temp[,3]
ht$rep <- temp[,4]


ht$spp[ht$spp == "*SPITOM"] <- "SPITOM"
ht$spp[ht$spp == "BETPOPX"] <- "BETPOP"


unique(ht$Year) # 2019, 2021, 2022
htData <- ht[complete.cases(ht$Height),]

htData <- subset(htData, pop != "X")
unique(htData$pop)
# 2019
htData2019 <- subset(htData, Year == "2019")
htData2019$sppFact <- as.numeric(as.factor(htData2019$spp))
htData2019$popFact <- as.numeric(as.factor(htData2019$pop))
datalistHt <- with(htData2019,
                        list(y = Height,
                             sp = sppFact,
                             pop = popFact,
                             N = nrow(htData2019),
                             n_sp = length(unique(htData2019$spp)),
                             n_pop = length(unique(htData2019$pop))
                        )
)

htData2021 <- subset(htData, Year == "2021")
htData2021$sppFact <- as.numeric(as.factor(htData2021$spp))
htData2021$popFact <- as.numeric(as.factor(htData2021$pop))

# 2022
htData2022 <- subset(htData, Year == "2022")
htData2022$sppFact <- as.numeric(as.factor(htData2022$spp))
htData2022$popFact <- as.numeric(as.factor(htData2022$pop))
datalistHt <- with(htData2022, 
                   list(y = Height,  
                        sp = sppFact,
                        pop = popFact,
                        N = nrow(htData2022),
                        n_sp = length(unique(htData2022$spp)),
                        n_pop = length(unique(htData2022$pop))
                   )
)

#########################################################################
# Ht full model:
htFull <- htData 
htFull$sppFact <- as.numeric(as.factor(htFull$spp))
htFull$popFact <- as.numeric(as.factor(htFull$pop))
datalistHt <- with(htFull, 
                   list(y = Height,  
                        sp = sppFact,
                        pop = popFact,
                        N = nrow(htFull),
                        n_sp = length(unique(htFull$spp)),
                        n_pop = length(unique(htFull$pop))
                   )
)


mdlHtFull = stan('stan/trait_3levelwpop.stan',
                 data = datalistHt, iter = 4000, warmup= 3000,
                 chains = 4,
                 #control=list(adapt_delta=0.999,max_treedepth = 15)
                 #adapt_delta=0.999,
                 
)

save(mdlHtFull, file="output/mdlHtFull.Rda")

#load("analyses/output/mdlHt2019.Rda")
sumHtFull <- summary(mdlHtFull)$summary


range(sumHtFull[, "n_eff"])
range(sumHtFull[, "Rhat"])

plot(mdlHtFull, pars = c("mu_grand", "mu_sp", "mu_pop"))
#########################################################################
# Ht full year 13 spp model:
htData2022 <- subset(htData, Year == "2022")
spp22 <- unique(htData2022$spp)
htFull14 <- htData[htData$spp %in% spp22,] 
htFull14$sppFact <- as.numeric(as.factor(htFull14$spp))
htFull14$popFact <- as.numeric(as.factor(htFull14$pop))

datalistHt14 <- with(htFull14, 
                   list(y = Height,  
                        sp = sppFact,
                        pop = popFact,
                        N = nrow(htFull14),
                        n_sp = length(unique(htFull14$spp)),
                        n_pop = length(unique(htFull14$pop))
                   )
)

mdlHtFull14 = stan('stan/trait_3levelwpop.stan',
                 data = datalistHt14, iter = 5000, warmup= 4000,
                 chains=4,
                 control=list(adapt_delta=0.999,max_treedepth = 15)
                 #adapt_delta=0.999,
                 
)

save(mdlHtFull14, file="output/mdlHtFull14.Rda")
# 
# #load("analyses/output/mdlHt2019.Rda")
sumHt14 <- summary(mdlHtFull14)$summary

ssm <-  as.shinystan(mdlHt)
launch_shinystan(ssm)
# 
# range(sumHt19[, "n_eff"])
# range(sumHt19[, "Rhat"])

#########################################################################
# read in leaf trait data:
leaf <- read.csv("input/SLA_LDMC_trait.csv")

colnames(leaf)[colnames(leaf) == "species"] <- "spp"
colnames(leaf)[colnames(leaf) == "site"] <- "pop"
colnames(leaf)[colnames(leaf) == "ind"] <- "indiv"

######### SLA #################################
slaData <- leaf[complete.cases(leaf$SLA),]
unique(slaData$year)

slaData <- subset(slaData, pop != "XX")

#############################################################
slaFull <- slaData
slaFull$sppFact <- as.numeric(as.factor(slaFull$spp))
slaFull$popFact <- as.numeric(as.factor(slaFull$pop))
datalistSlaFull <- with(slaFull, 
                    list(y = SLA,  
                         sp = sppFact,
                         pop = popFact,
                         N = nrow(slaFull),
                         n_sp = length(unique(slaFull$spp)),
                         n_pop = length(unique(slaFull$pop)),
                         prior_mu_grand_mu = 0,
                         prior_mu_grand_sigma = 20,
                         prior_sigma_sp_mu = 0,
                         prior_sigma_sp_sigma = 10,
                         prior_sigma_pop_mu = 0,
                         prior_sigma_pop_sigma = 10,
                         prior_sigma_y_mu = 0,
                         prior_sigma_y_sigma = 10
                    )
)


mdlSLAFull = stan('stan/trait_3levelwpop_ldmc.stan',
                   data = datalistSlaFull, iter = 8000, warmup= 7000,
                   chains=4,
                   control=list(adapt_delta=0.999,max_treedepth = 15)
                   #adapt_delta=0.999,
                   
)

save(mdlSLAFull, file="output/mdlSLAFull.Rda")

######################################
slaData2022 <- subset(slaData, year == "2022")
spp22 <- unique(slaData2022$spp)

slaFull14 <- slaFull[slaFull$spp %in% spp22,] 
slaFull14$sppFact <- as.numeric(as.factor(slaFull14$spp))
slaFull14$popFact <- as.numeric(as.factor(slaFull14$pop))
datalistSlaFull <- with(slaFull14, 
                        list(y = SLA,  
                             sp = sppFact,
                             pop = popFact,
                             N = nrow(slaFull14),
                             n_sp = length(unique(slaFull14$spp)),
                             n_pop = length(unique(slaFull14$pop)),
                             prior_mu_grand_mu = 0,
                             prior_mu_grand_sigma = 20,
                             prior_sigma_sp_mu = 0,
                             prior_sigma_sp_sigma = 10,
                             prior_sigma_pop_mu = 0,
                             prior_sigma_pop_sigma = 10,
                             prior_sigma_y_mu = 0,
                             prior_sigma_y_sigma = 10
                        )
)


mdlSLAFull14 = stan('stan/trait_3levelwpop_ldmc.stan',
                  data = datalistSlaFull, iter = 8000, warmup= 7000,
                  chains=4,
                  control=list(adapt_delta=0.999, max_treedepth = 15)
)
save(mdlSLAFull14, file="output/mdlSLAFull14.Rda")
# runs with no issues

# ############### LDMC ##############################
ldmcData <- leaf[complete.cases(leaf$LDMC),]
unique(ldmcData$year)
ldmcData <- subset(ldmcData, pop != "XX")
unique(ldmcData$pop)
# 2019
ldmcData2019 <- subset(ldmcData, year == "2019")

#2022
ldmcData2022 <- subset(ldmcData, year == "2022")

#############################################################
# Full ldmc mdl:
ldmcDataFull <- ldmcData
ldmcDataFull$sppFact <- as.numeric(as.factor(ldmcDataFull$spp)) #18 spp
ldmcDataFull$popFact <- as.numeric(as.factor(ldmcDataFull$pop))

datalistLdmc <- with(ldmcDataFull, 
                     list(y = LDMC,  
                          sp = sppFact,
                          pop = popFact,
                          N = nrow(ldmcDataFull),
                          n_sp = length(unique(ldmcDataFull$spp)),
                          n_pop = length(unique(ldmcDataFull$pop)),
                          prior_mu_grand_mu = 300,
                          prior_mu_grand_sigma = 50,
                          prior_sigma_sp_mu = 50,
                          prior_sigma_sp_sigma = 20,
                          prior_sigma_pop_mu = 0,
                          prior_sigma_pop_sigma = 10,
                          prior_sigma_y_mu = 50,
                          prior_sigma_y_sigma = 10
                     )
)

mdlLDMC = stan('stan/trait_3levelwpop_ldmc.stan',
               data = datalistLdmc, iter = 8000, warmup = 7000, chains = 4,
               control = list(max_treedepth = 15,adapt_delta = 0.99))

ldmcData2022 <- subset(ldmcData, year == "2022")
spp22 <- unique(ldmcData2022$spp)
ldmcDataFull14 <- ldmcData[ldmcData$spp %in% spp22, ]
ldmcDataFull14$sppFact <- as.numeric(as.factor(ldmcDataFull14$spp))
ldmcDataFull14$popFact <- as.numeric(as.factor(ldmcDataFull14$pop))
datalistLdmc14<- with(ldmcDataFull14, 
                     list(y = LDMC,  
                          sp = sppFact,
                          pop = popFact,
                          N = nrow(ldmcDataFull14),
                          n_sp = length(unique(ldmcDataFull14$spp)),
                          n_pop = length(unique(ldmcDataFull14$pop)),
                          prior_mu_grand_mu = 300,
                          prior_mu_grand_sigma = 50,
                          prior_sigma_sp_mu = 50,
                          prior_sigma_sp_sigma = 20,
                          prior_sigma_pop_mu = 0,
                          prior_sigma_pop_sigma = 10,
                          prior_sigma_y_mu = 50,
                          prior_sigma_y_sigma = 10
                     )
)


mdlLDMCFull14 = stan('stan/trait_3levelwpop_ldmc.stan',
                   data = datalistLdmc14, iter = 8000, warmup = 7000,
                   chains=4, control=list(adapt_delta = 0.99, max_treedepth = 15)
)
 save(mdlLDMCFull14, file="output/mdlLDMCFull14.Rda")
################# SSD ############################################

# read in stem trait data: only one year of data
stem <- read.csv("input/SSD_trait.csv")

ssdData <- stem[complete.cases(stem$ssd),]
ssdData <- subset(ssdData, site != "XX")

ssdData$sppFact <- as.numeric(as.factor(ssdData$species))
ssdData$popFact <- as.numeric(as.factor(ssdData$site))

datalistSSD <- with(ssdData, 
     list(y = ssd,  
          sp = sppFact,
          pop = popFact,
          N = nrow(ssdData),
          n_sp = length(unique(ssdData$spp)),
          n_pop = length(unique(ssdData$pop)),
          prior_mu_grand_mu = 0,
          prior_mu_grand_sigma = 1,
          prior_sigma_sp_mu = 0,
          prior_sigma_sp_sigma = 1,
          prior_sigma_pop_mu = 0,
          prior_sigma_pop_sigma = 1,
          prior_sigma_y_mu = 0,
          prior_sigma_y_sigma = 10
     )
)

mdlSSD = stan('stan/trait_3levelwpop.stan',
              data = datalistSSD, iter = 4000, warmup= 3000,
              chains=4,
              control=list(adapt_delta=0.9999,
              max_treedepth = 19)
)

mdlSSDSmall = stan('stan/trait_3levelwpop_ldmc.stan',
              data = datalistSSD, iter = 4000, warmup= 3000,
              chains=4,
              control=list(adapt_delta=0.99,
                           max_treedepth = 15)
)
# 

#############################################
# calculate relative growth rates:
rgr <- merge(htData2019, htData2021, by = c("Plot", "spp", "pop", "indiv"), all = T)
rgr <- rgr[,c("Plot", "spp", "pop", "indiv", "ID.x", "ID.y", "Height.x", "Height.y")]
names(rgr) <- c("Plot", "spp", "pop", "indiv", "ID19", "ID21", "ht19", "ht21")

rgr <- merge(rgr, htData2022, by = c("Plot", "spp", "pop", "indiv"), all = T)
rgr <- rgr[,c("Plot", "spp", "pop", "indiv", "ID19", "ID21","ID", "ht19", "ht21", "Height")]
names(rgr) <- c("Plot", "spp", "pop", "indiv", "ID19", "ID21", "ID22","ht19", "ht21", "ht22")

rgr$rgr21 <- rgr$ht21 - rgr$ht19
rgr$rgr22 <- rgr$ht22 - rgr$ht21
rgr$rgrtot <- (rgr$ht22 - rgr$ht19)/3
rgr$grtot <- (rgr$ht22 - rgr$ht19)

# there are a concerning number of rgr that are very negative

# 1. for now removing all the really negative values:
rgr21 <- subset(rgr, rgr21 >-0.05)
rgr$rgr21[which(rgr$rgr21 < 0)] <- "0"

rgr22 <- subset(rgr, rgr22 > -0.05)

tot <- subset(rgr, rgrtot > -0.05)
small <- subset(tot, rgrtot < 0)
rest <- subset(tot, rgrtot > 0)

small$rgrtot <- 0

relGrow <- rbind(small, rest)

relGrow$sppFact <- as.numeric(as.factor(relGrow$spp))
relGrow$popFact <- as.numeric(as.factor(relGrow$pop))

datalistRGR <- with(relGrow, 
                     list(y = rgrtot,  
                          sp = sppFact,
                          pop = popFact,
                          N = nrow(relGrow),
                          n_sp = length(unique(relGrow$spp)),
                          n_pop = length(unique(relGrow$pop))
                     )
)


pdf("postPrior.pdf", height = 10, width = 10)
par(mfrow = c(2,2))
h1 <- hist(rnorm(1000, 0,20), col=rgb(0,0,1,1/4))
hist(post$mu_grand,  col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000,0,10), col=rgb(0,0,1,1/4))
hist(post$sigma_sp, add = T,  col=rgb(1,0,1,1/4))

h1 <- hist(rnorm(1000, 0,10), col=rgb(0,0,1,1/4))
hist(post$sigma_pop,  col=rgb(1,0,1,1/4), add =T)

h1 <- hist(rnorm(1000, 0,10), col=rgb(0,0,1,1/4))
hist(post$sigma_y, col=rgb(1,0,1,1/4), add = T)

hist(post$mu_grand,  col=rgb(1,0,1,1/4))

###########################################
# Compare the 14 spp model to the full models and see if different
load("output/mdlLDMCFull.Rda")
load("output/mdlLDMCFull14.Rda")

load("output/mdlSLAFull.Rda")
load("output/mdlSLAFull14.Rda")

