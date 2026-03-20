# Wildchrokie model of previous years condition on following year's growth
# CRD 19 March 2025

# Goal: check if the condition of the previous year on current year's growth

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)

# Load library 
library(ggplot2)
library(rstan)
library(wesanderson)

if (length(grep("christophe_rouleau-desrochers", getwd())) > 0) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if (length(grep("lizzie", getwd())) > 0) {
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
} else  {
  setwd("/home/crouleau/wildchrokie/analyses")
}

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)
# my function to extract parameters
source('rcode/utilExtractParam.R')

emp <- read.csv("output/empiricalDataMAIN.csv")
rw <- read.csv("output/wildchrokieRingWidth.csv")
gdd <- read.csv("/Users/christophe_rouleau-desrochers/github/coringtreespotters/analyses/output/gddByYear.csv")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# calculate gdd by year between 1st april (91) and 1st august (213)  
gddsub <- subset(gdd, doy >90 & doy < 214)
agg <- aggregate(GDD_5 ~ year, gddsub, FUN = min)
gddsub$mingdd <- agg$GDD_5[match(gddsub$year, agg$year)]
gddsub$diff <- gddsub$GDD_5 - gddsub$mingdd
agg2 <- aggregate(diff ~ year, gddsub, FUN = max)

years <- unique(rw$year)
rw$yeardiff <- NA

for (i in years) {
  rw$yeardiff[rw$year == i] <- rw$year[rw$year == i] - 1
}

# add gdd of previous year indexed by yeardiff
rw$gdd <- agg2$diff[match(rw$yeardiff, agg2$year)]

rw$lengthMM <- rw$lengthCM * 10

# transform my groups to numeric values
rw$site_num <- match(rw$site, unique(rw$site))
rw$spp_num <- match(rw$spp, unique(rw$spp))
rw$treeid_num <- match(rw$treeid, unique(rw$treeid))

# transform data in vectors for GDD
rw <- rw[!is.na(rw$lengthMM) & !is.na(rw$gdd),]
rw <- subset(rw, site != "XX" & !spp %in% c("QUERNA", "QUEROB"))
y <- rw$lengthMM
N <- nrow(rw)
Nspp <- length(unique(rw$spp_num))
Nsite <- length(unique(rw$site_num))
site <- as.numeric(as.character(rw$site_num))
species <- as.numeric(as.character(rw$spp_num))
treeid <- as.numeric(rw$treeid_num)
Ntreeid <- length(unique(treeid))
gdd <- rw$gdd

# Fit model GDD
rstan_options(auto_write = TRUE)
gddmodel <- stan_model("stan/twolevelhierint.stan")
fit <- sampling(gddmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gdd"),
                   warmup = 1000, iter=2000, chains=4)
saveRDS(fit, "output/stanOutput/fitGrowthPreviousYear")
