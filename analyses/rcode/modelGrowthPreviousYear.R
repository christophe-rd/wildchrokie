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
rw$gddpreviousyr <- agg2$diff[match(rw$yeardiff, agg2$year)]
rw$gddcurrentyr <- agg2$diff[match(rw$year, agg2$year)]

rw$lengthMM <- rw$lengthCM * 10

# transform my groups to numeric values
rw$site_num <- match(rw$site, unique(rw$site))
rw$spp_num <- match(rw$spp, unique(rw$spp))
rw$treeid_num <- match(rw$treeid, unique(rw$treeid))

# transform data in vectors for GDD
rw <- rw[!is.na(rw$lengthMM) & !is.na(rw$gddcurrentyr) & !is.na(rw$gddpreviousyr),]
rw <- subset(rw, site != "XX" & !spp %in% c("QUERNA", "QUEROB"))
y <- rw$lengthMM
N <- nrow(rw)
Nspp <- length(unique(rw$spp_num))
Nsite <- length(unique(rw$site_num))
site <- as.numeric(as.character(rw$site_num))
species <- as.numeric(as.character(rw$spp_num))
treeid <- as.numeric(rw$treeid_num)
Ntreeid <- length(unique(treeid))
gdd	<- rw$gddcurrentyr
gddyr <- rw$gddpreviousyr

# Fit model GDD
rstan_options(auto_write = TRUE)
gddmodel <- stan_model("stan/modelGrowthPreviousYear.stan")
fit <- sampling(gddmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gdd", "gddyr"),
                   warmup = 1000, iter=2000, chains=4)
saveRDS(fit, "output/stanOutput/fitGrowthPreviousYear")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
samples <- util$extract_expectand_vals(fit)
jpeg(
  filename = "figures/troubleShootingGrowthModel/retrodictiveCheckHistPrvsYr.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
util$plot_hist_quantiles(samples, "y_rep", 
                         -5, # lower x axis limit
                         15, # upper x axis limit
                         0.5, # binning
                         baseline_values = y,
                         xlab = "Ring width (mm)")
dev.off()


##### Plot posterior vs priors for gdd fit #####
# pdf(file = "figures/empiricalData/gddModelPriorVSPosteriorPrvsYr.pdf", width = 8, height = 10)
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
treeid_df <- df_fit[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
site_df <- df_fit[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

jpeg("figures/empiricalData/gddModelPriorVSPosteriorPrvsYr.jpeg", 
     width =2400, height = 3600, res =300)
pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 3))

# a
plot(density(df_fit[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fit[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fit[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fit[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fit[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fit[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fit[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fit[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fit[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bspyr
plot(density(df_fit[, "bspyr_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bspPreviousYr", 
     xlab = "bspPreviousYr", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()
