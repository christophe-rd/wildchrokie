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
source('rcode/tools.R')

emp <- read.csv("output/empiricalDataMAIN.csv")
rw <- read.csv("output/wildchrokieRingWidth.csv")
gdd <- read.csv("/Users/christophe_rouleau-desrochers/github/coringtreespotters/analyses/output/gddByYear.csv")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# calculate gdd by year between 1st May (121) and 1st august (244)  
gddsub <- subset(gdd, doy >121 & doy < 244)
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
y <- log(rw$lengthMM)
N <- nrow(rw)
Nspp <- length(unique(rw$spp_num))
Nsite <- length(unique(rw$site_num))
site <- as.numeric(as.character(rw$site_num))
species <- as.numeric(as.character(rw$spp_num))
treeid <- as.numeric(rw$treeid_num)
Ntreeid <- length(unique(treeid))
gdd	<- rw$gddcurrentyr / 200
gddyr <- rw$gddpreviousyr / 200

# Fit model GDD 
gddmodel <- stan_model("stan/modelGrowthPreviousYear.stan")
fit <- sampling(gddmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gdd", "gddyr"),
                   warmup = 1000, iter = 2000, chains = 4)
saveRDS(fit, "output/stanOutput/fitGrowthPreviousYear")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# fit <- readRDS("output/stanOutput/fitGrowthPreviousYear")
samples <- util$extract_expectand_vals(fit)
jpeg(
  filename = "figures/growthPreviousYearModel/retrodictiveCheckHistPrvsYr.jpeg", 
  width = 2400, height = 2400, res = 300          
)
util$plot_hist_quantiles(samples, "y_rep", 
                         -2, # lower x axis limit
                         5, # upper x axis limit
                         0.2, # binning
                         baseline_values = y,
                         xlab = "log(ring width")
dev.off()


##### Plot posterior vs priors for gdd fit #####
# pdf(file = "figures/empiricalData/gddModelPriorVSPosteriorPrvsYr.pdf", width = 8, height = 10)
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns) & !grepl("yr", columns)]]
bsppyr_df <- df_fit[, columns[grepl("bspyr", columns)]]
treeid_df <- df_fit[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
site_df <- df_fit[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(bsppyr_df) <- 1:ncol(bsppyr_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

jpeg("figures/growthPreviousYearModel/gddModelPriorVSPosteriorPrvsYr.jpeg", 
     width =2400, height = 3600, res =300)
pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 3))

# a
plot(density(df_fit[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.1))
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
     xlab = "aspp", xlim = c(-60, 60), ylim = c(0, 0.1))
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
  lines(density(bsppyr_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Compare bspp vs bsppyr ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Open device
# posterior summaries
bspp_df2_current <- extract_params(df_fit, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
bspp_df2_current <- subset(bspp_df2_current, spp %in% bspp_df2_current$spp[!grepl("yr", bspp_df2_current$spp)])
bspp_df2_previous <- extract_params(df_fit, "bspyr", "fit_bspp", "spp", "bspyr\\[(\\d+)\\]")

jpeg("figures/growthPreviousYearModel/bsppCurrentVSpreviousYR.jpeg", width = 6, height = 9, units = "in", res = 300)
par(mfrow = c(2,1))
n_spp <- length(unique(bspp_df2_current$spp))
y_pos <- rev(1:n_spp)

# Current year
plot(bspp_df2_current$fit_bspp, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_spp + 0.5), 
     xlab = "slope current year", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_current$fit_bspp_per5,  y_pos, bspp_df2_current$fit_bspp_per95, y_pos,
         col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_current$fit_bspp_per25, y_pos, bspp_df2_current$fit_bspp_per75, y_pos,
         col = wccolslatbi, lwd = 3)
mtext("Current year", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: Previous year
plot(bspp_df2_previous$fit_bspp, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_spp + 0.5),
     xlab = "slope previous year", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = FALSE,      
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_previous$fit_bspp_per5,  y_pos, bspp_df2_previous$fit_bspp_per95, y_pos,
         col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_previous$fit_bspp_per25, y_pos, bspp_df2_previous$fit_bspp_per75, y_pos,
         col = wccolslatbi, lwd = 3)
mtext("Previous year", side = 3, adj = 0, font = 2, cex = 0.9)
dev.off()


