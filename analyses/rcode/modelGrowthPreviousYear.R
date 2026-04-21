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

runmodels <- F
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
data <- list(
  y = log(rw$lengthMM),
  N = nrow(rw),
  Nspp = length(unique(rw$spp_num)),
  Nsite = length(unique(rw$site_num)),
  site = as.numeric(as.character(rw$site_num)),
  species = as.numeric(as.character(rw$spp_num)),
  treeid = as.numeric(rw$treeid_num),
  Ntreeid = length(unique(as.numeric(rw$treeid_num))),
  gdd = rw$gddcurrentyr / 200,
  gddyr = rw$gddpreviousyr / 200
  # gdd = (rw$gddcurrentyr - mean(rw$gddcurrentyr)) / sd(rw$gddcurrentyr),
  # gddyr = (rw$gddpreviousyr - mean(rw$gddpreviousyr)) / sd(rw$gddpreviousyr)
)

# Fit model GDD 
if(runmodels) {
gddmodel <- stan_model("stan/modelGrowthPreviousYear.stan")
fit <- sampling(gddmodel, data = data, iter = 2000, chains = 4)
saveRDS(fit, "output/stanOutput/fitGrowthPreviousYear")
diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)
}

library(ggplot2)
ggplot(rw) + 
  geom_smooth(aes(x = gddcurrentyr, y = lengthMM)) + 
  geom_point(aes(x = gddcurrentyr, y = lengthMM)) 
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

# discs by species
jpeg(
  filename = "figures/growthPreviousYearModel/retrodictiveDiskSpp.jpeg",
  width = 3600, height = 2000, res = 300)
par(mfrow = c(1,data$Nspp))
for (s in unique(data$species)) { # s = 1
  idxs <- which(data$species == s)
  util$plot_disc_pushforward_quantiles(samples,
                                       paste0("y_rep[", idxs, "]"),
                                       baseline_values = data$y[idxs],
                                       ylab = "Leafout",
                                       main = paste("Spp", s))
}
dev.off()
# discs by year
par(mfrow = c(1,data$Nyear))
for (y in unique(data$year)) { # s = 1
  idxs <- which(data$year == y)
  util$plot_disc_pushforward_quantiles(samples,
                                       paste0("y_rep[", idxs, "]"),
                                       baseline_values = data$y[idxs],
                                       ylab = "Leafout",
                                       main = paste("Year", y))
}


##### Plot posterior vs priors for gdd fit #####
# pdf(file = "figures/empiricalData/gddModelPriorVSPosteriorPrvsYr.pdf", width = 8, height = 10)
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns) & !grepl("yr", columns)]]
bsppyr_df <- df_fit[, columns[grepl("bspyr", columns)]]
# treeid_df <- df_fit[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
site_df <- df_fit[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(bsppyr_df) <- 1:ncol(bsppyr_df)
# colnames(treeid_df) <- 1:ncol(treeid_df)
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
plot(bspp_df2_current$mean, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_spp + 0.5), 
     xlab = "slope current year", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_current$p5,  y_pos, bspp_df2_current$p95, y_pos,
         col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_current$p25, y_pos, bspp_df2_current$p75, y_pos,
         col = wccolslatbi, lwd = 3)
mtext("Current year", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: Previous year
plot(bspp_df2_previous$mean, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_spp + 0.5),
     xlab = "slope previous year", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = FALSE,      
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_previous$p5,  y_pos, bspp_df2_previous$p95, y_pos,
         col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_previous$p25, y_pos, bspp_df2_previous$p75, y_pos,
         col = wccolslatbi, lwd = 3)
mtext("Previous year", side = 3, adj = 0, font = 2, cex = 0.9)
dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Simulated data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
set.seed(124)
Nspp <- 10
Nsite <- 4
N_tree_per_spp <- sample(10:50, Nspp, replace = TRUE)
Nyear <- 4
Ntreeid <- sum(N_tree_per_spp)
N <- Ntreeid * Nyear

# index vectors
tree_species_idxs <- rep(1:Nspp, times = N_tree_per_spp)
tree_site_idxs <- rep(sample(1:Nsite, Ntreeid, replace = TRUE))
treeid <- rep(1:Ntreeid, each = Nyear)
species <- tree_species_idxs[treeid]
site <- tree_site_idxs[treeid]
year <- rep(1:Nyear, times = Ntreeid)

# sigma params
sigma_y <- 1
sigma_aspp <- 3
sigma_asite <- 1
sigma_bsp <- 2
sigma_bspyr <- 2

# simulate raw gdd by year, then standardize across all obs
gddcurrentyr_yr <- rnorm(Nyear, 1200, 150)
gddpreviousyr_yr <- rnorm(Nyear, 1200, 150)
gddcurrentyr_raw <- gddcurrentyr_yr[year]
gddpreviousyr_raw <- gddpreviousyr_yr[year]
gdd <- (gddcurrentyr_raw - mean(gddcurrentyr_raw)) / sd(gddcurrentyr_raw)
gddyr <- (gddpreviousyr_raw - mean(gddpreviousyr_raw)) / sd(gddpreviousyr_raw)

# sim parameters
a <- 2
aspp <- rnorm(Nspp, 0, sigma_aspp)
asite <- rnorm(Nsite, 0, sigma_asite)
bsp <- rnorm(Nspp, 0, sigma_bsp)
bspyr <- rnorm(Nspp, 0, sigma_bspyr)

# df
sim <- data.frame(
  treeid = treeid,
  species = species,
  site = site,
  year = year,
  gdd = gdd,
  gddyr = gddyr
)
sim$a <- a
sim$aspp <- aspp[sim$species]
sim$asite <- asite[sim$site]
sim$bsp <- bsp[sim$species]
sim$bspyr <- bspyr[sim$species]

mu <- sim$a + sim$aspp + sim$asite + sim$bsp * sim$gdd + sim$bspyr * sim$gddyr
sim$y <- rnorm(N, mu, sigma_y)

# data list
simd <- list(
  y = sim$y,
  N = nrow(sim),
  Nspp = Nspp,
  Nsite = Nsite,
  species = sim$species,
  site = sim$site,
  treeid = sim$treeid,
  Ntreeid = Ntreeid,
  gdd = sim$gdd,
  gddyr = sim$gddyr
)

# fit sim data
gddmodel <- stan_model("stan/modelGrowthPreviousYear.stan")
fitsim <- sampling(gddmodel, data = simd, iter = 2000, chains = 4)
saveRDS(fitsim, "output/stanOutput/simGrowthPreviousYear")
diagnostics <- util$extract_hmc_diagnostics(fitsim) 
util$check_all_hmc_diagnostics(diagnostics)

##### Prior vs posterior #####
pdf(file = "figures/growthPreviousYearModel/simPrvsYr.pdf", width = 10, height = 10)
pal <- wes_palette("AsteroidCity1")[3:4]
par(mfrow = c(2, 3))

df_fitsim <- as.data.frame(fitsim)

columns <- colnames(df_fitsim)[!grepl("prior", colnames(df_fitsim))]
aspp_df <- df_fitsim[, columns[grepl("aspp", columns)]]
asite_df <- df_fitsim[, columns[grepl("asite", columns)]]
bsp_df <- df_fitsim[, columns[grepl("^bsp\\[", columns)]]
bspyr_df <- df_fitsim[, columns[grepl("^bspyr\\[", columns)]]

# a
plot(density(df_fitsim[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.1))
lines(density(df_fitsim[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitsim[, "sigma_y_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_sigma_y",
     xlab = "sigma_y", ylim = c(0, 2))
lines(density(df_fitsim[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitsim[, "aspp_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_aspp",
     xlab = "aspp", xlim = c(-50, 50), ylim = c(0, 0.2))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitsim[, "asite_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_asite",
     xlab = "asite", xlim = c(-30, 30), ylim = c(0, 0.1))
for (col in colnames(asite_df)) {
  lines(density(asite_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitsim[, "bsp_prior"]),
     col = pal[1], lwd = 2, xlim = c(-30, 30),
     main = "priorVSposterior_bsp",
     xlab = "bsp", ylim = c(0, 1))
for (col in colnames(bsp_df)) {
  lines(density(bsp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bspyr
plot(density(df_fitsim[, "bspyr_prior"]),
     col = pal[1], lwd = 2, xlim = c(-30, 30),
     main = "priorVSposterior_bspyr",
     xlab = "bspyr", ylim = c(0, 1))
for (col in colnames(bspyr_df)) {
  lines(density(bspyr_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()


##### Check parameter recovery #####
df_fitsim <- as.data.frame(fitsim)

aspp_df2 <- extract_params(df_fitsim, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
asite_df2 <- extract_params(df_fitsim, "asite", "fit_asite", "site", "asite\\[(\\d+)\\]")
bsp_df2 <- extract_params(df_fitsim, "bsp", "fit_bsp", "spp", "bsp\\[(\\d+)\\]")
bsp_df2 <- bsp_df2[!grepl("yr", bsp_df2$spp),]
bspyr_df2 <- extract_params(df_fitsim, "bspyr", "fit_bspyr", "spp", "bspyr\\[(\\d+)\\]")

aspp_df2$sim_aspp <- sim$asp[match(aspp_df2$spp, sim$species)]
asite_df2$sim_asite <- sim$asite[match(asite_df2$site, sim$site)]
bsp_df2$sim_bsp <- sim$bsp[match(bsp_df2$spp, sim$species)]
bspyr_df2$sim_bspyr <- sim$bspyr[match(bsp_df2$spp, sim$species)]

jpeg("figures/growthPreviousYearModel/simVSfit.jpeg", width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 4))

# aspp
plot(aspp_df2$sim_aspp, aspp_df2$fit_aspp,
     xlab = "sim", ylab = "fit", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2$fit_aspp_per5, aspp_df2$fit_aspp_per95)),
     xlim = range(aspp_df2$sim_aspp))
arrows(x0 = aspp_df2$sim_aspp, y0 = aspp_df2$fit_aspp_per5,
       x1 = aspp_df2$sim_aspp, y1 = aspp_df2$fit_aspp_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2$sim_aspp, aspp_df2$fit_aspp, pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# asite
plot(asite_df2$sim_asite, asite_df2$fit_asite,
     xlab = "sim", ylab = "fit", main = "asite", type = "n", frame = FALSE,
     ylim = range(c(asite_df2$fit_asite_per5, asite_df2$fit_asite_per95)),
     xlim = range(asite_df2$sim_asite))
arrows(x0 = asite_df2$sim_asite, y0 = asite_df2$fit_asite_per5,
       x1 = asite_df2$sim_asite, y1 = asite_df2$fit_asite_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(asite_df2$sim_asite, asite_df2$fit_asite, pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# bsp
plot(bsp_df2$sim_bsp, bsp_df2$fit_bsp,
     xlab = "sim", ylab = "fit", main = "bsp", type = "n", frame = FALSE,
     ylim = range(c(bsp_df2$fit_bsp_per5, bsp_df2$fit_bsp_per95)),
     xlim = range(bsp_df2$sim_bsp))
arrows(x0 = bsp_df2$sim_bsp, y0 = bsp_df2$fit_bsp_per5,
       x1 = bsp_df2$sim_bsp, y1 = bsp_df2$fit_bsp_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bsp_df2$sim_bsp, bsp_df2$fit_bsp, pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# bspyr
plot(bspyr_df2$sim_bspyr, bspyr_df2$fit_bspyr,
     xlab = "sim", ylab = "fit", main = "bspyr", type = "n", frame = FALSE,
     ylim = range(c(bspyr_df2$fit_bspyr_per5, bspyr_df2$fit_bspyr_per95)),
     xlim = range(bspyr_df2$sim_bspyr))
arrows(x0 = bspyr_df2$sim_bspyr, y0 = bspyr_df2$fit_bspyr_per5,
       x1 = bspyr_df2$sim_bspyr, y1 = bspyr_df2$fit_bspyr_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspyr_df2$sim_bspyr, bspyr_df2$fit_bspyr, pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

dev.off()
