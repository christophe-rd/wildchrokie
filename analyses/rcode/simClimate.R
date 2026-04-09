# Wildchrokie model
# CRD 8 April 2026
# We want to figure out if fitting the nested model returns the same thing

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
library(patchwork)

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
source('rcode/utilExtractParam.R')

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Simulate data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
set.seed(124)

Nspp <- 10
N_tree_per_spp <- sample(1:15, Nspp)
Nyear <- 3
Ntreeid <- sum(N_tree_per_spp)
N <- Ntreeid * Nyear

tree_species_idxs <- rep(1:Nspp, times = N_tree_per_spp)

# index
treeid <- rep(1:Ntreeid, each = Nyear)
species <- tree_species_idxs[treeid]
year <- rep(1:Nyear, times = Ntreeid)

sigma_y <- 4
sigma_aspp <- 3
sigma_bspp <- 1.5
sigma_ayear <- 3.5

# sim climate variable
climpredictoryr <- rnorm(Nyear, 10, 0.8)
climpredictor <- climpredictoryr[year]

# sim intercepts
a <- 2
aspp <- rnorm(Nspp, 0, sigma_aspp)
bspp <- rnorm(Nspp, 0, sigma_bspp)
ayear <- rnorm(Nyear, 0, sigma_ayear)

# df
sim <- data.frame( 
  treeid = treeid,
  species = species,
  year = year)

sim$a <- 2
sim$aspp <- aspp[sim$species]
sim$bspp <- bspp[sim$species]
sim$ayear <- ayear[sim$year]

y <- sim$a + sim$aspp + sim$ayear + (sim$bspp * climpredictor[year])

y <- rnorm(N, y, sigma_y)

hist(y)
fit <- stan("stan/TSclimatePredictors.stan", 
                  data = c("N","y",
                         "Nspp","species",
                         "Nyear", "year",
                         "climpredictor"),
                  warmup = 1000, iter = 2000, chains=4)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot posterior vs priors for gdd fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
pdf(file = "figures/climate/simClimVar.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(2, 2))
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
ayear_df <- df_fit[, columns[grepl("ayear", columns)]]

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
     xlab = "aspp", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fit[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-50, 50), ylim = c(0, 0.1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fit[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-10, 10),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Recover and plot parameters ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Nested #####
df_fit <- as.data.frame(fit)

# posterior summaries
sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
ayear_df2 <- extract_params(df_fit, "ayear", "fit_ayear", 
                            "year", "ayear\\[(\\d+)\\]")
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
bspp_df2   <- extract_params(df_fit, "bsp", "fit_bsp", 
                             "spp", "bsp\\[(\\d+)\\]")

# add sim coef
aspp_df2$sim_aspp <- aspp
ayear_df2$sim_ayear <- ayear
bspp_df2$sim_bspp <- bspp

jpeg("figures/climate/simVSfit.jpeg", width = 9, height = 6, units = "in", res = 300)
par(mfrow = c(1, 3))

# aspp
plot(aspp_df2$sim_aspp, aspp_df2$fit_aspp,
     xlab = "sim", ylab = "fit", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2$fit_aspp_per5, aspp_df2$fit_aspp_per95)),
     xlim = range(aspp_df2$sim_aspp))
arrows(x0 = aspp_df2$sim_aspp, y0 = aspp_df2$fit_aspp_per5,
       x1 = aspp_df2$sim_aspp, y1 = aspp_df2$fit_aspp_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2$sim_aspp, aspp_df2$fit_aspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# ayear
plot(ayear_df2$sim_ayear, ayear_df2$fit_ayear,
     xlab = "sim", ylab = "fit", main = "ayear", type = "n", frame = FALSE,
     ylim = range(c(ayear_df2$fit_ayear_per5, ayear_df2$fit_ayear_per95)),
     xlim = range(ayear_df2$sim_ayear))
arrows(x0 = ayear_df2$sim_ayear, y0 = ayear_df2$fit_ayear_per5,
       x1 = ayear_df2$sim_ayear, y1 = ayear_df2$fit_ayear_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(ayear_df2$sim_ayear, ayear_df2$fit_ayear,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# bspp
plot(bspp_df2$sim_bspp, bspp_df2$fit_bsp,
     xlab = "sim", ylab = "fit", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2$fit_bsp_per5, bspp_df2$fit_bsp_per95)),
     xlim = range(bspp_df2$sim_bspp))
arrows(x0 = bspp_df2$sim_bspp, y0 = bspp_df2$fit_bsp_per5,
       x1 = bspp_df2$sim_bspp, y1 = bspp_df2$fit_bsp_per95,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2$sim_bspp, bspp_df2$fit_bsp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)
dev.off()


