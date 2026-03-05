# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
# options(mc.cores = parallel::detectCores())
options(digits = 3)
# quartz()

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(shinystan)
library(wesanderson)
library(patchwork)
library(rstanarm)

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

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
emp <- read.csv("output/empiricalDataMAIN.csv")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit model GDD
emp2 <- emp[!is.na(emp$pgsGDD5),]

# transform my groups to numeric values
emp2$site_num <- match(emp2$site, unique(emp2$site))
emp2$spp_num <- match(emp2$spp, unique(emp2$spp))
emp2$treeid_num <- match(emp2$treeid, unique(emp2$treeid))

# transform data in vectors for GDD
y <- emp2$lengthCM*10 # ring width in mm
N <- nrow(emp2)
Nspp <- length(unique(emp2$spp_num))
Nsite <- length(unique(emp2$site_num))
site <- as.numeric(as.character(emp2$site_num))
species <- as.numeric(as.character(emp2$spp_num))
treeid <- as.numeric(emp2$treeid_num)
Ntreeid <- length(unique(treeid))

# different response variables
gdd <- emp2$pgsGDD5/200
gsl <- as.numeric(emp2$pgsGSL/10)
sos <- emp2$leafout/10
eos <- emp2$budset/10

# Fit model GDD
rstan_options(auto_write = TRUE)
gddmodel <- stan_model("stan/twolevelhierint.stan")
fitgdd <- sampling(gddmodel, data = c("N","y",
                                "Nspp","species",
                                "Nsite", "site", 
                                "Ntreeid", "treeid", 
                                "gdd"),
                warmup = 1000, iter=2000, 
                chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDD")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fitgdd) 
util$check_all_hmc_diagnostics(diagnostics)

# Fit model GSL
rstan_options(auto_write = TRUE)
gslmodel <- stan_model("stan/modelGrowthGSL.stan")
fitgsl <- sampling(gslmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site", 
                                  "Ntreeid", "treeid", 
                                  "gsl"),
                warmup = 1000, iter = 2000, 
                chains = 4)
saveRDS(fitgsl, "output/stanOutput/fitGrowthGSL")

# Fit model SOS
rstan_options(auto_write = TRUE)
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsos <- sampling(sosmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site",
                                  "Ntreeid", "treeid",
                                  "sos"),
                warmup = 1000, iter = 2000,
                chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOS")

# Fit model EOS
rstan_options(auto_write = TRUE)
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteos <- sampling(eosmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site",
                                  "Ntreeid", "treeid",
                                  "eos"),
                warmup = 1000, iter = 2000,
                chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOS")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# jpeg("figures/pairs.jpg", width = 5000, height = 5000, 
     # units = "px", res = 300)

# fit@model_pars

# pairs(fit, pars = c("a", "b",
#                     "sigma_atreeid",
#                     "sigma_y"))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GDD fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/empiricalData/gddModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitgdd[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fitgdd[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgdd[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fitgdd[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgdd[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitgdd[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgdd[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgdd[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgdd[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GSL fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitgsl <- as.data.frame(fitgsl)

# full posterior
columns <- colnames(df_fitgsl)[!grepl("prior", colnames(df_fitgsl))]
sigma_df <- df_fitgsl[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgsl[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgsl[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fitgsl[, columns[grepl("aspp", columns)]]
site_df <- df_fitgsl[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgsl, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgsl, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgsl, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgsl, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fitgsl, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for GSL fit #####
pdf(file = "figures/empiricalData/gslModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitgsl[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fitgsl[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgsl[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fitgsl[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgsl[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitgsl[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgsl[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgsl[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgsl[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot SOS fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitsos <- as.data.frame(fitsos)

# full posterior
columns <- colnames(df_fitsos)[!grepl("prior", colnames(df_fitsos))]
sigma_df <- df_fitsos[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitsos[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitsos[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fitsos[, columns[grepl("aspp", columns)]]
site_df <- df_fitsos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for sos fit #####
pdf(file = "figures/empiricalData/sosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitsos[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fitsos[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitsos[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fitsos[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitsos[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fitsos[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitsos[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitsos[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitsos[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot EOS fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fiteos <- as.data.frame(fiteos)

# full posterior
columns <- colnames(df_fiteos)[!grepl("prior", colnames(df_fiteos))]
sigma_df <- df_fiteos[, columns[grepl("sigma", columns)]]
bspp_df <- df_fiteos[, columns[grepl("bsp", columns)]]
treeid_df <- df_fiteos[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fiteos[, columns[grepl("aspp", columns)]]
site_df <- df_fiteos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for eos fit #####
pdf(file = "figures/empiricalData/eosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fiteos[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fiteos[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fiteos[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fiteos[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fiteos[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fiteos[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fiteos[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fiteos[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fiteos[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (FALSE){
  
samples <- util$extract_expectand_vals(fitgsl)
jpeg(
  filename = "figures/modelGrowthGDD/retrodictiveCheckHist.jpeg",
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Diagnostics ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Parameterization for treeid
if (FALSE) { 
  samples <- util$extract_expectand_vals(fit)
  
  # atreeid
  atreeid <- names(samples)[grepl("zatreeid", names(samples))]
  atreeid <- atreeid[!grepl("sigma", atreeid)]
  atreeid <- atreeid[sample(length(unique(atreeid)), 9)]
  # pdf("figures/troubleShootingGrowthModel/atreeidParameterization_only_atreeid_no_b.pdf", 
  #     width = 6, height = 18)
  jpeg("figures/atreeidParameterization.jpeg", 
       width = 2000, height = 3000,
       units = "px", res = 300)
  util$plot_div_pairs(atreeid, "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
  dev.off()
}

# other diagnostics?
if (FALSE) {
  
  diagnostics <- util$extract_hmc_diagnostics(fit_noncentered) 
  util$check_all_hmc_diagnostics(diagnostics)
  
  samples <- util$extract_expectand_vals(fit_noncentered)
  
  util$plot_div_pairs("zbsp[1]", "sigma_bsp", samples, diagnostics, transforms = list("sigma_bsp" = 1))
  
  util$plot_div_pairs("zaspp[1]", "sigma_aspp", samples, diagnostics, transforms = list("sigma_aspp" = 1))
  
  util$plot_div_pairs("zasite[1]", "sigma_asite", samples, diagnostics, transforms = list("sigma_asite" = 1))
  
  util$plot_div_pairs("atreeid[1]", "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
  
}

}


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Full data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit model GDD
empgdd <- emp[!is.na(emp$pgsGDD5),]
nrow(emp[!is.na(emp$pgsGDD5),])
nrow(emp[!is.na(emp$pgsGSL),])
nrow(emp[!is.na(emp$leafout),])
nrow(emp[!is.na(emp$budset),])

# transform my groups to numeric values
empgdd$site_num <- match(empgdd$site, unique(empgdd$site))
empgdd$spp_num <- match(empgdd$spp, unique(empgdd$spp))
empgdd$treeid_num <- match(empgdd$treeid, unique(empgdd$treeid))

# transform data in vectors for GDD
y <- empgdd$lengthCM*10 # ring width in mm
N <- nrow(empgdd)
Nspp <- length(unique(empgdd$spp_num))
Nsite <- length(unique(empgdd$site_num))
site <- as.numeric(as.character(empgdd$site_num))
species <- as.numeric(as.character(empgdd$spp_num))
treeid <- as.numeric(empgdd$treeid_num)
Ntreeid <- length(unique(treeid))
gdd <- empgdd$pgsGDD5/200

# Fit model GDD
rstan_options(auto_write = TRUE)
gddmodel <- stan_model("stan/twolevelhierint.stan")
fitgdd <- sampling(gddmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gdd"),
                   warmup = 1000, iter=2000, 
                   chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDDFull")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fitgdd) 
util$check_all_hmc_diagnostics(diagnostics)

# Fit model GSL
empgsl <- emp[!is.na(emp$pgsGSL),]

# transform my groups to numeric values
empgsl$site_num <- match(empgsl$site, unique(empgsl$site))
empgsl$spp_num <- match(empgsl$spp, unique(empgsl$spp))
empgsl$treeid_num <- match(empgsl$treeid, unique(empgsl$treeid))

# transform data in vectors for gsl
y <- empgsl$lengthCM*10 # ring width in mm
N <- nrow(empgsl)
Nspp <- length(unique(empgsl$spp_num))
Nsite <- length(unique(empgsl$site_num))
site <- as.numeric(as.character(empgsl$site_num))
species <- as.numeric(as.character(empgsl$spp_num))
treeid <- as.numeric(empgsl$treeid_num)
Ntreeid <- length(unique(treeid))
gsl <- empgsl$pgsGSL/10

rstan_options(auto_write = TRUE)
gslmodel <- stan_model("stan/modelGrowthGSL.stan")
fitgsl <- sampling(gslmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gsl"),
                   warmup = 1000, iter = 2000, 
                   chains = 4)
saveRDS(fitgsl, "output/stanOutput/fitGrowthGSLFull")

# Fit model SOS
empsos <- emp[!is.na(emp$leafout),]
nrow(emp2) - nrow(empsos)
nrow(empgdd)
setdiff(emp2$treeid, emp$treeid)
nrow(subset(emp2, treeid %in% setdiff(emp2$treeid, emp$treeid)))
nrow(empgsl)
nrow(emp)
# transform my groups to numeric values
empgsl$site_num <- match(empgsl$site, unique(empgsl$site))
empgsl$spp_num <- match(empgsl$spp, unique(empgsl$spp))
empgsl$treeid_num <- match(empgsl$treeid, unique(empgsl$treeid))

# transform data in vectors for gsl
y <- empgsl$lengthCM*10 # ring width in mm
N <- nrow(empgsl)
Nspp <- length(unique(empgsl$spp_num))
Nsite <- length(unique(empgsl$site_num))
site <- as.numeric(as.character(empgsl$site_num))
species <- as.numeric(as.character(empgsl$spp_num))
treeid <- as.numeric(empgsl$treeid_num)
Ntreeid <- length(unique(treeid))
gsl <- empgsl$pgsGSL/10
rstan_options(auto_write = TRUE)
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsos <- sampling(sosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "sos"),
                   warmup = 1000, iter = 2000,
                   chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOSFull")

# Fit model EOS
rstan_options(auto_write = TRUE)
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteos <- sampling(eosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "eos"),
                   warmup = 1000, iter = 2000,
                   chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOSFull")