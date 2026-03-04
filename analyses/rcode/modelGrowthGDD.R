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

# checks
nrow(emp[!is.na(emp$pgsGDDAVG),]) - nrow(emp[!is.na(emp$pgsGDD5),])


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
gdd <- emp2$pgsGDD5/200
# gdd <- emp$pgsGDD10/200

rstan_options(auto_write = TRUE)
wcmodel <- stan_model("stan/twolevelhierint.stan")
fit <- sampling(wcmodel, data = c("N","y",
                                "Nspp","species",
                                "Nsite", "site", 
                                "Ntreeid", "treeid", 
                                "gdd"),
                iter=1000, chains=4, cores=4)

saveRDS(fit, "output/stanOutput/fitGrowthGDD")

# Fit model GSL
emp2 <- emp[!is.na(emp$pgsGSL),]

# transform my groups to numeric values
emp2$site_num <- match(emp2$site, unique(emp2$site))
emp2$spp_num <- match(emp2$spp, unique(emp2$spp))
emp2$treeid_num <- match(emp2$treeid, unique(emp2$treeid))

y <- emp2$lengthCM*10 # ring width in mm
N <- nrow(emp2)
Nspp <- length(unique(emp2$spp_num))
Nsite <- length(unique(emp2$site_num))
site <- as.numeric(as.character(emp2$site_num))
species <- as.numeric(as.character(emp2$spp_num))
treeid <- as.numeric(emp2$treeid_num)
Ntreeid <- length(unique(treeid))
gsl <- as.numeric(emp2$pgsGSL/10)

rstan_options(auto_write = TRUE)
gslmodel <- stan_model("stan/modelGrowthGSL.stan")
fitgsl <- sampling(gslmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site", 
                                  "Ntreeid", "treeid", 
                                  "gsl"),
                warmup = 1000, iter = 2000, 
                chains = 4)

saveRDS(fit, "output/stanOutput/fitGrowthGSL")

# Fit model SOS
sos <- emp2$pgsGDD5/200
rstan_options(auto_write = TRUE)
wcmodel <- stan_model("stan/modelGrowthSOS.stan")
fit <- sampling(wcmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site", 
                                  "Ntreeid", "treeid", 
                                  "sos"),
                warmup = 1000, iter = 2000, 
                chains=4)
saveRDS(fit, "output/stanOutput/fitGrowthSOS")

# Fit model EOS
eos <- emp2$pgsGDD5/200
rstan_options(auto_write = TRUE)
wcmodel <- stan_model("stan/modelGrowthEOS.stan")
fit <- sampling(wcmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site", 
                                  "Ntreeid", "treeid", 
                                  "eos"),
                warmup = 1000, iter = 2000, 
                chains=4)
saveRDS(fit, "output/stanOutput/fitGrowthEOS")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# jpeg("figures/pairs.jpg", width = 5000, height = 5000, 
     # units = "px", res = 300)

# fit@model_pars

# pairs(fit, pars = c("a", "b",
#                     "sigma_atreeid",
#                     "sigma_y"))

# === === === === === === === === === === === === #
##### Recover parameters from the posterior ##### 
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)

# full posterior
columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
treeid_df <- df_fit[, grepl("treeid", columns) & 
                      !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
site_df <- df_fit[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(asite_df) <- 1:ncol(asite_df)

# posterior summaries
sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fit, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fit, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fit, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# create prior vectors
sigmaatreeidprior <- df_fit[, "sigma_atreeid_prior"]
sigmayprior <- df_fit[, "sigma_y_prior"]
bsppprior <- df_fit[, "bsp_prior"]
asppprior <- df_fit[, "aspp_prior"]
asiteprior <- df_fit[, "asite_prior"]

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot posterior vs priors for gdd fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
pdf(file = "figures/empiricalData/gddModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

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

dev.off()




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


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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
