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
                iter=4000, chains=4, cores=4)

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

sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fit, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fit, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z", treeid))
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fit, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

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

###### Plot a prior vs posterior ######
a_posterior <- df_fit[, colnames(df_fit) %in% "a"]

a_prior <- rnorm(1e4, 2, 3)

priora <- ggplot() +
  geom_density(data = data.frame(a = a_prior),
               aes(x = a, colour = "Prior at N(5,3)"),
               linewidth = 1) +
  geom_density(data = data.frame(value = a_posterior),
               aes(x = value, colour = "Posterior"),
               linewidth = 1) +
  labs(title = "priorVSposterior_a",
       x = "a", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
priora

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot sigmas prior vs posterior ######
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
sigma_long <- reshape(
  sigma_df,
  direction = "long",
  varying = list(names(sigma_df)),
  v.names = "value",
  timevar = "parameter",
  times = names(sigma_df),
  idvar = "draw"
)
sigma_long

sigma_long$prior <- NA
sigma_long$prior[which(sigma_long$parameter == "sigma_atreeid")] <- rnorm(8e3, 0, 0.5)
sigma_long$prior[which(sigma_long$parameter == "sigma_y")] <- rnorm(8e3, 0, 3)

priorsigmas <- ggplot(sigma_long) +
  geom_density(aes(x = prior, colour = "Prior sigma_atreeid  at N(0, 0.5)
Prior sigma_y at N(0,3)"),
               linewidth = 0.8) +
  geom_density(aes(x = value, colour = "Posterior"),
               linewidth = 0.8) +
  facet_wrap(~parameter) + 
  labs(title = "priorVSposterior_sigmas",
       x = "", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
priorsigmas

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot atreeid prior vs posterior ######
treeid_long <- reshape(
  treeid_df,
  direction = "long",
  varying = list(names(treeid_df)),
  v.names = "value",
  timevar = "treeid",
  times = names(treeid_df),
  idvar = "draw"
)
treeid_long

# simulate priors
hyperparameter_draws <- 8000
parameter_draws <- 1000
n_sigmatreeid <- 200

# set to prior values
sigmatreeid_vec <- abs(rnorm(n_sigmatreeid, 0, 0.5))

prior_treeid <- rep(NA, parameter_draws*length(sigmatreeid_vec))

for (i in 1: length(sigmatreeid_vec)) {
  prior_treeid[((i - 1)*parameter_draws + 1):(i*parameter_draws)] <- rnorm(parameter_draws, 0, sigmatreeid_vec[i])
}
prior_treeid

# sub of some treeids for plotting
subtreeid <- subset(treeid_long, treeid %in% sample(treeid_long$treeid, 5))

prioratreeid <- ggplot() +
  geom_density(data = data.frame(prior_treeid = prior_treeid),
               aes(x = prior_treeid, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = subtreeid,
               aes(x = value, colour = "Posterior"),
               linewidth = 0.8) +
  facet_wrap(~treeid) + 
  labs(title = "priorVSposterior_treeid",
       x = "treeid", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
prioratreeid

###### Plot aspp prior vs posterior ######
# convert posterior distribution to long format
aspp_long <- reshape(
  aspp_df,
  direction = "long",
  varying = list(names(aspp_df)),
  v.names = "value",
  timevar = "spp",
  times = names(aspp_df),
  idvar = "draw"
)
aspp_long

# aspp prior
aspp_prior <- rnorm(1e4, 0, 6)

prioraspp <- ggplot() +
  geom_density(data = data.frame(aspp_prior = aspp_prior),
               aes(x = aspp_prior, colour = "Prior at N(0, 6)"),
               linewidth = 0.8) +
  geom_density(data = aspp_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_aspp",
       x = "aspp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  xlim(c(-20, 20)) +
  theme_minimal()
prioraspp

###### Plot asite prior vs posterior ######
# convert posterior distribution to long format
asite_long <- reshape(
  site_df,
  direction = "long",
  varying = list(names(site_df)),
  v.names = "value",
  timevar = "site",
  times = names(site_df),
  idvar = "draw"
)
asite_long

asite_prior <- rnorm(1e4, 0, 2)

priorasite <- ggplot() +
  geom_density(data = data.frame(asite_prior = asite_prior),
               aes(x = asite_prior, colour = "Prior at N(0, 2)"),
               linewidth = 0.8) +
  geom_density(data = asite_long,
               aes(x = value, colour = "Posterior", group = site),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_asite",
       x = "asite", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  xlim(c(-6, 6)) +
  theme_minimal()
priorasite

###### Plot bsp prior vs posterior ######
# convert posterior distribution to long format
bsp_long <- reshape(
  bspp_df,
  direction = "long",
  varying = list(names(bspp_df)),
  v.names = "value",
  timevar = "spp",
  times = names(bspp_df),
  idvar = "draw"
)
bsp_long

# aspp prior
bsp_prior <- rnorm(1e4, 0, 0.5)

priorbsp <- ggplot() +
  geom_density(data = data.frame(bsp_prior = bsp_prior),
               aes(x = bsp_prior, colour = "Prior at N(0, 0.3)"),
               linewidth = 0.8) +
  geom_density(data = bsp_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_bsp",
       x = "bsp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
priorbsp

#  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
priorcombined <- (priorsigmas) / (priora) / (prioraspp) / (priorasite) / (priorbsp)
ggsave("figures/priorVSposteriorCombined.jpeg", priorcombined, width = 8, height = 12, units = "in", dpi = 300)
#  --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---



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
