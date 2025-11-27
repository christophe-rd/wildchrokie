
# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)
# quartz()

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(shinystan)
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

# === === === === === === === === === === === === === === === === 
#### Run model on empirical data ####
# === === === === === === === === === === === === === === === === 
emp <- read.csv("output/empiricalDataMAIN.csv")

# transform my groups to numeric values
emp$spp_num <- match(emp$spp, unique(emp$spp))

# transform data in vectors
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
Nspp <- length(unique(emp$spp_num))
species <- as.numeric(as.character(emp$spp_num))

# check that everything is fine
table(treeid,species)

rstan_options(auto_write = TRUE)
 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fit <- stan("stan/twolevelhierint_only_asp.stan", 
            data=c("N","y","Nspp","species", "gdd"),
            iter=4000, chains=4, cores=4)

# saveRDS(fit_only_asp, "output/stanOutput/fit_only_asp")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# jpeg("figures/pairs.jpg", width = 5000, height = 5000, 
     # units = "px", res = 300)

# fit@model_pars

# pairs(fit, pars = c("a", "b",
#                     "sigma_asp",
#                     "sigma_y"))

# === === === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a spp  ######
aspp_cols <- colnames(df_fit)[grepl("asp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("zasp", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("sigma", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("asp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_a_spp = numeric(ncol(aspp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per25 = NA,
  fit_a_spp_per75 = NA,
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_a_spp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_a_spp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_a_spp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_a_spp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_a_spp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Diagnostics ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Parameterization for asp

diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)

asp <- names(samples)[grepl("asp", names(samples))]
asp <- asp[!grepl("zasp", asp)]
asp <- asp[!grepl("sigma", asp)]

zasp <- names(samples)[grepl("zasp", names(samples))]

jpeg("figures/aspParameterization_only_asp.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
util$plot_div_pairs(asp, "sigma_asp", samples, diagnostics, transforms = list("sigma_asp" = 1))
dev.off()

jpeg("figures/zaspParameterization_only_asp.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
util$plot_div_pairs(zasp, "sigma_asp", samples, diagnostics, transforms = list("sigma_asp" = 1))
dev.off()

# === === === === === === === #
# Plot parameter recovery #####
# === === === === === === === #

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot sigmas ######
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
sigma_long$prior[which(sigma_long$parameter == "sigma_bsp")] <- rnorm(8e3, 0, 0.2)
sigma_long$prior[which(sigma_long$parameter == "sigma_asp")] <- rnorm(8e3, 0, 0.3)
sigma_long$prior[which(sigma_long$parameter == "sigma_asite")] <- rnorm(8e3, 0, 0.3)
sigma_long$prior[which(sigma_long$parameter == "sigma_atreeid")] <- rnorm(8e3, 0, 0.1)
sigma_long$prior[which(sigma_long$parameter == "sigma_y")] <- rnorm(8e3, 0, 0.1)

ggplot(sigma_long) +
  geom_density(aes(x = prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(aes(x = value, colour = "Posterior"),
               linewidth = 0.8) +
  facet_wrap(~parameter) + 
  labs(title = "priorVSposterior_sigmas",
       x = "", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot b spp ######
bspp_long <- reshape(
  bspp_df,
  direction = "long",
  varying = list(names(bspp_df)),
  v.names = "value",
  timevar = "spp",
  times = names(bspp_df),
  idvar = "draw"
)
bspp_long

# simulate priors
hyperparameter_draws <- 8000
parameter_draws <- 1000
n_sigma_bsp <- 200

# set to prior values
sigma_bsp_vec <- abs(rnorm(n_sigma_bsp, 0, 0.2))

prior_bsp <- rep(NA, parameter_draws*length(sigma_bsp_vec))

for (i in 1: length(sigma_bsp_vec)) {
  prior_bsp[((i - 1)*parameter_draws + 1):(i*parameter_draws)] <- rnorm(parameter_draws, 0, sigma_bsp_vec[i])
}
prior_bsp

ggplot() +
  geom_density(data = data.frame(prior_bsp = prior_bsp),
               aes(x = prior_bsp, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = bspp_long,
               aes(x = value, colour = "Posterior"),
               linewidth = 0.8) +
  facet_wrap(~spp) + 
  labs(title = "priorVSposterior_bsp",
       x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
bspp_long <- reshape(
  bspp_df,
  direction = "long",
  varying = list(names(bspp_df)),
  v.names = "value",
  timevar = "spp",
  times = names(bspp_df),
  idvar = "draw"
)
bspp_long

ggplot() +
  # geom_density(data = data.frame(sigma_bsp_draw = sigma_bsp_draw),
  #              aes(x = sigma_bsp_draw, colour = "Prior"),
  #              linewidth = 0.8) +
  geom_density(data = bspp_long,
               aes(x = value, colour = "Posterior"),
               linewidth = 0.8) +
  facet_wrap(~spp) + 
  labs(title = "priorVSposterior_bsp",
       x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot a spp ######

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######

