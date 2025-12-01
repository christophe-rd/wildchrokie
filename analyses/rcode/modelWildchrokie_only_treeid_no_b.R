
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
nrow(emp)
# transform my groups to numeric values
alnbetpop <- subset(emp, spp %in% c("ALNINC", "BETPOP"))

emp$treeid_num <- match(emp$treeid, unique(emp$treeid))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$site_num <- match(emp$site, unique(emp$site))

# transform data in vectors
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))
species <- as.numeric(as.character(emp$spp_num))
Nspp <- length(unique(species))
site <- as.numeric(as.character(emp$site_num))
Nsite<- length(unique(site))

# check that everything is fine
table(treeid, species)
table(site)
table(species)

rstan_options(auto_write = TRUE)
 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fit <- stan("stan/twolevelhierint_only_atreeid_no_b.stan", 
            data=c("N","y", 
                   "Ntreeid", "treeid",
                   "Nspp","species",
                   "Nsite","site"),
            iter=4000, chains=4, cores=4)


# saveRDS(fit, "output/stanOutput/fit_a_atreeid_ppONasp_all_spp")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# jpeg("figures/pairs.jpg", width = 5000, height = 5000, 
     # units = "px", res = 300)

# fit@model_pars

# pairs(fit, pars = c("a", "b",
#                     "sigma_asp",
#                     "sigma_y"))

# Diagnostics ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Parameterization for asp

diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectand_vals(fit)



# asp
asp <- names(samples)[grepl("asp", names(samples))]
asp <- asp[!grepl("sigma", asp)]

jpeg("figures/troubleShootingGrowthModel/aspParameterization.jpg", 
     width = 2000, height = 2000,
     units = "px", res = 300)
# util$plot_div_pairs(asp, "sigma_asp", samples, diagnostics, transforms = list("sigma_asp" = 1))
dev.off()

# # asite
# asite <- names(samples)[grepl("asite", names(samples))]
# asite <- asite[!grepl("sigma", asite)]
# 
# jpeg("figures/asiteParameterization_only_asp_asite_atreeid.jpg", width = 2000, height = 2000, 
#      units = "px", res = 300)
# util$plot_div_pairs(asite, "sigma_asite", samples, diagnostics, transforms = list("sigma_asite" = 1))
# dev.off()

# atreeid
atreeid <- names(samples)[grepl("zatreeid", names(samples))]
atreeid <- atreeid[!grepl("sigma", atreeid)]
atreeid <- atreeid[sample(length(unique(atreeid)), 9)]
# pdf("figures/troubleShootingGrowthModel/atreeidParameterization_only_atreeid_no_b.pdf", 
#     width = 6, height = 18)
jpeg("figures/troubleShootingGrowthModel/atreeidParameterization_only_atreeid_no_b.jpeg", 
     width = 2000, height = 3000,
     units = "px", res = 300)
util$plot_div_pairs(atreeid, "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
dev.off()


###### Plot treeid ######
df_fit <- as.data.frame(fit)

# Prior checks ####
# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]

# remove sigma_asp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_a_treeid = numeric(ncol(treeid_df)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_a_treeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_a_treeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# plot
atreeid_long <- reshape(
  treeid_df,
  direction = "long",
  varying = list(names(treeid_df)),
  v.names = "value",
  timevar = "spp",
  times = names(treeid_df),
  idvar = "draw"
)
atreeid_long

# simulate priors
hyperparameter_draws <- 1000
parameter_draws <- 1000
n_sigma_atreeid <- 200

# set to prior values
sigma_atreeid_vec <- abs(rnorm(n_sigma_atreeid, 0, 1))

prior_atreeid <- rep(NA, parameter_draws*length(sigma_atreeid_vec))

for (i in 1: length(sigma_atreeid_vec)) {
  prior_atreeid[((i - 1)*parameter_draws + 1):(i*parameter_draws)] <- rnorm(parameter_draws, 0, sigma_atreeid_vec[i])
}
prior_atreeid

ggplot() +
  geom_density(data = data.frame(prior_atreeid = prior_atreeid),
               aes(x = prior_atreeid, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = atreeid_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.1) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_atreeid",
       x = "atreeid", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

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

# asp prior
asp_prior <- rnorm(1e4, 0, 0.5)

ggplot() +
  geom_density(data = data.frame(asp_prior = asp_prior),
               aes(x = asp_prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = aspp_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_atreeid",
       x = "asp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

###### Recover a site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per25 = NA,
  fit_a_site_per75 = NA,
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2

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

# asite prior
asite_prior <- rnorm(1e4, 0, 0.5)

ggplot() +
  geom_density(data = data.frame(asite_prior = asite_prior),
               aes(x = asite_prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = asite_long,
               aes(x = value, colour = "Posterior", group = site),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_atreeid",
       x = "asite", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()



# a #####
df_fit <- as.data.frame(fit)

# Prior checks ####
# grab treeid 
a_posterior <- df_fit[, colnames(df_fit) %in% "a"]

a_prior <- rnorm(1e4, 5, 1)

ggplot() +
  geom_density(data = data.frame(a = a_prior),
               aes(x = a, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = data.frame(value = a_posterior),
               aes(x = value, colour = "Posterior"),
               linewidth = 0.2) +
  labs(title = "priorVSposterior_a",
       x = "a", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

