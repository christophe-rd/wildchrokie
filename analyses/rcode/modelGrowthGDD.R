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
# nrow(emp[!is.na(emp$pgsGDDAVG),])
emp <- emp[!is.na(emp$pgsGDD5),]

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors for GDD
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
x <- emp$pgsGDD5/200
# gdd <- emp$pgsGDD10/200
Nspp <- length(unique(emp$spp_num))
Nsite <- length(unique(emp$site_num))
site <- as.numeric(as.character(emp$site_num))
species <- as.numeric(as.character(emp$spp_num))
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))

# check that everything is fine
table(treeid,species)

rstan_options(auto_write = TRUE)
 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
wcmodel <- stan_model("stan/twolevelhierint.stan")
fit <- sampling(wcmodel, data = c("N","y",
                                "Nspp","species",
                                "Nsite", "site", 
                                "Ntreeid", "treeid", 
                                "x"),
            warmup = 1000, iter = 2000, 
            chains=4)

saveRDS(fit, "output/stanOutput/fitGrowthGDD")

if (fitlmer) {

# fit stanlmer to check differences
emp$gdd <- emp$pgsGDD/200
emp$y <- emp$lengthCM*10

emp$site_fac <- as.factor(emp$site_num)
emp$spp_fac <- as.factor(emp$spp_num)
emp$treeid_fac <- as.factor(emp$treeid_num)

unique(emp$spp_fac)
emp$gdd_c <- scale(emp$gdd, scale = FALSE)

fitlmer_partialpooling1 <- stan_lmer(
  y ~
    (1|site_fac) +
    (gdd | spp_fac:treeid_fac),
  data = emp,
  chains = 4,
  iter = 4000,
  cores = 4
)

fitlmer_partialpooling3 <- stan_lmer(
  y ~
    (1|site_fac) +
    (gdd | spp_fac),
  data = emp,
  chains = 4,
  iter = 4000,
  cores = 4
)

colnames(as.data.frame(fitlmer_partialpooling))[1:20]
colnames(as.data.frame(fitlmer_partialpooling))[80:300]
fitlmer$stanfit
pairs(fit, pars = c("a", "b",
                    "sigma_atreeid",
                    "sigma_y"))

}

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
# === === === === === === === #
##### Plot parameter recovery #####
# === === === === === === === #

###### Plot a prior vs posterior ######
a_posterior <- df_fit[, colnames(df_fit) %in% "a"]

a_prior <- rnorm(1e4, 5, 3)

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


# RECOVER FROM STAN_LMER ####
df_fit_lmer <- as.data.frame(fitlmer_partialpooling1)

coef(fitlmer_partialpooling1)$`spp_fac:treeid_fac`
ranef(fitlmer_partialpooling1)$`spp_fac:treeid_fac`

VarCorr(fitlmer_partialpooling1)$`spp_fac:treeid_fac`
print(colnames(as.matrix(fitlmer_partialpooling1)))
vcov(fitlmer_partialpooling1, correlation = TRUE)
vcov(fitlmer_partialpooling1)

cov_mat <- VarCorr(fitlmer_partialpooling1)$`spp_fac:treeid_fac`

fitlmer_partialpooling1$stan_function
fitlmer_partialpooling1$stanfit
decov(fitlmer_partialpooling1)
# aspp
fixef(fitlmer_partialpooling1)
ranef(fitlmer_partialpooling1)
mean(df_fit[,"a"])

treeidtest <- ranef(fitlmer_partialpooling1)["spp_fac:treeid_fac"]
treeid_lmer <- as.data.frame(treeidtest)

colnames(treeid_lmer) <- c("aspp_lmer", "bspp_lmer")
treeid_lmer$spp <- substr(rownames(treeid_lmer), 1,1)
treeid_lmer$treeid <- substr(rownames(treeid_lmer), 3,4)

# calculate average per species
aspp_aver <- aggregate(aspp_lmer ~ spp, treeid_lmer, FUN = mean)
aspp_quan <- aggregate(aspp_lmer ~ spp, treeid_lmer, FUN = quantile)
aspp_lmer <- data.frame(
  spp = aspp_aver$spp,
  mean_lmer = aspp_aver$aspp_lmer,
  per25_lmer = aspp_quan$aspp_lmer[, 1],
  per75_lmer = aspp_quan$aspp_lmer[, 2]
)

asppcomp <- merge(aspp_df2, aspp_lmer, by = "spp")

asppcompplot <- ggplot(asppcomp, aes(x = mean_lmer, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 0.5) +
  geom_point(color = "#046C9A", size = 3) +
  labs(x = "stanlmer", y = "rstan",title = "aspp") +
  theme_minimal()
asppcompplot

# bspp
bspp_aver <- aggregate(bspp_lmer ~ spp, treeid_lmer, FUN = mean)
bspp_quan <- aggregate(bspp_lmer ~ spp, treeid_lmer, FUN = quantile)
bspp_lmer <- data.frame(
  spp = bspp_aver$spp,
  mean_lmer = bspp_aver$bspp_lmer,
  per25_lmer = bspp_quan$bspp_lmer[, 1],
  per75_lmer = bspp_quan$bspp_lmer[, 2]
)

bsppcomp <- merge(bspp_df2, bspp_lmer, by = "spp")

bsppcompplot <- ggplot(bsppcomp, aes(x = mean_lmer, y = fit_bspp)) +
  geom_errorbar(aes(ymin = fit_bspp_per25, ymax = fit_bspp_per75),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 0.5) +
  geom_point(color = "#046C9A", size = 3) +
  labs(x = "stanlmer", y = "rstan", title = "bspp") +
  theme_minimal()
bsppcompplot

# asite
asitetest <- ranef(fitlmer_partialpooling1)["site_fac"]
site_lmer <- as.data.frame(asitetest)
colnames(site_lmer) <- c("asite_lmer")

sitebind <- cbind(site_df2, site_lmer)

asitelmer <- ggplot(sitebind, aes(x = asite_lmer, y = fit_asite)) +
  geom_errorbar(aes(ymin = fit_asite_per25, ymax = fit_asite_per75),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 0.5) +
  geom_point(color = "#046C9A", size = 3) +
  labs(x = "stanlmer", y = "rstan", title = "asite") +
  theme_minimal()
asitelmer

combined_plot <- (asppcompplot + bsppcompplot + asitelmer)
combined_plot
ggsave("figures/troubleShootingGrowthModel/combinedPlots_lmerVSStanM1.jpeg", combined_plot, width = 10, height = 8, units = "in", dpi = 300)


fitlmer_partialpooling$coefficients["spp_fac"]

sigma_cols <- colnames(df_fit)[
  grepl("igma", colnames(df_fit))
]

sigma_df_lmer <- df_fit[, colnames(df_fit) %in% sigma_cols][,1:2]

colnames(sigma_df_lmer) <- c("sigma_y", "sigma_atreeid")

# SIGMAS
sigma_df2_lmer <- data.frame(
  sigma = character(ncol(sigma_df_lmer)),
  mean = numeric(ncol(sigma_df_lmer)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)

for (i in 1:ncol(sigma_df_lmer)) {
  sigma_df2_lmer$sigma[i] <- colnames(sigma_df_lmer)[i]
  sigma_df2_lmer$mean[i] <- round(mean(sigma_df_lmer[[i]]),3)
  sigma_df2_lmer$per5[i] <- round(quantile(sigma_df_lmer[[i]], 0.05),3)
  sigma_df2_lmer$per25[i] <- round(quantile(sigma_df_lmer[[i]], 0.25),3)
  sigma_df2_lmer$per75[i] <- round(quantile(sigma_df_lmer[[i]], 0.75),3)
  sigma_df2_lmer$per95[i] <- round(quantile(sigma_df_lmer[[i]], 0.95),3)
}
sigma_df2_lmer

# BSPP
bspp_cols <- colnames(df_fit)[
  grepl("gdd:", colnames(df_fit))
]

bspp_df_lmer <- df_fit[, colnames(df_fit) %in% bspp_cols]

colnames(bspp_df_lmer) <- 1:4

bspp_df2_lmer <- data.frame(
  spp = character(ncol(bspp_df_lmer)),
  fit_bspp = numeric(ncol(bspp_df_lmer)),  
  fit_bspp_per5 = NA, 
  fit_bspp_per25 = NA,
  fit_bspp_per75 = NA,
  fit_bspp_per95 = NA
)

for (i in 1:ncol(bspp_df_lmer)) {
  bspp_df2_lmer$spp[i] <- colnames(bspp_df_lmer)[i]
  bspp_df2_lmer$fit_bspp[i] <- round(mean(bspp_df_lmer[[i]]),3)
  bspp_df2_lmer$fit_bspp_per5[i] <- round(quantile(bspp_df_lmer[[i]], 0.05),3)
  bspp_df2_lmer$fit_bspp_per25[i] <- round(quantile(bspp_df_lmer[[i]], 0.25),3)
  bspp_df2_lmer$fit_bspp_per75[i] <- round(quantile(bspp_df_lmer[[i]], 0.75),3)
  bspp_df2_lmer$fit_bspp_per95[i] <- round(quantile(bspp_df_lmer[[i]], 0.95),3)
}
bspp_df2_lmer

# ASPP
aspp_cols <- colnames(df_fit)[
  !grepl(":gdd", colnames(df_fit))
]
aspp_cols <- aspp_cols[
  grepl("spp_", aspp_cols) & !grepl("Intercept", aspp_cols)
]

aspp_df_lmer <- df_fit[, colnames(df_fit) %in% aspp_cols]
colnames(aspp_df_lmer) <- 1:4

aspp_df2_lmer <- data.frame(
  spp = character(ncol(aspp_df_lmer)),
  fit_aspp = numeric(ncol(aspp_df_lmer)),  
  fit_aspp_per5 = NA, 
  fit_aspp_per25 = NA,
  fit_aspp_per75 = NA,
  fit_aspp_per95 = NA
)

for (i in 1:ncol(aspp_df_lmer)) {
  aspp_df2_lmer$spp[i] <- colnames(aspp_df_lmer)[i]
  aspp_df2_lmer$fit_aspp[i] <- round(mean(aspp_df_lmer[[i]]),3)
  aspp_df2_lmer$fit_aspp_per5[i] <- round(quantile(aspp_df_lmer[[i]], 0.05),3)
  aspp_df2_lmer$fit_aspp_per25[i] <- round(quantile(aspp_df_lmer[[i]], 0.25),3)
  aspp_df2_lmer$fit_aspp_per75[i] <- round(quantile(aspp_df_lmer[[i]], 0.75),3)
  aspp_df2_lmer$fit_aspp_per95[i] <- round(quantile(aspp_df_lmer[[i]], 0.95),3)
}
aspp_df2_lmer

# SITE
site_cols <- colnames(df_fit)[grepl("site_fac", colnames(df_fit))]
site_cols <- site_cols[!grepl("Sigma", site_cols)]
site_df_lmer <- df_fit[, colnames(df_fit) %in% site_cols]

colnames(site_df_lmer) <- 1:4

site_df2_lmer <- data.frame(
  site = character(ncol(site_df_lmer)),
  fit_asite = numeric(ncol(site_df_lmer)),  
  fit_asite_per5 = NA, 
  fit_asite_per25 = NA,
  fit_asite_per75 = NA,
  fit_asite_per95 = NA
)

for (i in 1:ncol(site_df_lmer)) {
  site_df2_lmer$site[i] <- colnames(site_df_lmer)[i]
  site_df2_lmer$fit_asite[i] <- round(mean(site_df_lmer[[i]]),3)
  site_df2_lmer$fit_asite_per5[i] <- round(quantile(site_df_lmer[[i]], 0.05),3)
  site_df2_lmer$fit_asite_per25[i] <- round(quantile(site_df_lmer[[i]], 0.25),3)
  site_df2_lmer$fit_asite_per75[i] <- round(quantile(site_df_lmer[[i]], 0.75),3)
  site_df2_lmer$fit_asite_per95[i] <- round(quantile(site_df_lmer[[i]], 0.95),3)
}
site_df2_lmer

# TREEID
treeid_cols <- colnames(df_fit)[grepl("treeid", colnames(df_fit))]
treeid_cols <- treeid_cols[!grepl("Sigma", treeid_cols)]
treeid_df_lmer <- df_fit[, colnames(df_fit) %in% treeid_cols]

colnames(treeid_df_lmer) <- unique(emp$treeid_num)

treeid_df2_lmer <- data.frame(
  treeid = character(ncol(treeid_df_lmer)),
  fit_atreeid = numeric(ncol(treeid_df_lmer)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)

for (i in 1:ncol(treeid_df_lmer)) {
  treeid_df2_lmer$treeid[i] <- colnames(treeid_df_lmer)[i]
  treeid_df2_lmer$fit_atreeid[i] <- round(mean(treeid_df_lmer[[i]]),3)
  treeid_df2_lmer$fit_atreeid_per5[i] <- round(quantile(treeid_df_lmer[[i]], 0.05),3)
  treeid_df2_lmer$fit_atreeid_per25[i] <- round(quantile(treeid_df_lmer[[i]], 0.25),3)
  treeid_df2_lmer$fit_atreeid_per75[i] <- round(quantile(treeid_df_lmer[[i]], 0.75),3)
  treeid_df2_lmer$fit_atreeid_per95[i] <- round(quantile(treeid_df_lmer[[i]], 0.95),3)
}
treeid_df2_lmer

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Plot two model comparisons #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot sigmas ######
colnames(sigma_df2_lmer)[2:ncol(sigma_df2_lmer)] <- paste(colnames(sigma_df2_lmer)[2:ncol(sigma_df2_lmer)], "lmer", sep = "_")

sigmaforplot <- merge(sigma_df2, sigma_df2_lmer, by = "sigma")

sigma_stanXstanlmer <- ggplot(sigmaforplot, aes(x = mean, y = mean_lmer)) +
  geom_errorbar(aes(xmin = per25, xmax = per75),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_errorbar(aes(ymin = per25_lmer, ymax = per75_lmer),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 0.5) +
  geom_point(color = "#046C9A", size = 3) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "") +
  theme_minimal()
sigma_stanXstanlmer
ggsave("figures/sigma_stanXstanlmer.jpeg", sigma_stanXstanlmer, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot b spp ######
colnames(bspp_df2_lmer)[2:ncol(bspp_df2_lmer)] <- paste(colnames(bspp_df2_lmer)[2:ncol(bspp_df2_lmer)], "lmer", sep = "_")

bsppforplot <- merge(bspp_df2, bspp_df2_lmer, by = "spp")

bspp_stanXstanlmer <- ggplot(bsppforplot, aes(x = fit_bspp, y = fit_bspp_lmer)) +
  geom_errorbar(aes(xmin = fit_bspp_per25, xmax = fit_bspp_per75), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=1) +
  geom_errorbar(aes(ymin = fit_bspp_per25_lmer, ymax = fit_bspp_per75_lmer), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_point(color = "#046C9A", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 0.5) +
  labs(x = "sim bspp", y = "fit bspp", title = "") +
  theme_minimal()
bspp_stanXstanlmer
# ggsave!
ggsave("figures/bspp_stanXstanlmer.jpeg", bspp_stanXstanlmer, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
colnames(treeid_df2_lmer)[2:ncol(treeid_df2_lmer)] <- paste(colnames(treeid_df2_lmer)[2:ncol(treeid_df2_lmer)], "lmer", sep = "_")

treeidforplot <- merge(treeid_df2, treeid_df2_lmer, by = "treeid")

# plot treeid
atreeid_stanXstanlmer <- ggplot(treeidforplot, aes(x = fit_atreeid, y = fit_atreeid_lmer)) +
  geom_errorbar(aes(xmin = fit_atreeid_per25, xmax = fit_atreeid_per75), 
                width = 0, linewidth = 0.7, color = "darkgray", alpha= 0.7) +
  geom_errorbar(aes(ymin = fit_atreeid_per25_lmer, ymax = fit_atreeid_per75_lmer), 
                width = 0, linewidth = 0.7, color = "darkgray", alpha = 0.7) +
  geom_point(color = "#046C9A", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 0.5) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "") +
  theme_minimal()
atreeid_stanXstanlmer
# ggsave!
ggsave("figures/atreeid_stanXstanlmer.jpeg", atreeid_stanXstanlmer, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot a spp ######
colnames(aspp_df2_lmer)[2:ncol(aspp_df2_lmer)] <- paste(colnames(aspp_df2_lmer)[2:ncol(aspp_df2_lmer)], "lmer", sep = "_")

asppforplot <- merge(aspp_df2, aspp_df2_lmer, by = "spp")

aspp_stanXstanlmer <- ggplot(asppforplot, aes(x = fit_aspp, y = fit_aspp_lmer)) +
  geom_errorbar(aes(xmin = fit_aspp_per5, xmax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25_lmer, ymax = fit_aspp_per75_lmer), 
                width = 0, linewidth = 0.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 0.5) +
  labs(x = "sim aspp", y = "fit aspp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
aspp_stanXstanlmer
# ggsave!
ggsave("figures/aspp_stanXstanlmer.jpeg", aspp_stanXstanlmer, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######
colnames(site_df2_lmer)[2:ncol(site_df2_lmer)] <- paste(colnames(site_df2_lmer)[2:ncol(site_df2_lmer)], "lmer", sep = "_")

asiteforplot <- merge(site_df2, site_df2_lmer, by = "site")

asite_stanXstanlmer <- ggplot(asiteforplot, aes(x = fit_asite, y = fit_asite_lmer)) +
  geom_errorbar(aes(xmin = fit_asite_per5, xmax = fit_asite_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_asite_per25_lmer, ymax = fit_asite_per75_lmer), 
                width = 0, linewidth = 0.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 0.5) +
  labs(x = "sim asite", y = "fit asite", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
asite_stanXstanlmer
# ggsave!
ggsave("figures/asite_stanXstanlmer.jpeg", asite_stanXstanlmer, width = 6, height = 6, units = "in", dpi = 300)

###### Combine plots  ######
combined_plot <- (atreeid_stanXstanlmer) /
  (sigma_stanXstanlmer + bspp_stanXstanlmer ) /
  (aspp_stanXstanlmer + asite_stanXstanlmer)
combined_plot
ggsave("figures/troubleShootingGrowthModel/combinedPlots.jpeg", combined_plot, width = 10, height = 8, units = "in", dpi = 300)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
samples <- util$extract_expectand_vals(fit)
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
