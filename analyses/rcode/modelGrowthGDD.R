
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

runSimData <- FALSE
if (runSimData) {
# === === === === === === === === === === === === === === === === 
#### SIMULATED DATA ####
# === === === === === === === === === === === === === === === ===

# set parameters
set.seed(124)
a <- 15
b <- 4
sigma_y <- 0.1
sigma_aspp <- 0.5 # This is pretty low, but I guess you think your species are closely related and will be similar?
sigma_atreeid <- 0.15
sigma_asite <- 0.3
sigma_bspp <- 0.25

n_site <- 10 # number of sites
n_spp <- 10 # number of species
n_perspp <- 10 # number of individuals per species
n_treeid <- n_perspp * n_spp * n_site # number of treeid
n_meas <- 5 # repeated measurements per id
N <- n_treeid * n_meas # total number of measurements
N

# get replicated treeid
treeid <- rep(1:n_treeid, each = n_meas)
# non replicated treeid
treeidnonrep <- rep(rep(1:n_perspp, times = n_spp), times = n_site)
# replicated spp
spp <- rep(rep(rep(1:n_spp, each = n_perspp), times = n_site), each = n_meas) 
# non replicated spp
spp_nonrep <- rep(rep(1:n_spp, each = n_perspp), each = n_site) 
# replicated site
site <- rep(rep(rep(1:n_site, each = n_spp), each = n_perspp), each = n_meas)
# non replicated site
site_nonrep <- rep(rep(1:n_site, each = n_spp), each = n_perspp)
# quick check 
table(treeidnonrep, site_nonrep)

simcoef <- data.frame(
  site = site,
  spp = spp,
  treeid = treeid
)

# get intercept values for each species
aspp <- rnorm(n_spp, 0, sigma_aspp)
asite <- rnorm(n_site, 0, sigma_asite)
atreeid <- rnorm(n_treeid, 0, sigma_atreeid)

# get slope values for each speciess
bspp <- rnorm(n_spp, 0, sigma_bspp)

# Add my parameters to the df
simcoef$atreeid <- atreeid[treeid]
simcoef$asite <- asite[simcoef$site]
simcoef$aspp <- aspp[simcoef$spp]
simcoef$bspp <- bspp[simcoef$spp]

# add the rest of the boring stuff 
simcoef$a <- a
simcoef$b <- b
simcoef$sigma_y <- sigma_y
simcoef$sigma_atreeid <- sigma_atreeid
simcoef$sigma_aspp <- sigma_aspp
simcoef$sigma_asite <- sigma_asite
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gdd <- rnorm(N, 1800, 100)
simcoef$gddcons <- simcoef$gdd/200

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$asite + 
  simcoef$aspp + 
  simcoef$atreeid + 
  simcoef$a +
  (simcoef$b*simcoef$gddcons) + 
  (simcoef$bspp*simcoef$gddcons)+
  simcoef$error

# prepare grouping factors
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
simcoef$treeid <- factor(simcoef$treeid)

# === === === === === #
##### Run model #####
# === === === === === #
y <- simcoef$ringwidth
N <- nrow(simcoef)
gdd <- simcoef$gddcons
Nspp <- length(unique(simcoef$spp))
Nsite <- length(unique(simcoef$site))
site <- as.numeric(as.character(simcoef$site))
species <- as.numeric(as.character(simcoef$spp))
treeid <- treeid
Ntreeid <- length(unique(treeid))
table(treeid)

rstan_options(auto_write = TRUE)

fit <- stan("stan/twolevelhierint.stan", 
                    data=c("N","y",
                           "Nspp","species",
                           "Nsite", "site", 
                           "Ntreeid", "treeid", 
                           "gdd"),
                    iter=4000, chains=4, cores=4)

# saveRDS(fit, "output/stanOutput/GDDleafout/fit")
# fit <- readRDS("output/stanOutput/GDDleafout/fit")

fit_with_b <- stan("stan/twolevelhierint_only_inter_bsp.stan", 
                   data=c("N","y", 
                          "Ntreeid", "treeid",
                          "Nspp","species",
                          "Nsite","site",
                          "gdd"),
                   iter=4000, chains=4, cores=4)
# === === === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)

if (FALSE){

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
unique(colnames(df_fit))
sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]

sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]

sigma_df2 <- data.frame(
  sigma = character(ncol(sigma_df)),
  mean = numeric(ncol(sigma_df)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)
sigma_df2

for (i in 1:ncol(sigma_df)) { # i = 1
  sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
  sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
  sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.05), 3)
  sigma_df2$per25[i] <- round(quantile(sigma_df[[i]], probs = 0.25), 3)
  sigma_df2$per75[i] <- round(quantile(sigma_df[[i]], probs = 0.75), 3)
  sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.95), 3)
}

sigma_df2$sim_sigma <- c(
  sigma_bspp,
                         sigma_aspp, 
                         sigma_asite, 
                         sigma_atreeid, 
                         sigma_y)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover b spp ######
bspp_cols <- colnames(df_fit)[grepl("bsp", colnames(df_fit))]
# remove sigma_aspp for now
bspp_cols <- bspp_cols[2:length(bspp_cols)]

bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# change their names
colnames(bspp_df) <- sub("bsp\\[(\\d+)\\]", "\\1", colnames(bspp_df))
#empty spp df
bspp_df2 <- data.frame(
  spp = character(ncol(bspp_df)),
  fit_bspp = numeric(ncol(bspp_df)),  
  fit_bspp_per5 = NA, 
  fit_bspp_per25 = NA,
  fit_bspp_per75 = NA,
  fit_bspp_per95 = NA
)
for (i in 1:ncol(bspp_df)) { # i = 1
  bspp_df2$spp[i] <- colnames(bspp_df)[i]         
  bspp_df2$fit_bspp[i] <- round(mean(bspp_df[[i]]),3)  
  bspp_df2$fit_bspp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
  bspp_df2$fit_bspp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
  bspp_df2$fit_bspp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
  bspp_df2$fit_bspp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
}



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
# remove sigma_aspp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_atreeid = numeric(ncol(treeid_df)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_atreeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_atreeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_atreeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_atreeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_atreeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a spp  ######
aspp_cols <- colnames(df_fit)[grepl("aspp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("zaspp", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("sigma", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("aspp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_aspp = numeric(ncol(aspp_df)),  
  fit_aspp_per5 = NA, 
  fit_aspp_per25 = NA,
  fit_aspp_per75 = NA,
  fit_aspp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_aspp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_aspp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_aspp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_aspp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_aspp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]
# remove sigma_asp for now
site_cols <- site_cols[2:length(site_cols)]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_asite = numeric(ncol(site_df)),  
  fit_asite_per5 = NA, 
  fit_asite_per25 = NA,
  fit_asite_per75 = NA,
  fit_asite_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_asite[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_asite_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_asite_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_asite_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_asite_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2

# === === === === === === === #
##### Plot parameter recovery #####
# === === === === === === === #

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot sigmas ######
sigma_simXfit_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_errorbar(aes(ymin = per25, ymax = per75),
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 1) +
  geom_errorbar(aes(ymin = per5, ymax = per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  geom_point(color = "#046C9A", size = 3) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "") +
  theme_minimal()
sigma_simXfit_plot
ggsave("figures/sigma_simXfit_plot.jpeg", sigma_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot b spp ######
bspptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "bspp")], 
  bspp_df2[!duplicated(bspp_df2$spp), 
           c("spp", "fit_bspp", "fit_bspp_per25", "fit_bspp_per75", "fit_bspp_per5", "fit_bspp_per95")], 
  by = "spp"
)
bspptoplot

bspp_simXfit_plot <- ggplot(bspptoplot, aes(x = bspp, y = fit_bspp)) +
  geom_errorbar(aes(ymin = fit_bspp_per5, ymax = fit_bspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_errorbar(aes(ymin = fit_bspp_per25, ymax = fit_bspp_per75), 
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 0.7) +
  geom_point(color = "#046C9A", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim bspp", y = "fit bspp", title = "") +
  theme_minimal()
bspp_simXfit_plot
# ggsave!
ggsave("figures/bspp_simXfit_plot2.jpeg", bspp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
# add sim to fit treeid df
treeidtoplot <- merge(
  simcoef[!duplicated(simcoef$treeid), 
          c("treeid", "atreeid")], 
  treeid_df2[!duplicated(treeid_df2$treeid), 
             c("treeid", "fit_atreeid", 
               "fit_atreeid_per5", 
               "fit_atreeid_per25",
               "fit_atreeid_per75",
               "fit_atreeid_per95")], 
  by = "treeid"
)
treeidtoplot
# plot treeid
atreeid_simXfit_plot <- ggplot(treeidtoplot, aes(x = atreeid, y = fit_atreeid)) +
  geom_errorbar(aes(ymin = fit_atreeid_per5, ymax = fit_atreeid_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "") +
  theme_minimal()
atreeid_simXfit_plot
# ggsave!
ggsave("figures/atreeid_simXfit_plot.jpeg", atreeid_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot a spp ######
aspptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "aspp")], 
  aspp_df2[!duplicated(aspp_df2$spp), 
           c("spp", "fit_aspp", 
             "fit_aspp_per5", 
             "fit_aspp_per25", 
             "fit_aspp_per75", 
             "fit_aspp_per95")], 
  by = "spp"
)
aspptoplot

aspp_simXfit_plot <- ggplot(aspptoplot, aes(x = aspp, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per5, ymax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim aspp", y = "fit aspp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
aspp_simXfit_plot
# ggsave!
ggsave("figures/aspp_simXfit_plot.jpeg", aspp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######
sitetoplot <- merge(
  simcoef[!duplicated(simcoef$site), 
          c("site", "asite")], 
  site_df2[!duplicated(site_df2$site), 
           c("site", "fit_asite", 
             "fit_asite_per5", 
             "fit_asite_per25", 
             "fit_asite_per75", 
             "fit_asite_per95")], 
  by = "site"
)
sitetoplot

asite_simXfit_plot <- ggplot(sitetoplot, aes(x = asite, y = fit_asite)) +
  geom_errorbar(aes(ymin = fit_asite_per5, ymax = fit_asite_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_asite_per25, ymax = fit_asite_per75), 
                width = 0, linewidth = 1.5, color = "darkgray", alpha=0.9) +
  geom_point(color = "#046C9A", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim asite", y = "fit asite", title = "") +
  theme_minimal()
asite_simXfit_plot
# ggsave!
ggsave("figures/asite_simXfit_plot.jpeg", asite_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

###### Combine plots  ######
combined_plot <- (atreeid_simXfit_plot) /
  (sigma_simXfit_plot + bspp_simXfit_plot ) /
  (aspp_simXfit_plot + asite_simXfit_plot)
combined_plot
ggsave("figures/combinedPlots.jpeg", combined_plot, width = 10, height = 8, units = "in", dpi = 300)

#  === === === === === === === === === === === === === === === === 

}

##### Diagnostics #####
# === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
###### Priors VS Posterior ######
# === === === === === === === === === === === === === === === === 
# prior predictive checks. Simulating prior values from the values set in the model block
hyperparameter_draws <- 8000
parameter_draws <- 1000

###### Priors sigma_bsp ######
sigma_bsp_draw <- abs(rnorm(hyperparameter_draws, 0, 0.2))   
ggplot() +
  geom_density(data = data.frame(sigma_bsp_draw = sigma_bsp_draw),
               aes(x = sigma_bsp_draw, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = sigma_df,
               aes(x = sigma_bsp, colour = "Posterior"),
               linewidth = 0.8) +
  labs(title = "priorVSposterior_bsp",
       x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
ggsave("figures/priorsPredictiveChecks/priorVSposterior_sigma_bsp.jpeg", width = 10, height = 8, units = "in", dpi = 300)

sigma_aspp_draw <- abs(rnorm(draws, 0, 0.5))
sigma_asite_draw <- abs(rnorm(draws, 0, 0.5))
sigma_atree_draw <- abs(rnorm(draws, 0, 0.05))
sigma_y_draw <- abs(rnorm(draws, 0, 5))

######Priors bsp ######
n_sigma_bsp <- 200

# set to prior values
sigma_bsp_vec <- abs(rnorm(n_sigma_bsp, 0, 0.2))

prior_bsp <- rep(NA, parameter_draws*length(sigma_bsp_vec))

for (i in 1: length(sigma_bsp_vec)) {
  prior_bsp[((i - 1)*parameter_draws + 1):(i*parameter_draws)] <- rnorm(parameter_draws, 0, sigma_bsp_vec[i])
}
prior_bsp

# Get the posterior distribution
bspp_df3 <- bspp_df

bspp_df3$draw <- rownames(bspp_df3)

colnames(bspp_df3) <- c(paste0("spp", 1:10), "draw") 
bspp_df3

long_post_bspp <- reshape(
  bspp_df3,
  direction = "long",
  varying = paste0("spp", 1:10),
  v.names = "post_bsp",
  idvar = "draw",
  timevar = "spp"
)
long_post_bspp

ggplot() +
  geom_density(data = data.frame(prior_bsp = prior_bsp),
               aes(x = prior_bsp, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = long_post_bspp,
               aes(x = post_bsp, colour = "Posterior", group = spp),
               linewidth = 0.3) +
  labs(title = "priorVSposterior_bsp",
       x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
ggsave("figures/priorsPredictiveChecks/priorVSposterior_bsp.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# now add row for prior_bsp
prior_bsp <- rnorm(nrow(sigma_df), 0, sigma_df$prior_sigma_bsp)

# convert each parameter to long format
sigma_long_bsp <- data.frame(
  value  = c(sigma_df$post_sigma_bsp, sigma_df$prior_sigma_bsp),
  source = rep(c("post_sigma_bsp", "prior_sigma_bsp"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_bsp, aes(x = value, color = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

# Sigma aspp
sigma_long_aspp <- data.frame(
  value  = c(sigma_df$post_sigma_aspp, sigma_df$prior_sigma_aspp),
  source = rep(c("post_sigma_aspp", "prior_sigma_aspp"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_aspp, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# Sigma asite
sigma_long_asite <- data.frame(
  value  = c(sigma_df$post_sigma_asite, sigma_df$prior_sigma_asite),
  source = rep(c("post_sigma_asite", "prior_sigma_asite"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_asite, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# Sigma atreeid
sigma_long_atreeid <- data.frame(
  value  = c(sigma_df$post_sigma_atreeid, sigma_df$prior_sigma_atreeid),
  source = rep(c("post_sigma_atreeid", "prior_sigma_atreeid"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_atreeid, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# sigma_y
sigma_long_atreeid <- data.frame(
  value  = c(sigma_df$post_sigma_atreeid, sigma_df$prior_sigma_y),
  source = rep(c("post_sigma_atreeid", "prior_sigma_atreeid"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_atreeid, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

}

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
emp <- read.csv("output/empiricalDataMAIN.csv")

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
gdd <- emp$pgsGDD/200
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
fit <- stan("stan/twolevelhierint.stan", 
            data=c("N","y",
                   "Nspp","species",
                   "Nsite", "site", 
                   "Ntreeid", "treeid", 
                   "gdd"),
            iter=4000, chains=4, cores=4)

# fit stanlmer to check differences
emp$gdd <- emp$pgsGDD/200
emp$y <- emp$lengthCM*10

emp$site_fac <- as.factor(emp$site_num)
emp$spp_fac <- as.factor(emp$spp_num)
emp$treeid_fac <- as.factor(emp$treeid_num)

unique(emp$spp_fac)
emp$gdd_c <- scale(emp$gdd, scale = FALSE)

fitlmer <- stan_lmer(
  y ~ 
    1 +                           
    (1|site_fac) +                    
    (1|spp_fac) +                     
    gdd:spp_fac +                 
    (1 | spp_fac:treeid_fac),     
  data = emp,
  prior = normal(0, 1),
  prior_intercept = normal(5, 6),
  chains = 4,
  adapt_delta = 0.99,
  iter = 4000,
  cores = 4 
)

fitlmer <- stan_lmer(
  y ~ 
    1 +                           
    site_fac +                    
    spp_fac +                     
    gdd:spp_fac +                 
    (1 | spp_fac:treeid_fac),     
  data = emp,
  prior = normal(0, 1),
  prior_intercept = normal(5, 6),
  chains = 4,
  adapt_delta = 0.99,
  iter = 4000,
  cores = 4 
)

colnames(as.data.frame(fitlmer))[1:20]
fitlmer$stanfit
pairs(fit, pars = c("a", "b",
                    "sigma_atreeid",
                    "sigma_y"))

saveRDS(fit, "output/stanOutput/fit")
# check warnings
diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# jpeg("figures/pairs.jpg", width = 5000, height = 5000, 
     # units = "px", res = 300)

fit@model_pars

# pairs(fit, pars = c("a", "b",
#                     "sigma_atreeid",
#                     "sigma_y"))

# === === === === === === === === === === === === #
##### Recover parameters from the posterior ##### 
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
unique(colnames(df_fit))
sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]

sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]

sigma_df2 <- data.frame(
  sigma = character(ncol(sigma_df)),
  mean = numeric(ncol(sigma_df)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)
sigma_df2

for (i in 1:ncol(sigma_df)) { # i = 1
  sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
  sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
  sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.05), 3)
  sigma_df2$per25[i] <- round(quantile(sigma_df[[i]], probs = 0.25), 3)
  sigma_df2$per75[i] <- round(quantile(sigma_df[[i]], probs = 0.75), 3)
  sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.95), 3)
}
sigma_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover b spp ######
bspp_cols <- colnames(df_fit)[grepl("bsp", colnames(df_fit))]
# remove sigma_bspp for now
# bspp_cols <- bspp_cols[2:length(bspp_cols)]

bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# change their names
colnames(bspp_df) <- sub("bsp\\[(\\d+)\\]", "\\1", colnames(bspp_df))
#empty spp df
bspp_df2 <- data.frame(
  spp = character(ncol(bspp_df)),
  fit_bspp = numeric(ncol(bspp_df)),  
  fit_bspp_per5 = NA, 
  fit_bspp_per25 = NA,
  fit_bspp_per75 = NA,
  fit_bspp_per95 = NA
)
for (i in 1:ncol(bspp_df)) { # i = 1
  bspp_df2$spp[i] <- colnames(bspp_df)[i]         
  bspp_df2$fit_bspp[i] <- round(mean(bspp_df[[i]]),3)  
  bspp_df2$fit_bspp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
  bspp_df2$fit_bspp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
  bspp_df2$fit_bspp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
  bspp_df2$fit_bspp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
}
bspp_df2


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]
treeid_cols <- treeid_cols[!grepl("sigma", treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]

# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_atreeid = numeric(ncol(treeid_df)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_atreeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_atreeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_atreeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_atreeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_atreeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a spp  ######
aspp_cols <- colnames(df_fit)[grepl("aspp", colnames(df_fit))]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("aspp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_aspp = numeric(ncol(aspp_df)),  
  fit_aspp_per5 = NA, 
  fit_aspp_per25 = NA,
  fit_aspp_per75 = NA,
  fit_aspp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_aspp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_aspp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_aspp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_aspp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_aspp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_asite = numeric(ncol(site_df)),  
  fit_asite_per5 = NA, 
  fit_asite_per25 = NA,
  fit_asite_per75 = NA,
  fit_asite_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_asite[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_asite_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_asite_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_asite_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_asite_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2

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
bsp_prior <- rnorm(1e4, 0, 0.4)

priorbsp <- ggplot() +
  geom_density(data = data.frame(bsp_prior = bsp_prior),
               aes(x = bsp_prior, colour = "Prior at N(0, 0.4)"),
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
df_fit <- as.data.frame(fitlmer)
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

