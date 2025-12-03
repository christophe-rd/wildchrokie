
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

gdd <- emp$pgsGDD/200

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

fit_with_b <- stan("stan/twolevelhierint_only_atreeid_no_b.stan", 
            data=c("N","y", 
                   "Ntreeid", "treeid",
                   "Nspp","species",
                   "Nsite","site",
                   "gdd"),
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


# Recover parameters ####
###### treeid ######
df_fit <- as.data.frame(fit)

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

###### aspp  ######
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

###### asite ######
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

# Prior checks ####
##### atreeid #####
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

##### asp #####
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


##### asite #####
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


##### a #####
df_fit <- as.data.frame(fit)

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

# Model comparison ####
##### recover the model with bsp first #####
###### treeid ######
df_fit_with_b <- as.data.frame(fit_with_b)

# grab treeid 
treeid_cols <- colnames(df_fit_with_b)[grepl("atreeid", colnames(df_fit_with_b))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]

# remove sigma_asp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

treeid_df <- df_fit_with_b[, colnames(df_fit_with_b) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2_with_b <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_a_treeid = numeric(ncol(treeid_df)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2_with_b$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2_with_b$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2_with_b$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2_with_b$fit_a_treeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2_with_b$fit_a_treeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2_with_b$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df

###### aspp  ######
aspp_cols <- colnames(df_fit_with_b)[grepl("asp", colnames(df_fit_with_b))]
aspp_cols <- aspp_cols[!grepl("zasp", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("sigma", aspp_cols)]

aspp_df <- df_fit_with_b[, colnames(df_fit_with_b) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("asp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2_with_b <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_a_spp = numeric(ncol(aspp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per25 = NA,
  fit_a_spp_per75 = NA,
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2_with_b$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2_with_b$fit_a_spp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2_with_b$fit_a_spp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2_with_b$fit_a_spp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2_with_b$fit_a_spp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2_with_b$fit_a_spp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2_with_b

###### asite ######
site_cols <- colnames(df_fit_with_b)[grepl("asite", colnames(df_fit_with_b))]

site_df <- df_fit_with_b[, colnames(df_fit_with_b) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2_with_b <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per25 = NA,
  fit_a_site_per75 = NA,
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2_with_b$site[i] <- colnames(site_df)[i]         
  site_df2_with_b$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2_with_b$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2_with_b$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2_with_b$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2_with_b$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2_with_b

##### Plots #####
###### asp ######
aspp_df2$spp <- as.numeric(aspp_df2$spp)
aspp_df2$spp_name <- emp$spp[match(aspp_df2$spp, emp$spp_num)]

aspp_df2_with_b$spp <- as.numeric(aspp_df2_with_b$spp)
aspp_df2_with_b$spp_name <- emp$spp[match(aspp_df2_with_b$spp, emp$spp_num)]

aspp_df2$model <- "intercept only"
aspp_df2_with_b$model <- "intercept and b"

aspp_df_binded <- rbind(aspp_df2, aspp_df2_with_b)

pos <- position_jitter(height = 0.15)

##### aspp intercepts means #####
ggplot(aspp_df_binded, aes(x = fit_a_spp, y = spp_name, color = model)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_spp_per5, 
                     xmax = fit_a_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_spp_per25, 
                     xmax = fit_a_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Model parameters") +
  facet_wrap(~ model, nrow =2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_bw() +
  scale_y_discrete(limits = rev)  
ggsave("figures/interceptVSinterceptslope_spp.jpeg", width = 8, height = 6, 
       units = "in", dpi = 300)


###### asite ######
site_df2$site <- as.numeric(site_df2$site)
site_df2$site_name <- emp$site[match(site_df2$site, emp$site_num)]

site_df2_with_b$site <- as.numeric(site_df2_with_b$site)
site_df2_with_b$site_name <- emp$site[match(site_df2_with_b$site, emp$site_num)]

site_df2$model <- "intercept only"
site_df2_with_b$model <- "intercept and b"

asite_df_binded <- rbind(site_df2, site_df2_with_b)

##### asite intercepts means #####
ggplot(asite_df_binded, aes(x = fit_a_site, y = site_name, color = model)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_site_per5, 
                     xmax = fit_a_site_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_site_per25, 
                     xmax = fit_a_site_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Site", x = "Site intercept values", color = "Model parameters") +
  facet_wrap(~ model, nrow =2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_bw() +
  scale_y_discrete(limits = rev)  
ggsave("figures/interceptVSinterceptslope_site.jpeg", width = 8, height = 6, 
       units = "in", dpi = 300)
