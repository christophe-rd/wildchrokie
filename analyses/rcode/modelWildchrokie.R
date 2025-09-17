# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 3)
quartz()

# Load library 
library(rstanarm)
library(ggplot2)
library(rstan)
library(shinystan)
library("wesanderson")

fitalpha = FALSE
fitbeta = TRUE

if(length(grep("christophe_rouleau-desrochers", getwd()) > 0)) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if(length(grep("lizzie", getwd())) > 0){
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
}

# === === === === === === === === === === === === === === === === 
#### Step 1. Come up with a model ####
# === === === === === === === === === === === === === === === === 
# As discussed with Victor, we will use a linear model to predict growth from GDD. 

# === === === === === === === === === === === === === === === === 
#### Step 2. Simulate data ####
# === === === === === === === === === === === === === === === === 
# set parameters

# set parameters
set.seed(124)
a <- 1.5
b <- 0.4
sigma_y <- 0.2
sigma_a_spp <- 0.3 # This is pretty low, but I guess you think your species are closely related and will be similar?
sigma_a_treeid <- 0.5
sigma_a_site <- 0.1
sigma_b_spp <- 0.2

n_site <- 4 # number of sites
n_spp <- 10 # number of species
n_perspp <- 10 # number of individuals per species
n_treeid <- n_perspp * n_spp * n_site # number of treeid
n_meas <- 5 # repeated measurements per id
N <- n_treeid * n_meas # total number of measurements


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
a_spp <- rnorm(n_spp, 0, sigma_a_spp)
a_site <- rnorm(n_site, 0, sigma_a_site)
a_treeid <- rnorm(n_treeid, 0, sigma_a_treeid)

# get slope values for each species
b_spp <- rnorm(n_spp, 0, sigma_b_spp)

# Add my parameters to the df
simcoef$a_treeid <- a_treeid[treeid]
simcoef$a_site <- a_site[simcoef$site]
simcoef$a_spp <- a_spp[simcoef$spp]
simcoef$b_spp <- b_spp[simcoef$spp]

# add the rest of the boring stuff 
simcoef$a <- a
simcoef$b <- b
simcoef$sigma_y <- sigma_y
simcoef$sigma_a_treeid <- sigma_a_treeid
simcoef$sigma_a_spp <- sigma_a_spp
simcoef$sigma_a_site <- sigma_a_site
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gddcons <- rnorm(N, 1800, 100)/200

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$a_site + 
  simcoef$a_spp + 
  simcoef$a_treeid + 
  simcoef$a + 
  (simcoef$b*simcoef$gddcons) + 
  (simcoef$b_spp*simcoef$gddcons)+
  simcoef$error

# prepare grouping factors
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
simcoef$treeid <- factor(simcoef$treeid)

# === === === === === #
##### Run models #####
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

fit <- rstan::stan("stan/twolevelhierint.stan", 
                      data=c("N","y","Nspp","species","Nsite", "site", "Ntreeid", "treeid", "gdd"),
                      iter=4000, chains=4, cores=4)  

summary(fit)$summary

# launch_shinystan(fit2)

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
  per95 = NA
)
sigma_df2
for (i in 1:ncol(sigma_df)) { # i = 1
  sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
  sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
  sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.055), 3)
  sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.945), 3)
}

sigma_df2$sim_sigma <- c(sigma_b_spp, sigma_a_spp, sigma_a_site, sigma_a_treeid, sigma_y)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
if(fitbeta){
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
    fit_b_spp = numeric(ncol(bspp_df)),  
    fit_b_spp_per5 = NA, 
    fit_b_spp_per95 = NA
  )
  bspp_df2
  for (i in 1:ncol(bspp_df)) { # i = 1
    bspp_df2$spp[i] <- colnames(bspp_df)[i]         
    bspp_df2$fit_b_spp[i] <- round(mean(bspp_df[[i]]),3)  
    bspp_df2$fit_b_spp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.055), 3)
    bspp_df2$fit_b_spp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.945), 3)
  }
  bspp_df2
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
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
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.055), 3)
  treeid_df2$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.945), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover spp  ######
spp_cols <- colnames(df_fit)[grepl("asp", colnames(df_fit))]
# remove sigma_asp for now
spp_cols <- spp_cols[2:length(spp_cols)]

spp_df <- df_fit[, colnames(df_fit) %in% spp_cols]
# change their names
colnames(spp_df) <- sub("asp\\[(\\d+)\\]", "\\1", colnames(spp_df))
#empty spp df
spp_df2 <- data.frame(
  spp = character(ncol(spp_df)),
  fit_a_spp = numeric(ncol(spp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per95 = NA
)
spp_df2
for (i in 1:ncol(spp_df)) { # i = 1
  spp_df2$spp[i] <- colnames(spp_df)[i]         
  spp_df2$fit_a_spp[i] <- round(mean(spp_df[[i]]),3)  
  spp_df2$fit_a_spp_per5[i] <- round(quantile(spp_df[[i]], probs = 0.055), 3)
  spp_df2$fit_a_spp_per95[i] <- round(quantile(spp_df[[i]], probs = 0.945), 3)
}
spp_df2


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]
# remove sigma_asp for now
site_cols <- site_cols[2:length(site_cols)]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.055), 3)
  site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.945), 3)
}
site_df2

# === === === === === === === #
# Plot parameter recovery #####
# === === === === === === === #

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot sigmas ######
sigma_simXfit_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_point(color = "#046C9A", size = 3) +
  geom_errorbar(aes(ymin = per5, ymax = per95),
                width = 0, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "fit vs sim sigmas") +
  theme_minimal()
sigma_simXfit_plot
ggsave("figures/sigma_simXfit_plot.jpeg", sigma_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot b spp ######
bspptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "b_spp")], 
  bspp_df2[!duplicated(bspp_df2$spp), 
           c("spp", "fit_b_spp", "fit_b_spp_per5", "fit_b_spp_per95")], 
  by = "spp"
)
bspptoplot

b_spp_simXfit_plot <- ggplot(bspptoplot, aes(x = b_spp, y = fit_b_spp)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_b_spp_per5, ymax = fit_b_spp_per95), width = 0, color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim b_spp", y = "fit b_spp", title = "") +
  theme_minimal()
b_spp_simXfit_plot
# ggsave!
ggsave("figures/b_spp_simXfit_plot.jpeg", b_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
# add sim to fit treeid df
treeidtoplot <- merge(
  simcoef[!duplicated(simcoef$treeid), 
          c("treeid", "a_treeid")], 
  treeid_df2[!duplicated(treeid_df2$treeid), 
         c("treeid", "fit_a_treeid", "fit_a_treeid_per5", "fit_a_treeid_per95")], 
  by = "treeid"
)
treeidtoplot
# plot treeid
a_treeid_simXfit_plot <- ggplot(treeidtoplot, aes(x = a_treeid, y = fit_a_treeid)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_treeid_per5, ymax = fit_a_treeid_per95), width = 0, color = "darkgray", alpha=0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_treeid", y = "fit a_treeid", title = "") +
  theme_minimal()
a_treeid_simXfit_plot
# ggsave!
ggsave("figures/a_treeid_simXfit_plot.jpeg", a_treeid_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot spp ######
spptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "a_spp")], 
  spp_df2[!duplicated(spp_df2$spp), 
         c("spp", "fit_a_spp", "fit_a_spp_per5", "fit_a_spp_per95")], 
  by = "spp"
  )
spptoplot

a_spp_simXfit_plot <- ggplot(spptoplot, aes(x = a_spp, y = fit_a_spp)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_spp_per5, ymax = fit_a_spp_per95), width = 0, color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_spp", y = "fit a_spp", title = "") +
  theme_minimal()
a_spp_simXfit_plot
# ggsave!
ggsave("figures/a_spp_simXfit_plot.jpeg", a_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######
sitetoplot <- merge(
  simcoef[!duplicated(simcoef$site), 
          c("site", "a_site")], 
  site_df2[!duplicated(site_df2$site), 
         c("site", "fit_a_site", "fit_a_site_per5", "fit_a_site_per95")], 
  by = "site"
)
sitetoplot

a_site_simXfit_plot <- ggplot(sitetoplot, aes(x = a_site, y = fit_a_site)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_site_per5, ymax = fit_a_site_per95), width = 0, color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_site", y = "fit a_site", title = "") +
  theme_minimal()
a_site_simXfit_plot
# ggsave!
ggsave("figures/a_site_simXfit_plot.jpeg", a_site_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)


# === === === === === === === === === === === === === === === === 
#### Step 3. Set your priors ####
# === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
#### Step 4. Run model on empirical data ####
# === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
#### Step 5. Perform retrodictive checks using the model fit to your empiral data ####
# === === === === === === === === === === === === === === === ===


