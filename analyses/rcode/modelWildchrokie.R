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
library(arm)
library(RColorBrewer)
library(shinystan)
library(rethinking)
library("wesanderson")

runmodels <- FALSE
runoldcode <- FALSE

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")


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
sigma_y <- 0.2
sigma_a_spp <- 0.3
sigma_a_ids <- 0.15
sigma_a_site <- 0.1

n_site <- 4 # number of sites
n_spp <- 10 # number of species
n_perspp <- 50 # number of individuals per species
n_ids <- n_perspp * n_spp * n_site # number of ids
n_meas <- 5 # repeated measurements per id
N <- n_ids * n_meas # total number of measurements


# get replicated ids
ids <- rep(rep(rep(1:n_perspp, times = n_spp), times = n_site), each = n_meas)
# non replicated ids
idsnonrep <- rep(rep(1:n_perspp, times = n_spp), times = n_site)
# replicated spp
spp <- rep(rep(rep(1:n_spp, each = n_perspp), times = n_site), each = n_meas) 
# non replicated spp
spp_nonrep <- rep(rep(1:n_spp, each = n_perspp), each = n_site) 
# replicated site
site <- rep(rep(rep(1:n_site, each = n_spp), each = n_perspp), each = n_meas)
# non replicated site
site_nonrep <- rep(rep(1:n_site, each = n_spp), each = n_perspp)

simcoef <- data.frame(
  site = site,
  spp = spp,
  ids = ids
)

# get 50 intercept values for each species
a_spp <- rnorm(n_spp, 0, sigma_a_spp)

a_site <- rnorm(n_site, 0, sigma_a_site)

# Option 2
a_ids <- rnorm(n_ids, 0, sigma_a_ids)

# Add my parameters to the df
simcoef$a_site <- a_site[simcoef$site]
simcoef$a_spp <- a_spp[simcoef$spp]

### add index to allow ids to be input at the right place
id_index <- rep(1:n_ids, each = n_meas)

simcoef$a_ids <- a_ids[id_index]

# add the rest of the boring stuff 
simcoef$a <- 1.5
simcoef$b <- 0.4
simcoef$sigma_y <- sigma_y
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gddcons <- rnorm(N, 1800, 100)/200

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$a_site + 
  simcoef$a_spp + 
  simcoef$a_ids + 
  simcoef$a + 
  (simcoef$b*simcoef$gddcons) + 
  simcoef$error

# prepare grouping factors for stan_lmer (ensure ids are unique within spp)
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
simcoef$ids <- factor(simcoef$ids)
simcoef$ids_unique <- paste(simcoef$site, simcoef$spp, simcoef$ids, sep = "_")

# === === === === === #
##### Run models #####
# === === === === === #

###### Model nested on the intercept #######
fitnestedrun <- TRUE
if(fitnestedrun) {
  fitnested <- stan_lmer(
    ringwidth ~ 1 + gddcons + 
      (1|site) +
      (1|spp) +
      (1|ids_unique),
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
  saveRDS(fitnested, "output/fitnested")
}
fitnested

###### Model partial pooled on b ######
runmodeWithPartialPooledBeta <- FALSE
if(runmodeWithPartialPooledBeta) {
  fit <- stan_lmer(
    ringwidth ~ 1 + gddcons + 
      (1 | ids) +
      (1 | spp) +
      (1|spp/ids),  #(1 | spp) means that I am partial pooling for the spp interecept
    # to confirm: to partial pool on both the intercept AND the slope : (1 + gddcons|spp)
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
}


print(fit, digits=3)

# === === === === === === #
##### Parameter recovery #####
# === === === === === === #

# === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === #
fitnested <- readRDS("output/fitnested")
df_fit <- as.data.frame(fitnested)

# recover slope
colnames(df_fit)
# grab ids nested in spp
ids_cols <- colnames(df_fit)[grepl("ids_unique", colnames(df_fit))]
ids_cols <- ids_cols[1:length(ids_cols)-1]
ids_df <- df_fit[, colnames(df_fit) %in% ids_cols]
# change their names
colnames(ids_df) <- sub(".*ids_unique(.*)\\]$", "\\1", colnames(ids_df))
# empty ids dataframe
ids_df2 <- data.frame(
  ids_spp = character(ncol(ids_df)),
  fit_a_ids_spp = numeric(ncol(ids_df)),  
  fit_per5 = NA, 
  fit_per95 = NA,
  fit_sd = NA
)
for (i in 1:ncol(ids_df)) { # i = 1
  ids_df2$ids_spp[i] <- colnames(ids_df)[i]         
  ids_df2$fit_a_ids_spp[i] <- round(mean(ids_df[[i]]),3)  
  ids_df2$fit_per5[i] <- round(quantile(ids_df[[i]], probs = 0.055), 3)
  ids_df2$fit_per95[i] <- round(quantile(ids_df[[i]], probs = 0.945), 3)
  ids_df2$fit_sd[i] <- round(sd(ids_df[[i]]), 3)
}
ids_df2

# grab spp 
spp_cols <- colnames(df_fit)[grepl(" spp:", colnames(df_fit))]
spp_df <- df_fit[, colnames(df_fit) %in% spp_cols]
# change their names
colnames(spp_df) <- sub(".*spp:([0-9]+).*", "\\1", colnames(spp_df))
#empty spp df
spp_df2 <- data.frame(
  spp = character(ncol(spp_df)),
  fit_a_spp = numeric(ncol(spp_df)),  
  fit_per5 = NA, 
  fit_per95 = NA,
  fit_sd = NA
)
for (i in 1:ncol(spp_df)) { # i = 1
  spp_df2$spp[i] <- colnames(spp_df)[i]         
  spp_df2$fit_a_spp[i] <- round(mean(spp_df[[i]]),3)  
  spp_df2$fit_per5[i] <- round(quantile(spp_df[[i]], probs = 0.055), 3)
  spp_df2$fit_per95[i] <- round(quantile(spp_df[[i]], probs = 0.945), 3)
  spp_df2$fit_sd[i] <- round(sd(spp_df[[i]]), 3)
}
spp_df2

# === === === === === === === === === === === === === === === === 
# Plot old vs new way tp recover parameters #####
# === === === === === === === === === === === === === === === === 
colnames(a_spp_mergedwithranef) <- c("spp", "a_spp_lme4", "per5_lme4", "per95_lme4")
colnames(spp_df2) <- c("spp", "a_spp_loop", "per5_loop", "per95_loop")
recovComparison <- merge(a_spp_mergedwithranef, spp_df2, by = "spp")

recovComparison_plot <- ggplot(recovComparison, aes(x = a_spp_lme4, y = a_spp_loop)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(xmin = per5_lme4, xmax = per95_lme4), width = 0, color = "darkgray", alpha=0.5) +
  geom_errorbar(aes(ymin = per5_loop, ymax = per95_loop), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "a_spp lme4 functions", y = "a_spp manual posterior recovery", title = "") +
  theme_minimal()
# save ggplot!
ggsave("figures/recovComparison_plot.jpeg", recovComparison_plot, width = 6, height = 6, units = "in", dpi = 300)



simVSfit <- merge(simcoeftoplot, spp_df2, by = "spp")
simVSfit_plot <- ggplot(simVSfit, aes(x = sim_a_spp, y = a_spp_loop)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5_loop, ymax = per95_loop), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "a_spp sim", y = "a_spp fit", title = "") +
  theme_minimal()
# save ggplot!
simVSfit_plot
ggsave("figures/simVSfit_plot.jpeg", simVSfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# === === === === === === === === === === === === === === === === 
#### Step 3. Set your priors ####
# === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
#### Step 4. Run model on empirical data ####
# === === === === === === === === === === === === === === === === 
# read GDD data
gdd <- read.csv("output/gddData.csv")
gdd18 <- subset(gdd, year == "2018")

# sum from DOY 100 to DOY 250
test <- subset(gdd18, doy>100)
test2 <- subset(test, doy<250)
sumgdd <- sum(test2$GDD_10)
# === === === === === === === === === === === === === === === === 
#### Step 5. Perform retrodictive checks using the model fit to your empiral data ####
# === === === === === === === === === === === === === === === ===


