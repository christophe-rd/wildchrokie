# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 200) 

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)
library(dplyr)

runmodels <- FALSE
runoldcode <- FALSE
codetofix <- FALSE

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

# === === === === === === === === === === === === === === === === 
#### Step 1. Come up with a model ####
# === === === === === === === === === === === === === === === === 
# As discussed with Victor, we will use a linear model to predict growth from GDD. 

# First, I'll try with the simplest version of the model, using whole year GDD instead of the growing degree days between leaf out and budset

gdd <- read.csv("output/gddData.csv")

# Simulate ring width for each year, from 2017 to 2023 for 4 species

### assign more digestible names and years to help me understand what the hell im doing
spp <- c("alninc","betpap","betpop","betall")
plot <- c("wm1","wm2","wm3","hf1","hf2","hf3","pg1","pg2","pg3","sh1","sh2","shg3")
yr <- 2016:2023
rep <- 1:3

# set species specific slopes to depend on gdd
spp_slopes <- c(
  alninc = 0.00002, 
  betpap = 0.00005, 
  betpop = 0.00008, 
  betall = 0.00004
  )

#make a dataframe for ring width
wdat <- data.frame(
  yr = rep(yr, each = length(plot) * length(spp) * length(rep)),
  plot = rep(rep(plot, each = length(spp) * length(rep)), times = length(yr)),
  spp = rep(rep(spp, each = length(rep)), times = length(yr) * length(plot)),
  rep = rep(rep, times = length(yr) * length(plot) * length(spp))
)

# simulate gdd 
gdd_sim <- data.frame(
  year = yr,
  gdd = round(rnorm(Nyr, 180000, 10000))  # setting a very approximate cumulative mean gdd between doy 100 and 300
)

# merge wdat and gdd
wdat <- merge(wdat, gdd_sim, by.x = "yr", by.y = "year", all.x = TRUE)

# paste plot and spp to give unique ids
wdat$treeid <- paste(wdat$spp, wdat$plot, sep = "_")

# the grand mean of ring width: baseline ring width when all other effects are zero.
mu_grand <- 0.04 

# Simulate random effect to set unique intercepts for each parameter
wdat$yr_effect <- rnorm(nrow(wdat), 0, 0.05)[match(wdat$yr, yr)] # with mean of 0 and SD of 0.05
wdat$tree_effect <- rnorm(nrow(wdat), 0, 0.02)[match(wdat$plot, plot)] # SD of 0.02 across trees
wdat$spp_effect <- rnorm(nrow(wdat), 0, 0.06)[match(wdat$spp, spp)]

# Simulate ring width (including GDD effect) aka yhat?
wdat$ringwidth <- mu_grand + 
  wdat$yr_effect + 
  wdat$tree_effect + 
  wdat$spp_effect +
  spp_slopes[wdat$spp] * wdat$gdd +  # gdd effect varies by spp
  rnorm(nrow(wdat), 0, 0.5) #overall error

# run models
if(runmodels){
  fit <- stan_lmer(
      ringwidth ~ spp * gdd + (1 | treeid)+ (1 | yr),  
    data = wdat,
    chains = 4,
    iter = 2000,
    core=4
  )
print(fit, digits=6)
}


# code I need to fix:
if (codetofix){
  
  shinystan::launch_shinystan(fit)
  
  # verify how well I am returning my parameters:
  # simulated data parameters:
  mu.grand
  sim_spp_effects <- unique(wdat[, c("spp", "mu.spp")])
  sim_tree_effects <- unique(wdat[, c("treeid", "mu.tree")])
  sim_yr_effects <- unique(wdat[, c("yr", "mu.yr")])
  sim_intercept <- mu.grand + sim_spp_effects$mu.spp[sim_spp_effects$spp == "alninc"]
  sim_resid_sd <- sqrt(w.var)
  
  
  # pull model parameters
  # sim coef relative to alninc
  true_spp_coefs <- sim_spp_effects$mu.spp - sim_spp_effects$mu.spp[sim_spp_effects$spp == "alninc"]
  names(true_spp_coefs) <- paste0("spp", sim_spp_effects$spp)
  
  # Compare with model estimates
  data.frame(
    Parameter = c("(Intercept)", names(sim_spp_effects)[-1]), 
    Estimated = fit$coefficients[1:4],  # Assuming (Intercept), sppbetall, etc.
    True = c(sim_intercept, sim_spp_effects[-1])
  )
}


# === === === === === === === === === === === === === === === === 
#### Step 2. Simulate data ####
# === === === === === === === === === === === === === === === === 
# fit1 <- stan_glm(y~x, data=fake)
print(fit1, digits=2) # first row shows estimated intercept with its uncertainty. second row shows the estimated slope with its uncertainty. 
# in the auxiliary parameters, it shows the residual standard deviation sigma with its uncertainty

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


# old code
if (runoldcode) {
  Nspp <- length(spp) # number of species
  Nplot <- length(plot) # number of plot per species
  Nyr <- length(yr) # number of growth years per tree
  Nrep <- length(rep) # number of measurements per tree
  # First making a data frame for the growth ring data
  Nw <- Nplot*Nyr*Nrep*Nspp # number of measurements per species
}
