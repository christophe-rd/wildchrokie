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
fit1 <- fit1 <- fit1 <- FALSE
fit2 <- FALSE

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")


# === === === === === === === === === === === === === === === === 
#### Step 1. Come up with a model ####
# === === === === === === === === === === === === === === === === 
# As discussed with Victor, we will use a linear model to predict growth from GDD. 

# First, I'll try with the simplest version of the model, using whole year GDD instead of the growing degree days between leaf out and budset

gdd <- read.csv("output/gddData.csv")

# Simulate ring width for each year, from 2017 to 2023 for 4 species

# Start with Alnus incana. Let's go with a sample size of 100 cookies that have ring width data for each year

### assign more digestible names and years to help me understand what the hell im doing
spp <- c("alninc","betpap","betpop","betall")
plot <- c("wm1","wm2","wm3","hf1","hf2","hf3","pg1","pg2","pg3","sh1","sh2","shg3")
yr <- 2016:2023
rep <- 1:3

Nspp <- length(spp) # number of species
Nplot <- length(plot) # number of plot per species
Nyr <- length(yr) # number of growth years per tree
Nrep <- length(rep) # number of measurements per tree

# First making a data frame for the growth ring data
Nw <- Nplot*Nyr*Nrep*Nspp # number of measurements per species
Nw

#make a dataframe for ring width
wdat <- data.frame(matrix(NA, Nw, 1))
wdat$yr <- rep(c(yr), each=Nrep)
wdat$plot <- rep(c(plot), each=Nyr*Nrep)
wdat$spp <- rep(c(spp), each=Nyr*Nrep*Nplot)

# paste plot and spp to give unique ids
wdat$treeid <- paste(wdat$spp, wdat$plot, sep = "_")

mu.grand <- 0.4 # the grand mean of ring width

# now generating the year effect with data
sigma.yr <- 0.05  # 
mu.yr <- rnorm(Nyr, 0, sigma.yr) #intercept for each study
wdat$mu.yr <- rep(mu.yr, each = Nrep) # generate data for ea study

# now generating the tree cookie effect within a same species
sigma.tree <- 0.02 # setting it low for now
mu.tree <- rnorm(Nplot, 0, sigma.tree) 
wdat$mu.tree <- rep(mu.tree, each=Nyr*Nrep) # add width data for each species
wdat

#now generating the species effect
sigma.spp <- sigma.yr*1.2 # seting it for now as twice of the year effect 
mu.spp <- rnorm(Nspp, 0, sigma.spp)
wdat$mu.spp <- rep(mu.spp, each=Nyr*Nrep*Nplot) # add width data for each species

# general variance
w.var <- 0.5 #will have to change
wdat$w.er <- rnorm(Nw, 1, w.var)

# generate yhat aka (I guess) yhat?
wdat$ringwidth <- mu.grand + wdat$mu.yr +wdat$mu.tree + wdat$mu.spp + wdat$w.er
wdat

wdat$spp <- factor(wdat$spp)
wdat$plot  <- factor(wdat$plot)
wdat$treeid <- factor(wdat$treeid)
wdat$yr <- factor(wdat$yr)

fit <- stan_lmer(ringwidth ~ spp + # main effects
                   (1 | yr) + (1 | treeid), # random effects
                 data = wdat,
                 chains = 2, 
                 cores = 2, 
                 iter = 4000, 
                 warmup=2000,
                 control = list(adapt_delta = 0.999))

print(fit, digits=2)

mu_gdd <- 94000 # mean of GDD
a <- 10 # a for alpha: intercept whatever, I have no clue what to set it at right now
b <- 0.0001 # b for beta: slope whatever. Assuming, that for each DGG, growth ring rises by that value
sigma_gdd <- 10 # standard deviation

# Find GDD. For now I'll set it as it starts on DOY 100 and ends on DOY 250. 
x <- rnorm(n, mu_gdd, sigma_gdd)
# hist(x)
# get ring width
y <- a +b*x + sigma_gdd*rnorm(n)
hist(y)
plot(y~x)

fake <- data.frame(x,y)

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