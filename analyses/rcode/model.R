# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)

fit1 <- FALSE
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
n <- 100
nspp <- 1
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
fit1 <- stan_glm(y~x, data=fake)
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