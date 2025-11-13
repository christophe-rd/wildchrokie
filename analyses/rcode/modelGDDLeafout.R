# Wildchrokie model -- leaf out prediction
# CRD 13 November 2025
# Estimate GDD until leafout (or budburst) using an intercept only model with species, ID and provenance.

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

# simulate data
set.seed(124)
a <- 150
sigma_y <- 6
sigma_asp <- 8
sigma_atreeid <- 2
sigma_asite <- 4

n_site <- 10 # number of sites
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

sim <- data.frame(
  site = site,
  spp = spp,
  treeid = treeid
)

# get intercept values for each parameter
asp <- rnorm(n_spp, 0, sigma_asp)
asite <- rnorm(n_site, 0, sigma_asite)
atreeid <- rnorm(n_treeid, 0, sigma_atreeid)

# Add my parameters to the df
sim$atreeid <- atreeid[treeid]
sim$asite <- asite[sim$site]
sim$asp <- asp[sim$spp]

# add the rest
sim$a <- a
sim$sigma_y <- sigma_y

sim$sigma_atreeid <- sigma_atreeid
sim$sigma_asp <- sigma_asp
sim$sigma_asite <- sigma_asite
sim$error <- rnorm(N, 0, sigma_y)

# adding both options of tree rings
sim$gddleafout <- 
  sim$asite + 
  sim$asp + 
  sim$atreeid +
  sim$a +
  sim$error


hist(sim$gddleafout)
emp <- read.csv("output/empiricalDataMAIN.csv")

gdd <- read.csv("output/gddByYear.csv")

gdd$year <- as.factor(gdd$year)
ggplot(gdd)  +
  geom_point(aes(x = doy, y = GDD_10, color = year)) + 
  geom_vline(xintercept = mean(emp$leafout))
