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

# empirical data
emp <- read.csv("output/empiricalDataMAIN.csv")

# gdd data
gdd <- read.csv("output/gddByYear.csv")


# simulate data gdd at leafout
set.seed(124)
a <- 150
sigma_y <- 6
sigma_asp <- 8
sigma_atreeid <- 2
sigma_asite <- 4

n_site <- 4 # number of sites
n_spp <- 4 # number of species
n_perspp <- 5 # number of individuals per species
n_treeid <- n_perspp * n_spp * n_site # number of treeid
n_meas <- 3 # repeated measurements per id
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
sim

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
sim

hist(sim$gddleafout)
 
sim$spp <- as.factor(sim$spp)
sim$sppname <- paste("spp", sim$spp)
sim$site <- as.factor(sim$site)
sim$sitename <- paste("site", sim$site)

sim$a_asite <- sim$a + sim$asite
sim$a_asp <- sim$a + sim$asp

# plot sim data
ggplot(sim) +
  geom_hline(aes(yintercept = gddleafout, color = sitename), alpha = 0.7) +
  geom_hline(aes(yintercept = a_asite), linewidth = 0.9, alpha = 0.8) +
  geom_hline(aes(yintercept = a_asp), linewidth = 0.9, linetype = 2) +
  facet_wrap(sitename~sppname) +
  labs(y = "gdd", title = "gdd at leafout") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()
ggsave("figures/sim_gddLeafout.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# plot gdd
gdd$year <- as.factor(gdd$year)
ggplot(gdd)  +
  geom_point(aes(x = doy, y = GDD_10, color = year)) + 
  geom_vline(xintercept = mean(emp$leafout))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Run model ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
y <- sim$gddleafout
N <- nrow(sim) 
Nspp <- length(unique(sim$spp))
species <- as.numeric(as.character(sim$spp))
Nsite <- length(unique(sim$site))
site <- as.numeric(as.character(sim$site))
Ntreeid <- length(unique(sim$treeid))
treeid <- treeid


table(treeid, species)

fit <- stan("stan/twolevelhierint.stan", 
            data=c("N","y","Nspp","species","Nsite", "site", "Ntreeid", "treeid"),
            iter=4000, chains=4, cores=4)


