# 2 September
# By Christophe RD
# Figure out nesting a_ids wihtin species

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 

# Load library 
library(rstanarm)
library(ggplot2)
library(wesanderson)

runmodels <- FALSE
runoldcode <- FALSE

# === === === === === #
# Simulate data ####
# === === === === === #

# set parameters
set.seed(124)
a <- 1.5
sigma_y <- 0.2
sigma_a_spp <- 0.6
sigma_a_ids <- 0.1

n_spp <- 10 # number of species
n_perspp <- 50 # number of individuals per species
n_ids <- n_perspp * n_spp # number of ids
n_meas <- 15 # repeated measurements per id
N <- n_ids * n_meas # total number of measurements

# set spp names and ids 
spp <- 1:n_spp
ids <- rep(rep(1:n_perspp, times = n_spp), each = n_meas)
idsnonrep <- rep(rep(1:n_perspp, times = n_spp))
tree_spp <- rep(rep(1:n_spp, each = n_perspp), each = n_meas) # 1500 obs
tree_spp_nonrep <- rep(rep(1:n_spp, each = n_perspp)) # 500 obs

# get 50 intercept values for each species
a_spp_values <- rnorm(n_spp, 0, sigma_a_spp)

# === === === === === #
# Option 1 #
# === === === === === #
a_ids_spp_values <- rnorm(n_ids, a_spp_values[tree_spp_nonrep], sigma_a_ids)

# === === === === === #
# Option 2 #
# === === === === === #
a_ids <- rnorm(n_ids, 0, sigma_a_ids)

testdf <- data.frame(
  tree_spp_nonrep = tree_spp_nonrep,
  idsnonrep = idsnonrep,
  a_spp_values = a_spp_values[tree_spp_nonrep],
  a_ids = a_ids,
  a_ids_spp_values = a_ids_spp_values,
  a_ids_spp_values2 = a_spp_values[tree_spp_nonrep] + a_ids
)
testdf


a_ids_spp1 <- ggplot(testdf) +
  geom_vline(aes(xintercept = a_ids_spp_values), 
             color = "red", alpha = 0.2) +
  geom_vline(aes(xintercept = a_spp_values), 
             color = "blue", alpha = 0.8) +
  facet_wrap(~tree_spp_nonrep) +
  theme_minimal()
ggsave("figures/a_ids_spp1.jpeg", a_ids_spp1, width = 8, height = 6)

a_ids_spp2 <- ggplot(testdf) +
  geom_vline(aes(xintercept = a_ids_spp_values2), 
             color = "red", alpha = 0.2) +
  geom_vline(aes(xintercept = a_spp_values), 
             color = "blue", alpha = 0.8) +
  facet_wrap(~tree_spp_nonrep) +
  theme_minimal()
ggsave("figures/a_ids_spp2.jpeg", a_ids_spp2, width = 8, height = 6)

