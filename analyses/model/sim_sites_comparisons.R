# 5 September 2025
# By Christophe RD
# Comparing likelihood methods as discussed in issue #4

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 3)

# Load library 
library(rstanarm)
library(ggplot2)
library(wesanderson)

runmodels <- FALSE
runoldcode <- FALSE

# === === === === === === === === === === === #
# Simulate generic data for both options ####
# === === === === === === === === === === === #

# set parameters
set.seed(124)
a <- 1.5
sigma_y <- 0.2
sigma_a_spp <- 0.3
sigma_a_ids <- 0.15
sigma_a_site <- 0.1

n_site <- 4 # number of sites
n_spp <- 50 # number of species
n_perspp <- 10 # number of individuals per species
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

# Option 1 
a_ids_spp <- rnorm(N, a_spp[spp_nonrep], sigma_a_ids)

# Option 2
a_ids <- rnorm(n_ids, 0, sigma_a_ids)

# Add my parameters to the df
simcoef$a_site <- a_site[simcoef$site]
simcoef$a_spp <- a_spp[simcoef$spp]

### add index to allow ids to be input at the right place
id_index <- rep(1:n_ids, each = n_meas)

simcoef$a_ids_opt1 <- a_ids_spp[id_index]
simcoef$a_ids_opt2 <- a_ids[id_index]

# add the rest of the boring stuff 
simcoef$a <- 1.5
simcoef$b <- 0.4
simcoef$sigma_y <- sigma_y
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gddcons <- rnorm(N, 1800, 100)/200

# adding both options of tree rings
simcoef$ringwidth_opt1 <- simcoef$a + 
  simcoef$a_ids_opt1 + simcoef$a_site + 
  (simcoef$b*simcoef$gddcons) + simcoef$error
simcoef$ringwidth_opt2 <- simcoef$a + 
  simcoef$a_spp + simcoef$a_ids_opt2 + simcoef$a_site + 
  (simcoef$b*simcoef$gddcons) + simcoef$error

# === === == #
# Plots ####
# === === == #

# scatter plot
ggplot(simcoef, aes(ringwidth_opt1, ringwidth_opt2, color = factor(site))) +
  geom_point(alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "Option 1", y = "Option 2", color = "Site") +
  facet_wrap(~site) +
  scale_color_manual(values = wes_palettes$Moonrise3) +
  theme_minimal()
ggsave("figures/simulated_ringwidth_scatter_plot.jpeg", width = 8, height = 6)

# prep for density plot
simcoef_long <- cbind(
  simcoef,
  stack(simcoef[c("ringwidth_opt1", "ringwidth_opt2")])
)
names(simcoef_long)[(ncol(simcoef_long)-1):ncol(simcoef_long)] <- c("ringwidth", "method")

ggplot(simcoef_long, aes(ringwidth, color = method)) +
  geom_density(size = 1.5) +  
  labs(x = "Ring Width", y = "Density") +
  scale_color_manual(values = wes_palettes$FantasticFox1[c(3, 5)]) +
  theme_minimal()
# save as jpeg!
ggsave("figures/simulated_ringwidth_density_plot.jpeg", width = 8, height = 6)

