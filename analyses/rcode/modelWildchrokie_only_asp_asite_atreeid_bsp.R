
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

# transform my groups to numeric values
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$site_num <- match(emp$site, unique(emp$site))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
Nspp <- length(unique(emp$spp_num))
species <- as.numeric(as.character(emp$spp_num))
Nsite <- length(unique(emp$site_num))
site <- as.numeric(as.character(emp$site_num))
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))

# check that everything is fine
table(species)
emp$count <- 1
emp2 <- emp[!duplicated(emp$treeid),]
aggregate(count ~ spp, FUN = sum, data = emp2)

agg <- aggregate(count ~ treeid + spp, FUN = sum, data = emp)
str(agg)
hist(emp$treeid)

ggplot(agg) +
  geom_histogram(aes(x = count), bins = 3) +
  facet_wrap(~spp) + 
  theme_minimal() + 
  labs(x = "number of observations")
ggsave("figures/empiricalData/countSpp.jpeg", width = 6, height = 6, units = "in", dpi = 300)

rstan_options(auto_write = TRUE)
 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fit <- stan("stan/twolevelhierint_only_asp_asite_atreeid_bsp.stan", 
            data=c("N","y","Nspp","species","Nsite", 
                   "site", "Ntreeid", "treeid", "gdd"),
            iter=4000, chains=4, cores=4)

# saveRDS(fit_only_asp, "output/stanOutput/fit_only_asp")

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

# asp
asp <- names(samples)[grepl("asp", names(samples))]
asp <- asp[!grepl("sigma", asp)]

jpeg("figures/aspParameterization_only_asp_asite_atreeid_bsp.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
util$plot_div_pairs(asp, "sigma_asp", samples, diagnostics, transforms = list("sigma_asp" = 1))
dev.off()

# asite
asite <- names(samples)[grepl("asite", names(samples))]
asite <- asite[!grepl("sigma", asite)]

jpeg("figures/asiteParameterization_only_asp_asite_atreeid_bsp.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
util$plot_div_pairs(asite, "sigma_asite", samples, diagnostics, transforms = list("sigma_asite" = 1))
dev.off()

# atreeid
atreeid <- names(samples)[grepl("atreeid", names(samples))]
atreeid <- atreeid[!grepl("sigma", atreeid)]
atreeid <- atreeid[sample(length(unique(atreeid)), 21)]
pdf("figures/atreeidParameterization_only_asp_asite_atreeid_bsp.pdf", width = 6, height = 18)
util$plot_div_pairs(atreeid, "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
dev.off()

# bsp
bsp <- names(samples)[grepl("bsp", names(samples))]
bsp <- bsp[!grepl("sigma", bsp)]

jpeg("figures/bspParameterization_only_asp_asite_atreeid_bsp.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
util$plot_div_pairs(bsp, "sigma_bsp", samples, diagnostics, transforms = list("sigma_bsp" = 1))
dev.off()
