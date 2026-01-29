# Wildchrokie model
# CRD 29 January 2026

# larger priors on atreeid

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
library(rstanarm)

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

emp <- read.csv("output/empiricalDataMAIN.csv")

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
gdd <- emp$pgsGDD/200
Nspp <- length(unique(emp$spp_num))
Nsite <- length(unique(emp$site_num))
site <- as.numeric(as.character(emp$site_num))
species <- as.numeric(as.character(emp$spp_num))
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))

# check that everything is fine
table(treeid,species)

rstan_options(auto_write = TRUE)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Regular prior on sigma_atreeid
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fit <- stan("stan/twolevelhierint.stan", 
            data=c("N","y",
                   "Nspp","species",
                   "Nsite", "site", 
                   "Ntreeid", "treeid", 
                   "gdd"),
            iter=4000, chains=4, cores=4)
df_fit <- as.data.frame(fit)
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Larger prior on sigma_atreeid
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fit_LP <- stan("stan/twolevelhierint_largerPriors.stan", 
            data=c("N","y",
                   "Nspp","species",
                   "Nsite", "site", 
                   "Ntreeid", "treeid", 
                   "gdd"),
            iter=4000, chains=4, cores=4)

df_fit_LP <- as.data.frame(fit_LP)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]
sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]

sigma_cols_LP <- colnames(df_fit_LP)[grepl("sigma", colnames(df_fit_LP))]
sigma_df_LP <- df_fit_LP[, colnames(df_fit_LP) %in% sigma_cols_LP]

sigma_long <- reshape(
  sigma_df,
  direction = "long",
  varying = list(names(sigma_df)),
  v.names = "value",
  timevar = "parameter",
  times = names(sigma_df),
  idvar = "draw"
)

sigma_long_LP <- reshape(
  sigma_df_LP,
  direction = "long",
  varying = list(names(sigma_df_LP)),
  v.names = "value",
  timevar = "parameter",
  times = names(sigma_df_LP),
  idvar = "draw"
)

sigma_long$model <- "normal prior"
sigma_long_LP$model <- "3x prior"

sigma_long$LP <- sigma_long_LP$value
# binded <- rbind(sigma_long, sigma_long_LP)

# remove sigma y
sigma_long2 <- subset(sigma_long, parameter == "sigma_atreeid")

# sigma_long$prior <- NA
# sigma_long$prior[which(sigma_long$parameter == "sigma_atreeid")] <- rnorm(8e3, 0, 0.5)
# sigma_long$prior[which(sigma_long$parameter == "sigma_y")] <- rnorm(8e3, 0, 3)

postComparisonLargerPriors <- ggplot(sigma_long2) +
  geom_density(aes(x = value, colour = "Posterior with prior N(0, 0.5)"),
               linewidth = 0.8) +
  geom_density(aes(x = LP, colour = "Posterior with prior N(0, 1.5)"),
               linewidth = 0.8) +
  labs(title = "Posterior comparison with larger priors on sigma_atreeid",
       x = "", y = "Density", color = "Curve") +
  scale_color_manual(values = c("#C52E19", "#E58601")) +
  xlim(-2, 2) +
  theme_minimal()
postComparisonLargerPriors
ggsave("figures/troubleShootingGrowthModel/postComparisonLargerPriors_atreeid.jpeg", postComparisonLargerPriors, width = 8, height = 6, units = "in", dpi = 300)
