# compare gdd at leafout model with and without bound
# crd 28 January 2026

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

# directories
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

# RUN MODEL ON EMPIRICAL DATA ####
# empirical data
emp <- read.csv("output/empiricalDataNORing.csv")

no_naleafout <- emp[!is.na(emp$leafoutGDD),]

# give numeric ids to my groups 
no_naleafout$site_num <- match(no_naleafout$site, unique(no_naleafout$site))
no_naleafout$spp_num <- match(no_naleafout$spp, unique(no_naleafout$spp))
no_naleafout$treeid_num <- match(no_naleafout$treeid, unique(no_naleafout$treeid))

# set a value to scale down the gdd
scale <- 20
y <- no_naleafout$leafoutGDD/scale
N <- nrow(no_naleafout) 
species <- as.numeric(as.character(no_naleafout$spp_num))
Nspp <- length(unique(no_naleafout$spp_num))
site <- as.numeric(as.character(no_naleafout$site_num))
Nsite <- length(unique(no_naleafout$site_num))
treeid <- as.numeric(no_naleafout$treeid_num)
Ntreeid <- length(unique(treeid))

fitNoBound <- stan("stan/modelGDDatLeafout.stan", 
            data=c("N","y",
                   "Nspp","species",
                   "Nsite","site",
                   "Ntreeid", "treeid"),
            iter=4000, chains=4, cores=4)

fitWithBound <- stan("stan/modelGDDatLeafout_bound.stan", 
                   data=c("N","y",
                          "Nspp","species",
                          "Nsite","site",
                          "Ntreeid", "treeid"),
                   iter=4000, chains=4, cores=4)


samplesNoBound <- util$extract_expectand_vals(fitNoBound)
samplesWithBound <- util$extract_expectand_vals(fitWithBound)

# back convert the object samples to the original gdd scale
y2 <- y*scale
samplesNoBound2 <- lapply(samplesNoBound, function(x) x*scale)
samplesWithBound2 <- lapply(samplesWithBound, function(x) x*scale)

yrepNamesNoBound <- grep("^y_rep\\[", names(samplesNoBound2), value = TRUE)
yrepNamesWithBound <- grep("^y_rep\\[", names(samplesWithBound2), value = TRUE)

yrepNoBound <- unlist(samplesNoBound2[yrepNamesNoBound])
yrepWithBound <- unlist(samplesWithBound2[yrepNamesWithBound])

nobound <- data.frame(
  model = "nobound",
  vals = yrepNoBound)

bound <- data.frame(
  model = "withbound",
  vals = yrepWithBound)

binded <- rbind(nobound, bound)

ggplot(binded, aes(x = vals, color = factor(model))) +
  geom_histogram(
    binwidth = 10,
    position = "identity",
    alpha = 0,
    linewidth = 1
  ) +
  geom_vline(xintercept = 100) +
  scale_color_manual(values = wes_palette("AsteroidCity1")[c(1,3)]) +
  # scale_x_continuous(breaks = seq(0, 700, by = 100)) +
  xlim(-50, 400) +
  labs(color = "Model with bound at 120 or without bound", title = "Generated quantities") +
  theme_minimal()
ggsave("figures/troubleShootingLeafoutGDDmodel/genQuantitiesBoundVSnoBound.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# just the posterior of a
fitNoBound_df <- as.data.frame(fitNoBound)
fitWithBound_df <- as.data.frame(fitWithBound)

aNoBound <- data.frame(
  a = fitNoBound_df["a"],
  model = "nobound")

aWithBound <- data.frame(
  a = fitWithBound_df["a"],
  model = "withbound")

abinded <- rbind(aNoBound, aWithBound)
abinded$a <- abinded$a* scale


ggplot(abinded, aes(x = a, color = factor(model))) +
  geom_histogram(
    binwidth = 10,
    position = "identity",
    alpha = 0,
    linewidth = 1
  ) +
  geom_vline(xintercept = 100) +
  scale_color_manual(values = wes_palette("AsteroidCity1")[c(1,3)]) +
  xlim(-50, 400) +
  # scale_x_continuous(breaks = seq(-50, 700, by = 20)) +
  labs(color = "Model with bound at 120 or without bound", title = "Posterior distributions of a") +
  theme_minimal()
ggsave("figures/troubleShootingLeafoutGDDmodel/post_a_BoundVSnoBound.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
samplesNoBound <- util$extract_expectand_vals(fitNoBound)
samplesWithBound <- util$extract_expectand_vals(fitWithBound)

util$plot_hist_quantiles(samplesNoBound2, "y_rep", 
                         -50, # lower x axis limit
                         700, # upper x axis limit
                         20, # binning
                         baseline_values = y2,
                         xlab = "gdd at leafout")

util$plot_hist_quantiles(samplesWithBound2, "y_rep", 
                         -50, # lower x axis limit
                         700, # upper x axis limit
                         20, # binning
                         baseline_values = y2,
                         xlab = "gdd at leafout")
