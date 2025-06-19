# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)
library(RColorBrewer)
runmodels <- FALSE
runoldcode <- FALSE

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

# === === === === === === === === === === === === === === === === 
#### Step 1. Come up with a model ####
# === === === === === === === === === === === === === === === === 
# As discussed with Victor, we will use a linear model to predict growth from GDD. 

# === === === === === === === === === === === === === === === === 
#### Step 2. Simulate data ####
# === === === === === === === === === === === === === === === === 
# set parameters
set.seed(124)
a <- 1.5
b <- 0.4
sigma_y <- 0.3

n_perspp <- 10
n_spp <- 50
n_ids <- n_perspp * n_spp
rep <- 3
N <- n_ids * rep

# set spp names and ids 
spp <- paste0(rep("spp", each = n_spp), 1:n_spp)
tree_ids <- paste0(rep(spp, each = n_perspp), "_", 1:n_perspp)  # 100 unique tree IDs
ids <- rep(tree_ids, each = rep)  # repeat each tree ID 3 times
tree_spp <- rep(rep(spp, each = n_perspp), each = rep)  # matching species for each ID

# partial pooling
sigma_ids <- 0.8 / 2.57
sigma_spp <- 1 / 2.57
a_ids_values <- rnorm(n_ids, 0, sigma_ids)
a_spp_values <- rnorm(n_spp, 0, sigma_spp)

# match to observations
a_ids <- rep(a_ids_values, each = rep)  # repeat each a_ids 3 times
a_spp <- rep(a_spp_values, each = n_perspp * rep)  # repeat each a_spp 75 times

# gdd and devide by constant
gdd <- round(rnorm(N, 1800, 100))
gddcons <- gdd / 200

# error
error <- rnorm(N, 0, sigma_y)

# calculate ring width
ringwidth <- a + a_ids + a_spp + b * gddcons + error
# set df
simcoef <- data.frame(
  ids = ids,
  spp = tree_spp,
  gddcons = gddcons,
  b = b,
  a = a,
  a_ids = a_ids,
  a_spp = a_spp,
  sigma_y = sigma_y,
  ringwidth = ringwidth
)
# === === === === === === === === #
##### Check simulated data ######
# === === === === === === === === #

# select 15 spp randomly out of the spp column
spp_to_plot <- sample(unique(simcoef$spp), 16)

subtoplot <- subset(simcoef, spp %in% spp_to_plot)
subtoplot$intercepttemp <- subtoplot$a_ids+subtoplot$a_spp+a

# ring width X gdd cons by spp with intercept
ringXgddcons <- ggplot(subtoplot, aes(gddcons, ringwidth)) +
  geom_point() +
  geom_abline(
    aes(intercept = intercepttemp, slope = 0.4),
    data = subtoplot,
    color = "red", linetype = "solid", alpha =0.5
  ) +
  facet_wrap(~ spp) +
  theme_minimal()
# save gg plot
ringXgddcons
ggsave("figures/ringXgddcons.jpeg", ringXgddcons, width = 8, height = 6)

#a_ids X a_spp
a_idsXa_spp <- ggplot(subtoplot) +
  geom_point(aes(a_spp, a_ids)) +
  geom_abline(
    aes(intercept = a_ids, slope = 0.4),
    data = subtoplot,
    color = "red", linetype = "dashed", alpha =0.5
  ) +
  geom_abline(aes(intercept = a_spp, slope = 0.4),
              data = subtoplot,
              color = "blue", linetype = "solid", alpha = 0.5
              )+
  facet_wrap(~spp)+
  # Optional: plot actual points if available 
  # labs(y = "", x = "") +
  theme_minimal()
a_idsXa_spp
# save 
ggsave("figures/a_idsXa_spp.jpeg", a_idsXa_spp, width = 8, height = 6)

# run models
runmodels <- FALSE
if(runmodels) {
  fit <- stan_lmer(
    ringwidth ~ gddcons + (1 | ids) + (1 | spp),  
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
}

print(fit, digits=3)


# coef: medians are used for point estimates. the sum of the random and fixed effects coefficients for each explanatory variable for each level of each grouping factor.
# se:  The se function returns standard errors based on mad. See the Uncertainty estimates section in print.stanmvreg for more details.
# Uncertainty estimates (MAD_SD): The standard deviations reported (labeled MAD_SD in the print output) are computed from the same set of draws described above and are proportional to the median absolute deviation (mad) from the median. Compared to the raw posterior standard deviation, the MAD_SD will be more robust for long-tailed distributions. These are the same as the values returned by se.

# === === === === === === #
#### Using ranef and fixef ####
# === === === === === === #
fitef <- ranef(fit)
# === === === === === === #
##### Recover a_ids #####
# === === === === === === #
# Access the 'ids' data frame
ids_df <- fitef$ids

# Now extract the tree IDs and intercepts
a_idswithranef <- data.frame(
  ids = rownames(ids_df),
  fit_a_ids = ids_df[["(Intercept)"]]
)
messyinter <- as.data.frame(posterior_interval(fit))
messyinter$messyids <- rownames(messyinter)
a_ids_messyinter <- subset(messyinter, grepl("ids", messyids))
a_ids_messyinter$ids <- sub(".*ids:([^]]+)]", "\\1", a_ids_messyinter$messyids)
# remove non necessary columns
a_ids_messyinter <- a_ids_messyinter[, c("ids", "5%", "95%")]
# renames 5% and 95%
colnames(a_ids_messyinter) <- c("ids", "per5", "per95")
# merge both df by ids
a_ids_mergedwithranef <- merge(a_idswithranef, a_ids_messyinter, by = "ids")
# add simulation data and merge!
simcoeftoplot2 <- simcoef[, c("ids", "a_ids")]
colnames(simcoeftoplot2) <- c("ids", "sim_a_ids")
a_ids_mergedwithranef <- merge(simcoeftoplot2, a_ids_mergedwithranef, by = "ids")

plot_a_ids_mergedwithranef<- ggplot(a_ids_mergedwithranef, aes(x = sim_a_ids, y = fit_a_ids)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
  theme_minimal()
plot_a_ids_mergedwithranef
ggsave("figures/a_ids_mergedwithranef.jpeg", plot_a_ids_mergedwithsum, width = 8, height = 6)

# === === === === === === #
##### Recover a_spp #####
# === === === === === === #
spp_df <- fitef$spp

# Now extract the tree IDs and intercepts
a_sppwithranef <- data.frame(
  spp = rownames(spp_df),
  fit_a_spp = spp_df[["(Intercept)"]]
)
# recover only conf intervals spp from previously created df
a_spp_messyinter <- subset(messyinter, grepl("spp:", messyids))

a_spp_messyinter$spp <- sub(".*spp:(spp[0-9]+)]", "\\1", a_spp_messyinter$messyids)
# remove unecessary columns
a_spp_messyinter <- a_spp_messyinter[, c("spp", "5%", "95%")]
# change 5 and 95% names
colnames(a_spp_messyinter) <- c("spp", "per5", "per95")
# merge!
a_spp_mergedwithranef <- merge(a_sppwithranef, a_spp_messyinter, by = "spp")
# get sim data ready to merge
simcoeftoplot <- simcoef[, c("spp", "a_spp")]
colnames(simcoeftoplot) <- c("spp", "sim_a_spp")
simcoeftoplot <- simcoeftoplot[!duplicated(simcoeftoplot),]
# now merge
a_spp_mergedwithranef <- merge(simcoeftoplot, a_spp_mergedwithranef, by = "spp")


plot_a_spp_mergedwithranef <- ggplot(a_spp_mergedwithranef, aes(x = sim_a_spp, y = fit_a_spp)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
  theme_minimal()
plot_a_spp_mergedwithranef
ggsave("figures/a_spp_mergedwithranef.jpeg", plot_a_spp_mergedwithranef, width = 8, height = 6)

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

