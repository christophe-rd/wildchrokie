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

n_perspp <- 25
n_spp <- 20
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
plot(ringwidth~gddcons, data=simcoef)

# run models
runmodels <- TRUE
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

#### Trying different functions for parameter recovery ####
posterior_interval(fit) # nice!
coefficients(fit)
fitted.values(fit)
summary(fit, pars = "alpha")

plot(fit)
as.data.frame(fit)

ranef(fit)$spp[, "(Intercept)"] 

# coef: medians are used for point estimates. the sum of the random and fixed effects coefficients for each explanatory variable for each level of each grouping factor.
# se:  The se function returns standard errors based on mad. See the Uncertainty estimates section in print.stanmvreg for more details.
# Uncertainty estimates (MAD_SD): The standard deviations reported (labeled MAD_SD in the print output) are computed from the same set of draws described above and are proportional to the median absolute deviation (mad) from the median. Compared to the raw posterior standard deviation, the MAD_SD will be more robust for long-tailed distributions. These are the same as the values returned by se.
# === === === === === === === === === === === === === === === === === === === === == 
#### Recover model coef using the coef function probably recovering the median #####
# === === === === === === === === === === === === === === === === === === === === == 

# === === === === === === #
##### start with IDS #####
# === === === === === === #
fitcoef <- as.data.frame(coef(fit)$ids)
fitcoef$a_ids <- fitcoef$`(Intercept)`- fixef(fit)[1] 
# fitcoef$a_spp <- fitcoef

# remove column that is overall intercept + a_ids
fitcoef <- fitcoef[,-1]
fitcoef$ids <- row.names(fitcoef)
fitcoef$a <- fixef(fit)[1]
fitcoef$b <-  fixef(fit)[2]
posterior <- as.data.frame(fit)
sigma_draws <- posterior$sigma
fitcoef$sigma_y <- mean(sigma_draws)

# rename and reorganize!
# add species colum by grep ids
fitcoef$spp <- sub("\\d+", "", fitcoef$ids)  # extract species from ids

# now work to grab intercept estimates for spp!
sppoutput <- as.data.frame(coef(fit)$spp)
sppoutput$a_spp <- sppoutput$`(Intercept)`-fixef(fit)
sppoutput$spp <- rownames(sppoutput)

fitcoef$a_spp <- sppoutput[match(fitcoef$spp, sppoutput$spp), "a_spp"]

colnames(fitcoef) <- c("gddcons", "a_ids", "ids", "a", "b", "sigma_y", "spp", "a_spp")
fitcoef <- fitcoef[, c("ids", "spp", "b", "a", "a_ids", "a_spp", "sigma_y")]
fitcoef$coefsource <- "model"

# prep sim for merge!
simcoef <- simcoef[, c("ids", "spp", "b", "a", "a_ids", "a_spp", "sigma_y")]
simcoef$coefsource <- "simulated"

# reorganize to make a xy plot
simcoef[!duplicated(simcoef$a_ids), c(1,5)]
fitcoef[,c(1,4)]
merged_a_idswithcoef <- merge(simcoef[!duplicated(simcoef$a_ids), c(1,5)], fitcoef[,c(1,5)], by="ids")
colnames(merged_a_idswithcoef) <- c("ids", "sim_a_ids", "fit_a_ids")

plot_merged_a_idswithcoef <- ggplot(merged_a_idswithcoef, aes(x = fit_a_ids, y = sim_a_ids)) +
  geom_point(color = "blue", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
  theme_minimal()
plot_merged_a_idswithcoef
ggsave("figures/a_ids_with_coef.jpeg", plot_merged_a_idswithcoef, width = 8, height = 6)

# === === === === === === #
##### now for spp #####
# === === === === === === #
a_sppwithcoef <- as.data.frame(coef(fit)$spp)
a_sppwithcoef$a_spp_subtracted <- a_sppwithcoef$`(Intercept)`-fixef(fit)[1]
a_sppwithcoef$spp <- row.names(a_sppwithcoef)
# bring simcoef and remove duplicated cols
sppsimcoef <- simcoef[!duplicated(simcoef$spp), c(2,6)]
colnames(sppsimcoef) <- c("spp", "sim_a_spp")
merged_a_sppwithcoef <- merge(sppsimcoef[, c("spp", "sim_a_spp")], a_sppwithcoef[, c("spp", "a_spp_subtracted")], by = "spp")

plot_merged_a_sppwithcoef <- ggplot(merged_a_sppwithcoef, aes(x = sim_a_spp, y = a_spp_subtracted)) +
  geom_point(color = "blue", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_spp", y = "Model a_spp", title = "") +
  theme_minimal()
ggsave("figures/a_spp_with_coef.jpeg", plot_merged_a_sppwithcoef, width = 8, height = 6)
plot_merged_a_sppwithcoef

# === === === === === === === === === === === === === === === === === === === === == 
#### Try with the summary function -- using mean ####
# === === === === === === === === === === === === === === === === === === === === == 
fitsum <- as.data.frame(summary(fit)) 
fitsum$interceptType <- rownames(fitsum)
rownames(fitsum) <- NULL

# === === === === === === #
##### Recover a_ids #####
# === === === === === === #
# Subset for ids by using grep from the $interceptType
a_ids <- subset(fitsum, grepl("ids", interceptType))
a_ids$ids <- sub(".*ids:([^]]+)]", "\\1", a_ids$interceptType)
a_ids <- a_ids[, c(ncol(a_ids), 1:c(ncol(a_ids)-2))]
a_ids$spp <- sub("\\d+", "", a_ids$ids)  

## select columns that wil be in the plot
a_idstoplot <- a_ids[,c("ids", "mean", "10%", "90%")]
colnames(a_idstoplot) <- c("ids", "a_ids", "per10", "per90")
simcoeftoplot <- simcoef[, c("ids", "a_ids")]
simcoeftoplot <- simcoeftoplot[!duplicated(simcoeftoplot),]
a_ids_mergedwithsum <- merge(simcoeftoplot, a_idstoplot, by = "ids")
colnames(a_ids_mergedwithsum) <- c("ids", "sim_a_ids", "fit_a_ids", "per10", "per90")

plot_a_ids_mergedwithsum <- ggplot(a_ids_mergedwithsum, aes(x = sim_a_ids, y = fit_a_ids)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per10, ymax = per90), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
  theme_minimal()
plot_a_ids_mergedwithsum
ggsave("figures/a_ids_mergedwithsum.jpeg", plot_a_ids_mergedwithsum, width = 8, height = 6)

# === === === === === === #
##### Recover a_spp #####
# === === === === === === #
# Subset for spp by using grep from the $interceptType
a_spp <- subset(fitsum, grepl("spp:", interceptType))
a_spp$spp <- sub(".*spp:(spp[0-9]+)]", "\\1", a_spp$interceptType)
a_spp <- a_spp[, c(ncol(a_spp), 1:c(ncol(a_spp)-2))]

## select columns that wil be in the plot
a_spptoplot <- a_spp[1:20,c("spp", "mean", "10%", "90%")] # for now I don't select the sigma row
colnames(a_spptoplot) <- c("spp", "fit_a_spp", "fit_per10", "fit_per90%")
simcoeftoplot <- simcoef[, c("ids", "spp", "a_spp")]
colnames(simcoeftoplot) <- c("ids", "spp", "sim_a_spp")
simcoeftoplot <- simcoeftoplot[!duplicated(simcoeftoplot$ids),]

#make copy
sppcoefwithsumm <- simcoeftoplot

sppcoefwithsumm$fit_a_spp <- a_spptoplot$fit_a_spp[match(sppcoefwithsumm$spp, a_spptoplot$spp)]
sppcoefwithsumm$fit_per10 <- a_spptoplot$fit_per10[match(sppcoefwithsumm$spp, a_spptoplot$spp)]
sppcoefwithsumm$fit_per90 <- a_spptoplot$fit_per90[match(sppcoefwithsumm$spp, a_spptoplot$spp)]

sppcoefwithsumm <- sppcoefwithsumm[!duplicated(sppcoefwithsumm$spp),]


plot_sppcoefwithsumm <- ggplot(sppcoefwithsumm, aes(x = sim_a_spp, y = fit_a_spp)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = fit_per10, ymax = fit_per90), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_spp", y = "Model a_spp", title = "") +
  theme_minimal()
plot_sppcoefwithsumm
ggsave("figures/a_spp_mergedwithsum.jpeg", plot_sppcoefwithsumm, width = 8, height = 6)

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
  geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.2) +
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

# === === === === === === #
#### Compare methods ####
# === === === === === === #
a_spp_mergedwithranef
sppcoefwithsumm
# merge by ids, but keep only both df fit_a_spp
a_spp_mergedwithcoef <- merge(sppcoefwithsumm[, c("spp", "fit_a_spp")], a_spp_mergedwithranef[, c("spp", "fit_a_spp")], by = "spp")

plottocompare <- ggplot(a_spp_mergedwithcoef, aes(x = fit_a_spp.x, y = fit_a_spp.y)) +
  geom_point(color = "blue", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) 
plottocompare
# save plot to compare summary vs ranef functions
ggsave("figures/compare_spp_coef.jpeg", plottocompare, width = 8, height = 6)

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

