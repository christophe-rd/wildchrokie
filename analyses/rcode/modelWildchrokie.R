# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
quartz()

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)
library(RColorBrewer)
library(shinystan)
library(rethinking)
library("wesanderson")

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
sigma_y <- 0.2

n_perspp <- 10
n_spp <- 50
n_ids <- n_perspp * n_spp
n_meas <- 3
N <- n_ids * n_meas

# set spp names and ids 
spp <- 1:n_spp
ids <- rep(rep(1:n_perspp, times = n_spp), each = n_meas)
tree_spp <- rep(rep(1:n_spp, each = n_perspp), each = n_meas) # 1500 obs
tree_spp_nonrep <- rep(rep(1:n_spp, each = n_perspp)) # 500 obs

# === === === === #
# partial pooling
sigma_a_spp <- 0.3
sigma_a_ids <- 0.15

# get 50 intercept values for each species
a_spp_values <- rnorm(n_spp, 0, sigma_a_spp)

# combine species-level mean + id deviation at the non-repeated id level
a_ids_spp_values <- rnorm(n_ids, a_spp_values[tree_spp_nonrep], sigma_a_ids)
# a_ids_spp_values <- a_ids_spp_values[tree_spp_nonrep] - a_spp_values[tree_spp_nonrep]

a_ids_obs <- rep(a_ids_spp_values, each = n_meas)

# gdd and devide by constant
gdd <- round(rnorm(N, 1800, 100))
gddcons <- gdd / 200
# centering (optional)
# gddcons <- gddcons - mean(gddcons)

# error
error <- rnorm(N, 0, sigma_y)

# match a_ids to ids and spp (keeps your merge structure)
df1 <- data.frame(
  ids_nonrep  = rep(1:n_perspp, times = n_spp),
  spp = tree_spp_nonrep,
  a_spp = a_spp_values[tree_spp_nonrep],
  a_ids_spp_values = a_ids_spp_values
)

df2 <- data.frame(
  ids = ids,
  spp = tree_spp
)

simcoef <- merge(df2, df1, by.x = c("ids", "spp"), by.y = c("ids_nonrep", "spp"), all.x = TRUE)
simcoef <- simcoef[order(simcoef$spp, simcoef$ids), ]

simcoef$a <- a
simcoef$b <- b
simcoef$gddcons <- gddcons
simcoef$sigma_y <- sigma_y
simcoef$error <- error

# calculate ring width
simcoef$ringwidth <- simcoef$a + 
  simcoef$a_ids_spp_values + 
  # simcoef$a_spp + 
  (simcoef$b * simcoef$gddcons) + 
  error

# prepare grouping factors for stan_lmer (ensure ids are unique within spp)
simcoef$spp <- factor(simcoef$spp)
# simcoef$ids <- factor(paste0(simcoef$spp, "_", as.character(simcoef$ids)))


# loop through cols and round to 3 decimal points
for (i in 3:ncol(simcoef)) {
  simcoef[, i] <- round(simcoef[, i], 3)
}

# === === === === === === === === #
##### Run model #####
# === === === === === === === === #

###### Model nested on the intercept #######
fitnested <- TRUE
if(fitnested) {
  fitnested <- stan_lmer(
    ringwidth ~ 1 + gddcons + 
      (1|spp/ids), # this estimates a spp level intercept AND ids nested within spp
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
  # saveRDS(fitnested, "output/fitnested")
}
fitnested

###### Model partial pooled on b ######
runmodeWithPartialPooledBeta <- FALSE
if(runmodeWithPartialPooledBeta) {
  fit <- stan_lmer(
    ringwidth ~ 1 + gddcons + 
      (1 | ids) +
      (1 | spp) +
      (1|spp/ids),  #(1 | spp) means that I am partial pooling for the spp interecept
    # to confirm: to partial pool on both the intercept AND the slope : (1 + gddcons|spp)
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
}


print(fit, digits=3)

# === === === === === === #
##### Parameter recovery #####
# === === === === === === #

# === === === === === === #
###### Recover a_ids ######
# === === === === === === #
# Access the 'ids' data frame
fitnested_ranef <- ranef(fitnested)
ids_df <- fitnested_ranef$ids

# Now extract the tree IDs and intercepts
a_idswithranef <- data.frame(
  ids_spp = rownames(ids_df),
  fit_a_ids_spp = ids_df[["(Intercept)"]]
)
# clean ids col
a_idswithranef$ids <- sub("([^:]+):.*", "\\1", a_idswithranef$ids_spp)
a_idswithranef$spp <- sub(".*:(.*)", "\\1", a_idswithranef$ids_spp)

messyinter <- as.data.frame(posterior_interval(fitnested))
messyinter$messyids <- rownames(messyinter)
a_ids_spp_messyinter <- subset(messyinter, grepl("ids:spp", messyids))

# extract spp and ids names
a_ids_spp_messyinter$ids <- sub(".*ids:spp:([0-9]+):.*", "\\1", a_ids_spp_messyinter$messyids)
a_ids_spp_messyinter$spp <- sub(".*:([0-9]+)\\]$", "\\1", a_ids_spp_messyinter$messyids)

# remove non necessary columns
a_ids_spp_messyinter <- a_ids_spp_messyinter[, c("ids", "spp", "5%", "95%")]
# renames 5% and 95%
colnames(a_ids_spp_messyinter) <- c("ids", "spp", "per5", "per95")
# merge both df by ids
a_ids_mergedwithranef <- merge(a_idswithranef, a_ids_spp_messyinter, by = c("ids", "spp"))
# add simulation data and merge!
simcoeftoplot2 <- simcoef[, c("ids", "spp", "a_ids_spp_values", "a_spp")]
colnames(simcoeftoplot2) <- c("ids", "spp", "sim_a_ids_spp", "a_spp")

a_ids_mergedwithranef <- merge(simcoeftoplot2, a_ids_mergedwithranef, by = c("ids", "spp"))


# plot!
plot_a_ids_mergedwithranef_nested_a <- 
  ggplot(a_ids_mergedwithranef, aes(x = sim_a_ids_spp, y = fit_a_ids_spp)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.05) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
  theme_minimal()
plot_a_ids_mergedwithranef_nested_a

# === === === === === === #
###### Recover a_spp ######
# === === === === === === #
spp_df <- fitnested_ranef$spp

# Now extract the tree IDs and intercepts
a_sppwithranef <- data.frame(
  spp = rownames(spp_df),
  fit_a_spp = spp_df[["(Intercept)"]]
)
# recover only conf intervals spp from previously created df
a_spp_messyinter <- subset(messyinter, grepl("Intercept) spp:", messyids))
a_spp_messyinter$spp <- sub(".*spp:([0-9]+)]", "\\1", a_spp_messyinter$messyids)

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
a_spp_mergedwithranef2 <- merge(simcoeftoplot, a_spp_mergedwithranef, by = "spp")

# plot!
plot_a_spp_mergedwithranef_nested_a <- ggplot(a_spp_mergedwithranef2, aes(x = sim_a_spp, y = fit_a_spp)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_spp", y = "Model a_spp", title = "") +
  theme_minimal()
plot_a_spp_mergedwithranef_nested_a


# === === === === === === === === === === #
##### New way to recover parameters #####
# === === === === === === === === === === #
df_fit <- as.data.frame(fitnested)

# recover slope
colnames(df_fit)
# grab ids nested in spp
ids_cols <- colnames(df_fit)[grepl("ids:spp:", colnames(df_fit))]
ids_cols <- ids_cols[1:length(ids_cols)-1]
ids_df <- df_fit[, colnames(df_fit) %in% ids_cols]
# change their names
colnames(ids_df) <- sub(".*ids:spp:(.*)\\]$", "\\1", colnames(ids_df))
# empty ids dataframe
ids_df2 <- data.frame(
  ids_spp = character(ncol(ids_df)),
  fit_a_ids_spp = numeric(ncol(ids_df)),  
  fit_per5 = NA, 
  fit_per95 = NA,
  fit_sd = NA
)
for (i in 1:ncol(ids_df)) { # i = 1
  ids_df2$ids_spp[i] <- colnames(ids_df)[i]         
  ids_df2$fit_a_ids_spp[i] <- round(mean(ids_df[[i]]),3)  
  ids_df2$fit_per5[i] <- round(quantile(ids_df[[i]], probs = 0.055), 3)
  ids_df2$fit_per95[i] <- round(quantile(ids_df[[i]], probs = 0.945), 3)
  ids_df2$fit_sd[i] <- round(sd(ids_df[[i]]), 3)
}
ids_df2

# grab spp 
spp_cols <- colnames(df_fit)[grepl(" spp:", colnames(df_fit))]
spp_df <- df_fit[, colnames(df_fit) %in% spp_cols]
# change their names
colnames(spp_df) <- sub(".*spp:([0-9]+).*", "\\1", colnames(spp_df))
#empty spp df
spp_df2 <- data.frame(
  spp = character(ncol(spp_df)),
  fit_a_spp = numeric(ncol(spp_df)),  
  fit_per5 = NA, 
  fit_per95 = NA,
  fit_sd = NA
)
for (i in 1:ncol(spp_df)) { # i = 1
  spp_df2$spp[i] <- colnames(spp_df)[i]         
  spp_df2$fit_a_spp[i] <- round(mean(spp_df[[i]]),3)  
  spp_df2$fit_per5[i] <- round(quantile(spp_df[[i]], probs = 0.055), 3)
  spp_df2$fit_per95[i] <- round(quantile(spp_df[[i]], probs = 0.945), 3)
  spp_df2$fit_sd[i] <- round(sd(spp_df[[i]]), 3)
}
spp_df2

# === === === === === === === === === === === === === === === === 
# Plot old vs new way tp recover parameters #####
# === === === === === === === === === === === === === === === === 
colnames(a_spp_mergedwithranef) <- c("spp", "a_spp_lme4", "per5_lme4", "per95_lme4")
colnames(spp_df2) <- c("spp", "a_spp_loop", "per5_loop", "per95_loop")
recovComparison <- merge(a_spp_mergedwithranef, spp_df2, by = "spp")

recovComparison_plot <- ggplot(recovComparison, aes(x = a_spp_lme4, y = a_spp_loop)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(xmin = per5_lme4, xmax = per95_lme4), width = 0, color = "darkgray", alpha=0.5) +
  geom_errorbar(aes(ymin = per5_loop, ymax = per95_loop), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "a_spp lme4 functions", y = "a_spp manual posterior recovery", title = "") +
  theme_minimal()
# save ggplot!
ggsave("figures/recovComparison_plot.jpeg", recovComparison_plot, width = 6, height = 6, units = "in", dpi = 300)



simVSfit <- merge(simcoeftoplot, spp_df2, by = "spp")
simVSfit_plot <- ggplot(simVSfit, aes(x = sim_a_spp, y = a_spp_loop)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5_loop, ymax = per95_loop), width = 0, color = "darkgray", alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "a_spp sim", y = "a_spp fit", title = "") +
  theme_minimal()
# save ggplot!
simVSfit_plot
ggsave("figures/simVSfit_plot.jpeg", simVSfit_plot, width = 6, height = 6, units = "in", dpi = 300)

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


