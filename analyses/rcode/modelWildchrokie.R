# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
# rm(list=ls()) 
# options(stringsAsFactors = FALSE)
options(max.print = 200) 

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
set.seed(123)
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

#figures
makeplot <- FALSE
if (makeplot) {
  # make long
  coef_long <- coefbind %>%
    select(ids, a_ids, a, coefsource) %>%
    pivot_longer(cols = c(a_ids, a), names_to = "type", values_to = "value")
  
  box <- ggplot(coef_long, aes(x = type, y = value, color = type)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1.5) +
    facet_wrap(~coefsource)+ 
    scale_color_manual(values = c("a_ids" = "#1b9e77",  
                                  "a" = "#d95f02")) +   
    labs(title = "sim vs model intercept",
         x = "tntercept type",
         y = "Value") +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave("figures/paramrecovery_box.jpeg", box, width = 8, height = 6)
  
}

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


# old code
if (runoldcode) {
  Nspp <- length(spp) # number of species
  Nplot <- length(plot) # number of plot per species
  Nyr <- length(yr) # number of growth years per tree
  Nrep <- length(rep) # number of measurements per tree
  # First making a data frame for the growth ring data
  Nw <- Nplot*Nyr*Nrep*Nspp # number of measurements per species
  
  wdat <- data.frame(
    yr = rep(yr, each = length(plot) * length(spp) * length(rep)),
    plot = rep(rep(plot, each = length(spp) * length(rep)), times = length(yr)),
    spp = rep(rep(spp, each = length(rep)), times = length(yr) * length(plot)),
    rep = rep(rep, times = length(yr) * length(plot) * length(spp))
  )
  
  # get posterior mean using ranef (random effects...). I guess it shows how much each tree's intercept deviates from the overall mean
  ranef_estimates <- ranef(fit)$ids
  ranef_estimates
  # overall intercept and slope 
  fixef(fit)
  
  # get a df with my 4 parameters estimated from the model
  param <- c("slope", "intercept", "sigma_y", "sigma_ids" )
  # slope
  fit_b <- fixef(fit)["gddcons"]
  fit_b_uncer <- 0.044
  # sigma_y
  fit_sigma_y <- 0.310
  fit_sigma_y_uncer <- 0.015
  # sigma_ids
  fit_sigma_ids <- 0.500
  fit_sigma_ids_uncer <- NA
  # a 
  fit_a <- 1.277
  fit_a_uncer <- 0.399
  
  # model ouput means
  fit_means <- c(fit_b, fit_a, fit_sigma_y, fit_sigma_ids)
  # model output sd
  fit_sd <- c(fit_b_uncer, fit_a_uncer, fit_sigma_y_uncer, fit_sigma_ids_uncer)
  
  fit_estimate <- data.frame(
    parameter= param,
    assumed_value = sim_param,
    estimate = fit_means,
    uncertainty = fit_sd
  )
  
  recovery <- ggplot(fit_estimate) +
    geom_point(aes(x = parameter, y = estimate, color = "Model estimate"), size = 3) +
    geom_errorbar(aes(x = parameter,
                      ymin = estimate - uncertainty,
                      ymax = estimate + uncertainty,
                      color = "Model estimate"),
                  width = 0,
                  na.rm = TRUE) +
    geom_point(aes(x = parameter, y = assumed_value, color = "Simulated"), size = 3, shape = 17) +
    labs(title = "",
         x = "Parameter",
         y = "Model estimate",
         color = "Type") +  # Legend title
    scale_color_manual(values = c("Model estimate" = "blue", "Simulated" = "orange")) +
    theme_minimal(base_size = 14)
  recovery
  ggsave("figures/model_parameter_recovery.jpeg", recovery, width = 8, height = 6)
  
  shinystan::launch_shinystan(fit)
  
  # verify how well I am returning my parameters:
  # simulated data parameters:
  mu.grand
  sim_spp_effects <- unique(wdat[, c("spp", "mu.spp")])
  sim_tree_effects <- unique(wdat[, c("treeid", "mu.tree")])
  sim_yr_effects <- unique(wdat[, c("yr", "mu.yr")])
  sim_intercept <- mu.grand + sim_spp_effects$mu.spp[sim_spp_effects$spp == "alninc"]
  sim_resid_sd <- sqrt(w.var)
  
  
  # pull model parameters
  # sim coef relative to alninc
  true_spp_coefs <- sim_spp_effects$mu.spp - sim_spp_effects$mu.spp[sim_spp_effects$spp == "alninc"]
  names(true_spp_coefs) <- paste0("spp", sim_spp_effects$spp)
  
  # Compare with model estimates
  data.frame(
    Parameter = c("(Intercept)", names(sim_spp_effects)[-1]), 
    Estimated = fit$coefficients[1:4],  # Assuming (Intercept), sppbetall, etc.
    True = c(sim_intercept, sim_spp_effects[-1])
  )
  
  
  #make a dataframe for ring width
  wdat <- data.frame(
    yr = rep(yr, each = length(plot) * length(spp) * length(rep)),
    plot = rep(rep(plot, each = length(spp) * length(rep)), times = length(yr)),
    spp = rep(rep(spp, each = length(rep)), times = length(yr) * length(plot)),
    rep = rep(rep, times = length(yr) * length(plot) * length(spp))
  )
  
  #
  
  # merge wdat and gdd
  wdat <- merge(wdat, gdd_sim, by.x = "yr", by.y = "year", all.x = TRUE)
  
  # paste plot and spp to give unique ids
  wdat$treeid <- paste(wdat$spp, wdat$plot, sep = "_")
  
  # the grand mean of ring width: baseline ring width when all other effects are zero.
  mu_grand <- 0.04 
  
  # Simulate random effect to set unique intercepts for each parameter
  wdat$yr_effect <- rnorm(nrow(wdat), 0, 0.05)[match(wdat$yr, yr)] # with mean of 0 and SD of 0.05
  wdat$tree_effect <- rnorm(nrow(wdat), 0, 0.02)[match(wdat$plot, plot)] # SD of 0.02 across trees
  wdat$spp_effect <- rnorm(nrow(wdat), 0, 0.06)[match(wdat$spp, spp)]
  
  # Simulate ring width (including GDD effect) aka yhat?
  wdat$ringwidth <- mu_grand + 
    wdat$yr_effect + 
    wdat$tree_effect + 
    wdat$spp_effect +
    spp_slopes[wdat$spp] * wdat$gdd +  # gdd effect varies by spp
    rnorm(nrow(wdat), 0, 0.5) #overall error
}
