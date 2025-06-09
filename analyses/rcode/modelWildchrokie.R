# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
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
a <- 1.5
b <- 0.4
sigma_y <- 0.3 # standard deviation

n_ids <- 100
rep <- 3
N <- n_ids * rep

# partial pool for each tree id
sigma_ids <- 1.5/2.57 # I specify how much the data varies across individuals. 1.96, j'aurais 95% de mes valeurs qui seraient entre 1.5 et -1.5. Quantile de la loi normale
a_ids <- rnorm(n_ids, 0, sigma_ids)

# set ids
ids <- rep(paste0("t", 1:n_ids), each = rep)

# just trees
trees <- paste0("t", 1:n_ids)

# growing degree days
gdd <- round(rnorm(N, 1800, 100))

# divide by a constant to adjust the scale
gddcons <- gdd / 200

# overall error
error <- rnorm(N, 0, sigma_y) 

# match a_ids to tree so its the same for each ids replicate
a_ids <- a_ids[match(ids, trees)] 

# calculate ring width
ringwidth <- a + a_ids +  b * gddcons + error

# create df
simcoef <- data.frame(
  ids = ids,
  gddcons = gddcons,
  b = b,
  a = a,
  a_ids = a_ids,
  sigma_y = sigma_y,
  ringwidth = ringwidth
)
plot(ringwidth~gddcons, data=simcoef)

# keep parameters in the following vector
sim_param <- c(b,a, sigma_y, sigma_ids)

# run models
if(runmodels){
  fit <- stan_lmer(
    ringwidth ~ gddcons + (1 | ids),  
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
}
print(fit, digits=3)
summary(fit)
# recover model coefs
fitcoef <- as.data.frame(coef(fit)$ids)
fitcoef$a_ids <- fitcoef$`(Intercept)`- fixef(fit)[1]
# remove column that is overall intercept + a_ids
fitcoef <- fitcoef[,-1]
fitcoef$ids <- row.names(fitcoef)
fitcoef$a <- fixef(fit)[1]
fitcoef$b <-  fixef(fit)[2]
posterior <- as.data.frame(fit)
sigma_draws <- posterior$sigma
fitcoef$sigma_y <- mean(sigma_draws)

# rename and reorganize!
colnames(fitcoef) <- c("gddcons", "a_ids", "ids", "a", "b", "sigma_y")
fitcoef <- fitcoef[, c("ids", "b", "a", "a_ids", "sigma_y" )]
fitcoef$coefsource <- "model"

# prep sim for merge!
simcoef <- simcoef[, c("ids", "b", "a", "a_ids", "sigma_y")]
simcoef$coefsource <- "simulated"

# bind by row
coefbind <- rbind(fitcoef, simcoef)

#double check length of a_ids
length(unique(fitcoef$a_ids))
length(unique(simcoef$a_ids))

# reorganize to make a xy plot
simcoef[!duplicated(simcoef$a_ids), c(1,4)]
fitcoef[,c(1,4)]
merged <- merge(simcoef[!duplicated(simcoef$a_ids), c(1,4)], fitcoef[,c(1,4)], by="ids")
colnames(merged) <- c("ids", "sim_a_ids", "fit_a_ids")

# xy plot with 0,1 abline
jpeg("figures/xyplot_intercept.jpeg", width = 1600, height = 1200,  quality = 95,res = 150)
# Create the plot
plot(merged$fit_a_ids, merged$sim_a_ids,
     xlab = "Model a_ids", ylab = "Simulated a_ids",
     main = "Model vs Simulated a_ids", pch = 19, col = "blue")

# Add 1:1 reference line
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Close the device to save the file
dev.off()

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
