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
sigma_a_ids <- 0.2

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

# calculate ring width
simcoef$ringwidth <- simcoef$a + 
  simcoef$a_ids_spp_values + 
  simcoef$a_spp + 
  (simcoef$b * simcoef$gddcons) + 
  error

# prepare grouping factors for stan_lmer (ensure ids are unique within spp)
simcoef$spp <- factor(simcoef$spp)
simcoef$ids <- factor(paste0(simcoef$spp, "_", as.character(simcoef$ids)))


# loop through cols and round to 3 decimal points
for (i in 3:ncol(simcoef)) {
  simcoef[, i] <- round(simcoef[, i], 3)
}

# === === === === === === === === #
##### Run model #####
# === === === === === === === === #
###### Model non nested ######
runmodelnonested <- FALSE # this is only for the intercepts, so I need to re-create old df
if(runmodelnonested) {
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
  
  fitnonested <- stan_lmer(
    ringwidth ~ 1 + gddcons + 
      (1 | ids) +
      (1 | spp),
      # (1|spp/ids),  #(1 | spp) means that I am partial pooling for the spp interecept
    # to confirm: to partial pool on both the intercept AND the slope : (1 + gddcons|spp)
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
}

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


# coef: medians are used for point estimates. the sum of the random and fixed effects coefficients for each explanatory variable for each level of each grouping factor.
# se:  The se function returns standard errors based on mad. See the Uncertainty estimates section in print.stanmvreg for more details.
# Uncertainty estimates (MAD_SD): The standard deviations reported (labeled MAD_SD in the print output) are computed from the same set of draws described above and are proportional to the median absolute deviation (mad) from the median. Compared to the raw posterior standard deviation, the MAD_SD will be more robust for long-tailed distributions. These are the same as the values returned by se.
# random efects (ranef): ceofficients that vary by group (see gelman and hill p.259)

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
# extract the "1_1" part
a_ids_spp_messyinter$ids <- sub(".*ids:spp:([0-9]+_[0-9]+):.*", "\\1", a_ids_spp_messyinter$messyids)
a_ids_spp_messyinter$spp <- sub(".*ids:spp:([0-9]+)_.*", "\\1", a_ids_spp_messyinter$messyids)

# remove non necessary columns
a_ids_spp_messyinter <- a_ids_spp_messyinter[, c("ids", "5%", "95%")]
# renames 5% and 95%
colnames(a_ids_spp_messyinter) <- c("ids", "per5", "per95")
# merge both df by ids
a_ids_mergedwithranef <- merge(a_idswithranef, a_ids_spp_messyinter, by = c("ids"))
# add simulation data and merge!
simcoeftoplot2 <- simcoef[, c("ids", "a_ids_spp_values")]
colnames(simcoeftoplot2) <- c("ids", "sim_a_ids_spp")
a_ids_mergedwithranef <- merge(simcoeftoplot2, a_ids_mergedwithranef, by = c("ids"))

# plot!
plot_a_ids_mergedwithranef_nested_a <- ggplot(a_ids_mergedwithranef, 
                                              aes(x = sim_a_ids_spp, y = fit_a_ids_spp)) +
  geom_point(color = "blue", size = 2) +
  geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
  theme_minimal()
plot_a_ids_mergedwithranef_nested_a
ggsave("figures/a_ids_mergedwithranef_nested_a.jpeg", plot_a_ids_mergedwithranef_nested_a, width = 8, height = 6)

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
ggsave("figures/plot_a_spp_mergedwithranef_nested_a.jpeg", plot_a_spp_mergedwithranef_nested_a, width = 8, height = 6)

# === === === === === === === === === === #
#### Add model intercept to sim data ####
# === === === === === === === === === === #
simcoef$asppfull <- simcoef$a + simcoef$a_spp
# add fitnonested data for that one
a_spp_mergedwithranef$asppfull_fit <- a_spp_mergedwithranef$fit_a_spp+fixef(fitnonested)[1]
# merge with simcoef
intercept_fit_sim <- merge(simcoef, a_spp_mergedwithranef[, c("spp", "asppfull_fit")], by = "spp")

# select 15 spp randomly out of the spp column
spp_to_plot <- sample(unique(intercept_fit_sim$spp), 50)
subtoplot2 <- subset(intercept_fit_sim, spp %in% spp_to_plot)



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

# Plot simulated data ####
# ------------------------------------------------------------------------------
# open jpeg device
jpeg(filename = "figures/a_ids_spp_values_sigbyspp.jpg", width = 8, height = 6, units = "in", res = 300)
plot(x_vals, y_vals, type = "l", lwd = 2, col = "black",
     xlab = "Value", ylab = "Density",
     main = "a_spp_values with sigma unique per spp")
abline(v = a_ids_spp_values, col = "red", lwd = 0.2)
abline(v = a_spp_values, col = "blue", lwd = 1)
dev.off()

# ------------------------------------------------------------------------------
jpeg(filename = "figures/a_ids_spp_values_uniquesig.jpg", width = 8, height = 6, units = "in", res = 300)
plot(x_vals, y_vals, type = "l", lwd = 2, col = "black",
     xlab = "Value", ylab = "Density",
     main = "a_spp_values with a unique sigma value")
abline(v = a_ids_spp_values, col = "red", lwd = 0.2)
abline(v = a_spp_values, col = "blue", lwd = 1)
dev.off()

simcoefvisal <- expand.grid(
  ids = unique(ids),
  spp = unique(tree_spp_num)
)

# attach your 500 values
simcoefvisal$a_ids_spp_values <- a_ids_spp_values
simcoefvisal$a_spp_values <- a_spp_values[simcoefvisal$spp]

simcoefvisal$spp <- as.factor(simcoefvisal$spp)

nestedbyspp_uniquesigma <- ggplot(simcoefvisal) +
  geom_vline(aes(xintercept = a_ids_spp_values), 
             color = "red", alpha = 0.3) +
  geom_vline(aes(xintercept = a_spp_values), 
             color = "blue", alpha = 0.8) +
  facet_wrap(~spp) +
  theme_minimal()
nestedbyspp_uniquesigma
# save!
ggsave("figures/nestedbyspp_uniquesigma.jpeg", nestedbyspp_uniquesigma, width = 8, height = 6)

# === === === === === === === === #
##### Check simulated data ######
# === === === === === === === === #
simcoef$intercepttemp <- simcoef$a_ids + simcoef$a_spp + simcoef$a
simcoef$asppfull <- simcoef$a + simcoef$a_spp

# select 15 spp randomly out of the spp column
spp_to_plot <- sample(unique(simcoef$spp), 16)
subtoplot <- subset(simcoef, spp %in% spp_to_plot)

# try to plot the values exclusively of a_ids_spp
ggplot(simcoef)+
  geom_point(aes(x=gddcons, y=a_ids, color = spp)) +
  facet_wrap(~spp)


# ring width X gdd cons by spp with intercept
ringXgddcons <- ggplot(subtoplot, aes(gddcons, ringwidth)) +
  geom_point(aes(color = "sim data")) +
  geom_abline(
    aes(intercept = intercepttemp, slope = 0.4, color = "sim a + a_spp"),
    data = subtoplot,
    linetype = "solid", alpha = 0.5
  ) +
  facet_wrap(~ spp) +
  
  scale_color_manual(
    name = "Legend",
    values = c(
      "sim data" = "black",
      "sim a + a_spp" = "red"
    )
  ) +
  theme_minimal()
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

# check species level intercept:
# ok so that's one option
diff_intercept_comparison <- ggplot(simcoef, aes(x = spp)) +
  geom_point(aes(y = intercepttemp, color = "a + a_spp + a_ids"), alpha = 0.7) +
  geom_point(aes(y = asppfull, color = "a + a_spp"), size = 2) +
  geom_hline(aes(yintercept = a, color = "a")) +
  scale_color_manual(
    name = "",
    values = c(
      "a + a_spp" = "black",
      "a" = "red",
      "a + a_spp + a_ids" = "darkgrey"
    )
  ) +
  labs(
    title = "",
    x = "spp",
    y = "values"
  ) +
  theme_minimal() +  # apply minimal theme first
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
diff_intercept_comparison 
ggsave("figures/diff_intercept_comparison.jpeg", diff_intercept_comparison, width = 10, height = 6)

# OLD CODE ####
# Non nested model #####
if(runmodelnonested) {
  fitnonested_ranef <- ranef(fitnonested)
  
  # === === === === === === #
  ###### Recover a_ids ######
  # === === === === === === #
  # Access the 'ids' data frame
  ids_df <- fitnonested_ranef$ids
  
  # Now extract the tree IDs and intercepts
  a_idswithranef <- data.frame(
    ids = rownames(ids_df),
    fit_a_ids = ids_df[["(Intercept)"]]
  )
  messyinter <- as.data.frame(posterior_interval(fitnonested))
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
  ggsave("figures/a_ids_mergedwithranef.jpeg", plot_a_ids_mergedwithranef, width = 8, height = 6)
  
  # === === === === === === #
  ###### Recover a_spp ######
  # === === === === === === #
  spp_df <- fitnonested_ranef$spp
  
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
  a_spp_mergedwithranef2 <- merge(simcoeftoplot, a_spp_mergedwithranef, by = "spp")
  
  # plot!
  plot_a_spp_mergedwithranef <- ggplot(a_spp_mergedwithranef2, aes(x = sim_a_spp, y = fit_a_spp)) +
    geom_point(color = "blue", size = 2) +
    geom_errorbar(aes(ymin = per5, ymax = per95), width = 0, color = "darkgray", alpha=0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
    labs(x = "Simulated a_ids", y = "Model a_ids", title = "") +
    theme_minimal()
  plot_a_spp_mergedwithranef
  ggsave("figures/a_spp_mergedwithranef.jpeg", plot_a_spp_mergedwithranef, width = 8, height = 6)
  
  # === === === === === === === === === === #
  #### Add model intercept to sim data ####
  # === === === === === === === === === === #
  simcoef$asppfull <- simcoef$a + simcoef$a_spp
  # add fitnonested data for that one
  a_spp_mergedwithranef$asppfull_fit <- a_spp_mergedwithranef$fit_a_spp+fixef(fitnonested)[1]
  # merge with simcoef
  intercept_fit_sim <- merge(simcoef, a_spp_mergedwithranef[, c("spp", "asppfull_fit")], by = "spp")
  
  # select 15 spp randomly out of the spp column
  spp_to_plot <- sample(unique(intercept_fit_sim$spp), 50)
  subtoplot2 <- subset(intercept_fit_sim, spp %in% spp_to_plot)
  
  # ring width X gdd cons by spp with intercept
  ringXgddcons2 <- ggplot(subtoplot2, aes(gddcons, ringwidth)) +
    # Add black dots with a legend
    geom_point(aes(color = "sim data")) +
    # Red line with custom label
    geom_abline(
      aes(intercept = asppfull, slope = 0.4, color = "sim a + a_spp"),
      data = subtoplot2,
      linetype = "solid", alpha = 0.5
    ) +
    # Blue line with custom label
    geom_abline(
      aes(intercept = asppfull_fit, slope = 0.4, color = "fitnonested a + a_spp"),
      data = subtoplot2,
      linetype = "solid", alpha = 0.5
    ) +
    scale_color_manual(
      name = "Legend",
      values = c(
        "sim data" = "black",
        "sim a + a_spp" = "red",
        "fitnonested a + a_spp" = "blue"
      )
    ) +
    facet_wrap(~ spp) +
    theme_minimal()
  # Show plot
  ringXgddcons2
  ggsave("figures/ringXgddcons_simANDfit2.jpeg", ringXgddcons2, width = 12, height = 8)
  
  
  # do the same graph, but only for a_spp to check if the overall intecept recovery is the problem.
  intercept_fit_sim2 <- merge(simcoef, a_spp_mergedwithranef[, c("spp", "fit_a_spp")], by = "spp")
  
  ringXgddcons3 <- ggplot(intercept_fit_sim2, aes(gddcons, ringwidthnooverall_a)) +
    geom_point(aes(color = "sim data")) +
    geom_abline(
      aes(intercept = a_spp, slope = 0.4, color = "sim a + a_spp"),
      data = intercept_fit_sim2,
      linetype = "solid", alpha = 0.5
    ) +
    geom_abline(
      aes(intercept = fit_a_spp, slope = 0.4, color = "fitnonested a + a_spp"),
      data = intercept_fit_sim2,
      linetype = "solid", alpha = 0.5
    ) +
    scale_color_manual(
      name = "Legend",
      values = c(
        "sim data" = "black",
        "sim a + a_spp" = "red",
        "fitnonested a + a_spp" = "blue"
      )
    ) +
    facet_wrap(~ spp) +
    theme_minimal()
  # Show plot
  ringXgddcons3
  ggsave("figures/ringXgddcons_simANDfit3.jpeg", ringXgddcons3, width = 12, height = 8)
  # I can visually see that 15 fit are > sim
  # vs 11 sim > fit
  # trying to understand why my model underestimates my overall interecept
}
