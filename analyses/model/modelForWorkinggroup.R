# 1 September
# By Christophe RD
# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 

# Load library 
library(rstanarm)
library(ggplot2)
library(wesanderson)

runmodels <- FALSE
runoldcode <- FALSE

# === === === === ===#
# Simulate data ####
# === === === === ===#

# set parameters
set.seed(124)
a <- 1.5
b <- 0.4
sigma_y <- 0.2
sigma_a_spp <- 0.3
sigma_a_ids <- 0.15

n_spp <- 50 # number of species
n_perspp <- 10 # number of individuals per species
n_ids <- n_perspp * n_spp # number of ids
n_meas <- 15 # repeated measurements per id
N <- n_ids * n_meas # total number of measurements

# set spp names and ids 
spp <- 1:n_spp
ids <- rep(rep(1:n_perspp, times = n_spp), each = n_meas)
tree_spp <- rep(rep(1:n_spp, each = n_perspp), each = n_meas) # 1500 obs
tree_spp_nonrep <- rep(rep(1:n_spp, each = n_perspp)) # 500 obs

# get 50 intercept values for each species
a_spp_values <- rnorm(n_spp, 0, sigma_a_spp)

# get intercept values for ids nested in spp
a_ids_spp_values <- rnorm(n_ids, a_spp_values[tree_spp_nonrep], sigma_a_ids)

# gdd and devide by constant
gdd <- round(rnorm(N, 1800, 100))
gddcons <- gdd / 200

# error
error <- rnorm(N, 0, sigma_y)

# get index of each measurement
id_index <- rep(1:n_ids, each = n_meas)

# build coefficient df
simcoef <- data.frame(
  ids  = ids,                    
  spp  = tree_spp,               
  a_spp = a_spp_values[tree_spp],
  a_ids_spp = a_ids_spp_values[id_index]        
)

# order it
simcoef <- simcoef[order(simcoef$spp, simcoef$ids), ]

# add other general coefficients
simcoef$a <- a
simcoef$b <- b
simcoef$gddcons <- gddcons
simcoef$sigma_y <- sigma_y
simcoef$error <- error

# --- --- --- --- --- --- --- --- --- --- --- ---
# calculate ring width option 1
simcoef$ringwidth <- simcoef$a + 
  simcoef$a_ids_spp + 
  # simcoef$a_spp + 
  (simcoef$b * simcoef$gddcons) + 
  error
# --- --- --- --- --- --- --- --- --- --- --- ---

# --- --- --- --- --- --- --- --- --- --- --- ---
# calculate ring width option 2
simcoef$ringwidth2 <- simcoef$a + 
  simcoef$a_ids_spp + 
  simcoef$a_spp +
  (simcoef$b * simcoef$gddcons) + 
  error
# --- --- --- --- --- --- --- --- --- --- --- ---

# round to 3 decimal points
for (i in 3:ncol(simcoef)) {
  simcoef[, i] <- round(simcoef[, i], 3)
}

# === === === === === === #
# Plot simulated data ####
# === === === === === === #
simcoef$fullintercept <- simcoef$a + 
  simcoef$a_ids_spp +
  simcoef$a_spp +
  simcoef$error

simcoef$fullintercept1 <- simcoef$a + 
  simcoef$a_ids_spp 

simcoef$fullintercept2 <- simcoef$a + 
  simcoef$a_ids_spp + 
  simcoef$a_spp 

# --- --- --- --- --- --- --- --- --- --- --- ---
fit1 <- stan_lmer(
    ringwidth ~ 1 + gddcons + (1|spp/ids),
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
)
fit2 <- stan_lmer(
    fullintercept ~ (1|spp/ids),
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
)

# === === === === === #
# Recover fit 2 ####
# === === === === === #
fitnested_ranef <- ranef(fit2)

ids_df <- fitnested_ranef$ids

a_idswithranef <- data.frame(
  ids = rownames(ids_df),
  fit_a_ids = ids_df[["(Intercept)"]]
)
# clean ids col
a_idswithranef$spp <- sub(".*:([0-9]+).*", "\\1", a_idswithranef$ids)
a_idswithranef$ids <- sub("^([0-9]+):.*", "\\1", a_idswithranef$ids)

# recover a_spp
spp_df <- fitnested_ranef$spp
# Now extract the tree IDs and intercepts
a_sppwithranef <- data.frame(
  spp = rownames(spp_df),
  fit_a_spp = spp_df[["(Intercept)"]]
)

# merge both dfs
a_merged <- merge(a_idswithranef, a_sppwithranef, by = c( "spp"))

# overall intercept 
a_merged$fit_a <- fit2$coefficients["(Intercept)"]

a_merged$fit_a_total <- a_merged$fit_a + a_merged$fit_a_spp + a_merged$fit_a_ids

# merge with sim data
a_recovermerged <- merge(a_merged, simcoef, by = c("spp", "ids"))

a_recovermerged$ids <- as.numeric(as.character(a_recovermerged$ids))
a_recovermerged$spp <- as.numeric(as.character(a_recovermerged$spp))

subforplot <- subset(a_recovermerged, spp %in% 
                       sample(unique(a_recovermerged$spp), 5) )

# Plot summed intercepts
histgridsim <- ggplot(subforplot, aes(x = fullintercept, fill = factor(spp))) +
  geom_histogram(bins = 30) +
  geom_vline(aes(xintercept = fullintercept2), linetype = "solid", size = 0.5) + 
  geom_vline(aes(xintercept = fit_a_total), linetype = "dashed", size = 0.5) + 
  facet_grid(ids ~ spp) +
  scale_fill_manual(values = wes_palettes$Moonrise3) +
  theme_minimal()
histgridsim

# === === === === === #
# Recover fit 1 ####
# === === === === === #
# Access the 'ids' data frame
fitnested_ranef <- ranef(fit1)
ids_df <- fitnested_ranef$ids

# Now extract the tree IDs and intercepts
a_idswithranef <- data.frame(
  ids_spp = rownames(ids_df),
  fit_a_ids_spp = ids_df[["(Intercept)"]]
)
# clean ids col
a_idswithranef$ids <- sub("([^:]+):.*", "\\1", a_idswithranef$ids_spp)
a_idswithranef$spp <- sub(".*:(.*)", "\\1", a_idswithranef$ids_spp)

messyinter <- as.data.frame(posterior_interval(fit1))
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
simcoeftoplot2 <- simcoef[, c("ids", "spp", "a_ids_spp", "a_spp")]
colnames(simcoeftoplot2) <- c("ids", "spp", "sim_a_ids_spp", "a_spp")
# --- --- ---  --- --- ---  --- --- ---  --- --- ---  --- --- --- 
simcoeftoplot2$test <- simcoeftoplot2$a_spp+simcoeftoplot2$sim_a_ids_spp
# --- --- ---  --- --- ---  --- --- ---  --- --- ---  --- --- --- 
a_ids_mergedwithranef <- merge(simcoeftoplot2, a_ids_mergedwithranef, by = c("ids", "spp"))

# --- --- ---  --- --- ---  --- --- ---  --- --- ---  --- --- --- 

# plot!
plot_a_ids_mergedwithranef_nested_a <- ggplot(a_ids_mergedwithranef, 
                                              aes(x = test, y = fit_a_ids_spp)) +
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





ppred <- posterior_predict(fit1)  
simcoef$fit_ringwidth <- colMeans(ppred)

# ringwidthrecovery <- 
  ggplot(simcoef,aes(ringwidth, fit_ringwidth)) +
  geom_point(color = "blue", size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(y = "fit_ringwidth", x = "sim_ringwidth") +
  theme_minimal()
