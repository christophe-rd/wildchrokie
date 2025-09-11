# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 3)
quartz()

# Load library 
library(rstanarm)
library(ggplot2)
library(rethinking)
library("wesanderson")

runmodels <- FALSE
runoldcode <- FALSE

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")


# === === === === === === === === === === === === === === === === 
#### Step 1. Come up with a model ####
# === === === === === === === === === === === === === === === === 
# As discussed with Victor, we will use a linear model to predict growth from GDD. 

# === === === === === === === === === === === === === === === === 
#### Step 2. Simulate data ####
# === === === === === === === === === === === === === === === === 
# set parameters

# set parameters
set.seed(124)
a <- 1.5
sigma_y <- 0.2
sigma_a_spp <- 0.3 # This is pretty low, but I guess you think your species are closely related and will be similar?
sigma_a_ids <- 0.8
sigma_a_site <- 0.1

n_site <- 4 # number of sites
n_spp <- 10 # number of species
n_perspp <- 10 # number of individuals per species
n_ids <- n_perspp * n_spp * n_site # number of ids
n_meas <- 5 # repeated measurements per id
N <- n_ids * n_meas # total number of measurements


# get replicated ids
ids <- rep(rep(rep(1:n_perspp, times = n_spp), times = n_site), each = n_meas)
# non replicated ids
idsnonrep <- rep(rep(1:n_perspp, times = n_spp), times = n_site)
# replicated spp
spp <- rep(rep(rep(1:n_spp, each = n_perspp), times = n_site), each = n_meas) 
# non replicated spp
spp_nonrep <- rep(rep(1:n_spp, each = n_perspp), each = n_site) 
# replicated site
site <- rep(rep(rep(1:n_site, each = n_spp), each = n_perspp), each = n_meas)
# non replicated site
site_nonrep <- rep(rep(1:n_site, each = n_spp), each = n_perspp)
# quick check 
table(idsnonrep, site_nonrep)

simcoef <- data.frame(
  site = site,
  spp = spp,
  ids = ids
)

# get 50 intercept values for each species
a_spp <- rnorm(n_spp, 0, sigma_a_spp)

a_site <- rnorm(n_site, 0, sigma_a_site)

# Option 2
a_ids <- rnorm(n_ids, 0, sigma_a_ids)

# Add my parameters to the df
simcoef$a_site <- a_site[simcoef$site]
simcoef$a_spp <- a_spp[simcoef$spp]

### add index to allow ids to be input at the right place
id_index <- rep(1:n_ids, each = n_meas)

simcoef$a_ids <- a_ids[id_index]

# add the rest of the boring stuff 
simcoef$a <- 1.5
simcoef$b <- 0.4
simcoef$sigma_y <- sigma_y
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gddcons <- rnorm(N, 1800, 100)/200

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$a_site + 
  simcoef$a_spp + 
  simcoef$a_ids + 
  simcoef$a + 
  (simcoef$b*simcoef$gddcons) + 
  simcoef$error

# prepare grouping factors for stan_lmer (ensure ids are unique within spp)
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
simcoef$ids <- factor(simcoef$ids)
simcoef$ids_uni <- paste(simcoef$site, simcoef$spp, simcoef$ids, sep = "_")

# === === === === === #
##### Run models #####
# === === === === === #

###### Model nested on the intercept #######
fitnestedrun <- FALSE
if(fitnestedrun) {
  fitnested <- stan_lmer(
    ringwidth ~ 1 + gddcons + 
      (1|site) +
      (1|spp) +
      (1|ids_uni),
    data = simcoef,
    chains = 4,
    iter = 4000,
    core=4
  )
  saveRDS(fitnested, "output/fitnested")
}
fitnested

y <- simcoef$ringwidth
N <- nrow(simcoef)
gdd <- simcoef$gddcons
Nspp <- length(unique(simcoef$spp))
Nsite <- length(unique(simcoef$site))
site <- as.numeric(as.character(simcoef$site))
species <- as.numeric(as.character(simcoef$sp))
treeid <- rep(1:length(unique(simcoef$ids_uni)), each=5)
Ntreeid <- length(unique(treeid))
table(treeid)

library(rstan)
fit <- stan("stan/twolevelhierint.stan", 
  data=c("N","y","Nspp","species","Nsite", "site", "Ntreeid", "treeid", "gdd"), 
  iter=4000, chains=4, cores=4)
summary(fit)$summary
fitpost <- extract(fit)

# === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === #
fitnested <- readRDS("output/fitnested")
df_fit <- as.data.frame(fitnested)

# recover slope
colnames(df_fit)
# grab ids nested in spp
ids_cols <- colnames(df_fit)[grepl("ids_uni", colnames(df_fit))]
ids_cols <- ids_cols[1:length(ids_cols)-1]
ids_df <- df_fit[, colnames(df_fit) %in% ids_cols]
# change their names
colnames(ids_df) <- sub(".*ids_uni:([^]]+)\\]$", "\\1", colnames(ids_df))
# empty ids dataframe
ids_df2 <- data.frame(
  ids_uni = character(ncol(ids_df)),
  fit_a_ids = numeric(ncol(ids_df)),  
  fit_a_ids_per5 = NA, 
  fit_a_ids_per95 = NA
)
for (i in 1:ncol(ids_df)) { # i = 1
  ids_df2$ids_uni[i] <- colnames(ids_df)[i]         
  ids_df2$fit_a_ids[i] <- round(mean(ids_df[[i]]),3)  
  ids_df2$fit_a_ids_per5[i] <- round(quantile(ids_df[[i]], probs = 0.055), 3)
  ids_df2$fit_a_ids_per95[i] <- round(quantile(ids_df[[i]], probs = 0.945), 3)
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
  fit_a_spp_per5 = NA, 
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(spp_df)) { # i = 1
  spp_df2$spp[i] <- colnames(spp_df)[i]         
  spp_df2$fit_a_spp[i] <- round(mean(spp_df[[i]]),3)  
  spp_df2$fit_a_spp_per5[i] <- round(quantile(spp_df[[i]], probs = 0.055), 3)
  spp_df2$fit_a_spp_per95[i] <- round(quantile(spp_df[[i]], probs = 0.945), 3)
}
spp_df2

# grab site
site_cols <- colnames(df_fit)[grepl(" site:", colnames(df_fit))]
site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub(".*site:([0-9]+).*", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.055), 3)
  site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.945), 3)
}
site_df2

# clean and re-arrange in a single df!
# create ids, spp and site cols
tmp <- do.call(rbind, strsplit(ids_df2$ids_uni, "_"))
# re-add the new cols
ids_df2$site <- as.numeric(tmp[,1])
ids_df2$spp  <- as.numeric(tmp[,2])
ids_df2$ids  <- as.numeric(tmp[,3])

# merge spp and ids
merge1 <- merge(spp_df2, ids_df2, by = "spp")
merge2 <- merge(merge1, site_df2, by = "site")


# === === === === === === === #
# Plot parameter recovery #####
# === === === === === === === #

# Start with ids ######
# merge simcoef and model
idstoplot <- merge(
  simcoef[!duplicated(simcoef$ids_uni), 
          c("ids_uni", "a_ids")], 
                   
  merge2[!duplicated(merge2$ids_uni), 
         c("ids_uni", "fit_a_ids", "fit_a_ids_per5", "fit_a_ids_per95")], 
  by = "ids_uni"
  )

# plot ids
a_ids_simXfit_plot <- ggplot(idstoplot, aes(x = a_ids, y = fit_a_ids)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_ids_per5, ymax = fit_a_ids_per95), width = 0, color = "darkgray", alpha=0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_ids", y = "fit a_ids", title = "") +
  theme_minimal()
a_ids_simXfit_plot
# ggsave!
ggsave("figures/a_ids_simXfit_plot.jpeg", a_ids_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# === === === === === === === === === === === === === === === === === === === ===

# plot spp
spptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "a_spp")], 
  merge2[!duplicated(merge2$spp), 
         c("spp", "fit_a_spp", "fit_a_spp_per5", "fit_a_spp_per95")], 
  by = "spp"
  )
spptoplot

a_spp_simXfit_plot <- ggplot(spptoplot, aes(x = a_spp, y = fit_a_spp)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_spp_per5, ymax = fit_a_spp_per95), width = 0, color = "darkgray", alpha=0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_spp", y = "fit a_spp", title = "") +
  theme_minimal()
a_spp_simXfit_plot
# ggsave!
ggsave("figures/a_spp_simXfit_plot.jpeg", a_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# === === === === === === === === === === === === === === === === === === === ===
sitetoplot <- merge(
  simcoef[!duplicated(simcoef$site), 
          c("site", "a_site")], 
  merge2[!duplicated(merge2$site), 
         c("site", "fit_a_site", "fit_a_site_per5", "fit_a_site_per95")], 
  by = "site"
)
sitetoplot

a_site_simXfit_plot <- ggplot(sitetoplot, aes(x = a_site, y = fit_a_site)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_site_per5, ymax = fit_a_site_per95), width = 0, color = "darkgray", alpha=0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_site", y = "fit a_site", title = "") +
  theme_minimal()
a_site_simXfit_plot
# ggsave!
ggsave("figures/a_site_simXfit_plot.jpeg", a_site_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)


# === === === === === === === === === === === === === === === === 
#### Step 3. Set your priors ####
# === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
#### Step 4. Run model on empirical data ####
# === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
#### Step 5. Perform retrodictive checks using the model fit to your empiral data ####
# === === === === === === === === === === === === === === === ===


