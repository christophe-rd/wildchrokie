# Wildchrokie model -- leaf out prediction
# CRD 13 November 2025
# Estimate GDD until leafout (or budburst) using an intercept only model with species, ID and provenance.

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

# run models
runsimdata <- FALSE

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


if (runsimdata) {
# Simulate data gdd at leafout ####
set.seed(124)
a <- 150
sigma_y <- 3
sigma_aspp <- 6
sigma_atreeid <- 2
sigma_asite <- 4

n_site <- 10 # number of sites
n_spp <- 10 # number of species
n_perspp <- 5 # number of individuals per species
n_treeid <- n_perspp * n_spp * n_site # number of treeid
n_meas <- 10 # repeated measurements per id
N <- n_treeid * n_meas # total number of measurements

# get replicated treeid
treeid <- rep(1:n_treeid, each = n_meas)
# non replicated treeid
treeidnonrep <- rep(rep(1:n_perspp, times = n_spp), times = n_site)
# replicated spp
spp <- rep(rep(rep(1:n_spp, each = n_perspp), times = n_site), each = n_meas) 
# non replicated spp
spp_nonrep <- rep(rep(1:n_spp, each = n_perspp), each = n_site) 
# replicated site
site <- rep(rep(rep(1:n_site, each = n_spp), each = n_perspp), each = n_meas)
# non replicated site
site_nonrep <- rep(rep(1:n_site, each = n_spp), each = n_perspp)

sim <- data.frame(
  site = site,
  spp = spp,
  treeid = treeid
)
sim

# get intercept values for each parameter
aspp <- rnorm(n_spp, 0, sigma_aspp)
asite <- rnorm(n_site, 0, sigma_asite)
atreeid <- rnorm(n_treeid, 0, sigma_atreeid)

# Add my parameters to the df
sim$atreeid <- atreeid[treeid]
sim$asite <- asite[sim$site]
sim$aspp <- aspp[sim$spp]

# add the rest
sim$a <- a
sim$sigma_y <- sigma_y

sim$sigma_atreeid <- sigma_atreeid
sim$sigma_aspp <- sigma_aspp
sim$sigma_asite <- sigma_asite
sim$error <- rnorm(N, 0, sigma_y)

# adding both options of tree rings
sim$gddleafout <- 
  sim$asite + 
  sim$aspp + 
  sim$atreeid +
  sim$a +
  sim$error
sim

hist(sim$gddleafout)
 
sim$spp <- as.factor(sim$spp)
sim$sppname <- paste("spp", sim$spp)
sim$site <- as.factor(sim$site)
sim$sitename <- paste("site", sim$site)

sim$a_asite <- sim$a + sim$asite
sim$a_aspp <- sim$a + sim$aspp

# plot sim data
ggplot(sim) +
  geom_histogram(aes(gddleafout, color = sitename, fill = sitename), 
                 binwidth = 2) +
  geom_vline(aes(xintercept = a_asite, linetype = "Site mean"),
             linewidth = 0.9, alpha = 0.8, color = "black") +
  geom_vline(aes(xintercept = a_aspp, linetype = "Species mean"),
             linewidth = 0.9, color = "black") +
  geom_vline(aes(xintercept = a, linetype = "Grand mean"),
             linewidth = 1.2, color = "black") +
  facet_wrap(sitename ~ sppname) +
  labs(y = "",
       title = "gdd at leafout",
       linetype = "line type") +
  # scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  scale_linetype_manual(values = c(
    "Grand mean" = "solid",
    "Site mean" = "dotted",
    "Species mean" = "dashed"
  )) +
  theme_minimal() +
  theme(
    legend.key.height = unit(4, "lines")
  )
ggsave("figures/gddLeafout_simData/sim_gddLeafout4x4.jpeg", width = 8, height = 6, units = "in", dpi = 300)

ggplot(sim) +
  geom_histogram(aes(gddleafout, color = sitename, fill = sitename), 
                 binwidth = 2) +
  # geom_vline(aes(xintercept = a_asite, linetype = "Site mean"),
  #            linewidth = 0.9, alpha = 0.8, color = "black") +
  geom_vline(aes(xintercept = a_aspp, linetype = "Species mean"),
             linewidth = 0.9, color = "black") +
  geom_vline(aes(xintercept = a, linetype = "Grand mean"),
             linewidth = 1.2, color = "black") +
  facet_wrap(~sppname) +
  labs(y = "",
       title = "gdd at leafout",
       linetype = "line type") +
  # scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  scale_linetype_manual(values = c(
    "Grand mean" = "solid",
    "Site mean" = "dotted",
    "Species mean" = "dashed"
  )) +
  theme_minimal() +
  theme(
    legend.key.height = unit(4, "lines")
  )
ggsave("figures/gddLeafout_simData/sim_gddLeafoutSppFacet.jpeg", width = 8, height = 6, units = "in", dpi = 300)

ggplot(sim) +
  geom_density(aes(gddleafout, color = sitename), 
                 linewidth = 1) +
  # geom_vline(aes(xintercept = a_asite, linetype = "Site mean"),
  #            linewidth = 0.9, alpha = 0.8, color = "black") +
  geom_vline(aes(xintercept = a_aspp, linetype = "Species mean"),
             linewidth = 0.9, color = "black") +
  geom_vline(aes(xintercept = a, linetype = "Grand mean"),
             linewidth = 1.2, color = "black") +
  facet_wrap(~sppname) +
  labs(y = "",
       title = "gdd at leafout",
       linetype = "line type") +
  # scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  scale_linetype_manual(values = c(
    "Grand mean" = "solid",
    "Site mean" = "dotted",
    "Species mean" = "dashed"
  )) +
  theme_minimal() +
  theme(
    legend.key.height = unit(4, "lines")
  )
ggsave("figures/gddLeafout_simData/sim_gddLeafoutSppFacet_density.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# plot gdd
# gdd$year <- as.factor(gdd$year)
# ggplot(gdd)  +
#   geom_point(aes(x = doy, y = GDD_10, color = year)) + 
#   geom_vline(xintercept = mean(emp$leafout))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Run model on sim data #####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
y <- sim$gddleafout
N <- nrow(sim) 
Nspp <- length(unique(sim$spp))
species <- as.numeric(as.character(sim$spp))
Nsite <- length(unique(sim$site))
site <- as.numeric(as.character(sim$site))
Ntreeid <- length(unique(sim$treeid))
treeid <- treeid

N 
table(treeid, species)

if(FALSE){
fit <- stan("stan/modelGDDatLeafout.stan", 
            data=c("N","y","Nspp","species","Nsite", "site", "Ntreeid", "treeid"),
            iter=4000, chains=4, cores=4)
writeRDS(fit, "output/stanOutput/gddLeafout_simData/fit")
}

fit@model_pars
fit@inits
pairs(fit, pars = c("sigma_aspp", "aspp"))
pairs(fit, pars = c("sigma_asite", "asite"))


# === === === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]

sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]

sigma_df2 <- data.frame(
  sigma = character(ncol(sigma_df)),
  mean = numeric(ncol(sigma_df)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)
sigma_df2

for (i in 1:ncol(sigma_df)) { # i = 1
  sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
  sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
  sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.05), 3)
  sigma_df2$per25[i] <- round(quantile(sigma_df[[i]], probs = 0.25), 3)
  sigma_df2$per75[i] <- round(quantile(sigma_df[[i]], probs = 0.75), 3)
  sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.95), 3)
}

sigma_df2$sim_sigma <- c(
  sigma_aspp, 
  sigma_asite, 
  sigma_atreeid, 
  sigma_y)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######
# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
# remove sigma_asp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_atreeid = numeric(ncol(treeid_df)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_atreeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_atreeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_atreeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_atreeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_atreeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover aspp  ######
aspp_cols <- colnames(df_fit)[grepl("aspp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("zaspp", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("sigma", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("aspp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_aspp = numeric(ncol(aspp_df)),  
  fit_aspp_per5 = NA, 
  fit_aspp_per25 = NA,
  fit_aspp_per75 = NA,
  fit_aspp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_aspp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_aspp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_aspp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_aspp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_aspp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]
# remove sigma_aspp for now
site_cols <- site_cols[2:length(site_cols)]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_asite = numeric(ncol(site_df)),  
  fit_asite_per5 = NA, 
  fit_asite_per25 = NA,
  fit_asite_per75 = NA,
  fit_asite_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_asite[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_asite_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_asite_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_asite_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_asite_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2

# === === === === === === === #
# Plot parameter recovery #####
# === === === === === === === #

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot sigmas ######
sigma_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_errorbar(aes(ymin = per25, ymax = per75),
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 1) +
  geom_errorbar(aes(ymin = per5, ymax = per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  geom_point(color = "#046C9A", size = 3) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "") +
  theme_minimal()
sigma_plot
ggsave("figures/gddLeafout_simData/sigma_simXfit_plot.jpeg", sigma_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
# add sim to fit treeid df
treeidtoplot <- merge(
  sim[!duplicated(sim$treeid), 
          c("treeid", "atreeid")], 
  treeid_df2[!duplicated(treeid_df2$treeid), 
             c("treeid", "fit_atreeid", 
               "fit_atreeid_per5", 
               "fit_atreeid_per25",
               "fit_atreeid_per75",
               "fit_atreeid_per95")], 
  by = "treeid"
)
treeidtoplot
# plot treeid
atreeid_plot <- ggplot(treeidtoplot, aes(x = atreeid, y = fit_atreeid)) +
  geom_errorbar(aes(ymin = fit_atreeid_per5, ymax = fit_atreeid_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "") +
  theme_minimal()
atreeid_plot
# ggsave!
ggsave("figures/gddLeafout_simData/atreeid_simXfit_plot.jpeg", atreeid_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot aspp ######
aspptoplot <- merge(
  sim[!duplicated(sim$spp), 
          c("spp", "aspp")], 
  aspp_df2[!duplicated(aspp_df2$spp), 
           c("spp", "fit_aspp", 
             "fit_aspp_per5", 
             "fit_aspp_per25", 
             "fit_aspp_per75", 
             "fit_aspp_per95")], 
  by = "spp"
)
aspptoplot

aspp_plot <- ggplot(aspptoplot, aes(x = aspp, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per5, ymax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim aspp", y = "fit aspp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
aspp_plot
# ggsave!
ggsave("figures/gddLeafout_simData/aspp_simXfit_plot.jpeg", aspp_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######
sitetoplot <- merge(
  sim[!duplicated(sim$site), 
          c("site", "asite")], 
  site_df2[!duplicated(site_df2$site), 
           c("site", "fit_asite", 
             "fit_asite_per5", 
             "fit_asite_per25", 
             "fit_asite_per75", 
             "fit_asite_per95")], 
  by = "site"
)
sitetoplot

asite_plot <- ggplot(sitetoplot, aes(x = asite, y = fit_asite)) +
  geom_errorbar(aes(ymin = fit_asite_per5, ymax = fit_asite_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_asite_per25, ymax = fit_asite_per75), 
                width = 0, linewidth = 1.5, color = "darkgray", alpha=0.9) +
  geom_point(color = "#046C9A", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim asite", y = "fit asite", title = "") +
  theme_minimal()
asite_plot
# ggsave!
ggsave("figures/gddLeafout_simData/asite_simXfit_plot.jpeg", asite_plot, width = 6, height = 6, units = "in", dpi = 300)

##### Combine plots #####
combined_plot <- (atreeid_plot + sigma_plot) /
  (aspp_plot + asite_plot)
combined_plot
ggsave("figures/gddLeafout_simData/combinedPlots.jpeg", combined_plot, width = 10, height = 8, units = "in", dpi = 300)

##### Parameterization diagnostics #####

# Parameterization for treeid
treeid_df3 <- treeid_df

# subset of 30 treeid
treeid_df3 <- treeid_df[, sample(length(unique(treeid)), 40)]

colnames(treeid_df3) <- paste("treeid", colnames(treeid_df3), sep = "")
sigmaXtreeid <- cbind(sigma_df, treeid_df3)

predictors <- colnames(treeid_df3)

# jpeg("figures/gddLeafout_simData/treeidParameterization.jpg", width = 4000, height = 8000, 
#      units = "px", res = 300)
pdf("figures/gddLeafout_simData/treeidParameterization.pdf", width = 12, height = 48)
par(mfrow = c(10,5)) # unique(treeid)

for (p in predictors) {
  plot(
    sigmaXtreeid[[p]],
    log(sigmaXtreeid$sigma_atreeid),
    xlab = p,
    ylab = "log(sigmatreeid)",
    pch = 16,
    col = adjustcolor("#B40F20", alpha.f = 0.09)
  )
}
dev.off()

# Parameterization for aspp
aspp_df3 <- aspp_df
colnames(aspp_df3) <- paste("aspp", colnames(aspp_df3), sep = "")
sigmaXaspp <- cbind(sigma_df, aspp_df3)

predictors <- colnames(aspp_df3)

jpeg("figures/gddLeafout_simData/asppParameterization.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
par(mfrow = c(3, 3))

for (p in predictors) {
  plot(
    sigmaXaspp[[p]],
    log(sigmaXaspp$sigma_aspp),
    xlab = p,
    ylab = "log(sigma_aspp)",
    pch = 16,
    col = adjustcolor("#B40F20", alpha.f = 0.09)
  )
}
dev.off()

# Parameterization for site
site_df3 <- site_df
colnames(site_df3) <- paste("site", colnames(site_df3), sep = "")
sigmaXsite <- cbind(sigma_df, site_df3)

predictors <- colnames(site_df3)

jpeg("figures/gddLeafout_simData/siteParameterization.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
par(mfrow = c(2, 2))

for (p in predictors) {
  
  plot(
    sigmaXsite[[p]],
    log(sigmaXsite$sigma_asite),
    xlab = p,
    ylab = "log(sigmasite)",
    pch = 16,
    col = adjustcolor("#B40F20", alpha.f = 0.09)
  )
}
dev.off()

}

# RUN MODEL ON EMPIRICAL DATA ####
# empirical data
emp <- read.csv("output/empiricalDataNORing.csv")
gdd <- read.csv("output/gddByYear.csv")

no_naleafout <- emp[!is.na(emp$leafoutGDD),]
nrow(no_naleafout)

# some checks
checks <- subset(no_naleafout, leafoutGDD > 300)
checks # there are two entries for very late leafout dates...

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

table(treeid, species)
table(treeid)
table(species)


fit <- stan("stan/modelGDDatLeafout.stan", 
            data=c("N","y",
                   "Nspp","species",
                   "Nsite","site",
                   "Ntreeid", "treeid"),
            iter=4000, chains=4, cores=4)

saveRDS(fit, "output/stanOutput/gddLeafout_empData_fit")

# Diagnostics for partial pooling ####
if (FALSE) {
diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fit)

# atreeid
atreeid <- names(samples)[grepl("zatreeid", names(samples))]
atreeid <- atreeid[!grepl("sigma", atreeid)]
atreeid <- atreeid[sample(length(unique(atreeid)), 21)]
pdf("figures/gddLeafout_empData/atreeidParameterization.pdf", width = 6, height = 18)
util$plot_div_pairs(atreeid, "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
dev.off()
}

# === === === === === === === === === === === === === === === === === === === ==
# Recover parameters ####
# === === === === === === === === === === === === === === === === === === === ==
df_fit <- as.data.frame(fit)

# get values to original scale
for (i in 1:ncol(df_fit)){
  df_fit[[i]] <- df_fit[[i]] * scale
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
unique(colnames(df_fit))
sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]

sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]

sigma_df2 <- data.frame(
  sigma = character(ncol(sigma_df)),
  mean = numeric(ncol(sigma_df)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)
sigma_df2

for (i in 1:ncol(sigma_df)) { # i = 1
  sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
  sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
  sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.05), 3)
  sigma_df2$per25[i] <- round(quantile(sigma_df[[i]], probs = 0.25), 3)
  sigma_df2$per75[i] <- round(quantile(sigma_df[[i]], probs = 0.75), 3)
  sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.95), 3)
}
sigma_df2


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######
# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]
# remove sigma_aspp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_atreeid = numeric(ncol(treeid_df)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_atreeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_atreeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_atreeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_atreeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_atreeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover aspp  ######
aspp_cols <- colnames(df_fit)[grepl("aspp", colnames(df_fit))]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("aspp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_aspp = numeric(ncol(aspp_df)),  
  fit_aspp_per5 = NA, 
  fit_aspp_per25 = NA,
  fit_aspp_per75 = NA,
  fit_aspp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_aspp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_aspp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_aspp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_aspp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_aspp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_asite = numeric(ncol(site_df)),  
  fit_asite_per5 = NA, 
  fit_asite_per25 = NA,
  fit_asite_per75 = NA,
  fit_asite_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_asite[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_asite_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_asite_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_asite_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_asite_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2

df_fit <- as.data.frame(fit)

# Prior checks ####
###### a ######
a_posterior <- df_fit[, colnames(df_fit) %in% "a"]

a_prior <- rnorm(1e4, 8, 2)

mu8 <- ggplot() +
  geom_density(data = data.frame(a = a_prior*20),
               aes(x = a, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = data.frame(value = a_posterior*20),
               aes(x = value, colour = "Posterior"),
               linewidth = 0.8) +
  labs(title = "prior at mu = 160",
       x = "a", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
mu8
ggsave("figures/gddLeafout_empData/a_prior_mu8.jpeg", 
       mu8,
       width = 8, height = 6, units = "in", dpi = 300)

###### aspp ######

# convert posterior distribution to long format
aspp_long <- reshape(
  aspp_df,
  direction = "long",
  varying = list(names(aspp_df)),
  v.names = "value",
  timevar = "spp",
  times = names(aspp_df),
  idvar = "draw"
)
aspp_long

# aspp prior
aspp_prior <- rnorm(1e4, 0, 1)

ggplot() +
  geom_density(data = data.frame(aspp_prior = aspp_prior),
               aes(x = aspp_prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = aspp_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_atreeid",
       x = "aspp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()


###### asite ######

# convert posterior distribution to long format
asite_long <- reshape(
  site_df,
  direction = "long",
  varying = list(names(site_df)),
  v.names = "value",
  timevar = "site",
  times = names(site_df),
  idvar = "draw"
)
asite_long

# asite prior
asite_prior <- rnorm(1e4, 0, 0.5)

ggplot() +
  geom_density(data = data.frame(asite_prior = asite_prior),
               aes(x = asite_prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = asite_long,
               aes(x = value, colour = "Posterior", group = site),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_atreeid",
       x = "asite", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()


###### atreeid ######
# plot
atreeid_long <- reshape(
  treeid_df,
  direction = "long",
  varying = list(names(treeid_df)),
  v.names = "value",
  timevar = "spp",
  times = names(treeid_df),
  idvar = "draw"
)
atreeid_long

# sample to plot
vec <- sample(treeid, 20)
sub_atreeid <- subset(atreeid_long, spp %in% vec)

# simulate priors
hyperparameter_draws <- 1000
parameter_draws <- 1000
n_sigma_atreeid <- 200

# set to prior values
sigma_atreeid_vec <- abs(rnorm(n_sigma_atreeid, 0, 0.5))

prior_atreeid <- rep(NA, parameter_draws*length(sigma_atreeid_vec))

for (i in 1: length(sigma_atreeid_vec)) {
  prior_atreeid[((i - 1)*parameter_draws + 1):(i*parameter_draws)] <- rnorm(parameter_draws, 0, sigma_atreeid_vec[i])
}
prior_atreeid

ggplot() +
  geom_density(data = data.frame(prior_atreeid = prior_atreeid),
               aes(x = prior_atreeid, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = sub_atreeid,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.1) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_atreeid",
       x = "atreeid", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Generate point estimates ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# create an empty data frame in which I have more years than I currently have
length(unique(no_naleafout$treeid_num))
dfpoints <- data.frame(
  year = rep(2017:2024, times = 105),
  treeid_num = rep(unique(no_naleafout$treeid_num), each = 8)
)
test <- merge(dfpoints, 
              no_naleafout[, c("year", "treeid_num", "spp_num", "site_num", "leafoutGDD")],   
              by = c("year", "treeid_num"),
              all.x = TRUE)

# fill in the remaining spp and site by matching to of df
test$site_num <- no_naleafout$site_num[match(test$treeid_num, no_naleafout$treeid_num)]
test$spp_num <- no_naleafout$spp_num[match(test$treeid_num, no_naleafout$treeid_num)]

# add parameter mean posterior estimates 
test$a <- mean(df_fit[,"a"])*20
test$aspp <- aspp_df2$fit_aspp[match(test$spp_num, aspp_df2$spp)]
test$asite <- site_df2$fit_asite[match(test$site_num, site_df2$site)]
test$atreeid <- treeid_df2$fit_atreeid[match(test$treeid_num, treeid_df2$treeid)]
test$estimates <- test$a + test$aspp + test$asite + test$atreeid

# quickly plot results to compare to observed data
forplot <- test[!is.na(test$leafoutGDD),]
par(mfrow=c(1,2))
hist(forplot$leafoutGDD)
hist(forplot$estimates) # at least it's consistent with my generated quantities bloc

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
samples <- util$extract_expectand_vals(fit)

# back convert the object samples to the original gdd scale
samples2 <- lapply(samples, function(x) x*scale)
y2 <- y*scale

jpeg(
  filename = "figures/gddLeafout_empData/retrodictiveHist_lowerBound160.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
util$plot_hist_quantiles(samples2, "y_rep", 
                         -50, # lower x axis limit
                         700, # upper x axis limit
                         20, # binning
                         baseline_values = y2,
                         xlab = "gdd at leafout")
dev.off()

