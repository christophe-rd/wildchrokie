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

if (length(grep("christophe_rouleau-desrochers", getwd())) > 0) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if (length(grep("lizzie", getwd())) > 0) {
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
} else  {
  setwd("/home/crouleau/wildchrokie/analyses")
}

# empirical data
emp <- read.csv("output/empiricalDataMAIN.csv")

# gdd data
gdd <- read.csv("output/gddByYear.csv")


# Simulate data gdd at leafout ####
set.seed(124)
a <- 150
sigma_y <- 3
sigma_asp <- 6
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
asp <- rnorm(n_spp, 0, sigma_asp)
asite <- rnorm(n_site, 0, sigma_asite)
atreeid <- rnorm(n_treeid, 0, sigma_atreeid)

# Add my parameters to the df
sim$atreeid <- atreeid[treeid]
sim$asite <- asite[sim$site]
sim$asp <- asp[sim$spp]

# add the rest
sim$a <- a
sim$sigma_y <- sigma_y

sim$sigma_atreeid <- sigma_atreeid
sim$sigma_asp <- sigma_asp
sim$sigma_asite <- sigma_asite
sim$error <- rnorm(N, 0, sigma_y)

# adding both options of tree rings
sim$gddleafout <- 
  sim$asite + 
  sim$asp + 
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
sim$a_asp <- sim$a + sim$asp

# plot sim data
ggplot(sim) +
  geom_histogram(aes(gddleafout, color = sitename, fill = sitename), 
                 binwidth = 2) +
  geom_vline(aes(xintercept = a_asite, linetype = "Site mean"),
             linewidth = 0.9, alpha = 0.8, color = "black") +
  geom_vline(aes(xintercept = a_asp, linetype = "Species mean"),
             linewidth = 0.9, color = "black") +
  geom_vline(aes(xintercept = a, linetype = "Grand mean"),
             linewidth = 1.2, color = "black") +
  facet_wrap(sitename ~ sppname) +
  labs(y = "",
       title = "gdd at leafout",
       linetype = "line type") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
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
  geom_vline(aes(xintercept = a_asp, linetype = "Species mean"),
             linewidth = 0.9, color = "black") +
  geom_vline(aes(xintercept = a, linetype = "Grand mean"),
             linewidth = 1.2, color = "black") +
  facet_wrap(~sppname) +
  labs(y = "",
       title = "gdd at leafout",
       linetype = "line type") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
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
  geom_vline(aes(xintercept = a_asp, linetype = "Species mean"),
             linewidth = 0.9, color = "black") +
  geom_vline(aes(xintercept = a, linetype = "Grand mean"),
             linewidth = 1.2, color = "black") +
  facet_wrap(~sppname) +
  labs(y = "",
       title = "gdd at leafout",
       linetype = "line type") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
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
gdd$year <- as.factor(gdd$year)
ggplot(gdd)  +
  geom_point(aes(x = doy, y = GDD_10, color = year)) + 
  geom_vline(xintercept = mean(emp$leafout))

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


table(treeid, species)

if(FALSE){
fit <- stan("stan/modelGDDatLeafout.stan", 
            data=c("N","y","Nspp","species","Nsite", "site", "Ntreeid", "treeid"),
            iter=4000, chains=4, cores=4)

}

fit@model_pars
fit@inits
pairs(fit, pars = c("sigma_asp", "asp"))
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
  sigma_asp, 
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
###### Recover a spp  ######
asp_cols <- colnames(df_fit)[grepl("asp", colnames(df_fit))]
asp_cols <- asp_cols[!grepl("zasp", asp_cols)]
asp_cols <- asp_cols[!grepl("sigma", asp_cols)]

asp_df <- df_fit[, colnames(df_fit) %in% asp_cols]
# change their names
colnames(asp_df) <- sub("asp\\[(\\d+)\\]", "\\1", colnames(asp_df))
#empty asp df
asp_df2 <- data.frame(
  spp = character(ncol(asp_df)),
  fit_asp = numeric(ncol(asp_df)),  
  fit_asp_per5 = NA, 
  fit_asp_per25 = NA,
  fit_asp_per75 = NA,
  fit_asp_per95 = NA
)
for (i in 1:ncol(asp_df)) { # i = 1
  asp_df2$spp[i] <- colnames(asp_df)[i]         
  asp_df2$fit_asp[i] <- round(mean(asp_df[[i]]),3)  
  asp_df2$fit_asp_per5[i] <- round(quantile(asp_df[[i]], probs = 0.05), 3)
  asp_df2$fit_asp_per25[i] <- round(quantile(asp_df[[i]], probs = 0.25), 3)
  asp_df2$fit_asp_per75[i] <- round(quantile(asp_df[[i]], probs = 0.75), 3)
  asp_df2$fit_asp_per95[i] <- round(quantile(asp_df[[i]], probs = 0.95), 3)
}
asp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a site ######
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]
# remove sigma_asp for now
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
###### Plot a spp ######
asptoplot <- merge(
  sim[!duplicated(sim$spp), 
          c("spp", "asp")], 
  asp_df2[!duplicated(asp_df2$spp), 
           c("spp", "fit_asp", 
             "fit_asp_per5", 
             "fit_asp_per25", 
             "fit_asp_per75", 
             "fit_asp_per95")], 
  by = "spp"
)
asptoplot

asp_plot <- ggplot(asptoplot, aes(x = asp, y = fit_asp)) +
  geom_errorbar(aes(ymin = fit_asp_per5, ymax = fit_asp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_asp_per25, ymax = fit_asp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim asp", y = "fit asp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
asp_plot
# ggsave!
ggsave("figures/gddLeafout_simData/asp_simXfit_plot.jpeg", asp_plot, width = 6, height = 6, units = "in", dpi = 300)

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
  (asp_plot + asite_plot)
combined_plot
ggsave("figures/gddLeafout_simData/combinedPlots.jpeg", combined_plot, width = 10, height = 8, units = "in", dpi = 300)

# Parameterization diagnostics ####

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
    ylab = "log(sigma_treeid)",
    pch = 16,
    col = adjustcolor("#B40F20", alpha.f = 0.09)
  )
}
dev.off()

# Parameterization for asp
asp_df3 <- asp_df
colnames(asp_df3) <- paste("asp", colnames(asp_df3), sep = "")
sigmaXasp <- cbind(sigma_df, asp_df3)

predictors <- colnames(asp_df3)

jpeg("figures/gddLeafout_simData/aspParameterization.jpg", width = 2000, height = 2000, 
     units = "px", res = 300)
par(mfrow = c(3, 3))

for (p in predictors) {
  plot(
    sigmaXasp[[p]],
    log(sigmaXasp$sigma_asp),
    xlab = p,
    ylab = "log(sigma_asp)",
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
    ylab = "log(sigma_site)",
    pch = 16,
    col = adjustcolor("#B40F20", alpha.f = 0.09)
  )
}
dev.off()
