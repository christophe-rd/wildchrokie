# Wildchrokie model -- leaf out prediction
# CRD 17 December 2025
# Related to modelGDDLeafout script where it estimates the gdd at leafout. 
# This is an additional code to plot those results 

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


# RUN MODEL ON EMPIRICAL DATA ####
# empirical data
emp <- read.csv("output/empiricalDataNORing.csv")
gdd <- read.csv("output/gddByYear.csv")
fit <- readRDS("output/stanOutput/gddLeafout_empData_fit") 
fit_largerPriors <- readRDS("output/stanOutput/gddLeafout_empData_fit_largerPriors") 

no_naleafout <- emp[!is.na(emp$leafoutGDD),]

# give numeric ids to my groups 
no_naleafout$site_num <- match(no_naleafout$site, unique(no_naleafout$site))
no_naleafout$spp_num <- match(no_naleafout$spp, unique(no_naleafout$spp))
no_naleafout$treeid_num <- match(no_naleafout$treeid, unique(no_naleafout$treeid))

scale <- 20

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


# Plots ####
aspp_df2$spp <- as.numeric(aspp_df2$spp)
aspp_df2$spp_name <- no_naleafout$spp[match(aspp_df2$spp, no_naleafout$spp_num)]

##### aspp intercepts means #####
asppEstimates_plot <- ggplot(aspp_df2, aes(x = fit_aspp, y = spp_name, color = spp_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_aspp_per5, 
                     xmax = fit_aspp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_aspp_per25, 
                     xmax = fit_aspp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species", title = "Priors unchanged")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  # xlim(-80, 80) +
  theme_bw() +
  scale_y_discrete(limits = rev)   
asppEstimates_plot
ggsave("figures/gddLeafout_empData/sppEstimates.jpeg", asppEstimates_plot, width = 8, height = 6, units = "in", dpi = 300)

## site intercept means ####

site_df2$site <- as.numeric(site_df2$site)
site_df2$site_name <- no_naleafout$site[match(site_df2$site, no_naleafout$site_num)]

siteEstimates_plot <- ggplot(site_df2, aes(x = fit_asite, y = site_name, color = site_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_asite_per5, 
                     xmax = fit_asite_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_asite_per25, 
                     xmax = fit_asite_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Site", x = "Site intercept values", color = "Site", title = "Priors unchanged")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  # xlim(-60, 60) +
  theme_bw() +
  scale_y_discrete(limits = rev)   
siteEstimates_plot
ggsave("figures/gddLeafout_empData/siteEstimates.jpeg", siteEstimates_plot, width = 8, height = 6, units = "in", dpi = 300)

# Mean plots with atreeid ####
treeid_df2$treeid <- as.numeric(treeid_df2$treeid)
treeid_df2$treeid_name <- no_naleafout$treeid[match(treeid_df2$treeid, no_naleafout$treeid_num)]

# now do the same, but for species
treeid_df2$spp <- no_naleafout$spp[match(treeid_df2$treeid, no_naleafout$treeid_num)]

# same for site
treeid_df2$site <- no_naleafout$site[match(treeid_df2$treeid, no_naleafout$treeid_num)]


# quick check that I didn't mess anything up
un <- no_naleafout[!duplicated(no_naleafout$treeid),]
table(un$spp)
table(treeid_df2$spp)

sub <- subset(no_naleafout, select = c("treeid_num", "spp_num", "site_num"))
sub <- sub[!duplicated(sub$treeid_num),]
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# sum atreeid, aspp, asite
aspp_df3 <- as.matrix(aspp_df)
site_df3 <- as.matrix(site_df)
treeid_df3 <- as.matrix(treeid_df)

fullatreeid <-
  treeid_df3[, sub$treeid_num] +
  aspp_df3[, sub$spp_num] +
  site_df3[, sub$site_num] 
fullatreeid

fullatreeid2 <- as.data.frame(fullatreeid)
# get posterior means and quantiles

# empty treeid dataframe
treeid_df4 <- data.frame(
  treeid = character(ncol(fullatreeid2)),
  fit_atreeid = numeric(ncol(fullatreeid2)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(fullatreeid2)) { # i = 1
  treeid_df4$treeid[i] <- colnames(fullatreeid2)[i]         
  treeid_df4$fit_atreeid[i] <- round(mean(fullatreeid2[[i]]),3)  
  treeid_df4$fit_atreeid_per5[i] <- round(quantile(fullatreeid2[[i]], probs = 0.05), 3)
  treeid_df4$fit_atreeid_per25[i] <- round(quantile(fullatreeid2[[i]], probs = 0.25), 3)
  treeid_df4$fit_atreeid_per75[i] <- round(quantile(fullatreeid2[[i]], probs = 0.75), 3)
  treeid_df4$fit_atreeid_per95[i] <- round(quantile(fullatreeid2[[i]], probs = 0.95), 3)
}
treeid_df4

# get the og treeid names, spp and site back:
treeid_df4$treeid <- as.numeric(treeid_df4$treeid)
treeid_df4$treeid_name <- no_naleafout$treeid[match(treeid_df4$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4$spp_name <- no_naleafout$spp[match(treeid_df4$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4$site_name <- no_naleafout$site[match(treeid_df4$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4


if (TRUE){
# Recover parameters for larger priors ####
# === === === === === === === === === === === === === === === === === === === ==
df_fit_largerPriors <- as.data.frame(fit_largerPriors)

# get values to original scale
for (i in 1:ncol(df_fit_largerPriors)){
  df_fit_largerPriors[[i]] <- df_fit_largerPriors[[i]] * scale
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
unique(colnames(df_fit_largerPriors))
sigma_cols <- colnames(df_fit_largerPriors)[grepl("sigma", colnames(df_fit_largerPriors))]

sigma_df_largerPriors <- df_fit_largerPriors[, colnames(df_fit_largerPriors) %in% sigma_cols]

sigma_df_largerPriors2 <- data.frame(
  sigma = character(ncol(sigma_df_largerPriors)),
  mean = numeric(ncol(sigma_df_largerPriors)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)
sigma_df_largerPriors2

for (i in 1:ncol(sigma_df_largerPriors)) { # i = 1
  sigma_df_largerPriors2$sigma[i] <- colnames(sigma_df_largerPriors)[i]         
  sigma_df_largerPriors2$mean[i] <- round(mean(sigma_df_largerPriors[[i]]),3)  
  sigma_df_largerPriors2$per5[i] <- round(quantile(sigma_df_largerPriors[[i]], probs = 0.05), 3)
  sigma_df_largerPriors2$per25[i] <- round(quantile(sigma_df_largerPriors[[i]], probs = 0.25), 3)
  sigma_df_largerPriors2$per75[i] <- round(quantile(sigma_df_largerPriors[[i]], probs = 0.75), 3)
  sigma_df_largerPriors2$per95[i] <- round(quantile(sigma_df_largerPriors[[i]], probs = 0.95), 3)
}
sigma_df_largerPriors2


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######
# grab treeid 
treeid_cols <- colnames(df_fit_largerPriors)[grepl("atreeid", colnames(df_fit_largerPriors))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]
# remove sigma_aspp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

treeid_df_largerPriors <- df_fit_largerPriors[, colnames(df_fit_largerPriors) %in% treeid_cols]
# change their names
colnames(treeid_df_largerPriors) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df_largerPriors))
# empty treeid dataframe
treeid_df_largerPriors2 <- data.frame(
  treeid = character(ncol(treeid_df_largerPriors)),
  fit_atreeid = numeric(ncol(treeid_df_largerPriors)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(treeid_df_largerPriors)) { # i = 1
  treeid_df_largerPriors2$treeid[i] <- colnames(treeid_df_largerPriors)[i]         
  treeid_df_largerPriors2$fit_atreeid[i] <- round(mean(treeid_df_largerPriors[[i]]),3)  
  treeid_df_largerPriors2$fit_atreeid_per5[i] <- round(quantile(treeid_df_largerPriors[[i]], probs = 0.05), 3)
  treeid_df_largerPriors2$fit_atreeid_per25[i] <- round(quantile(treeid_df_largerPriors[[i]], probs = 0.25), 3)
  treeid_df_largerPriors2$fit_atreeid_per75[i] <- round(quantile(treeid_df_largerPriors[[i]], probs = 0.75), 3)
  treeid_df_largerPriors2$fit_atreeid_per95[i] <- round(quantile(treeid_df_largerPriors[[i]], probs = 0.95), 3)
}
treeid_df_largerPriors2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover aspp  ######
aspp_cols <- colnames(df_fit_largerPriors)[grepl("aspp", colnames(df_fit_largerPriors))]

aspp_df_largerPriors <- df_fit_largerPriors[, colnames(df_fit_largerPriors) %in% aspp_cols]
# change their names
colnames(aspp_df_largerPriors) <- sub("aspp\\[(\\d+)\\]", "\\1", colnames(aspp_df_largerPriors))
#empty aspp df
aspp_df_largerPriors2 <- data.frame(
  spp = character(ncol(aspp_df_largerPriors)),
  fit_aspp = numeric(ncol(aspp_df_largerPriors)),  
  fit_aspp_per5 = NA, 
  fit_aspp_per25 = NA,
  fit_aspp_per75 = NA,
  fit_aspp_per95 = NA
)
for (i in 1:ncol(aspp_df_largerPriors)) { # i = 1
  aspp_df_largerPriors2$spp[i] <- colnames(aspp_df_largerPriors)[i]         
  aspp_df_largerPriors2$fit_aspp[i] <- round(mean(aspp_df_largerPriors[[i]]),3)  
  aspp_df_largerPriors2$fit_aspp_per5[i] <- round(quantile(aspp_df_largerPriors[[i]], probs = 0.05), 3)
  aspp_df_largerPriors2$fit_aspp_per25[i] <- round(quantile(aspp_df_largerPriors[[i]], probs = 0.25), 3)
  aspp_df_largerPriors2$fit_aspp_per75[i] <- round(quantile(aspp_df_largerPriors[[i]], probs = 0.75), 3)
  aspp_df_largerPriors2$fit_aspp_per95[i] <- round(quantile(aspp_df_largerPriors[[i]], probs = 0.95), 3)
}
aspp_df_largerPriors2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a site ######
site_cols <- colnames(df_fit_largerPriors)[grepl("asite", colnames(df_fit_largerPriors))]

site_df_largerPriors <- df_fit_largerPriors[, colnames(df_fit_largerPriors) %in% site_cols]
# change their names
colnames(site_df_largerPriors) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df_largerPriors))
# empty site df
site_df_largerPriors2 <- data.frame(
  site = character(ncol(site_df_largerPriors)),
  fit_asite = numeric(ncol(site_df_largerPriors)),  
  fit_asite_per5 = NA, 
  fit_asite_per25 = NA,
  fit_asite_per75 = NA,
  fit_asite_per95 = NA
)
for (i in 1:ncol(site_df_largerPriors)) { # i = 1
  site_df_largerPriors2$site[i] <- colnames(site_df_largerPriors)[i]         
  site_df_largerPriors2$fit_asite[i] <- round(mean(site_df_largerPriors[[i]]),3)  
  site_df_largerPriors2$fit_asite_per5[i] <- round(quantile(site_df_largerPriors[[i]], probs = 0.05), 3)
  site_df_largerPriors2$fit_asite_per25[i] <- round(quantile(site_df_largerPriors[[i]], probs = 0.25), 3)
  site_df_largerPriors2$fit_asite_per75[i] <- round(quantile(site_df_largerPriors[[i]], probs = 0.75), 3)
  site_df_largerPriors2$fit_asite_per95[i] <- round(quantile(site_df_largerPriors[[i]], probs = 0.95), 3)
}
site_df_largerPriors2


# Plots ####
aspp_df_largerPriors2$spp <- as.numeric(aspp_df_largerPriors2$spp)
aspp_df_largerPriors2$spp_name <- no_naleafout$spp[match(aspp_df_largerPriors2$spp, no_naleafout$spp_num)]

##### aspp intercepts means #####
aspp_df_largerPrior_plot <- ggplot(aspp_df_largerPriors2, aes(x = fit_aspp, y = spp_name, color = spp_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_aspp_per5, 
                     xmax = fit_aspp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_aspp_per25, 
                     xmax = fit_aspp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species", title = "Priors with 3x sigmas")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  xlim(-80, 80) +
  theme_bw() +
  scale_y_discrete(limits = rev)   
ggsave("figures/gddLeafout_empData/sppEstimates_largerPriors.jpeg", aspp_df_largerPrior_plot, width = 8, height = 6, units = "in", dpi = 300)

## site intercept means ####

site_df_largerPriors2$site <- as.numeric(site_df_largerPriors2$site)
site_df_largerPriors2$site_name <- no_naleafout$site[match(site_df_largerPriors2$site, no_naleafout$site_num)]

site_df_largerPrior_plot <- ggplot(site_df_largerPriors2, aes(x = fit_asite, y = site_name, color = site_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_asite_per5, 
                     xmax = fit_asite_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_asite_per25, 
                     xmax = fit_asite_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Site", x = "Site intercept values", color = "Site", title = "Priors with 3x sigmas")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  xlim(-60, 60) +
  theme_bw() +
  scale_y_discrete(limits = rev)   
ggsave("figures/gddLeafout_empData/siteEstimates_largerPriors.jpeg", site_df_largerPrior_plot, width = 8, height = 6, units = "in", dpi = 300)

asppCombined <- (asppEstimates_plot) /
  (aspp_df_largerPrior_plot) 
asppCombined
ggsave("figures/gddLeafout_empData/asppCombined.jpeg", asppCombined, width = 6, height = 8, units = "in", dpi = 300)

sitesCombined <- (siteEstimates_plot) /
  (site_df_largerPrior_plot) 
sitesCombined
ggsave("figures/gddLeafout_empData/sitesCombined.jpeg", sitesCombined, width = 6, height = 8, units = "in", dpi = 300)




treeid_df_largerPriors2$treeid <- as.numeric(treeid_df_largerPriors2$treeid)
treeid_df_largerPriors2$treeid_name <- no_naleafout$treeid[match(treeid_df_largerPriors2$treeid, no_naleafout$treeid_num)]

# now do the same, but for sppecies
treeid_df_largerPriors2$spp <- no_naleafout$spp[match(treeid_df_largerPriors2$treeid, no_naleafout$treeid_num)]

# same for site
treeid_df_largerPriors2$site <- no_naleafout$site[match(treeid_df_largerPriors2$treeid, no_naleafout$treeid_num)]


# quick check that I didn't mess anything up
un <- no_naleafout[!duplicated(no_naleafout$treeid),]
table(un$spp)
table(treeid_df_largerPriors2$spp)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# sum atreeid, aspp, asite
aspp_df3_largerPriors <- as.matrix(aspp_df_largerPriors)
site_df3_largerPriors <- as.matrix(site_df_largerPriors)
treeid_df3_largerPriors <- as.matrix(treeid_df_largerPriors)

fullatreeid_largerPriors <-
  treeid_df3_largerPriors[, sub$treeid_num] +
  aspp_df3_largerPriors[, sub$spp_num] +
  site_df3_largerPriors[, sub$site_num] 
fullatreeid_largerPriors

fullatreeid_largerPriors2 <- as.data.frame(fullatreeid_largerPriors)
# get posterior means and quantiles

# empty treeid dataframe
treeid_df4_largerPriors <- data.frame(
  treeid = character(ncol(fullatreeid_largerPriors2)),
  fit_atreeid = numeric(ncol(fullatreeid_largerPriors2)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(fullatreeid_largerPriors2)) { # i = 1
  treeid_df4_largerPriors$treeid[i] <- colnames(fullatreeid_largerPriors2)[i]         
  treeid_df4_largerPriors$fit_atreeid[i] <- round(mean(fullatreeid_largerPriors2[[i]]),3)  
  treeid_df4_largerPriors$fit_atreeid_per5[i] <- round(quantile(fullatreeid_largerPriors2[[i]], probs = 0.05), 3)
  treeid_df4_largerPriors$fit_atreeid_per25[i] <- round(quantile(fullatreeid_largerPriors2[[i]], probs = 0.25), 3)
  treeid_df4_largerPriors$fit_atreeid_per75[i] <- round(quantile(fullatreeid_largerPriors2[[i]], probs = 0.75), 3)
  treeid_df4_largerPriors$fit_atreeid_per95[i] <- round(quantile(fullatreeid_largerPriors2[[i]], probs = 0.95), 3)
}
treeid_df4_largerPriors

# get the og treeid names, spp and site back:
treeid_df4_largerPriors$treeid <- as.numeric(treeid_df4_largerPriors$treeid)
treeid_df4_largerPriors$treeid_name <- no_naleafout$treeid[match(treeid_df4_largerPriors$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4_largerPriors$spp_name <- no_naleafout$spp[match(treeid_df4_largerPriors$treeid,
                                              no_naleafout$treeid_num)]
treeid_df4_largerPriors$site_name <- no_naleafout$site[match(treeid_df4_largerPriors$treeid,
                                                no_naleafout$treeid_num)]
treeid_df4_largerPriors

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Prior vs posterior ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### asp #####
# Regular prior
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
aspp_prior <- rnorm(1e4, 0, 1)*20

asp_priorunchanged <- ggplot() +
  geom_density(data = data.frame(aspp_prior = aspp_prior),
               aes(x = aspp_prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = aspp_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_asp back 
       converted to original scale",
       x = "aspp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  xlim(-200, 200)+
  theme_minimal()
asp_priorunchanged

# asp prior 3x
aspp_long_LP <- reshape(
  aspp_df_largerPriors,
  direction = "long",
  varying = list(names(aspp_df_largerPriors)),
  v.names = "value",
  timevar = "spp",
  times = names(aspp_df_largerPriors),
  idvar = "draw"
)
aspp_long_LP

# aspp prior
aspp_prior_LP <- rnorm(1e4, 0, 3)*20

asp_priorunchanged_LP <- ggplot() +
  geom_density(data = data.frame(aspp_prior_LP = aspp_prior_LP),
               aes(x = aspp_prior_LP, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = aspp_long_LP,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_asp back 
       converted to original scale with 3x priors",
       x = "aspp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  xlim(-200, 200) +
  theme_minimal()
asp_priorunchanged_LP

asp_priorPosterior <- (asp_priorunchanged)/
  asp_priorunchanged_LP
ggsave("figures/gddLeafout_empData/asp_priorPosterior.jpeg", asp_priorPosterior, width = 6, height = 10, units = "in", dpi = 300)

##### asite #####
# Regular prior
asite_long <- reshape(
  site_df,
  direction = "long",
  varying = list(names(site_df)),
  v.names = "value",
  timevar = "spp",
  times = names(site_df),
  idvar = "draw"
)
asite_long

# asite prior
asite_prior <- rnorm(1e4, 0, 1)*20

asite_priorunchanged <- ggplot() +
  geom_density(data = data.frame(asite_prior = asite_prior),
               aes(x = asite_prior, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = asite_long,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_asp back 
       converted to original scale",
       x = "asite", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  xlim(-200, 200)+
  theme_minimal()
asite_priorunchanged

# asite prior 3x
asite_long_LP <- reshape(
  site_df_largerPriors,
  direction = "long",
  varying = list(names(site_df_largerPriors)),
  v.names = "value",
  timevar = "spp",
  times = names(site_df_largerPriors),
  idvar = "draw"
)
asite_long_LP

# asite prior
asite_prior_LP <- rnorm(1e4, 0, 3)*20

asite_priorunchanged_LP <- ggplot() +
  geom_density(data = data.frame(asite_prior_LP = asite_prior_LP),
               aes(x = asite_prior_LP, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = asite_long_LP,
               aes(x = value, colour = "Posterior", group = spp),
               linewidth = 0.5) +
  # facet_wrap(~spp) + 
  labs(title = "priorVSposterior_asp back 
       converted to original scale with 3x priors",
       x = "asite", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  xlim(-200, 200) +
  theme_minimal()
asite_priorunchanged_LP

asite_priorPosterior <- (asite_priorunchanged)/
  asite_priorunchanged_LP
ggsave("figures/gddLeafout_empData/asite_priorPosterior.jpeg", asite_priorPosterior, width = 6, height = 10, units = "in", dpi = 300)

# Prior differences mean comparisons ####
colnames(treeid_df_largerPriors2)[2:ncol(treeid_df_largerPriors2)] <-
  paste(colnames(treeid_df_largerPriors2)[2:ncol(treeid_df_largerPriors2)], "LP", sep = "_")

treeidtoplot <- merge(treeid_df2, treeid_df_largerPriors2,
  by = "treeid"
)

treeidtoplot
# plot treeid
atreeid_plot_priorComp <- ggplot(treeidtoplot, aes(x = fit_atreeid, y = fit_atreeid_LP)) +
  geom_errorbar(aes(xmin = fit_atreeid_per25, xmax = fit_atreeid_per75), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_errorbar(aes(ymin = fit_atreeid_per25_LP, ymax = fit_atreeid_per75_LP), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "priors unchanged", y = "priors 3x larger", title = "atreeid estimates") +
  theme_minimal()
atreeid_plot_priorComp

# aspp
aspp_df_largerPriors23 <- aspp_df_largerPriors2
colnames(aspp_df_largerPriors23)[2:ncol(aspp_df_largerPriors23)] <-
  paste(colnames(aspp_df_largerPriors23)[2:ncol(aspp_df_largerPriors23)], "LP", sep = "_")

aspptoplot <- merge(aspp_df2, aspp_df_largerPriors23,
                      by = "spp"
)
# plot aspp
aspp_plot_priorComp <- ggplot(aspptoplot, aes(x = fit_aspp, y = fit_aspp_LP)) +
  geom_errorbar(aes(xmin = fit_aspp_per25, xmax = fit_aspp_per75), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_errorbar(aes(ymin = fit_aspp_per25_LP, ymax = fit_aspp_per75_LP), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "priors unchanged", y = "priors 3x larger", title = "aspp estimates") +
  theme_minimal()
aspp_plot_priorComp

# asite
site_df_largerPriors2 <- site_df_largerPriors2
colnames(site_df_largerPriors2)[2:ncol(site_df_largerPriors2)] <-
  paste(colnames(site_df_largerPriors2)[2:ncol(site_df_largerPriors2)], "LP", sep = "_")

asitetoplot <- merge(site_df2, site_df_largerPriors2,
                   by = "site"
)

# plot asite
asite_plot_priorComp <- ggplot(asitetoplot, aes(x = fit_asite, y = fit_asite_LP)) +
  geom_errorbar(aes(xmin = fit_asite_per25, xmax = fit_asite_per75), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_errorbar(aes(ymin = fit_asite_per25_LP, ymax = fit_asite_per75_LP), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "priors unchanged", y = "priors 3x larger", title = "asite estimates") +
  theme_minimal()
asite_plot_priorComp

prior11plot <- (atreeid_plot_priorComp + aspp_plot_priorComp + asite_plot_priorComp)
prior11plot
ggsave("figures/gddLeafout_empData/priorComp11plots.jpeg", prior11plot, width = 10, height = 6, units = "in", dpi = 300)


}


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# open device
jpeg(
  filename = "figures/gddLeafout_empData/meanPlot_treeidBYspp.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
par(mar = c(
  4, 
  6, 
  4, 
  5))

# define a gap between species clusters
gap <- 3

# y positions
treeid_df4$y_pos <- NA
current_y <- 1

species_order <- c(
  "ALNINC", 
  "BETALL", 
  "BETPAP", 
  "BETPOP")
# site_order <- c("SH", "GR", "WM", "HF")
site_order <- c(
  "HF",
  "WM",
  "GR", 
  "SH")

# col
my_colors <- c(
  ALNINC = wes_palette("AsteroidCity1")[1],
  BETALL = wes_palette("AsteroidCity1")[2],
  BETPAP = wes_palette("AsteroidCity1")[3],
  BETPOP = wes_palette("AsteroidCity1")[4]
)
# shapes for sites
my_shapes <- c(
  HF = 19,
  WM = 18,
  GR = 15,
  SH = 17
)


treeid_df4$spp  <- factor(treeid_df4$spp, levels = species_order)
treeid_df4$site <- factor(treeid_df4$site, levels = site_order)

treeid_df4 <- treeid_df4[
  order(treeid_df4$spp, treeid_df4$site, treeid_df4$treeid),
]

treeid_df4$y_pos <- seq_len(nrow(treeid_df4))

for(spp in species_order){
  idx <- which(treeid_df4$spp == spp)
  n <- length(idx)
  
  # assign sequential positions for this species
  treeid_df4$y_pos[idx] <- current_y:(current_y + n - 1)
  
  # move cursor down with a gap before next species cluster
  current_y <- current_y + n + gap
}

treeid_df4$y_pos

# Set up empty plot
plot(
  NA, NA,
  # xlim = range(c(treeid_df4$fit_atreeid_per5-15,
  #                treeid_df4$fit_atreeid_per95+15)),
  xlim = range(c(-100, 100)),
  ylim = c(1, max(treeid_df4$y_pos) + 0.5),
  xlab = "Estimate of GDD change from the grand mean",
  ylab = "",
  yaxt = "n"  
)


# --- Add horizontal error bars (5–95%) ---
segments(
  x0 = treeid_df4$fit_atreeid_per5,
  x1 = treeid_df4$fit_atreeid_per95,
  y0 = treeid_df4$y_pos,
  col = adjustcolor(my_colors[treeid_df4$spp], alpha.f = 0.4),
  lwd = 1
)

# --- Add thicker horizontal error bars (25–75%) ---
segments(
  x0 = treeid_df4$fit_atreeid_per25,
  x1 = treeid_df4$fit_atreeid_per75,
  y0 = treeid_df4$y_pos,
  col = adjustcolor(my_colors[treeid_df4$spp], alpha.f = 0.4),
  lwd = 1.5
)

# --- Add the points ---
points(
  treeid_df4$fit_atreeid,
  treeid_df4$y_pos,
  cex = 0.8,
  pch = my_shapes[treeid_df4$site],
  col = adjustcolor(my_colors[treeid_df4$spp], alpha.f = 0.4)
)

# aspp_df2$spp <- aspp_df2$spp_name

spp_y <- tapply(treeid_df4$y_pos, treeid_df4$spp, mean)

aspp_df2$y_pos <- spp_y[aspp_df2$spp]


segments(
  x0 = aspp_df2$fit_aspp_per5,
  x1 = aspp_df2$fit_aspp_per95,
  y0 = aspp_df2$y_pos,
  col = adjustcolor(my_colors[aspp_df2$spp], alpha.f = 0.9),
  lwd = 2
)

segments(
  x0 = aspp_df2$fit_aspp_per25,
  x1 = aspp_df2$fit_aspp_per75,
  y0 = aspp_df2$y_pos,
  col = my_colors[aspp_df2$spp],
  lwd = 4
)
points(
  aspp_df2$fit_aspp,
  aspp_df2$y_pos,
  pch = 21,
  bg  = my_colors[aspp_df2$spp],
  col = "black",
  cex = 1.5
)


# --- Add vertical line at 0 ---
abline(v = 0, lty = 2)

# --- Add custom y-axis labels (reverse order if needed) ---
axis(
  side = 2,
  at = treeid_df4$y_pos,
  labels = treeid_df4$treeid_name,
  cex.axis = 0.5,
  las = 1
)
# spp mean
spp_y <- tapply(treeid_df4$y_pos, treeid_df4$spp, mean)
site_y <- tapply(treeid_df4$y_pos, treeid_df4$site, max)

## order species by mean y descending (top of plot first)
species_legend_order <- names(sort(spp_y, decreasing = TRUE))
site_legend_order <- names(sort(site_y, decreasing = FALSE))

## species legend (colors matched by name)
legend(
  x = max(treeid_df4$fit_atreeid_per95) - 6,
  y = max(treeid_df4$y_pos) - 2,
  legend = species_legend_order,
  col = my_colors[species_legend_order],    # index so colors match
  pch = 16,
  pt.cex = 1.2,
  title = "Species",
  bty = "n"
)

site_legend_order <- c("SH", "GR", "WM", "HF")
# site legend
legend(
  x = max(treeid_df4$fit_atreeid_per95) - 6,
  y = max(treeid_df4$y_pos) - 25,
  legend = site_legend_order,
  pch = my_shapes[site_legend_order],
  pt.cex = 1.2,
  title = "Sites",
  bty = "n"
)

dev.off()

# === === === ===  === === === ===  === === === ===  === === === ===  === === === ===
# === === === ===  === === === ===  === === === ===  === === === ===  === === === ===
# === === === ===  === === === ===  === === === ===  === === === ===  === === === ===
# LARGER priors below
# open device
jpeg(
  filename = "figures/gddLeafout_empData/meanPlot_treeidBYspp_largerPriors.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
par(mar = c(
  4, 
  6, 
  4, 
  5))

# define a gap between species clusters
gap <- 3

# y positions
treeid_df4_largerPriors$y_pos <- NA
current_y <- 1

species_order <- c(
  "ALNINC", 
  "BETALL", 
  "BETPAP", 
  "BETPOP")
# site_order <- c("SH", "GR", "WM", "HF")
site_order <- c(
  "HF",
  "WM",
  "GR", 
  "SH")

# col
my_colors <- c(
  ALNINC = wes_palette("AsteroidCity1")[1],
  BETALL = wes_palette("AsteroidCity1")[2],
  BETPAP = wes_palette("AsteroidCity1")[3],
  BETPOP = wes_palette("AsteroidCity1")[4]
)
# shapes for sites
my_shapes <- c(
  HF = 19,
  WM = 18,
  GR = 15,
  SH = 17
)


treeid_df4_largerPriors$spp  <- factor(treeid_df4_largerPriors$spp, levels = species_order)
treeid_df4_largerPriors$site <- factor(treeid_df4_largerPriors$site, levels = site_order)

treeid_df4_largerPriors <- treeid_df4_largerPriors[
  order(treeid_df4_largerPriors$spp, treeid_df4_largerPriors$site, treeid_df4_largerPriors$treeid),
]

treeid_df4_largerPriors$y_pos <- seq_len(nrow(treeid_df4_largerPriors))

for(spp in species_order){
  idx <- which(treeid_df4_largerPriors$spp == spp)
  n <- length(idx)
  
  # assign sequential positions for this species
  treeid_df4_largerPriors$y_pos[idx] <- current_y:(current_y + n - 1)
  
  # move cursor down with a gap before next species cluster
  current_y <- current_y + n + gap
}

treeid_df4_largerPriors$y_pos

# Set up empty plot
plot(
  NA, NA,
  # xlim = range(c(treeid_df4_largerPriors$fit_atreeid_per5-15,
  #                treeid_df4_largerPriors$fit_atreeid_per95+15)),
  xlim = range(c(-100, 100)),
  ylim = c(1, max(treeid_df4_largerPriors$y_pos) + 0.5),
  xlab = "Estimate of GDD change from the grand mean",
  ylab = "",
  yaxt = "n"  
)


# --- Add horizontal error bars (5–95%) ---
segments(
  x0 = treeid_df4_largerPriors$fit_atreeid_per5,
  x1 = treeid_df4_largerPriors$fit_atreeid_per95,
  y0 = treeid_df4_largerPriors$y_pos,
  col = adjustcolor(my_colors[treeid_df4_largerPriors$spp], alpha.f = 0.4),
  lwd = 1
)

# --- Add thicker horizontal error bars (25–75%) ---
segments(
  x0 = treeid_df4_largerPriors$fit_atreeid_per25,
  x1 = treeid_df4_largerPriors$fit_atreeid_per75,
  y0 = treeid_df4_largerPriors$y_pos,
  col = adjustcolor(my_colors[treeid_df4_largerPriors$spp], alpha.f = 0.4),
  lwd = 1.5
)

# --- Add the points ---
points(
  treeid_df4_largerPriors$fit_atreeid,
  treeid_df4_largerPriors$y_pos,
  cex = 0.8,
  pch = my_shapes[treeid_df4_largerPriors$site],
  col = adjustcolor(my_colors[treeid_df4_largerPriors$spp], alpha.f = 0.4)
)

# aspp_df_largerPriors2$spp <- aspp_df_largerPriors2$spp_name

spp_y <- tapply(treeid_df4_largerPriors$y_pos, treeid_df4_largerPriors$spp, mean)

aspp_df_largerPriors2$y_pos <- spp_y[aspp_df_largerPriors2$spp]


segments(
  x0 = aspp_df_largerPriors2$fit_aspp_per5,
  x1 = aspp_df_largerPriors2$fit_aspp_per95,
  y0 = aspp_df_largerPriors2$y_pos,
  col = adjustcolor(my_colors[aspp_df_largerPriors2$spp], alpha.f = 0.9),
  lwd = 2
)

segments(
  x0 = aspp_df_largerPriors2$fit_aspp_per25,
  x1 = aspp_df_largerPriors2$fit_aspp_per75,
  y0 = aspp_df_largerPriors2$y_pos,
  col = my_colors[aspp_df_largerPriors2$spp],
  lwd = 4
)
points(
  aspp_df_largerPriors2$fit_aspp,
  aspp_df_largerPriors2$y_pos,
  pch = 21,
  bg  = my_colors[aspp_df_largerPriors2$spp],
  col = "black",
  cex = 1.5
)


# --- Add vertical line at 0 ---
abline(v = 0, lty = 2)

# --- Add custom y-axis labels (reverse order if needed) ---
axis(
  side = 2,
  at = treeid_df4_largerPriors$y_pos,
  labels = treeid_df4_largerPriors$treeid_name,
  cex.axis = 0.5,
  las = 1
)
# spp mean
spp_y <- tapply(treeid_df4_largerPriors$y_pos, treeid_df4_largerPriors$spp, mean)
site_y <- tapply(treeid_df4_largerPriors$y_pos, treeid_df4_largerPriors$site, max)

## order species by mean y descending (top of plot first)
species_legend_order <- names(sort(spp_y, decreasing = TRUE))
site_legend_order <- names(sort(site_y, decreasing = FALSE))

## species legend (colors matched by name)
legend(
  x = max(treeid_df4_largerPriors$fit_atreeid_per95) - 6,
  y = max(treeid_df4_largerPriors$y_pos) - 2,
  legend = species_legend_order,
  col = my_colors[species_legend_order],    # index so colors match
  pch = 16,
  pt.cex = 1.2,
  title = "Species",
  bty = "n"
)

site_legend_order <- c("SH", "GR", "WM", "HF")
# site legend
legend(
  x = max(treeid_df4_largerPriors$fit_atreeid_per95) - 6,
  y = max(treeid_df4_largerPriors$y_pos) - 25,
  legend = site_legend_order,
  pch = my_shapes[site_legend_order],
  pt.cex = 1.2,
  title = "Sites",
  bty = "n"
)

dev.off()


  
# Mis stuff ####
aspp_df2$percentchange <- aspp_df2$fit_aspp/mean(df_fit[,"a"])
  
##### Leaf out vs GDD #####
aspp_df23 <- aspp_df2
colnames(aspp_df23) <- c("spp_num", "fit_aspp", "fit_aspp_per5", "fit_aspp_per25", "fit_aspp_per75",  "fit_aspp_per95", "spp_name"
)
asppchecks <- merge(no_naleafout, aspp_df23, by = "spp_num")
asppchecks$a_aspp <- 157 + asppchecks$fit_aspp + site_df2$fit_asite
ggplot(asppchecks, aes(x = leafout)) +
  # geom_errorbar(aes(xmin = fit_aspp_per25, xmax = fit_aspp_per75), 
  #               width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  # geom_errorbar(aes(ymin = fit_aspp_per25_LP, ymax = fit_aspp_per75_LP), 
  #               width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_histogram(aes(color = spp_name, fill = spp_name), binwidth = 4, alpha = 0.7) +
  facet_wrap(~spp_name, nrow = 4, ncol = 3) +
  geom_vline(aes(xintercept = a_aspp)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  # labs(x = "priors unchanged", y = "priors 3x larger", title = "aspp estimates") +
  theme_minimal()



# how many days species leaf out different
ag <- aggregate(leafout ~ spp + year, no_naleafout,  mean)

### check how many gdd around those days
vec <- c(min(ag$leafout) : max(ag$leafout))
period <- subset(gdd, doy %in% vec)

base_temp <- 10

# average across years
periodmean <- aggregate(GDD_10 ~ doy, period, FUN = mean)
periodmean$GDD_10_diff <- NA

# calculates gdd increment per day
for (i in 2:nrow(periodmean)) {
  periodmean$GDD_10_diff[i] <- periodmean$GDD_10[i] - periodmean$GDD_10[i - 1]
}

# averages gdd increment per day
meangddperiod <- mean(periodmean$GDD_10_diff, na.rm = TRUE)

(max(ag$leafout) - min(ag$leafout)) * meangddperiod 

(max(ag$leafout) - min(ag$leafout)) * meangddperiod / 2

quantile(rnorm(1e4, 0, 50), prob = 0.05)
hist(rnorm(1e4, 0, 8)/20)
hist(rnorm(1e4, 0, 8/20))
# across 3 years, for the period of doy 125 to doy 135, there is approximately 4 GDD per day, so my prior need to span at least 22 GDD 


ggplot(no_naleafout) +
  geom_point(aes(x = leafout, y = leafoutGDD)) +
  facet_wrap(~year)
plot(no_naleafout$leafout)

##### Empirical data distribution #####
ggplot(no_naleafout, aes(x = leafoutGDD)) +
  geom_histogram(binwidth = 20, alpha = 0.7) +
  facet_wrap(~spp, nrow = 2, ncol = 2) +
  # geom_vline(aes(xintercept = a_aspp)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  # labs(x = "priors unchanged", y = "priors 3x larger", title = "aspp estimates") +
  theme_classic()
ggsave("figures/gddLeafout_empData/histSppgddleafoutggplot.jpeg", width = 8, height = 6, units = "in", dpi = 300)

par(mfrow=c(1, 2))

leafout_wide <- reshape(
  no_naleafout[, c("year","treeid","spp","leafoutGDD")],
  idvar = c("year","treeid"),
  timevar = "spp",
  direction = "wide"
)
par(mfrow=c(1, 1))
util$plot_line_hists(data$x, data$y, -6, 7, 0.5, xlab="")

ainc <- leafout_wide$leafoutGDD.ALNINC
ball <- leafout_wide$leafoutGDD.BETALL
bpap <- leafout_wide$leafoutGDD.BETPAP
bpop <- leafout_wide$leafoutGDD.BETPOP

jpeg(
  filename = "figures/gddLeafout_empData/histSppgddleafoutmike.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
par(mfrow=c(2, 2))
util$plot_line_hist(ainc, 0, 750, 20, xlab="gdd at leafout", main="ALNINC")
util$plot_line_hist(ball, 0, 750, 20, xlab="gdd at leafout", main="BETALL")
util$plot_line_hist(bpap, 0, 750, 20, xlab="gdd at leafout", main="BETPAP")
util$plot_line_hist(bpop, 0, 750, 20, xlab="gdd at leafout", main="BETPOP")
dev.off()
