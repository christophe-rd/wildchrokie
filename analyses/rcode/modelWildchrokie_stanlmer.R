# Wildchrokie model
# CRD 29 October  2025
# temporary model fit for presentation purpose 

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)
# quartz()

# Load library 
library(ggplot2)
library(rstanarm)
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

# sim data
a <- 1.5
b <- 0.4
sigma_y <- 0.1
sigma_a_spp <- 0.3
sigma_a_treeid <- 0.15
sigma_a_site <- 0.3
sigma_b_spp <- 0.25

n_site <- 5 # number of sites
n_spp <- 20 # number of species
n_perspp <- 5 # number of individuals per species
n_treeid <- n_perspp * n_spp * n_site # number of treeid
n_meas <- 5 # repeated measurements per id
N <- n_treeid * n_meas # total number of measurements
N

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
# quick check 
table(treeidnonrep, site_nonrep)

simcoef <- data.frame(
  site = site,
  spp = spp,
  treeid = treeid
)

# get intercept values for each species
a_spp <- rnorm(n_spp, 0, sigma_a_spp)
a_site <- rnorm(n_site, 0, sigma_a_site)
a_treeid <- rnorm(n_treeid, 0, sigma_a_treeid)

# get slope values for each speciess
b_spp <- rnorm(n_spp, 0, sigma_b_spp)

# Add my parameters to the df
simcoef$a_treeid <- a_treeid[treeid]
simcoef$a_site <- a_site[simcoef$site]
simcoef$a_spp <- a_spp[simcoef$spp]
simcoef$b_spp <- b_spp[simcoef$spp]

# add the rest of the boring stuff 
simcoef$a <- a
simcoef$b <- b
simcoef$sigma_y <- sigma_y
simcoef$sigma_a_treeid <- sigma_a_treeid
simcoef$sigma_a_spp <- sigma_a_spp
simcoef$sigma_a_site <- sigma_a_site
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gdd <- rnorm(N, 1800, 100)
simcoef$gddcons <- simcoef$gdd/200

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$a_site + 
  simcoef$a_spp + 
  simcoef$a_treeid + 
  # simcoef$a + 
  (simcoef$b*simcoef$gddcons) + 
  # (simcoef$b_spp*simcoef$gddcons)+
  simcoef$error

# prepare grouping factors
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
simcoef$treeid <- factor(simcoef$treeid)

# run model
fit <- stan_lmer(
  ringwidth ~ 1 + gddcons + 
    (1|site) + (1|spp) + (1|treeid),
  data = simcoef,
  chains = 4,
  iter = 4000,
  core=4
)

df_fit <- as.data.frame(fit)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover b spp ######
# bspp_cols <- colnames(df_fit)[grepl("bsp", colnames(df_fit))]
# # remove sigma_aspp for now
# bspp_cols <- bspp_cols[2:length(bspp_cols)]
# 
# bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# # change their names
# colnames(bspp_df) <- sub("bsp\\[(\\d+)\\]", "\\1", colnames(bspp_df))
# #empty spp df
# bspp_df2 <- data.frame(
#   spp = character(ncol(bspp_df)),
#   fit_b_spp = numeric(ncol(bspp_df)),  
#   fit_b_spp_per5 = NA, 
#   fit_b_spp_per25 = NA,
#   fit_b_spp_per75 = NA,
#   fit_b_spp_per95 = NA
# )
# for (i in 1:ncol(bspp_df)) { # i = 1
#   bspp_df2$spp[i] <- colnames(bspp_df)[i]         
#   bspp_df2$fit_b_spp[i] <- round(mean(bspp_df[[i]]),3)  
#   bspp_df2$fit_b_spp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
#   bspp_df2$fit_b_spp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
#   bspp_df2$fit_b_spp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
#   bspp_df2$fit_b_spp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
# }



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover treeid ######

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("treeid", colnames(df_fit))]
# remove sigma_asp for now
treeid_cols <- treeid_cols[1:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub(".*treeid:(\\d+).*", "\\1", colnames(treeid_df))

# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_a_treeid = numeric(ncol(treeid_df)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_a_treeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_a_treeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover a spp  ######
aspp_cols <- colnames(df_fit)[grepl("spp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("Sigma", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub(".*spp:(\\d+).*", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_a_spp = numeric(ncol(aspp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per25 = NA,
  fit_a_spp_per75 = NA,
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_a_spp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_a_spp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_a_spp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_a_spp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_a_spp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover site ######
site_cols <- colnames(df_fit)[grepl("site", colnames(df_fit))]
# remove for now
site_cols <- site_cols[1:length(site_cols)]
site_cols <- site_cols[!grepl("Sigma", site_cols)]
site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub(".*site:(\\d+).*", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per25 = NA,
  fit_a_site_per75 = NA,
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2

# === === === === === === === #
# Plot parameter recovery #####
# === === === === === === === #

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot b spp ######
# bspptoplot <- merge(
#   simcoef[!duplicated(simcoef$spp), 
#           c("spp", "b_spp")], 
#   bspp_df2[!duplicated(bspp_df2$spp), 
#            c("spp", "fit_b_spp", "fit_b_spp_per25", "fit_b_spp_per75", "fit_b_spp_per5", "fit_b_spp_per95")], 
#   by = "spp"
# )
# bspptoplot
# 
# b_spp_simXfit_plot <- ggplot(bspptoplot, aes(x = b_spp, y = fit_b_spp)) +
#   geom_errorbar(aes(ymin = fit_b_spp_per5, ymax = fit_b_spp_per95), 
#                 width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
#   geom_errorbar(aes(ymin = fit_b_spp_per25, ymax = fit_b_spp_per75), 
#                 width = 0, linewidth = 1.5, color = "darkgray", alpha = 0.7) +
#   geom_point(color = "#046C9A", size = 2) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
#   labs(x = "sim b_spp", y = "fit b_spp", title = "") +
#   theme_minimal()
# b_spp_simXfit_plot
# # ggsave!
# ggsave("figures/lmer_b_spp_simXfit_plot2.jpeg", b_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
# add sim to fit treeid df
treeidtoplot <- merge(
  simcoef[!duplicated(simcoef$treeid), 
          c("treeid", "a_treeid")], 
  treeid_df2[!duplicated(treeid_df2$treeid), 
             c("treeid", "fit_a_treeid", 
               "fit_a_treeid_per5", 
               "fit_a_treeid_per25",
               "fit_a_treeid_per75",
               "fit_a_treeid_per95")], 
  by = "treeid"
)
treeidtoplot
# plot treeid
a_treeid_simXfit_plot <- ggplot(treeidtoplot, aes(x = a_treeid, y = fit_a_treeid)) +
  geom_errorbar(aes(ymin = fit_a_treeid_per5, ymax = fit_a_treeid_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_treeid", y = "fit a_treeid", title = "") +
  theme_minimal()
a_treeid_simXfit_plot
# ggsave!
ggsave("figures/lmer_a_treeid_simXfit_plot.jpeg", a_treeid_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot a spp ######
aspptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "a_spp")], 
  aspp_df2[!duplicated(aspp_df2$spp), 
           c("spp", "fit_a_spp", 
             "fit_a_spp_per5", 
             "fit_a_spp_per25", 
             "fit_a_spp_per75", 
             "fit_a_spp_per95")], 
  by = "spp"
)
aspptoplot

a_spp_simXfit_plot <- ggplot(aspptoplot, aes(x = a_spp, y = fit_a_spp)) +
  geom_errorbar(aes(ymin = fit_a_spp_per5, ymax = fit_a_spp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_a_spp_per25, ymax = fit_a_spp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_spp", y = "fit a_spp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
a_spp_simXfit_plot
# ggsave!
ggsave("figures/lmer_a_spp_simXfit_plot.jpeg", a_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######
sitetoplot <- merge(
  simcoef[!duplicated(simcoef$site), 
          c("site", "a_site")], 
  site_df2[!duplicated(site_df2$site), 
           c("site", "fit_a_site", 
             "fit_a_site_per5", 
             "fit_a_site_per25", 
             "fit_a_site_per75", 
             "fit_a_site_per95")], 
  by = "site"
)
sitetoplot

a_site_simXfit_plot <- ggplot(sitetoplot, aes(x = a_site, y = fit_a_site)) +
  geom_errorbar(aes(ymin = fit_a_site_per5, ymax = fit_a_site_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_a_site_per25, ymax = fit_a_site_per75), 
                width = 0, linewidth = 1.5, color = "darkgray", alpha=0.9) +
  geom_point(color = "#046C9A", size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_site", y = "fit a_site", title = "") +
  theme_minimal()
a_site_simXfit_plot
# ggsave!
ggsave("figures/lmer_a_site_simXfit_plot.jpeg", a_site_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)



