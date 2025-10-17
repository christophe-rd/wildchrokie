
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
library(ggplot2)
library(rstan)
library(shinystan)
library(wesanderson)

if(length(grep("christophe_rouleau-desrochers", getwd()) > 0)) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if(length(grep("lizzie", getwd())) > 0){
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
}

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
b <- 0.4
sigma_y <- 0.2
sigma_a_spp <- 0.3 # This is pretty low, but I guess you think your species are closely related and will be similar?
sigma_a_treeid <- 0.5
sigma_a_site <- 0.1
sigma_b_spp <- 0.25

n_site <- 4 # number of sites
n_spp <- 10 # number of species
n_perspp <- 10 # number of individuals per species
n_treeid <- n_perspp * n_spp * n_site # number of treeid
n_meas <- 5 # repeated measurements per id
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

# get slope values for each species
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
simcoef$gddcons <- simcoef$gdd - mean(simcoef$gdd)

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$a_site + 
  simcoef$a_spp + 
  simcoef$a_treeid + 
  simcoef$a + 
  (simcoef$b*simcoef$gddcons) + 
  (simcoef$b_spp*simcoef$gddcons)+
  simcoef$error

# prepare grouping factors
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
simcoef$treeid <- factor(simcoef$treeid)

# === === === === === #
##### Run models #####
# === === === === === #
y <- simcoef$ringwidth
N <- nrow(simcoef)
gdd <- simcoef$gddcons
Nspp <- length(unique(simcoef$spp))
Nsite <- length(unique(simcoef$site))
site <- as.numeric(as.character(simcoef$site))
species <- as.numeric(as.character(simcoef$spp))
treeid <- treeid
Ntreeid <- length(unique(treeid))
table(treeid)

fit <- rstan::stan("stan/twolevelhierint.stan", 
                    data=c("N","y","Nspp","species","Nsite", "site", "Ntreeid", "treeid", "gdd"),
                    iter=4000, chains=4, cores=4,
                    control = list(max_treedepth = 15))  

# summary(fit)$summary
# saveRDS(fit, "output/fit")
launch_shinystan(fit)


# === === === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)

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

sigma_df2$sim_sigma <- c(sigma_b_spp, sigma_a_spp, sigma_a_site, sigma_a_treeid, sigma_y)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover b spp ######
bspp_cols <- colnames(df_fit)[grepl("bsp", colnames(df_fit))]
# remove sigma_aspp for now
bspp_cols <- bspp_cols[2:length(bspp_cols)]

bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# change their names
colnames(bspp_df) <- sub("bsp\\[(\\d+)\\]", "\\1", colnames(bspp_df))
#empty spp df
bspp_df2 <- data.frame(
  spp = character(ncol(bspp_df)),
  fit_b_spp = numeric(ncol(bspp_df)),  
  fit_b_spp_per5 = NA, 
  fit_b_spp_per25 = NA,
  fit_b_spp_per75 = NA,
  fit_b_spp_per95 = NA
)
for (i in 1:ncol(bspp_df)) { # i = 1
  bspp_df2$spp[i] <- colnames(bspp_df)[i]         
  bspp_df2$fit_b_spp[i] <- round(mean(bspp_df[[i]]),3)  
  bspp_df2$fit_b_spp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
  bspp_df2$fit_b_spp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
  bspp_df2$fit_b_spp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
  bspp_df2$fit_b_spp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
}



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
aspp_cols <- colnames(df_fit)[grepl("asp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("zasp", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("sigma", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("asp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
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
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]
# remove sigma_asp for now
site_cols <- site_cols[2:length(site_cols)]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
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
###### Plot sigmas ######
sigma_simXfit_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_point(color = "#046C9A", size = 3) +
  geom_errorbar(aes(ymin = per5, ymax = per95),
                width = 0, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "fit vs sim sigmas") +
  theme_minimal()
sigma_simXfit_plot
ggsave("figures/sigma_simXfit_plot.jpeg", sigma_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot b spp ######
bspptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "b_spp")], 
  bspp_df2[!duplicated(bspp_df2$spp), 
           c("spp", "fit_b_spp", "fit_b_spp_per5", "fit_b_spp_per95")], 
  by = "spp"
)
bspptoplot

b_spp_simXfit_plot <- ggplot(bspptoplot, aes(x = b_spp, y = fit_b_spp)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_b_spp_per5, ymax = fit_b_spp_per95), width = 0, color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim b_spp", y = "fit b_spp", title = "") +
  theme_minimal()
b_spp_simXfit_plot
# ggsave!
ggsave("figures/b_spp_simXfit_plot2.jpeg", b_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot treeid ######
# add sim to fit treeid df
treeidtoplot <- merge(
  simcoef[!duplicated(simcoef$treeid), 
          c("treeid", "a_treeid")], 
  treeid_df2[!duplicated(treeid_df2$treeid), 
             c("treeid", "fit_a_treeid", "fit_a_treeid_per5", "fit_a_treeid_per95")], 
  by = "treeid"
)
treeidtoplot
# plot treeid
a_treeid_simXfit_plot <- ggplot(treeidtoplot, aes(x = a_treeid, y = fit_a_treeid)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_treeid_per5, ymax = fit_a_treeid_per95), width = 0, color = "darkgray", alpha=0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_treeid", y = "fit a_treeid", title = "") +
  theme_minimal()
a_treeid_simXfit_plot
# ggsave!
ggsave("figures/a_treeid_simXfit_plot.jpeg", a_treeid_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot a spp ######
aspptoplot <- merge(
  simcoef[!duplicated(simcoef$spp), 
          c("spp", "a_spp")], 
  aspp_df2[!duplicated(aspp_df2$spp), 
           c("spp", "fit_a_spp", "fit_a_spp_per5", "fit_a_spp_per95")], 
  by = "spp"
)
aspptoplot

a_spp_simXfit_plot <- ggplot(aspptoplot, aes(x = a_spp, y = fit_a_spp)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_spp_per5, ymax = fit_a_spp_per95), width = 0, color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_spp", y = "fit a_spp", title = "") +
  theme_minimal()
a_spp_simXfit_plot
# ggsave!
ggsave("figures/a_spp_simXfit_plot.jpeg", a_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Plot site ######
sitetoplot <- merge(
  simcoef[!duplicated(simcoef$site), 
          c("site", "a_site")], 
  site_df2[!duplicated(site_df2$site), 
           c("site", "fit_a_site", "fit_a_site_per5", "fit_a_site_per95")], 
  by = "site"
)
sitetoplot

a_site_simXfit_plot <- ggplot(sitetoplot, aes(x = a_site, y = fit_a_site)) +
  geom_point(color = "#046C9A", size = 2) +
  geom_errorbar(aes(ymin = fit_a_site_per5, ymax = fit_a_site_per95), width = 0, color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim a_site", y = "fit a_site", title = "") +
  theme_minimal()
a_site_simXfit_plot
# ggsave!
ggsave("figures/a_site_simXfit_plot.jpeg", a_site_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)


# === === === === === === === === === === === === === === === === 
#### Step 3. Look at my priors####
# === === === === === === === === === === === === === === === === 
n <- 1e4
prior_a <- rnorm(n, 2, 4)
hist(prior_a)
prior_b <- rnorm(n, 0.4, 1)
hist(prior_b)
prior_bsp <- rnorm(n, 0, 0.1)
hist(prior_bsp)
prior_zasp <-  rnorm(n, 0, 1)
hist(prior_zasp)
prior_asite <-  rnorm(n, 0, 0.5)
hist(prior_asite)
prior_sigma_bsp <-  rnorm(n, 0, 0.2)
hist(prior_sigma_bsp)
prior_sigma_asp <-  rnorm(n, 0, 0.5)
hist(prior_sigma_asp)
prior_sigma_asite <-  rnorm(n, 0, 0.06 )
hist(prior_sigma_asite)
prior_sigma_atreeid <-  rnorm(n, 0, 0.05)
hist(prior_sigma_atreeid)

##### Priors VS Posterior #####
# prior predictive checks. Simulating prior values from the values set in the model block
Ndraws <- 1000
colnames(sigma_df) <- paste("post", colnames(sigma_df), sep = "_")

###### sigmas ######
sigma_bsp_draw <- abs(rnorm(Ndraws, 0, 0.2))   
sigma_asp_draw <- abs(rnorm(Ndraws, 0, 0.5))
sigma_asite_draw <- abs(rnorm(Ndraws, 0, 0.5))
sigma_atree_draw <- abs(rnorm(Ndraws, 0, 0.05))
sigma_y_draw <- abs(rnorm(Ndraws, 0, 5))

inds <- 1:Ndraws

bsp_list <- vector("list", Ndraws)   
asp_list <- vector("list", Ndraws)  
asite_list <- vector("list", Ndraws)
atree_list <- vector("list", Ndraws)

for (i in seq_along(inds)) {
  idx <- inds[i]
  bsp_list[[i]]   <- rnorm(n_spp, 0, sigma_bsp_draw[idx])
  asp_list[[i]]   <- rnorm(n_spp, 0, sigma_asp_draw[idx]) 
  asite_list[[i]] <- rnorm(n_site, 0, sigma_asite_draw[idx])
  atree_list[[i]] <- rnorm(n_treeid, 0, sigma_atree_draw[idx])
}

# Flatten bsp/asp into data.frames for plotting
bsp_df <- do.call(rbind, lapply(1:Ndraws, function(i) {
  data.frame(draw = i, sigma_bsp = sigma_bsp_draw[inds[i]],
             species = 1:Nspp, bsp = bsp_list[[i]])
}))
asp_df <- do.call(rbind, lapply(1:Ndraws, function(i) {
  data.frame(draw = i, sigma_asp = sigma_asp_draw[inds[i]],
             species = 1:Nspp, asp = asp_list[[i]])
}))

ggplot(sigma_long_bsp, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()



# now add row for prior_bsp
prior_bsp <- rnorm(nrow(sigma_df), 0, sigma_df$prior_sigma_bsp)

# convert each parameter to long format
sigma_long_bsp <- data.frame(
  value  = c(sigma_df$post_sigma_bsp, sigma_df$prior_sigma_bsp),
  source = rep(c("post_sigma_bsp", "prior_sigma_bsp"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_bsp, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# Sigma asp
sigma_long_asp <- data.frame(
  value  = c(sigma_df$post_sigma_asp, sigma_df$prior_sigma_asp),
  source = rep(c("post_sigma_asp", "prior_sigma_asp"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_asp, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# Sigma asite
sigma_long_asite <- data.frame(
  value  = c(sigma_df$post_sigma_asite, sigma_df$prior_sigma_asite),
  source = rep(c("post_sigma_asite", "prior_sigma_asite"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_asite, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# Sigma atreeid
sigma_long_atreeid <- data.frame(
  value  = c(sigma_df$post_sigma_atreeid, sigma_df$prior_sigma_atreeid),
  source = rep(c("post_sigma_atreeid", "prior_sigma_atreeid"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_atreeid, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

# sigma_y
sigma_long_atreeid <- data.frame(
  value  = c(sigma_df$post_sigma_atreeid, sigma_df$prior_sigma_y),
  source = rep(c("post_sigma_atreeid", "prior_sigma_atreeid"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_atreeid, aes(x = value, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")[3:4])+
  theme_minimal()

#  === === === === === === === === === === === === === === === === 
#### Diagnostics ####
# === === === === === === === === === === === === === === === === 

# Non-centered parameterization #####

# bspp
# y: log sigma for all my partial pooled parameters 
# x: sigma of every parameter (e.g. spp1, spp2), need to check all of them individually.
# tree id 

sigma_df
sim_df
# aspp

# site

# === === === === === === === === === === === === === === === === 
#### Run model on empirical data ####
# === === === === === === === === === === === === === === === === 
