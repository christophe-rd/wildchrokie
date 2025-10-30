# Wildchrokie model
# CRD 30 October 2025
# removing tree id to trouble shoot the model

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
sigma_y <- 0.1
sigma_a_spp <- 0.5 # This is pretty low, but I guess you think your species are closely related and will be similar?
sigma_a_treeid <- 0.15
sigma_a_site <- 0.3
sigma_b_spp <- 0.25

n_site <- 10 # number of sites
n_spp <- 10 # number of species
n_perspp <- 10 # number of individuals per species
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
# a_treeid <- rnorm(n_treeid, 0, sigma_a_treeid)

# get slope values for each speciess
b_spp <- rnorm(n_spp, 0, sigma_b_spp)

# Add my parameters to the df
# simcoef$a_treeid <- a_treeid[treeid]
simcoef$a_site <- a_site[simcoef$site]
simcoef$a_spp <- a_spp[simcoef$spp]
simcoef$b_spp <- b_spp[simcoef$spp]

# add the rest of the boring stuff 
simcoef$a <- a
simcoef$b <- b
simcoef$sigma_y <- sigma_y
# simcoef$sigma_a_treeid <- sigma_a_treeid
simcoef$sigma_a_spp <- sigma_a_spp
simcoef$sigma_a_site <- sigma_a_site
simcoef$error <- rnorm(N, 0, sigma_y)
simcoef$gdd <- rnorm(N, 1800, 100)
simcoef$gddcons <- simcoef$gdd/200

# adding both options of tree rings
simcoef$ringwidth <- 
  simcoef$a_site + 
  simcoef$a_spp + 
  # simcoef$a_treeid + 
  simcoef$a +
  (simcoef$b*simcoef$gddcons) + 
  (simcoef$b_spp*simcoef$gddcons)+
  simcoef$error

# prepare grouping factors
simcoef$site <- factor(simcoef$site) 
simcoef$spp <- factor(simcoef$spp)
# simcoef$treeid <- factor(simcoef$treeid)

# === === === === === #
##### Run model #####
# === === === === === #
y <- simcoef$ringwidth
N <- nrow(simcoef)
gdd <- simcoef$gddcons
Nspp <- length(unique(simcoef$spp))
Nsite <- length(unique(simcoef$site))
site <- as.numeric(as.character(simcoef$site))
species <- as.numeric(as.character(simcoef$spp))
# treeid <- treeid
# Ntreeid <- length(unique(treeid))
table(treeid)

rstan_options(auto_write = TRUE)

fit <- stan("stan/twolevelhierint_notreeid.stan", 
                    data=c("N","y","Nspp","species","Nsite", "site", 
                           # "Ntreeid", "treeid", 
                           "gdd"),
                    iter=4000, chains=4, cores=4)

saveRDS(fit, "output/stanOutput/fit_no_a")
fit <- readRDS("output/stanOutput/fit_withbspp")
run_fit_noSite <- FALSE

if (FALSE) {

# === === === === === === === === === === === === #
##### Recover parameters from the posterior #####
# === === === === === === === === === === === === #
df_fit <- as.data.frame(fit)


if (FALSE){
  write.csv(df_fit, "output/df_fit.csv")
  
  df_fit <- read.csv("output/df_fit.csv")
  
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

sigma_df2$sim_sigma <- c(
  sigma_b_spp,
                         sigma_a_spp, 
                         sigma_a_site, 
                         sigma_a_treeid, 
                         sigma_y)


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
# temporary removal of treeid
sigma_simXfit_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_errorbar(aes(ymin = per25, ymax = per75),
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 1) +
  geom_errorbar(aes(ymin = per5, ymax = per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  geom_point(color = "#046C9A", size = 3) +
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
           c("spp", "fit_b_spp", "fit_b_spp_per25", "fit_b_spp_per75", "fit_b_spp_per5", "fit_b_spp_per95")], 
  by = "spp"
)
bspptoplot

b_spp_simXfit_plot <- ggplot(bspptoplot, aes(x = b_spp, y = fit_b_spp)) +
  geom_errorbar(aes(ymin = fit_b_spp_per5, ymax = fit_b_spp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_errorbar(aes(ymin = fit_b_spp_per25, ymax = fit_b_spp_per75), 
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 0.7) +
  geom_point(color = "#046C9A", size = 2) +
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
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.1) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.1) +
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
ggsave("figures/a_spp_simXfit_plot.jpeg", a_spp_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)

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
ggsave("figures/a_site_simXfit_plot.jpeg", a_site_simXfit_plot, width = 6, height = 6, units = "in", dpi = 300)


# === === === === === === === === === === === === === === === === 
#### Look at my priors####
# === === === === === === === === === === === === === === === === 

##### Priors VS Posterior #####
# prior predictive checks. Simulating prior values from the values set in the model block
draws <- 8000

vec <- rep(NA, draws*length(sigma_b_spp))

###### sigmas ######
sigma_bsp_draw <- abs(rnorm(Ndraws, 0, 0.3))   
sigma_asp_draw <- abs(rnorm(Ndraws, 0, 0.5))
sigma_asite_draw <- abs(rnorm(Ndraws, 0, 0.5))
sigma_atree_draw <- abs(rnorm(Ndraws, 0, 0.05))
sigma_y_draw <- abs(rnorm(Ndraws, 0, 5))

inds <- 1:Ndraws
bsp_list <- vector("list", Ndraws)   
asp_list <- vector("list", Ndraws)  
asite_list <- vector("list", Ndraws)
atree_list <- vector("list", Ndraws)

for (i in seq_along(inds)) { # i =1
  idx <- inds[i]
  bsp_list[[i]]   <- rnorm(n_spp, 0, sigma_bsp_draw[idx])
  asp_list[[i]]   <- rnorm(n_spp, 0, sigma_asp_draw[idx]) 
  asite_list[[i]] <- rnorm(n_site, 0, sigma_asite_draw[idx])
  atree_list[[i]] <- rnorm(n_treeid, 0, sigma_atree_draw[idx])
}


for (i in seq_along(inds)) { # i =1
  idx <- inds[i]
  bsp_list[[i]]   <- rnorm(n_spp, 0, sigma_bsp_draw[idx])
}
  
# Flatten bsp/asp into data.frames for plotting
prior_bsp_df <- do.call(rbind, lapply(1:Ndraws, function(i) {
  data.frame(draw = i, sigma_bsp = sigma_bsp_draw[inds[i]],
             spp = 1:Nspp, prior_bsp = bsp_list[[i]])
}))

# copy of bspp_df
bspp_df3 <- bspp_df

bspp_df3$draw <- rownames(bspp_df3)

colnames(bspp_df3)[1:20] <- paste0("spp", 1:20)

long_post_bspp <- reshape(
  bspp_df3,
  direction = "long",
  varying = paste0("spp", 1:20),
  v.names = "post_bsp",
  idvar = "draw",
  timevar = "spp"
)

# join in the posterior estimates
joined_bsp <- cbind(prior_bsp_df, long_post_bspp)

# get a subset of 10000 randomly selected rows
sub <- joined_bsp[sample(0:160000, 1e4, replace = TRUE), 
                  c(1,3,4,7)
                  ]
sub

ggplot(sub) +
  geom_density(aes(x = prior_bsp, colour = "Prior", group = spp),
               linewidth = 0.3) +
  geom_density(aes(x = post_bsp, colour = "Posterior", group = spp),
               linewidth = 0.3) +
  labs(title = "priorVSposterior_bsp", x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
ggsave("figures/priorVSposterior_bsp.jpeg", width = 8, height = 6, units = "in", dpi = 300)


asp_df <- do.call(rbind, lapply(1:Ndraws, function(i) {
  data.frame(draw = i, sigma_asp = sigma_asp_draw[inds[i]],
             species = 1:Nspp, asp = asp_list[[i]])
}))

sigma_bsp_draw <- rnorm(Ndraws, 0, 0.3)
sigma_df$prior_sigma_bsp <- sigma_bsp_draw

ggplot(sigma_df) +
  geom_density(aes(x = prior_sigma_bsp, colour = "Prior"),
               linewidth = 0.3) +
  geom_density(aes(x = post_sigma_bsp, colour = "Posterior"),
               linewidth = 0.3) +
  labs(title = "priorVSposterior_sigma_bsp", x = "sigma_bsp", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
ggsave("figures/priorVSposterior_sigma_bsp.jpeg", width = 8, height = 6, units = "in", dpi = 300)




# now add row for prior_bsp
prior_bsp <- rnorm(nrow(sigma_df), 0, sigma_df$prior_sigma_bsp)

# convert each parameter to long format
sigma_long_bsp <- data.frame(
  value  = c(sigma_df$post_sigma_bsp, sigma_df$prior_sigma_bsp),
  source = rep(c("post_sigma_bsp", "prior_sigma_bsp"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_bsp, aes(x = value, color = source)) +
  geom_density(alpha = 0.3) +
  labs(color = "Parameter", fill = "Parameter") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
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

# get asp df and add prefix to spp
colnames(aspp_df) <- paste("asp", colnames(aspp_df), sep = "")

sigmaXasp <- cbind(sigma_df, aspp_df)

plot(log(sigmaXasp$sigma_asp) ~ sigmaXasp$asp1)

asp1 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp1, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp2 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp2, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp3 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp3, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp4 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp4, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp5 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp5, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp6 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp6, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp7 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp7, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp8 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp8, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp9 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp9, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()

asp10 <- ggplot(sigmaXasp) + 
  geom_point(aes(asp10, log(sigma_asp)), alpha = 0.03, color = "#B40F20") + 
  theme_minimal()


combined_plot <- (asp1 + asp2 + asp3 + asp4 + asp5) /
  (asp6 + asp7 + asp8 + asp9 + asp10)

combined_plot

# save combined plot
ggsave("figures/asp_parameterization.jpeg", combined_plot, width = 12, height = 6, units = "in", dpi = 300)
ggsave("figures/asp_parameterization.jpeg", combined_plot,
       width = 12, height = 6, units = "in", dpi = 300, device = "jpeg")

# site

# get asp df and add prefix to spp
colnames(site_df) <- paste("site", colnames(site_df), sep = "")

sigmaXasite <- cbind(sigma_df, aspp_df)

site1 <- ggplot(sigmaXasite) + 
  geom_point(aes(site1, log(sigma_asite)), alpha = 0.03, color = "#273046") + 
  theme_minimal()

site2 <- ggplot(sigmaXasite) + 
  geom_point(aes(site2, log(sigma_asite)), alpha = 0.03, color = "#273046") + 
  theme_minimal()

site3 <- ggplot(sigmaXasite) + 
  geom_point(aes(site3, log(sigma_asite)), alpha = 0.03, color = "#273046") + 
  theme_minimal()

site4 <- ggplot(sigmaXasite) + 
  geom_point(aes(site4, log(sigma_asite)), alpha = 0.03, color = "#273046") + 
  theme_minimal()

combined_plot <- (site1 + site2 + site3 + site4)

combined_plot

# save combined plot
ggsave("figures/asite_parameterization.jpeg", combined_plot,
       width = 12, height = 6, units = "in", dpi = 300, device = "jpeg")
}
# === === === === === === === === === === === === === === === === 
#### Run model on empirical data ####
# === === === === === === === === === === === === === === === === 

