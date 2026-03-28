# Wildchrokie model
# CRD 27 February 2026
# We want to figure out if fitting the nested model returns the same thing

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)

# Load library 
library(ggplot2)
library(rstan)
library(wesanderson)
library(patchwork)

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
source('rcode/utilExtractParam.R')

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Simulate data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
set.seed(124)

Nspp <- 20 
N_tree_per_spp <- 15
N_obs_per_tree <- 12
Ntreeid <- Nspp*N_tree_per_spp
N <- Ntreeid * N_obs_per_tree

tree_species_idxs <- rep(1:Nspp, times = N_tree_per_spp)

# index
treeid <- rep(1:Ntreeid, each = N_obs_per_tree)
species <- tree_species_idxs[treeid]

sigma_y <- 3
sigma_atreeid <- 1.5
sigma_aspp <- 2

aspp <- rnorm(Nspp, 0, sigma_aspp)

# Nested 
alpha_tree_species <- rnorm(Ntreeid, aspp[tree_species_idxs], sigma_atreeid)
alpha_tree <- alpha_tree_species - aspp[tree_species_idxs]

y <- rnorm(N, alpha_tree_species[treeid], sigma_y)

fitnested <- stan("stan/twolevelhierint_nested.stan", 
                  data=c("N","y",
                         "Nspp","species",
                         "Ntreeid", "treeid"),
                  warmup = 1000, iter = 2000, chains=4)

# Non-nested
atreeid_nonested <- rnorm(Ntreeid, 0, sigma_atreeid)
atreeid_aspp_nonnested <- aspp[tree_species_idxs] + atreeid_nonested

y <- rnorm(N, atreeid_aspp_nonnested[treeid], sigma_y)

fitnonnested <- stan("stan/twolevelhierint_nested.stan", 
                  data=c("N","y",
                         "Nspp","species",
                         "Ntreeid", "treeid"),
                  warmup = 1000, iter = 2000, chains=4)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Recover and plot parameters ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Nested #####
df_fit <- as.data.frame(fitnested)

# posterior summaries
sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
treeid_df2 <- extract_params(df_fit, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")

# add sim coef
treeid_df2$sim_treeid <- alpha_tree
aspp_df2$sim_aspp <- aspp


# Plot treeid 
treeidnested <- ggplot(treeid_df2, aes(x = sim_treeid, y = fit_atreeid)) +
  geom_errorbar(aes(ymin = fit_atreeid_per5, ymax = fit_atreeid_per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 1, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "atreeid nested") +
  theme_minimal()

# Plot aspp 
asppnested <- ggplot(aspp_df2, aes(x = sim_aspp, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per5, ymax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim aspp", y = "fit aspp", title = "aspp nested") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()

##### Non-nested #####
df_fitnonnested <- as.data.frame(fitnonnested)

# posterior summaries
sigma_df2_nonnested  <- extract_params(df_fitnonnested, "sigma", "mean", "sigma")
treeid_df2_nonnested <- extract_params(df_fitnonnested, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_nonnested <- subset(treeid_df2_nonnested, !grepl("z|sigma", treeid))
aspp_df2_nonnested   <- extract_params(df_fitnonnested, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")

# add sim coef
treeid_df2_nonnested$sim_treeid <- atreeid_nonested
aspp_df2_nonnested$sim_aspp <- aspp


# Plot treeid 
treeidnonnested <- ggplot(treeid_df2_nonnested, aes(x = sim_treeid, y = fit_atreeid)) +
  geom_errorbar(aes(ymin = fit_atreeid_per5, ymax = fit_atreeid_per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 1, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "atreeid non-nested") +
  theme_minimal()

# Plot aspp 
asppnonnested <- ggplot(aspp_df2_nonnested, aes(x = sim_aspp, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per5, ymax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim aspp", y = "fit aspp", title = "aspp non-nested") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()

combined_plot <- (treeidnested + asppnested) /
(treeidnonnested + asppnonnested)
combined_plot
ggsave("figures/simdata/combinedPrmRecoveryNested.jpeg", combined_plot, 
       width = 8, height = 8, units = "in", dpi = 300)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
samples <- util$extract_expectand_vals(fit)
jpeg(
  filename = "figures/modelGrowthGDD/retrodictiveCheckHist.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
util$plot_hist_quantiles(samples, "y_rep", 
                         -5, # lower x axis limit
                         15, # upper x axis limit
                         0.5, # binning
                         baseline_values = y,
                         xlab = "Ring width (mm)")
dev.off()
