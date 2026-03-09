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

# set parameters
n_spp <- 15 # number of species
n_perspp <- 10 # number of individuals per species
n_treeid <- n_perspp * n_spp # number of treeid
n_meas <- 12 # repeated measurements per id
N <- n_treeid * n_meas # total number of measurements

a <- 5
sigma_y <- 1
sigma_atreeid <- 0.5
sigma_aspp <- 2
error <- rnorm(N, 0, sigma_y)

# get replicated treeid
treeid <- rep(1:n_treeid, each = n_meas)

# replicated spp
spp <- rep(rep(1:n_spp, each = n_perspp), each = n_meas) 

# get intercept values for each species
aspp <- rnorm(n_spp, 0, sigma_aspp)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Simulate data in the nested way #####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
simnest <- data.frame(
  treeid = treeid,
  spp = spp
)

atreeid <- rnorm(N, aspp[spp], sigma_atreeid)

# Add my parameters to the df
simnest$aspp <- aspp[simnest$spp]
simnest$atreeid <- atreeid[simnest$treeid]
simnest$atreeid <- simnest$atreeid - simnest$aspp

# add the rest
simnest$a <- a
simnest$sigma_y <- sigma_y
simnest$sigma_atreeid <- sigma_atreeid
simnest$sigma_aspp <- sigma_aspp
simnest$error <- error

# adding both options of tree rings
simnest$ringwidth <- 
  simnest$a +
  simnest$atreeid + 
  simnest$aspp +
  simnest$error

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Simulate data in the additive way #####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
simadd <- data.frame(
  treeid = treeid,
  spp = spp
)

# get intercept values for each species
# aspp <- rnorm(n_spp, 0, sigma_aspp)
atreeid <- rnorm(N, 0, sigma_atreeid)

# Add my parameters to the df
simadd$atreeid <- atreeid[simnest$treeid]
simadd$aspp <- aspp[simadd$spp]

# add the rest
simadd$a <- a
simadd$sigma_y <- sigma_y
simadd$sigma_atreeid <- sigma_atreeid
simadd$sigma_aspp <- sigma_aspp
simadd$error <- error

# adding both options of tree rings
simadd$ringwidth <- 
  simadd$a +
  simadd$atreeid +
  simadd$aspp + 
  simadd$error

jpeg(
  filename = "figures/simdata/simNestedVSAdditive11plot.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)
plot(simnest$ringwidth ~ simadd$ringwidth, xlab = "rw additive", ylab = "rw nested")
# plot(simnest$aspp ~ simadd$aspp, xlab = "rw additive", ylab = "rw nested")
abline(0, 1)
dev.off()

pdf(file = "figures/simData/simNestedVSAdditive.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
for (i in unique(simnest$spp)) { # i = 3
  nest_sub <- simnest[simnest$spp == i, ]
  add_sub  <- simadd[simadd$spp == i, ]
  
  plot(0, type="n",
       xlim=c(0, 10), xlab="RW",
       ylim=c(0, 100), ylab="Counts",
       main = paste("Species", i))
  
  util$plot_line_hist(nest_sub$ringwidth, 0, 25, 0.6, col="black", 
                      add=TRUE)
  util$plot_line_hist(add_sub$ringwidth, 0, 25, 0.6, col=util$c_mid_teal, 
                      add=TRUE)
  
  # abline(v = sum(unique(nest_sub$a), unique(nest_sub$aspp)), col = "black")
  abline(v = sum(unique(add_sub$a), unique(add_sub$aspp)), col = util$c_mid_teal)
  legend("topright",
         legend = c("Nested", "Additive"),
         col    = c("black", util$c_mid_teal),
         lty    = 1,
         lwd    = 2,
         bty    = "n")
}
dev.off()


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Run model ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

y <- simnest$ringwidth
N <- nrow(simnest)
Nspp <- length(unique(simnest$spp))
species <- as.numeric(as.character(simnest$spp))
treeid <- treeid
Ntreeid <- length(unique(treeid))
table(treeid)

rstan_options(auto_write = TRUE)

fitnested <- stan("stan/twolevelhierint_nested.stan", 
                    data=c("N","y",
                           "Nspp","species",
                           "Ntreeid", "treeid"),
                    warmup = 1000, iter = 2000, chains=4)

diagnostics <- util$extract_hmc_diagnostics(fitnested) 
util$check_all_hmc_diagnostics(diagnostics)

saveRDS(fitnested, "output/stanOutput/fitSimNested")
# fit <- readRDS("output/stanOutput/fitSimNested") 

# Fit additive
y <- simadd$ringwidth
fitadd <- stan("stan/twolevelhierint_nested.stan", 
                  data=c("N","y",
                         "Nspp","species",
                         "Ntreeid", "treeid"),
                  warmup = 1000, iter = 2000, chains=4)

diagnostics <- util$extract_hmc_diagnostics(fitadd) 
util$check_all_hmc_diagnostics(diagnostics)

saveRDS(fitadd, "output/stanOutput/fitSimAdditive")
# fit <- readRDS("output/stanOutput/fitSimNested") 

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Recover and plot parameters ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Nested #####
df_fit <- as.data.frame(fitnested)

# full posterior
columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
treeid_df <- df_fit[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]

# change colnames
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)

# posterior summaries
sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
treeid_df2 <- extract_params(df_fit, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")

# Plot parameter recovery
sigma_df2$sim_sigma <- c(sigma_atreeid, sigma_y)
aspp_df2$sim_aspp <- simnest$aspp[match(aspp_df2$spp, simnest$spp)]
treeid_df2$sim_atreeid <- simnest$atreeid[match(treeid_df2$treeid, simnest$treeid)]

# Plot sigmas
sigma_simXfit_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_errorbar(aes(ymin = mean_per25, ymax = mean_per75),
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 1) +
  geom_errorbar(aes(ymin = mean_per5, ymax = mean_per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  geom_point(color = "#046C9A", size = 3) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "") +
  theme_minimal()
sigma_simXfit_plot

# Plot treeid 
atreeid_simXfit_plot <- ggplot(treeid_df2, aes(x = sim_atreeid, y = fit_atreeid)) +
  geom_errorbar(aes(ymin = fit_atreeid_per5, ymax = fit_atreeid_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "") +
  theme_minimal()
atreeid_simXfit_plot

# Plot aspp 
aspp_simXfit_plot <- ggplot(aspp_df2, aes(x = sim_aspp, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per5, ymax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim aspp", y = "fit aspp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
aspp_simXfit_plot


# Combine plots 
combined_plot <- (atreeid_simXfit_plot + sigma_simXfit_plot + aspp_simXfit_plot )
combined_plot
ggsave("figures/simdata/combinedPrmRecoveryNested.jpeg", combined_plot, width = 12, height = 4, units = "in", dpi = 300)

##### Additive #####
df_fit <- as.data.frame(fitadd)

# full posterior
columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
treeid_df <- df_fit[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]

# change colnames
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)

# posterior summaries
sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
treeid_df2 <- extract_params(df_fit, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")

# Plot parameter recovery
sigma_df2$sim_sigma <- c(sigma_atreeid, sigma_y)
aspp_df2$sim_aspp <- simadd$aspp[match(aspp_df2$spp, simadd$spp)]
treeid_df2$sim_atreeid <- simadd$atreeid[match(treeid_df2$treeid, simadd$treeid)]

# Plot sigmas
sigma_simXfit_plot <- ggplot(sigma_df2, aes(x = sim_sigma, y = mean)) +
  geom_errorbar(aes(ymin = mean_per25, ymax = mean_per75),
                width = 0, linewidth = 1.5, color = "darkgray", alpha = 1) +
  geom_errorbar(aes(ymin = mean_per5, ymax = mean_per95),
                width = 0, linewidth = 0.5, color = "darkgray", alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "#B40F20", linewidth = 1) +
  geom_point(color = "#046C9A", size = 3) +
  ggrepel::geom_text_repel(aes(label = sigma), size = 3) +
  labs(x = "sim sigma", y = "fit sigma",
       title = "") +
  theme_minimal()
sigma_simXfit_plot

# Plot treeid 
atreeid_simXfit_plot <- ggplot(treeid_df2, aes(x = sim_atreeid, y = fit_atreeid)) +
  geom_errorbar(aes(ymin = fit_atreeid_per5, ymax = fit_atreeid_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.7) +
  geom_point(color = "#046C9A", size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim atreeid", y = "fit atreeid", title = "") +
  theme_minimal()
atreeid_simXfit_plot

# Plot aspp 
aspp_simXfit_plot <- ggplot(aspp_df2, aes(x = sim_aspp, y = fit_aspp)) +
  geom_errorbar(aes(ymin = fit_aspp_per5, ymax = fit_aspp_per95), 
                width = 0, linewidth = 0.5, color = "darkgray", alpha=0.9) +
  geom_errorbar(aes(ymin = fit_aspp_per25, ymax = fit_aspp_per75), 
                width = 0, linewidth = 1.5,  color = "darkgray", alpha=0.9) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B40F20", linewidth = 1) +
  labs(x = "sim aspp", y = "fit aspp", title = "") +
  geom_point(color = "#046C9A", size = 2) +
  theme_minimal()
aspp_simXfit_plot


# Combine plots 
combined_plot <- (atreeid_simXfit_plot + sigma_simXfit_plot + aspp_simXfit_plot )
combined_plot
ggsave("figures/simdata/combinedPrmRecoveryAdditive.jpeg", combined_plot, width = 12, height = 4, units = "in", dpi = 300)

#  === === === === === === === === === === === === === === === === 

# === === === === === === === === === === === === === === === === 
###### Priors VS Posterior ######
# === === === === === === === === === === === === === === === === 
# prior predictive checks. Simulating prior values from the values set in the model block
hyperparameter_draws <- 8000
parameter_draws <- 1000

###### Priors sigma_bsp ######
sigma_bsp_draw <- abs(rnorm(hyperparameter_draws, 0, 0.2))   
ggplot() +
  geom_density(data = data.frame(sigma_bsp_draw = sigma_bsp_draw),
               aes(x = sigma_bsp_draw, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = sigma_df,
               aes(x = sigma_bsp, colour = "Posterior"),
               linewidth = 0.8) +
  labs(title = "priorVSposterior_bsp",
       x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
ggsave("figures/priorsPredictiveChecks/priorVSposterior_sigma_bsp.jpeg", width = 10, height = 8, units = "in", dpi = 300)

sigma_aspp_draw <- abs(rnorm(draws, 0, 0.5))
sigma_asite_draw <- abs(rnorm(draws, 0, 0.5))
sigma_atree_draw <- abs(rnorm(draws, 0, 0.05))
sigma_y_draw <- abs(rnorm(draws, 0, 5))

######Priors bsp ######
n_sigma_bsp <- 200

# set to prior values
sigma_bsp_vec <- abs(rnorm(n_sigma_bsp, 0, 0.2))

prior_bsp <- rep(NA, parameter_draws*length(sigma_bsp_vec))

for (i in 1: length(sigma_bsp_vec)) {
  prior_bsp[((i - 1)*parameter_draws + 1):(i*parameter_draws)] <- rnorm(parameter_draws, 0, sigma_bsp_vec[i])
}
prior_bsp

# Get the posterior distribution
bspp_df3 <- bspp_df

bspp_df3$draw <- rownames(bspp_df3)

colnames(bspp_df3) <- c(paste0("spp", 1:10), "draw") 
bspp_df3

long_post_bspp <- reshape(
  bspp_df3,
  direction = "long",
  varying = paste0("spp", 1:10),
  v.names = "post_bsp",
  idvar = "draw",
  timevar = "spp"
)
long_post_bspp

ggplot() +
  geom_density(data = data.frame(prior_bsp = prior_bsp),
               aes(x = prior_bsp, colour = "Prior"),
               linewidth = 0.8) +
  geom_density(data = long_post_bspp,
               aes(x = post_bsp, colour = "Posterior", group = spp),
               linewidth = 0.3) +
  labs(title = "priorVSposterior_bsp",
       x = "BSP", y = "Density", color = "Curve") +
  scale_color_manual(values = wes_palette("AsteroidCity1")[3:4]) +
  theme_minimal()
ggsave("figures/priorsPredictiveChecks/priorVSposterior_bsp.jpeg", width = 8, height = 6, units = "in", dpi = 300)

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

# Sigma aspp
sigma_long_aspp <- data.frame(
  value  = c(sigma_df$post_sigma_aspp, sigma_df$prior_sigma_aspp),
  source = rep(c("post_sigma_aspp", "prior_sigma_aspp"), 
               each = nrow(sigma_df))
)

ggplot(sigma_long_aspp, aes(x = value, color = source, fill = source)) +
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
