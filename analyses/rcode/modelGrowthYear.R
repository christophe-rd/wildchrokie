# Wildchrokie model
# CRD 30 March 2026
# compare a year model as the predictor with the gdd model and cross-validate to figure out which one is the most predictive.

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)

# Load library 
library(loo)
library(rstan)

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
# my function to extract parameters
source('rcode/tools.R')

emp <- read.csv("output/empiricalDataMAIN.csv")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit model GDD
emp <- emp[!is.na(emp$pgsGDD5),]

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))
emp$year_num <- as.integer(as.factor(emp$year))
emp$lengthMM <- emp$lengthCM*10

# transform data in vectors for GDD
data <- list( 
  y = log(emp$lengthMM),
  N = nrow(emp),
  Nspp = length(unique(emp$spp_num)),
  Nsite = length(unique(emp$site_num)),
  site = as.numeric(as.character(emp$site_num)),
  species = as.numeric(as.character(emp$spp_num)),
  treeid = as.numeric(emp$treeid_num),
  Ntreeid = length(unique(emp$treeid)),
  Nyear = length(unique(emp$year)),
  year = emp$year_num)

# Fit model year with other parameters
yearmodel <- stan_model("stan/modelGrowthYear.stan")
fityear <- sampling(yearmodel, data = data,
                warmup = 1000, iter=2000, chains=4)
saveRDS(fityear, "output/stanOutput/fitGrowthYear")

# Fit model year without other parameters
yearmodelonly <- stan_model("stan/modelGrowthOnlyYear.stan")
fityearonly <- sampling(yearmodelonly, data = data,
                    warmup = 1000, iter=2000, chains=4)
saveRDS(fityearonly, "output/stanOutput/fitGrowthOnlyYear")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fityearonly) 
util$check_all_hmc_diagnostics(diagnostics)

# loo
log_lik_1 <- extract_log_lik(fityear, merge_chains = FALSE)

# relative effective sample sizes
r_eff <- relative_eff(exp(log_lik_1)) 

loo_1 <- loo(log_lik_1, r_eff = r_eff)

plot(loo_1)
# which observations are problematic
pareto_k <- loo_1$diagnostics$pareto_k
bad_obs <- which(pareto_k > 0.7)
bad_obs

emp[bad_obs, ]

# comparison 
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD")
log_lik_gdd <- extract_log_lik(fitgdd, merge_chains = FALSE)
r_eff_gdd <- relative_eff(exp(log_lik_gdd)) 
loo_gdd <- loo(log_lik_gdd, r_eff = r_eff_gdd)

comp <- loo_compare(loo_1, loo_gdd)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot Year fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fityear <- as.data.frame(fityear)

# full posterior
columns <- colnames(df_fityear)[!grepl("prior", colnames(df_fityear))]
sigma_df <- df_fityear[, columns[grepl("sigma", columns)]]
treeid_df <- df_fityear[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fityear[, columns[grepl("aspp", columns)]]
site_df <- df_fityear[, columns[grepl("asite", columns)]]
ayear_df <- df_fityear[, columns[grepl("ayear", columns)]]

# change colnames
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
sigma_df2  <- extract_params(df_fityear, "sigma", "mean", "sigma")
treeid_df2 <- extract_params(df_fityear, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fityear, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fityear, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")
year_df2   <- extract_params(df_fityear, "ayear", "fit_ayear", 
                             "year", "ayear\\[(\\d+)\\]")

##### Plot posterior vs priors for year fit #####
pdf(file = "figures/growthYearModel/yearModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fityear[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fityear[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fityear[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fityear[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fityear[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2.5))
lines(density(df_fityear[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fityear[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.2))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fityear[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fityear[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", ylim = c(0, 0.6))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot Year fit ONLY ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fityearonly <- as.data.frame(fityearonly)

# full posterior
columns <- colnames(df_fityearonly)[!grepl("prior", colnames(df_fityearonly))]
sigma_df <- df_fityearonly[, columns[grepl("sigma", columns)]]
ayear_df <- df_fityearonly[, columns[grepl("ayear", columns)]]

# change colnames
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
year_df2_only <- extract_params(df_fityearonly, "ayear", "fit_ayear", 
                             "year", "ayear\\[(\\d+)\\]")
year_df2_only$year_name <- emp$year[match(year_df2_only$year, emp$year_num)]

##### Plot posterior vs priors for year fit #####
pdf(file = "figures/growthYearModel/yearOnlyModelPriorVSPosterior.pdf", width = 6, height = 9)

par(mfrow = c(2, 1))
pal <- wes_palette("AsteroidCity1")[3:4]

# a
plot(density(df_fityearonly[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fityearonly[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fityearonly[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", ylim = c(0, 0.6))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

colsyear <- c("#dd5129", "#0f7ba2", "#43b284")
# plots years
jpeg("figures/growthYearModel/muYears.jpeg", width = 6, height = 6, 
     units = "in", res = 300)
par(mfrow = c(1,1))
n_year <- length(unique(year_df2_only$year))
y_pos <- 1:n_year

# Current year
plot(year_df2_only$fit_ayear, y_pos,
     xlim = c(-5, 5), ylim = c(0.5, n_year + 0.5), 
     xlab = "year intercept", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = colsyear, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(year_df2_only$fit_ayear_per5, y_pos, year_df2_only$fit_ayear_per95, y_pos,
         col = colsyear, lwd = 1.5)
segments(year_df2_only$fit_ayear_per25, y_pos, year_df2_only$fit_ayear_per75, y_pos,
         col = colsyear, lwd = 3)
axis(2, at = y_pos, labels = year_df2_only$year_name, las = 1)
dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (FALSE){
  
samples <- util$extract_expectand_vals(fityear)
jpeg(
  filename = "figures/growthYearModel/retrodictiveCheckHist.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
util$plot_hist_quantiles(samples, "y_rep", 
                         -5, # lower x axis limit
                         15, # upper x axis limit
                         0.5, # binning
                         baseline_values = data$y,
                         xlab = "Ring width (mm)")
dev.off()
}
