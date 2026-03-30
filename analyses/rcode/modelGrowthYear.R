# Wildchrokie model
# CRD 23 April 2025
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
source('rcode/utilExtractParam.R')

emp <- read.csv("output/empiricalDataMAIN.csv")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit model GDD
# emp <- emp[!is.na(emp$year),]

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

emp$lengthMM <- emp$lengthCM*10

# transform data in vectors for GDD
data <- list( 
  y = emp$lengthMM,
  N = nrow(emp),
  Nspp = length(unique(emp$spp_num)),
  Nsite = length(unique(emp$site_num)),
  site = as.numeric(as.character(emp$site_num)),
  species = as.numeric(as.character(emp$spp_num)),
  treeid = as.numeric(emp$treeid_num),
  Ntreeid = length(unique(treeid)),
  Nyear = length(unique(emp$year)),
  year = as.integer(as.factor(emp$year)))

# Fit model GDD
rstan_options(auto_write = TRUE)

yearmodel <- stan_model("stan/modelGrowthYear.stan")
fityear <- sampling(yearmodel, data = data,
                warmup = 1000, iter=2000, chains=4)
saveRDS(fityear, "output/stanOutput/fitGrowthYear")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fityear) 
util$check_all_hmc_diagnostics(diagnostics)

# loo
extract_log_lik(yearmodel, merge_chains = FALSE)
loo(yearmodel)
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


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FULL DATA ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
emp <- read.csv("output/empiricalDataMAIN.csv")
# Fit model SOS --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
empsos <- emp[!is.na(emp$leafout),]

# transform my groups to numeric values
empsos$site_num <- match(empsos$site, unique(empsos$site))
empsos$spp_num <- match(empsos$spp, unique(empsos$spp))
empsos$treeid_num <- match(empsos$treeid, unique(empsos$treeid))

# transform data in vectors for gsl
y <- empsos$lengthCM*10 # ring width in mm
N <- nrow(empsos)
Nspp <- length(unique(empsos$spp_num))
Nsite <- length(unique(empsos$site_num))
site <- as.numeric(as.character(empsos$site_num))
species <- as.numeric(as.character(empsos$spp_num))
treeid <- as.numeric(empsos$treeid_num)
Ntreeid <- length(unique(treeid))
sos <- empsos$leafout / 5

rstan_options(auto_write = TRUE)
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsosfull <- sampling(sosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "sos"),
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOSFull")

# Fit model EOS
empeos <- emp[!is.na(emp$budset),]

# transform my groups to numeric values
empeos$site_num <- match(empeos$site, unique(empeos$site))
empeos$spp_num <- match(empeos$spp, unique(empeos$spp))
empeos$treeid_num <- match(empeos$treeid, unique(empeos$treeid))

# transform data in vectors for gsl
y <- empeos$lengthCM*10 # ring width in mm
N <- nrow(empeos)
Nspp <- length(unique(empeos$spp_num))
Nsite <- length(unique(empeos$site_num))
site <- as.numeric(as.character(empeos$site_num))
species <- as.numeric(as.character(empeos$spp_num))
treeid <- as.numeric(empeos$treeid_num)
Ntreeid <- length(unique(treeid))
eos <- empeos$budset/10

rstan_options(auto_write = TRUE)
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteosfull <- sampling(eosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "eos"),
                   warmup = 1000, iter = 2000,
                   chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOSFull")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Recover and plot parameters SOS restricted vs full #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
df_fitsos <- as.data.frame(fitsos)
# posterior summaries
sigma_df2_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("z|sigma", treeid))
aspp_df2_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("prior", treeid))
site_df2_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# SOS full
df_fitsos <- as.data.frame(fitsosfull)
# posterior summaries
sigma_df2_full_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_full_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_full_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_full_sos <- subset(treeid_df2_full_sos, !grepl("z|sigma", treeid))
aspp_df2_full_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_full_sos <- subset(treeid_df2_full_sos, !grepl("prior", treeid))
site_df2_full_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# Recover and plot parameters EOS restricted vs full 
df_fiteos <- as.data.frame(fiteos)
# posterior summaries
sigma_df2_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("z|sigma", treeid))
aspp_df2_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("prior", treeid))
site_df2_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# EOS full
df_fiteos <- as.data.frame(fiteosfull)
# posterior summaries
sigma_df2_full_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_full_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_full_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_full_eos <- subset(treeid_df2_full_eos, !grepl("z|sigma", treeid))
aspp_df2_full_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_full_eos <- subset(treeid_df2_full_eos, !grepl("prior", treeid))
site_df2_full_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# Open device
jpeg("figures/empiricalData/FullVSRestricted.jpeg", width = 9, height = 6, units = "in", res = 300)
par(mfrow = c(2,3), oma = c(0, 2, 0, 0))

plot(sigma_df2_sos$mean, sigma_df2_full_sos$mean,
     xlab = "restricted", ylab = "full", main = "", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2_full_sos$mean_per25, sigma_df2_full_sos$mean_per75)),
     xlim = range(c(sigma_df2_sos$mean_per25, sigma_df2_sos$mean_per75)))
arrows(x0 = sigma_df2_sos$mean, y0 = sigma_df2_full_sos$mean_per25,
       x1 = sigma_df2_sos$mean, y1 = sigma_df2_full_sos$mean_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_sos$mean_per25, y0 = sigma_df2_full_sos$mean,
       x1 = sigma_df2_sos$mean_per75, y1 = sigma_df2_full_sos$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_sos$mean, sigma_df2_full_sos$mean,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)
points(sigma_df2_sos$mean, sigma_df2_full_sos$mean, pch = 16, col = "#046C9A", cex = 1.5)
text(sigma_df2_sos$mean_per75, sigma_df2_full_sos$mean_per25, labels = sigma_df2_sos$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2_sos$fit_bspp, bspp_df2_full_sos$fit_bspp,
     xlab = "restricted", ylab = "full", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2_full_sos$fit_bspp_per25, bspp_df2_full_sos$fit_bspp_per75)),
     xlim = range(c(bspp_df2_sos$fit_bspp_per25, bspp_df2_sos$fit_bspp_per75)))
arrows(x0 = bspp_df2_sos$fit_bspp, y0 = bspp_df2_full_sos$fit_bspp_per25,
       x1 = bspp_df2_sos$fit_bspp, y1 = bspp_df2_full_sos$fit_bspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2_sos$fit_bspp_per25, y0 = bspp_df2_full_sos$fit_bspp,
       x1 = bspp_df2_sos$fit_bspp_per75, y1 = bspp_df2_full_sos$fit_bspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2_sos$fit_bspp, bspp_df2_full_sos$fit_bspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# aspp
plot(aspp_df2_sos$fit_aspp, aspp_df2_full_sos$fit_aspp,
     xlab = "restricted", ylab = "full", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2_full_sos$fit_aspp_per25, aspp_df2_full_sos$fit_aspp_per75)),
     xlim = range(c(aspp_df2_sos$fit_aspp_per25, aspp_df2_sos$fit_aspp_per75)))
arrows(x0 = aspp_df2_sos$fit_aspp, y0 = aspp_df2_full_sos$fit_aspp_per25,
       x1 = aspp_df2_sos$fit_aspp, y1 = aspp_df2_full_sos$fit_aspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2_sos$fit_aspp_per25, y0 = aspp_df2_full_sos$fit_aspp,
       x1 = aspp_df2_sos$fit_aspp_per75, y1 = aspp_df2_full_sos$fit_aspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2_sos$fit_aspp, aspp_df2_full_sos$fit_aspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)


# add label
mtext("a)", side = 2, outer = TRUE, at = 0.95, font = 2, las = 1, line = 0.5)

# EOS
plot(sigma_df2_eos$mean, sigma_df2_full_eos$mean,
     xlab = "restricted", ylab = "full", main = "", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2_full_eos$mean_per25, sigma_df2_full_eos$mean_per75)),
     xlim = range(c(sigma_df2_eos$mean_per25, sigma_df2_eos$mean_per75)))
arrows(x0 = sigma_df2_eos$mean, y0 = sigma_df2_full_eos$mean_per25,
       x1 = sigma_df2_eos$mean, y1 = sigma_df2_full_eos$mean_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_eos$mean_per25, y0 = sigma_df2_full_eos$mean,
       x1 = sigma_df2_eos$mean_per75, y1 = sigma_df2_full_eos$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_eos$mean, sigma_df2_full_eos$mean,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)
points(sigma_df2_eos$mean, sigma_df2_full_eos$mean, pch = 16, col = "#046C9A", cex = 1.5)
text(sigma_df2_eos$mean_per75, sigma_df2_full_eos$mean_per25, labels = sigma_df2_eos$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2_eos$fit_bspp, bspp_df2_full_eos$fit_bspp,
     xlab = "restricted", ylab = "full", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2_full_eos$fit_bspp_per25, bspp_df2_full_eos$fit_bspp_per75)),
     xlim = range(c(bspp_df2_eos$fit_bspp_per25, bspp_df2_eos$fit_bspp_per75)))
arrows(x0 = bspp_df2_eos$fit_bspp, y0 = bspp_df2_full_eos$fit_bspp_per25,
       x1 = bspp_df2_eos$fit_bspp, y1 = bspp_df2_full_eos$fit_bspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2_eos$fit_bspp_per25, y0 = bspp_df2_full_eos$fit_bspp,
       x1 = bspp_df2_eos$fit_bspp_per75, y1 = bspp_df2_full_eos$fit_bspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2_eos$fit_bspp, bspp_df2_full_eos$fit_bspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# aspp
plot(aspp_df2_eos$fit_aspp, aspp_df2_full_eos$fit_aspp,
     xlab = "restricted", ylab = "full", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2_full_eos$fit_aspp_per25, aspp_df2_full_eos$fit_aspp_per75)),
     xlim = range(c(aspp_df2_eos$fit_aspp_per25, aspp_df2_eos$fit_aspp_per75)))
arrows(x0 = aspp_df2_eos$fit_aspp, y0 = aspp_df2_full_eos$fit_aspp_per25,
       x1 = aspp_df2_eos$fit_aspp, y1 = aspp_df2_full_eos$fit_aspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2_eos$fit_aspp_per25, y0 = aspp_df2_full_eos$fit_aspp,
       x1 = aspp_df2_eos$fit_aspp_per75, y1 = aspp_df2_full_eos$fit_aspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2_eos$fit_aspp, aspp_df2_full_eos$fit_aspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# add label
mtext("b)", side = 2, outer = TRUE, at = 0.42, font = 2, las = 1, line = 0.5)

dev.off()

}