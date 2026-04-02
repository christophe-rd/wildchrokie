# Wildchrokie model
# CRD 23 April 2025
# Started in Boston shortly after completing field work for the tree spotters

# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)

# Load library 
library(rstan)
library(wesanderson)

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

# flags
fitmodels <- FALSE
fitmodelsZscored <- FALSE

emp <- read.csv("output/empiricalDataMAIN.csv")

# change ring width
emp$lengthMM <- emp$lengthCM*10
emp$loglength <- log(emp$lengthMM)

empfullsos <- emp[!is.na(emp$leafout),]
empfulleos <- emp[!is.na(emp$budset),]

gddyr <- read.csv("output/gddByYear.csv")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit model GDD
emp <- emp[!is.na(emp$pgsGDD5),]

# scale gdd to how many gdd are in 10 average spring days
temp<- subset(gddyr, doy <151 & doy > 120)
temp$mingddperiod <- ave(temp$GDD_5, temp$year, FUN = min)
temp$gdddiff <- temp$GDD_5 - temp$mingddperiod

temp <- temp[order(temp$year, temp$doy), ]

temp$bin10 <- ave(temp$doy, temp$year, FUN = function(x) ceiling((x - min(x) + 1) / 10))
gdd_10day <- aggregate(gdddiff ~ year + bin10, data = temp, max)
gddscale <- mean(gdd_10day$gdddiff)

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors for GDD
y <- emp$loglength # ring width in mm
N <- nrow(emp)
Nspp <- length(unique(emp$spp_num))
Nsite <- length(unique(emp$site_num))
site <- as.numeric(as.character(emp$site_num))
species <- as.numeric(as.character(emp$spp_num))
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))

gdd <- emp$pgsGDD5 / gddscale
gsl <- emp$pgsGSL / 10
sos <- emp$leafout / 5
eos <- emp$budset / 10

if (fitmodels){
# Fit model GDD
gddmodel <- stan_model("stan/modelGrowthGDD.stan")
fitgdd <- sampling(gddmodel, data = c("N","y",
                                "Nspp","species",
                                "Nsite", "site", 
                                "Ntreeid", "treeid", 
                                "gdd"),
                warmup = 1000, iter=2000, chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDD")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fitgdd) 
util$check_all_hmc_diagnostics(diagnostics)

# Fit model GSL
gslmodel <- stan_model("stan/modelGrowthGSL.stan")
fitgsl <- sampling(gslmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site", 
                                  "Ntreeid", "treeid", 
                                  "gsl"),
                warmup = 1000, iter = 2000, chains = 4)
saveRDS(fitgsl, "output/stanOutput/fitGrowthGSL")

# Fit model SOS
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsos <- sampling(sosmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site",
                                  "Ntreeid", "treeid",
                                  "sos"),
                warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOS")

# Fit model EOS
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteos <- sampling(eosmodel, data = c("N","y",
                                  "Nspp","species",
                                  "Nsite", "site",
                                  "Ntreeid", "treeid",
                                  "eos"),
                warmup = 1000, iter = 2000, chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOS")


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GDD fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gddModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitgdd[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitgdd[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgdd[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitgdd[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgdd[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitgdd[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgdd[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgdd[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgdd[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GSL fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitgsl <- as.data.frame(fitgsl)

# full posterior
columns <- colnames(df_fitgsl)[!grepl("prior", colnames(df_fitgsl))]
sigma_df <- df_fitgsl[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgsl[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgsl[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fitgsl[, columns[grepl("aspp", columns)]]
site_df <- df_fitgsl[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgsl, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgsl, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgsl, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgsl, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fitgsl, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for GSL fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gslModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitgsl[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitgsl[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgsl[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitgsl[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgsl[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitgsl[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgsl[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgsl[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgsl[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot SOS fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitsos <- as.data.frame(fitsos)

# full posterior
columns <- colnames(df_fitsos)[!grepl("prior", colnames(df_fitsos))]
sigma_df <- df_fitsos[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitsos[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitsos[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fitsos[, columns[grepl("aspp", columns)]]
site_df <- df_fitsos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("z|sigma", treeid))
aspp_df2_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("prior", treeid))
site_df2_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for sos fit #####
pdf(file = "figures/growthModelsMain/diagnostics/sosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitsos[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitsos[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitsos[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitsos[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitsos[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitsos[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitsos[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitsos[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitsos[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot EOS fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fiteos <- as.data.frame(fiteos)

# full posterior
columns <- colnames(df_fiteos)[!grepl("prior", colnames(df_fiteos))]
sigma_df <- df_fiteos[, columns[grepl("sigma", columns)]]
bspp_df <- df_fiteos[, columns[grepl("bsp", columns)]]
treeid_df <- df_fiteos[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fiteos[, columns[grepl("aspp", columns)]]
site_df <- df_fiteos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("z|sigma", treeid))
aspp_df2_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("prior", treeid))
site_df2_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for eos fit #####
pdf(file = "figures/growthModelsMain/diagnostics/eosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fiteos[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fiteos[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fiteos[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fiteos[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fiteos[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fiteos[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fiteos[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fiteos[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fiteos[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (FALSE){
  
samples <- util$extract_expectand_vals(fitgsl)
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Diagnostics ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Parameterization for treeid
if (FALSE) { 
  samples <- util$extract_expectand_vals(fit)
  
  # atreeid
  atreeid <- names(samples)[grepl("zatreeid", names(samples))]
  atreeid <- atreeid[!grepl("sigma", atreeid)]
  atreeid <- atreeid[sample(length(unique(atreeid)), 9)]
  # pdf("figures/troubleShootingGrowthModel/atreeidParameterization_only_atreeid_no_b.pdf", 
  #     width = 6, height = 18)
  jpeg("figures/atreeidParameterization.jpeg", 
       width = 2000, height = 3000,
       units = "px", res = 300)
  util$plot_div_pairs(atreeid, "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
  dev.off()
}

# other diagnostics?
if (FALSE) {
  
  diagnostics <- util$extract_hmc_diagnostics(fit_noncentered) 
  util$check_all_hmc_diagnostics(diagnostics)
  
  samples <- util$extract_expectand_vals(fit_noncentered)
  
  util$plot_div_pairs("zbsp[1]", "sigma_bsp", samples, diagnostics, transforms = list("sigma_bsp" = 1))
  
  util$plot_div_pairs("zaspp[1]", "sigma_aspp", samples, diagnostics, transforms = list("sigma_aspp" = 1))
  
  util$plot_div_pairs("zasite[1]", "sigma_asite", samples, diagnostics, transforms = list("sigma_asite" = 1))
  
  util$plot_div_pairs("atreeid[1]", "sigma_atreeid", samples, diagnostics, transforms = list("sigma_atreeid" = 1))
  
}

}


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FULL DATA ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit model SOS
# transform my groups to numeric values
empfullsos$site_num <- match(empfullsos$site, unique(empfullsos$site))
empfullsos$spp_num <- match(empfullsos$spp, unique(empfullsos$spp))
empfullsos$treeid_num <- match(empfullsos$treeid, unique(empfullsos$treeid))

# transform data in vectors for gsl
y <- empfullsos$loglength # ring width in mm
N <- nrow(empfullsos)
Nspp <- length(unique(empfullsos$spp_num))
Nsite <- length(unique(empfullsos$site_num))
site <- as.numeric(as.character(empfullsos$site_num))
species <- as.numeric(as.character(empfullsos$spp_num))
treeid <- as.numeric(empfullsos$treeid_num)
Ntreeid <- length(unique(treeid))
sos <- empfullsos$leafout / 5


sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsosfull <- sampling(sosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "sos"),
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOSFull")

# Fit model EOS
# transform my groups to numeric values
empfulleos$site_num <- match(empfulleos$site, unique(empfulleos$site))
empfulleos$spp_num <- match(empfulleos$spp, unique(empfulleos$spp))
empfulleos$treeid_num <- match(empfulleos$treeid, unique(empfulleos$treeid))

# transform data in vectors for gsl
y <- empfulleos$loglength # ring width in mm
N <- nrow(empfulleos)
Nspp <- length(unique(empfulleos$spp_num))
Nsite <- length(unique(empfulleos$site_num))
site <- as.numeric(as.character(empfulleos$site_num))
species <- as.numeric(as.character(empfulleos$spp_num))
treeid <- as.numeric(empfulleos$treeid_num)
Ntreeid <- length(unique(treeid))
eos <- empfulleos$budset/10

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
jpeg("figures/growthModelsMain/FullVSRestricted.jpeg", width = 9, height = 6, units = "in", res = 300)
par(mfrow = c(2,3), oma = c(0, 2, 0, 0))

plot(sigma_df2_sos$mean, sigma_df2_full_sos$mean,
     xlab = "restricted", ylab = "full", main = "sigmas", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2_full_sos$mean_per25, sigma_df2_full_sos$mean_per75)),
     xlim = range(c(sigma_df2_sos$mean_per25, sigma_df2_sos$mean_per75)))
arrows(x0 = sigma_df2_sos$mean, y0 = sigma_df2_full_sos$mean_per25,
       x1 = sigma_df2_sos$mean, y1 = sigma_df2_full_sos$mean_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_sos$mean_per25, y0 = sigma_df2_full_sos$mean,
       x1 = sigma_df2_sos$mean_per75, y1 = sigma_df2_full_sos$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_sos$mean, sigma_df2_full_sos$mean,
       pch = 16, col = "#79ad41", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
points(sigma_df2_sos$mean, sigma_df2_full_sos$mean, pch = 16, col = "#79ad41", cex = 1.5)
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
       pch = 16, col = "#79ad41", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

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
       pch = 16, col = "#79ad41", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)


# add label
mtext("a)", side = 2, outer = TRUE, at = 0.95, font = 2, las = 1, line = 0.5)

# EOS
plot(sigma_df2_eos$mean, sigma_df2_full_eos$mean,
     xlab = "restricted", ylab = "full", main = "sigmas", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2_full_eos$mean_per25, sigma_df2_full_eos$mean_per75)),
     xlim = range(c(sigma_df2_eos$mean_per25, sigma_df2_eos$mean_per75)))
arrows(x0 = sigma_df2_eos$mean, y0 = sigma_df2_full_eos$mean_per25,
       x1 = sigma_df2_eos$mean, y1 = sigma_df2_full_eos$mean_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_eos$mean_per25, y0 = sigma_df2_full_eos$mean,
       x1 = sigma_df2_eos$mean_per75, y1 = sigma_df2_full_eos$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_eos$mean, sigma_df2_full_eos$mean,
       pch = 16, col = "#e67424", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
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
       pch = 16, col = "#e67424", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

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
       pch = 16, col = "#e67424", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# add label
mtext("b)", side = 2, outer = TRUE, at = 0.42, font = 2, las = 1, line = 0.5)

dev.off()
}

# === === === === === === === === === === === === === === === === === === === #
# === === === === === === === === === === === === === === === === === === === #

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Z-SCORED ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (fitmodelsZscored) {
# different response variables Z scored
gdd <- (emp$pgsGDD5 - mean(emp$pgsGDD5)) / sd(emp$pgsGDD5)
gsl <- (emp$pgsGSL - mean(emp$pgsGSL)) / sd(emp$pgsGSL)
sos <- (emp$leafout - mean(emp$leafout)) / sd(emp$leafout)
eos <- (emp$budset - mean(emp$budset)) / sd(emp$budset)

gddmodel <- stan_model("stan/modelGrowthGDD.stan")
fitgdd <- sampling(gddmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gdd"),
                   warmup = 1000, iter=2000, chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDDZscored")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fitgdd) 
util$check_all_hmc_diagnostics(diagnostics)

# Fit model GSL
gslmodel <- stan_model("stan/modelGrowthGSL.stan")
fitgsl <- sampling(gslmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site", 
                                      "Ntreeid", "treeid", 
                                      "gsl"),
                   warmup = 1000, iter = 2000, chains = 4)
saveRDS(fitgsl, "output/stanOutput/fitGrowthGSLZscored")

# Fit model SOS
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsos <- sampling(sosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "sos"),
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOSZscored")

# Fit model EOS
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteos <- sampling(eosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "eos"),
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOSZscored")


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GDD fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/growthModelsMain/diagnostics/ZgddModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitgdd[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitgdd[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgdd[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitgdd[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgdd[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitgdd[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgdd[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgdd[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgdd[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GSL fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitgsl <- as.data.frame(fitgsl)

# full posterior
columns <- colnames(df_fitgsl)[!grepl("prior", colnames(df_fitgsl))]
sigma_df <- df_fitgsl[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgsl[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgsl[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fitgsl[, columns[grepl("aspp", columns)]]
site_df <- df_fitgsl[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgsl, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgsl, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgsl, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgsl, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fitgsl, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for GSL fit #####
pdf(file = "figures/growthModelsMain/diagnostics/ZgslModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitgsl[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitgsl[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgsl[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitgsl[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgsl[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitgsl[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgsl[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgsl[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgsl[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot SOS fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fitsos <- as.data.frame(fitsos)

# full posterior
columns <- colnames(df_fitsos)[!grepl("prior", colnames(df_fitsos))]
sigma_df <- df_fitsos[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitsos[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitsos[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fitsos[, columns[grepl("aspp", columns)]]
site_df <- df_fitsos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("z|sigma", treeid))
aspp_df2_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("prior", treeid))
site_df2_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for sos fit #####
pdf(file = "figures/growthModelsMain/diagnostics/ZsosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fitsos[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitsos[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitsos[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitsos[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitsos[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitsos[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitsos[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitsos[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitsos[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot EOS fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fiteos <- as.data.frame(fiteos)

# full posterior
columns <- colnames(df_fiteos)[!grepl("prior", colnames(df_fiteos))]
sigma_df <- df_fiteos[, columns[grepl("sigma", columns)]]
bspp_df <- df_fiteos[, columns[grepl("bsp", columns)]]
treeid_df <- df_fiteos[, grepl("treeid", columns) & 
                         !grepl("z|sigma", columns)]
aspp_df <- df_fiteos[, columns[grepl("aspp", columns)]]
site_df <- df_fiteos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("z|sigma", treeid))
aspp_df2_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("prior", treeid))
site_df2_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

##### Plot posterior vs priors for eos fit #####
pdf(file = "figures/growthModelsMain/diagnostics/ZeosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fiteos[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fiteos[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fiteos[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fiteos[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fiteos[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fiteos[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fiteos[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fiteos[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fiteos[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

}