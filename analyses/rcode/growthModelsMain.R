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

# flags
fitmodels <- F
fitmodelfull <- F
fitmodelsZscored <- F

emp <- read.csv("output/empiricalDataMAIN.csv")

# change ring width
emp$lengthMM <- emp$lengthCM*10
emp$loglength <- log(emp$lengthMM)

empfullsos <- emp[!is.na(emp$leafout),]
empfulleos <- emp[!is.na(emp$budset),]

gddyr <- read.csv("output/gddByYear.csv")
nrow(empfullsos)
nrow(empfulleos)

lineplotseqlength <- 10
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Most restricted amount of data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
emp <- emp[!is.na(emp$pgsGDD5),]

# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))
emp$year_num <- match(emp$year, unique(emp$year))

# order by tree id
treeid_spp_site <- unique(emp[, c("treeid_num", "spp_num", "site_num",
                                  "treeid", "spp", "site", 
                                  "latbi")])

treeid_spp_site_ordered <- treeid_spp_site[order(treeid_spp_site$treeid_num), ]

# scale gdd to how many gdd are in 10 average spring days
temp<- subset(gddyr, doy <151 & doy > 120)
temp$mingddperiod <- ave(temp$GDD_5, temp$year, FUN = min)
temp$gdddiff <- temp$GDD_5 - temp$mingddperiod

temp <- temp[order(temp$year, temp$doy), ]

temp$bin10 <- ave(temp$doy, temp$year, FUN = function(x) ceiling((x - min(x) + 1) / 7))
gdd_7day <- aggregate(gdddiff ~ year + bin10, data = temp, max)
wcgddscale <- mean(gdd_7day$gdddiff)

gddseq <-  seq(min(emp$pgsGDD5), max(emp$pgsGDD5), length.out = lineplotseqlength)

# data list for gdd
dgdd <- list(
  y = emp$loglength,
  N = nrow(emp),
  Nspp = length(unique(emp$spp_num)),
  Nsite = length(unique(emp$site_num)),
  site = as.numeric(as.character(emp$site_num)),
  species = as.numeric(as.character(emp$spp_num)),
  treeid = as.numeric(emp$treeid_num),
  Ntreeid = length(unique(as.numeric(emp$treeid_num))),
  year = as.numeric(emp$year_num),
  Nyear = length(unique(emp$year_num)),
  treeid_species = treeid_spp_site_ordered$spp_num,
  treeid_site = treeid_spp_site_ordered$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered$spp_num)),
  gdd = emp$pgsGDD5 / wcgddscale,
  gddseq = gddseq,
  wcgddscale = wcgddscale,
  Ngddseq = length(gddseq)
)
dgdd

# Set model GSL data
gslscale <- 7
gsl <- emp$pgsGSL / gslscale
gslseq <-  seq(min(emp$pgsGSL), max(emp$pgsGSL), length.out = lineplotseqlength)

# data list for GSL
dgsl <- list(
  y = emp$loglength,
  N = nrow(emp),
  Nspp = length(unique(emp$spp_num)),
  Nsite = length(unique(emp$site_num)),
  site = as.numeric(as.character(emp$site_num)),
  species = as.numeric(as.character(emp$spp_num)),
  treeid = as.numeric(emp$treeid_num),
  Ntreeid = length(unique(as.numeric(emp$treeid_num))),
  year = as.numeric(emp$year_num),
  Nyear = length(unique(emp$year_num)),
  treeid_species = treeid_spp_site_ordered$spp_num,
  treeid_site = treeid_spp_site_ordered$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered$spp_num)),
  gsl = emp$pgsGSL / gslscale,
  gslseq = gslseq,
  gslscale = gslscale,
  Ngslseq = length(gslseq)
)

sosscale <- 7
sos <- emp$leafout / sosscale
sosseq <-  seq(min(emp$leafout), max(emp$leafout), length.out = lineplotseqlength)

# data list for sos
dsos <- list(
  y = emp$loglength,
  N = nrow(emp),
  Nspp = length(unique(emp$spp_num)),
  Nsite = length(unique(emp$site_num)),
  site = as.numeric(as.character(emp$site_num)),
  species = as.numeric(as.character(emp$spp_num)),
  treeid = as.numeric(emp$treeid_num),
  Ntreeid = length(unique(as.numeric(emp$treeid_num))),
  year = as.numeric(emp$year_num),
  Nyear = length(unique(emp$year_num)),
  treeid_species = treeid_spp_site_ordered$spp_num,
  treeid_site = treeid_spp_site_ordered$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered$spp_num)),
  sos = emp$leafout / sosscale,
  sosseq = sosseq,
  sosscale = sosscale,
  Nsosseq = length(sosseq)
)

eosscale <- 7
eos <- emp$budset / eosscale
eosseq <-  seq(min(emp$budset), max(emp$budset), length.out = lineplotseqlength)

# data list for eos
deos <- list(
  y = emp$loglength,
  N = nrow(emp),
  Nspp = length(unique(emp$spp_num)),
  Nsite = length(unique(emp$site_num)),
  site = as.numeric(as.character(emp$site_num)),
  species = as.numeric(as.character(emp$spp_num)),
  treeid = as.numeric(emp$treeid_num),
  Ntreeid = length(unique(as.numeric(emp$treeid_num))),
  year = as.numeric(emp$year_num),
  Nyear = length(unique(emp$year_num)),
  treeid_species = treeid_spp_site_ordered$spp_num,
  treeid_site = treeid_spp_site_ordered$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered$spp_num)),
  eos = emp$budset / eosscale,
  eosseq = eosseq,
  eosscale = eosscale,
  Neosseq = length(eosseq)
)

if (fitmodels){
# Fit model GDD
gddmodel <- stan_model("stan/modelGrowthGDD.stan")
fitgdd <- sampling(gddmodel, data = dgdd,
                warmup = 1500, iter = 3000, chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDD")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fitgdd) 
util$check_all_hmc_diagnostics(diagnostics)

# fit gsl
gslmodel <- stan_model("stan/modelGrowthGSL.stan")
fitgsl <- sampling(gslmodel, data = dgsl,
                warmup = 1000, iter = 2000, chains = 4)
saveRDS(fitgsl, "output/stanOutput/fitGrowthGSL")

# Fit model SOS
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsos <- sampling(sosmodel, data = dsos,
                warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOS")

# Fit model EOS
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteos <- sampling(eosmodel, data = deos,
                warmup = 1500, iter = 3000, chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOS")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GDD fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD")
fitgsl <- readRDS("output/stanOutput/fitGrowthGSL")
fitsos <- readRDS("output/stanOutput/fitGrowthSOS")
fiteos <- readRDS("output/stanOutput/fitGrowthEOS")

##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma|slope|full", columns)]
aspp_df <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df <- df_fitgdd[, columns[grepl("asite", columns)]]
ayear_df <- df_fitgdd[, columns[grepl("ayear", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

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
site_df2 <- subset(site_df2, !grepl("z|sigma", site))
ayear_df2  <- extract_params(df_fitgdd, "ayear", "fit_ayear", 
                             "year", "ayear\\[(\\d+)\\]")
ayear_df2 <- subset(ayear_df2, !grepl("mean", year))

# save csvs
write.csv(sigma_df2,  "output/GM_GDDparam_sigma.csv",  row.names = FALSE)
write.csv(bspp_df2,   "output/GM_GDDparam_bspp.csv",   row.names = FALSE)
write.csv(treeid_df2, "output/GM_GDDparam_treeid.csv", row.names = FALSE)
write.csv(aspp_df2,   "output/GM_GDDparam_aspp.csv",   row.names = FALSE)
write.csv(site_df2,   "output/GM_GDDparam_site.csv",   row.names = FALSE)
write.csv(ayear_df2,  "output/GM_GDDparam_ayear.csv",  row.names = FALSE)

##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gddModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(4, 2))

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

# sigma_asite
plot(density(df_fitgdd[, "sigma_asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_asite", 
     xlab = "sigma_asite", ylim = c(0,4))
lines(density(df_fitgdd[, "sigma_asite"]), col = pal[2], lwd = 2)
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
     ylim = c(0, 0.3))
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

# ayear
plot(density(df_fitgdd[, "ayear_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_ayear",
     xlab = "ayear", xlim = c(-3, 3), ylim = c(0, 1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
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
ayear_df <- df_fitgsl[, columns[grepl("ayear", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgsl, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgsl, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgsl, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgsl, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fitgsl, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")
site_df2 <- subset(site_df2, !grepl("z|sigma", site))
ayear_df2  <- extract_params(df_fitgsl, "ayear", "fit_ayear", "year", "ayear\\[(\\d+)\\]")
ayear_df2 <- subset(ayear_df2, !grepl("mean", year))

# save csvs
write.csv(sigma_df2,  "output/GM_GSLparam_sigma.csv",  row.names = FALSE)
write.csv(bspp_df2,   "output/GM_GSLparam_bspp.csv",   row.names = FALSE)
write.csv(treeid_df2, "output/GM_GSLparam_treeid.csv", row.names = FALSE)
write.csv(aspp_df2,   "output/GM_GSLparam_aspp.csv",   row.names = FALSE)
write.csv(site_df2,   "output/GM_GSLparam_site.csv",   row.names = FALSE)
write.csv(ayear_df2,  "output/GM_GSLparam_ayear.csv",  row.names = FALSE)

##### Plot posterior vs priors for GSL fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gslModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(4, 2))

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

# sigma_asite
plot(density(df_fitgsl[, "sigma_asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_asite", 
     xlab = "sigma_asite", ylim = c(0,4))
lines(density(df_fitgsl[, "sigma_asite"]), col = pal[2], lwd = 2)
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
     ylim = c(0, 0.3))
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

# ayear
plot(density(df_fitgsl[, "ayear_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_ayear",
     xlab = "ayear", xlim = c(-3, 3), ylim = c(0, 1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
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
ayear_df <- df_fitsos[, columns[grepl("ayear", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
sigma_df2_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("z|sigma", treeid))
aspp_df2_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("prior", treeid))
site_df2_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")
site_df2_sos <- subset(site_df2_sos, !grepl("z|sigma", site))
ayear_df2_sos  <- extract_params(df_fitsos, "ayear", "fit_ayear", "year", "ayear\\[(\\d+)\\]")
ayear_df2_sos <- subset(ayear_df2_sos, !grepl("mean", year))

# save csvs
write.csv(sigma_df2_sos,  "output/GM_SOSparam_sigma.csv",  row.names = FALSE)
write.csv(bspp_df2_sos,   "output/GM_SOSparam_bspp.csv",   row.names = FALSE)
write.csv(treeid_df2_sos, "output/GM_SOSparam_treeid.csv", row.names = FALSE)
write.csv(aspp_df2_sos,   "output/GM_SOSparam_aspp.csv",   row.names = FALSE)
write.csv(site_df2_sos,   "output/GM_SOSparam_site.csv",   row.names = FALSE)
write.csv(ayear_df2_sos,  "output/GM_SOSparam_ayear.csv",  row.names = FALSE)

##### Plot posterior vs priors for sos fit #####
pdf(file = "figures/growthModelsMain/diagnostics/sosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(4, 2))

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

# sigma_asite
plot(density(df_fitsos[, "sigma_asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_asite", 
     xlab = "sigma_asite", ylim = c(0,4))
lines(density(df_fitsos[, "sigma_asite"]), col = pal[2], lwd = 2)
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
     ylim = c(0, 0.3))
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

# ayear
plot(density(df_fitsos[, "ayear_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_ayear",
     xlab = "ayear", xlim = c(-3, 3), ylim = c(0, 1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
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
ayear_df <- df_fiteos[, columns[grepl("ayear", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
sigma_df2_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("z|sigma", treeid))
aspp_df2_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("prior", treeid))
site_df2_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")
site_df2_eos <- subset(site_df2_eos, !grepl("z|sigma", site))
ayear_df2_eos  <- extract_params(df_fiteos, "ayear", "fit_ayear", "year", "ayear\\[(\\d+)\\]")
ayear_df2_eos <- subset(ayear_df2_eos, !grepl("mean", year))

# save csvs
write.csv(sigma_df2_eos,  "output/GM_EOSparam_sigma.csv",  row.names = FALSE)
write.csv(bspp_df2_eos,   "output/GM_EOSparam_bspp.csv",   row.names = FALSE)
write.csv(treeid_df2_eos, "output/GM_EOSparam_treeid.csv", row.names = FALSE)
write.csv(aspp_df2_eos,   "output/GM_EOSparam_aspp.csv",   row.names = FALSE)
write.csv(site_df2_eos,   "output/GM_EOSparam_site.csv",   row.names = FALSE)
write.csv(ayear_df2_eos,  "output/GM_EOSparam_ayear.csv",  row.names = FALSE)

##### Plot posterior vs priors for eos fit #####
pdf(file = "figures/growthModelsMain/diagnostics/eosModelPriorVSPosterior.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(4, 2))

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

# sigma_asite
plot(density(df_fiteos[, "sigma_asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_asite", 
     xlab = "sigma_asite", ylim = c(0,4))
lines(density(df_fiteos[, "sigma_asite"]), col = pal[2], lwd = 2)
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
     ylim = c(0, 0.3))
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

# ayear
plot(density(df_fiteos[, "ayear_prior"]),
     col = pal[1], lwd = 2,
     main = "priorVSposterior_ayear",
     xlab = "ayear", xlim = c(-3, 3), ylim = c(0, 1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Diagnostics ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Parameterization for treeid
if (FALSE) { 
  samples <- util$extract_expectand_vals(fitgdd)
  
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
  diagnostics_gdd <- util$extract_hmc_diagnostics(fitgdd) 
  diagnostics_gsl <- util$extract_hmc_diagnostics(fitgsl) 
  diagnostics_sos <- util$extract_hmc_diagnostics(fitsos) 
  diagnostics_eos <- util$extract_hmc_diagnostics(fiteos) 
  
  samples_gdd <- util$extract_expectand_vals(fitgdd)
  samples_gsl <- util$extract_expectand_vals(fitgsl)
  samples_sos <- util$extract_expectand_vals(fitsos)
  samples_eos <- util$extract_expectand_vals(fiteos)
  
  util$plot_div_pairs("zatreeid[1]", "sigma_atreeid", samples_gdd, diagnostics_gdd, transforms = list("sigma_atreeid" = 1))
  
  # check aspp
  aspp <- paste0("aspp[", 1:4, "]")
  util$plot_div_pairs(aspp, aspp, samples_gdd, diagnostics_gdd)
  util$plot_div_pairs(aspp, aspp, samples_gsl, diagnostics_gsl)
  util$plot_div_pairs(aspp, aspp, samples_sos, diagnostics_sos)
  util$plot_div_pairs(aspp, aspp, samples_eos, diagnostics_eos)
  
  # check asite vs a
  asite <- paste0("asite[", 1:4, "]")
  util$plot_div_pairs(asite, "a", samples_gdd, diagnostics_gdd)
  util$plot_div_pairs(asite, "a", samples_gsl, diagnostics_gsl)
  util$plot_div_pairs(asite, "a", samples_sos, diagnostics_sos)
  util$plot_div_pairs(asite, "a", samples_eos, diagnostics_eos)
  
  # check asite vs aspp
  par(mfrow = c(4,4))
  pdf("figures/growthModelsMain/diagnostics/pairsSiteVSaspp.pdf", 
      width = 6, height = 9)
  util$plot_div_pairs(asite, aspp, samples_gdd, diagnostics_gdd)
  dev.off()
  util$plot_div_pairs(asite, aspp, samples_gsl, diagnostics_gsl)
  util$plot_div_pairs(asite, aspp, samples_sos, diagnostics_sos)
  util$plot_div_pairs(asite, aspp, samples_eos, diagnostics_eos)
  dev.off()  
  
  # check bspp
  bspp <- paste0("bsp[", 1:4, "]")
  util$plot_div_pairs(bspp, bspp, samples_gdd, diagnostics_gdd)
  util$plot_div_pairs(bspp, bspp, samples_gsl, diagnostics_gsl)
  util$plot_div_pairs(bspp, bspp, samples_sos, diagnostics_sos)
  util$plot_div_pairs(bspp, bspp, samples_eos, diagnostics_eos)
  
  # check bsppyr
  bsppyr <- paste0("bspyr[", 1:4, "]")
  util$plot_div_pairs(bsppyr, bsppyr, samples, diagnostics)

  
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
jpeg(filename = "figures/growthModelsMain/diagnostics/retrodictiveCheckHist.jpeg",
  width = 2400, height = 2400, res = 300)
util$plot_hist_quantiles(samples_gdd, "y_rep", 
                         -2, # lower x axis limit
                         4, # upper x axis limit
                         0.3, # binning
                         baseline_values = dgdd$y,
                         xlab = "Ring width (mm)")
dev.off()


# Hist by species
# GDD
jpeg(filename = "figures/growthModelsMain/diagnostics/retrodictiveHistGDD.jpeg",
  width = 3600, height = 2000, res = 300)
par(mfrow = c(2, dgdd$Nspp/2))
for (s in unique(dgdd$species)) { # s = 2
  idxs <- which(dgdd$species == s)
  samples_sub <- samples_gdd[grep(paste0("^y_rep\\[(", paste(idxs, collapse="|"), ")\\]$"), names(samples_gdd))]
  util$plot_hist_quantiles(samples_sub,
                           "y_rep",
                           -2,
                           4,
                           0.4,
                           baseline_values = dgdd$y[idxs],
                           xlab = "log(ring width)",
                           main = unique(emp$latbi[which(emp$spp_num == s)]))
}
dev.off()

# GSL
jpeg(filename = "figures/growthModelsMain/diagnostics/retrodictiveHistGSL.jpeg",
  width = 3600, height = 2000, res = 300)
par(mfrow = c(2, dgsl$Nspp/2))
for (s in unique(dgsl$species)) { # s = 2
  idxs <- which(dgsl$species == s)
  samples_sub <- samples_gsl[grep(paste0("^y_rep\\[(", paste(idxs, collapse="|"), ")\\]$"), names(samples_gsl))]
  util$plot_hist_quantiles(samples_sub,
                           "y_rep",
                           -2,
                           4,
                           0.3,
                           baseline_values = dgsl$y[idxs],
                           xlab = "log(ring width)",
                           main = unique(emp$latbi[which(emp$spp_num == s)]))
}
dev.off()


# SOS
jpeg(filename = "figures/growthModelsMain/diagnostics/retrodictiveHistSOS.jpeg",
  width = 3600, height = 2000, res = 300)
par(mfrow = c(2, dsos$Nspp/2))
for (s in unique(dsos$species)) { # s = 2
  idxs <- which(dsos$species == s)
  samples_sub <- samples_sos[grep(paste0("^y_rep\\[(", paste(idxs, collapse="|"), ")\\]$"), names(samples_sos))]
  util$plot_hist_quantiles(samples_sub,
                           "y_rep",
                           -2,
                           4,
                           0.3,
                           baseline_values = dsos$y[idxs],
                           xlab = "log(ring width)",
                           main = unique(emp$latbi[which(emp$spp_num == s)]))
}
dev.off()

# EOS
jpeg(filename = "figures/growthModelsMain/diagnostics/retrodictiveHistEOS.jpeg",
  width = 3600, height = 2000, res = 300)
par(mfrow = c(2, deos$Nspp/2))
for (s in unique(deos$species)) { # s = 2
  idxs <- which(deos$species == s)
  samples_sub <- samples_eos[grep(paste0("^y_rep\\[(", paste(idxs, collapse="|"), ")\\]$"), names(samples_eos))]
  util$plot_hist_quantiles(samples_sub,
                           "y_rep",
                           -2,
                           4,
                           0.3,
                           baseline_values = deos$y[idxs],
                           xlab = "log(ring width)",
                           main = unique(emp$latbi[which(emp$spp_num == s)]))
}
dev.off()

# Site GDD
jpeg(filename = "figures/growthModelsMain/diagnostics/retrodictiveHistGDDSite.jpeg",
     width = 3600, height = 2000, res = 300)
par(mfrow = c(2, dgdd$Nsite/2))
for (s in unique(dgdd$site)) { # s = 2
  idxs <- which(dgdd$site == s)
  samples_sub <- samples_gdd[grep(paste0("^y_rep\\[(", paste(idxs, collapse="|"), ")\\]$"), names(samples_gdd))]
  util$plot_hist_quantiles(samples_sub,
                           "y_rep",
                           -2,
                           4,
                           0.4,
                           baseline_values = dgdd$y[idxs],
                           xlab = "log(ring width)",
                           main = unique(emp$site[which(emp$site_num == s)]))
}
dev.off()
}


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FULL DATA ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (fitmodelfull) {
# Fit model SOS
# transform my groups to numeric values
empfullsos$site_num <- match(empfullsos$site, unique(empfullsos$site))
empfullsos$spp_num <- match(empfullsos$spp, unique(empfullsos$spp))
empfullsos$treeid_num <- match(empfullsos$treeid, unique(empfullsos$treeid))

# order by tree id
treeid_spp_site_sos <- unique(empfullsos[, c("treeid_num", "spp_num", "site_num",
                                  "treeid", "spp", "site", "latbi")])

treeid_spp_site_ordered_sos <- treeid_spp_site_sos[order(treeid_spp_site_sos$treeid_num), ]

# transform data in vectors for gsl
dsosfull <- list(
  y = empfullsos$loglength,
  N = nrow(empfullsos),
  Nspp = length(unique(empfullsos$spp_num)),
  Nsite = length(unique(empfullsos$site_num)),
  site = as.numeric(as.character(empfullsos$site_num)),
  species = as.numeric(as.character(empfullsos$spp_num)),
  treeid = as.numeric(empfullsos$treeid_num),
  Ntreeid = length(unique(empfullsos$treeid_num)),
  treeid_species = treeid_spp_site_ordered_sos$spp_num,
  treeid_site = treeid_spp_site_ordered_sos$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered_sos$spp_num)),
  sos = empfullsos$leafout / sosscale,
  sosseq = sosseq,
  sosscale = sosscale,
  Nsosseq = length(sosseq)
)

sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsosfull <- sampling(sosmodel, data = dsosfull,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsosfull, "output/stanOutput/fitGrowthSOSFull")

# Fit model EOS
# transform my groups to numeric values
empfulleos$site_num <- match(empfulleos$site, unique(empfulleos$site))
empfulleos$spp_num <- match(empfulleos$spp, unique(empfulleos$spp))
empfulleos$treeid_num <- match(empfulleos$treeid, unique(empfulleos$treeid))

# order by tree id
treeid_spp_site_eos <- unique(empfulleos[, c("treeid_num", "spp_num", "site_num",
                                             "treeid", "spp", "site", "latbi")])

treeid_spp_site_ordered_eos <- treeid_spp_site_eos[order(treeid_spp_site_eos$treeid_num), ]


# transform data in vectors for gsl
deosfull <- list(
  y = empfulleos$loglength,
  N = nrow(empfulleos),
  Nspp = length(unique(empfulleos$spp_num)),
  Nsite = length(unique(empfulleos$site_num)),
  site = as.numeric(as.character(empfulleos$site_num)),
  species = as.numeric(as.character(empfulleos$spp_num)),
  treeid = as.numeric(empfulleos$treeid_num),
  Ntreeid = length(unique(as.numeric(empfulleos$treeid_num))),
  treeid_species = treeid_spp_site_ordered_eos$spp_num,
  treeid_site = treeid_spp_site_ordered_eos$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered_eos$spp_num)),
  eos = empfulleos$budset / eosscale,
  eosseq = eosseq,
  eosscale = eosscale,
  Neosseq = length(eosseq)
)
eosmodel <- stan_model("stan/modelGrowthEOS.stan")
fiteosfull <- sampling(eosmodel, data = deosfull,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fiteosfull, "output/stanOutput/fitGrowthEOSFull")

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
       pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
points(sigma_df2_sos$mean, sigma_df2_full_sos$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
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
       pch = 16, col = "#0a6a3c", cex = 1.5)
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
       pch = 16, col = "#0a6a3c", cex = 1.5)
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
       pch = 16, col = "#d39822", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
points(sigma_df2_eos$mean, sigma_df2_full_eos$mean, pch = 16, col = "#d39822", cex = 1.5)
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
       pch = 16, col = "#d39822", cex = 1.5)
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
       pch = 16, col = "#d39822", cex = 1.5)
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
genericmodel <- stan_model("stan/modelGrowth_z.stan")

# Fit model GDD
gddz <- (emp$pgsGDD5 - mean(emp$pgsGDD5)) / sd(emp$pgsGDD5)
dgddz <- dgdd[1:10]
dgddz$covariate <- gddz

fitgdd <- sampling(genericmodel, data = dgddz,
                   warmup = 1000, iter=2000, chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDDZscored")

# Fit model GSL
gslz <- (emp$pgsGSL - mean(emp$pgsGSL)) / sd(emp$pgsGSL)
dgslz <- dgdd[1:10]
dgslz$covariate <- gslz

fitgsl <- sampling(genericmodel, data = dgslz,
                   warmup = 1000, iter = 2000, chains = 4)
saveRDS(fitgsl, "output/stanOutput/fitGrowthGSLZscored")

# Fit model SOS
sosz <- (emp$leafout - mean(emp$leafout)) / sd(emp$leafout)
dsosz <- dgdd[1:10]
dsosz$covariate <- sosz

fitsos <- sampling(genericmodel, data = dsosz,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOSZscored")

# Fit model EOS
eosz <- (emp$budset - mean(emp$budset)) / sd(emp$budset)
deosz <- dgdd[1:10]
deosz$covariate <- eosz

fiteos <- sampling(genericmodel, data = deosz,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fiteos, "output/stanOutput/fitGrowthEOSZscored")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Diagnostics #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
diagnostics <- util$extract_hmc_diagnostics(fitgsl) 
util$check_all_hmc_diagnostics(diagnostics)
samples <- util$extract_expectand_vals(fitgsl)

# asite
asite <- names(samples)[grepl("asite", names(samples))]
asite <- asite[!grepl("sigma|prior", asite)]
# pdf("figures/troubleShootingGrowthModel/atreeidParameterization_only_atreeid_no_b.pdf", 
#     width = 6, height = 18)
jpeg("figures/growthModelsMain/zscored/asiteParameterization.jpeg", 
     width = 2000, height = 3000,
     units = "px", res = 300)
util$plot_div_pairs(asite, "sigma_asite", samples, diagnostics, transforms = list("sigma_asite" = 1))
dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot GDD fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgdd <- readRDS("output/stanOutput/fitGrowthGDDZscored")
fitgsl <- readRDS("output/stanOutput/fitGrowthGSLZscored")
fitsos <- readRDS("output/stanOutput/fitGrowthSOSZscored")
fiteos <- readRDS("output/stanOutput/fitGrowthEOSZscored")

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
site_df2 <- subset(site_df2, !grepl("z|sigma", site))

# save csvs
write.csv(sigma_df2,  "output/GM_GDDparam_Z_sigma.csv",  row.names = FALSE)
write.csv(bspp_df2,   "output/GM_GDDparam_Z_bspp.csv",   row.names = FALSE)
write.csv(treeid_df2, "output/GM_GDDparam_Z_treeid.csv", row.names = FALSE)
write.csv(aspp_df2,   "output/GM_GDDparam_Z_aspp.csv",   row.names = FALSE)
write.csv(site_df2,   "output/GM_GDDparam_Z_site.csv",   row.names = FALSE)


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
site_df2 <- subset(site_df2, !grepl("z|sigma", site))

# save csvs
write.csv(sigma_df2,  "output/GM_GSLparam_Z_sigma.csv",  row.names = FALSE)
write.csv(bspp_df2,   "output/GM_GSLparam_Z_bspp.csv",   row.names = FALSE)
write.csv(treeid_df2, "output/GM_GSLparam_Z_treeid.csv", row.names = FALSE)
write.csv(aspp_df2,   "output/GM_GSLparam_Z_aspp.csv",   row.names = FALSE)
write.csv(site_df2,   "output/GM_GSLparam_Z_site.csv",   row.names = FALSE)

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
site_df2_sos <- subset(site_df2_sos, !grepl("z|sigma", site))

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
site_df2_eos <- subset(site_df2_eos, !grepl("z|sigma", site))

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




# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# asite partial pooling comparison ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (fitmodels) {
emp$year_num <- match(emp$year, unique(emp$year))
dgdd$year <- as.numeric(as.character(emp$year_num))
dgdd$Nyear <- length(unique(dgdd$year))

gddmodelpp <- stan_model("stan/modelGrowthGDD_PPsite.stan")
fitgddppsite <- sampling(gddmodelpp, data = dgdd,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitgddppsite, "output/stanOutput/fitGrowthGDD_PPsite")
fitgddppsite <- readRDS("output/stanOutput/fitGrowthGDD_PPsite")

##### Recover parameters #####
df_fitgddpp <- as.data.frame(fitgddppsite)

# full posterior
columns <- colnames(df_fitgddpp)[!grepl("prior", colnames(df_fitgddpp))]
sigma_df <- df_fitgddpp[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgddpp[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgddpp[, grepl("treeid", columns) & !grepl("z|sigma|slope|full", columns)]
aspp_df <- df_fitgddpp[, columns[grepl("aspp", columns)]]
site_df <- df_fitgddpp[, columns[grepl("asite", columns)]]
ayear_df <- df_fitgddpp[, columns[grepl("ayear", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgddpp, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgddpp, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgddpp, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgddpp, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fitgddpp, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")
site_df2 <- subset(site_df2, !grepl("z|sigma", site))
ayear_df2   <- extract_params(df_fitgddpp, "ayear", "fit_ayear", 
                             "year", "ayear\\[(\\d+)\\]")
##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gddModelPriorVSPosterior_PPsite.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 3))

# a
plot(density(df_fitgddpp[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitgddpp[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgddpp[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitgddpp[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgddpp[, "sigma_asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_asite", 
     xlab = "sigma_asite", ylim = c(0,4))
lines(density(df_fitgddpp[, "sigma_asite"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgddpp[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitgddpp[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgddpp[, "aspp_prior"]), 
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
plot(density(df_fitgddpp[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgddpp[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

plot(density(df_fitgddpp[, "atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     # main = "priorVSposterior_asite", 
     xlab = "atreeid", xlim = c(-6, 6), ylim = c(0, 5))
for (col in colnames(treeid_df)) {
  lines(density(treeid_df[, col]), col = pal[2], lwd = 0.4)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# Recover fitgdd without partial pooling
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD")

##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df_noPP <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df_noPP <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df_noPP <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma|slope|full", columns)]
aspp_df_noPP <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df_noPP <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df_noPP) <- 1:ncol(bspp_df_noPP)
colnames(treeid_df_noPP) <- 1:ncol(treeid_df_noPP)
colnames(aspp_df_noPP) <- 1:ncol(aspp_df_noPP)
colnames(site_df_noPP) <- 1:ncol(site_df_noPP)

sigma_df2_noPP  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2_noPP   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2_noPP <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_noPP <- subset(treeid_df2_noPP, !grepl("z|sigma", treeid))
aspp_df2_noPP   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2_noPP   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


##### Compare model output with and without partial pooling #####
# Open device
jpeg("figures/growthModelsMain/sitePPvsNoPP.jpeg", width = 9, height = 6, units = "in", res = 300)
par(mfrow = c(2,3), oma = c(0, 2, 0, 0))

# sigma
sigma_df3 <- subset(sigma_df2, sigma != "sigma_asite")
plot(sigma_df2_noPP$mean, sigma_df3$mean,
     xlab = "no partial pooling on asite", ylab = "partial pooling on asite", main = "sigmas", type = "n", frame = FALSE,
     ylim = range(c(sigma_df3$p25, sigma_df3$p75)),
     xlim = range(c(sigma_df2_noPP$p25, sigma_df2_noPP$p75+0.2)))
arrows(x0 = sigma_df2_noPP$mean, y0 = sigma_df3$p25,
       x1 = sigma_df2_noPP$mean, y1 = sigma_df3$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_noPP$p25, y0 = sigma_df3$mean,
       x1 = sigma_df2_noPP$p75, y1 = sigma_df3$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_noPP$mean, sigma_df3$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
points(sigma_df2_noPP$mean, sigma_df3$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
text(sigma_df2_noPP$p75, sigma_df3$p25, labels = sigma_df2_noPP$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2_noPP$mean, bspp_df2$mean,
     xlab = "no partial pooling on asite", ylab = "partial pooling on asite", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2$p25, bspp_df2$p75)),
     xlim = range(c(bspp_df2_noPP$p25, bspp_df2_noPP$p75)))
arrows(x0 = bspp_df2_noPP$mean, y0 = bspp_df2$p25,
       x1 = bspp_df2_noPP$mean, y1 = bspp_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2_noPP$p25, y0 = bspp_df2$mean,
       x1 = bspp_df2_noPP$p75, y1 = bspp_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2_noPP$mean, bspp_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# aspp
plot(aspp_df2_noPP$mean, aspp_df2$mean,
     xlab = "no partial pooling on asite", ylab = "partial pooling on asite", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2$p25, aspp_df2$p75)),
     xlim = range(c(aspp_df2_noPP$p25, aspp_df2_noPP$p75)))
arrows(x0 = aspp_df2_noPP$mean, y0 = aspp_df2$p25,
       x1 = aspp_df2_noPP$mean, y1 = aspp_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2_noPP$p25, y0 = aspp_df2$mean,
       x1 = aspp_df2_noPP$p75, y1 = aspp_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2_noPP$mean, aspp_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# asite
plot(site_df2_noPP$mean, site_df2$mean,
     xlab = "no partial pooling on asite", ylab = "partial pooling on asite", main = "asite", type = "n", frame = FALSE,
     ylim = range(c(site_df2$p25, site_df2$p75)),
     xlim = range(c(site_df2_noPP$p25, site_df2_noPP$p75)))
arrows(x0 = site_df2_noPP$mean, y0 = site_df2$p25,
       x1 = site_df2_noPP$mean, y1 = site_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = site_df2_noPP$p25, y0 = site_df2$mean,
       x1 = site_df2_noPP$p75, y1 = site_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(site_df2_noPP$mean, site_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)


# atreeid
plot(treeid_df2_noPP$mean, treeid_df2$mean,
     xlab = "no partial pooling on asite", ylab = "partial pooling on asite", main = "atreeid", type = "n", frame = FALSE,
     ylim = range(c(treeid_df2$p25, treeid_df2$p75)),
     xlim = range(c(treeid_df2_noPP$p25, treeid_df2_noPP$p75)))
arrows(x0 = treeid_df2_noPP$mean, y0 = treeid_df2$p25,
       x1 = treeid_df2_noPP$mean, y1 = treeid_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1, col = "darkgray")
arrows(x0 = treeid_df2_noPP$p25, y0 = treeid_df2$mean,
       x1 = treeid_df2_noPP$p75, y1 = treeid_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1, col = "darkgray")
points(treeid_df2_noPP$mean, treeid_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
dev.off()

}
 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Add ayear to model ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (fitmodels) {
emp$year_num <- match(emp$year, unique(emp$year))
dgdd$year <- as.numeric(as.character(emp$year_num))
dgdd$Nyear <- length(unique(dgdd$year))

gddmodelayear <- stan_model("stan/modelGrowthGDD_ayear.stan")
fitgddayear <- sampling(gddmodelayear, data = dgdd,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitgddayear, "output/stanOutput/fitGrowthGDD_ayear")
# fitgddayear <- readRDS("output/stanOutput/fitGrowthGDD_ayear")


##### Recover parameters #####
df_fitgddayear <- as.data.frame(fitgddayear)

# full posterior
columns <- colnames(df_fitgddayear)[!grepl("prior", colnames(df_fitgddayear))]
sigma_df <- df_fitgddayear[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgddayear[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgddayear[, grepl("treeid", columns) & !grepl("z|sigma|slope|full", columns)]
aspp_df <- df_fitgddayear[, columns[grepl("aspp", columns)]]
site_df <- df_fitgddayear[, columns[grepl("asite", columns)]]
ayear_df <- df_fitgddayear[, columns[grepl("year", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)
colnames(ayear_df) <- 1:ncol(ayear_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgddayear, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgddayear, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgddayear, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgddayear, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fitgddayear, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")
ayear_df2   <- extract_params(df_fitgddayear, "ayear", "fit_a_year", 
                             "year", "ayear\\[(\\d+)\\]")

##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gddModelPriorVSPosterior_ayear.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 3))

# a
plot(density(df_fitgddayear[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0, 1))
lines(density(df_fitgddayear[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fitgddayear[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,4))
lines(density(df_fitgddayear[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fitgddayear[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0, 4))
lines(density(df_fitgddayear[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fitgddayear[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fitgddayear[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", 
     # xlim = c(-5, 5), 
     ylim = c(0, 1))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fitgddayear[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 1))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fitgddayear[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 5))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

plot(density(df_fitgddayear[, "atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     # main = "priorVSposterior_asite", 
     xlab = "atreeid", xlim = c(-6, 6), ylim = c(0, 5))
for (col in colnames(treeid_df)) {
  lines(density(treeid_df[, col]), col = pal[2], lwd = 0.4)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# Recover fitgdd without partial pooling
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD")

##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df_noayr <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df_noayr <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df_noayr <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma|slope|full", columns)]
aspp_df_noayr <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df_noayr <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df_noayr) <- 1:ncol(bspp_df_noayr)
colnames(treeid_df_noayr) <- 1:ncol(treeid_df_noayr)
colnames(aspp_df_noayr) <- 1:ncol(aspp_df_noayr)
colnames(site_df_noayr) <- 1:ncol(site_df_noayr)

sigma_df2_noayr  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2_noayr   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2_noayr <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_noayr <- subset(treeid_df2_noayr, !grepl("z|sigma", treeid))
aspp_df2_noayr   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2_noayr   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")

# quick mu plot for bspp
par(mfrow = c(2,1))
y_pos <- rev(1:4)
plot(bspp_df2$mean, y_pos,
     xlim = c(-0.5, 0.6), ylim = c(0.5, 4 + 0.5),
     xlab = "log(ring width) change per 10 spring days GDD", ylab = "",
     main = "with ayear",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2$p5,  y_pos, bspp_df2$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2$p25, y_pos, bspp_df2$p75, y_pos, col = wccolslatbi, lwd = 3)


plot(bspp_df2_noayr$mean, y_pos,
     xlim = c(-0.5, 0.6), ylim = c(0.5, 4 + 0.5),
     xlab = "log(ring width) change per 10 spring days GDD", ylab = "",
     main = "without ayear",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_noayr$p5,  y_pos, bspp_df2_noayr$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_noayr$p25, y_pos, bspp_df2_noayr$p75, y_pos, col = wccolslatbi, lwd = 3)

##### Compare model output with and without partial pooling #####
# Open device
jpeg("figures/growthModelsMain/sitePPvsnoayr.jpeg", width = 9, height = 6, units = "in", res = 300)
par(mfrow = c(2,3), oma = c(0, 2, 0, 0))

# sigma
plot(sigma_df2_noayr$mean, sigma_df2$mean,
     xlab = "no year intercept", ylab = "with year intercept", main = "sigmas", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2$p25, sigma_df2$p75)),
     xlim = range(c(sigma_df2_noayr$p25, sigma_df2_noayr$p75+0.2)))
arrows(x0 = sigma_df2_noayr$mean, y0 = sigma_df2$p25,
       x1 = sigma_df2_noayr$mean, y1 = sigma_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_noayr$p25, y0 = sigma_df2$mean,
       x1 = sigma_df2_noayr$p75, y1 = sigma_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_noayr$mean, sigma_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
points(sigma_df2_noayr$mean, sigma_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
text(sigma_df2_noayr$p75, sigma_df2$p25, labels = sigma_df2_noayr$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2_noayr$mean, bspp_df2$mean,
     xlab = "no year intercept", ylab = "with year intercept", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2$p25, bspp_df2$p75)),
     xlim = range(c(bspp_df2_noayr$p25, bspp_df2_noayr$p75)))
arrows(x0 = bspp_df2_noayr$mean, y0 = bspp_df2$p25,
       x1 = bspp_df2_noayr$mean, y1 = bspp_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2_noayr$p25, y0 = bspp_df2$mean,
       x1 = bspp_df2_noayr$p75, y1 = bspp_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2_noayr$mean, bspp_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# aspp
plot(aspp_df2_noayr$mean, aspp_df2$mean,
     xlab = "no year intercept", ylab = "with year intercept", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2$p25, aspp_df2$p75)),
     xlim = range(c(aspp_df2_noayr$p25, aspp_df2_noayr$p75)))
arrows(x0 = aspp_df2_noayr$mean, y0 = aspp_df2$p25,
       x1 = aspp_df2_noayr$mean, y1 = aspp_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2_noayr$p25, y0 = aspp_df2$mean,
       x1 = aspp_df2_noayr$p75, y1 = aspp_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2_noayr$mean, aspp_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# asite
plot(site_df2_noayr$mean, site_df2$mean,
     xlab = "no year intercept", ylab = "with year intercept", main = "asite", type = "n", frame = FALSE,
     ylim = range(c(site_df2$p25, site_df2$p75)),
     xlim = range(c(site_df2_noayr$p25, site_df2_noayr$p75)))
arrows(x0 = site_df2_noayr$mean, y0 = site_df2$p25,
       x1 = site_df2_noayr$mean, y1 = site_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = site_df2_noayr$p25, y0 = site_df2$mean,
       x1 = site_df2_noayr$p75, y1 = site_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(site_df2_noayr$mean, site_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)


# atreeid
plot(treeid_df2_noayr$mean, treeid_df2$mean,
     xlab = "no year intercept", ylab = "with year intercept", main = "atreeid", type = "n", frame = FALSE,
     ylim = range(c(treeid_df2$p25, treeid_df2$p75)),
     xlim = range(c(treeid_df2_noayr$p25, treeid_df2_noayr$p75)))
arrows(x0 = treeid_df2_noayr$mean, y0 = treeid_df2$p25,
       x1 = treeid_df2_noayr$mean, y1 = treeid_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1, col = "darkgray")
arrows(x0 = treeid_df2_noayr$p25, y0 = treeid_df2$mean,
       x1 = treeid_df2_noayr$p75, y1 = treeid_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1, col = "darkgray")
points(treeid_df2_noayr$mean, treeid_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
dev.off()

# plots years
colsyear <- c("#dd5129", "#0f7ba2", "#43b284")
jpeg("figures/growthModelsMain/muayear.jpeg", width = 6, height = 6, 
     units = "in", res = 300)
par(mfrow = c(1,1))
ayear_df2$year_name <- emp$year[match(ayear_df2$year, emp$year_num)]
n_year <- length(unique(ayear_df2$year))
y_pos <- 1:n_year

plot(ayear_df2$mean, y_pos,
     xlim = c(-5, 5), ylim = c(0.5, n_year + 0.5), 
     xlab = "year intercept", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = colsyear, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(ayear_df2$p5, y_pos, ayear_df2$p95, y_pos,
         col = colsyear, lwd = 1.5)
segments(ayear_df2$p25, y_pos, ayear_df2$p75, y_pos,
         col = colsyear, lwd = 3)
axis(2, at = y_pos, labels = ayear_df2$year_name, las = 1)
dev.off()

}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Look at GDD > 30 degrees celsius ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (FALSE){
subset(gddyr, meanTempC >30)[, 1:4]
check <- subset(gddyr, maxTempC >30)
nrow(check)
gddyr$max30 <- ifelse(gddyr$maxTempC > 30, 1, 0)
gddyr$yrdoy <- paste(gddyr$year, gddyr$doy, sep = "_")

emp$yrdoy <- paste(emp$year, emp$doy, sep = "_")

emp$idyr <- paste(emp$treeid, emp$year, sep = "_")
emp$idyr_num <- match(emp$idyr, unique(emp$idyr))

idvec <- unique(emp$idyr_num)

emp$n30 <- NA
for (i in seq_len(nrow(emp))) { # i = 5
  e <- emp[emp$idyr_num == i,]
  g <- gddyr[which(gddyr$year == e$year & 
                     gddyr$doy >= e$leafout &
                     gddyr$doy <= e$budset),]
  sum <- sum(g$max30 == 1)
  emp$n30[i] <- sum
}

dgdd$gddabv <- as.numeric(emp$n30) / 4

# Fit model GDD with slope on temp above 30 C
gddmodelabv <- stan_model("stan/modelGrowthGDD_bMax30.stan")
fitgdd <- sampling(gddmodelabv, data = dgdd,
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDD_abv30")
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD_abv30")

##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgdd[, columns[grepl("bsp", columns)]]
bsppabv_df <- df_fitgdd[, columns[grepl("bspabv", columns)]]
treeid_df <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(bsppabv_df) <- 1:ncol(bsppabv_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
bspp_df2 <- subset(bspp_df2, !grepl("bsp", spp))
bsppabv_df2   <- extract_params(df_fitgdd, "bspabv", "fit_bspp", 
                             "spp", "bspabv\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/growthModelsMain/diagnostics/gddModelPriorVSPosterior_abv30.pdf", width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 3))

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

# bspabv
plot(density(df_fitgdd[, "bspabv_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bspabv30", 
     xlab = "bspabv30", ylim = c(0, 5))
for (col in colnames(bsppabv_df)) {
  lines(density(bsppabv_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

par(mfrow = c(1,1))
y_pos <- rev(1:4)
n_spp <- 4
plot(bsppabv_df2$mean, y_pos,
     xlim = c(-0.5, 0.6), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change per 4 days above 30C", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bsppabv_df2$p5,  y_pos, bsppabv_df2$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bsppabv_df2$p25, y_pos, bsppabv_df2$p75, y_pos, col = wccolslatbi, lwd = 3)

# Plot empirical data

jpeg(file = "figures/empiricalData/nMax30VariationSppYr.jpeg",
     width = 2400, height = 2000, res = 300)

# mu plot dimensions and stuff
species_order <- c(
  "Alnus incana", 
  "Betula alleghaniensis", 
  "Betula papyrifera", 
  "Betula populifolia")

gap <- 2
years <- c(2018, 2019, 2020)
n_sp <- length(species_order)


aggabv30 <- aggregate(n30 ~ latbi + year, emp, FUN = mean)
quants <- aggregate(n30 ~ latbi + year, emp, FUN = function(x) quantile(x, c(0.05, 0.25, 0.75, 0.95)))
aggabv30$p5  <- quants$n30[,1]
aggabv30$p25 <- quants$n30[,2]
aggabv30$p75 <- quants$n30[,3]
aggabv30$p95 <- quants$n30[,4]

total_rows <- nrow(aggabv30) + (length(years) - 1) * gap

current_y <- total_rows
aggabv30$y_pos <- NA

aggabv30 <- aggabv30[order(aggabv30$latbi),]
aggabv30$spp_num <- as.integer(as.factor(aggabv30$latbi))

for(s in unique(aggabv30$spp_num)){ # s = 2
  idx <- which(aggabv30$spp_num == s)
  aggabv30$y_pos[idx] <- current_y:(current_y - length(idx) + 1)
  current_y <- current_y - length(idx) - gap
}
aggabv30
par(mar = c(4,6,4,2))

plot(aggabv30$n30, aggabv30$y_pos,
     xlim = c(25, 60), 
     ylim = c(min(aggabv30$y_pos), max(aggabv30$y_pos) + 0.5),
     xlab = "anomalized gsl (days)", ylab = "",
     yaxt = "n",
     pch = 16, cex = 2, col = wccolslatbi[aggabv30$latbi], frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aggabv30$p5,  aggabv30$y_pos, aggabv30$p95, aggabv30$y_pos, 
         col = wccolslatbi[aggabv30$latbi], lwd = 1.5)
segments(aggabv30$p25, aggabv30$y_pos, aggabv30$p75, aggabv30$y_pos, 
         col = wccolslatbi[aggabv30$latbi], lwd = 3)
abline(v = 0, lty = 2)

# custom y axis label
ylabel <- aggregate(y_pos ~ latbi + year, aggabv30, mean)
ylabel
aggabv30$ylabel <- ylabel$y_pos[match(aggabv30$year, ylabel$year)]
axis(
  side = 2,
  at = aggabv30$y_pos,
  labels = aggabv30$year,
  cex.axis = 1,
  las = 1
)

# # add n per year and species
# sum <- aggregate(anomgsl ~ latbi + year, emp, function(x) length(x))
# aggabv30$count <- sum$anomgsl

legend("right",
       legend = sapply(unique(aggabv30$latbi), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = wccolslatbi,
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)
dev.off()                 


##### Compare model output with and without slope on max temp #####
# Recover fitgdd without partial pooling
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD")

##### Recover parameters #####
df_fitgdd <- as.data.frame(fitgdd)

# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df_noayr <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df_nobspabv <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df_nobspabv <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma|slope|full", columns)]
aspp_df_nobspabv <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df_nobspabv <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df_nobspabv) <- 1:ncol(bspp_df_nobspabv)
colnames(treeid_df_nobspabv) <- 1:ncol(treeid_df_nobspabv)
colnames(aspp_df_nobspabv) <- 1:ncol(aspp_df_nobspabv)
colnames(site_df_nobspabv) <- 1:ncol(site_df_nobspabv)

sigma_df2_nobspabv  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2_nobspabv   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                                   "spp", "bsp\\[(\\d+)\\]")
treeid_df2_nobspabv <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                                   "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_nobspabv <- subset(treeid_df2_nobspabv, !grepl("z|sigma", treeid))
aspp_df2_nobspabv   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                                   "spp", "aspp\\[(\\d+)\\]")
site_df2_nobspabv   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                                   "site", "asite\\[(\\d+)\\]")

# Open device
jpeg("figures/growthModelsMain/bsppVSbsppabv.jpeg", width = 9, height = 6, units = "in", res = 300)
par(mfrow = c(2,3), oma = c(0, 2, 0, 0))

# sigma
plot(sigma_df2_nobspabv$mean, sigma_df2$mean,
     xlab = "no bsp on max temp", ylab = "with bspp on max temp", main = "sigmas", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2$p25, sigma_df2$p75)),
     xlim = range(c(sigma_df2_nobspabv$p25, sigma_df2_nobspabv$p75+0.2)))
arrows(x0 = sigma_df2_nobspabv$mean, y0 = sigma_df2$p25,
       x1 = sigma_df2_nobspabv$mean, y1 = sigma_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2_nobspabv$p25, y0 = sigma_df2$mean,
       x1 = sigma_df2_nobspabv$p75, y1 = sigma_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(sigma_df2_nobspabv$mean, sigma_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
points(sigma_df2_nobspabv$mean, sigma_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
text(sigma_df2_nobspabv$p75, sigma_df2$p25, labels = sigma_df2_nobspabv$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2_nobspabv$mean, bspp_df2$mean,
     xlab = "no bsp on max temp", ylab = "with bspp on max temp", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2$p25, bspp_df2$p75)),
     xlim = range(c(bspp_df2_nobspabv$p25, bspp_df2_nobspabv$p75)))
arrows(x0 = bspp_df2_nobspabv$mean, y0 = bspp_df2$p25,
       x1 = bspp_df2_nobspabv$mean, y1 = bspp_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2_nobspabv$p25, y0 = bspp_df2$mean,
       x1 = bspp_df2_nobspabv$p75, y1 = bspp_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2_nobspabv$mean, bspp_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# aspp
plot(aspp_df2_nobspabv$mean, aspp_df2$mean,
     xlab = "no bsp on max temp", ylab = "with bspp on max temp", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2$p25, aspp_df2$p75)),
     xlim = range(c(aspp_df2_nobspabv$p25, aspp_df2_nobspabv$p75)))
arrows(x0 = aspp_df2_nobspabv$mean, y0 = aspp_df2$p25,
       x1 = aspp_df2_nobspabv$mean, y1 = aspp_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2_nobspabv$p25, y0 = aspp_df2$mean,
       x1 = aspp_df2_nobspabv$p75, y1 = aspp_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2_nobspabv$mean, aspp_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)

# asite
plot(site_df2_nobspabv$mean, site_df2$mean,
     xlab = "no bsp on max temp", ylab = "with bspp on max temp", main = "asite", type = "n", frame = FALSE,
     ylim = range(c(site_df2$p25, site_df2$p75)),
     xlim = range(c(site_df2_nobspabv$p25, site_df2_nobspabv$p75)))
arrows(x0 = site_df2_nobspabv$mean, y0 = site_df2$p25,
       x1 = site_df2_nobspabv$mean, y1 = site_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = site_df2_nobspabv$p25, y0 = site_df2$mean,
       x1 = site_df2_nobspabv$p75, y1 = site_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(site_df2_nobspabv$mean, site_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)


# atreeid
plot(treeid_df2_nobspabv$mean, treeid_df2$mean,
     xlab = "no bsp on max temp", ylab = "with bspp on max temp", main = "atreeid", type = "n", frame = FALSE,
     ylim = range(c(treeid_df2$p25, treeid_df2$p75)),
     xlim = range(c(treeid_df2_nobspabv$p25, treeid_df2_nobspabv$p75)))
arrows(x0 = treeid_df2_nobspabv$mean, y0 = treeid_df2$p25,
       x1 = treeid_df2_nobspabv$mean, y1 = treeid_df2$p75,
       angle = 90, code = 3, length = 0, lwd = 1, col = "darkgray")
arrows(x0 = treeid_df2_nobspabv$p25, y0 = treeid_df2$mean,
       x1 = treeid_df2_nobspabv$p75, y1 = treeid_df2$mean,
       angle = 90, code = 3, length = 0, lwd = 1, col = "darkgray")
points(treeid_df2_nobspabv$mean, treeid_df2$mean, pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
dev.off()
}