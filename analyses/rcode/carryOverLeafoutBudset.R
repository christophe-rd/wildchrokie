# Wildchrokie looking at carryover effect of leafout on budset
# CRD 26 March 2026
# look at evidence of a leafout carryover with budset. 

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

emp <- read.csv("output/empiricalDataMAIN.csv")

# Fit model GDD
emp2 <- emp[!is.na(emp$pgsGDD5),]

# transform my groups to numeric values
emp2$site_num <- match(emp2$site, unique(emp2$site))
emp2$spp_num <- match(emp2$spp, unique(emp2$spp))
emp2$treeid_num <- match(emp2$treeid, unique(emp2$treeid))

data <- list(
  y = emp2$budset / 10,
  N = nrow(emp2),
  Nspp = length(unique(emp2$spp_num)),
  Nsite = length(unique(emp2$site_num)),
  site = as.numeric(as.character(emp2$site_num)),
  species = as.numeric(as.character(emp2$spp_num)),
  treeid = as.numeric(emp2$treeid_num),
  Ntreeid = length(unique(as.numeric(emp2$treeid_num))),
  leafout = emp2$leafout / 10
)
# Fit model GDD
rstan_options(auto_write = TRUE)
comodel <- stan_model("stan/carryOverLeafoutBudset.stan")
fit <- sampling(comodel, data = data,
                warmup = 1000, iter = 2000, chains=4)
saveRDS(fit, "output/stanOutput/fitCarryOver")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##### Recover parameters #####
df_fit <- as.data.frame(fit)

# full posterior
columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
treeid_df <- df_fit[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
site_df <- df_fit[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)

# posterior summaries
sigma_df2  <- extract_params(df_fit, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fit, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fit, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fit, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2   <- extract_params(df_fit, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


##### Plot posterior vs priors for gdd fit #####
pdf(file = "figures/carryOverModel/carryOverModelPriorVSPosterior.pdf", 
    width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))

# a
plot(density(df_fit[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fit[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_atreeid
plot(density(df_fit[, "sigma_atreeid_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_atreeid", 
     xlab = "sigma_atreeid", ylim = c(0,2))
lines(density(df_fit[, "sigma_atreeid"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fit[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fit[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fit[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fit[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-6, 6), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fit[, "bsp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Retrodictive checks ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if (FALSE){
  
samples <- util$extract_expectand_vals(fit)
jpeg(
  filename = "figures/carryOverModel/retrodictiveCheckHist.jpeg",
  width = 2400,      
  height = 2400,
  res = 300          
)
util$plot_hist_quantiles(samples, "y_rep", 
                         20, # lower x axis limit
                         40, # upper x axis limit
                         0.5, # binning
                         baseline_values = y,
                         xlab = "Budset dates")
dev.off()

for (i in unique(data$treeid))
t <- 60
idxs <- which(data$treeid==t)
util$plot_disc_pushforward_quantiles(samples, 
                                     paste0("y_rep[",idxs,"]"),
                                     baseline_values = data$y[idxs], ylab = "Budset", xticklabs = 1:length(idxs))
}
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
sos <- empsos$leafout/10

rstan_options(auto_write = TRUE)
sosmodel <- stan_model("stan/modelGrowthSOS.stan")
fitsosfull <- sampling(sosmodel, data = c("N","y",
                                      "Nspp","species",
                                      "Nsite", "site",
                                      "Ntreeid", "treeid",
                                      "sos"),
                   warmup = 1000, iter = 2000, chains=4)
saveRDS(fitsos, "output/stanOutput/fitGrowthSOSFull")

# Fit model EOS --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
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
# SOS restricted
df_fitsos <- as.data.frame(fitsos)

# posterior summaries
sigma_df2  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# SOS full
df_fitsos <- as.data.frame(fitsosfull)

# posterior summaries
sigma_df2_full  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_full   <- extract_params(df_fitsos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_full <- extract_params(df_fitsos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_full <- subset(treeid_df2_full, !grepl("z|sigma", treeid))
aspp_df2_full   <- extract_params(df_fitsos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_full <- subset(treeid_df2_full, !grepl("prior", treeid))
site_df2_full   <- extract_params(df_fitsos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# Open device
jpeg("figures/empiricalData/SOSFullVSRestricted.jpeg", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))

plot(sigma_df2$mean, sigma_df2_full$mean,
     xlab = "restricted", ylab = "full", main = "", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2_full$mean_per25, sigma_df2_full$mean_per75)),
     xlim = range(c(sigma_df2$mean_per25, sigma_df2$mean_per75)))
arrows(x0 = sigma_df2$mean, y0 = sigma_df2_full$mean_per25,
       x1 = sigma_df2$mean, y1 = sigma_df2_full$mean_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2$mean_per25, y0 = sigma_df2_full$mean,
       x1 = sigma_df2$mean_per75, y1 = sigma_df2_full$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")

points(sigma_df2$mean, sigma_df2_full$mean,
       pch = 16, col = "#046C9A", cex = 1.5)

abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)
points(sigma_df2$mean, sigma_df2_full$mean, pch = 16, col = "#046C9A", cex = 1.5)
text(sigma_df2$mean_per75, sigma_df2_full$mean_per25, labels = sigma_df2$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2$fit_bspp, bspp_df2_full$fit_bspp,
     xlab = "restricted", ylab = "full", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2_full$fit_bspp_per25, bspp_df2_full$fit_bspp_per75)),
     xlim = range(c(bspp_df2$fit_bspp_per25, bspp_df2$fit_bspp_per75)))
arrows(x0 = bspp_df2$fit_bspp, y0 = bspp_df2_full$fit_bspp_per25,
       x1 = bspp_df2$fit_bspp, y1 = bspp_df2_full$fit_bspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2$fit_bspp_per25, y0 = bspp_df2_full$fit_bspp,
       x1 = bspp_df2$fit_bspp_per75, y1 = bspp_df2_full$fit_bspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2$fit_bspp, bspp_df2_full$fit_bspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# aspp
plot(aspp_df2$fit_aspp, aspp_df2_full$fit_aspp,
     xlab = "restricted", ylab = "full", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2_full$fit_aspp_per25, aspp_df2_full$fit_aspp_per75)),
     xlim = range(c(aspp_df2$fit_aspp_per25, aspp_df2$fit_aspp_per75)))
arrows(x0 = aspp_df2$fit_aspp, y0 = aspp_df2_full$fit_aspp_per25,
       x1 = aspp_df2$fit_aspp, y1 = aspp_df2_full$fit_aspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2$fit_aspp_per25, y0 = aspp_df2_full$fit_aspp,
       x1 = aspp_df2$fit_aspp_per75, y1 = aspp_df2_full$fit_aspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2$fit_aspp, aspp_df2_full$fit_aspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

###### Plot treeid ######
treeid_df2_full <- subset(treeid_df2_full, !treeid %in% setdiff(treeid_df2_full$treeid, treeid_df2$treeid))
plot(treeid_df2$fit_atreeid, treeid_df2_full$fit_atreeid,
     xlab = "restricted", ylab = "full", main = "atreeid", type = "n", frame = FALSE,
     ylim = range(c(treeid_df2_full$fit_atreeid_per25, treeid_df2_full$fit_atreeid_per75)),
     xlim = range(c(treeid_df2$fit_atreeid_per25, treeid_df2$fit_atreeid_per75)))
arrows(x0 = treeid_df2$fit_atreeid, y0 = treeid_df2_full$fit_atreeid_per25,
       x1 = treeid_df2$fit_atreeid, y1 = treeid_df2_full$fit_atreeid_per75,
       angle = 90, code = 3, length = 0, lwd = 0.8, col = "darkgray")
arrows(x0 = treeid_df2$fit_atreeid_per25, y0 = treeid_df2_full$fit_atreeid,
       x1 = treeid_df2$fit_atreeid_per75, y1 = treeid_df2_full$fit_atreeid,
       angle = 90, code = 3, length = 0, lwd = 0.8, col = "darkgray")
points(treeid_df2$fit_atreeid, treeid_df2_full$fit_atreeid,
       pch = 16, col = "#046C9A", cex = 0.8)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

dev.off()

# check which tree ids differ and if it's the result 
treeid_df2$fit_atreeid_full <- treeid_df2_full$fit_atreeid[match(treeid_df2$treeid, treeid_df2_full$treeid)]

treeid_df2$meandiff <- treeid_df2$fit_atreeid - treeid_df2$fit_atreeid_full
min(treeid_df2$meandiff)
max(treeid_df2$meandiff)
treeiddiffs <- c(head(treeid_df2[order(treeid_df2$meandiff), ], 3)$treeid, 
                 head(treeid_df2[order(-treeid_df2$meandiff), ], 3)$treeid)
fulleosleafout <- aggregate(leafout ~ treeid_num, empsos, FUN = length)
restreosleafout <- aggregate(leafout ~ treeid_num, emp2, FUN = length)
# checks
sum(fulleosleafout$leafout)
sum(restreosleafout$leafout)

fulleosleafout$restr <- restreosleafout$leafout[match(fulleosleafout$treeid_num, 
                                                      restreosleafout$treeid_num)]

checkleafout <- subset(fulleosleafout, treeid_num %in% treeiddiffs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Recover and plot parameters EOS restricted vs full #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# EOS restricted
df_fiteos <- as.data.frame(fiteos)

# posterior summaries
sigma_df2  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2 <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2 <- subset(treeid_df2, !grepl("prior", treeid))
site_df2   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# EOS full
df_fiteos <- as.data.frame(fiteosfull)

# posterior summaries
sigma_df2_full  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_full   <- extract_params(df_fiteos, "bsp", "fit_bspp", "spp", "bsp\\[(\\d+)\\]")
treeid_df2_full <- extract_params(df_fiteos, "atreeid", "fit_atreeid", "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_full <- subset(treeid_df2_full, !grepl("z|sigma", treeid))
aspp_df2_full   <- extract_params(df_fiteos, "aspp", "fit_aspp", "spp", "aspp\\[(\\d+)\\]")
treeid_df2_full <- subset(treeid_df2_full, !grepl("prior", treeid))
site_df2_full   <- extract_params(df_fiteos, "asite", "fit_a_site", "site", "asite\\[(\\d+)\\]")

# Open device
jpeg("figures/empiricalData/EOSFullVSRestricted.jpeg", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))

plot(sigma_df2$mean, sigma_df2_full$mean,
     xlab = "restricted", ylab = "full", main = "", type = "n", frame = FALSE,
     ylim = range(c(sigma_df2_full$mean_per25, sigma_df2_full$mean_per75)),
     xlim = range(c(sigma_df2$mean_per25, sigma_df2$mean_per75)))
arrows(x0 = sigma_df2$mean, y0 = sigma_df2_full$mean_per25,
       x1 = sigma_df2$mean, y1 = sigma_df2_full$mean_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = sigma_df2$mean_per25, y0 = sigma_df2_full$mean,
       x1 = sigma_df2$mean_per75, y1 = sigma_df2_full$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")

points(sigma_df2$mean, sigma_df2_full$mean,
       pch = 16, col = "#046C9A", cex = 1.5)

abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)
points(sigma_df2$mean, sigma_df2_full$mean, pch = 16, col = "#046C9A", cex = 1.5)
text(sigma_df2$mean_per75, sigma_df2_full$mean_per25, labels = sigma_df2$sigma, pos = c(3,3), cex = 0.75)

# bspp
plot(bspp_df2$fit_bspp, bspp_df2_full$fit_bspp,
     xlab = "restricted", ylab = "full", main = "bspp", type = "n", frame = FALSE,
     ylim = range(c(bspp_df2_full$fit_bspp_per25, bspp_df2_full$fit_bspp_per75)),
     xlim = range(c(bspp_df2$fit_bspp_per25, bspp_df2$fit_bspp_per75)))
arrows(x0 = bspp_df2$fit_bspp, y0 = bspp_df2_full$fit_bspp_per25,
       x1 = bspp_df2$fit_bspp, y1 = bspp_df2_full$fit_bspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = bspp_df2$fit_bspp_per25, y0 = bspp_df2_full$fit_bspp,
       x1 = bspp_df2$fit_bspp_per75, y1 = bspp_df2_full$fit_bspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(bspp_df2$fit_bspp, bspp_df2_full$fit_bspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

# aspp
plot(aspp_df2$fit_aspp, aspp_df2_full$fit_aspp,
     xlab = "restricted", ylab = "full", main = "aspp", type = "n", frame = FALSE,
     ylim = range(c(aspp_df2_full$fit_aspp_per25, aspp_df2_full$fit_aspp_per75)),
     xlim = range(c(aspp_df2$fit_aspp_per25, aspp_df2$fit_aspp_per75)))
arrows(x0 = aspp_df2$fit_aspp, y0 = aspp_df2_full$fit_aspp_per25,
       x1 = aspp_df2$fit_aspp, y1 = aspp_df2_full$fit_aspp_per75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = aspp_df2$fit_aspp_per25, y0 = aspp_df2_full$fit_aspp,
       x1 = aspp_df2$fit_aspp_per75, y1 = aspp_df2_full$fit_aspp,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(aspp_df2$fit_aspp, aspp_df2_full$fit_aspp,
       pch = 16, col = "#046C9A", cex = 1.5)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

###### Plot treeid ######
treeid_df2_full <- subset(treeid_df2_full, !treeid %in% setdiff(treeid_df2_full$treeid, treeid_df2$treeid))
plot(treeid_df2$fit_atreeid, treeid_df2_full$fit_atreeid,
     xlab = "restricted", ylab = "full", main = "atreeid", type = "n", frame = FALSE,
     ylim = range(c(treeid_df2_full$fit_atreeid_per25, treeid_df2_full$fit_atreeid_per75)),
     xlim = range(c(treeid_df2$fit_atreeid_per25, treeid_df2$fit_atreeid_per75)))
arrows(x0 = treeid_df2$fit_atreeid, y0 = treeid_df2_full$fit_atreeid_per25,
       x1 = treeid_df2$fit_atreeid, y1 = treeid_df2_full$fit_atreeid_per75,
       angle = 90, code = 3, length = 0, lwd = 0.8, col = "darkgray")
arrows(x0 = treeid_df2$fit_atreeid_per25, y0 = treeid_df2_full$fit_atreeid,
       x1 = treeid_df2$fit_atreeid_per75, y1 = treeid_df2_full$fit_atreeid,
       angle = 90, code = 3, length = 0, lwd = 0.8, col = "darkgray")
points(treeid_df2$fit_atreeid, treeid_df2_full$fit_atreeid,
       pch = 16, col = "#046C9A", cex = 0.8)
abline(0, 1, lty = 2, col = "#B40F20", lwd = 2)

dev.off()