# Wildchrokie model
# CRD 21 April 2026
# add a baseline on predictors to compare mostly, site

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
fitmodels <- FALSE

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
nrow(emp)
# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# order by tree id
treeid_spp_site <- unique(emp[, c("treeid_num", "spp_num", "site_num",
                                  "treeid", "spp", "site", "latbi")])

treeid_spp_site_ordered <- treeid_spp_site[order(treeid_spp_site$treeid_num), ]

# scale gdd to how many gdd are in 10 average spring days
temp<- subset(gddyr, doy <151 & doy > 120)
temp$mingddperiod <- ave(temp$GDD_5, temp$year, FUN = min)
temp$gdddiff <- temp$GDD_5 - temp$mingddperiod

temp <- temp[order(temp$year, temp$doy), ]

temp$bin10 <- ave(temp$doy, temp$year, FUN = function(x) ceiling((x - min(x) + 1) / 10))
gdd_10day <- aggregate(gdddiff ~ year + bin10, data = temp, max)
wcgddscale <- mean(gdd_10day$gdddiff)

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
  treeid_species = treeid_spp_site_ordered$spp_num,
  treeid_site = treeid_spp_site_ordered$site_num,
  Ntreeid_per_spp = as.integer(table(treeid_spp_site_ordered$spp_num)),
  gdd = (emp$pgsGDD5 - 1800)/wcgddscale,
  gddseq = gddseq,
  wcgddscale = wcgddscale,
  Ngddseq = length(gddseq)
)
dgdd

# Fit model GDD
gddmodel <- stan_model("stan/modelGrowthGDD.stan")
fitgdd <- sampling(gddmodel, data = dgdd,
                warmup = 1000, iter = 2000, chains=4)
saveRDS(fitgdd, "output/stanOutput/fitGrowthGDDbaseline")

# check warnings
diagnostics <- util$extract_hmc_diagnostics(fitgdd) 
util$check_all_hmc_diagnostics(diagnostics)

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
pdf(file = "figures/growthModelsMain/diagnostics/gddPvsPbaseline.pdf", width = 8, height = 10)

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

dev.off()



# 1:1 line plot for asite between the model with and without a baseline
# recover no baseline site
fitgddnobase <- readRDS("output/stanOutput/fitGrowthGDD")
df_fitgdd_nobase <- as.data.frame(fitgddnobase)

sigma_df2_nobase  <- extract_params(df_fitgdd_nobase, "sigma", "mean", "sigma")
bspp_df2_nobase   <- extract_params(df_fitgdd_nobase, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2_nobase <- extract_params(df_fitgdd_nobase, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_nobase <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2_nobase   <- extract_params(df_fitgdd_nobase, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2_nobase <- extract_params(df_fitgdd_nobase, "asite", "fit_a_site", 
                                  "site", "asite\\[(\\d+)\\]")

jpeg("figures/growthModelsMain/baselineVSnoneAsite.jpeg", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(1,1))

# asite
plot(site_df2$mean, site_df2_nobase$mean,
     xlab = "with baseline", ylab = "no baselin", main = "asite", type = "n", frame = FALSE,
     ylim = range(c(site_df2_nobase$p25, site_df2_nobase$p75)),
     xlim = range(c(site_df2$p25, site_df2$p75)))
arrows(x0 = site_df2$mean, y0 = site_df2_nobase$p25,
       x1 = site_df2$mean, y1 = site_df2_nobase$p75,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
arrows(x0 = site_df2$p25, y0 = site_df2_nobase$mean,
       x1 = site_df2$p75, y1 = site_df2_nobase$mean,
       angle = 90, code = 3, length = 0, lwd = 1.5, col = "darkgray")
points(site_df2$mean, site_df2_nobase$mean,
       pch = 16, col = "#0a6a3c", cex = 1.5)
abline(0, 1, lty = 2, col = "black", lwd = 2)
dev.off()

# just look at the intercepts
aspp_df2$a <- mean(df_fitgdd[,"a"])
aspp_df2_nobase$a <- mean(df_fitgdd_nobase[,"a"])
aspp_df2$a_aspp <- aspp_df2$mean + aspp_df2$a
aspp_df2_nobase$a_aspp <- aspp_df2_nobase$mean + aspp_df2_nobase$a

aspp_df2$bspp <- bspp_df2$mean[match(aspp_df2$spp, bspp_df2$spp)]
aspp_df2_nobase$bspp <- bspp_df2_nobase$mean[match(aspp_df2_nobase$spp, bspp_df2$spp)]

# spp names
aspp_df2$spp_name <- emp$latbi[match(aspp_df2$spp, emp$spp_num)]
aspp_df2_nobase$spp_name <- emp$latbi[match(aspp_df2_nobase$spp, emp$spp_num)]

jpeg("figures/growthModelsMain/baselineVSnoneSlopes.jpeg", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(2,1))
# no baseline
plot(x = emp$pgsGDD5, y = dgdd$y,
     xlim = c(0, max(emp$pgsGDD5)), ylim = c(-4, 3),
     xlab = "gdd with no baseline", ylab = "log(ring width)",
     main = "gdd model without a baseline")
sppvenum <- as.numeric(aspp_df2_nobase$spp)
for (i in sppvenum) { # i = 3
  d <- aspp_df2_nobase[aspp_df2_nobase$spp == i,]
  
  spp_name <- aspp_df2_nobase$spp_name[i]
  line_col <- wccolslatbi[spp_name]
  
  abline(a = d$a_aspp, b = d$bspp / wcgddscale, col = line_col)  
  abline(a = d$a, b = 0, lty = 2)
  # abline(v = 1800, lty = 1, lwd = 0.5)

}
             
# with a baseline
plot(x = emp$pgsGDD5 - 1800, y = dgdd$y,
     xlim = c(-1000, 1000), ylim = c(-4, 3),
     xlab = "gdd with a baseline", ylab = "log(ring width)",
     main = "gdd model with a baseline")
sppvenum <- as.numeric(aspp_df2$spp)

for (i in sppvenum) { # i = 3
  d <- aspp_df2[aspp_df2$spp == i,]
  
  spp_name <- aspp_df2$spp_name[i]
  line_col <- wccolslatbi[spp_name]
  
  abline(a = d$a_aspp, b = d$bspp / wcgddscale, col = line_col)  
  abline(a = d$a, b = 0, lty = 2)
  abline(0,1, lty = 1, lwd = 0.5)
  text(x = 350, y = -3, "baseline = 1800 gdd")
}
dev.off()
