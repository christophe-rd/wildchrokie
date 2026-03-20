# Wildchrokie model
# CRD 14 December 2025

# Goal: Plot model output because modelGrowthGDD is becoming too long and messy

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)
options(max.print = 150)
options(mc.cores = parallel::detectCores())
options(digits = 3)

# Load library 
library(ggplot2)
library(rstan)
library(future)
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
# my function to extract parameters
source('rcode/utilExtractParam.R')

# flags
makeplots <- FALSE
interceptmuplots <- TRUE
# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
climatesum <- read.csv("output/climateSummariesYear.csv")
gddyr <- read.csv("output/gddByYear.csv")
emp <- read.csv("output/empiricalDataMAIN.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

commonNames <- c(
  "Alnus incana"          = "Grey alder",
  "Betula alleghaniensis" = "Yellow birch",
  "Betula papyrifera"     = "Paper birch",
  "Betula populifolia"    = "Gray birch"
)

emp$commonName <- commonNames[emp$latbi]
emp$lengthMM <- emp$lengthCM*10
emp <- emp[!is.na(emp$pgsGDD5),]
# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors
y <- emp$lengthMM # ring width in mm
N <- nrow(emp)
gdd <- emp$pgsGDD5/200
Nspp <- length(unique(emp$spp_num))
Nsite <- length(unique(emp$site_num))
site <- as.numeric(as.character(emp$site_num))
species <- as.numeric(as.character(emp$spp_num))
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GDD posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgdd <- readRDS("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/stanOutput/fitGrowthGDD")

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GSL posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgsl <- readRDS("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/stanOutput/fitGrowthGSL")

df_fitgsl <- as.data.frame(fitgsl)

# full posterior
columns <- colnames(df_fitgsl)[!grepl("prior", colnames(df_fitgsl))]
sigma_df_gsl <- df_fitgsl[, columns[grepl("sigma", columns)]]
bspp_df_gsl <- df_fitgsl[, columns[grepl("bsp", columns)]]
treeid_df_gsl <- df_fitgsl[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df_gsl <- df_fitgsl[, columns[grepl("aspp", columns)]]
site_df_gsl <- df_fitgsl[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df_gsl) <- 1:ncol(bspp_df)
colnames(treeid_df_gsl) <- 1:ncol(treeid_df)
colnames(aspp_df_gsl) <- 1:ncol(aspp_df)
colnames(site_df_gsl) <- 1:ncol(site_df)

# posterior summaries
sigma_df2_gsl  <- extract_params(df_fitgsl, "sigma", "mean", "sigma")
bspp_df2_gsl   <- extract_params(df_fitgsl, "bsp", "fit_bspp", 
                                 "spp", "bsp\\[(\\d+)\\]")
treeid_df2_gsl <- extract_params(df_fitgsl, "atreeid", "fit_atreeid", 
                                 "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_gsl <- subset(treeid_df2, !grepl("z|sigma", treeid))
aspp_df2_gsl   <- extract_params(df_fitgsl, "aspp", "fit_aspp", 
                                 "spp", "aspp\\[(\\d+)\\]")
site_df2_gsl   <- extract_params(df_fitgsl, "asite", "fit_a_site", 
                                 "site", "asite\\[(\\d+)\\]")

treeid_df2_gsl$treeid <- as.numeric(treeid_df2_gsl$treeid)
treeid_df2_gsl$treeid_name <- emp$treeid[match(treeid_df2_gsl$treeid, emp$treeid_num)]
bspp_df2_gsl$spp_name <- emp$latbi[match(bspp_df2_gsl$spp, emp$spp_num)]
site_df2_gsl$site_name <- emp$site[match(site_df2_gsl$site, emp$site_num)]
aspp_df2_gsl$spp_name <- emp$latbi[match(aspp_df2_gsl$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# SOS posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitsos <- readRDS("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/stanOutput/fitGrowthSOS")

df_fitsos <- as.data.frame(fitsos)

# full posterior
columns <- colnames(df_fitsos)[!grepl("prior", colnames(df_fitsos))]
sigma_df_sos <- df_fitsos[, columns[grepl("sigma", columns)]]
bspp_df_sos <- df_fitsos[, columns[grepl("bsp", columns)]]
treeid_df_sos <- df_fitsos[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df_sos <- df_fitsos[, columns[grepl("aspp", columns)]]
site_df_sos <- df_fitsos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df_sos) <- 1:ncol(bspp_df_sos)
colnames(treeid_df_sos) <- 1:ncol(treeid_df_sos)
colnames(aspp_df_sos) <- 1:ncol(aspp_df_sos)
colnames(site_df_sos) <- 1:ncol(site_df_sos)

# posterior summaries
sigma_df2_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", 
                                 "spp", "bsp\\[(\\d+)\\]")
treeid_df2_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", 
                                 "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_sos <- subset(treeid_df2_sos, !grepl("z|sigma", treeid))
aspp_df2_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", 
                                 "spp", "aspp\\[(\\d+)\\]")
site_df2_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", 
                                 "site", "asite\\[(\\d+)\\]")

treeid_df2_sos$treeid <- as.numeric(treeid_df2_sos$treeid)
treeid_df2_sos$treeid_name <- emp$treeid[match(treeid_df2_sos$treeid, emp$treeid_num)]
bspp_df2_sos$spp_name <- emp$latbi[match(bspp_df2_sos$spp, emp$spp_num)]
site_df2_sos$site_name <- emp$site[match(site_df2_sos$site, emp$site_num)]
aspp_df2_sos$spp_name <- emp$latbi[match(aspp_df2_sos$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# EOS posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fiteos <- readRDS("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/stanOutput/fitGrowthEOS")

df_fiteos <- as.data.frame(fiteos)

# full posterior
columns <- colnames(df_fiteos)[!grepl("prior", colnames(df_fiteos))]
sigma_df_eos <- df_fiteos[, columns[grepl("sigma", columns)]]
bspp_df_eos <- df_fiteos[, columns[grepl("bsp", columns)]]
treeid_df_eos <- df_fiteos[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df_eos <- df_fiteos[, columns[grepl("aspp", columns)]]
site_df_eos <- df_fiteos[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df_eos) <- 1:ncol(bspp_df_eos)
colnames(treeid_df_eos) <- 1:ncol(treeid_df_eos)
colnames(aspp_df_eos) <- 1:ncol(aspp_df_eos)
colnames(site_df_eos) <- 1:ncol(site_df_eos)

# posterior summaries
sigma_df2_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", 
                                 "spp", "bsp\\[(\\d+)\\]")
treeid_df2_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", 
                                 "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_eos <- subset(treeid_df2_eos, !grepl("z|sigma", treeid))
aspp_df2_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", 
                                 "spp", "aspp\\[(\\d+)\\]")
site_df2_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", 
                                 "site", "asite\\[(\\d+)\\]")

treeid_df2_eos$treeid <- as.numeric(treeid_df2_eos$treeid)
treeid_df2_eos$treeid_name <- emp$treeid[match(treeid_df2_eos$treeid, emp$treeid_num)]
bspp_df2_eos$spp_name <- emp$latbi[match(bspp_df2_eos$spp, emp$spp_num)]
site_df2_eos$site_name <- emp$site[match(site_df2_eos$site, emp$site_num)]
aspp_df2_eos$spp_name <- emp$latbi[match(aspp_df2_eos$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot lines with quantiles ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
treeid_df2$treeid <- as.numeric(treeid_df2$treeid)
treeid_df2$treeid_name <- emp$treeid[match(treeid_df2$treeid, emp$treeid_num)]
bspp_df2$spp_name <- emp$latbi[match(bspp_df2$spp, emp$spp_num)]
site_df2$site_name <- emp$site[match(site_df2$site, emp$site_num)]
aspp_df2$spp_name <- emp$latbi[match(aspp_df2$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Define objects used throught the models ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
sitefull <- c(
  "GR" = "Dartmouth College (NH)",
  "HF" = "Harvard Forest (MA)",
  "SH" = "St-Hyppolyte (Qc)",
  "WM" = "White Mountains (NH)"
)

sppcols <- c(wes_palette("AsteroidCity1"))[1:4]

subyvec <- vector()
for (i in 1:length(unique(emp$treeid_num))) {
  subyvec[i] <- paste("atreeid", "[",i,"]", sep = "")  
}
subyvec

# get the spp and site identities for each tree id
treeid_spp_site <- unique(emp[, c("treeid_num", "spp_num", "site_num",
                                  "treeid", "spp", "site", "latbi")])

atreeidsub <- subset(df_fitgdd, select = subyvec)
colnames(atreeidsub) <- 1:length(subyvec)

# get a vector for each treeid for each species
spp1vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 1]
spp2vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 2]
spp3vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 3]
spp4vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 4]

spp_list <- list(
  "1" = spp1vec,
  "2" = spp2vec,
  "3" = spp3vec,
  "4" = spp4vec
)

sppvecnum <- 1:4
sppvecname <- unique(treeid_spp_site$latbi)

n_spp <- nrow(bspp_df2)
n_site <- nrow(site_df2)
y_pos <- 1:n_spp 

sitecolors <- c(wes_palette("Darjeeling1"))[1:4]

if(makeplots) {
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GDD: Prep posterior reconstruction #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# start by filling a df with treeid intercepts only
# the spp values for each tree id
treeid_aspp <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fitgdd)))
colnames(treeid_aspp) <- colnames(atreeidsub)

for (i in seq_len(ncol(treeid_aspp))) { # i = 1
  tree_id <- as.integer(colnames(treeid_aspp)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_aspp[, i] <- aspp_df[, spp_id]
}
treeid_aspp

# the site values for each tree id
treeid_asite <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fitgdd)))
colnames(treeid_asite) <- colnames(atreeidsub)

for (i in seq_len(ncol(treeid_asite))) { # i = 1
  tree_id <- as.integer(colnames(treeid_asite)[i])
  site_id <- treeid_spp_site$site_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_asite[, i] <- site_df[, site_id]
}
treeid_asite

# recover a
treeid_a <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fitgdd)))
colnames(treeid_a) <- colnames(atreeidsub)

for (i in seq_len(ncol(treeid_a))) { # i = 1
  treeid_a[, i] <- df_fitgdd[, "a"]
}

# sum all 3 dfs together to get the full intercept for each treeid
fullintercept <-
  treeid_a + 
  atreeidsub +
  treeid_aspp +
  treeid_asite
fullintercept

# now get the slope for each treeid
treeid_bspp <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fitgdd)))
colnames(treeid_bspp) <- colnames(atreeidsub)

# back convert the slopes to their original scales
bspp_df4 <- bspp_df
for (i in 1:ncol(bspp_df4)){
  bspp_df4[[i]] <- bspp_df4[[i]] / 200
}

for (i in seq_len(ncol(treeid_bspp))) { # i = 30
  tree_id <- as.integer(colnames(treeid_bspp)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_bspp[, i] <- bspp_df4[, spp_id]
}
treeid_bspp

treeidvecnum <- 1:ncol(fullintercept)
treeidvecname <- treeid_spp_site$treeid
x <- seq(min(emp$pgsGDD5), max(emp$pgsGDD5), length.out = 100)  
y_post_list <- list()  # store posterior predictions in a list where each tree id gets matrixad

# below I create a list where each row is the posterior estimate for each value of gdd (so the first row correspond to the model estimate for the first gdd value stored in x) and each column is the iteration (from 1 to 8000)
for (i in seq_along(treeidvecnum)) { # i = 1
  tree_col <- as.character(treeidvecnum[i]) 

  y_post <- sapply(1:nrow(df_fitgdd), function(f) {
    rnorm(length(x), 
          fullintercept[f, tree_col] + treeid_bspp[f, tree_col] * x,
          sigma_df$sigma_y[f])
  })
  y_post_list[[tree_col]] <- y_post
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GDD: per treeid, facet #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# PDF output
pdf(file = "figures/empiricalData/growthModelSlopesperTreeid.pdf", width = 10, height = 8)

# Layout: 2 rows × 2 columns per page
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Loop over trees again to plot each tree individually
for (i in seq_along(treeidvecnum)) { # i = 1
  tree_col <- as.character(treeidvecnum[i])
  tree_col_name <- as.character(treeidvecname[i])
  y_post <- y_post_list[[tree_col]]
  
  # color line by spp
  tree_id_num <- as.integer(tree_col)
  
  # index the dots per treeid
  emp_treeid <- emp[emp$treeid_num == tree_id_num, ]
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  # empty plot first
  plot(emp$pgsGDD5, y, type = "n", 
       ylim = range(c(emp_treeid$lengthCM * 10, y_low, y_high), na.rm = TRUE),
       xlab = "Primary growing season GDD", ylab = "Ring width (mm)",
       main = tree_col_name) # set the name for each plot
  
  spp_id <- treeid_spp_site$spp_num[
    match(tree_id_num, treeid_spp_site$treeid_num)
  ]
  line_col <- sppcols[spp_id]
  
  # shaded interval
  polygon(c(x, rev(x)), 
          c(y_low, # lower interval
            rev(y_high)), # high interval
          col = adjustcolor(line_col, alpha.f = 0.3), 
          border = NA)
  
  # mean line
  lines(x, y_mean,
        col = line_col,
        lwd = 2)
  
  points(
    emp_treeid$pgsGDD5,
    emp_treeid$lengthCM * 10,
    pch = 16,
    cex = 2,
    col = line_col)
}
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GDD: per Spp, facet #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
mean_post_list <- list()
for (i in seq_along(treeidvecnum)) {
  tree_col <- as.character(treeidvecnum[i])
  mean_post_list[[tree_col]] <- sapply(1:nrow(df_fitgdd), function(f) {
    fullintercept[f, tree_col] + treeid_bspp[f, tree_col] * x
    # no sigma_y yet
  })
}

# average the mean predictions across trees within species
spp_mean_list <- lapply(spp_list, function(tree_vec) {
  Reduce("+", mean_post_list[as.character(tree_vec)]) / length(tree_vec)
})

# re-simulate sigma on the averaged mean
spp_post_list <- lapply(spp_mean_list, function(mean_mat) {
  sapply(1:nrow(df_fitgdd), function(f) {
    rnorm(length(x), mean_mat[, f], sigma_df$sigma_y[f])
  })
})

x <- seq(min(emp$pgsGDD5), max(emp$pgsGDD5), length.out = 100)   
sppcols <- c(wes_palette("AsteroidCity1"))[1:4]

# jpeg output
jpeg(
  filename = "figures/empiricalData/growthModelSlopesperSppFacet.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)
# Layout: 2 rows × 2 columns per page
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  
  # define spp num
  spp_num <- as.integer(spp_column)
  
  # subset empirical data correctly
  emp_spp <- emp[emp$spp_num == spp_num, ]
  
  spp_column <- as.character(sppvecnum[i]) 
  y_post <- spp_post_list[[spp_column]]
  
  # summaries
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  # species-specific ylim
  ylim_spp <- range(c(emp_spp$lengthCM * 10, y_low, y_high), na.rm = TRUE)
  
  plot(emp_spp$pgsGDD5, emp_spp$lengthCM * 10,
       type = "n",
       ylim = ylim_spp,
       xlab = "Primary growing season GDD",
       ylab = "Ring width (mm)",
       main = spp_column_name,
       frame = FALSE)
  
  # color
  line_col <- sppcols[spp_num]
  
  polygon(
    c(x, rev(x)),
    c(y_low, rev(y_high)),
    col = adjustcolor(line_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(x, y_mean, col = line_col, lwd = 2)
  
  points(
    emp_spp$pgsGDD5,
    emp_spp$lengthCM * 10,
    pch = 16,
    cex = 1,
    col = line_col
  )
}

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### GDD: per spp, non facetted #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# PDF output
jpeg(
  filename = "figures/empiricalData/growthModelSlopesperSpp.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)

# Layout: 2 rows × 2 columns per page
par(mar = c(4, 4, 2, 1))

# below I create a list where each row is the posterior estimate for each value of gdd (so the first row correspond to the model estimate for the first gdd value stored in x) and each column is the iteration (from 1 to 8000)

plot(emp$pgsGDD5, y, type = "n", 
     ylim = range(min(emp$lengthCM*10), max(emp$lengthCM*10)), 
     xlab = "Primary growing season GDD", ylab = "Ring width (mm)",
     main = "species growth responses")

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  y_post <- spp_post_list[[spp_column]]
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  # shaded interval
  polygon(c(x, rev(x)), 
          c(y_low, # lower interval
            rev(y_high)), # high interval
          col = adjustcolor(line_col, alpha.f = 0.2), 
          border = NA)
  
  # mean line
  lines(x, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]

  points(
    emp_spp$pgsGDD5,
    emp_spp$lengthCM * 10,
    pch = 16,
    cex = 1,
    col = line_col)
  
  legend(
    "topleft",
    legend = sppvecname,
    col = sppcols[as.integer(sppvecnum)],
    lwd = 2,
    pch = 16,
    bty = "n",
    cex = 1.5
  )
  
}
dev.off()


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GSL: prep posterior reconstruction ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
subyvec <- vector()
for (i in 1:length(unique(emp$treeid_num))) {
  subyvec[i] <- paste("atreeid", "[",i,"]", sep = "")  
}
subyvec

atreeidsub_gsl <- subset(df_fitgsl, select = subyvec)

colnames(atreeidsub_gsl) <- 1:length(subyvec)

# the spp values for each tree id
treeid_aspp_gsl <- data.frame(matrix(ncol = ncol(atreeidsub_gsl), nrow = nrow(df_fitgsl)))
colnames(treeid_aspp_gsl) <- colnames(atreeidsub_gsl)

for (i in seq_len(ncol(treeid_aspp_gsl))) { # i = 1
  tree_id <- as.integer(colnames(treeid_aspp_gsl)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_aspp_gsl[, i] <- aspp_df_gsl[, spp_id]
}
treeid_aspp_gsl

# the site values for each tree id
treeid_asite_gsl <- data.frame(matrix(ncol = ncol(atreeidsub_gsl), nrow = nrow(df_fitgsl)))
colnames(treeid_asite_gsl) <- colnames(atreeidsub_gsl)

for (i in seq_len(ncol(treeid_asite_gsl))) { # i = 1
  tree_id <- as.integer(colnames(treeid_asite_gsl)[i])
  site_id <- treeid_spp_site$site_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_asite_gsl[, i] <- site_df_gsl[, site_id]
}
treeid_asite_gsl

# recover a
treeid_a_gsl <- data.frame(matrix(ncol = ncol(atreeidsub_gsl), nrow = nrow(df_fitgsl)))
colnames(treeid_a_gsl) <- colnames(atreeidsub_gsl)

for (i in seq_len(ncol(treeid_a_gsl))) { # i = 1
  treeid_a_gsl[, i] <- df_fitgsl[, "a"]
}

# sum all 3 dfs together to get the full intercept for each treeid
fullintercept_gsl <-
  treeid_a_gsl + 
  atreeidsub_gsl +
  treeid_aspp_gsl +
  treeid_asite_gsl
fullintercept_gsl

# now get the slope for each treeid
treeid_bspp_gsl <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fitgsl)))
colnames(treeid_bspp_gsl) <- colnames(atreeidsub)

# back convert the slopes to their original scales
bspp_df4_gsl <- bspp_df_gsl
for (i in 1:ncol(bspp_df4_gsl)){
  bspp_df4_gsl[[i]] <- bspp_df4_gsl[[i]] / 10
}

for (i in seq_len(ncol(treeid_bspp_gsl))) { # i = 30
  tree_id <- as.integer(colnames(treeid_bspp_gsl)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_bspp_gsl[, i] <- bspp_df4_gsl[, spp_id]
}
treeid_bspp_gsl

treeidvecnum <- 1:ncol(fullintercept_gsl)
treeidvecname <- treeid_spp_site$treeid
x <- seq(min(emp$pgsGSL), max(emp$pgsGSL), length.out = 100)  

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GSL: per Spp, facet #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
mean_post_list_gsl <- list()
for (i in seq_along(treeidvecnum)) {
  tree_col <- as.character(treeidvecnum[i])
  mean_post_list_gsl[[tree_col]] <- sapply(1:nrow(df_fitgsl), function(f) {
    fullintercept_gsl[f, tree_col] + treeid_bspp_gsl[f, tree_col] * x
    # no sigma_y yet
  })
}

# average the mean predictions across trees within species
spp_mean_list_gsl <- lapply(spp_list, function(tree_vec) {
  Reduce("+", mean_post_list_gsl[as.character(tree_vec)]) / length(tree_vec)
})

# re-simulate sigma on the averaged mean
spp_post_list_gsl <- lapply(spp_mean_list_gsl, function(mean_mat) {
  sapply(1:nrow(df_fitgsl), function(f) {
    rnorm(length(x), mean_mat[, f], sigma_df_gsl$sigma_y[f])
  })
})

# jpeg output
jpeg(
  filename = "figures/empiricalData/growthModelSlopesperSppFacetGSL.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)
# Layout: 2 rows × 2 columns per page
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  
  # define spp num
  spp_num <- as.integer(spp_column)
  
  # spp color
  line_col <- sppcols[spp_num]
  
  # subset empirical data correctly
  emp_spp <- emp[emp$spp_num == spp_num, ]
  
  spp_column <- as.character(sppvecnum[i]) 
  y_post_gsl <- spp_post_list_gsl[[spp_column]]
  
  # summaries
  y_mean_gsl <- apply(y_post_gsl, 1, mean)
  y_low_gsl  <- apply(y_post_gsl, 1, quantile, 0.25)
  y_high_gsl <- apply(y_post_gsl, 1, quantile, 0.75)
  
  # species-specific ylim
  # ylim_spp <- range(c(emp_spp$lengthCM * 10, y_low, y_high), na.rm = TRUE)
  
  plot(emp_spp$pgsGSL, emp_spp$lengthCM * 10,
       type = "n",
       ylim = c(0,14),
       xlab = "Primary growing season GSL",
       ylab = "Ring width (mm)",
       main = spp_column_name,
       frame = FALSE)
  
  polygon(
    c(x, rev(x)),
    c(y_low_gsl, rev(y_high_gsl)),
    col = adjustcolor(line_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(x, y_mean_gsl, col = line_col, lwd = 2)
  
  points(
    emp_spp$pgsGSL,
    emp_spp$lengthCM * 10,
    pch = 16,
    cex = 1,
    col = line_col
  )
}

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Mu plots #####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
###### asp ######
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
jpeg("figures/empiricalData/aspp_mean_plot.jpeg", width = 8, height = 6, units = "in", res = 300)

par(mar = c(5, 10, 2, 2)) 

plot(aspp_df2$fit_aspp, y_pos,
     xlim = range(c(aspp_df2$fit_aspp_per5, aspp_df2$fit_aspp_per95)),
     ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)",
     ylab = "",
     yaxt = "n",
     pch = 16,
     cex = 2,
     col = sppcols,
     frame.plot = FALSE)

# color labels
for (i in seq_along(y_pos)) {
  axis(2, at = y_pos[i],
       labels = aspp_df2$spp_name[i],
       las = 2,
       col.axis = sppcols[i],
       tick = FALSE,
       cex.axis = 1)
}

# error bars and dashed line
segments(aspp_df2$fit_aspp_per5,  y_pos,
         aspp_df2$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2$fit_aspp_per25, y_pos,
         aspp_df2$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)

abline(v = 0, lty = 2, col = "black")

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
###### asite ######
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
emp$sitefull <- sitefull[emp$site]
site_df2$sitefull <- sitefull[site_df2$site_name]

jpeg("figures/empiricalData/asite_mean_plot.jpeg", width = 8, height = 6, units = "in", res = 300)

par(mar = c(5, 10, 2, 2)) 

plot(site_df2$fit_a_site, y_pos,
     xlim = range(c(site_df2$fit_a_site_per5, site_df2$fit_a_site_per95)),
     ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)",
     ylab = "",
     yaxt = "n",
     pch = 16,
     cex = 2,
     col = sitecolors,
     frame.plot = FALSE)

# color labels
for (i in seq_along(y_pos)) {
  axis(2, at = y_pos[i],
       labels = site_df2$sitefull[i],
       las = 2,
       col.axis = sitecolors[i],
       tick = FALSE,
       cex.axis = 1)
}

# error bars and dashed line
segments(site_df2$fit_a_site_per5,  y_pos,
         site_df2$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2$fit_a_site_per25, y_pos,
         site_df2$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)

abline(v = 0, lty = 2, col = "black")

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
###### bsp ######
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

jpeg("figures/empiricalData/bspp_mean_plot.jpeg", width = 8, height = 6, units = "in", res = 300)

par(mar = c(5, 10, 2, 2)) 

plot(bspp_df2$fit_bspp, y_pos,
     xlim = range(c(bspp_df2$fit_bspp_per5, bspp_df2$fit_bspp_per95)),
     ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width (mm) change/200 GDD",
     ylab = "",
     yaxt = "n",
     pch = 16,
     cex = 2,
     col = sppcols,
     frame.plot = FALSE)

# color labels
for (i in seq_along(y_pos)) {
  axis(2, at = y_pos[i],
       labels = bspp_df2$spp_name[i],
       las = 2,
       col.axis = sppcols[i],
       tick = FALSE,
       cex.axis = 1)
}

# error bars and dashed line
segments(bspp_df2$fit_bspp_per5,  y_pos,
         bspp_df2$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2$fit_bspp_per25, y_pos,
         bspp_df2$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

abline(v = 0, lty = 2, col = "black")

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### full treeid mu plots #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
  
# Mean plots with atreeid ####

# now do the same, but for species
treeid_df2$spp <- emp$spp[match(treeid_df2$treeid, emp$treeid_num)]

# same for site
treeid_df2$site <- emp$site[match(treeid_df2$treeid, emp$treeid_num)]

# quick check that I didn't mess anything up
un <- emp[!duplicated(emp$treeid),]
table(un$spp)
table(treeid_df2$spp)

sub <- subset(emp, select = c("treeid_num", "spp_num", "site_num"))
sub <- sub[!duplicated(sub$treeid_num),]

# recalculate the full intercept without the grand mean
fullintercept2 <-
  atreeidsub +
  treeid_aspp +
  treeid_asite
fullintercept2

# get posterior means and quantiles

# empty treeid dataframe
treeid_df4 <- data.frame(
  treeid = character(ncol(fullintercept2)),
  fit_atreeid = numeric(ncol(fullintercept2)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(fullintercept2)) { # i = 1
  treeid_df4$treeid[i] <- colnames(fullintercept2)[i]         
  treeid_df4$fit_atreeid[i] <- round(mean(fullintercept2[[i]]),3)  
  treeid_df4$fit_atreeid_per5[i] <- round(quantile(fullintercept2[[i]], probs = 0.05), 3)
  treeid_df4$fit_atreeid_per25[i] <- round(quantile(fullintercept2[[i]], probs = 0.25), 3)
  treeid_df4$fit_atreeid_per75[i] <- round(quantile(fullintercept2[[i]], probs = 0.75), 3)
  treeid_df4$fit_atreeid_per95[i] <- round(quantile(fullintercept2[[i]], probs = 0.95), 3)
}
treeid_df4

# get the og treeid names, spp and site back:
treeid_df4$treeid <- as.numeric(treeid_df4$treeid)
treeid_df4$treeid_name <- emp$treeid[match(treeid_df4$treeid,
                                                    emp$treeid_num)]
treeid_df4$spp_name <- emp$latbi[match(treeid_df4$treeid,
                                              emp$treeid_num)]
treeid_df4$spp_num <- emp$spp_num[match(treeid_df4$treeid,
                                       emp$treeid_num)]
treeid_df4$site_name <- emp$site[match(treeid_df4$treeid,
                                                emp$treeid_num)]
treeid_df4$site_num <- emp$site_num[match(treeid_df4$treeid,
                                       emp$treeid_num)]

# Prep for the figure

# define a gap between species clusters
gap <- 3

# y positions
# treeid_df4$y_pos <- NA
current_y <- 1

species_order <- c(
  "Alnus incana", 
  "Betula alleghaniensis", 
  "Betula papyrifera", 
  "Betula populifolia")

site_order <- c(
  "HF",
  "WM",
  "GR", 
  "SH")

# col
my_colors <- c(
  "Alnus incana" = wes_palette("AsteroidCity1")[1],
  "Betula alleghaniensis" = wes_palette("AsteroidCity1")[2],
  "Betula papyrifera" = wes_palette("AsteroidCity1")[3],
  "Betula populifolia" = wes_palette("AsteroidCity1")[4]
)
# shapes for sites
my_shapes <- c(
  HF = 19,
  WM = 18,
  GR = 15,
  SH = 17
)


# open device
pdf(
  file = "figures/empiricalData/meanPlotGrowthGDD_treeidBYspp.pdf",
  width = 8,  
  height = 8
)
par(mar = c(
  4, # lower margin 
  6, 
  4, # upper margin  
  6)) # right margin

treeid_df4$spp_name  <- factor(treeid_df4$spp_name, levels = species_order)
treeid_df4$site_num <- factor(treeid_df4$site_num, levels = site_order)

treeid_df4 <- treeid_df4[
  order(treeid_df4$spp_name, treeid_df4$site_num, treeid_df4$treeid),
]

treeid_df4$y_pos <- seq_len(nrow(treeid_df4))

total_rows <- nrow(treeid_df4) + (length(species_order) - 1) * gap
current_y <- total_rows 

for(sp in species_order){ # sp = "Alnus incana"
  idx <- which(treeid_df4$spp_name == sp)
  n <- length(idx)
  
  # assign sequential positions for this species
  treeid_df4$y_pos[idx] <- current_y:(current_y - n + 1)
  
  # move cursor down with a gap before next species cluster
  current_y <- current_y - n - gap
}

# Set up empty plot
plot(
  NA, NA,
  xlim = range(c(treeid_df4$fit_atreeid_per5-2,
                 treeid_df4$fit_atreeid_per95+2)),
  ylim = c(0.5, max(treeid_df4$y_pos) + 0.5),
  xlab = "treeid intercept values",
  ylab = "",
  yaxt = "n",
  bty = "l"
)


# --- Add horizontal error bars (5–95%) ---
segments(
  x0 = treeid_df4$fit_atreeid_per5,
  x1 = treeid_df4$fit_atreeid_per95,
  y0 = treeid_df4$y_pos,
  col = adjustcolor(my_colors[treeid_df4$spp_name], alpha.f = 0.7),
  lwd = 1
)

# --- Add thicker horizontal error bars (25–75%) ---
segments(
  x0 = treeid_df4$fit_atreeid_per25,
  x1 = treeid_df4$fit_atreeid_per75,
  y0 = treeid_df4$y_pos,
  col = adjustcolor(my_colors[treeid_df4$spp_name], alpha.f = 0.7),
  lwd = 1.5
)

# --- Add the points ---
points(
  treeid_df4$fit_atreeid,
  treeid_df4$y_pos,
  cex = 0.8,
  pch = my_shapes[treeid_df4$site_name],
  col = adjustcolor(my_colors[treeid_df4$spp_name], alpha.f = 0.7)
)

spp_y_top <- tapply(treeid_df4$y_pos, treeid_df4$spp_name, max)
aspp_df2$y_pos <- spp_y_top[aspp_df2$spp_name] + 1

segments(
  x0 = aspp_df2$fit_aspp_per5,
  x1 = aspp_df2$fit_aspp_per95,
  y0 = aspp_df2$y_pos,
  col = adjustcolor(my_colors[aspp_df2$spp_name], alpha.f = 0.9),
  lwd = 2
)

segments(
  x0 = aspp_df2$fit_aspp_per25,
  x1 = aspp_df2$fit_aspp_per75,
  y0 = aspp_df2$y_pos,
  col = my_colors[aspp_df2$spp_name],
  lwd = 3
)
points(
  aspp_df2$fit_aspp,
  aspp_df2$y_pos,
  pch = 16,
  bg  = my_colors[aspp_df2$spp_name],
  col = my_colors[aspp_df2$spp_name],
  cex = 1.5
)

# --- Add vertical line at 0 ---
abline(v = 0, lty = 2)

# --- Add custom y-axis labels (reverse order if needed) ---
axis(
  side = 2,
  at = treeid_df4$y_pos,
  labels = treeid_df4$treeid_name,
  cex.axis = 0.5,
  las = 1
)
# spp_name mean
spp_y <- tapply(treeid_df4$y_pos, treeid_df4$spp_name, mean)
site_y <- tapply(treeid_df4$y_pos, treeid_df4$site_num, max)

## order species by mean y descending (top of plot first)
species_legend_order <- names(sort(spp_y, decreasing = TRUE))
site_legend_order <- names(sort(site_y, decreasing = FALSE))

## species legend (colors matched by name)
legend(
  x = max(treeid_df4$fit_atreeid_per95) - 5,
  y = max(treeid_df4$y_pos) + 1,
  legend = species_legend_order,
  col = my_colors[species_legend_order],    # index so colors match
  pch = 16,
  pt.cex = 1.2,
  title = "Species",
  bty = "n"
)

site_legend_order <- c("SH", "GR", "WM", "HF")
# site_num legen
legend(
  x = max(treeid_df4$fit_atreeid_per95) - 2,
  y = max(treeid_df4$y_pos) - 15,
  legend = site_legend_order,
  pch = my_shapes[site_legend_order],
  pt.cex = 1.2,
  title = "Sites",
  bty = "n"
)
dev.off()


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GDD mu plots Together ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
jpeg(file = "figures/empiricalData/muGDD.jpeg", 
     width = 3000, height = 1000, res = 300)
par(mfrow = c(1, 3))
par(mar = c(5, 10, 2, 2)) 

##### asp #####
plot(aspp_df2$fit_aspp, y_pos,
     xlim = range(c(aspp_df2$fit_aspp_per5, aspp_df2$fit_aspp_per95)),
     ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)",
     ylab = "",
     yaxt = "n",
     pch = 16,
     cex = 2,
     col = sppcols,
     frame.plot = FALSE)

# color labels
for (i in seq_along(y_pos)) {
  axis(2, at = y_pos[i],
       labels = aspp_df2$spp_name[i],
       las = 2,
       col.axis = sppcols[i],
       tick = FALSE,
       cex.axis = 1)
}

# error bars and dashed line
segments(aspp_df2$fit_aspp_per5,  y_pos,
         aspp_df2$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2$fit_aspp_per25, y_pos,
         aspp_df2$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)

abline(v = 0, lty = 2, col = "black")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### asite #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

sitefull <- c(
  "GR" = "Dartmouth College (NH)",
  "HF" = "Harvard Forest (MA)",
  "SH" = "St-Hyppolyte (Qc)",
  "WM" = "White Mountains (NH)"
)

emp$sitefull <- sitefull[emp$site]
site_df2$sitefull <- sitefull[site_df2$site_name]

par(mar = c(5, 10, 2, 2)) 

plot(site_df2$fit_a_site, y_pos,
     xlim = range(c(site_df2$fit_a_site_per5, site_df2$fit_a_site_per95)),
     ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)",
     ylab = "",
     yaxt = "n",
     pch = 16,
     cex = 2,
     col = sitecolors,
     frame.plot = FALSE)

# color labels
for (i in seq_along(y_pos)) {
  axis(2, at = y_pos[i],
       labels = site_df2$sitefull[i],
       las = 2,
       col.axis = sitecolors[i],
       tick = FALSE,
       cex.axis = 1)
}

# error bars and dashed line
segments(site_df2$fit_a_site_per5,  y_pos,
         site_df2$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2$fit_a_site_per25, y_pos,
         site_df2$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)

abline(v = 0, lty = 2, col = "black")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### bsp #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = c(5, 10, 2, 2)) 

plot(bspp_df2$fit_bspp, y_pos,
     xlim = range(c(bspp_df2$fit_bspp_per5, bspp_df2$fit_bspp_per95)),
     ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width (mm) change/200 GDD",
     ylab = "",
     yaxt = "n",
     pch = 16,
     cex = 2,
     col = sppcols,
     frame.plot = FALSE)

# color labels
for (i in seq_along(y_pos)) {
  axis(2, at = y_pos[i],
       labels = bspp_df2$spp_name[i],
       las = 2,
       col.axis = sppcols[i],
       tick = FALSE,
       cex.axis = 1)
}

# error bars and dashed line
segments(bspp_df2$fit_bspp_per5,  y_pos,
         bspp_df2$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2$fit_bspp_per25, y_pos,
         bspp_df2$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

abline(v = 0, lty = 2, col = "black")
dev.off()

}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Phenology carry-over ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
plot(budset ~ leafout, emp)
lm <- lmer(budset ~ leafout + (1 | year), data = emp)
sum <- summary(lm)
new_data <- data.frame(leafout = x_seq)
pred <- predict(lm, newdata = new_data, re.form = NA)  
lines(x_seq, pred, col = "black", lwd = 2)
