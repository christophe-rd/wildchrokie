# Wildchrokie model
# CRD 19 March 2026

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
emp <- read.csv("output/empiricalDataMAIN.csv")
rw <- read.csv("output/wildchrokieRingWidth.csv")
gdd <- read.csv("/Users/christophe_rouleau-desrochers/github/coringtreespotters/analyses/output/gddByYear.csv")

commonNames <- c(
  "ALNINC" = "Alnus incana",           
  "BETALL" = "Betula alleghaniensis",
  "BETPAP" = "Betula papyrifera",      
  "BETPOP" = "Betula populifolia"     
)

emp$commonName <- commonNames[emp$latbi]
# calculate gdd by year between 1st april (91) and 1st august (213)  
gddsub <- subset(gdd, doy >90 & doy < 214)
agg <- aggregate(GDD_5 ~ year, gddsub, FUN = min)
gddsub$mingdd <- agg$GDD_5[match(gddsub$year, agg$year)]
gddsub$diff <- gddsub$GDD_5 - gddsub$mingdd
agg2 <- aggregate(diff ~ year, gddsub, FUN = max)

years <- unique(rw$year)
rw$yeardiff <- NA

for (i in years) {
  rw$yeardiff[rw$year == i] <- rw$year[rw$year == i] - 1
}

# add gdd of previous year indexed by yeardiff
rw$gdd <- agg2$diff[match(rw$yeardiff, agg2$year)]

rw$lengthMM <- rw$lengthCM * 10

# transform my groups to numeric values
rw$site_num <- match(rw$site, unique(rw$site))
rw$spp_num <- match(rw$spp, unique(rw$spp))
rw$treeid_num <- match(rw$treeid, unique(rw$treeid))

# transform data in vectors for GDD
rw <- rw[!is.na(rw$lengthMM) & !is.na(rw$gdd),]
rw <- subset(rw, site != "XX" & !spp %in% c("QUERNA", "QUEROB"))
y <- rw$lengthMM
N <- nrow(rw)
Nspp <- length(unique(rw$spp_num))
Nsite <- length(unique(rw$site_num))
site <- as.numeric(as.character(rw$site_num))
species <- as.numeric(as.character(rw$spp_num))
treeid <- as.numeric(rw$treeid_num)
Ntreeid <- length(unique(treeid))
gdd <- rw$gdd

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GDD posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fit <- readRDS("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/stanOutput/fitGrowthPreviousYear")

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot lines with quantiles ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
treeid_df2$treeid <- as.numeric(treeid_df2$treeid)
treeid_df2$treeid_name <- rw$treeid[match(treeid_df2$treeid, rw$treeid_num)]
bspp_df2$spp_name <- rw$latbi[match(bspp_df2$spp, rw$spp_num)]
site_df2$site_name <- rw$site[match(site_df2$site, rw$site_num)]
aspp_df2$spp_name <- rw$latbi[match(aspp_df2$spp, rw$spp_num)]

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
for (i in 1:length(unique(rw$treeid_num))) {
  subyvec[i] <- paste("atreeid", "[",i,"]", sep = "")  
}
subyvec

# get the spp and site identities for each tree id
treeid_spp_site <- unique(rw[, c("treeid_num", "spp_num", "site_num",
                                  "treeid", "spp", "site", "latbi")])

atreeidsub <- subset(df_fit, select = subyvec)
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
treeid_aspp <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
colnames(treeid_aspp) <- colnames(atreeidsub)

for (i in seq_len(ncol(treeid_aspp))) { # i = 1
  tree_id <- as.integer(colnames(treeid_aspp)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_aspp[, i] <- aspp_df[, spp_id]
}
treeid_aspp

# the site values for each tree id
treeid_asite <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
colnames(treeid_asite) <- colnames(atreeidsub)

for (i in seq_len(ncol(treeid_asite))) { # i = 1
  tree_id <- as.integer(colnames(treeid_asite)[i])
  site_id <- treeid_spp_site$site_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_asite[, i] <- site_df[, site_id]
}
treeid_asite

# recover a
treeid_a <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
colnames(treeid_a) <- colnames(atreeidsub)

for (i in seq_len(ncol(treeid_a))) { # i = 1
  treeid_a[, i] <- df_fit[, "a"]
}

# sum all 3 dfs together to get the full intercept for each treeid
fullintercept <-
  treeid_a + 
  atreeidsub +
  treeid_aspp +
  treeid_asite
fullintercept

# now get the slope for each treeid
treeid_bspp <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
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
x <- seq(min(rw$pgsGDD5), max(rw$pgsGDD5), length.out = 100)  
y_post_list <- list()  # store posterior predictions in a list where each tree id gets matrixad

# below I create a list where each row is the posterior estimate for each value of gdd (so the first row correspond to the model estimate for the first gdd value stored in x) and each column is the iteration (from 1 to 8000)
for (i in seq_along(treeidvecnum)) { # i = 1
  tree_col <- as.character(treeidvecnum[i]) 

  y_post <- sapply(1:nrow(df_fit), function(f) {
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
pdf(file = "figures/rwiricalData/growthModelSlopesperTreeid.pdf", width = 10, height = 8)

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
  rw_treeid <- rw[rw$treeid_num == tree_id_num, ]
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  # rwty plot first
  plot(rw$pgsGDD5, y, type = "n", 
       ylim = range(c(rw_treeid$lengthCM * 10, y_low, y_high), na.rm = TRUE),
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
    rw_treeid$pgsGDD5,
    rw_treeid$lengthCM * 10,
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
  mean_post_list[[tree_col]] <- sapply(1:nrow(df_fit), function(f) {
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
  sapply(1:nrow(df_fit), function(f) {
    rnorm(length(x), mean_mat[, f], sigma_df$sigma_y[f])
  })
})

x <- seq(min(rw$pgsGDD5), max(rw$pgsGDD5), length.out = 100)   
sppcols <- c(wes_palette("AsteroidCity1"))[1:4]

# jpeg output
jpeg(
  filename = "figures/rwiricalData/growthModelSlopesperSppFacet.jpeg",
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
  
  # subset rwirical data correctly
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



