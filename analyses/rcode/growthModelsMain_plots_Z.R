# Wildchrokie model
# CRD 2 April 2026

# plotting z-scored data

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

# source model code
source("rcode/growthModelsMain.R")

# flags
makeplots <- TRUE

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
climatesum <- read.csv("output/climateSummariesYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

commonNames <- c(
  "Alnus incana"          = "Grey alder",
  "Betula alleghaniensis" = "Yellow birch",
  "Betula papyrifera"     = "Paper birch",
  "Betula populifolia"    = "Gray birch"
)

emp$commonName <- commonNames[emp$latbi]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GDD posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgdd <- readRDS("output/stanOutput/fitGrowthGDDZscored")

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
sigma_df2_z  <- extract_params(df_fitgdd, "sigma", "mean", "sigma")
bspp_df2_z   <- extract_params(df_fitgdd, "bsp", "fit_bspp", 
                             "spp", "bsp\\[(\\d+)\\]")
treeid_df2_z <- extract_params(df_fitgdd, "atreeid", "fit_atreeid", 
                             "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_z <- subset(treeid_df2_z, !grepl("z|sigma", treeid))
aspp_df2_z   <- extract_params(df_fitgdd, "aspp", "fit_aspp", 
                             "spp", "aspp\\[(\\d+)\\]")
site_df2_z   <- extract_params(df_fitgdd, "asite", "fit_a_site", 
                             "site", "asite\\[(\\d+)\\]")


treeid_df2_z$treeid_name <- emp$treeid[match(treeid_df2_z$treeid, emp$treeid_num)]
bspp_df2_z$spp_name <- emp$latbi[match(bspp_df2_z$spp, emp$spp_num)]
site_df2_z$site_name <- emp$site[match(site_df2_z$site, emp$site_num)]
aspp_df2_z$spp_name <- emp$latbi[match(aspp_df2_z$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GSL posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgsl <- readRDS("output/stanOutput/fitGrowthGSLZscored")

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
sigma_df2_z_gsl  <- extract_params(df_fitgsl, "sigma", "mean", "sigma")
bspp_df2_z_gsl   <- extract_params(df_fitgsl, "bsp", "fit_bspp", 
                                 "spp", "bsp\\[(\\d+)\\]")
treeid_df2_z_gsl <- extract_params(df_fitgsl, "atreeid", "fit_atreeid", 
                                 "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_z_gsl <- subset(treeid_df2_z, !grepl("z|sigma", treeid))
aspp_df2_z_gsl   <- extract_params(df_fitgsl, "aspp", "fit_aspp", 
                                 "spp", "aspp\\[(\\d+)\\]")
site_df2_z_gsl   <- extract_params(df_fitgsl, "asite", "fit_a_site", 
                                 "site", "asite\\[(\\d+)\\]")

treeid_df2_z_gsl$treeid <- as.numeric(treeid_df2_z_gsl$treeid)
treeid_df2_z_gsl$treeid_name <- emp$treeid[match(treeid_df2_z_gsl$treeid, emp$treeid_num)]
bspp_df2_z_gsl$spp_name <- emp$latbi[match(bspp_df2_z_gsl$spp, emp$spp_num)]
site_df2_z_gsl$site_name <- emp$site[match(site_df2_z_gsl$site, emp$site_num)]
aspp_df2_z_gsl$spp_name <- emp$latbi[match(aspp_df2_z_gsl$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# SOS posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitsos <- readRDS("output/stanOutput/fitGrowthSOSZscored")

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
sigma_df2_z_sos  <- extract_params(df_fitsos, "sigma", "mean", "sigma")
bspp_df2_z_sos   <- extract_params(df_fitsos, "bsp", "fit_bspp", 
                                 "spp", "bsp\\[(\\d+)\\]")
treeid_df2_z_sos <- extract_params(df_fitsos, "atreeid", "fit_atreeid", 
                                 "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_z_sos <- subset(treeid_df2_z_sos, !grepl("z|sigma", treeid))
aspp_df2_z_sos   <- extract_params(df_fitsos, "aspp", "fit_aspp", 
                                 "spp", "aspp\\[(\\d+)\\]")
site_df2_z_sos   <- extract_params(df_fitsos, "asite", "fit_a_site", 
                                 "site", "asite\\[(\\d+)\\]")

treeid_df2_z_sos$treeid <- as.numeric(treeid_df2_z_sos$treeid)
treeid_df2_z_sos$treeid_name <- emp$treeid[match(treeid_df2_z_sos$treeid, emp$treeid_num)]
bspp_df2_z_sos$spp_name <- emp$latbi[match(bspp_df2_z_sos$spp, emp$spp_num)]
site_df2_z_sos$site_name <- emp$site[match(site_df2_z_sos$site, emp$site_num)]
aspp_df2_z_sos$spp_name <- emp$latbi[match(aspp_df2_z_sos$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# EOS posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fiteos <- readRDS("output/stanOutput/fitGrowthEOSZscored")

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
sigma_df2_z_eos  <- extract_params(df_fiteos, "sigma", "mean", "sigma")
bspp_df2_z_eos   <- extract_params(df_fiteos, "bsp", "fit_bspp", 
                                 "spp", "bsp\\[(\\d+)\\]")
treeid_df2_z_eos <- extract_params(df_fiteos, "atreeid", "fit_atreeid", 
                                 "treeid", "atreeid\\[(\\d+)\\]")
treeid_df2_z_eos <- subset(treeid_df2_z_eos, !grepl("z|sigma", treeid))
aspp_df2_z_eos   <- extract_params(df_fiteos, "aspp", "fit_aspp", 
                                 "spp", "aspp\\[(\\d+)\\]")
site_df2_z_eos   <- extract_params(df_fiteos, "asite", "fit_a_site", 
                                 "site", "asite\\[(\\d+)\\]")

treeid_df2_z_eos$treeid <- as.numeric(treeid_df2_z_eos$treeid)
treeid_df2_z_eos$treeid_name <- emp$treeid[match(treeid_df2_z_eos$treeid, emp$treeid_num)]
bspp_df2_z_eos$spp_name <- emp$latbi[match(bspp_df2_z_eos$spp, emp$spp_num)]
site_df2_z_eos$site_name <- emp$site[match(site_df2_z_eos$site, emp$site_num)]
aspp_df2_z_eos$spp_name <- emp$latbi[match(aspp_df2_z_eos$spp, emp$spp_num)]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Mu plots####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if(makeplots) {
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Define objects used throught the models #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
sitefull <- c(
  "GR" = "Dartmouth College (NH)",
  "HF" = "Harvard Forest (MA)",
  "SH" = "St-Hyppolyte (Qc)",
  "WM" = "White Mountains (NH)"
)

site_order <- c(
  "HF",
  "WM",
  "GR", 
  "SH")

locations <- data.frame(
  name = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hyppolyte (Qc)"),
  shortnames = c("HF", "WM", "GR", "SH"),
  lon = c(-72.2, -71, -70.66, -74.01),
  lat = c(42.6, 44.1, 44.9, 45.9)
)

# shapes for sites
my_shapes <- c( HF = 19, WM = 18, GR = 15, SH = 17)
sppcols <- c(wes_palette("AsteroidCity1"))[1:4]

subyvec <- vector()
for (i in 1:length(unique(emp$treeid_num))) {
  subyvec[i] <- paste("atreeid", "[",i,"]", sep = "")  
}
subyvec

# get the spp and site identities for each tree id
treeid_spp_site <- unique(emp[, c("treeid_num", "spp_num", "site_num",
                                  "treeid", "spp", "site", "latbi")])

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

n_spp <- nrow(bspp_df2_z)
n_site <- nrow(site_df2_z)
y_pos <- rev(1:n_spp)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### bspp ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/growthModelsMain/zscored/muALLbspp.jpeg",
     width = 1800, height = 2200, res = 300)

layout(matrix(c(
  1, 5,
  2, 5,
  3, 5,
  4, 5
), nrow = 4, byrow = TRUE),
widths = c(0.7, 0.4))

# Row 1: GDD
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_z$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "GDD standardized effect size", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_z$fit_bspp_per5,  y_pos, bspp_df2_z$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_z$fit_bspp_per25, y_pos, bspp_df2_z$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: GSL
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_z_gsl$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "GSL standardized effect size", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_z_gsl$fit_bspp_per5,  y_pos, bspp_df2_z_gsl$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_z_gsl$fit_bspp_per25, y_pos, bspp_df2_z_gsl$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3: SOS
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_z_sos$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "SOS standardized effect size", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_z_sos$fit_bspp_per5,  y_pos, bspp_df2_z_sos$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_z_sos$fit_bspp_per25, y_pos, bspp_df2_z_sos$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4: EOS
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_z_eos$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "EOS standardized effect size", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_z_eos$fit_bspp_per5,  y_pos, bspp_df2_z_eos$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_z_eos$fit_bspp_per25, y_pos, bspp_df2_z_eos$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("End of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Slot 5: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("center",
       legend = sapply(unique(bspp_df2_z$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(sppcols),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### aspp ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/growthModelsMain/zscored/muALLaspp.jpeg",
     width = 1800, height = 2200, res = 300)

layout(matrix(c(
  1, 5,
  2, 5,
  3, 5,
  4, 5
), nrow = 4, byrow = TRUE),
widths = c(0.7, 0.4))

# Row 1: GDD
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_z$fit_aspp, y_pos,
     xlim = c(-15, 15), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_z$fit_aspp_per5,  y_pos, aspp_df2_z$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_z$fit_aspp_per25, y_pos, aspp_df2_z$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: GSL
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_z_gsl$fit_aspp, y_pos,
     xlim = c(-15, 15), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_z_gsl$fit_aspp_per5,  y_pos, aspp_df2_z_gsl$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_z_gsl$fit_aspp_per25, y_pos, aspp_df2_z_gsl$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3: SOS
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_z_sos$fit_aspp, y_pos,
     xlim = c(-15, 15),ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "", 
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_z_sos$fit_aspp_per5,  y_pos, aspp_df2_z_sos$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_z_sos$fit_aspp_per25, y_pos, aspp_df2_z_sos$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4: EOS
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_z_eos$fit_aspp, y_pos,
     xlim = c(-15, 15), ylim = c(0.5, n_spp + 0.5), 
     xlab = "Ring width intercept values (mm)", ylab = "", 
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_z_eos$fit_aspp_per5,  y_pos, aspp_df2_z_eos$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_z_eos$fit_aspp_per25, y_pos, aspp_df2_z_eos$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("End of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Slot 5: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("center",
       legend = sapply(unique(aspp_df2_z$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(sppcols),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### asite ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/growthModelsMain/zscored/muALLasite.jpeg",
     width = 1800, height = 2200, res = 300)

layout(matrix(c(
  1, 5,
  2, 5,
  3, 5,
  4, 5
), nrow = 4, byrow = TRUE),
widths = c(0.7, 0.4))

y_pos <- match(site_df2_z$site_name, site_order)
site_color_map <- setNames(c(wes_palette("Darjeeling1"))[1:4], site_order)
sitecolors <- site_color_map[site_df2_z$site_name]

site_df2_z$lat <- locations$lat[match(site_df2_z$site_name, locations$shortnames)]
site_df2_z_gsl$lat <- locations$lat[match(site_df2_z_gsl$site_name, locations$shortnames)]
site_df2_z_sos$lat <- locations$lat[match(site_df2_z_sos$site_name, locations$shortnames)]
site_df2_z_eos$lat <- locations$lat[match(site_df2_z_eos$site_name, locations$shortnames)]

lat_labels <- locations$lat[match(site_order, locations$shortnames)]


# Row 1: GDD
par(mar = c(5, 8, 2, 2))
plot(site_df2_z$fit_a_site, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "Latitude",
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE, panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE,      
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(site_df2_z$fit_a_site_per5,  y_pos, site_df2_z$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_z$fit_a_site_per25, y_pos, site_df2_z$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)

mtext("Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: GSL
par(mar = c(5, 8, 2, 2))
plot(site_df2_z_gsl$fit_a_site, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "Latitude",
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE, panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE,      
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(site_df2_z_gsl$fit_a_site_per5,  y_pos, site_df2_z_gsl$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_z_gsl$fit_a_site_per25, y_pos, site_df2_z_gsl$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)

mtext("Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3: SOS
par(mar = c(5, 8, 2, 2))
plot(site_df2_z_sos$fit_a_site, y_pos,
     xlim = c(-2, 2),ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "Latitude", 
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE, panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE,      
     panel.first = abline(v = 0, lty = 2, col = "black"),
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(site_df2_z_sos$fit_a_site_per5,  y_pos, site_df2_z_sos$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_z_sos$fit_a_site_per25, y_pos, site_df2_z_sos$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)

mtext("Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4: EOS
par(mar = c(5, 8, 2, 2))
plot(site_df2_z_eos$fit_a_site, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5), 
     xlab = "Ring width intercept values (mm)", ylab = "Latitude", 
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE, panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE,      
     panel.first = abline(v = 0, lty = 2, col = "black"),
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(site_df2_z_eos$fit_a_site_per5,  y_pos, site_df2_z_eos$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_z_eos$fit_a_site_per25, y_pos, site_df2_z_eos$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)

mtext("End of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Slot 5: Site legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("center",
       legend = rev(locations$name[match(site_order, locations$shortnames)]),
       col    = rev(site_color_map[site_order]),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Site", title.font = 2)
dev.off()

}
