# Wildchrokie model
# CRD 14 December 2025

# Goal: Plot model output because modelGrowthGDD is becoming too long and messy
  
if (length(grep("christophe_rouleau-desrochers", getwd())) > 0) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if (length(grep("lizzie", getwd())) > 0) {
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
} else  {
  setwd("/home/crouleau/wildchrokie/analyses")
}

# source model code
source("rcode/growthModelsMain.R")

library(ggplot2)

# flags
makeplots <- TRUE
# interceptmuplots <- TRUE

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

# Load parameter summaries generated in growthModelsMain.R ####
sigma_df2  <- read.csv("output/GM_GDDparam_sigma.csv")
bspp_df2   <- read.csv("output/GM_GDDparam_bspp.csv")
treeid_df2 <- read.csv("output/GM_GDDparam_treeid.csv")
aspp_df2   <- read.csv("output/GM_GDDparam_aspp.csv")
site_df2   <- read.csv("output/GM_GDDparam_site.csv")

treeid_df2$treeid <- as.numeric(treeid_df2$treeid)  
treeid_df2$treeid_name <- emp$treeid[match(treeid_df2$treeid, emp$treeid_num)]
bspp_df2$spp_name <- emp$latbi[match(bspp_df2$spp, emp$spp_num)]
site_df2$site_name <- emp$site[match(site_df2$site, emp$site_num)]
aspp_df2$spp_name <- emp$latbi[match(aspp_df2$spp, emp$spp_num)]

# GSL
sigma_df2_gsl  <- read.csv("output/GM_GSLparam_sigma.csv")
bspp_df2_gsl   <- read.csv("output/GM_GSLparam_bspp.csv")
treeid_df2_gsl <- read.csv("output/GM_GSLparam_treeid.csv")
aspp_df2_gsl   <- read.csv("output/GM_GSLparam_aspp.csv")
site_df2_gsl   <- read.csv("output/GM_GSLparam_site.csv")

treeid_df2_gsl$treeid <- as.numeric(treeid_df2_gsl$treeid)
treeid_df2_gsl$treeid_name <- emp$treeid[match(treeid_df2_gsl$treeid, emp$treeid_num)]
bspp_df2_gsl$spp_name <- emp$latbi[match(bspp_df2_gsl$spp, emp$spp_num)]
site_df2_gsl$site_name <- emp$site[match(site_df2_gsl$site, emp$site_num)]
aspp_df2_gsl$spp_name <- emp$latbi[match(aspp_df2_gsl$spp, emp$spp_num)]

# SOS 
sigma_df2_sos  <- read.csv("output/GM_SOSparam_sigma.csv")
bspp_df2_sos   <- read.csv("output/GM_SOSparam_bspp.csv")
treeid_df2_sos <- read.csv("output/GM_SOSparam_treeid.csv")
aspp_df2_sos   <- read.csv("output/GM_SOSparam_aspp.csv")
site_df2_sos   <- read.csv("output/GM_SOSparam_site.csv")

treeid_df2_sos$treeid <- as.numeric(treeid_df2_sos$treeid)
treeid_df2_sos$treeid_name <- emp$treeid[match(treeid_df2_sos$treeid, emp$treeid_num)]
bspp_df2_sos$spp_name <- emp$latbi[match(bspp_df2_sos$spp, emp$spp_num)]
site_df2_sos$site_name <- emp$site[match(site_df2_sos$site, emp$site_num)]
aspp_df2_sos$spp_name <- emp$latbi[match(aspp_df2_sos$spp, emp$spp_num)]

# EOS
sigma_df2_eos  <- read.csv("output/GM_EOSparam_sigma.csv")
bspp_df2_eos   <- read.csv("output/GM_EOSparam_bspp.csv")
treeid_df2_eos <- read.csv("output/GM_EOSparam_treeid.csv")
aspp_df2_eos   <- read.csv("output/GM_EOSparam_aspp.csv")
site_df2_eos   <- read.csv("output/GM_EOSparam_site.csv")

treeid_df2_eos$treeid <- as.numeric(treeid_df2_eos$treeid)
treeid_df2_eos$treeid_name <- emp$treeid[match(treeid_df2_eos$treeid, emp$treeid_num)]
bspp_df2_eos$spp_name <- emp$latbi[match(bspp_df2_eos$spp, emp$spp_num)]
site_df2_eos$site_name <- emp$site[match(site_df2_eos$site, emp$site_num)]
aspp_df2_eos$spp_name <- emp$latbi[match(aspp_df2_eos$spp, emp$spp_num)]


if(makeplots) {
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GDD posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgdd <- readRDS("output/stanOutput/fitGrowthGDD")

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# GSL posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitgsl <- readRDS("output/stanOutput/fitGrowthGSL")

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# SOS posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fitsos <- readRDS("output/stanOutput/fitGrowthSOS")

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# EOS posterior recovery ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
fiteos <- readRDS("output/stanOutput/fitGrowthEOS")

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot lines with quantiles ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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

treeidvecnum <- unique(treeid_spp_site$treeid_num)
treeidvecname <- treeid_spp_site$treeid

sppvecnum <- 1:4
sppvecname <- unique(treeid_spp_site$latbi)

n_spp <- nrow(bspp_df2)
n_site <- nrow(site_df2)
y_pos <- rev(1:n_spp)

# y axis for mean plots
ylimline <- c(-1, 3)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GDD: Prep posterior reconstruction #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# recover full intercepts for each tree id
fullintercept_cols <- grep("^fullintercept", colnames(df_fitgdd), value = TRUE)
fullintercept <- df_fitgdd[, fullintercept_cols]
colnames(fullintercept) <- 1:ncol(fullintercept)

# recover each slope
treeid_slope_cols <- grep("^treeid_slope", colnames(df_fitgdd), value = TRUE)
treeid_bspp <- df_fitgdd[, treeid_slope_cols] / wcgddscale
colnames(treeid_bspp) <- 1:ncol(treeid_bspp)

# recover sim ypred for each gdd X tree id for each iterations
y_post_array <- extract(fitgdd, "y_post")$y_post

gddseq <- dgdd$gddseq

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GDD: per treeid, facet #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# PDF output
pdf(file = "figures/growthModelsMain/growthModelSlopesperTreeid.pdf", width = 10, height = 8)

# Layout: 2 rows × 2 columns per page
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

for (i in seq_along(treeidvecnum)) { # i = 6
  tree_col <- as.character(treeidvecnum[i])
  tree_col_name <- as.character(treeidvecname[i])
  
  # new extraction: y_post_array is [n_draws, Ngddseq, Ntreeid]
  # slice out this tree, transpose to [Ngddseq, n_draws] to match old y_post shape
  y_post <- t(y_post_array[, , i])
  
  tree_id_num <- as.integer(tree_col)
  emp_treeid <- emp[emp$treeid_num == tree_id_num, ]
  
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  plot(emp_treeid$pgsGDD5, emp_treeid$loglength, type = "n", 
       ylim = range(c(emp_treeid$loglength, y_low, y_high), na.rm = TRUE),
       xlab = "Primary growing season GDD", ylab = "Ring width (mm)",
       main = tree_col_name)
  
  spp_id <- treeid_spp_site$spp_num[match(tree_id_num, treeid_spp_site$treeid_num)]
  line_col <- sppcols[spp_id]
  
  polygon(c(gddseq, rev(gddseq)), 
          c(y_low, rev(y_high)),
          col = adjustcolor(line_col, alpha.f = 0.3), 
          border = NA)
  
  lines(gddseq, y_mean, col = line_col, lwd = 2)
  
  points(emp_treeid$pgsGDD5, emp_treeid$loglength, pch = 16, cex = 2, col = line_col)
}
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GDD: per Spp, facet #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
spp_post_array <- extract(fitgdd, "spp_post")$spp_post
# dimensions: [n_draws, Ngddseq, Nspp]

# jpeg output
jpeg(
  filename = "figures/growthModelsMain/growthModelSlopesperSppFacet.jpeg",
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
  y_post <- t(spp_post_array[, , i])
  
  # summaries
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  # species-specific ylim
  ylim_spp <- range(c(emp_spp$loglength, y_low, y_high), na.rm = TRUE)
  
  plot(emp_spp$pgsGDD5, emp_spp$loglength,
       type = "n",
       ylim = ylimline,
       xlab = "Primary growing season GDD",
       ylab = "Ring width (mm) log",
            main = bquote(italic(.(spp_column_name))),
       frame = FALSE)
  
  # add panel letter 
  mtext(paste0("(", letters[i], ")"), 
        side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  # color
  line_col <- sppcols[spp_num]
  
  polygon(
    c(gddseq, rev(gddseq)),
    c(y_low, rev(y_high)),
    col = adjustcolor(line_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(gddseq, y_mean, col = line_col, lwd = 2)
  
  points(
    emp_spp$pgsGDD5,
    emp_spp$loglength,
    pch = my_shapes[emp_spp$site],
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
  filename = "figures/growthModelsMain/growthModelSlopesperSpp.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)

# Layout: 2 rows × 2 columns per page
par(mar = c(4, 4, 2, 1))

plot(emp$pgsGDD5, dgdd$y, type = "n", 
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Primary growing season GDD", ylab = "Ring width (mm)",
     main = "species growth responses")

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  y_post <- t(spp_post_array[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  # shaded interval
  polygon(c(gddseq, rev(gddseq)), 
          c(y_low, # lower interval
            rev(y_high)), # high interval
          col = adjustcolor(line_col, alpha.f = 0.2), 
          border = NA)
  
  # mean line
  lines(gddseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]

  points(
    emp_spp$pgsGDD5,
    emp_spp$loglength,
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
# recover full intercepts for each tree id
fullintercept_cols_gsl <- grep("^fullintercept", colnames(df_fitgsl), value = TRUE)
fullintercept_gsl <- df_fitgsl[, fullintercept_cols_gsl]
colnames(fullintercept_gsl) <- 1:ncol(fullintercept_gsl)

# recover each slope
treeid_slope_cols_gsl <- grep("^treeid_slope", colnames(df_fitgsl), value = TRUE)
treeid_bspp_gsl <- df_fitgsl[, treeid_slope_cols_gsl] / gslscale
colnames(treeid_bspp_gsl) <- 1:ncol(treeid_bspp_gsl)

# recover sim ypred for each gsl X tree id for each iterations
y_post_array_gsl <- extract(fitgsl, "y_post")$y_post

# posterior array for each species. This recovers the simulated full intercept + bspp*predictor, with sigma_y
spp_post_array_gsl <- extract(fitgsl, "spp_post")$spp_post

gslseq <- dgsl$gslseq

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### GSL: per Spp, facet #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# jpeg output
jpeg(
  filename = "figures/growthModelsMain/growthModelSlopesperSppFacetGSL.jpeg",
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
  
  y_post_gsl <- t(spp_post_array_gsl[, , i])
  
  # summaries
  y_mean_gsl <- apply(y_post_gsl, 1, mean)
  y_low_gsl  <- apply(y_post_gsl, 1, quantile, 0.25)
  y_high_gsl <- apply(y_post_gsl, 1, quantile, 0.75)
  
  # species-specific ylim
  # ylim_spp <- range(c(emp_spp$loglength, y_low, y_high), na.rm = TRUE)
  
  plot(emp_spp$pgsGSL, emp_spp$loglength,
       type = "n",
       ylim = ylimline,
       xlab = "Primary growing season GSL",
       ylab = "Ring width (mm)",
       main = bquote(italic(.(spp_column_name))),
       frame = FALSE)
  
  # add panel letter 
  mtext(paste0("(", letters[i], ")"), 
        side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  polygon(
    c(gslseq, rev(gslseq)),
    c(y_low_gsl, rev(y_high_gsl)),
    col = adjustcolor(line_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(gslseq, y_mean_gsl, col = line_col, lwd = 2)
  
  points(
    emp_spp$pgsGSL,
    emp_spp$loglength,
    pch = my_shapes[emp_spp$site],
    cex = 1,
    col = line_col
  )
}

dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# SOS: prep posterior reconstruction ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# recover full intercepts for each tree id
fullintercept_cols_sos <- grep("^fullintercept", colnames(df_fitsos), value = TRUE)
fullintercept_sos <- df_fitsos[, fullintercept_cols_sos]
colnames(fullintercept_sos) <- 1:ncol(fullintercept_sos)

# recover each slope
treeid_slope_cols_sos <- grep("^treeid_slope", colnames(df_fitsos), value = TRUE)
treeid_bspp_sos <- df_fitsos[, treeid_slope_cols_sos] / sosscale
colnames(treeid_bspp_sos) <- 1:ncol(treeid_bspp_sos)

# recover sim ypred for each sos X tree id for each iterations
y_post_array_sos <- extract(fitsos, "y_post")$y_post

# posterior array for each species. This recovers the simulated full intercept + bspp*predictor, with sigma_y
spp_post_array_sos <- extract(fitsos, "spp_post")$spp_post

sosseq <- dsos$sosseq

# jpeg output
jpeg(
  filename = "figures/growthModelsMain/growthModelSlopesperSppFacetSOS.jpeg",
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
  
  y_post_sos <- t(spp_post_array_sos[, , i])
  
  # summaries
  y_mean_sos <- apply(y_post_sos, 1, mean)
  y_low_sos  <- apply(y_post_sos, 1, quantile, 0.25)
  y_high_sos <- apply(y_post_sos, 1, quantile, 0.75)
  
  # species-specific ylim
  # ylim_spp <- range(c(emp_spp$loglength, y_low, y_high), na.rm = TRUE)
  
  plot(emp_spp$leafout, emp_spp$loglength,
       type = "n",
       ylim = ylimline,
       xlab = "Leafout day of year",
       ylab = "Ring width (mm)",
       main = bquote(italic(.(spp_column_name))),
       frame = FALSE)
  
  # add panel letter 
  mtext(paste0("(", letters[i], ")"), 
        side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  polygon(
    c(sosseq, rev(sosseq)),
    c(y_low_sos, rev(y_high_sos)),
    col = adjustcolor(line_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(sosseq, y_mean_sos, col = line_col, lwd = 2)
  
  points(
    emp_spp$leafout,
    emp_spp$loglength,
    pch = my_shapes[emp_spp$site],
    cex = 1,
    col = line_col
  )
}

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### SOS: per Spp, non-facetted #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# PDF output
jpeg(
  filename = "figures/growthModelsMain/growthModelSlopesperSppNoFacetSOS.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)

# Layout: 2 rows × 2 columns per page
par(mar = c(4, 4, 2, 1))

plot(emp$leafout, dsos$y, type = "n", frame = FALSE,
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Leafout day of year", ylab = "Ring width (mm)",
     main = "")

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  
  y_post_sos <- t(spp_post_array_sos[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post_sos, 1, mean)
  y_low  <- apply(y_post_sos, 1, quantile, 0.25)
  y_high <- apply(y_post_sos, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  # shaded interval
  polygon(c(sosseq, rev(sosseq)),
          c(y_low, # lower interval
            rev(y_high)), # high interval
          col = adjustcolor(line_col, alpha.f = 0.1),
          border = NA)
  
  # mean line
  lines(sosseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]
  
  # points(
  #   x = emp_spp$leafout, y = emp_spp$loglength,
  #   pch = 16, cex = 1, col = line_col)
  
  # legend(
  #   "topleft", legend = sppvecname,
  #   col = sppcols[as.integer(sppvecnum)], lwd = 2, pch = 16, bty = "n", cex = 1.5
  # )
  
}
dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# EOS: prep posterior reconstruction ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# recover full intercepts for each tree id
fullintercept_cols_eos <- grep("^fullintercept", colnames(df_fiteos), value = TRUE)
fullintercept_eos <- df_fiteos[, fullintercept_cols_eos]
colnames(fullintercept_eos) <- 1:ncol(fullintercept_eos)

# recover each slope
treeid_slope_cols_eos <- grep("^treeid_slope", colnames(df_fiteos), value = TRUE)
treeid_bspp_eos <- df_fiteos[, treeid_slope_cols_eos] / eosscale
colnames(treeid_bspp_eos) <- 1:ncol(treeid_bspp_eos)

# recover sim ypred for each eos X tree id for each iterations
y_post_array_eos <- extract(fiteos, "y_post")$y_post

# posterior array for each species. This recovers the simulated full intercept + bspp*predictor, with sigma_y
spp_post_array_eos <- extract(fiteos, "spp_post")$spp_post

eosseq <- deos$eosseq

# jpeg output
jpeg(
  filename = "figures/growthModelsMain/growthModelSlopesperSppFacetEOS.jpeg",
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
  
  y_post_eos <- t(spp_post_array_eos[, , i])
  
  # summaries
  y_mean_eos <- apply(y_post_eos, 1, mean)
  y_low_eos  <- apply(y_post_eos, 1, quantile, 0.25)
  y_high_eos <- apply(y_post_eos, 1, quantile, 0.75)
  
  # species-specific ylim
  # ylim_spp <- range(c(emp_spp$loglength, y_low, y_high), na.rm = TRUE)
  
  plot(emp_spp$budset, emp_spp$loglength,
       type = "n",
       ylim = ylimline,
       xlab = "Budset day of year",
       ylab = "Ring width (mm)",
       main = bquote(italic(.(spp_column_name))),
       frame = FALSE)
  
  # add panel letter 
  mtext(paste0("(", letters[i], ")"), 
        side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  polygon(
    c(eosseq, rev(eosseq)),
    c(y_low_eos, rev(y_high_eos)),
    col = adjustcolor(line_col, alpha.f = 0.3),
    border = NA
  )
  
  lines(eosseq, y_mean_eos, col = line_col, lwd = 2)
  
  points(
    emp_spp$budset,
    emp_spp$loglength,
    pch = my_shapes[emp_spp$site],
    cex = 1,
    col = line_col
  )
}

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### EOS: per Spp, non-facetted #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# PDF output
jpeg(
  filename = "figures/growthModelsMain/growthModelSlopesperSppNoFacetEOS.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)

# Layout: 2 rows × 2 columns per page
par(mar = c(4, 4, 2, 1))

plot(emp$budset, deos$y, type = "n", frame = FALSE,
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Budset day of year", ylab = "Ring width (mm)",
     main = "")

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  
  y_post_eos <- t(spp_post_array_eos[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post_eos, 1, mean)
  y_low  <- apply(y_post_eos, 1, quantile, 0.25)
  y_high <- apply(y_post_eos, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  polygon(c(eosseq, rev(eosseq)),
          c(y_low, rev(y_high)), # low and high interval
          col = adjustcolor(line_col, alpha.f = 0.1), border = NA)
  
  lines(eosseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]
  
  # points(
  #   x = emp_spp$budset, y = emp_spp$loglength,
  #   pch = 16, cex = 1, col = line_col)
  
  # legend(
  #   "topleft", legend = sppvecname,
  #   col = sppcols[as.integer(sppvecnum)], lwd = 2, pch = 16, bty = "n", cex = 1.5
  # )
  
}
dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Mu plots #####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
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

# get posterior means and quantiles
treeid_df4 <- data.frame(
  treeid = character(ncol(fullintercept)),
  fit_atreeid = numeric(ncol(fullintercept)),  
  fit_atreeid_per5 = NA, 
  fit_atreeid_per25 = NA,
  fit_atreeid_per75 = NA,
  fit_atreeid_per95 = NA
)
for (i in 1:ncol(fullintercept)) { # i = 1
  treeid_df4$treeid[i] <- colnames(fullintercept)[i]         
  treeid_df4$fit_atreeid[i] <- round(mean(fullintercept[[i]]),3)  
  treeid_df4$fit_atreeid_per5[i] <- round(quantile(fullintercept[[i]], probs = 0.05), 3)
  treeid_df4$fit_atreeid_per25[i] <- round(quantile(fullintercept[[i]], probs = 0.25), 3)
  treeid_df4$fit_atreeid_per75[i] <- round(quantile(fullintercept[[i]], probs = 0.75), 3)
  treeid_df4$fit_atreeid_per95[i] <- round(quantile(fullintercept[[i]], probs = 0.95), 3)
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

# col
my_colors <- c(
  "Alnus incana" = wes_palette("AsteroidCity1")[1],
  "Betula alleghaniensis" = wes_palette("AsteroidCity1")[2],
  "Betula papyrifera" = wes_palette("AsteroidCity1")[3],
  "Betula populifolia" = wes_palette("AsteroidCity1")[4]
)

# open device
pdf(
  file = "figures/growthModelsMain/meanPlotGrowthGDD_treeidBYspp.pdf",
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
# Combined mu plots bspp (GDD / GSL / SOS / EOS) ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### bspp ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/growthModelsMain/muALLbspp.jpeg",
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
plot(bspp_df2$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width (mm) change in averaged GDD of 10 spring days", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2$fit_bspp_per5,  y_pos, bspp_df2$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2$fit_bspp_per25, y_pos, bspp_df2$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)

mtext("Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: GSL
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_gsl$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width (mm) change per 10 days of GSL", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_gsl$fit_bspp_per5,  y_pos, bspp_df2_gsl$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_gsl$fit_bspp_per25, y_pos, bspp_df2_gsl$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3: SOS
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_sos$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width (mm) change per 5 days of leafout", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_sos$fit_bspp_per5,  y_pos, bspp_df2_sos$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_sos$fit_bspp_per25, y_pos, bspp_df2_sos$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4: EOS
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_eos$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width (mm) change per 10 days of budset", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_eos$fit_bspp_per5,  y_pos, bspp_df2_eos$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_eos$fit_bspp_per25, y_pos, bspp_df2_eos$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("End of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Slot 5: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("center",
       legend = sapply(unique(bspp_df2$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(sppcols),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### bspp with lines #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
jpeg(file = "figures/growthModelsMain/muALLbspp.jpeg",
     width = 2800, height = 2800, res = 300)

layout(matrix(c(
  1, 5, 9,
  2, 6, 9,
  3, 7, 9,
  4, 8, 9
), nrow = 4, byrow = TRUE),
widths = c(1.3, 1, 0.5))

# Row 1, Col 1, Slot 5 : GDD
par(mar = c(5, 8, 2, 2))
plot(bspp_df2$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Log ring width (mm) change in averaged GDD of 10 spring days", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2$fit_bspp_per5,  y_pos, bspp_df2$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2$fit_bspp_per25, y_pos, bspp_df2$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(a) Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2, Col 1, Slot 6 : GSL
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_gsl$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Log ring width (mm) change per 10 days of GSL", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_gsl$fit_bspp_per5,  y_pos, bspp_df2_gsl$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_gsl$fit_bspp_per25, y_pos, bspp_df2_gsl$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(b) Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3, Col 1, Slot 7 : SOS
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_sos$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Log ring width (mm) change per 5 days of leafout", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_sos$fit_bspp_per5,  y_pos, bspp_df2_sos$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_sos$fit_bspp_per25, y_pos, bspp_df2_sos$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(c) Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4, Col 1, Slot 8 : EOS
par(mar = c(5, 8, 2, 2))
plot(bspp_df2_eos$fit_bspp, y_pos,
     xlim = c(-0.8, 0.8), ylim = c(0.5, n_spp + 0.5),
     xlab = "Log ring width (mm) change per 10 days of budset", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(bspp_df2_eos$fit_bspp_per5,  y_pos, bspp_df2_eos$fit_bspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(bspp_df2_eos$fit_bspp_per25, y_pos, bspp_df2_eos$fit_bspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(d) End of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 1, Col 2, Slot 5 : GDD
par(mar = c(5, 8, 2, 2))
plot(emp$pgsGDD5, dgdd$y, type = "n", frame = FALSE,
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Growing season growing degree days (GDD)", ylab = "Log ring width (mm)",
     main = "")
mtext("(e)", side = 3, adj = 0, font = 2, cex = 0.9)

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  
  y_post <- t(spp_post_array[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  polygon(c(gddseq, rev(gddseq)),
          c(y_low, rev(y_high)), # low and high interval
          col = adjustcolor(line_col, alpha.f = 0.1), border = NA)
  
  lines(gddseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]
}
# Row 2, Col 2, Slot 6 : GSL
par(mar = c(5, 8, 2, 2))
plot(emp$pgsGSL, dgsl$y, type = "n", frame = FALSE,
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Growing season length (days)", ylab = "Log ring width (mm)",
     main = "")
mtext("(f)", side = 3, adj = 0, font = 2, cex = 0.9)

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  
  y_post <- t(spp_post_array_gsl[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  polygon(c(gslseq, rev(gslseq)),
          c(y_low, rev(y_high)), # low and high interval
          col = adjustcolor(line_col, alpha.f = 0.1), border = NA)
  
  lines(gslseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]
}

# Row 3, Col 2, Slot 7 : SOS
par(mar = c(5, 8, 2, 2))
plot(emp$leafout, dsos$y, type = "n", frame = FALSE,
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Leafout day of year", ylab = "Log ring width (mm)",
     main = "")
mtext("(g)", side = 3, adj = 0, font = 2, cex = 0.9)

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  y_post <- t(spp_post_array_sos[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  polygon(c(sosseq, rev(sosseq)),
          c(y_low, rev(y_high)), # low and high interval
          col = adjustcolor(line_col, alpha.f = 0.1), border = NA)
  
  lines(sosseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]
}

# Row 4, Col 2, Slot 8 : EOS
par(mar = c(5, 8, 2, 2))
plot(emp$budset, deos$y, type = "n", frame = FALSE,
     ylim = range(min(emp$loglength), max(emp$loglength)), 
     xlab = "Budset day of year", ylab = "Log ring width (mm)",
     main = "")
mtext("(h)", side = 3, adj = 0, font = 2, cex = 0.9)

# Loop over trees again to plot each tree individually
for (i in seq_along(sppvecnum)) { # i = 1
  spp_column <- as.character(sppvecnum[i])
  spp_column_name <- as.character(sppvecname[i])
  y_post <- t(spp_post_array_eos[, , i])
  
  # color line by spp
  spp_num <- as.integer(spp_column)
  
  # calculate mean and 50% credible interval (25%-75%)
  y_mean <- apply(y_post, 1, mean)
  y_low  <- apply(y_post, 1, quantile, 0.25)
  y_high <- apply(y_post, 1, quantile, 0.75)
  
  line_col <- sppcols[spp_num]
  
  polygon(c(eosseq, rev(eosseq)),
          c(y_low, rev(y_high)), # low and high interval
          col = adjustcolor(line_col, alpha.f = 0.1), border = NA)
  
  lines(eosseq, y_mean,
        col = line_col,
        lwd = 2)
  
  emp_spp <- emp[emp$spp_num == spp_num, ]
}

# Slot 9: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("center",
       legend = sapply(unique(bspp_df2$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(sppcols),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)

dev.off()


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### aspp ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/growthModelsMain/muALLaspp.jpeg",
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
plot(aspp_df2$fit_aspp, y_pos,
     xlim = c(-15, 15), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "",
     yaxt = "n", pch = 15, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2$fit_aspp_per5,  y_pos, aspp_df2$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2$fit_aspp_per25, y_pos, aspp_df2$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(a) Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: GSL
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_gsl$fit_aspp, y_pos,
     xlim = c(-15, 15), ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "",
     yaxt = "n", pch = 15, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_gsl$fit_aspp_per5,  y_pos, aspp_df2_gsl$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_gsl$fit_aspp_per25, y_pos, aspp_df2_gsl$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(b) Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3: SOS
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_sos$fit_aspp, y_pos,
     xlim = c(-15, 15),ylim = c(0.5, n_spp + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "", 
     yaxt = "n", pch = 15, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_sos$fit_aspp_per5,  y_pos, aspp_df2_sos$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_sos$fit_aspp_per25, y_pos, aspp_df2_sos$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(c) Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4: EOS
par(mar = c(5, 8, 2, 2))
plot(aspp_df2_eos$fit_aspp, y_pos,
     xlim = c(-15, 15), ylim = c(0.5, n_spp + 0.5), 
     xlab = "Ring width intercept values (mm)", ylab = "", 
     yaxt = "n", pch = 15, cex = 2, col = sppcols, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
segments(aspp_df2_eos$fit_aspp_per5,  y_pos, aspp_df2_eos$fit_aspp_per95, y_pos,
         col = sppcols, lwd = 1.5)
segments(aspp_df2_eos$fit_aspp_per25, y_pos, aspp_df2_eos$fit_aspp_per75, y_pos,
         col = sppcols, lwd = 3)
mtext("(d) End of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Slot 5: species legend
par(mar = c(1, 1, 1, 1))
plot.new()
legend("center",
       legend = sapply(unique(aspp_df2$spp_name), 
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(sppcols),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### asite ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/growthModelsMain/muALLasite.jpeg",
     width = 1800, height = 2200, res = 300)

layout(matrix(c(
  1, 5,
  2, 5,
  3, 5,
  4, 5
), nrow = 4, byrow = TRUE),
widths = c(0.7, 0.4))

y_pos <- match(site_df2$site_name, site_order)
site_color_map <- setNames(c(wes_palette("Darjeeling1"))[1:4], site_order)
sitecolors <- site_color_map[site_df2$site_name]

site_df2$lat <- locations$lat[match(site_df2$site_name, locations$shortnames)]
site_df2_gsl$lat <- locations$lat[match(site_df2_gsl$site_name, locations$shortnames)]
site_df2_sos$lat <- locations$lat[match(site_df2_sos$site_name, locations$shortnames)]
site_df2_eos$lat <- locations$lat[match(site_df2_eos$site_name, locations$shortnames)]

lat_labels <- locations$lat[match(site_order, locations$shortnames)]


# Row 1: GDD
par(mar = c(5, 8, 2, 2))
plot(site_df2$fit_a_site, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "Latitude",
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE)
segments(site_df2$fit_a_site_per5,  y_pos, site_df2$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2$fit_a_site_per25, y_pos, site_df2$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)
mtext("Growing degree days", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 2: GSL
par(mar = c(5, 8, 2, 2))
plot(site_df2_gsl$fit_a_site, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "Latitude",
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE)
segments(site_df2_gsl$fit_a_site_per5,  y_pos, site_df2_gsl$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_gsl$fit_a_site_per25, y_pos, site_df2_gsl$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)
mtext("Growing season length", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 3: SOS
par(mar = c(5, 8, 2, 2))
plot(site_df2_sos$fit_a_site, y_pos,
     xlim = c(-2, 2),ylim = c(0.5, n_site + 0.5),
     xlab = "Ring width intercept values (mm)", ylab = "Latitude", 
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE,
     panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE)
segments(site_df2_sos$fit_a_site_per5,  y_pos, site_df2_sos$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_sos$fit_a_site_per25, y_pos, site_df2_sos$fit_a_site_per75, y_pos,
         col = sitecolors, lwd = 3)
mtext("Start of season", side = 3, adj = 0, font = 2, cex = 0.9)

# Row 4: EOS
par(mar = c(5, 8, 2, 2))
plot(site_df2_eos$fit_a_site, y_pos,
     xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5), 
     xlab = "Ring width intercept values (mm)", ylab = "Latitude", 
     yaxt = "n", pch = 16, cex = 2, col = sitecolors, frame.plot = FALSE, 
     panel.first = abline(v = 0, lty = 2, col = "black"))
axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE)
segments(site_df2_eos$fit_a_site_per5,  y_pos, site_df2_eos$fit_a_site_per95, y_pos,
         col = sitecolors, lwd = 1.5)
segments(site_df2_eos$fit_a_site_per25, y_pos, site_df2_eos$fit_a_site_per75, y_pos,
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### asite with map ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)

# set full names of sites
locations <- data.frame(
  name       = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hyppolyte (Qc)"),
  shortnames = c("HF", "WM", "GR", "SH"),   
  Longitude  = c(-72.2,  -71.0, -70.66, -74.01),
  Latitude   = c( 42.55,  44.11,  44.92,  45.98)
)
site_order <- locations$shortnames[order(locations$Latitude)]     
site_color_map <- setNames(wes_palette("Darjeeling1")[1:4], site_order)

# order same as figure
locations2 <- locations[order(locations$Latitude), ]              
locations2$col <- wes_palette("Darjeeling1")[1:4]     

# site map 
n_site     <- length(site_order)
y_pos      <- match(site_df2$site_name, site_order)
sitecolors <- site_color_map[site_df2$site_name]
lat_labels <- locations$Latitude[match(site_order, locations$shortnames)]

world <- ne_countries(scale = "medium", returnclass = "sf")
lat_min <- 41; lat_max <- 48
lon_min <- -78; lon_max <- -63

special_point <- data.frame(
  name = "Arnold Arboretum of\nHarvard University (MA)",
  Longitude  = -71.13358611669867,
  Latitude  =  42.29601035316377
)

special_sf <- st_as_sf(special_point, coords = c("Longitude", "Latitude"), crs = 4326)
points_sf  <- st_as_sf(locations2, coords = c("Longitude", "Latitude"), crs = 4326)


map_plot <- ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = points_sf, color = locations2$col, size = 4) +
  geom_text(data = locations2,
            aes(x = Longitude, y = Latitude, label = name),
            nudge_y = 0.35, size = 4.5, fontface = "bold") +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  theme_minimal() +
  theme(
    strip.text        = element_blank(),
    legend.key.height = unit(1.5, "lines"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# mu plot for asite
forest_grob <- as_grob(function() {
  par(mar = c(7, 5, 5, 0.5))
  
  plot(site_df2$fit_a_site, y_pos,
       xlim = c(-2, 2), ylim = c(0.5, n_site + 0.5),
       xlab = "Ring width intercept values (mm)", ylab = "Latitude",
       yaxt = "n", pch = 16, cex = 2, col = sitecolors,
       frame.plot = TRUE,
       panel.first = abline(v = 0, lty = 2, col = "black"))
  axis(2, at = 1:n_site, labels = lat_labels, las = 2, tick = TRUE)
  segments(site_df2$fit_a_site_per5,  y_pos,
           site_df2$fit_a_site_per95, y_pos,
           col = sitecolors, lwd = 1.5)
  segments(site_df2$fit_a_site_per25, y_pos,
           site_df2$fit_a_site_per75, y_pos,
           col = sitecolors, lwd = 3)
})

combined <- plot_grid(forest_grob, map_plot, ncol = 2, rel_widths = c(0.4, 0.7))

combined_labeled <- ggdraw(combined) +
  draw_plot_label(
    label    = c("(a)", "(b)"),
    x        = c(0.098, 0.418),   # x position: left edge of each panel
    y        = c(0.84, 0.84),     # y position: top of figure
    size     = 14,
    fontface = "bold"
  )


ggsave("figures/growthModelsMain/asiteMap.pdf", combined_labeled, width = 10, height = 6)


}
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Phenology carry-over ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# plot(budset ~ leafout, emp)
# lm <- lmer(budset ~ leafout + (1 | year), data = emp)
# sum <- summary(lm)
# new_data <- data.frame(leafout = x_seq)
# pred <- predict(lm, newdata = new_data, re.form = NA)  
# lines(x_seq, pred, col = "black", lwd = 2)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Ring width summaries ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
rwmin <- aggregate(lengthMM ~ spp, emp, FUN = min)
rwmax <- aggregate(lengthMM ~ spp, emp, FUN = max)
rwsum <- merge(rwmin, rwmax, by = "spp")
colnames(rwsum) <- c("spp","min", "max")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Z-scored output ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### GDD posterior recovery #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
fitgdd <- readRDS("output/stanOutput/fitGrowthGDDZscored")

df_fitgdd <- as.data.frame(fitgdd)

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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### GSL posterior recovery #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
fitgsl <- readRDS("output/stanOutput/fitGrowthGSLZscored")

df_fitgsl <- as.data.frame(fitgsl)

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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### SOS posterior recovery #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
fitsos <- readRDS("output/stanOutput/fitGrowthSOSZscored")

df_fitsos <- as.data.frame(fitsos)


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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### EOS posterior recovery #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
fiteos <- readRDS("output/stanOutput/fitGrowthEOSZscored")

df_fiteos <- as.data.frame(fiteos)

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

# Add predictors and bind
bspp_df2_z$pred <- "GDD"
bspp_df2_z_gsl$pred <- "GSL"
bspp_df2_z_sos$pred <- "SOS"
bspp_df2_z_eos$pred <- "EOS"
bspp_z_binded <- rbind(bspp_df2_z, bspp_df2_z_gsl, bspp_df2_z_sos, bspp_df2_z_eos)

ainc_z <- subset(bspp_z_binded, spp_name %in% "Alnus incana")
ball_z <- subset(bspp_z_binded, spp_name %in% "Betula alleghaniensis")
bpap_z <- subset(bspp_z_binded, spp_name %in% "Betula papyrifera")
bpop_z <- subset(bspp_z_binded, spp_name %in% "Betula populifolia")

bspp_z_binded$fit_bspp_abs <- abs(bspp_z_binded$fit_bspp)

agg_z <- aggregate(fit_bspp_abs ~ spp_name, bspp_z_binded, function(f) abs(max(f)))

# Aggregate to get only the max effect size for each species

max_ES <- merge(agg_z, bspp_z_binded[, c("fit_bspp", "fit_bspp_per5", "fit_bspp_per95", "pred", "fit_bspp_abs", "spp_name")], 
                   by = c("spp_name", "fit_bspp_abs"))

max_ES$fit_bspp_abs <- NULL

max_ES <- max_ES[order(max_ES$pred),]
