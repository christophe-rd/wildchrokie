# Wildchrokie model
# CRD 14 December 2025

# Goal: Plot model output because modelGrowthGDD is becoming too long and messy

# housekeeping
# rm(list=ls())  
# options(stringsAsFactors = FALSE)
# options(max.print = 150) 
# options(mc.cores = parallel::detectCores())
# options(digits = 3)
cat("SCRIPT STARTED\n")
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

# flags
# makeplots <- FALSE
# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
emp <- read.csv("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/empiricalDataMAIN.csv")
emp <- emp[!is.na(emp$pgsGDD),]
# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# transform data in vectors
y <- emp$lengthCM*10 # ring width in mm
N <- nrow(emp)
gdd <- emp$pgsGDD/200
Nspp <- length(unique(emp$spp_num))
Nsite <- length(unique(emp$site_num))
site <- as.numeric(as.character(emp$site_num))
species <- as.numeric(as.character(emp$spp_num))
treeid <- as.numeric(emp$treeid_num)
Ntreeid <- length(unique(treeid))


# === === === === === === === === === === === === #
#### Recover parameters from the posterior ####
# === === === === === === === === === === === === #
fit <- readRDS("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/output/stanOutput/fitGrowthGDD")
df_fit <- as.data.frame(fit)
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Recover sigmas #####
# unique(colnames(df_fit))
# sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]
# 
# sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]
# 
# sigma_df2 <- data.frame(
#   sigma = character(ncol(sigma_df)),
#   mean = numeric(ncol(sigma_df)),  
#   per5 = NA, 
#   per25 = NA,
#   per75 = NA,
#   per95 = NA
# )
# sigma_df2
# 
# for (i in 1:ncol(sigma_df)) { # i = 1
#   sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
#   sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
#   sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.05), 3)
#   sigma_df2$per25[i] <- round(quantile(sigma_df[[i]], probs = 0.25), 3)
#   sigma_df2$per75[i] <- round(quantile(sigma_df[[i]], probs = 0.75), 3)
#   sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.95), 3)
# }
# sigma_df2
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Recover b spp #####
# bspp_cols <- colnames(df_fit)[grepl("bsp", colnames(df_fit))]
# # remove sigma_bspp for now
# # bspp_cols <- bspp_cols[2:length(bspp_cols)]
# 
# bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# # change their names
# colnames(bspp_df) <- sub("bsp\\[(\\d+)\\]", "\\1", colnames(bspp_df))
# #empty spp df
# bspp_df2 <- data.frame(
#   spp = character(ncol(bspp_df)),
#   fit_bspp = numeric(ncol(bspp_df)),  
#   fit_bspp_per5 = NA, 
#   fit_bspp_per25 = NA,
#   fit_bspp_per75 = NA,
#   fit_bspp_per95 = NA
# )
# for (i in 1:ncol(bspp_df)) { # i = 1
#   bspp_df2$spp[i] <- colnames(bspp_df)[i]         
#   bspp_df2$fit_bspp[i] <- round(mean(bspp_df[[i]]),3)  
#   bspp_df2$fit_bspp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
#   bspp_df2$fit_bspp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
#   bspp_df2$fit_bspp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
#   bspp_df2$fit_bspp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
# }
# bspp_df2
# 
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Recover treeid #####
# 
# # grab treeid 
# treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
# treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]
# treeid_cols <- treeid_cols[!grepl("sigma", treeid_cols)]
# 
# treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# 
# # change their names
# colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# # empty treeid dataframe
# treeid_df2 <- data.frame(
#   treeid = character(ncol(treeid_df)),
#   fit_atreeid = numeric(ncol(treeid_df)),  
#   fit_atreeid_per5 = NA, 
#   fit_atreeid_per25 = NA,
#   fit_atreeid_per75 = NA,
#   fit_atreeid_per95 = NA
# )
# for (i in 1:ncol(treeid_df)) { # i = 1
#   treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
#   treeid_df2$fit_atreeid[i] <- round(mean(treeid_df[[i]]),3)  
#   treeid_df2$fit_atreeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
#   treeid_df2$fit_atreeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
#   treeid_df2$fit_atreeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
#   treeid_df2$fit_atreeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
# }
# treeid_df2
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Recover a spp  #####
# aspp_cols <- colnames(df_fit)[grepl("aspp", colnames(df_fit))]
# 
# aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# # change their names
# colnames(aspp_df) <- sub("aspp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
# #empty aspp df
# aspp_df2 <- data.frame(
#   spp = character(ncol(aspp_df)),
#   fit_aspp = numeric(ncol(aspp_df)),  
#   fit_aspp_per5 = NA, 
#   fit_aspp_per25 = NA,
#   fit_aspp_per75 = NA,
#   fit_aspp_per95 = NA
# )
# for (i in 1:ncol(aspp_df)) { # i = 1
#   aspp_df2$spp[i] <- colnames(aspp_df)[i]         
#   aspp_df2$fit_aspp[i] <- round(mean(aspp_df[[i]]),3)  
#   aspp_df2$fit_aspp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
#   aspp_df2$fit_aspp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
#   aspp_df2$fit_aspp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
#   aspp_df2$fit_aspp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
# }
# aspp_df2
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Recover a site #####
# site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]
# 
# site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# # change their names
# colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# # empty site df
# site_df2 <- data.frame(
#   site = character(ncol(site_df)),
#   fit_a_site = numeric(ncol(site_df)),  
#   fit_a_site_per5 = NA, 
#   fit_a_site_per25 = NA,
#   fit_a_site_per75 = NA,
#   fit_a_site_per95 = NA
# )
# for (i in 1:ncol(site_df)) { # i = 1
#   site_df2$site[i] <- colnames(site_df)[i]         
#   site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
#   site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
#   site_df2$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
#   site_df2$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
#   site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
# }
# site_df2
# 
# # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# # Plot lines with quantiles ####
# # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# 
#  
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Per treeid #####
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# subyvec <- vector()
# for (i in 1:length(unique(emp$treeid_num))) {
#   subyvec[i] <- paste("atreeid", "[",i,"]", sep = "")  
# }
# subyvec
# 
# atreeidsub <- subset(df_fit, select = subyvec)
# 
# colnames(atreeidsub) <- 1:length(subyvec)
# 
# # start by filling a df with treeid intercepts only
# 
# # get the spp and site identities for each tree id
# treeid_spp_site <- unique(emp[, c("treeid_num", "spp_num", "site_num",
#                                   "treeid", "spp", "site", "latbi")])
# 
# # the spp values for each tree id
# treeid_aspp <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
# colnames(treeid_aspp) <- colnames(atreeidsub)
# 
# for (i in seq_len(ncol(treeid_aspp))) { # i = 1
#   tree_id <- as.integer(colnames(treeid_aspp)[i])
#   spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
#   treeid_aspp[, i] <- aspp_df[, spp_id]
# }
# treeid_aspp
# 
# # the site values for each tree id
# treeid_asite <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
# colnames(treeid_asite) <- colnames(atreeidsub)
# 
# for (i in seq_len(ncol(treeid_asite))) { # i = 1
#   tree_id <- as.integer(colnames(treeid_asite)[i])
#   site_id <- treeid_spp_site$site_num[match(tree_id, treeid_spp_site$treeid_num)]
#   treeid_asite[, i] <- site_df[, site_id]
# }
# treeid_asite
# 
# # recover a
# treeid_a <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
# colnames(treeid_a) <- colnames(atreeidsub)
# 
# for (i in seq_len(ncol(treeid_a))) { # i = 1
#   treeid_a[, i] <- df_fit[, "a"]
# }
# 
# # sum all 3 dfs together to get the full intercept for each treeid
# fullintercept <-
#   treeid_a + 
#   atreeidsub +
#   treeid_aspp +
#   treeid_asite
# fullintercept
# 
# # now get the slope for each treeid
# treeid_bspp <- data.frame(matrix(ncol = ncol(atreeidsub), nrow = nrow(df_fit)))
# colnames(treeid_bspp) <- colnames(atreeidsub)
# 
# # back convert the slopes to their original scales
# bspp_df4 <- bspp_df
# for (i in 1:ncol(bspp_df4)){
#   bspp_df4[[i]] <- bspp_df4[[i]] / 200
# }
# 
# for (i in seq_len(ncol(treeid_bspp))) { # i = 30
#   tree_id <- as.integer(colnames(treeid_bspp)[i])
#   spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
#   treeid_bspp[, i] <- bspp_df4[, spp_id]
# }
# treeid_bspp
# 
# treeidvecnum <- 1:ncol(fullintercept)
# treeidvecname <- treeid_spp_site$treeid
# x <- seq(min(emp$pgsGDD), max(emp$pgsGDD), length.out = 100)  
# y_post_list <- list()  # store posterior predictions in a list where each tree id gets matrixad
# 
# # below I create a list where each row is the posterior estimate for each value of gdd (so the first row correspond to the model estimate for the first gdd value stored in x) and each column is the iteration (from 1 to 8000)
# for (i in seq_along(treeidvecnum)) { # i = 1
#   tree_col <- as.character(treeidvecnum[i]) 
#   # TO CHANGE: get the 8000 samples back
#   y_post <- sapply(1:nrow(df_fit), function(f) {
#     rnorm(length(x), 
#           fullintercept[f, tree_col] + treeid_bspp[f, tree_col] * x,
#           sigma_df$sigma_y[f])
#   })
#   y_post_list[[tree_col]] <- y_post
# }
# 
# # PDF output
# sppcols <- c(wes_palette("AsteroidCity1"))[1:4]
# 
# pdf(file = "figures/empiricalData/growthModelSlopesperTreeid.pdf", width = 10, height = 8)
# 
# # Layout: 2 rows × 2 columns per page
# par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# 
# # Loop over trees again to plot each tree individually
# for (i in seq_along(treeidvecnum)) { # i = 1
#   tree_col <- as.character(treeidvecnum[i])
#   tree_col_name <- as.character(treeidvecname[i])
#   y_post <- y_post_list[[tree_col]]
#   
#   # color line by spp
#   tree_id_num <- as.integer(tree_col)
#   
#   # index the dots per treeid
#   emp_treeid <- emp[emp$treeid_num == tree_id_num, ]
#   
#   # calculate mean and 50% credible interval (25%-75%)
#   y_mean <- apply(y_post, 1, mean)
#   y_low  <- apply(y_post, 1, quantile, 0.25)
#   y_high <- apply(y_post, 1, quantile, 0.75)
#   
#   # empty plot first
#   plot(emp$pgsGDD, y, type = "n", 
#        ylim = range(c(emp_treeid$lengthCM * 10, y_low, y_high), na.rm = TRUE),
#        xlab = "Primary growing season GDD", ylab = "Ring width (mm)",
#        main = tree_col_name) # set the name for each plot
#   
#   spp_id <- treeid_spp_site$spp_num[
#     match(tree_id_num, treeid_spp_site$treeid_num)
#   ]
#   line_col <- sppcols[spp_id]
#   
#   # shaded interval
#   polygon(c(x, rev(x)), 
#           c(y_low, # lower interval
#             rev(y_high)), # high interval
#           col = adjustcolor(line_col, alpha.f = 0.3), 
#           border = NA)
#   
#   # mean line
#   lines(x, y_mean,
#         col = line_col,
#         lwd = 2)
#   
#   points(
#     emp_treeid$pgsGDD,
#     emp_treeid$lengthCM * 10,
#     pch = 16,
#     cex = 2,
#     col = line_col)
# }
# dev.off()
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### Per Spp, facet #####
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# # get a vector for each treeid for each species
# spp1vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 1]
# spp2vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 2]
# spp3vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 3]
# spp4vec <- treeid_spp_site$treeid_num[treeid_spp_site$spp_num == 4]
# 
# spp_list <- list(
#   "1" = spp1vec,
#   "2" = spp2vec,
#   "3" = spp3vec,
#   "4" = spp4vec
# )
# 
# mean_post_list <- list()
# for (i in seq_along(treeidvecnum)) {
#   tree_col <- as.character(treeidvecnum[i])
#   mean_post_list[[tree_col]] <- sapply(1:nrow(df_fit), function(f) {
#     fullintercept[f, tree_col] + treeid_bspp[f, tree_col] * x
#     # no sigma_y yet
#   })
# }
# 
# # average the mean predictions across trees within species
# spp_mean_list <- lapply(spp_list, function(tree_vec) {
#   Reduce("+", mean_post_list[as.character(tree_vec)]) / length(tree_vec)
# })
# 
# # re-simulate sigma on the averaged mean
# spp_post_list <- lapply(spp_mean_list, function(mean_mat) {
#   sapply(1:nrow(df_fit), function(f) {
#     rnorm(length(x), mean_mat[, f], sigma_df$sigma_y[f])
#   })
# })
# 
# sppvecnum <- 1:4
# sppvecname <- unique(treeid_spp_site$latbi)
# 
# x <- seq(min(emp$pgsGDD), max(emp$pgsGDD), length.out = 100)   
# sppcols <- c(wes_palette("AsteroidCity1"))[1:4]
# 
# # jpeg output
# jpeg(
#   filename = "figures/empiricalData/growthModelSlopesperSppFacet.jpeg",
#   width = 2400,      # wider image (pixels) → more horizontal room
#   height = 2400,
#   res = 300          # good print-quality resolution
# )
# # Layout: 2 rows × 2 columns per page
# par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# 
# # Loop over trees again to plot each tree individually
# for (i in seq_along(sppvecnum)) { # i = 1
#   
#   spp_column <- as.character(sppvecnum[i])
#   spp_column_name <- as.character(sppvecname[i])
#   
#   # define spp num
#   spp_num <- as.integer(spp_column)
#   
#   # subset empirical data correctly
#   emp_spp <- emp[emp$spp_num == spp_num, ]
#   
#   spp_column <- as.character(sppvecnum[i]) 
#   y_post <- spp_post_list[[spp_column]]
#   
#   # summaries
#   y_mean <- apply(y_post, 1, mean)
#   y_low  <- apply(y_post, 1, quantile, 0.25)
#   y_high <- apply(y_post, 1, quantile, 0.75)
#   
#   # species-specific ylim
#   ylim_spp <- range(c(emp_spp$lengthCM * 10, y_low, y_high), na.rm = TRUE)
#   
#   plot(emp_spp$pgsGDD, emp_spp$lengthCM * 10,
#        type = "n",
#        ylim = ylim_spp,
#        xlab = "Primary growing season GDD",
#        ylab = "Ring width (mm)",
#        main = spp_column_name)
#   
#   # color
#   line_col <- sppcols[spp_num]
#   
#   polygon(
#     c(x, rev(x)),
#     c(y_low, rev(y_high)),
#     col = adjustcolor(line_col, alpha.f = 0.3),
#     border = NA
#   )
#   
#   lines(x, y_mean, col = line_col, lwd = 2)
#   
#   points(
#     emp_spp$pgsGDD,
#     emp_spp$lengthCM * 10,
#     pch = 16,
#     cex = 1,
#     col = line_col
#   )
# }
# 
# dev.off()
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ##### Per spp, non facetted #####
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# # PDF output
# jpeg(
#   filename = "figures/empiricalData/growthModelSlopesperSpp.jpeg",
#   width = 2400,      # wider image (pixels) → more horizontal room
#   height = 2400,
#   res = 300          # good print-quality resolution
# )
# 
# # Layout: 2 rows × 2 columns per page
# par(mar = c(4, 4, 2, 1))
# 
# # below I create a list where each row is the posterior estimate for each value of gdd (so the first row correspond to the model estimate for the first gdd value stored in x) and each column is the iteration (from 1 to 8000)
# 
# plot(emp$pgsGDD, y, type = "n", 
#      ylim = range(min(emp$lengthCM*10), max(emp$lengthCM*10)), 
#      xlab = "Primary growing season GDD", ylab = "Ring width (mm)",
#      main = "species growth responses")
# 
# # Loop over trees again to plot each tree individually
# for (i in seq_along(sppvecnum)) { # i = 1
#   spp_column <- as.character(sppvecnum[i])
#   spp_column_name <- as.character(sppvecname[i])
#   y_post <- spp_post_list[[spp_column]]
#   
#   # color line by spp
#   spp_num <- as.integer(spp_column)
#   
#   # calculate mean and 50% credible interval (25%-75%)
#   y_mean <- apply(y_post, 1, mean)
#   y_low  <- apply(y_post, 1, quantile, 0.25)
#   y_high <- apply(y_post, 1, quantile, 0.75)
#   
#   line_col <- sppcols[spp_num]
#   
#   # shaded interval
#   polygon(c(x, rev(x)), 
#           c(y_low, # lower interval
#             rev(y_high)), # high interval
#           col = adjustcolor(line_col, alpha.f = 0.2), 
#           border = NA)
#   
#   # mean line
#   lines(x, y_mean,
#         col = line_col,
#         lwd = 2)
#   
#   emp_spp <- emp[emp$spp_num == spp_num, ]
# 
#   points(
#     emp_spp$pgsGDD,
#     emp_spp$lengthCM * 10,
#     pch = 16,
#     cex = 1,
#     col = line_col)
#   
#   legend(
#     "topleft",
#     legend = sppvecname,
#     col = sppcols[as.integer(sppvecnum)],
#     lwd = 2,
#     pch = 16,
#     bty = "n",
#     cex = 1.5
#   )
#   
# }
# dev.off()
# 
# 
# # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# # Mu plots #####
# # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# n_spp <- nrow(bspp_df2)
# n_site <- nrow(site_df2)
# y_pos <- 1:n_spp 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ###### asp ######
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# aspp_df2$spp_name <- emp$latbi[match(aspp_df2$spp, emp$spp_num)]
# 
# jpeg("figures/empiricalData/aspp_mean_plot.jpeg", width = 8, height = 6, units = "in", res = 300)
# 
# par(mar = c(5, 10, 2, 2)) 
# 
# plot(aspp_df2$fit_aspp, y_pos,
#      xlim = range(c(aspp_df2$fit_aspp_per5, aspp_df2$fit_aspp_per95)),
#      ylim = c(0.5, n_spp + 0.5),
#      xlab = "Ring width intercept values (mm)",
#      ylab = "",
#      yaxt = "n",
#      pch = 16,
#      cex = 2,
#      col = sppcols,
#      frame.plot = FALSE)
# 
# # color labels
# for (i in seq_along(y_pos)) {
#   axis(2, at = y_pos[i],
#        labels = aspp_df2$spp_name[i],
#        las = 2,
#        col.axis = sppcols[i],
#        tick = FALSE,
#        cex.axis = 1)
# }
# 
# # error bars and dashed line
# segments(aspp_df2$fit_aspp_per5,  y_pos,
#          aspp_df2$fit_aspp_per95, y_pos,
#          col = sppcols, lwd = 1.5)
# segments(aspp_df2$fit_aspp_per25, y_pos,
#          aspp_df2$fit_aspp_per75, y_pos,
#          col = sppcols, lwd = 3)
# 
# abline(v = 0, lty = 2, col = "black")
# 
# dev.off()
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ###### asite ######
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# site_df2$site_name <- emp$site[match(site_df2$site, emp$site_num)]
# sitecolors <- c(wes_palette("Darjeeling1"))[1:4]
# 
# jpeg("figures/empiricalData/asite_mean_plot.jpeg", width = 8, height = 6, units = "in", res = 300)
# 
# par(mar = c(5, 10, 2, 2)) 
# 
# plot(site_df2$fit_a_site, y_pos,
#      xlim = range(c(site_df2$fit_a_site_per5, site_df2$fit_a_site_per95)),
#      ylim = c(0.5, n_site + 0.5),
#      xlab = "Ring width intercept values (mm)",
#      ylab = "",
#      yaxt = "n",
#      pch = 16,
#      cex = 2,
#      col = sitecolors,
#      frame.plot = FALSE)
# 
# # color labels
# for (i in seq_along(y_pos)) {
#   axis(2, at = y_pos[i],
#        labels = site_df2$site_name[i],
#        las = 2,
#        col.axis = sitecolors[i],
#        tick = FALSE,
#        cex.axis = 1)
# }
# 
# # error bars and dashed line
# segments(site_df2$fit_a_site_per5,  y_pos,
#          site_df2$fit_a_site_per95, y_pos,
#          col = sitecolors, lwd = 1.5)
# segments(site_df2$fit_a_site_per25, y_pos,
#          site_df2$fit_a_site_per75, y_pos,
#          col = sitecolors, lwd = 3)
# 
# abline(v = 0, lty = 2, col = "black")
# 
# dev.off()
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ###### asp ######
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# bspp_df2$spp_name <- emp$latbi[match(bspp_df2$spp, emp$spp_num)]
# 
# jpeg("figures/empiricalData/bspp_mean_plot.jpeg", width = 8, height = 6, units = "in", res = 300)
# 
# par(mar = c(5, 10, 2, 2)) 
# 
# plot(bspp_df2$fit_bspp, y_pos,
#      xlim = range(c(bspp_df2$fit_bspp_per5, bspp_df2$fit_bspp_per95)),
#      ylim = c(0.5, n_spp + 0.5),
#      xlab = "Ring width (mm) change/200 GDD",
#      ylab = "",
#      yaxt = "n",
#      pch = 16,
#      cex = 2,
#      col = sppcols,
#      frame.plot = FALSE)
# 
# # color labels
# for (i in seq_along(y_pos)) {
#   axis(2, at = y_pos[i],
#        labels = bspp_df2$spp_name[i],
#        las = 2,
#        col.axis = sppcols[i],
#        tick = FALSE,
#        cex.axis = 1)
# }
# 
# # error bars and dashed line
# segments(bspp_df2$fit_bspp_per5,  y_pos,
#          bspp_df2$fit_bspp_per95, y_pos,
#          col = sppcols, lwd = 1.5)
# segments(bspp_df2$fit_bspp_per25, y_pos,
#          bspp_df2$fit_bspp_per75, y_pos,
#          col = sppcols, lwd = 3)
# 
# abline(v = 0, lty = 2, col = "black")
# 
# dev.off()
# 
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ##### full treeid mu plots #####
# # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# # Mean plots with atreeid ####
# treeid_df2$treeid <- as.numeric(treeid_df2$treeid)
# treeid_df2$treeid_name <- emp$treeid[match(treeid_df2$treeid, emp$treeid_num)]
# 
# # now do the same, but for species
# treeid_df2$spp <- emp$spp[match(treeid_df2$treeid, emp$treeid_num)]
# 
# # same for site
# treeid_df2$site <- emp$site[match(treeid_df2$treeid, emp$treeid_num)]
# 
# # quick check that I didn't mess anything up
# un <- emp[!duplicated(emp$treeid),]
# table(un$spp)
# table(treeid_df2$spp)
# 
# sub <- subset(emp, select = c("treeid_num", "spp_num", "site_num"))
# sub <- sub[!duplicated(sub$treeid_num),]
# 
# # recalculate the full intercept without the grand mean
# fullintercept2 <-
#   atreeidsub +
#   treeid_aspp +
#   treeid_asite
# fullintercept2
# 
# # get posterior means and quantiles
# 
# # empty treeid dataframe
# treeid_df4 <- data.frame(
#   treeid = character(ncol(fullintercept2)),
#   fit_atreeid = numeric(ncol(fullintercept2)),  
#   fit_atreeid_per5 = NA, 
#   fit_atreeid_per25 = NA,
#   fit_atreeid_per75 = NA,
#   fit_atreeid_per95 = NA
# )
# for (i in 1:ncol(fullintercept2)) { # i = 1
#   treeid_df4$treeid[i] <- colnames(fullintercept2)[i]         
#   treeid_df4$fit_atreeid[i] <- round(mean(fullintercept2[[i]]),3)  
#   treeid_df4$fit_atreeid_per5[i] <- round(quantile(fullintercept2[[i]], probs = 0.05), 3)
#   treeid_df4$fit_atreeid_per25[i] <- round(quantile(fullintercept2[[i]], probs = 0.25), 3)
#   treeid_df4$fit_atreeid_per75[i] <- round(quantile(fullintercept2[[i]], probs = 0.75), 3)
#   treeid_df4$fit_atreeid_per95[i] <- round(quantile(fullintercept2[[i]], probs = 0.95), 3)
# }
# treeid_df4
# 
# # get the og treeid names, spp and site back:
# treeid_df4$treeid <- as.numeric(treeid_df4$treeid)
# treeid_df4$treeid_name <- emp$treeid[match(treeid_df4$treeid,
#                                                     emp$treeid_num)]
# treeid_df4$spp_name <- emp$latbi[match(treeid_df4$treeid,
#                                               emp$treeid_num)]
# treeid_df4$spp_num <- emp$spp_num[match(treeid_df4$treeid,
#                                        emp$treeid_num)]
# treeid_df4$site_name <- emp$site[match(treeid_df4$treeid,
#                                                 emp$treeid_num)]
# treeid_df4$site_num <- emp$site_num[match(treeid_df4$treeid,
#                                        emp$treeid_num)]
# 
# # Prep for the figure
# 
# # define a gap between species clusters
# gap <- 3
# 
# # y positions
# # treeid_df4$y_pos <- NA
# current_y <- 1
# 
# species_order <- c(
#   "Alnus incana", 
#   "Betula alleghaniensis", 
#   "Betula papyrifera", 
#   "Betula populifolia")
# # site_order <- c("SH", "GR", "WM", "HF")
# site_order <- c(
#   "HF",
#   "WM",
#   "GR", 
#   "SH")
# 
# # col
# my_colors <- c(
#   "Alnus incana" = wes_palette("AsteroidCity1")[1],
#   "Betula alleghaniensis" = wes_palette("AsteroidCity1")[2],
#   "Betula papyrifera" = wes_palette("AsteroidCity1")[3],
#   "Betula populifolia" = wes_palette("AsteroidCity1")[4]
# )
# # shapes for sites
# my_shapes <- c(
#   HF = 19,
#   WM = 18,
#   GR = 15,
#   SH = 17
# )
# 
# 
# # open device
# pdf(
#   file = "figures/empiricalData/meanPlotGrowthGDD_treeidBYspp.pdf",
#   width = 8,  
#   height = 8
# )
# par(mar = c(
#   4, # lower margin 
#   6, 
#   4, # upper margin  
#   6)) # right margin
# 
# treeid_df4$spp_name  <- factor(treeid_df4$spp_name, levels = species_order)
# treeid_df4$site_num <- factor(treeid_df4$site_num, levels = site_order)
# 
# treeid_df4 <- treeid_df4[
#   order(treeid_df4$spp_name, treeid_df4$site_num, treeid_df4$treeid),
# ]
# 
# treeid_df4$y_pos <- seq_len(nrow(treeid_df4))
# 
# total_rows <- nrow(treeid_df4) + (length(species_order) - 1) * gap
# current_y <- total_rows 
# 
# for(sp in species_order){ # sp = "Alnus incana"
#   idx <- which(treeid_df4$spp_name == sp)
#   n <- length(idx)
#   
#   # assign sequential positions for this species
#   treeid_df4$y_pos[idx] <- current_y:(current_y - n + 1)
#   
#   # move cursor down with a gap before next species cluster
#   current_y <- current_y - n - gap
# }
# 
# # Set up empty plot
# plot(
#   NA, NA,
#   xlim = range(c(treeid_df4$fit_atreeid_per5-2,
#                  treeid_df4$fit_atreeid_per95+2)),
#   ylim = c(0.5, max(treeid_df4$y_pos) + 0.5),
#   xlab = "treeid intercept values",
#   ylab = "",
#   yaxt = "n",
#   bty = "l"
# )
# 
# 
# # --- Add horizontal error bars (5–95%) ---
# segments(
#   x0 = treeid_df4$fit_atreeid_per5,
#   x1 = treeid_df4$fit_atreeid_per95,
#   y0 = treeid_df4$y_pos,
#   col = adjustcolor(my_colors[treeid_df4$spp_name], alpha.f = 0.7),
#   lwd = 1
# )
# 
# # --- Add thicker horizontal error bars (25–75%) ---
# segments(
#   x0 = treeid_df4$fit_atreeid_per25,
#   x1 = treeid_df4$fit_atreeid_per75,
#   y0 = treeid_df4$y_pos,
#   col = adjustcolor(my_colors[treeid_df4$spp_name], alpha.f = 0.7),
#   lwd = 1.5
# )
# 
# # --- Add the points ---
# points(
#   treeid_df4$fit_atreeid,
#   treeid_df4$y_pos,
#   cex = 0.8,
#   pch = my_shapes[treeid_df4$site_name],
#   col = adjustcolor(my_colors[treeid_df4$spp_name], alpha.f = 0.7)
# )
# 
# spp_y_top <- tapply(treeid_df4$y_pos, treeid_df4$spp_name, max)
# aspp_df2$y_pos <- spp_y_top[aspp_df2$spp_name] + 1
# 
# segments(
#   x0 = aspp_df2$fit_aspp_per5,
#   x1 = aspp_df2$fit_aspp_per95,
#   y0 = aspp_df2$y_pos,
#   col = adjustcolor(my_colors[aspp_df2$spp_name], alpha.f = 0.9),
#   lwd = 2
# )
# 
# segments(
#   x0 = aspp_df2$fit_aspp_per25,
#   x1 = aspp_df2$fit_aspp_per75,
#   y0 = aspp_df2$y_pos,
#   col = my_colors[aspp_df2$spp_name],
#   lwd = 3
# )
# points(
#   aspp_df2$fit_aspp,
#   aspp_df2$y_pos,
#   pch = 16,
#   bg  = my_colors[aspp_df2$spp_name],
#   col = my_colors[aspp_df2$spp_name],
#   cex = 1.5
# )
# 
# # --- Add vertical line at 0 ---
# abline(v = 0, lty = 2)
# 
# # --- Add custom y-axis labels (reverse order if needed) ---
# axis(
#   side = 2,
#   at = treeid_df4$y_pos,
#   labels = treeid_df4$treeid_name,
#   cex.axis = 0.5,
#   las = 1
# )
# # spp_name mean
# spp_y <- tapply(treeid_df4$y_pos, treeid_df4$spp_name, mean)
# site_y <- tapply(treeid_df4$y_pos, treeid_df4$site_num, max)
# 
# ## order species by mean y descending (top of plot first)
# species_legend_order <- names(sort(spp_y, decreasing = TRUE))
# site_legend_order <- names(sort(site_y, decreasing = FALSE))
# 
# ## species legend (colors matched by name)
# legend(
#   x = max(treeid_df4$fit_atreeid_per95) - 5,
#   y = max(treeid_df4$y_pos) + 1,
#   legend = species_legend_order,
#   col = my_colors[species_legend_order],    # index so colors match
#   pch = 16,
#   pt.cex = 1.2,
#   title = "Species",
#   bty = "n"
# )
# 
# site_legend_order <- c("SH", "GR", "WM", "HF")
# # site_num legen
# legend(
#   x = max(treeid_df4$fit_atreeid_per95) - 2,
#   y = max(treeid_df4$y_pos) - 15,
#   legend = site_legend_order,
#   pch = my_shapes[site_legend_order],
#   pt.cex = 1.2,
#   title = "Sites",
#   bty = "n"
# )
# dev.off()

