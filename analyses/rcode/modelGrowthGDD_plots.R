# Wildchrokie model
# CRD 14 December 2025

# Goal: Plot model output because modelGrowthGDD is becoming too long and messy

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(mc.cores = parallel::detectCores())
options(digits = 3)
# quartz()

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(shinystan)
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

runSimData <- FALSE

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
emp <- read.csv("output/empiricalDataMAIN.csv")

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
fit <- readRDS("output/stanOutput/fit")
df_fit <- as.data.frame(fit)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Recover sigmas #####
unique(colnames(df_fit))
sigma_cols <- colnames(df_fit)[grepl("sigma", colnames(df_fit))]

sigma_df <- df_fit[, colnames(df_fit) %in% sigma_cols]

sigma_df2 <- data.frame(
  sigma = character(ncol(sigma_df)),
  mean = numeric(ncol(sigma_df)),  
  per5 = NA, 
  per25 = NA,
  per75 = NA,
  per95 = NA
)
sigma_df2

for (i in 1:ncol(sigma_df)) { # i = 1
  sigma_df2$sigma[i] <- colnames(sigma_df)[i]         
  sigma_df2$mean[i] <- round(mean(sigma_df[[i]]),3)  
  sigma_df2$per5[i] <- round(quantile(sigma_df[[i]], probs = 0.05), 3)
  sigma_df2$per25[i] <- round(quantile(sigma_df[[i]], probs = 0.25), 3)
  sigma_df2$per75[i] <- round(quantile(sigma_df[[i]], probs = 0.75), 3)
  sigma_df2$per95[i] <- round(quantile(sigma_df[[i]], probs = 0.95), 3)
}
sigma_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Recover b spp #####
bspp_cols <- colnames(df_fit)[grepl("bsp", colnames(df_fit))]
# remove sigma_bspp for now
# bspp_cols <- bspp_cols[2:length(bspp_cols)]

bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# change their names
colnames(bspp_df) <- sub("bsp\\[(\\d+)\\]", "\\1", colnames(bspp_df))
#empty spp df
bspp_df2 <- data.frame(
  spp = character(ncol(bspp_df)),
  fit_b_spp = numeric(ncol(bspp_df)),  
  fit_b_spp_per5 = NA, 
  fit_b_spp_per25 = NA,
  fit_b_spp_per75 = NA,
  fit_b_spp_per95 = NA
)
for (i in 1:ncol(bspp_df)) { # i = 1
  bspp_df2$spp[i] <- colnames(bspp_df)[i]         
  bspp_df2$fit_b_spp[i] <- round(mean(bspp_df[[i]]),3)  
  bspp_df2$fit_b_spp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
  bspp_df2$fit_b_spp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
  bspp_df2$fit_b_spp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
  bspp_df2$fit_b_spp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
}
bspp_df2


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Recover treeid #####

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]
treeid_cols <- treeid_cols[!grepl("sigma", treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]

# change their names
colnames(treeid_df) <- sub("atreeid\\[(\\d+)\\]", "\\1", colnames(treeid_df))
# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_a_treeid = numeric(ncol(treeid_df)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_a_treeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_a_treeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Recover a spp  #####
aspp_cols <- colnames(df_fit)[grepl("asp", colnames(df_fit))]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub("asp\\[(\\d+)\\]", "\\1", colnames(aspp_df))
#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_a_spp = numeric(ncol(aspp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per25 = NA,
  fit_a_spp_per75 = NA,
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_a_spp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_a_spp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_a_spp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_a_spp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_a_spp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Recover a site #####
site_cols <- colnames(df_fit)[grepl("asite", colnames(df_fit))]

site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub("asite\\[(\\d+)\\]", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per25 = NA,
  fit_a_site_per75 = NA,
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2




##### Plot lines #####
# Gdd on the x axis and growth on y ####
aspp_df2$a <- mean(df_fit[,"a"])
aspp_df2$a_asp <- aspp_df2$a + aspp_df2$fit_a_spp
aspp_df2$b <- mean(df_fit[,"b"])
aspp_df2$bsp <- bspp_df2$fit_b_spp
aspp_df2$b_bsp <- aspp_df2$b + bspp_df2$fit_b_spp

colnames(aspp_df2)[colnames(aspp_df2) == "spp"] <- "spp_num"

emp2 <- emp
emp2 <- merge(emp2, aspp_df2[, c("spp_num", "a", "b", "b_bsp", "a_asp")], 
              by = "spp_num")
# plot lines
ggplot(emp2) +
  geom_point(aes(x = pgsGDD/200, y = lengthCM*10, colour = spp)) +
  geom_abline(aes(intercept = a_asp, slope = b_bsp, colour = spp), 
              linewidth = 0.5) +
  labs(title = "", x = "pgsGDD", y = "ring width in mm") +
  scale_colour_manual(values = wes_palette("AsteroidCity1")) +
  # facet_wrap(~ spp) +
  theme_minimal()
ggsave("figures/empiricalData/empiricalData/slope_intercepts_varyingslopes.jpeg", 
       width = 8, height = 6, units = "in", dpi = 300)

##### mu plots #####
###### asp ######
aspp_df2$spp <- as.numeric(aspp_df2$spp)
aspp_df2$spp_name <- emp$spp[match(aspp_df2$spp, emp$spp_num)]

asp_mean_plot <- ggplot(aspp_df2, aes(x = fit_a_spp, y = spp_name, color = spp_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_spp_per5, 
                     xmax = fit_a_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_spp_per25, 
                     xmax = fit_a_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values") +
  # facet_wrap(~ model, nrow =2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev)  
asp_mean_plot
ggsave("figures/empiricalData/asp_mean_plot.jpeg",
       asp_mean_plot,
       width = 8, height = 6, 
       units = "in", dpi = 300)

###### asite ######
site_df2$site <- as.numeric(site_df2$site)
site_df2$site_name <- emp$site[match(site_df2$site, emp$site_num)]

site_mean_plot <- ggplot(site_df2, 
                         aes(x = fit_a_site, 
                             y = site_name, 
                             color = site_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_site_per5, 
                     xmax = fit_a_site_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_site_per25, 
                     xmax = fit_a_site_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Site", x = "Site intercept values") +
  # facet_wrap(~ model, nrow =2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev)  
site_mean_plot
ggsave("figures/empiricalData/asite_mean_plot.jpeg",
       site_mean_plot,
       width = 8, height = 6, 
       units = "in", dpi = 300)

###### bsp ###### 
bspp_df2$spp <- as.numeric(bspp_df2$spp)
bspp_df2$spp_name <- emp$spp[match(bspp_df2$spp, emp$spp_num)]

bsp_mean_plot <- ggplot(bspp_df2, 
                        aes(x = fit_b_spp, 
                            y = spp_name, 
                            color = spp_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_b_spp_per5, 
                     xmax = fit_b_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_b_spp_per25, 
                     xmax = fit_b_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  # scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species slope values") +
  # facet_wrap(~ model, nrow =2) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_minimal() +
  scale_y_discrete(limits = rev)  
bsp_mean_plot
ggsave("figures/empiricalData/bsp_mean_plot.jpeg",
       bsp_mean_plot,
       width = 8, height = 6, 
       units = "in", dpi = 300)


combined_plot <- (asp_mean_plot) /
                  (site_mean_plot) / 
                  (bsp_mean_plot)
combined_plot
ggsave("figures/empiricalData/combinedMeanPlots.jpeg", combined_plot, width = 6, height = 8, units = "in", dpi = 300)


# full treeid mu plots ####

# Mean plots with atreeid ####
treeid_df2$treeid <- as.numeric(treeid_df2$treeid)
treeid_df2$treeid_name <- emp$treeid[match(treeid_df2$treeid, emp$treeid_num)]

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
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# sum atreeid, asp, asite
aspp_df3 <- as.matrix(aspp_df)
site_df3 <- as.matrix(site_df)
treeid_df3 <- as.matrix(treeid_df)

fullatreeid <-
  treeid_df3[, sub$treeid_num] +
  aspp_df3[, sub$spp_num] +
  site_df3[, sub$site_num] 
fullatreeid

fullatreeid2 <- as.data.frame(fullatreeid)
# get posterior means and quantiles

# empty treeid dataframe
treeid_df4 <- data.frame(
  treeid = character(ncol(fullatreeid2)),
  fit_a_treeid = numeric(ncol(fullatreeid2)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(fullatreeid2)) { # i = 1
  treeid_df4$treeid[i] <- colnames(fullatreeid2)[i]         
  treeid_df4$fit_a_treeid[i] <- round(mean(fullatreeid2[[i]]),3)  
  treeid_df4$fit_a_treeid_per5[i] <- round(quantile(fullatreeid2[[i]], probs = 0.05), 3)
  treeid_df4$fit_a_treeid_per25[i] <- round(quantile(fullatreeid2[[i]], probs = 0.25), 3)
  treeid_df4$fit_a_treeid_per75[i] <- round(quantile(fullatreeid2[[i]], probs = 0.75), 3)
  treeid_df4$fit_a_treeid_per95[i] <- round(quantile(fullatreeid2[[i]], probs = 0.95), 3)
}
treeid_df4

# get the og treeid names, spp and site back:
treeid_df4$treeid <- as.numeric(treeid_df4$treeid)
treeid_df4$treeid_name <- emp$treeid[match(treeid_df4$treeid,
                                                    emp$treeid_num)]
treeid_df4$spp_name <- emp$spp[match(treeid_df4$treeid,
                                              emp$treeid_num)]
treeid_df4$site_name <- emp$site[match(treeid_df4$treeid,
                                                emp$treeid_num)]
treeid_df4

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# open device
jpeg(
  filename = "figures/empiricalData/meanPlotGrowthGDD_treeidBYspp.jpeg",
  width = 2400,      # wider image (pixels) → more horizontal room
  height = 2400,
  res = 300          # good print-quality resolution
)
par(mar = c(
  4, 
  6, 
  4, 
  5))

# define a gap between species clusters
gap <- 3

# y positions
treeid_df4$y_pos <- NA
current_y <- 1

species_order <- c(
  "ALNINC", 
  "BETALL", 
  "BETPAP", 
  "BETPOP")
# site_order <- c("SH", "GR", "WM", "HF")
site_order <- c(
  "HF",
  "WM",
  "GR", 
  "SH")

# col
my_colors <- c(
  ALNINC = wes_palette("AsteroidCity1")[1],
  BETALL = wes_palette("AsteroidCity1")[2],
  BETPAP = wes_palette("AsteroidCity1")[3],
  BETPOP = wes_palette("AsteroidCity1")[4]
)
# shapes for sites
my_shapes <- c(
  HF = 19,
  WM = 18,
  GR = 15,
  SH = 17
)


treeid_df4$spp  <- factor(treeid_df4$spp, levels = species_order)
treeid_df4$site <- factor(treeid_df4$site, levels = site_order)

treeid_df4 <- treeid_df4[
  order(treeid_df4$spp, treeid_df4$site, treeid_df4$treeid),
]

treeid_df4$y_pos <- seq_len(nrow(treeid_df4))

for(sp in species_order){
  idx <- which(treeid_df4$spp == sp)
  n <- length(idx)
  
  # assign sequential positions for this species
  treeid_df4$y_pos[idx] <- current_y:(current_y + n - 1)
  
  # move cursor down with a gap before next species cluster
  current_y <- current_y + n + gap
}

treeid_df4$y_pos

# Set up empty plot
plot(
  NA, NA,
  xlim = range(c(treeid_df4$fit_a_treeid_per5-0.5,
                 treeid_df4$fit_a_treeid_per95+0.5)),
  ylim = c(0.5, max(treeid_df4$y_pos) + 0.5),
  xlab = "treeid intercept values",
  ylab = "",
  yaxt = "n"  
)


# --- Add horizontal error bars (5–95%) ---
segments(
  x0 = treeid_df4$fit_a_treeid_per5,
  x1 = treeid_df4$fit_a_treeid_per95,
  y0 = treeid_df4$y_pos,
  col = adjustcolor(my_colors[treeid_df4$spp], alpha.f = 0.4),
  lwd = 1
)

# --- Add thicker horizontal error bars (25–75%) ---
segments(
  x0 = treeid_df4$fit_a_treeid_per25,
  x1 = treeid_df4$fit_a_treeid_per75,
  y0 = treeid_df4$y_pos,
  col = adjustcolor(my_colors[treeid_df4$spp], alpha.f = 0.4),
  lwd = 1.5
)

# --- Add the points ---
points(
  treeid_df4$fit_a_treeid,
  treeid_df4$y_pos,
  cex = 0.8,
  pch = my_shapes[treeid_df4$site],
  col = adjustcolor(my_colors[treeid_df4$spp], alpha.f = 0.4)
)

aspp_df2$spp <- aspp_df2$spp_name

spp_y <- tapply(treeid_df4$y_pos, treeid_df4$spp, mean)

aspp_df2$y_pos <- spp_y[aspp_df2$spp]


segments(
  x0 = aspp_df2$fit_a_spp_per5,
  x1 = aspp_df2$fit_a_spp_per95,
  y0 = aspp_df2$y_pos,
  col = adjustcolor(my_colors[aspp_df2$spp], alpha.f = 0.9),
  lwd = 2
)

segments(
  x0 = aspp_df2$fit_a_spp_per25,
  x1 = aspp_df2$fit_a_spp_per75,
  y0 = aspp_df2$y_pos,
  col = my_colors[aspp_df2$spp],
  lwd = 4
)
points(
  aspp_df2$fit_a_spp,
  aspp_df2$y_pos,
  pch = 21,
  bg  = my_colors[aspp_df2$spp],
  col = "black",
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
# spp mean
spp_y <- tapply(treeid_df4$y_pos, treeid_df4$spp, mean)
site_y <- tapply(treeid_df4$y_pos, treeid_df4$site, max)

## order species by mean y descending (top of plot first)
species_legend_order <- names(sort(spp_y, decreasing = TRUE))
site_legend_order <- names(sort(site_y, decreasing = FALSE))

## species legend (colors matched by name)
legend(
  x = max(treeid_df4$fit_a_treeid_per95) - 0.1,
  y = max(treeid_df4$y_pos) + 1,
  legend = species_legend_order,
  col = my_colors[species_legend_order],    # index so colors match
  pch = 16,
  pt.cex = 1.2,
  title = "Species",
  bty = "n"
)

site_legend_order <- c("SH", "GR", "WM", "HF")
# site legen
legend(
  x = max(treeid_df4$fit_a_treeid_per95) - 0.1,
  y = max(treeid_df4$y_pos) - 15,
  legend = site_legend_order,
  pch = my_shapes[site_legend_order],
  pt.cex = 1.2,
  title = "Sites",
  bty = "n"
)

dev.off()

