# Wildchrokie model -- leaf out prediction
# CRD 17 December 2025
# Related to modelGDDLeafout script where it estimates the gdd at leafout. 
# This is an additional code to plot those results 

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

# directories
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


# RUN MODEL ON EMPIRICAL DATA ####
# empirical data
emp <- read.csv("output/empiricalDataNORing.csv")
gdd <- read.csv("output/gddByYear.csv")
fit <- readRDS("output/stanOutput/gddLeafout_empData_fit") 

no_naleafout <- emp[!is.na(emp$leafoutGDD),]

# give numeric ids to my groups 
no_naleafout$site_num <- match(no_naleafout$site, unique(no_naleafout$site))
no_naleafout$spp_num <- match(no_naleafout$spp, unique(no_naleafout$spp))
no_naleafout$treeid_num <- match(no_naleafout$treeid, unique(no_naleafout$treeid))

scale <- 20

# Recover parameters ####
# === === === === === === === === === === === === === === === === === === === ==
df_fit <- as.data.frame(fit)

# get values to original scale
for (i in 1:ncol(df_fit)){
  df_fit[[i]] <- df_fit[[i]] * scale
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
###### Recover sigmas ######
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
###### Recover treeid ######
# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("atreeid", colnames(df_fit))]
treeid_cols <- treeid_cols[!grepl("zatreeid", treeid_cols)]
# remove sigma_asp for now
treeid_cols <- treeid_cols[2:length(treeid_cols)]

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
###### Recover a spp  ######
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
###### Recover a site ######
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

df_fit <- as.data.frame(fit)



# Plots ####
aspp_df2$spp <- as.numeric(aspp_df2$spp)
aspp_df2$spp_name <- no_naleafout$spp[match(aspp_df2$spp, no_naleafout$spp_num)]

##### aspp intercepts means #####
ggplot(aspp_df2, aes(x = fit_a_spp, y = spp_name, color = spp_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_spp_per5, 
                     xmax = fit_a_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_spp_per25, 
                     xmax = fit_a_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_bw() +
  scale_y_discrete(limits = rev)   
ggsave("figures/gddLeafout_simData/sppEstimates.jpeg", width = 8, height = 6, units = "in", dpi = 300)

## site intercept means ####

site_df2$site <- as.numeric(site_df2$site)
site_df2$site_name <- no_naleafout$site[match(site_df2$site, no_naleafout$site_num)]

ggplot(site_df2, aes(x = fit_a_site, y = site_name, color = site_name)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_site_per5, 
                     xmax = fit_a_site_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_site_per25, 
                     xmax = fit_a_site_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_bw() +
  scale_y_discrete(limits = rev)   
ggsave("figures/gddLeafout_simData/sppEstimates.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# Mean plots with atreeid ####
no_naleafout

treeid_df2$treeid <- as.numeric(treeid_df2$treeid)
treeid_df2$treeid_name <- no_naleafout$treeid[match(treeid_df2$treeid, no_naleafout$treeid_num)]

# now do the same, but for species
treeid_df2$spp <- no_naleafout$spp[match(treeid_df2$treeid, no_naleafout$treeid_num)]

# same for site
treeid_df2$site <- no_naleafout$site[match(treeid_df2$treeid, no_naleafout$treeid_num)]


# quick check that I didn't mess anything up
un <- no_naleafout[!duplicated(no_naleafout$treeid),]
table(un$spp)
table(treeid_df2$spp)

sub <- subset(no_naleafout, select = c("treeid_num", "spp_num", "site_num"))
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
treeid_df4$treeid_name <- no_naleafout$treeid[match(treeid_df4$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4$spp_name <- no_naleafout$spp[match(treeid_df4$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4$site_name <- no_naleafout$site[match(treeid_df4$treeid,
                                                    no_naleafout$treeid_num)]
treeid_df4

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# open device
jpeg(
  filename = "figures/gddLeafout_empData/meanPlot_treeidBYspp.jpeg",
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
  xlim = range(c(treeid_df4$fit_a_treeid_per5-15,
                 treeid_df4$fit_a_treeid_per95+15)),
  ylim = c(1, max(treeid_df4$y_pos) + 0.5),
  xlab = "Estimate of GDD change from the grand mean",
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
  x = max(treeid_df4$fit_a_treeid_per95) - 6,
  y = max(treeid_df4$y_pos) - 2,
  legend = species_legend_order,
  col = my_colors[species_legend_order],    # index so colors match
  pch = 16,
  pt.cex = 1.2,
  title = "Species",
  bty = "n"
)

site_legend_order <- c("SH", "GR", "WM", "HF")
# site legend
legend(
  x = max(treeid_df4$fit_a_treeid_per95) - 6,
  y = max(treeid_df4$y_pos) - 25,
  legend = site_legend_order,
  pch = my_shapes[site_legend_order],
  pt.cex = 1.2,
  title = "Sites",
  bty = "n"
)

dev.off()


