# Wildchrokie model
# CRD 28 October April 2025
# Start plotting empirical data for wildchrokie

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 3)
# quartz()

# Load library 
library(ggplot2)
library(rstan)
library(rstanarm)
library(future)
library(shinystan)
library(wesanderson)
library(patchwork)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

sim <- read.csv("output/simdata.csv")
emp <- read.csv("output/empiricalDataMAIN.csv")
gdd <- read.csv("output/gddByYear.csv")

# load fit
fit <- readRDS("output/stanOutput/fitEmpirical_stanlmer")

# Plot empirical data
ggplot(emp, aes(x= pgsGDD, y = lengthCM, group = spp, color = spp)) +
  # geom_violin() +
  geom_boxplot() +
  facet_wrap(~spp, ncol = 1, nrow = 4) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) + 
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()
ggsave("figures/empiricalData/sppBoxPlots.jpeg", width = 6, height = 6, units = "in", dpi = 300)

# let's just look at 2018
y18 <- subset(emp, year == "2018")

ggplot(y18, aes(x= pgsGDD, y = lengthCM, group = spp, color = spp)) +
  geom_boxplot(position = "dodge2") +
  facet_wrap(~year) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) + 
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()

ggplot(emp, aes(x = pgsGDD, y = lengthCM, 
                color = spp, 
                fill = spp)) +
  geom_point(size = 2, alpha = 0.7, aes(shape = as.factor(year))) +  # optional: show raw data points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  facet_wrap(~spp)+
  theme_minimal()
ggsave("figures/empiricalData/sppLinearRegressions_nofacet.jpeg", width = 6, height = 6, units = "in", dpi = 300)


# plot gdd
gddstats <- aggregate(GDD_10 ~ doy, FUN = mean, data = gdd)
colnames(gddstats)[2] <-  "mean"
gddstats2 <- aggregate(GDD_10 ~ doy, FUN = min, data = gdd)
colnames(gddstats2)[2] <-  "min"
gddstats3 <- aggregate(GDD_10 ~ doy, FUN = max, data = gdd)
colnames(gddstats3)[2] <-  "max"

gddstats <- merge(gddstats, gddstats2, by = "doy")
gddstats <- merge(gddstats, gddstats3, by = "doy")

ggplot(gddstats, aes(x = doy, y = mean)) + 
  geom_line() +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, color = NA) +
  geom_vline(xintercept = 5) +
  labs(title = "GDD accumulation")

test <- data.frame(
  # aspp = rnorm(8000, 0, 1),
  aspp1 = rnorm(8000, 0.2, 1),
  aspp2 = rnorm(8000, -0.5, 1),
  aspp3 = rnorm(8000, -1, 1),
  aspp4 = rnorm(8000, 1, 1)
)
test 
simsub <- sim[sample(sample(1:40000, 1e4)),]
simsub <- subset(sim, spp %in% 1:4)

# get some basic stats
simcoef_summary <- aggregate(ringwidth ~ spp, data = simsub, function(x) {
  mean <- mean(x)
  serror <- sd(x) / sqrt(length(x))
  lower <- mean - 1.96 * serror
  upper <- mean + 1.96 * serror
  c(mean = mean, serror = serror, lower = lower, upper = upper)
})

# split
simcoef_summary <- do.call(data.frame, simcoef_summary)
names(simcoef_summary) <- c("spp", "mean_ringwidth", "serror", "lower", "upper")

# mean gdd and merge!
gdd_mean <- aggregate(gddcons ~ spp, data = simcoef, mean)
simcoef_summary <- merge(simcoef_summary, gdd_mean, by = "spp")

ggplot(simcoef_summary, aes(x = gddcons, y = factor(spp))) +
  geom_vline(xintercept = 0, color = "black") +
  
  # geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "grey50") +
  geom_point(size = 3, color = "black") +
  theme_minimal() +
  labs(
    x = "GDD (centered, gddcons)",
    y = "Species",
    title = "ringwidth by spp "
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot test parameters
test_long <- data.frame(
  spp = rep(names(test), each = nrow(test)),
  a_spp = as.vector(as.matrix(test))
)
test_long

ggplot(test_long, aes(x = a_spp, y = spp, group = spp, color = spp)) +
  geom_boxplot(position = "dodge2", outlier.shape = NA) +
  geom_vline(xintercept = 0, color = "black") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()
ggsave("figures/aspp_boxPlot.jpeg", width = 6, height = 6, units = "in", dpi = 300)
  

# plot fit parameters ####
###### Recover treeid ######
df_fit <- as.data.frame(fit)

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("treeid", colnames(df_fit))]
# remove sigma_asp for now
treeid_cols <- treeid_cols[1:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub(".*treeid:([^]]+)\\].*", "\\1", colnames(treeid_df))

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
aspp_cols <- colnames(df_fit)[grepl("spp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("Sigma", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub(".*spp:([^]]+)\\].*", "\\1", colnames(aspp_df))

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
###### Recover site ######
site_cols <- colnames(df_fit)[grepl("prov", colnames(df_fit))]
# remove for now
site_cols <- site_cols[1:length(site_cols)]
site_cols <- site_cols[!grepl("Sigma", site_cols)]
site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub(".*prov:([^]]+)\\].*", "\\1", colnames(site_df))
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

 
# Map
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# --- Get the map data ---
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- Define bounding box for northeastern North America ---
# Adjust these coordinates as needed
lat_min <- 40
lat_max <- 52
lon_min <- -82
lon_max <- -55

# --- Create example points along a latitudinal gradient ---
# These are arbitrary example locations
locations <- data.frame(
  name = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hyppolyte (Qc)"),
  lon = c(-72.2, -71, -70.66, -74.01),
  lat = c(42.55, 44.11, 44.92, 45.98)
)

special_point <- data.frame(
  name = "Arnorld Arboretum of 
  Harvard University (MA)",
  lon = -71.13358611669867,
  lat = 42.29601035316377
)
special_sf <- st_as_sf(special_point, coords = c("lon", "lat"), crs = 4326)

# Convert to sf object
points_sf <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)

# --- Plot the map ---
ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = points_sf, color = "#0A9F9D", size = 3) +
  geom_sf(data = special_sf, color = "#E54E21", shape = 8, size = 5, stroke = 1.2) +
  geom_text(data = locations, aes(x = lon, y = lat, label = name),
            nudge_y = 0.4, nudge_x = 0, size = 3.5) +
  geom_text(data = special_point,
            aes(x = lon, y = lat, label = name),
            nudge_y = -0.5, nudge_x = 6, color = "#E54E21", fontface = "bold") +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Source populations of seeds and location of the common garden",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted")
  )
ggsave("figures/mapSourcePop.jpeg", width = 6, height = 4, units = "in", dpi = 300)

# Fitting empirical data with stan_lmer ####

