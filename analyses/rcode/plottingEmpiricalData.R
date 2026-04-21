# Wildchrokie model
# CRD 28 October April 2025
# Start plotting empirical data for wildchrokie

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 5)
# quartz()

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(wesanderson)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

sim <- read.csv("output/simdata.csv")
emp <- read.csv("output/empiricalDataMAIN.csv")
gdd <- read.csv("output/gddByYear.csv")

# add full species
emp$sppfull <- NA
emp$sppfull[which(emp$spp == "ALNINC")] <- "Alnus incana"
emp$sppfull[which(emp$spp == "BETPOP")] <- "Betula populifolia"
emp$sppfull[which(emp$spp == "BETPAP")] <- "Betula papyrifera"
emp$sppfull[which(emp$spp == "BETALL")] <- "Betula alleghaniensis"
emp$lengthMM <- emp$lengthCM*10

# color coded by number of frost free days
frostfree <- subset(gdd, minTempC > 0)

# count the nb of frost free days per year
countfrost <- frostfree %>% count(year)
colnames(countfrost) <- c("year", "countFrostFree")

emp <- merge(emp, countfrost, by = "year")

# associate cols
emp$colfrost <- NA
unique(emp$countFrostFree)
emp$colfrost[which(emp$year %in% c(2018, 2019))] <- "#81A88D"
emp$colfrost[which(emp$year %in% c(2020))] <- "#02401B"

emp$countFrostFree <- as.factor(emp$countFrostFree)
ggplot(emp, aes(x = pgsGDD5, y = lengthMM)) +
  geom_point(size = 2, alpha = 1,
             aes(shape = site,
                 color = countFrostFree, 
                 fill = countFrostFree)) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  
  facet_wrap(~sppfull) +
  labs(y = "Ring width (mm)", 
       x = "Growing degree days (GDD)", 
       color = "Number of frost free days",
       fill = "Number of frost free days",
       shape = "Site") +  
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 21)),
         fill = guide_legend(override.aes = list(shape = 21)))
ggsave("figures/empiricalData/sppLinearRegressions_pgsGDD_frostFreeDays.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Map ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# --- Get the map data ---
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- Define bounding box for northeastern North America ---
# Adjust these coordinates as needed
lat_min <- 41
lat_max <- 48
lon_min <- -78
lon_max <- -63

# --- Create example points along a latitudinal gradient ---
# These are arbitrary example locations
locations <- data.frame(
  name = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hyppolyte (Qc)"),
  lon = c(-72.2, -71, -70.66, -74.01),
  lat = c(42.55, 44.11, 44.92, 45.98)
)

locations2 <- locations[order(-locations$lat), ]

locations2$col <-  wes_palettes$Darjeeling1[1:4]

special_point <- data.frame(
  name = "Arnold Arboretum of 
  Harvard University (MA)",
  lon = -71.13358611669867,
  lat = 42.29601035316377
)
special_sf <- st_as_sf(special_point, coords = c("lon", "lat"), crs = 4326)

# Convert to sf object
points_sf <- st_as_sf(locations2, coords = c("lon", "lat"), crs = 4326)

# --- Plot the map ---
ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = points_sf, color = locations2$col, size = 4) +
  geom_sf(data = special_sf, color = "#E54E21", shape = 8, size = 6, stroke = 1.2) +
  geom_text(data = locations2, aes(x = lon, y = lat, label = name),
            nudge_y = 0.35, nudge_x = 0, size = 4.5, fontface = "bold") +
  geom_text(data = special_point,
            aes(x = lon, y = lat, label = name),
            nudge_y = 0.2, nudge_x = 2.5, color = "#E54E21", size = 5, fontface = "bold") +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  theme_minimal() +   
  theme(strip.text = element_blank(),                    
        legend.key.height = unit(1.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted")
  )
ggsave("figures/mapSourcePop.jpeg", width = 9, height = 6, units = "in", dpi = 300)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Time series phenological data ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
comb <- read.csv("output/uncleanedTimeseriesPheno.csv")
comb$yeardoy <- paste(comb$year, comb$doy, sep = "_")
comb2 <- comb[!duplicated(comb$yeardoy),]

minbb <- aggregate(budburst ~ year, emp, FUN = min)
meanbb <- aggregate(budburst ~ year, emp, FUN = mean)
maxbb <- aggregate(budburst ~ year, emp, FUN = max)
bb <- merge(minbb, meanbb, by = "year")
bb <- merge(bb, maxbb, by = "year")
colnames(bb) <- c("year", "min", "mean", "max")

jpeg(
  filename = "figures/empiricalData/phenoTimeseries.jpeg", 
  width = 2000, height = 2800, res = 400)

par(mfrow = c(length(unique(comb2$year)), 1))
for (yr in sort(unique(comb2$year))) { # i =1
  hist(comb2$doy[comb2$year == yr], breaks = seq(0, 366, by = 14),
       main = yr, xlab = "Day of year", xlim = c(0, 366), ylim = c(0,4),
       ylab = "Number of observations per 14 days")
  bbx <- bb[bb$year == yr,]
  segments(x0 = bbx$min, y0 = 0, y1 = 4, lwd = 2, lty = 2, col = "darkgreen")
  segments(x0 = bbx$mean, y0 = 0, y1 = 4, lwd = 2, lty = 2, col = "black")
  segments(x0 = bbx$max, y0 = 0, y1 = 4, lwd = 2, lty = 2, col = "orange2")
}
dev.off()

