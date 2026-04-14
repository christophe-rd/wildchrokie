# Wildchrokie (and eventually coringTreespotters) climate data with phenology
# CRD 9 April 2026

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)
options(max.print = 150)
options(digits = 3)

# Load library 
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

# flags
makeplots <- FALSE

emp <- read.csv("output/empiricalDataMAIN.csv")
climatesum <- read.csv("output/climateSummariesYear.csv")
climatesummonth <- read.csv("output/climateSummariesByMonth.csv")
gddyr <- read.csv("output/gddByYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

# emp <- empir[!is.na(empir$leafout),]

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate data #### 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
emp$leafout <- as.integer(emp$leafout)
emp$budset <- as.integer(emp$budset)

gddyr$yeardoy <- paste(gddyr$year, gddyr$doy, sep = "_")
emp$yeardoybudburst <- paste(emp$year, emp$budburst, sep = "_")
emp$yeardoyleafout <- paste(emp$year, emp$leafout, sep = "_")
emp$yeardoybudset <- paste(emp$year, emp$budset, sep = "_")

plot(x = gddyr$doy, y = gddyr$GDD_5,
     xlab = "", ylab = "gdd",
     pch = 16, frame = FALSE, cex = 0,
     # col = yearcolors[match(emp$year, years)],
     main = "")

# start of season average
lo <- aggregate(leafout ~ latbi, emp, FUN = mean)
bs <- aggregate(budset ~ latbi, emp, FUN = mean)
gslength <- merge(lo, bs, by = "latbi")

years <- unique(gddyr$year)
for (i in seq_along(years)) { # i = 1
  
  year_dat <- gddyr[gddyr$year == years[i], ]
  
  # lm_fit <- lm(leafout ~ winterPptLeafout, data = year_dat)
  # x_seq  <- seq(min(year_dat$winterPptLeafout, na.rm = TRUE), 
  #               max(year_dat$winterPptLeafout, na.rm = TRUE), length.out = 200)
  # pred   <- predict(lm_fit, newdata = data.frame(winterPptLeafout = x_seq))
  # 
  lines(year_dat$do, year_dat$GDD_5, 
        col = "black",
        # col = yearcolors[i],
        lwd = 2)
  spp <- unique(gslength$latbi)
  y_base <- 100
  y_step <- 50
  for (s in seq_along(spp)) { # i = 1
    gs <- gslength[gslength$latbi == spp[s],]
    y_pos <- y_base + (s-1) * y_step
    segments(x0 = gs$leafout, x1 = gs$budset, y0 = y_pos, y1 = y_pos)
  }
}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# kind of a heat map ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

years   <- unique(emp$year)
wcclimatesum <- subset(climatesum, year %in% years)
seasons <- unique(climatesum$period)

# check the mean temperature for each season
anom <- aggregate(tmeanmean ~ period, climatesum, FUN = mean)

temp <- reshape(wcclimatesum,
                idvar = "year",
                timevar = "period",
                direction = "wide",
                drop = c("pdsi", "tmeanmax", "tmeanmin", "ppt", "rad"))
rownames(temp) <- temp$year
temp$year <- NULL
colnames(temp) <- c("DJF", "JJA", "MAM", "SON")
temp$DJF <- temp$DJF - anom$tmeanmean[which(anom$period %in% "DJF")]
temp$JJA <- temp$JJA - anom$tmeanmean[which(anom$period %in% "JJA")]
temp$MAM <- temp$MAM - anom$tmeanmean[which(anom$period %in% "MAM")]
temp$SON <- temp$SON - anom$tmeanmean[which(anom$period %in% "SON")]

# PDSI
pdsi <- reshape(wcclimatesum,
                idvar = "year",
                timevar = "period",
                direction = "wide",
                drop = c("tmeanmean", "tmeanmax", "tmeanmin", "ppt", "rad"))
rownames(pdsi) <- pdsi$year
pdsi$year <- NULL
colnames(pdsi) <- c("DJF", "JJA", "MAM", "SON")

# tree ring width
fityr <- readRDS("output/stanOutput/fitGrowthOnlyYear")
growth <- colMeans(extract(fityr)$"ayear")
names(growth) <- years


# Define number of levels
n_levels <- 100

# temp palet
temp_pal <- colorRampPalette(c("#9cc184", "#192813"))(n_levels)

# pdsi palet
pdsi_pal <- colorRampPalette(c("#74c8c3", "#0a2e57"))(n_levels)

# growth palet (model parameter)
growth_pal <- colorRampPalette(c("#a9845b", "#0f252f"))(n_levels)

# define limits
t_min <- -1.5
t_max <- 1.5
p_min <- -4
p_max <- 4
g_min <- min(growth) 
g_max <- max(growth)

ny  <- length(years)
ns  <- length(seasons)

# cell dimensions (in plot units)
cw  <- 1.0   # full cell width  (season col)
ch  <- 1.0   # cell height      (year row)
sw  <- cw / 2  # subcell width (Temp or PDSI)
gw  <- 0.7   # growth col width
lw  <- 0.8   # year-label col width (left margin)
gap <- 0.15  # gap between season groups

total_w <- lw + ns * cw + (ns - 1) * gap + gap + gw
total_h <- ny * ch + 1

par(mar = c(0.5, 0.5, 0.5, 0.5))

jpeg(
  filename = "figures/climate/heatGS.jpeg", 
  width = 3000, height = 2400, res = 300)

# open plot
plot.new()

# define the coordinates of the box
plot.window(xlim = c(0, total_w), ylim = c(0, total_h))

# x coordinates
growth_x <- lw + ns * cw + (ns - 1) * gap + gap

for (i in 1:ny) { # i = 1
  for (j in 1:ns) { # i = 1
    # 1. coordinates
    x0 <- lw + (j - 1) * (cw + gap)
    y0 <- (ny - i) * ch + 0.6
    
    # 2. temperature cols left subcell
    # calculate where the valeus sit between min max
    t_val <- temp[i, j]
    t_idx <- round((t_val - t_min) / (t_max - t_min) * (n_levels - 1)) + 1
    tc <- temp_pal[t_idx]
    
    # 3. PDSI colors right subcell
    p_val <- pdsi[i, j]
    p_idx <- round((p_val - p_min) / (p_max - p_min) * (n_levels - 1)) + 1
    pc    <- pdsi_pal[p_idx]
    
    # 4. draw rectangles
    # temp
    rect(x0, y0, x0 + sw, y0 + ch, col = tc, border = "white", lwd = 1.5)
    # PDSI
    rect(x0 + sw, y0, x0 + cw, y0 + ch, col = pc, border = "white", lwd = 1.5)
    # outer border
    rect(x0, y0, x0 + cw, y0 + ch, col = NA, border = "grey30", lwd = 1.2)
  }
  
  # 5. growth Cell
  y0_g <- (ny - i) * ch + 0.6
  g_val <- growth[i]
  g_idx <- round((g_val - g_min) / (g_max - g_min) * (n_levels - 1)) + 1
  gc    <- growth_pal[g_idx]
  
  rect(growth_x, y0_g, growth_x + gw, y0_g + ch, col = gc, border = "grey30")
  text(growth_x + gw/2, y0_g + ch/2, sprintf("%.2f", g_val), cex=0.72, font=2)
}
#legend under the first season col
legend_x0 <- lw           # Start legend under the first season column
# vertical position
legend_y0 <- 0.05
# Width of each color bar
bar_w     <- 1.6 
# Height of each color bar
bar_h     <- 0.18        
# Number of color segments in the bar
n_seg     <- 50           
# Width of a single segment
seg_w     <- bar_w / n_seg 

# 1. temp legend 
for (k in 1:n_seg) {
  rect(legend_x0 + (k-1)*seg_w, legend_y0, legend_x0 + k*seg_w, legend_y0 + bar_h,
       col = temp_pal[round((k/n_seg)*(n_levels-1))+1], border = NA)
}
text(legend_x0 + bar_w/2, legend_y0 + bar_h + 0.07, "Temp (Manual)", cex=0.65, font=2)
text(legend_x0, legend_y0 - 0.09, "Colder", cex = 0.7, adj = 0.3)
text(legend_x0 + bar_w, legend_y0 - 0.09, "Warmer", cex = 0.7, adj = 0.7)

# 2. PDSI legend aligned with its col
px0 <- legend_x0 + bar_w + 0.35 
for (k in 1:n_seg) {
  idx <- round((k - 1) / (n_seg - 1) * (n_levels - 1)) + 1
  rect(px0 + (k-1)*seg_w, legend_y0, 
       px0 + k*seg_w, legend_y0 + bar_h, 
       col = pdsi_pal[idx], border = NA)
}
rect(px0, legend_y0, px0 + bar_w, legend_y0 + bar_h, border = "grey30")
text(px0 + bar_w/2, legend_y0 + bar_h + 0.07, "PDSI", cex=0.65, font=2)
text(px0, legend_y0 - 0.09, "Drier", cex = 0.7, adj = 0.3)
text(px0 + bar_w, legend_y0 - 0.09, "Wetter", cex = 0.7, adj = 0.7)

# 3. growth gegend
gx0_leg <- px0 + bar_w + 0.35
for (k in 1:n_seg) {
  idx <- round((k - 1) / (n_seg - 1) * (n_levels - 1)) + 1
  rect(gx0_leg + (k-1)*seg_w, legend_y0, 
       gx0_leg + k*seg_w, legend_y0 + bar_h, 
       col = growth_pal[idx], border = NA)
}
rect(gx0_leg, legend_y0, gx0_leg + bar_w, legend_y0 + bar_h, border = "grey30")
text(gx0_leg + bar_w/2, legend_y0 + bar_h + 0.07, "Growth (mm)", cex=0.65, font=2)
text(gx0_leg, legend_y0 - 0.09, "More growth", cex = 0.7, adj = 0.3)
text(gx0_leg + bar_w, legend_y0 - 0.09, "Less growth", cex = 0.7, adj = 0.7)

# season titles
for (j in 1:ns) {
  cx <- lw + (j - 1) * (cw + gap) + (cw / 2)
  ty <- total_h -0.2 # just a little above the season
  
  # season names
  text(cx, ty, seasons[j], cex = 0.85, font = 2, col = "grey15")
  
  # temp and pdsi markers just a little below the seasonsd
  text(cx + sw/2, total_h - 0.3, "PDSI", cex = 0.55, col = "grey40")
  text(cx - sw/2, total_h - 0.3, "Temp",    cex = 0.55, col = "grey40")
}

# growth title
text(growth_x + gw / 2, total_h - 0.22, "Growth\n(mm)", 
     cex = 0.8, font = 2, col = "grey15")

for (i in 1:ny) {
  # Calculate the vertical center of the current row
  # This matches the 'y0' logic in your rectangle loop
  y_center <- (ny - i) * ch + 0.6 + (ch / 2)
  
  # Use the 'lw' (left width) constant to position the text
  # adj = 1 right-aligns the text against the grid
  text(x = lw - 0.1, 
       y = y_center, 
       labels = years[i], 
       cex = 0.85, 
       font = 2, 
       col = "grey15", 
       adj = 1) 
  
}
dev.off()