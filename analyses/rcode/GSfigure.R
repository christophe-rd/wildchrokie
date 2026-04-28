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

source("rcode/tools.R")

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

# calculate daily gdd with caping to 0 with cold temp
gddyr$dgdd <- pmax(gddyr$meanTempC - 5, 0)
dgddagg <- aggregate(dgdd ~ doy, gddyr, FUN = mean)

plot(x = dgddagg$doy, y = dgddagg$dgdd,
     xlab = "", ylab = "gdd",
     pch = 16, frame = FALSE, cex = 0,
     # col = yearcolors[match(emp$year, years)],
     main = "")

# start of season average
lo <- aggregate(leafout ~ latbi, emp, FUN = mean)
bs <- aggregate(budset ~ latbi, emp, FUN = mean)
gslength <- merge(lo, bs, by = "latbi")


years <- unique(gddyr$year)

lines(dgddagg$doy, dgddagg$dgdd, col = "black", cex = 0.2)

for (i in seq_along(years)) { # i = 1
  
  year_dat <- gddyr[gddyr$year == years[i], ]
  
  # lm_fit <- lm(leafout ~ winterPptLeafout, data = year_dat)?>
  # x_seq  <- seq(min(year_dat$winterPptLeafout, na.rm = TRUE), 
  #               max(year_dat$winterPptLeafout, na.rm = TRUE), length.out = 200)
  # pred   <- predict(lm_fit, newdata = data.frame(winterPptLeafout = x_seq))
  # 
 # cumulated gdd 
    # lines(year_dat$do, year_dat$GDD_5, 
    #     col = "black",
    #     # col = yearcolors[i],
    #     lwd = 2)
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
# Conceptual figure ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
plot(x = dgddagg$doy, y = dgddagg$dgdd,
     xlab = "", ylab = "gdd",
     pch = 16, frame = FALSE, cex = 0,
     # col = yearcolors[match(emp$year, years)],
     main = "")

# start of season average
lo <- aggregate(leafout ~ latbi, emp, FUN = mean)
bs <- aggregate(budset ~ latbi, emp, FUN = mean)
gslength <- merge(lo, bs, by = "latbi")


years <- unique(gddyr$year)

lines(dgddagg$doy, dgddagg$dgdd, col = "black", lwd = 0.4)

abline(v = mean(lo$leafout))
text(x = mean(lo$leafout) - 15, y = 15, "SOS")
abline(v = mean(bs$budset))
text(x = mean(bs$budset) + 15, y = 15, "EOS")
arrows(x0 = mean(lo$leafout), x1 = mean(lo$leafout) -20,
       y0 = 12, y1 = 12)
for (i in seq_along(years)) { # i = 1
  
  year_dat <- gddyr[gddyr$year == years[i], ]
  
  # lm_fit <- lm(leafout ~ winterPptLeafout, data = year_dat)?>
  # x_seq  <- seq(min(year_dat$winterPptLeafout, na.rm = TRUE), 
  #               max(year_dat$winterPptLeafout, na.rm = TRUE), length.out = 200)
  # pred   <- predict(lm_fit, newdata = data.frame(winterPptLeafout = x_seq))
  # 
  # cumulated gdd 
  # lines(year_dat$do, year_dat$GDD_5, 
  #     col = "black",
  #     # col = yearcolors[i],
  #     lwd = 2)
  spp <- unique(gslength$latbi)
  y_base <- 5
  y_step <- 2
  for (s in seq_along(spp)) { # i = 1
    gs <- gslength[gslength$latbi == spp[s],]
    y_pos <- y_base + (s-1) * y_step
    segments(x0 = gs$leafout, x1 = gs$budset, y0 = y_pos, y1 = y_pos,
             col = wccolslatbi[gs$latbi])
  }
}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# kind of a heat map ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

years   <- unique(emp$year)
wcclimatesum <- subset(climatesum, year %in% years)

seasons <- c("DJF", "MAM", "JJA", "SON")

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
growth <- colMeans(rstan::extract(fityr)$"ayear")
names(growth) <- years

leafoutagg <- aggregate(leafout~ year, emp, FUN = mean)
leafoutagg$leafout <- leafoutagg$leafout - mean(emp$leafout, na.rm = TRUE)
leafout <- leafoutagg$leafout 
names(leafout) <- years

budsetagg <- aggregate(budset~ year, emp, FUN = mean)
budsetagg$budset <- budsetagg$budset - mean(emp$budset, na.rm = TRUE)
budset <- budsetagg$budset 
names(budset) <- years

gslagg <- aggregate(pgsGSL~ year, emp, FUN = mean)
gslagg$pgsGSL <- gslagg$pgsGSL - mean(emp$pgsGSL, na.rm = TRUE)
gsl <- gslagg$pgsGSL
names(gsl) <- years

# Define number of levels
n_levels <- 100

# temp palet
# https://colorbrewer2.org/#type=sequential&scheme=Blues&n=9
temp_pal <- colorRampPalette(c("#ffffcc",
                               "#ffeda0",
                               "#fed976",
                               "#feb24c",
                               "#fd8d3c",
                               "#fc4e2a",
                               "#e31a1c",
                               "#bd0026",
                               "#800026"))(n_levels)
# temp_pal <- colorRampPalette(c("#9cc184", "#192813"))(n_levels)

# pdsi palet
# pdsi_pal <- colorRampPalette(c("#74c8c3", "#0a2e57"))(n_levels)
pdsi_pal <- colorRampPalette(c(
  "#f7fbff",
  "#deebf7",
  "#c6dbef",
  "#9ecae1",
  "#6baed6",
  "#4292c6",
  "#2171b5",
  "#08519c",
  "#08306b"))(n_levels)

# growth palet (model parameter)
growth_pal <- colorRampPalette(c(
  "#dadaeb",
  "#bcbddc",
  "#9e9ac8",
  "#807dba",
  "#6a51a3",
  "#54278f"
))(n_levels)

# define limits
t_min <- -1.5
t_max <- 1.5
p_min <- -4
p_max <- 4
g_min <- min(growth) 
g_max <- max(growth)
lo_min <- min(leafout) 
lo_max <- max(leafout)
bs_min <- min(budset)  
bs_max <- max(budset)
gs_min <- min(gsl)     
gs_max <- max(gsl)

ny  <- length(years)
ns  <- length(seasons)

# cell dimensions (in plot units)
cw  <- 1.0   # full cell width  (season col)
ch  <- 1.0   # cell height      (year row)
sw  <- cw / 2  # subcell width (Temp or PDSI)
gw  <- 0.7   # growth col width
lw  <- 0.8   # year-label col width (left margin)
gap <- 0.15  # gap between season groups

pw <- 0.7 # pheno cell width (leafout, budset, gsl)
total_w <- lw + ns * cw + (ns - 1) * gap + gap + 3 * pw + gap + gw

# define horizontal positions 
pheno_x <- lw + ns * cw + (ns - 1) * gap + gap
growth_x <- pheno_x + 3 * pw + gap


pheno_pal <- colorRampPalette(c(
  "#f0f0f0",
  "#d9d9d9",
  "#bdbdbd",
  "#969696",
  "#737373",
  "#525252"
))(n_levels)

total_h <- ny * ch + 1

par(mar = c(0.5, 0.5, 0.5, 0.5))

# jpeg(
# filename = "figures/climate/heatGS.jpeg", 
# width = 3000, height = 2400, res = 300)

# open plot
plot.new()

# define the coordinates of the box
plot.window(xlim = c(0, total_w), ylim = c(0, total_h))

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
    # temp
    rect(x0, y0, x0 + sw, y0 + ch, col = tc, border = "white", lwd = 1.5)
    text(x0 + sw/2, y0 + ch/2, sprintf("%.1f", t_val), cex = 0.55, font = 2)
    # PDSI
    rect(x0 + sw, y0, x0 + cw, y0 + ch, col = pc, border = "white", lwd = 1.5)
    text(x0 + sw + sw/2, y0 + ch/2, sprintf("%.1f", p_val), cex = 0.55, font = 2)
  }
  
  # 5. Phenology cells
  y0_g <- (ny - i) * ch + 0.6
  lo_val <- leafout[i]
  lo_idx <- round((lo_val - lo_min) / (lo_max - lo_min) * (n_levels - 1)) + 1
  rect(pheno_x, y0_g, pheno_x + pw, y0_g + ch, col = pheno_pal[lo_idx], border = "grey30")
  text(pheno_x + pw/2, y0_g + ch/2, sprintf("%.0f", lo_val), cex = 0.62, font = 2)
  
  bs_val <- budset[i]
  bs_idx <- round((bs_val - bs_min) / (bs_max - bs_min) * (n_levels - 1)) + 1
  rect(pheno_x + pw, y0_g, pheno_x + 2*pw,     y0_g + ch, col = pheno_pal[bs_idx], border = "grey30")
  text(pheno_x + 1.5*pw, y0_g + ch/2, sprintf("%.0f", bs_val), cex = 0.62, font = 2)
  
  gs_val <- gsl[i]
  gs_idx <- round((gs_val - gs_min) / (gs_max - gs_min) * (n_levels - 1)) + 1
  rect(pheno_x + 2*pw, y0_g, pheno_x + 3*pw,     y0_g + ch, col = pheno_pal[gs_idx], border = "grey30")
  text(pheno_x + 2.5*pw, y0_g + ch/2, sprintf("%.0f", gs_val), cex = 0.62, font = 2)
  
  # 6. growth Cell
  g_val <- growth[i]
  g_idx <- round((g_val - g_min) / (g_max - g_min) * (n_levels - 1)) + 1
  gc <- growth_pal[g_idx]
  
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

# 4. pheno legend
phenox0_leg <- gx0_leg + bar_w + 0.35
for (k in 1:n_seg) {
  idx <- round((k - 1) / (n_seg - 1) * (n_levels - 1)) + 1
  rect(phenox0_leg + (k-1)*seg_w, legend_y0,
       phenox0_leg + k*seg_w, legend_y0 + bar_h,
       col = pheno_pal[idx], border = NA)
}
rect(phenox0_leg, legend_y0, phenox0_leg + bar_w, legend_y0 + bar_h, border = "grey30")
text(phenox0_leg + bar_w/2, legend_y0 + bar_h + 0.07, "Phenology (doy)/Season length", cex=0.65, font=2)
text(phenox0_leg, legend_y0 - 0.09, "Earlier/Shorter", cex = 0.7, adj = 0.3)
text(phenox0_leg + bar_w, legend_y0 - 0.09, "Later/Longer", cex = 0.7, adj = 0.7)


# season titles
for (j in 1:ns) {
  cx <- lw + (j - 1) * (cw + gap) + (cw / 2)
  ty <- total_h -0.2 # just a little above the season
  
  # season names
  text(cx, ty, seasons[j], cex = 0.85, font = 2, col = "grey15")
  
  # temp and pdsi markers just a little below the seasonsd
  text(cx + sw/2, total_h - 0.3, "PDSI", cex = 0.55, col = "grey40")
  text(cx - sw/2, total_h - 0.3, "Temp",    cex = 0.55, col = "grey40")
  
  # Pheno titles
  text(pheno_x + pw/2,   total_h - 0.22, "Leafout\n(anomalized)", cex = 0.8, font = 2, col = "grey15")
  text(pheno_x + 1.5*pw, total_h - 0.22, "Budset\n(anomalized)",  cex = 0.8, font = 2, col = "grey15")
  text(pheno_x + 2.5*pw, total_h - 0.22, "GSL\n(anomalized)",    cex = 0.8, font = 2, col = "grey15")
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