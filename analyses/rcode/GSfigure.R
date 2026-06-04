# Wildchrokie (and eventually coringTreespotters) climate data with phenology
# CRD 9 April 2026

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)
options(max.print = 150)
options(digits = 3)

# Load library 
library(future)
library(patchwork) 
library(rsvg)
library(shape)
library(pollen)

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
# gddyr <- read.csv("output/gddByYear.csv")
gddyr <- read.csv("/Users/christophe_rouleau-desrochers/github/coringtreespotters/analyses/output/gddByYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

logan <- read.csv("/Users/christophe_rouleau-desrochers/github/coringtreespotters/analyses/output/loganLongTermCleaned.csv")

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

# plot(x = dgddagg$doy, y = dgddagg$dgdd,
#      xlab = "", ylab = "gdd",
#      pch = 16, frame = FALSE, cex = 0,
#      # col = yearcolors[match(emp$year, years)],
#      main = "")

# start of season average
lo <- aggregate(leafout ~ latbi, emp, FUN = mean)
bs <- aggregate(budset ~ latbi, emp, FUN = mean)
gslength <- merge(lo, bs, by = "latbi")


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Conceptual figure ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Set curves and stuff
# pre cc
colpre <- "#04a3bd"
colcc <- "#a00e00"

colspring <- "#247d3f"
colfall <- "#da7901"

axissize <- 1.2
labsize <- 1.5

# assign sos and eos values
ccsos <- 110
cceos <- 260
presos <- 130
preeos <- 250

# Logistic curves (Panel 2 still uses these)
doy_seq <- 30:330
gdd_pre <- 2500 / (1 + exp(-0.025 * (doy_seq - 172)))
gdd_cc  <- 3000 / (1 + exp(-0.025 * (doy_seq - 140)))

# calendar days
ticks <- seq(min(doy_seq), max(doy_seq), by = 30)
dates <- format(as.Date(ticks, origin = "2023-01-01"), "%d %b")

myxlimp3 <- c(min(doy_seq), max(doy_seq))
mylwd <- 3

# Panel margins
p1 <- c(0, 5, 0, 2)
p2 <- c(3, 5, 0, 2)
p3 <- c(0, 5, 0, 2)

# matrix heights
matheights <- c(1, 2.8, 2.4)

# ylim logistic
ylimlogis <- c(0, 3400)

# Real data from logan airport
# GDD for logistic curves
logan <- subset(logan, doy >= min(doy_seq) & doy <= max(doy_seq))
logan$GDD_5 <- NA

# Get unique years
years <- unique(logan$year)

# Loop through each year
for (y in years) {
  # Find rows for this year
  year_rows <- which(logan$year == y)
  
  # Calculate GDD for this year only
  logan$GDD_5[year_rows] <- gdd(tmax = logan$maxTempC[year_rows], 
                                tmin = logan$minTempC[year_rows], 
                                tbase = 5, 
                                type = "B")
}

baselineperiod <- subset(logan, year >1940 & year < 1981)
baselinemean <- mean(baselineperiod$meanTempC)

prewarm <- subset(logan, year > 1954 & year < 1976)
# prewarm$meanTempC <- prewarm$meanTempC - baselinemean
poswarm <- subset(logan, year > 2004 & year < 2026)
# poswarm$meanTempC <- poswarm$meanTempC - baselinemean

mean_pre <- aggregate(meanTempC ~ doy, data = prewarm, FUN = mean, na.rm = TRUE)
mean_pos  <- aggregate(meanTempC ~ doy, data = poswarm,  FUN = mean, na.rm = TRUE)

# Loess smoothing (span controls wiggliness: 0.3 = more wiggly, 0.5 = smoother)
loess_pre <- loess(meanTempC ~ doy, data = mean_pre, span = 0.4)
loess_pos  <- loess(meanTempC ~ doy, data = mean_pos,  span = 0.4)

smooth_pre <- predict(loess_pre, newdata = data.frame(doy = doy_seq))
smooth_cc  <- predict(loess_pos,  newdata = data.frame(doy = doy_seq))

# aggregate for GDD logistic curve
mean_pre_gdd <- aggregate(GDD_5 ~ doy, data = prewarm, FUN = mean, na.rm = TRUE)
mean_pos_gdd  <- aggregate(GDD_5 ~ doy, data = poswarm,  FUN = mean, na.rm = TRUE)

# Shaded area: below threshold on the PRE curve
threshold <- 5
x_poly <- c(doy_seq, rev(doy_seq))
y_poly <- c(pmin(smooth_pre, threshold), rev(rep(0, length(doy_seq))))
"#FFF2FF"
# y-axis limit based on real data
ylim_temp <- c(0,max(smooth_cc))

# rasterized pictograms
img_thermom <- rsvg::rsvg("figures/pictogramsLeaves/thermometer.svg")
img_calenda <- rsvg::rsvg("figures/pictogramsLeaves/calendar.svg")
img_leafout <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicLeafout.svg")
img_budset  <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicBudset.svg")

# Plot! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
jpeg("figures/climate/gsconceptualfig.jpeg", width = 11, height = 8, units = "in", res = 400)
layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# Panel 1 --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))

# CC season arrow
arrow_y <- 0.5
shaft_h <- 0.15
head_h  <- 0.15
x_start <- ccsos
x_end   <- cceos
x_neck_l <- x_start + 10
x_neck_r <- x_end - 10

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colcc, alpha.f = 0.7), border = NA)

# pictograms width and height
img_w <- 23
img_h <- 0.6

# pictograms scaler
smll <- 4.3
norm <- 2

rasterImage(img_calenda, 
            x_start + 125 - img_w/smll, 
            arrow_y - img_h/smll, 
            x_start + 125 + img_w/smll, 
            arrow_y + img_h/smll)
rasterImage(img_leafout, 
            x_start - 12 - img_w/norm, 
            arrow_y - img_h/norm, 
            x_start - 12 + img_w/norm, 
            arrow_y + img_h/norm)
rasterImage(img_budset,  
            x_end + 12 - img_w/norm,    
            arrow_y - img_h/norm, 
            x_end + 12 + img_w/norm,    
            arrow_y + img_h/norm)

text(x = ccsos + (cceos - ccsos)/2, y = arrow_y,
     "Longer calendar season", col = "black", cex = 1.9)

# Pre season arrow
arrow_y  <- 0.2
shaft_h  <- 0.06
head_h   <- 0.06
x_start  <- presos
x_end    <- preeos
x_neck_l <- x_start + 10
x_neck_r <- x_end - 10

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colpre, alpha.f = 0.4), border = NA)

text(x = ccsos + (cceos - ccsos)/2, y = arrow_y,
     "Pre climate change calendar season", col = "black", cex = 1)

# Panel 2: Temperature curves --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p1)
plot(doy_seq, smooth_pre, type = "n",
     xaxt = "n", ylim = ylim_temp,
     xlab = "",
     ylab = expression(paste("Temperature (", degree, "C)")), frame = FALSE,
     cex.axis = axissize, cex.lab = labsize)

# Shade area below threshold under pre curve
polygon(x_poly, y_poly, col = adjustcolor("grey", alpha.f = 0.6), border = NA)

lines(doy_seq, smooth_pre, lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.4))
lines(doy_seq, smooth_cc,  lwd = mylwd, col = colcc)

# # GS delimitations
# segments(x0 = ccsos, y0 = 0, y1 = smooth_cc[which.min(abs(doy_seq - ccsos))],  lwd = 1.5, lty = 2)
# segments(x0 = cceos, y0 = 0, y1 = smooth_cc[which.min(abs(doy_seq - cceos))],  lwd = 1.5, lty = 2)
# 
# # Pre-CC boundaries (lighter)
# segments(x0 = presos, y0 = 0, y1 = smooth_pre[which.min(abs(doy_seq - presos))], lwd = 0.3, lty = 2)
# segments(x0 = preeos, y0 = 0, y1 = smooth_pre[which.min(abs(doy_seq - preeos))], lwd = 0.3, lty = 2)

# GS delimitations
segments(x0 = ccsos, y0 = 0, y1 = 30,  lwd = 1.5, lty = 2)
segments(x0 = cceos, y0 = 0, y1 = 30,  lwd = 1.5, lty = 2)

# Pre-CC boundaries (lighter)
segments(x0 = presos, y0 = 0, y1 = 30, lwd = 0.3, lty = 2)
segments(x0 = preeos, y0 = 0, y1 = 30, lwd = 0.3, lty = 2)



# Phenology trend arrows
Arrows(x0 = ccsos + 20, y0 = 5, x1 = ccsos + 5, y1 = 5,
       arr.type = "triangle", arr.width = 0.3, lwd = 2, col = colspring)
Arrows(x0 = cceos - 10, y0 = 5, x1 = cceos - 2, y1 = 5,
       arr.type = "triangle", arr.width = 0.2, arr.lwd = 0.5, arr.length = 0.2, lwd = 2, col = colfall)

text(x = ccsos + 30, y = 7, "Earlier spring", col = colspring, cex = 1.9)
text(x = cceos - 20, y = 7, "Later fall",     col = colfall,   cex = 1.4)

# Segments that shows cooler temperature early in the season
Arrows(x0 = ccsos + 30, x1 = ccsos + 30,
       y0 = smooth_pre[which.min(abs(doy_seq - presos))], 
       y1 = smooth_pre[which.min(abs(doy_seq - ccsos))], 
       lwd = 1, lty = 1, col = "black", arr.type = "T", code = 3)
text(x = ccsos + 64, 
     y = mean(c(smooth_pre[which.min(abs(doy_seq - presos))], 
                smooth_pre[which.min(abs(doy_seq - ccsos))])) , 
     "Cooler first days of growth", col = "black",   cex = 1.4)

# a horizontal line with whiskers
# arrows(x0 = x - err, y0 = y,
#        x1 = x + err, y1 = y,
#        code = 3,        # arrowheads on both ends
#        angle = 90,      # flat = whisker caps
#        length = 0.05)   # cap width in inches

# Panel 3: GDD curves --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)
plot(doy_seq, gdd_cc, ylim = ylimlogis,
     type = "n", lwd = 1.2,
     xlab = "", ylab = "Accumulated GDD",
     xaxt = "n", 
     frame = FALSE,
     col = adjustcolor(colpre, alpha.f = 0.4),
     main = "", cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)


mean_pre_gdd <- subset(mean_pre_gdd, doy <= max(doy_seq))
mean_pos_gdd <- subset(mean_pos_gdd, doy <= max(doy_seq))

lines(doy_seq, mean_pre_gdd$GDD_5, type = "l", lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.4))
lines(doy_seq, mean_pos_gdd$GDD_5,  type = "l", lwd = mylwd, col = adjustcolor(colcc))

text(x = 188, y = 2750, "Warmer thermal season", col = "black", cex = 1.9)
img_w <- 23
img_h <- 3400 * 0.3
smll <- 4.3
norm <- 2
rasterImage(img_thermom, 
            x_start + 105 - img_w/smll, 
            2700 - img_h/smll, 
            x_start + 105 + img_w/smll, 
            2700 + img_h/smll)

Arrows(x0 = cceos, y0 = gdd_pre[200] + 200, x1 = 245, y1 = gdd_cc[200]-150,
       arr.type = "triangle", arr.width = 0.3, lwd = 2, col = colcc)

# Pre-CC boundaries (lighter)
segments(x0 = presos, y0 = -2, y1 = 3000, lwd = 0.3, lty = 2)
segments(x0 = preeos, y0 = -2, y1 = 3000, lwd = 0.3, lty = 2)

segments(x0 = ccsos, y0 = -100, y1 = 3000, lwd = 1.5, lty = 2)
segments(x0 = cceos, y0 = -100, y1 = 3000, lwd = 1.5, lty = 2)

dev.off()


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
years <- unique(gddyr$year)

# lines(dgddagg$doy, dgddagg$dgdd, col = "black", cex = 0.2)

# for (i in seq_along(years)) { # i = 1
#   
#   year_dat <- gddyr[gddyr$year == years[i], ]
#   
#   # lm_fit <- lm(leafout ~ winterPptLeafout, data = year_dat)?>
#   # x_seq  <- seq(min(year_dat$winterPptLeafout, na.rm = TRUE), 
#   #               max(year_dat$winterPptLeafout, na.rm = TRUE), length.out = 200)
#   # pred   <- predict(lm_fit, newdata = data.frame(winterPptLeafout = x_seq))
#   # 
#  # cumulated gdd 
#     # lines(year_dat$do, year_dat$GDD_5, 
#     #     col = "black",
#     #     # col = yearcolors[i],
#     #     lwd = 2)
#   spp <- unique(gslength$latbi)
#   y_base <- 100
#   y_step <- 50
#   for (s in seq_along(spp)) { # i = 1
#     gs <- gslength[gslength$latbi == spp[s],]
#     y_pos <- y_base + (s-1) * y_step
#     segments(x0 = gs$leafout, x1 = gs$budset, y0 = y_pos, y1 = y_pos)
#   }
# }

 # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# kind of a heat map ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
if(makeplots){
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
}