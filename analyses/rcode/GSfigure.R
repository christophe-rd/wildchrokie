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

# years for pre climate change and post climate change
preyr <- 1955:1975
ccyr <- 2005:2025

# Logistic curves (Panel 2 still uses these)
doy_seq <- 30:330
gdd_pre <- 2500 / (1 + exp(-0.025 * (doy_seq - 172)))
gdd_cc  <- 3000 / (1 + exp(-0.025 * (doy_seq - 140)))

# calendar days
ticks <- seq(min(doy_seq), max(doy_seq), by = 30)
dates <- format(as.Date(ticks, origin = "2023-01-01"), "%d %b")

myxlimp3 <- c(min(doy_seq), max(doy_seq))
mylwd <- 3
mysmalltxt <- 1.2
mylargetxt <- 1.9

# Panel margins
p1 <- c(0, 5, 0, 2)
p2 <- c(3, 5, 2, 2)
p3 <- c(0, 5, 0, 2)

# matrix heights
matheights <- c(1, 2.8, 2.4)

# ylim logistic
ylimlogis <- c(0, 3400)

# Real data from logan airport
prelogan <- subset(logan, year %in% preyr & doy >= presos & doy <= preeos)
cclogan <- subset(logan, year %in% ccyr & doy >= ccsos & doy <= cceos)

# GDD for logistic curves for both periods
prelogan$GDD_5 <- NA
cclogan$GDD_5 <- NA

# Start with pre climate change
years <- unique(prelogan$year)

# Loop through each year
for (y in years) {
  # Find rows for this year
  year_rows <- which(prelogan$year == y)
  
  # Calculate GDD for this year only
  prelogan$GDD_5[year_rows] <- gdd(tmax = prelogan$maxTempC[year_rows], 
                                   tmin = prelogan$minTempC[year_rows], 
                                   tbase = 5, 
                                   type = "B")
}
mean_pre_gdd <- aggregate(GDD_5 ~ doy, data = prelogan, FUN = mean, na.rm = TRUE)

# Then post climate change
years <- unique(cclogan$year)

# Loop through each year
for (y in years) {
  # Find rows for this year
  year_rows <- which(cclogan$year == y)
  
  # Calculate GDD for this year only
  cclogan$GDD_5[year_rows] <- gdd(tmax = cclogan$maxTempC[year_rows], 
                                  tmin = cclogan$minTempC[year_rows], 
                                  tbase = 5, 
                                  type = "B")
}
mean_cc_gdd  <- aggregate(GDD_5 ~ doy, data = cclogan,  FUN = mean, na.rm = TRUE)

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
arrow_y <- 0.3
cap_h <- 0.15
x_left <- ccsos
x_right <- cceos
segments(x0 = x_left, x1 = x_right, y0 = arrow_y, lwd = mylwd, col = colcc)
segments(x0 = x_left,  y0 = arrow_y - cap_h, y1 = arrow_y + cap_h, lwd = mylwd, col = colcc)
segments(x0 = x_right, y0 = arrow_y - cap_h, y1 = arrow_y + cap_h, lwd = mylwd, col = colcc)

# pictograms width and height
img_w <- 23
img_h <- 0.6
# pictograms scaler
smll <- 2
norm <- 2
rasterImage(img_leafout,
            x_left - 13 - img_w/norm,
            arrow_y - img_h/norm,
            x_left - 13 + img_w/norm,
            arrow_y + img_h/norm)
rasterImage(img_budset,
            x_right + 13 - img_w/norm,
            arrow_y - img_h/norm,
            x_right + 13 + img_w/norm,
            arrow_y + img_h/norm)
text(x = ccsos + (cceos - ccsos)/2 - 10 , y = arrow_y + 0.1,
     "Longer calendar season", col = "black", cex = mylargetxt)
arrow_y <- 0.4
rasterImage(img_calenda,
            x_left + 115 - img_w/smll,
            arrow_y - img_h/smll,
            x_left + 115 + img_w/smll,
            arrow_y + img_h/smll)

# Pre season arrow
arrow_y <- 0.1
cap_h <- 0.06
x_left <- presos
x_right <- preeos
segments(x0 = x_left, x1 = x_right, y0 = arrow_y, lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.8))
segments(x0 = x_left,  y0 = arrow_y - cap_h, y1 = arrow_y + cap_h, lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.8))
segments(x0 = x_right, y0 = arrow_y - cap_h, y1 = arrow_y + cap_h, lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.8))

# text(x = ccsos + (cceos - ccsos)/2, y = 0.25,
#      "Pre climate change calendar season", col = "black", cex = mysmalltxt)

# Panel 2: Temperature curves --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p1)
plot(doy_seq, smooth_pre, type = "n",
     xaxt = "n", ylim = ylim_temp,
     xlab = "", bty = "l",
     ylab = expression(paste("Temperature (", degree, "C)")), 
     # frame = FALSE,
     cex.axis = axissize, cex.lab = labsize)
mtext("(a)", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.3)

# Shade area below threshold under pre curve
polygon(x_poly, y_poly, col = adjustcolor("grey", alpha.f = 0.6), border = NA)

lines(doy_seq, smooth_pre, lwd = 1, col = adjustcolor(colpre, alpha.f = 0.7))
lines(doy_seq, smooth_cc,  lwd = 1, col = colcc)

# thicker lines within their respective leafout and budset timings
doypre <- presos : preeos
idx_pre <- which(doy_seq %in% doypre)
lines(doy_seq[idx_pre], smooth_pre[idx_pre], 
      lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.8))

doycc <- ccsos:cceos
idx_cc <- which(doy_seq %in% doycc)
lines(doy_seq[idx_cc], smooth_cc[idx_cc], lwd = mylwd, col = colcc)

# Cooler first days of growth delimitations
segments(x0 = presos, x1 = presos - 30, y0 = smooth_pre[doy_seq %in% presos],  
         lwd = 0.8, lty = 3, col = colpre)

segments(x0 = ccsos, x1 = presos - 30, y0 = smooth_cc[doy_seq %in% ccsos],  
         lwd = 0.8, lty = 3, col = colcc)

# Segments that shows cooler temperature early in the season
Arrows(x0 = presos - 35, x1 = presos - 35,
       y0 = smooth_cc[doy_seq %in% ccsos], 
       y1 = smooth_pre[doy_seq %in% presos], 
       lwd = 1, lty = 1, col = "black", arr.type = "T", code = 3)

text(x = presos - 72, 
     y = mean(c(smooth_pre[doy_seq %in% presos], 
                smooth_cc[doy_seq %in% ccsos])),
     "Early SOS = cooler first days of growth",
     col = "black", cex = mysmalltxt)

# End-of-season delimitations
segments(x0 = preeos, x1 = preeos - 10,
         y0 = smooth_pre[doy_seq %in% preeos],  
         lwd = 0.8, lty = 3, col = colpre)

segments(x0 = cceos, x1 = preeos - 10,
         y0 = smooth_cc[doy_seq %in% cceos],  
         lwd = 0.8, lty = 3, col = colcc)

# Difference at end of season
Arrows(x0 = preeos - 15, x1 = preeos - 15,
       y0 = smooth_cc[doy_seq %in% cceos], 
       y1 = smooth_pre[doy_seq %in% preeos], 
       lwd = 1, lty = 1, col = "black", arr.type = "T", code = 3)

text(x = preeos - 50, 
     y = mean(c(smooth_pre[doy_seq %in% preeos], 
                smooth_cc[doy_seq %in% cceos]))* 0.98,
     "Late EOS =\nlittle temperature difference",
     col = "black", cex = mysmalltxt, adj = 0)

# pre: start and end of thick segment
r <- 1.5  
filledcircle(r1 = r, mid = c(presos,  smooth_pre[doy_seq %in% presos]),
             col = adjustcolor(colpre, alpha.f = 1), lcol = NA)
filledcircle(r1 = r, mid = c(preeos,  smooth_pre[doy_seq %in% preeos]),
             col = adjustcolor(colpre, alpha.f = 1), lcol = NA)

# cc: start and end of thick segment
filledcircle(r1 = r, mid = c(ccsos, smooth_cc[doy_seq %in% ccsos]),
             col = colcc, lcol = NA)
filledcircle(r1 = r, mid = c(cceos, smooth_cc[doy_seq %in% cceos]),
             col = colcc, lcol = NA)

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

text(x = ccsos + 30, y = 7, "Earlier spring", col = colspring, cex = mylargetxt)
text(x = cceos - 20, y = 7, "Later fall",     col = colfall,   cex = mylargetxt)

# legend(x = ccsos - 80, y = 25, 
#        legend = c("Pre climate change",
#                   "Post climate change"),
#        bty = "o", lwd = 3, cex = 1.2,
#        col = c(colpre, colcc),
#        title = "Curves")

text(x = mean(mean_pre_gdd$doy) - 30 , y = max(smooth_pre) -9,
     "Pre climate change", col = adjustcolor(colpre, alpha.f = 1), 
     cex = mysmalltxt)
text(x = mean(mean_pre_gdd$doy) - 30 , y = max(smooth_cc) + -1,
     "Post climate change", col = colcc, 
     cex = mysmalltxt)

# Panel 3: GDD curves --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)
plot(doy_seq, gdd_cc, ylim = range(mean_cc_gdd$GDD_5),
     type = "n", lwd = 1.2,
     xlab = "", ylab = "Accumulated GDD",
     xaxt = "n", bty = "l",
     # frame = FALSE,
     col = adjustcolor(colpre, alpha.f = 0.4),
     main = "", cex.axis = axissize, cex.lab = labsize)
mtext("(b)", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.3)
axis(1, at = ticks, labels = dates, cex.axis = axissize)


# mean_pre_gdd <- subset(mean_pre_gdd, doy <= max(doy_seq))
# mean_cc_gdd <- subset(mean_cc_gdd, doy <= max(doy_seq))

lines(mean_pre_gdd$doy, mean_pre_gdd$GDD_5, type = "l", lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.8))
lines(mean_cc_gdd$doy, mean_cc_gdd$GDD_5,  type = "l", lwd = mylwd, col = adjustcolor(colcc))

# pre
filledcircle(r1 = r, mid = c(mean_pre_gdd$doy[1],  mean_pre_gdd$GDD_5[1]),
             col = adjustcolor(colpre, alpha.f = 1), lcol = NA)
filledcircle(r1 = r, mid = c(mean_pre_gdd$doy[nrow(mean_pre_gdd)], mean_pre_gdd$GDD_5[nrow(mean_pre_gdd)]),
             col = adjustcolor(colpre, alpha.f = 1), lcol = NA)

# cc
filledcircle(r1 = r, mid = c(mean_cc_gdd$doy[1], mean_cc_gdd$GDD_5[1]),
             col = adjustcolor(colcc), lcol = NA)
filledcircle(r1 = r, mid = c(mean_cc_gdd$doy[nrow(mean_cc_gdd)], mean_cc_gdd$GDD_5[nrow(mean_cc_gdd)]),
             col = adjustcolor(colcc), lcol = NA)

# Pre-CC boundaries (lighter)
segments(x0 = presos, y0 = -2, y1 = 3000, lwd = 0.3, lty = 2)
segments(x0 = preeos, y0 = -2, y1 = 3000, lwd = 0.3, lty = 2)

segments(x0 = ccsos, y0 = -100, y1 = 3000, lwd = 1.5, lty = 2)
segments(x0 = cceos, y0 = -100, y1 = 3000, lwd = 1.5, lty = 2)

# Early season gdd
# Segments that shows cooler temperature early in the season
segments(x0 = presos, x1 = presos - 30, 
         y0 = mean_pre_gdd$GDD_5[mean_pre_gdd$doy %in% presos],  
         lwd = 0.8, lty = 3, col = colpre)

segments(x0 = presos, x1 = presos - 30, 
         y0 = mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% presos],  
         lwd = 0.8, lty = 3, col = colcc)

Arrows(x0 = presos - 35, x1 = presos - 35,
       y0 = min(mean_cc_gdd$GDD_5), 
       y1 = mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% presos], 
       lwd = 1, lty = 1, col = "black", arr.type = "T", code = 3)

text(x = presos - 72,
     y = (mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% presos] +
            min(mean_cc_gdd$GDD_5)) / 2, 
     "Early SOS = little effect on GDD", 
     col = "black", cex = mysmalltxt)

# Late season gdd
segments(x0 = preeos, x1 = cceos - 30,
         y0 = mean_pre_gdd$GDD_5[mean_pre_gdd$doy %in% preeos],
         lwd = 0.8, lty = 3, col = colpre)

segments(x0 = cceos, x1 = cceos - 30,
         y0 = mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% cceos],
         lwd = 0.8, lty = 3, col = colcc)

Arrows(x0 = cceos - 35, x1 = cceos - 35,
       y0 = mean_pre_gdd$GDD_5[mean_pre_gdd$doy %in% preeos], 
       y1 = mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% cceos], 
       lwd = 1, lty = 1, col = "black", arr.type = "T", code = 3)

text(x = preeos - 65,
     y = (mean_pre_gdd$GDD_5[mean_pre_gdd$doy %in% preeos] +
            mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% cceos]) / 2, 
     "Late EOS = bigger effect on GDD", 
     col = "black", cex = mysmalltxt)

# Polygon for warmer thermal season
x_arrow <- cceos + 5

shaft_w <- 2
head_w  <- 4
head_l  <- 150

y_start <- mean_pre_gdd$GDD_5[mean_pre_gdd$doy %in% preeos]
y_end   <- mean_cc_gdd$GDD_5[mean_cc_gdd$doy %in% cceos]

direction <- sign(y_end - y_start)
y_neck <- y_end - direction * head_l

polygon(
  x = c(x_arrow - shaft_w,
        x_arrow - shaft_w,
        x_arrow - head_w,
        x_arrow,
        x_arrow + head_w,
        x_arrow + shaft_w,
        x_arrow + shaft_w),
  y = c(y_start,
        y_neck,
        y_neck,
        y_end,
        y_neck,
        y_neck,
        y_start),
  col = adjustcolor(colcc, alpha.f = 1),
  border = NA
)

text(x = x_arrow + 5, y = y_start + 200, 
     "Warmer thermal \nseason", col = "black", cex = mylargetxt, adj = 0)
img_w <- 23
img_h <- 3400 * 0.2
smll <- 4.3
norm <- 2

rasterImage(img_thermom, 
            x_arrow + 65 - img_w/norm, 
            y_end- 300 - img_h/norm, 
            x_arrow + 65 + img_w/norm, 
            y_end - 300 + img_h/norm)

text(x = mean(mean_pre_gdd$doy) - 30 , 
     y = 200,
     "Pre climate change", col = adjustcolor(colpre, alpha.f = 1), 
     cex = mysmalltxt, adj = 0)
text(x = mean(mean_pre_gdd$doy) - 30, 
     y = 1250,
     "Post climate change", col = colcc, cex = mysmalltxt, adj = 0)
 
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
