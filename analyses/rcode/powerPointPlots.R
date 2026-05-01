# plots for power point
# Wildchrokie model
# CRD 20 April 2026

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
library(rsvg)
library(shape)

# flags
makeplots <- TRUE
runzscore <- F
# interceptmuplots <- TRUE
par(family = "Helvetica")

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
climatesum <- read.csv("output/climateSummariesYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

emp$latbi[which(emp$latbi %in% "Alnus incana")] <- "A. incana"
emp$latbi[which(emp$latbi %in% "Betula alleghaniensis")] <- "B. alleghaniensis"
emp$latbi[which(emp$latbi %in% "Betula papyrifera")] <- "B. papyrifera"
emp$latbi[which(emp$latbi %in% "Betula populifolia")] <- "B. populifolia"

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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### Define objects used throught the models #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
sitefull <- c(
  "GR" = "Dartmouth College (NH)",
  "HF" = "Harvard Forest (MA)",
  "SH" = "St-Hippolyte (Qc)",
  "WM" = "White Mountains (NH)"
)

locations <- data.frame(
  name       = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hippolyte (Qc)"),
  shortnames = c("HF", "WM", "GR", "SH"),   
  Longitude  = c(-72.2,  -71.0, -70.66, -74.01),
  Latitude   = c( 42.55,  44.11,  44.92,  45.98)
)

site_order <- locations$shortnames[order(locations$Latitude)]
y_pos_site <- match(site_df2$site_name, site_order)

# shapes for sites
my_shapes <- c(HF = 19, WM = 18, GR = 15, SH = 17)

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

species_order <- c(
  "A. incana", 
  "B. alleghaniensis", 
  "B. papyrifera", 
  "B. populifolia")

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Combined mu plots bspp (GDD / GSL / SOS / EOS) ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
custommar <- c(4, 4, 2, 1.2)
img_thermom <- rsvg::rsvg("figures/pictogramsLeaves/thermometer.svg")
img_calenda <- rsvg::rsvg("figures/pictogramsLeaves/calendar.svg")
img_leafout <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicLeafout.svg")
img_budset  <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicBudset.svg")
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### bspp SOS ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# pdf("figures/powerPoint/muSOS.pdf", width = 7.5, height = 6)
jpeg(file = "figures/powerPoint/muSOS.jpeg",
     width = 2400, height = 1800, res = 400)

mumar <- c(4, 1, 4, 1)

par(mfrow = c(1,1))
plot(bspp_df2_sos$mean, y_pos,
     xlim = c(-0.3, 0.4), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change per 5 days of leafout", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2_sos$p5,  y_pos, bspp_df2_sos$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_sos$p25, y_pos, bspp_df2_sos$p75, y_pos, col = wccolslatbi, lwd = 3)
# mtext("(c) Start of season", adj = 0, side = 3, line = 2.5, font = 2, cex = 0.9)
# arrows(x0 = -0.05, y0 = n_spp + 0.85, x1 = -0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
# text(-0.18, n_spp + 0.85, "Larger/Earlier", pos = 3, xpd = TRUE, cex = 1.3)
# arrows(x0 = 0.05, y0 = n_spp + 0.85, x1 = 0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
# text(0.18, n_spp + 0.85, "Smaller/Later", pos = 3, xpd = TRUE, cex = 1.3)
# usr <- par("usr")
# rasterImage(img_leafout, usr[1], usr[4] - diff(usr[3:4]) * 0.35, usr[1] + diff(usr[1:2]) * 0.25, usr[4])

legend("right",
       legend = sapply(unique(bspp_df2$spp_name),
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(wccolslatbi),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### bspp EOS ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/powerPoint/muEOS.jpeg", width = 2400, height = 1800, res = 400)
par(mfrow = c(1,1))
plot(bspp_df2_eos$mean, y_pos,
     xlim = c(-0.3, 0.4), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change per 10 days of budset", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2_eos$p5,  y_pos, bspp_df2_eos$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_eos$p25, y_pos, bspp_df2_eos$p75, y_pos, col = wccolslatbi, lwd = 3)
# mtext("(d) End of season", adj = 0, side = 3, line = 2.5, font = 2, cex = 0.9)
# arrows(x0 = -0.05, y0 = n_spp + 0.85, x1 = -0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
# text(-0.18, n_spp + 0.85, "Larger/Earlier", pos = 3, xpd = TRUE, cex = 1.3)
# arrows(x0 = 0.05, y0 = n_spp + 0.85, x1 = 0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
# text(0.18, n_spp + 0.85, "Smaller/Later", pos = 3, xpd = TRUE, cex = 1.3)
# usr <- par("usr")
# rasterImage(img_budset, usr[1], usr[4] - diff(usr[3:4]) * 0.35, usr[1] + diff(usr[1:2]) * 0.25, usr[4])
legend("right",
       legend = sapply(unique(bspp_df2$spp_name),
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(wccolslatbi),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### bspp GDD ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/powerPoint/muGDD.jpeg", width = 2400, height = 1800, res = 400)
par(mfrow = c(1,1), mar = c(6, 7, 2, 4))
plot(bspp_df2$mean, y_pos,
     xlim = c(-0.2, 0.6), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change with longer thermal seasons (GDD)", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, lwd = 0.9, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2$p5,  y_pos, bspp_df2$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2$p25, y_pos, bspp_df2$p75, y_pos, col = wccolslatbi, lwd = 3)

# y-axis species labels
spp_labels <- sapply(unique(bspp_df2_gsl$spp_name),
                     function(x) parse(text = paste0("italic('", x, "')")))
axis(2, at = y_pos, labels = spp_labels, las = 1, tick = TRUE, cex.axis = 0.9)

dev.off()



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### bspp GSL ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/powerPoint/muGSL.jpeg", width = 2400, height = 1800, res = 400)
par(mfrow = c(1,1), mar = c(6, 7, 2, 4))  # widen left margin for species names
plot(bspp_df2_gsl$mean, y_pos,
     xlim = c(-0.2, 0.6), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change with longer calendar season", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, lwd = 0.9, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2_gsl$p5,  y_pos, bspp_df2_gsl$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_gsl$p25, y_pos, bspp_df2_gsl$p75, y_pos, col = wccolslatbi, lwd = 3)

# y-axis species labels
spp_labels <- sapply(unique(bspp_df2$spp_name),
                     function(x) parse(text = paste0("italic('", x, "')")))
axis(2, at = y_pos, labels = spp_labels, las = 1, tick = TRUE, cex.axis = 0.9)

dev.off()



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Conceptual figure broken down ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# pdf("figures/climate/gsconceptualfig.pdf", width = 8, height = 8)


##### Set common things that will be repeated #####
# cols
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
peak <- 185


# Logistic curves
doy <- 1:330
gdd_pre <- 2500 / (1 + exp(-0.03 * (doy - 172)))
gdd_cc  <- 3000 / (1 + exp(-0.04 * (doy - 140)))

# Normal curves
x <- as.integer(seq(0, 330, length.out = 330))
xs <- seq(0, peak, length.out = 330)
xf <- seq(peak, 330 , length.out = 330)

y1 <- dnorm(x, mean = peak, sd = 45)
ys <- dnorm(xs, mean = peak, sd = 55)
yf <- dnorm(xf, mean = peak, sd = 50)

# Scale y2 so both curves peak at the same height
scale_to <- 25

y1_scaled <- y1  * (scale_to / max(y1))
ys_scaled2 <- ys  * (scale_to / max(ys))
yf_scaled2 <- yf  * (scale_to / max(yf))
y2_scaled <- c(ys_scaled2, yf_scaled2)
length(ys_scaled2) + length(yf_scaled2)
# calendar days
ticks <- seq(0, 330, by = 30)  # wherever you want ticks
dates <- format(as.Date(ticks, origin = "2023-01-01"), "%d %b")

myxlimp3 <- c(min(doy), max(doy))

mylwd <- 3

# Panel margins
p1 <- c(0, 5, 2, 2)
p2 <- c(0, 5, 0, 2)
p3 <- c(0, 5, 0, 2)

# matrix heights 
matheights <- c(2.4, 2.8, 1)

# ylim logistic
ylimlogis <- c(0, 3400)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### 1. Basic curve with just temp and gdd #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg("figures/powerPoint/concept1.jpeg", width = 10, height = 8, units = "in", res = 400)
layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# Panel 1
par(mar = p1)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 3))

# Panel 2 with just pre industrial
par(mar = p2)
plot(x, y1_scaled, type = "l", lwd = mylwd, col = colpre, xaxt = "n",
     ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE, 
     cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)

# abline for 5C
segments(x0 = -11, x1 = 126, y0 = 11, y1 = 11, lty = 2)
Arrows(x0 = 120, y0 = 20, x1 = 120, y1 = 12.8, arr.type = "triangle", 
       arr.width = 0.4, lwd	= 2)
text(x = 110, y = 22, labels = expression("16"*degree*"C, GDD = 11"), 
     cex = 1.5, col = "black")

# Panel 3
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))
dev.off()


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### 2. Pre-CC accumulated GDD, with arrows of how many is accumulated daily #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg("figures/powerPoint/concept2.jpeg", width = 10, height = 8, units = "in", res = 400)
layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# Panel 1 --- ---  --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p1)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 3))

# Panel 2 : with just pre industrial --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)

plot(x, y1_scaled, type = "l", lwd = mylwd, col = colpre, xaxt = "n",
     ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE, 
     cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)

# Early spring
Arrows(x0 = 82, y0 = 7, x1 = 82, y1 = 3.5, arr.type = "triangle", 
       arr.width = 0.3, lwd	= 2)
segments(x0 = -12, y0 = 2, x1 = 85, y1 = 2, lty = 2, lwd = 1.2)
text(x = 55, y = 8.5, labels = "Early spring: 2 daily GDD", 
     cex = 1.5, col = "black")

# Mid summer
Arrows(x0 = 150, y0 = 25, x1 = 150, y1 = 21.5, arr.type = "triangle", 
       arr.width = 0.3, lwd	= 2)
segments(x0 = -12, y0 = 20, x1 = 152, y1 = 20, lty = 2, lwd = 1.2)
text(x = 130, y = 26, labels = "Mid summer: 20 daily GDD", 
     cex = 1.5, col = "black")

# Panel 3  --- ---  --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### 3. Pre-CC accumulated GDD, with arrows of how many is accumulated overall #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg("figures/powerPoint/concept3.jpeg", width = 10, height = 8, units = "in", res = 400)
layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# Panel 1 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p1)
# Pre CC
plot(doy, gdd_cc, ylim = ylimlogis,
     type = "n", lwd = mylwd,
     xlab = "", ylab = "Accumulated GDD",
     xaxt = "n", frame = FALSE, col = adjustcolor(colpre), main = "",
     cex.axis = axissize, cex.lab = labsize)

# draw logistic curves pre CC
lines(doy, gdd_pre, type = "l", lwd = mylwd, col = adjustcolor(colpre))

# Early spring
Arrows(x0 = 79, y0 = gdd_pre[79] + 600, x1 = 79, y1 = gdd_pre[79] + 200,
       arr.type = "triangle", arr.width = 0.3, lwd = mylwd)
segments(x0 = -12, y0 = gdd_pre[79], x1 = 79, y1 = gdd_pre[79], lty = 2, lwd = 1.2)
text(x = 57, y = gdd_pre[79] + 700, labels = "Early spring: 50 GDD",
     cex = 1.5, col = "black")

# Mid summer
Arrows(x0 = 158, y0 = gdd_pre[160] + 650, x1 = 158, y1 = gdd_pre[160] + 250,
       arr.type = "triangle", arr.width = 0.3, lwd = mylwd)
segments(x0 = -12, y0 = gdd_pre[160], x1 = 157, y1 = gdd_pre[160], lty = 2, lwd = 1.2)
text(x = 130, y = gdd_pre[160] + 800, labels = "Mid summer: 1000 GDD",
     cex = 1.5, col = "black")

# Panel 2 : with just pre industrial --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)

plot(x, y1_scaled, type = "l", lwd = mylwd, col = colpre, xaxt = "n",
     ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE, 
     cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)

# Early spring
Arrows(x0 = 82, y0 = 7, x1 = 82, y1 = 3.5, arr.type = "triangle", 
       arr.width = 0.3, lwd	= 2)
segments(x0 = -12, y0 = 2, x1 = 85, y1 = 2, lty = 2, lwd = 1.2)
text(x = 55, y = 8.5, labels = "Early spring: 2 daily GDD", 
     cex = 1.5, col = "black")

# Mid summer
Arrows(x0 = 150, y0 = 25, x1 = 150, y1 = 21.5, arr.type = "triangle", 
       arr.width = 0.3, lwd	= 2)
segments(x0 = -12, y0 = 20, x1 = 152, y1 = 20, lty = 2, lwd = 1.2)
text(x = 130, y = 26, labels = "Mid summer: 20 daily GDD", 
     cex = 1.5, col = "black")

# Panel 3 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### 4. GS arrow #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg("figures/powerPoint/concept4.jpeg", width = 10, height = 8, units = "in", res = 400)
layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# Panel 1 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p1)
# Pre CC
plot(doy, gdd_cc, ylim = ylimlogis,
     type = "n", lwd = mylwd,
     xlab = "", ylab = "Accumulated GDD",
     xaxt = "n", frame = FALSE, col = adjustcolor(colpre), main = "",
     cex.axis = axissize, cex.lab = labsize)

# draw logistic curves pre CC
lines(doy, gdd_pre, type = "l", lwd = mylwd, col = adjustcolor(colpre))

# Panel 2 : with just pre industrial --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)
plot(x, y1_scaled, type = "l", lwd = mylwd, col = colpre, xaxt = "n",
     ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE, 
     cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)

# gs boundaries
segments(x0 = preeos, x1 = preeos, y0 = -2, y1 = y1_scaled[presos] - 4, lwd = 1, lty = 2)
segments(x0 = presos, x1 = presos, y0 = -2, y1 = y1_scaled[preeos] + 2.8, lwd = 1, lty = 2)

# Panel 3 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))

# gsl arrow line pre CC
arrow_y <- 0.5
shaft_h <- 0.15
head_h  <- 0.15
x_start <- presos
x_end <- preeos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colpre, alpha.f = 0.7),
  border = NA)

# 0a6a3c
img_leafout <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicLeafout.svg")
img_budset  <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicBudset.svg")
usr <- par("usr")

img_w <- 23  # width in plot units, tune as needed
img_h <- 0.6  # height in plot units, tune as needed

rasterImage(img_leafout, x_start - 12 - img_w/2, arrow_y - img_h/2, x_start -12 + img_w/2, arrow_y + img_h/2)
rasterImage(img_budset,  x_end + 7 - img_w/2,   arrow_y - img_h/2, x_end + 7 + img_w/2,   arrow_y + img_h/2)

text(x = presos + (preeos - presos)/2, y = arrow_y, 
     "Growing season length", col = "black", cex = 1.9)

dev.off()

 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### 5. GS extension with Climate change WITHOUT Logistic CC #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg("figures/powerPoint/concept5.jpeg", width = 10, height = 8, units = "in", res = 400)

layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# gdd_pre <- 2500 / (1 + exp(-0.03 * (doy - 172)))
# gdd_cc  <- 3200 / (1 + exp(-0.03 * (doy - 140)))

# Panel 1 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p1)
# Pre CC
plot(doy, gdd_cc, ylim = ylimlogis,
     type = "n", lwd = 1.2,
     xlab = "", ylab = "Accumulated GDD",
     xaxt = "n", frame = FALSE, col = adjustcolor(colpre, alpha.f = 0.4), main = "",
     cex.axis = axissize, cex.lab = labsize)

# draw logistic curves pre CC
lines(doy, gdd_pre, type = "l", lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.4))

# Panel 2 : with just pre industrial --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)
plot(x, y1_scaled, type = "l", lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.4), 
     xaxt = "n", ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE, 
     cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)

lines(xs, ys_scaled2, lwd = mylwd, col = colcc)
lines(xf, yf_scaled2, lwd = mylwd, col = colcc)

# gs boundaries
segments(x0 = ccsos, x1 = ccsos, y0 = -2, y1 = 9.2, lwd = 1, lty = 2)
segments(x0 = cceos, x1 = cceos, y0 = -2, y1 = 8, lwd = 1, lty = 2)

# Panel 3 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))

# gsl arrow line CC RED
arrow_y <- 0.5
shaft_h <- 0.15
head_h  <- 0.15
x_start <- ccsos
x_end <- cceos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colcc, alpha.f = 0.7),
  border = NA)

# 0a6a3c
img_leafout <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicLeafout.svg")
img_budset  <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicBudset.svg")
usr <- par("usr")

img_w <- 23  # width in plot units, tune as needed
img_h <- 0.6  # height in plot units, tune as needed

rasterImage(img_leafout, x_start - 12 - img_w/2, arrow_y - img_h/2, x_start -12 + img_w/2, arrow_y + img_h/2)
rasterImage(img_budset,  x_end + 7 - img_w/2,   arrow_y - img_h/2, x_end + 7 + img_w/2,   arrow_y + img_h/2)

text(x = ccsos + (cceos - ccsos)/2, y = arrow_y, 
     "Longer calendar season", col = "black", cex = 1.9)

# gsl arrow line pre CC BLUE
arrow_y <- 0.2
shaft_h <- 0.06
head_h  <- 0.06
x_start <- presos
x_end <- preeos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colpre, alpha.f = 0.4),
  border = NA)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### 6. GS extension with Climate change WITH Logistic CC #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg("figures/powerPoint/concept6.jpeg", width = 10, height = 8, units = "in", res = 400)

layout(matrix(c(1, 2, 3), nrow = 3), heights = matheights)

# gdd_pre <- 2500 / (1 + exp(-0.03 * (doy - 172)))
# gdd_cc  <- 3200 / (1 + exp(-0.03 * (doy - 140)))

# Panel 1 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p1)
# Pre CC
plot(doy, gdd_cc, ylim = ylimlogis,
     type = "n", lwd = 1.2,
     xlab = "", ylab = "Accumulated GDD",
     xaxt = "n", frame = FALSE, col = adjustcolor(colpre, alpha.f = 0.4), main = "",
     cex.axis = axissize, cex.lab = labsize)

# draw logistic curves pre CC
lines(doy, gdd_pre, type = "l", lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.4))
# draw logistic curves CC
lines(doy, gdd_cc, type = "l", lwd = mylwd, col = adjustcolor(colcc))

text(x = 280, y = max(gdd_cc) + 150, "Warmer thermal season", col = "black", cex = 1.9)

Arrows(x0 = 300, y0 = max(gdd_pre), x1 = 300, y1 = max(gdd_cc) - 200, arr.type = "triangle",
       arr.width = 0.3, lwd	= 2, col = colcc)

# GS delimitations
# abline(v = ccsos + 0.3, lwd = 1, lty = 2)
# abline(v = cceos + 0.3, lwd = 1, lty = 2)
segments(x0 = ccsos, x1 = ccsos, y0 = -5, y1 = 680, lwd = 1, lty = 2)
segments(x0 = cceos, x1 = cceos, y0 = -5, y1 = 2900, lwd = 1, lty = 2)

# Panel 2 : with just pre industrial --- --- --- --- --- --- --- --- --- --- ---
par(mar = p2)
plot(x, y1_scaled, type = "l", lwd = mylwd, col = adjustcolor(colpre, alpha.f = 0.4), 
     xaxt = "n", ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE, 
     cex.axis = axissize, cex.lab = labsize)
axis(1, at = ticks, labels = dates, cex.axis = axissize)

lines(xs, ys_scaled2, lwd = mylwd, col = colcc)
lines(xf, yf_scaled2, lwd = mylwd, col = colcc)

# GS delimitations
abline(v = ccsos, lwd = 1, lty = 2)
abline(v = cceos, lwd = 1, lty = 2)

# Panel 3 --- --- --- --- --- --- --- --- --- --- ------ --- --- --- --- --- --- 
par(mar = p3)
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))

# gsl arrow line pre CC
arrow_y <- 0.5
shaft_h <- 0.15
head_h  <- 0.15
x_start <- ccsos
x_end <- cceos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colcc, alpha.f = 0.7),
  border = NA)

# 0a6a3c
img_leafout <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicLeafout.svg")
img_budset  <- rsvg::rsvg("figures/pictogramsLeaves/bepaPicBudset.svg")
usr <- par("usr")

img_w <- 23  # width in plot units, tune as needed
img_h <- 0.6  # height in plot units, tune as needed

rasterImage(img_leafout, x_start - 12 - img_w/2, arrow_y - img_h/2, x_start -12 + img_w/2, arrow_y + img_h/2)
rasterImage(img_budset,  x_end + 7 - img_w/2,   arrow_y - img_h/2, x_end + 7 + img_w/2,   arrow_y + img_h/2)

text(x = ccsos + (cceos - ccsos)/2, y = arrow_y, 
     "Longer calendar season", col = "black", cex = 1.9)

# gsl arrow line pre CC
arrow_y <- 0.2
shaft_h <- 0.06
head_h  <- 0.06
x_start <- presos
x_end <- preeos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colpre, alpha.f = 0.4),
  border = NA)

dev.off()

# === === === === === === === === === === === === === === === === === === === ===


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Slot 2: Season length arrows
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = c(0, 4, 0, 2))
plot.new()
plot.window(xlim = myxlimp3, ylim = c(0, 1))

# gsl arrow line pre CC
arrow_y <- 0.5
shaft_h <- 0.08
head_h  <- 0.08
x_start <- presos
x_end <- preeos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colpre, alpha.f = 0.7),
  border = NA
)
#0a6a3c
text(x = presos + (preeos - presos)/2, y = arrow_y, "GSL pre CC", col = "black", cex = 1)

# gsl arrow line CC
arrow_y <- 0.3
shaft_h <- 0.08
head_h  <- 0.08
x_start <- ccsos
x_end <- cceos
x_neck_l <- x_start + 10  
x_neck_r <- x_end - 10    

polygon(
  x = c(x_start, x_neck_l, x_neck_l, x_neck_r, x_neck_r, x_end, x_neck_r, x_neck_r, x_neck_l, x_neck_l),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y - head_h, arrow_y, arrow_y + head_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colcc, alpha.f = 0.7),
  border = NA
)
#0a6a3c
text(x = ccsos + (cceos - ccsos)/2, y = arrow_y, "GSL post CC", col = "black", cex = 1)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Slot 3: Daily GDD 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
par(mar = c(4, 4, 0, 2))

plot(x, y1_scaled, type = "l", lwd = mylwd, col = colpre,
     ylim = c(0, scale_to * 1.1),
     xlab = "Day of year", ylab = "Daily GDD", frame = FALSE)
lines(xs, ys_scaled2, lwd = mylwd, col = colcc)
lines(xf, yf_scaled2, lwd = mylwd, col = colcc)

# abline for 5C
abline(a = 5, b = 0, lty = 2)
text(x = 37, y = 5 + 2, labels = "GDD base temperature", 
     cex = 1.5, col = "black") 

# Spring arrow
arrow_y <- 30
shaft_h <- 1.2
head_h  <- 2
x_start <- ccsos
x_neck  <- ccsos + 12

polygon(
  x = c(x_start, x_neck, x_neck, x_start + 35, x_start + 35, x_neck, x_neck),
  y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
  col = adjustcolor(colspring, alpha.f = 0.7),
  border = NA
)

text(x = x_start - 24, y = arrow_y + 3, labels = "Earlier spring",
     adj = c(0, 0.5), cex = 1.5, col = "black") 

# # Autumn arrow
# arrow_y <- 30
# shaft_h <- 0.5
# head_h  <- 0.9
# x_start <- cceos
# x_neck  <- cceos - 12
# 
# polygon(
#   x = c(x_start, x_neck, x_neck, x_start - 20, x_start - 20, x_neck, x_neck),
#   y = c(arrow_y, arrow_y - head_h, arrow_y - shaft_h, arrow_y - shaft_h, arrow_y + shaft_h, arrow_y + shaft_h, arrow_y + head_h),
#   col = adjustcolor("#d39822", alpha.f = 0.7), border = F)
# 
# text(x = x_start - 12, y = arrow_y + 2, labels = "Later autumn",
#      adj = c(0, 0.5), cex = 0.8, col = "black")
# dev.off()
dev.off()
