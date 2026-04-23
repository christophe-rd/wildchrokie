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
# flags
makeplots <- TRUE
runzscore <- F
# interceptmuplots <- TRUE

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
climatesum <- read.csv("output/climateSummariesYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

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
  "SH" = "St-Hyppolyte (Qc)",
  "WM" = "White Mountains (NH)"
)

locations <- data.frame(
  name       = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hyppolyte (Qc)"),
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
  "Alnus incana", 
  "Betula alleghaniensis", 
  "Betula papyrifera", 
  "Betula populifolia")

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
     width = 2400, height = 1800, res = 300)

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
arrows(x0 = -0.05, y0 = n_spp + 0.85, x1 = -0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(-0.18, n_spp + 0.85, "Larger/Earlier", pos = 3, xpd = TRUE, cex = 1.3)
arrows(x0 = 0.05, y0 = n_spp + 0.85, x1 = 0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(0.18, n_spp + 0.85, "Smaller/Later", pos = 3, xpd = TRUE, cex = 1.3)
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
jpeg(file = "figures/powerPoint/muEOS.jpeg", width = 2400, height = 1800, res = 300)
par(mfrow = c(1,1))
plot(bspp_df2_eos$mean, y_pos,
     xlim = c(-0.3, 0.4), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change per 10 days of budset", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2_eos$p5,  y_pos, bspp_df2_eos$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_eos$p25, y_pos, bspp_df2_eos$p75, y_pos, col = wccolslatbi, lwd = 3)
# mtext("(d) End of season", adj = 0, side = 3, line = 2.5, font = 2, cex = 0.9)
arrows(x0 = -0.05, y0 = n_spp + 0.85, x1 = -0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(-0.18, n_spp + 0.85, "Larger/Earlier", pos = 3, xpd = TRUE, cex = 1.3)
arrows(x0 = 0.05, y0 = n_spp + 0.85, x1 = 0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(0.18, n_spp + 0.85, "Smaller/Later", pos = 3, xpd = TRUE, cex = 1.3)
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
jpeg(file = "figures/powerPoint/muGDD.jpeg", width = 2400, height = 1800, res = 300)
par(mfrow = c(1,1))
plot(bspp_df2$mean, y_pos,
     xlim = c(-0.2, 1), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change in 10 spring days GDD", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2$p5,  y_pos, bspp_df2$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2$p25, y_pos, bspp_df2$p75, y_pos, col = wccolslatbi, lwd = 3)
# mtext("(d) End of season", adj = 0, side = 3, line = 2.5, font = 2, cex = 0.9)
arrows(x0 = -0.05, y0 = n_spp + 0.85, x1 = -0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(-0.18, n_spp + 0.85, "Smaller/Cooler", pos = 3, xpd = TRUE, cex = 1.3)
arrows(x0 = 0.05, y0 = n_spp + 0.85, x1 = 0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(0.18, n_spp + 0.85, "Larger/Warmer", pos = 3, xpd = TRUE, cex = 1.3)
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
##### bspp GSL ##### 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(file = "figures/powerPoint/muGSL.jpeg", width = 2400, height = 1800, res = 300)
par(mfrow = c(1,1))
plot(bspp_df2_gsl$mean, y_pos,
     xlim = c(-0.2, 1), ylim = c(0.5, n_spp + 0.5),
     xlab = "log(ring width) change per 10 days of GSL", ylab = "",
     yaxt = "n", pch = 16, cex = 2, col = wccolslatbi, frame.plot = TRUE,
     panel.first = abline(v = 0, lty = 2, col = "black"), cex.axis = 1.2, cex.lab = 1.2)
segments(bspp_df2_gsl$p5,  y_pos, bspp_df2_gsl$p95, y_pos, col = wccolslatbi, lwd = 1.5)
segments(bspp_df2_gsl$p25, y_pos, bspp_df2_gsl$p75, y_pos, col = wccolslatbi, lwd = 3)
# mtext("(d) End of season", adj = 0, side = 3, line = 2.5, font = 2, cex = 0.9)
arrows(x0 = -0.05, y0 = n_spp + 0.85, x1 = -0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(-0.18, n_spp + 0.85, "Smaller/Shorter", pos = 3, xpd = TRUE, cex = 1.3)
arrows(x0 = 0.05, y0 = n_spp + 0.85, x1 = 0.3, y1 = n_spp + 0.85, length = 0.1, xpd = TRUE)
text(0.18, n_spp + 0.85, "Larger/Longer", pos = 3, xpd = TRUE, cex = 1.3)
# usr <- par("usr")
# rasterImage(img_budset, usr[1], usr[4] - diff(usr[3:4]) * 0.35, usr[1] + diff(usr[1:2]) * 0.25, usr[4])
legend("right",
       legend = sapply(unique(bspp_df2$spp_name),
                       function(x) parse(text = paste0("italic('", x, "')"))),
       col    = unique(wccolslatbi),
       pch    = 16, pt.cex = 1.5, bty = "n", cex = 1.2,
       title  = "Species", title.font = 2)
dev.off()

