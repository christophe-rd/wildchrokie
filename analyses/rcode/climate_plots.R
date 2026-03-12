# Wildchrokie (and eventually coringTreespotters) climate data with phenology
# CRD 12 March 2026

# housekeeping
rm(list=ls())
options(stringsAsFactors = FALSE)
options(max.print = 150)
options(digits = 3)

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(wesanderson)
library(patchwork) 
library(cartography)

if (length(grep("christophe_rouleau-desrochers", getwd())) > 0) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if (length(grep("lizzie", getwd())) > 0) {
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
} else  {
  setwd("/home/crouleau/wildchrokie/analyses")
}

# flags
makeplots <- TRUE

# === === === === === === === === === === === === === === === === 
# EMPIRICAL DATA ####
# === === === === === === === === === === === === === === === === 
emp <- read.csv("output/empiricalDataMAIN.csv")
climatesum <- read.csv("output/climateSummariesYear.csv")
gddyr <- read.csv("output/gddByYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

commonNames <- c(
  "Alnus incana"          = "Grey alder",
  "Betula alleghaniensis" = "Yellow birch",
  "Betula papyrifera"     = "Paper birch",
  "Betula populifolia"    = "Gray birch"
)

emp$commonName <- commonNames[emp$latbi]

emp <- emp[!is.na(emp$pgsGDD5),]
# transform my groups to numeric values
emp$site_num <- match(emp$site, unique(emp$site))
emp$spp_num <- match(emp$spp, unique(emp$spp))
emp$treeid_num <- match(emp$treeid, unique(emp$treeid))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate data #### 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
emp4 <- emp[!is.na(emp$leafout),]
emp4$leafout <- as.integer(emp4$leafout)
emp4$budset <- as.integer(emp4$budset)

gddyr$yeardoy <- paste(gddyr$year, gddyr$doy, sep = "_")
emp4$yeardoybudburst <- paste(emp4$year, emp4$budburst, sep = "_")
emp4$yeardoyleafout <- paste(emp4$year, emp4$leafout, sep = "_")
emp4$yeardoybudset <- paste(emp4$year, emp4$budset, sep = "_")

yearcolors <- c("#931e18", "#da7901", "#247d3f")
par(mfrow = c(1,1))
years <- sort(unique(emp4$year))

# leafout vs gdd
emp4$gddLeafout <- gddyr$GDD_5[match(emp4$yeardoyleafout, gddyr$yeardoy)]

plot(emp4$gddLeafout, emp4$leafout,
     xlab = "gddLeafout", ylab = "leafout",
     pch = 16, 
     col = yearcolors[match(emp4$year, years)],
     main = "leafout X gdd at leafout")

for (i in seq_along(years)) {
  year_dat <- emp4[emp4$year == years[i], ]
  
  lm_fit <- lm(leafout ~ gddLeafout, data = year_dat)
  x_seq  <- seq(min(year_dat$gddLeafout, na.rm = TRUE), 
                max(year_dat$gddLeafout, na.rm = TRUE), length.out = 200)
  pred   <- predict(lm_fit, newdata = data.frame(gddLeafout = x_seq))
  
  lines(x_seq, pred, 
        col = yearcolors[i],
        lwd = 2)
}

legend("bottomright",
       legend = years, 
       col= yearcolors, pch = 16, lty = 1, lwd = 2,
       title  = "Year")

# max temp MAM vs leafout 
emp4$tmeanmax_MAM <- climatesum$tmeanmax_MAM[match(emp4$year, climatesum$year)]

jpeg(
  filename = "figures/climate/leafoutXtmeanmax.jpeg", 
  width = 2400, height = 2400, res = 300)
plot(emp4$tmeanmax_MAM, emp4$leafout,
     xlab = "tmeanmax_MAM", ylab = "leafout",
     pch = 16, frame = FALSE,
     col = yearcolors[match(emp4$year, years)],
     main = "leafout X temperature mean max for March/April/May")

lm_fit <- lm(leafout ~ tmeanmax_MAM, data = emp4)
x_seq  <- seq(min(emp4$tmeanmax_MAM, na.rm = TRUE), 
              max(emp4$tmeanmax_MAM, na.rm = TRUE), length.out = 200)
pred   <- predict(lm_fit, newdata = data.frame(tmeanmax_MAM = x_seq))

lines(x_seq, pred, 
      col = "black",
      lwd = 2)

legend("bottomright",
       legend = years, 
       col= yearcolors, pch = 16, lty = 1, lwd = 2,
       title  = "Year")
dev.off()

}

# min temp MAM vs leafout 
emp4$tmeanmin_MAM <- climatesum$tmeanmin_MAM[match(emp4$year, climatesum$year)]

plot(emp4$tmeanmin_MAM, emp4$leafout,
     xlab = "tmeanmin_MAM", ylab = "leafout",
     pch = 16, 
     col = yearcolors[match(emp4$year, years)],
     main = "leafout X gdd at leafout")

lm_fit <- lm(leafout ~ tmeanmin_MAM, data = emp4)
x_seq  <- seq(min(emp4$tmeanmin_MAM, na.rm = TRUE), 
              max(emp4$tmeanmin_MAM, na.rm = TRUE), length.out = 200)
pred   <- predict(lm_fit, newdata = data.frame(tmeanmin_MAM = x_seq))

lines(x_seq, pred, 
      col = "black",
      lwd = 2)

# number of frost free days before leafout
weldhillclim$frostFreeDays <- ave(weldhillclim$minTempC, weldhillclim$year, FUN = function(x)
{sapply(seq_along(x), function(i) sum(x[seq_len(i-1)] < 0, na.rm = TRUE))
})

weldhillclim$yeardoy <- paste(weldhillclim$year, gddyr$doy, sep = "_")
emp5 <- emp4[!is.na(emp4$budburst),]
emp5$frostFbudburst <- weldhillclim$frostFreeDays[match(emp5$yeardoybudburst, 
                                                        weldhillclim$yeardoy)]
unique(emp5$budburst)
plot(emp5$frostFbudburst, emp5$budburst,
     xlab = "number of frost free days at budburst", ylab = "budburst",
     pch = 16, 
     col = yearcolors[match(emp5$year, years)],
     main = "leafout X frost free days at budburst")

for (i in seq_along(years)) {
  year_dat <- emp5[emp5$year == years[i], ]  # <-- years[i] not i
  
  lm_fit <- lm(budburst ~ frostFbudburst, data = year_dat)
  x_seq  <- seq(min(year_dat$frostFbudburst, na.rm = TRUE), 
                max(year_dat$frostFbudburst, na.rm = TRUE), length.out = 200)
  pred   <- predict(lm_fit, newdata = data.frame(frostFbudburst = x_seq))
  
  lines(x_seq, pred, 
        col = yearcolors[i],
        lwd = 2)
}

# precipitation at leafout
emp4$ppmMM <- weldhillclim$ppm[match(emp4$yeardoyleafout, weldhillclim$yeardoy)]

jpeg(
  filename = "figures/climate/leafoutXfrostFdays.jpeg", 
  width = 2400, height = 2400, res = 300
)
plot(emp4$budburst, emp4$leafout,
     xlab = "number of frost free days at budburst", ylab = "leafout",
     pch = 16, frame = FALSE,
     col = yearcolors[match(emp4$year, years)],
     main = "leafout X frost free days at budburst")

for (i in seq_along(years)) {
  year_dat <- emp4[emp4$year == years[i], ]
  
  lm_fit <- lm(leafout ~ budburst, data = year_dat)
  x_seq  <- seq(min(year_dat$budburst, na.rm = TRUE), 
                max(year_dat$budburst, na.rm = TRUE), length.out = 200)
  pred   <- predict(lm_fit, newdata = data.frame(budburst = x_seq))
  
  lines(x_seq, pred, 
        col = yearcolors[i],
        lwd = 2)
  legend("bottomright",
         legend = years, 
         col= yearcolors, pch = 16, lty = 1, lwd = 2,
         title  = "Year")
}
dev.off()


# 
# ggplot(emp) +
#   geom_point(aes(x = gddLeafout, y = leafout)) +
#   facet_wrap(~year) + theme_minimal()
# ggsave("figures/climate/leafoutGDD.jpeg", width = 8, height = 6, units = "in", dpi = 300)
# 
# # budset vs pdsi mam
# emp$pdsiMAM <- climatesum$pdsi_MAM[match(emp$year, climatesum$year)]
# ggplot(emp) +
#   geom_point(aes(x = pdsiMAM, y = budset, color = year)) +
#   # facet_wrap(~year) + 
#   theme_minimal()
# 
# # budset vs pdsi jja
# emp4$pdsiJJA <- climatesum$pdsi_JJA[match(emp4$year, climatesum$year)]
# 
# plot(emp4$pdsiJJA, emp4$budset,
#      xlab = "pdsiJJA", ylab = "budset",
#      pch = 16, 
#      col = yearcolors[match(emp4$year, years)],
#      main = "budset X summer PDSI")
# 
# lm_fit <- lm(budset ~ pdsiJJA, data = emp4)
# x_seq  <- seq(min(emp4$pdsiJJA, na.rm = TRUE), 
#               max(emp4$pdsiJJA, na.rm = TRUE), length.out = 200)
# pred   <- predict(lm_fit, newdata = data.frame(pdsiJJA = x_seq))
# 
# lines(x_seq, pred, 
#       col = "black",
#       lwd = 2)
# 
# 
# emp$sppyear <- paste(emp$spp, emp$year, sep = "_")
# emp$lengthMM <- emp$lengthCM*10
# empclim <- aggregate(lengthMM ~ sppyear, emp, mean)
# q25 <- aggregate(lengthMM ~ sppyear, emp, function(x) quantile(x, 0.25))
# q75 <- aggregate(lengthMM ~ sppyear, emp, function(x) quantile(x, 0.75))
# 
# empclim$q25 <- q25$lengthMM[match(empclim$sppyear, q25$sppyear)]
# empclim$q75 <- q75$lengthMM[match(empclim$sppyear, q75$sppyear)]
# 
# empclim$spp <- substr(empclim$sppyear, 1,6)
# empclim$year <- substr(empclim$sppyear, 8,11)
# 
# n_spp <- length(unique(empclim$spp))
# n_year <- length(unique(empclim$year))
# y_pos <- 1:n_year 
# 
# 
# empclim$pdsimam <- climatesum$pdsi_MAM[match(empclim$year, climatesum$year)]
# 
# pal <- colorRampPalette(MetBrewer::MetPalettes$VanGogh3[[1]])
# year_cols <- setNames(pal(n_year), unique(empclim$year))
# 
# year_pdsi <- unique(empclim[, c("year", "pdsimam")])
# year_pdsi <- year_pdsi[order(year_pdsi$pdsimam), ]
# 
# green_pal <- carto.pal(pal1 = "blue.pal", n1 = n_year)
# year_cols <- setNames(green_pal, year_pdsi$year)


if (makeplots){
  jpeg(
    filename = "figures/climate/rwPDSI.jpeg",
    width = 3600,      # wider image (pixels) → more horizontal room
    height = 2400,
    res = 300          # good print-quality resolution
  )
  par(mfrow = c(1, 4))
  y_pos <- rev(1:n_year)
  
  for (i in unique(empclim$spp)) { # i = "ALNINC"
    sub <- empclim[empclim$spp == i, ]
    
    plot(sub$lengthMM, y_pos,
         xlim = range(c(sub$q25-0.5, sub$q75+0.5)),
         ylim = c(0.5, n_year + 0.5),
         xlab = "Ring width (mm)",
         ylab = "",
         yaxt = "n",
         pch = 16,
         cex = 2,
         col = year_cols[sub$year],
         frame.plot = FALSE,
         main = i)
    
    axis(2, at = y_pos, labels = sub$year, las = 2, tick = FALSE, cex.axis = 1)
    
    abline(v = mean(sub$lengthMM), lty = 2)
    # error bars and dashed line
    segments(sub$q25, y_pos,
             sub$q75, y_pos,
             col = year_cols[sub$year],
             lwd = 3)
  }
  legend("bottomright",
         legend = paste(year_pdsi$year, round(year_pdsi$pdsimam, 2), sep = " PDSI: "),
         col    = year_cols[year_pdsi$year],
         pch    = 16,
         bty    = "n",
         cex    = 1,
         title  = "Year (PDSI MarchAprilMay)")
  dev.off()
}
