# Wildchrokie (and eventually coringTreespotters) climate data with phenology
# CRD 12 March 2026

# housekeeping
# rm(list=ls())
# options(stringsAsFactors = FALSE)
# options(max.print = 150)
# options(digits = 3)

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

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
years      <- sort(unique(emp$year))
firststeps <- colorRampPalette(c("#9cc184", "#192813"))(length(years))
emp$anomleafout <- emp$leafout - mean(emp$leafout)
emp$anombudset <- emp$budset - mean(emp$budset)

# leafout vs gdd
emp4$gddLeafout <- gddyr$GDD_5[match(emp4$yeardoyleafout, gddyr$yeardoy)]
if (makeplots){
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
     pch = 16, frame = FALSE,
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
  legend("bottomright",
         legend = years, 
         col= yearcolors, pch = 16, lty = 1, lwd = 2,
         title  = "Year")
}

# precipitation at leafout
emp4$winterPptLeafout <- mapply(function(leafout_doy, obs_year) {
  # takes the previous year accumulation of ppt in december
  sub <- weldhillclim[(weldhillclim$year == obs_year - 1 & weldhillclim$doy >= 335) |
                        # then going into the current year condition
                        (weldhillclim$year == obs_year & weldhillclim$doy <= leafout_doy), ]
  sum(sub$pptMM, na.rm = TRUE) # sum the ppt over our period of interest
}, emp4$leafout, emp4$year) # apply the function to each of those 2 arguments

plot(emp4$winterPptLeafout, emp4$leafout,
     xlab = "precipitation accumulation (mm) at leafout", ylab = "leafout",
     pch = 16, frame = FALSE,
     col = yearcolors[match(emp4$year, years)],
     main = "leafout X precipitation accumulation (mm) at leafout")

for (i in seq_along(years)) { # i = 2018
  year_dat <- emp4[emp4$year == years[i], ]
  
  lm_fit <- lm(leafout ~ winterPptLeafout, data = year_dat)
  x_seq  <- seq(min(year_dat$winterPptLeafout, na.rm = TRUE), 
                max(year_dat$winterPptLeafout, na.rm = TRUE), length.out = 200)
  pred   <- predict(lm_fit, newdata = data.frame(winterPptLeafout = x_seq))
  
  lines(x_seq, pred, 
        col = yearcolors[i],
        lwd = 2)
  legend("bottomright",
         legend = years, 
         col= yearcolors, pch = 16, lty = 1, lwd = 2,
         title  = "Year")
}

}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate summaries ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# common objects across budset and leafout
clim_vars  <- c("PDSI", 
                "TempMeanMax", 
                "TempMeanMean", 
                "TempMeanMin", 
                "Precip"      , 
                "Radiation")
emp_clim <- merge(emp, climatesum, by = "year", all.x = TRUE)

colnames(emp_clim)[which(colnames(emp_clim) %in% "pdsi")] <- "PDSI"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmax")] <- "TempMeanMax"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmean")] <- "TempMeanMean"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmin")] <- "TempMeanMin"
colnames(emp_clim)[which(colnames(emp_clim) %in% "ppt")] <- "Precip"
colnames(emp_clim)[which(colnames(emp_clim) %in% "rad")] <- "Radiation"

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Leafout ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
if (makeplots) {

jpeg(
  filename = "figures/climate/climSumLeafout.jpeg", 
  width = 2400, height = 3600, res = 300)

periods    <- c("DJF", "MAM")

par(mfrow = c(6, 2), 
    mar = c(4, 4, 2, 1),   
    oma = c(0, 0, 4, 8))

for (i in seq_along(clim_vars)) { # i = "tmeanmin"
  for (j in seq_along(periods)) { # j = "MAM"
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_clim[emp_clim$period == p & !is.na(emp_clim[[var]]) & 
                      !is.na(emp_clim$anomleafout), ]
    
    plot(dat[[var]], dat$anomleafout,
         xlab = var, ylab  = "leafout",
         ylim = c(min(emp$anomleafout), max(emp$anomleafout)),
         pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
         main = "")
    
    abline(h = 0, lty = 2, col = "gray50")
    
    # Add column headers 
    if (i == 1) {
      mtext(p, side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2)
    }
    
    if (nrow(dat) > 1) {
      tmp    <- data.frame(x = dat[[var]], y = dat$anomleafout)
      lm_fit <- lm(y ~ scale(x), data = tmp)
      sum <- summary(lm_fit)
      significance <- ifelse(sum$coefficients[2,4]<0.05, "signif", "nonsignif")
      x_seq  <- seq(min(tmp$x, na.rm = TRUE), max(tmp$x, na.rm = TRUE), 
                    length.out = 200)
      pred   <- predict(lm_fit, newdata = data.frame(x = x_seq))
      lines(x_seq, pred, col = "black", lwd = 2)
      slope <- round(coef(lm_fit)[2], 2)
      mtext(paste0("β = ", slope, ifelse(significance == "signif", " *", "")), 
            side = 3, line = -2, cex = 0.7, adj = 0.95)
    }
  }
}


# Legend in outer right margin
par(xpd = NA)
legend(x = par("usr")[2] + 2, y = mean(par("usr")[3:4]),
       legend = years, col = firststeps, pch = 16, lty = 1, lwd = 2,
       title = "Year", bty = "y", xjust = 0, yjust = -7)
dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Budset ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
jpeg(
  filename = "figures/climate/climSumBudset.jpeg", 
  width = 2400, height = 3600, res = 300)

periods    <- c("MAM", "JJA", "SON")

par(mfrow = c(6, 3), 
    mar = c(4, 4, 2, 1),   
    oma = c(0, 0, 4, 8))

# test linear model
test <- subset(emp_clim, period %in% "MAM")
lm(anomleafout ~ TempMeanMax + spp + year, data = test)

for (i in seq_along(clim_vars)) { # i = "tmeanmin"
  for (j in seq_along(periods)) { # j = "MAM"
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_clim[emp_clim$period == p & !is.na(emp_clim[[var]]) & 
                      !is.na(emp_clim$anombudset), ]
    
    plot(dat[[var]], dat$anombudset,
         xlab = var, ylab  = "budset",
         ylim = c(min(emp$anombudset), max(emp$anombudset)),
         pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
         main = "")
    
    abline(h = 0, lty = 2, col = "gray50")
    
    # Add column headers 
    if (i == 1) {
      mtext(p, side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2)
    }
    
    if (nrow(dat) > 1) {
      tmp    <- data.frame(x = dat[[var]], y = dat$anombudset)
      lm_fit <- lm(y ~ scale(x), data = tmp)
      sum <- summary(lm_fit)
      significance <- ifelse(sum$coefficients[2,4]<0.001, "signif", "nonsignif")
      x_seq  <- seq(min(tmp$x, na.rm = TRUE), max(tmp$x, na.rm = TRUE), 
                    length.out = 200)
      pred   <- predict(lm_fit, newdata = data.frame(x = x_seq))
      lines(x_seq, pred, col = "black", lwd = 2)
      slope <- round(coef(lm_fit)[2], 2)
      mtext(paste0("β = ", slope, ifelse(significance == "signif", " *", "")), 
            side = 3, line = -2, cex = 0.7, adj = 0.95)
    }
  }
}

# Legend in outer right margin
par(xpd = NA)
legend(x = par("usr")[2] + 2, y = mean(par("usr")[3:4]),
       legend = years, col = firststeps, pch = 16, lty = 1, lwd = 2,
       title = "Year", bty = "y", xjust = 0, yjust = -7)
dev.off()


# if (makeplots){
#   n_year <- length(unique(emp$year))
#   jpeg(
#     filename = "figures/climate/rwPDSI.jpeg",
#     width = 3600,      # wider image (pixels) → more horizontal room
#     height = 2400,
#     res = 300          # good print-quality resolution
#   )
#   par(mfrow = c(1, 4))
#   y_pos <- rev(1:n_year)
#   
#   for (i in unique(empclim$spp)) { # i = "ALNINC"
#     sub <- empclim[empclim$spp == i, ]
#     
#     plot(sub$lengthMM, y_pos,
#          xlim = range(c(sub$q25-0.5, sub$q75+0.5)),
#          ylim = c(0.5, n_year + 0.5),
#          xlab = "Ring width (mm)",
#          ylab = "",
#          yaxt = "n",
#          pch = 16,
#          cex = 2,
#          col = year_cols[sub$year],
#          frame.plot = FALSE,
#          main = i)
#     
#     axis(2, at = y_pos, labels = sub$year, las = 2, tick = FALSE, cex.axis = 1)
#     
#     abline(v = mean(sub$lengthMM), lty = 2)
#     # error bars and dashed line
#     segments(sub$q25, y_pos,
#              sub$q75, y_pos,
#              col = year_cols[sub$year],
#              lwd = 3)
#   }
#   legend("bottomright",
#          legend = paste(year_pdsi$year, round(year_pdsi$pdsimam, 2), sep = " PDSI: "),
#          col    = year_cols[year_pdsi$year],
#          pch    = 16,
#          bty    = "n",
#          cex    = 1,
#          title  = "Year (PDSI MarchAprilMay)")
#   dev.off()
# }

}

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Phenology ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
leafoutbyyr <- aggregate(leafout ~ year + latbi, emp4, FUN = mean)
leafoutbyyr$leafout <- round(leafoutbyyr$leafout, 2)
budsetbyyr <- aggregate(budset ~ year + latbi, emp4, FUN = mean)
budsetbyyr$budset <- round(budsetbyyr$budset, 2)
