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
