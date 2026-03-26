# Wildchrokie (and eventually coringTreespotters) climate data with phenology
# CRD 12 March 2026

# housekeeping
# rm(list=ls())
# options(stringsAsFactors = FALSE)
options(max.print = 150)
options(digits = 3)

# Load library 
library(ggplot2)
library(rstan)
library(future)
library(wesanderson)
library(patchwork) 
library(lme4)
library(lmerTest)  

if (length(grep("christophe_rouleau-desrochers", getwd())) > 0) {
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} else if (length(grep("lizzie", getwd())) > 0) {
  setwd("/Users/lizzie/Documents/git/projects/others/christophe/wildchrokie/analyses")
} else  {
  setwd("/home/crouleau/wildchrokie/analyses")
}

# flags
makeplots <- FALSE

empir <- read.csv("output/empiricalDataMAIN.csv")
climatesum <- read.csv("output/climateSummariesYear.csv")
gddyr <- read.csv("output/gddByYear.csv")
weldhillclim <- read.csv("output/weldhillClimateCleaned.csv")

commonNames <- c(
  "Alnus incana"          = "Grey alder",
  "Betula alleghaniensis" = "Yellow birch",
  "Betula papyrifera"     = "Paper birch",
  "Betula populifolia"    = "Gray birch"
)

empir$commonName <- commonNames[empir$latbi]

empir <- empir[!is.na(empir$pgsGDD5),]
# transform my groups to numeric values
empir$site_num <- match(empir$site, unique(empir$site))
empir$spp_num <- match(empir$spp, unique(empir$spp))
empir$treeid_num <- match(empir$treeid, unique(empir$treeid))
empir$year_num <- match(empir$year, unique(empir$year))

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate data #### 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
emp4 <- empir[!is.na(empir$leafout),]
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
years      <- sort(unique(empir$year))
firststeps <- colorRampPalette(c("#9cc184", "#192813"))(length(years))
empir$anomleafout <- empir$leafout - mean(empir$leafout)
empir$anombudset <- empir$budset - mean(empir$budset)

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

if (makeplots) {
  
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

}
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Climate summaries ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# common objects across budset and leafout
clim_vars  <- c("PDSI", 
                "TempMeanMax", 
                "TempMeanMean", 
                "TempMeanMin", 
                "Precip")

emp_clim <- merge(empir, climatesum, by = "year", all.x = TRUE)

colnames(emp_clim)[which(colnames(emp_clim) %in% "pdsi")] <- "PDSI"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmax")] <- "TempMeanMax"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmean")] <- "TempMeanMean"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmin")] <- "TempMeanMin"
colnames(emp_clim)[which(colnames(emp_clim) %in% "ppt")] <- "Precip"

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Leafout ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

  # Initialize results data frame
  climresultsleafout <- data.frame(
    predictor   = character(),
    period      = character(),
    slope       = numeric(),
    std_error   = numeric(),
    t_value     = numeric(),
    p_value     = numeric(),
    stringsAsFactors = FALSE
  )
  
  jpeg(
    filename = "figures/climate/climSumLeafout.jpeg", 
    width = 2400, height = 3600, res = 300)
  periods    <- c("DJF", "MAM")
  par(mfrow = c(5, 2), 
      mar = c(4, 4, 2, 1),   
      oma = c(0, 0, 4, 8))
  for (i in seq_along(clim_vars)) {
    for (j in seq_along(periods)) {
      
      p   <- periods[j]
      var <- clim_vars[i]
      
      dat <- emp_clim[emp_clim$period == p & !is.na(emp_clim[[var]]) & 
                        !is.na(emp_clim$anomleafout), ]
      
      plot(dat[[var]], dat$anomleafout,
           xlab = var, ylab  = "leafout",
           ylim = c(min(empir$anomleafout), max(empir$anomleafout)),
           pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
           main = "")
      
      if (i == 1) {
        mtext(p, side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2)
      }
      
      if (nrow(dat) > 1) {
        tmp    <- data.frame(x = dat[[var]], y = dat$anomleafout, year = dat$year)
        lm_fit <- lmer(y ~ scale(x) + (1 | year), data = tmp)
        sum <- summary(lm_fit)
        significance <- ifelse(sum$coefficients[2,4] < 0.05, "signif", "nonsignif")
        x_seq  <- seq(min(tmp$x, na.rm = TRUE), max(tmp$x, na.rm = TRUE), 
                      length.out = 200)
        pred <- predict(lm_fit, newdata = data.frame(x = x_seq, year = NA), re.form = NA)
        lines(x_seq, pred, col = "black", lwd = 2)
        slope <- round(fixef(lm_fit)[2], 2)
        mtext(paste0("β = ", slope), 
              side = 3, line = -2, cex = 1, adj = 0.8)
        if (significance == "signif") {
          mtext(" *", side = 3, line = -2, cex = 2, adj = 0.95)
        }
        
        # Collect results
        climresultsleafout <- rbind(climresultsleafout, data.frame(
          predictor    = var,
          period       = p,
          slope     = sum$coefficients[2, 1],
          std_error    = sum$coefficients[2, 2],
          t_value      = sum$coefficients[2, 4],
          p_value      = sum$coefficients[2, 5],
          stringsAsFactors = FALSE
        ))
        
      } else {
        # Still record the row but with NAs when insufficient data
        climresultsleafout <- rbind(climresultsleafout, data.frame(
          predictor    = var,
          period       = p,
          slope        = NA_real_,
          std_error    = NA_real_,
          t_value      = NA_real_,
          p_value      = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  par(xpd = NA)
  legend(x = par("usr")[2] + 2, y = mean(par("usr")[3:4]),
         legend = years, col = firststeps, pch = 16, lty = 1, lwd = 2,
         title = "Year", bty = "y", xjust = 0, yjust = -7)
  dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Budset ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
  climresultsbudset <- data.frame(
    predictor   = character(),
    period      = character(),
    slope       = numeric(),
    std_error   = numeric(),
    t_value     = numeric(),
    p_value     = numeric(),
    stringsAsFactors = FALSE
  )
  
  jpeg(
  filename = "figures/climate/climSumBudset.jpeg", 
  width = 2400, height = 3600, res = 300)

periods    <- c("MAM", "JJA", "SON")

par(mfrow = c(5, 3), 
    mar = c(4, 4, 2, 1),   
    oma = c(0, 0, 4, 8))

for (i in seq_along(clim_vars)) { # i = "tmeanmin"
  for (j in seq_along(periods)) { # j = "MAM"
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_clim[emp_clim$period == p & !is.na(emp_clim[[var]]) & 
                      !is.na(emp_clim$anombudset), ]
    
    plot(dat[[var]], dat$anombudset,
         xlab = var, ylab  = "budset",
         ylim = c(min(empir$anombudset), max(empir$anombudset)),
         pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
         main = "")
    
    abline(h = 0, lty = 2, col = "gray50")
    
    # Add column headers 
    if (i == 1) {
      mtext(p, side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2)
    }
    
    
    if (nrow(dat) > 1) {
      tmp    <- data.frame(x = dat[[var]], y = dat$anombudset, year = dat$year)
      lm_fit <- lmer(y ~ scale(x) + (1 | year), data = tmp)
      sum <- summary(lm_fit)
      significance <- ifelse(sum$coefficients[2,4]<0.05, "signif", "nonsignif")
      x_seq  <- seq(min(tmp$x, na.rm = TRUE), max(tmp$x, na.rm = TRUE), 
                    length.out = 200)
      pred <- predict(lm_fit, newdata = data.frame(x = x_seq, year = NA), re.form = NA)
      lines(x_seq, pred, col = "black", lwd = 2)
      slope <- round(fixef(lm_fit)[2], 2)
      mtext(paste0("β = ", slope), 
            side = 3, line = -2, cex = 1, adj = 0.8)
      if (significance == "signif") {
        mtext(" *", 
              side = 3, line = -2, cex = 2, adj = 0.95)
      }
      
      # Collect results
      climresultsbudset <- rbind(climresultsbudset, data.frame(
        predictor    = var,
        period       = p,
        slope     = sum$coefficients[2, 1],
        std_error    = sum$coefficients[2, 2],
        t_value      = sum$coefficients[2, 4],
        p_value      = sum$coefficients[2, 5],
        stringsAsFactors = FALSE
      ))
      
    } else {
      # Still record the row but with NAs when insufficient data
      climresultsbudset <- rbind(climresultsbudset, data.frame(
        predictor    = var,
        period       = p,
        slope        = NA_real_,
        std_error    = NA_real_,
        t_value      = NA_real_,
        p_value      = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Legend in outer right margin
par(xpd = NA)
legend(x = par("usr")[2] + 2, y = mean(par("usr")[3:4]),
       legend = years, col = firststeps, pch = 16, lty = 1, lwd = 2,
       title = "Year", bty = "y", xjust = 0, yjust = -7)
dev.off()

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# CoringTreespotters ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
empts <- read.csv("/Users/christophe_rouleau-desrochers/github/coringTreespotters/analyses/output/empiricalDataMAIN.csv")

empts$anomleafcolor <- empts$coloredLeaves - mean(empts$coloredLeaves)

emp_climts <- merge(empts, climatesum, by = "year", all.y = TRUE)

colnames(emp_climts)[which(colnames(emp_climts) %in% "pdsi")] <- "PDSI"
colnames(emp_climts)[which(colnames(emp_climts) %in% "tmeanmax")] <- "TempMeanMax"
colnames(emp_climts)[which(colnames(emp_climts) %in% "tmeanmean")] <- "TempMeanMean"
colnames(emp_climts)[which(colnames(emp_climts) %in% "tmeanmin")] <- "TempMeanMin"
colnames(emp_climts)[which(colnames(emp_climts) %in% "ppt")] <- "Precip"

years <- sort(unique(empts$year))
firststeps <- colorRampPalette(c("#9cc184", "#192813"))(length(years))

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Leafout ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
emp_climtslo <- emp_climts[!is.na(emp_climts$leafout),]
emp_climtslo$anomleafout <- emp_climtslo$leafout - mean(emp_climtslo$leafout)
jpeg("figures/climate/climSumLeafoutTS.jpeg",
     width = 2400, height = 3600, res = 300)
periods <- c("DJF", "MAM") 
par(mfrow = c(5, 2), 
    mar = c(4, 4, 2, 1),   
    oma = c(0, 0, 4, 8))

for (i in seq_along(clim_vars)) {
  for (j in seq_along(periods)) {
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_climtslo[emp_climtslo$period == p & !is.na(emp_climtslo[[var]]) & 
                        !is.na(emp_climtslo$anomleafout), ]
    
    plot(dat[[var]], dat$anomleafout,
         xlab = var, ylab  = "leafout",
         ylim = c(min(emp_climtslo$anomleafout), max(emp_climtslo$anomleafout)),
         pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
         main = "")
    
    if (i == 1) {
      mtext(p, side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2)
    }
    
    if (nrow(dat) > 1) {
      tmp    <- data.frame(x = dat[[var]], y = dat$anomleafout, year = dat$year)
      lm_fit <- lmer(y ~ x + (1 | year), data = tmp)
      sum <- summary(lm_fit)
      significance <- ifelse(sum$coefficients[2,4] < 0.05, "signif", "nonsignif")
      x_seq  <- seq(min(tmp$x, na.rm = TRUE), max(tmp$x, na.rm = TRUE), length.out = 200)
      pred <- predict(lm_fit, newdata = data.frame(x = x_seq, year = NA), re.form = NA)
      lines(x_seq, pred, col = "black", lwd = 2)
      slope <- round(fixef(lm_fit)[2], 2)
      mtext(paste0("β = ", slope), side = 3, line = -2, cex = 1, adj = 0.8)
      if (significance == "signif") {
        mtext(" *", side = 3, line = -2, cex = 2, adj = 0.95)
      }
    }
  }
}
par(xpd = NA)
legend(x = 0, y = 0,
       legend = years, col = firststeps, pch = 16, lty = 1, lwd = 2,
       title = "Year", bty = "y", xjust = -9, yjust = -3)

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### coloredLeaves ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
emp_climtscl <- emp_climts[!is.na(emp_climts$coloredLeaves),]
emp_climtscl$anomleafcolor <- emp_climtscl$coloredLeaves - 
  mean(emp_climtscl$coloredLeaves)

jpeg(
  filename = "figures/climate/climSumLeafColTS.jpeg", 
  width = 2400, height = 3600, res = 300)

periods    <- c("MAM", "JJA", "SON")

par(mfrow = c(5, 3), 
    mar = c(4, 4, 2, 1),   
    oma = c(0, 0, 4, 8))

for (i in seq_along(clim_vars)) { # i = "tmeanmin"
  for (j in seq_along(periods)) { # j = "MAM"
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_climtscl[emp_climtscl$period == p & !is.na(emp_climtscl[[var]]) & 
                          !is.na(emp_climtscl$anomleafcolor), ]
    
    plot(dat[[var]], dat$anomleafcolor,
         xlab = var, ylab  = "leaf color",
         ylim = c(min(emp_climtscl$anomleafcolor), 
                  max(emp_climtscl$anomleafcolor)+20),
         pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
         main = "")
    
    abline(h = 0, lty = 2, col = "gray50")
    
    # Add column headers 
    if (i == 1) {
      mtext(p, side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2)
    }
    
    
    if (nrow(dat) > 1) {
      tmp    <- data.frame(x = dat[[var]], y = dat$coloredLeaves, year = dat$year)
      lm_fit <- lmer(y ~ scale(x) + (1 | year), data = tmp)
      sum <- summary(lm_fit)
      significance <- ifelse(sum$coefficients[2,4]<0.05, "signif", "nonsignif")
      x_seq  <- seq(min(tmp$x, na.rm = TRUE), max(tmp$x, na.rm = TRUE), 
                    length.out = 200)
      pred <- predict(lm_fit, newdata = data.frame(x = x_seq, year = NA), re.form = NA)
      lines(x_seq, pred, col = "black", lwd = 2)
      slope <- round(fixef(lm_fit)[2], 2)
      mtext(paste0("β = ", slope), 
            side = 3, line = -2, cex = 1, adj = 0.8)
      if (significance == "signif") {
        mtext(" *", 
              side = 3, line = -2, cex = 2, adj = 0.95)
      }
    }
  }
}

# Legend in outer right margin
par(xpd = NA)
legend(x = par("usr")[2] + 2, y = mean(par("usr")[3:4]),
       legend = years, col = firststeps, pch = 16, lty = 1, lwd = 2,
       title = "Year", bty = "y", xjust = 0, yjust = -7)

dev.off()



# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Phenology ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
leafoutbyyr <- aggregate(leafout ~ year + latbi, emp4, FUN = mean)
leafoutbyyr$leafout <- round(leafoutbyyr$leafout, 2)
budsetbyyr <- aggregate(budset ~ year + latbi, emp4, FUN = mean)
budsetbyyr$budset <- round(budsetbyyr$budset, 2)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Fit the figures with stan ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
clim_vars <- c("TempMeanMax", "TempMeanMin", "TempMeanMean")
climvar <- clim_vars[1] 
period <- "MAM"

d <- emp_clim[emp_clim$period == period & !is.na(emp_clim$TempMeanMax) & 
                !is.na(emp_clim$anombudset), ]


# transform data in vectors for gsl
data <- list(
  y = d$leafout / 10,
  N = nrow(d),
  year = as.numeric(as.character(d$year_num)),
  species = as.numeric(as.character(d$spp_num)),
  site = as.numeric(as.character(d$site_num)),
  Nspp = length(unique(d$spp_num)),
  Nsite = length(unique(d$site_num)),
  Nyear = length(unique(d$year_num)),
  climpredictor = d$TempMeanMax
)

rstan_options(auto_write = TRUE)
climmodel <- stan_model("stan/climatePredictors.stan")
fit <- sampling(climmodel, data = data, 
                       warmup = 1000, iter = 2000, chains=4)

post_means <- summary(fit)$summary[, "mean"]

# Pull out what you need by name
a <- post_means["a"]
aspp <- post_means[grep("^aspp", names(post_means))]
asite <- post_means[grep("^asite", names(post_means))]
ayear <- post_means[grep("^ayear", names(post_means))]
bsp <- post_means[grep("^bsp",  names(post_means))]

x_vals <- unique(d$TempMeanMax)

# Set up empty plot
plot(NULL, xlim = range(x_vals), ylim = c(min(d$leafout/10), max(d$leafout/10)), 
     xlab = "x", ylab = "y", xaxt = "n")

cols <- c("firebrick", "steelblue", "forestgreen", "darkorange")

for (s in 1:4) { # i = 2
  intercept_s <- a + aspp[s]
  slope_s <- bsp[s]
  y_vals <- intercept_s + slope_s * x_vals
  lines(x_vals, y_vals, col = cols[s], lwd = 2)
  points(x_vals, y_vals, col = cols[s], pch = 16, cex = 1.2)
}

legend("topleft", legend = paste("spp", 1:4),
       col = cols, lwd = 2, pch = 16)


# add asite/ayear if needed for fitted values

# Example: plot per species (assuming x = some continuous predictor)
# bsp[1] = slope for spp 1, aspp[1] + a = intercept for spp 1

x_range <- seq(min(yourdata$x), max(yourdata$x), length.out = 100)

plot(yourdata$x, yourdata$y, col = yourdata$species, pch = 16,
     xlab = "x", ylab = "y")

for (i in seq_along(clim_vars)) { # i = "tmeanmin"
  for (j in seq_along(periods)) { # j = "MAM"
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_clim[emp_clim$period == p & !is.na(emp_clim[[var]]) & 
                      !is.na(emp_clim$anombudset), ]
    
    plot(dat[[var]], dat$anombudset,
         xlab = var, ylab  = "budset",
         ylim = c(min(empir$anombudset), max(empir$anombudset)),
         pch = 16, frame = FALSE, col = firststeps[match(dat$year, years)],
         main = "")
    
    abline(h = 0, lty = 2, col = "gray50")
  }
}