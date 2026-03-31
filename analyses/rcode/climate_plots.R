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

# precipitation at leafout
if (makeplots) {
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
clim_vars  <- c("TempMeanMax", "TempMeanMean","TempMeanMin")

emp_clim <- merge(empir, climatesum, by = "year", all.x = TRUE)

colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmax")] <- "TempMeanMax"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmean")] <- "TempMeanMean"
colnames(emp_clim)[which(colnames(emp_clim) %in% "tmeanmin")] <- "TempMeanMin"

# define plot objects and stuff
species_order <- c(
  "Alnus incana", 
  "Betula alleghaniensis", 
  "Betula papyrifera", 
  "Betula populifolia")

my_colors <- c(
  "Alnus incana" = wes_palette("AsteroidCity1")[1],
  "Betula alleghaniensis" = wes_palette("AsteroidCity1")[2],
  "Betula papyrifera" = wes_palette("AsteroidCity1")[3],
  "Betula populifolia" = wes_palette("AsteroidCity1")[4]
  
)

years      <- sort(unique(empir$year))
firststeps <- colorRampPalette(c("#9cc184", "#192813"))(length(years))

climmodel <- stan_model("stan/climatePredictors.stan")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
##### Leafout ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
clim_vars <- c("TempMeanMin", "TempMeanMean", "TempMeanMax")
periods <- c("DJF", "MAM")

jpeg(
  filename = "figures/climate/climSumLeafout.jpeg", 
  width = 2000, height = 3000, res = 300)

layout(matrix(c(1,2,7,
                3,4,7,
                5,6,7), nrow = 3, byrow = TRUE),
       widths = c(2, 2, 1))

# set df to recover the parameters
all_params <- data.frame()

for (i in seq_along(clim_vars)) { # i = 2
  for (j in seq_along(periods)) { # j = 1
    
    p   <- periods[j]
    var <- clim_vars[i]
    
    dat <- emp_clim[emp_clim$period == p & !is.na(emp_clim[[var]]), ]
    
    # transform data in vectors for gsl
    data <- list(
      y = dat$anomleafout,
      N = nrow(dat),
      year = as.numeric(as.character(dat$year_num)),
      species = as.numeric(as.character(dat$spp_num)),
      site = as.numeric(as.character(dat$site_num)),
      Nspp = length(unique(dat$spp_num)),
      Nsite = length(unique(dat$site_num)),
      Nyear = length(unique(dat$year_num)),
      climpredictor = dat[[var]]
    )
    
    # Fit models
    rstan_options(auto_write = TRUE)
    fit <- sampling(climmodel, data = data, 
                    warmup = 1000, iter = 2000, chains=4, refresh = 0)
    
    post_means <- summary(fit)$summary[, "mean"]
    
    # Extract full summary with quantiles
    post_summary <- summary(fit, probs = c(0.05, 0.95))$summary
    
    # Build dataframe for all parameters of interest
    param_indices <- c(
      "a",
      grep("^ayear(?!.*_prior)", rownames(post_summary), value = TRUE, perl = TRUE),
      grep("^asite(?!.*_prior)", rownames(post_summary), value = TRUE, perl = TRUE),
      grep("^aspp(?!.*_prior)",  rownames(post_summary), value = TRUE, perl = TRUE),
      grep("^bsp(?!.*_prior)",   rownames(post_summary), value = TRUE, perl = TRUE)
    )
    param_df <- data.frame(
      clim_var  = var,
      period    = p,
      parameter = param_indices,
      mean      = post_summary[param_indices, "mean"],
      q5        = post_summary[param_indices, "5%"],
      q95       = post_summary[param_indices, "95%"],
      row.names = NULL
    )
    
    all_params <- rbind(all_params, param_df)
    # pull what I need
    a <- post_means["a"]
    aspp <- post_means[grep("^aspp", names(post_means))]
    asite <- post_means[grep("^asite", names(post_means))]
    ayear <- post_means[grep("^ayear", names(post_means))]
    bsp <- post_means[grep("^bsp",  names(post_means))]
    
    x_vals <- sort(unique(dat[[var]]))
    
    # Set up empty plot
    x_range <- range(x_vals)
    x_pad <- diff(x_range) * 0.08  # 8% padding on each side
    
    plot(NULL, 
         xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
         ylim = c(min(dat$anomleafout), max(dat$anomleafout)), 
         xlab = var, ylab = "leafout", 
         main = p, frame = FALSE)
    abline(h = 0, lty = 2, col = "gray50")
    
    for (s in 1:data$Nspp) { # i = 2
      intercept_s <- a + aspp[s]
      slope_s <- bsp[s]
      y_vals <- intercept_s + slope_s * x_vals
      lines(x_vals, y_vals, col = my_colors[s], lwd = 2)
    }
    
    points(data$climpredictor, data$y,
           col = my_colors[data$species], pch = 16, cex = 0.8)
    
    # one year label per unique clim value, at top of each plot
    text(x_vals, rep(max(dat$anomleafout), length(x_vals)),
         labels = dat$year[match(x_vals, dat[[var]])],
         cex = 1, col = "black")
  }
}
par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
legend("center", legend = species_order,
       col = my_colors[species_order], pch = 16, lwd = 2,
       bty = "n", cex = 0.9, pt.cex = 1.5)


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
if (FALSE){

clim_vars <- c("TempMeanMax", "TempMeanMin", "TempMeanMean")
climvar <- clim_vars[1] 
period <- "DJF"

d <- emp_clim[emp_clim$period == period & !is.na(emp_clim$TempMeanMax) & 
                !is.na(emp_clim$anombudset), ]

# transform data in vectors for gsl
data <- list(
  y = d$anomleafout,
  N = nrow(d),
  year = as.numeric(as.character(d$year_num)),
  species = as.numeric(as.character(d$spp_num)),
  site = as.numeric(as.character(d$site_num)),
  Nspp = length(unique(d$spp_num)),
  Nsite = length(unique(d$site_num)),
  Nyear = length(unique(d$year_num)),
  climpredictor = d$TempMeanMax
)

# Fit models
rstan_options(auto_write = TRUE)
climmodel <- stan_model("stan/climatePredictors.stan")
fit <- sampling(climmodel, data = data, iter = 2000, chains=4, cores = 4)

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# Plot posterior vs priors for gdd fit ####
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
pdf(file = "figures/climate/climateModelPriorVSPosterior.pdf", 
    width = 8, height = 10)

pal <- wes_palette("AsteroidCity1")[3:4]

par(mfrow = c(3, 2))
df_fit <- as.data.frame(fit)

columns <- colnames(df_fit)[!grepl("prior", colnames(df_fit))]
sigma_df <- df_fit[, columns[grepl("sigma", columns)]]
aspp_df <- df_fit[, columns[grepl("aspp", columns)]]
bspp_df <- df_fit[, columns[grepl("bsp", columns)]]
site_df <- df_fit[, columns[grepl("asite", columns)]]
ayear_df <- df_fit[, columns[grepl("ayear", columns)]]

# a
plot(density(df_fit[, "a_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_a", 
     xlab = "a", ylim = c(0,0.5))
lines(density(df_fit[, "a"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# sigma_y
plot(density(df_fit[, "sigma_y_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_sigma_y", 
     xlab = "sigma_y", ylim = c(0,2))
lines(density(df_fit[, "sigma_y"]), col = pal[2], lwd = 2)
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# aspp
plot(density(df_fit[, "aspp_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_aspp", 
     xlab = "aspp", xlim = c(-20, 20), ylim = c(0, 0.15))
for (col in colnames(aspp_df)) {
  lines(density(aspp_df[, col]), col = pal[2], lwd = 1)
} 
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# asite
plot(density(df_fit[, "asite_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_asite", 
     xlab = "asite", xlim = c(-20, 20), ylim = c(0, 0.5))
for (col in colnames(site_df)) {
  lines(density(site_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# ayear
plot(density(df_fit[, "ayear_prior"]), 
     col = pal[1], lwd = 2, 
     main = "priorVSposterior_ayear", 
     xlab = "ayear", xlim = c(-20, 20), ylim = c(0, 0.5))
for (col in colnames(ayear_df)) {
  lines(density(ayear_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

# bsp
plot(density(df_fit[, "bsp_prior"]), 
     col = pal[1], lwd = 2, xlim = c(-15, 15),
     main = "priorVSposterior_bsp", 
     xlab = "bsp", ylim = c(0, 1.8))
for (col in colnames(bspp_df)) {
  lines(density(bspp_df[, col]), col = pal[2], lwd = 1)
}
legend("topright", legend = c("Prior", "Posterior"), col = pal, lwd = 2)

dev.off()
}