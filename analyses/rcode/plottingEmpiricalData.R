# Wildchrokie model
# CRD 28 October April 2025
# Start plotting empirical data for wildchrokie

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 5)
# quartz()

# Load library 
library(ggplot2)
library(rstan)
library(rstanarm)
library(future)
library(shinystan)
library(wesanderson)
library(patchwork)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

sim <- read.csv("output/simdata.csv")
emp <- read.csv("output/empiricalDataMAIN.csv")
gdd <- read.csv("output/gddByYear.csv")

# add full species
emp$sppfull <- NA
emp$sppfull[which(emp$spp == "ALNINC")] <- "Alnus incana"
emp$sppfull[which(emp$spp == "BETPOP")] <- "Betula populifolia"
emp$sppfull[which(emp$spp == "BETPAP")] <- "Betula papyrifera"
emp$sppfull[which(emp$spp == "BETALL")] <- "Betula alleghaniensis"

# load fit
fit <- readRDS("output/stanOutput/fitEmpirical_stanlmer")
fit_pgsNgrowingdays <- readRDS("output/stanOutput/fit_pgsNgrowingdays_Empirical_stanlmer")

# faceted
# order by spp
emp <- emp[order(emp$spp), ]
emp$lengthMM <- emp$lengthCM*10

emp$year <- as.factor(emp$year)


# At eco evo:
ggplot(emp, aes(x = pgsGDD, y = lengthMM, 
                color = sppfull, 
                fill = sppfull)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  facet_wrap(~sppfull) +
  labs(y = "Ring width (mm)", x = "Growing degree days (GDD)", color = "Tree Species") +
  theme_bw() +
  theme(legend.key.height = unit(1.5, "lines"),
        strip.text = element_text(face = "bold.italic", size = 10)) +
  guides(fill = "none", color = "none") 

# new symbols and stuff
ggplot(emp, aes(x = pgsGDD, y = lengthMM)) +
  geom_point(size = 2, alpha = 0.9,
             aes(shape = prov,
                 color = year, 
                 fill = year)) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +  
  facet_wrap(~sppfull) +
  labs(y = "Ring width (mm)", 
       x = "Growing degree days (GDD)", 
       color = "Year",
       fill = "Year",
       shape = "Site") +  
  theme_bw() +
  guides(color = guide_legend(override.aes = list(shape = 21)),
         fill = guide_legend(override.aes = list(shape = 21)))
ggsave("figures/empiricalData/sppLinearRegressions_pgsGDD2.jpeg", width = 8, height = 6, units = "in", dpi = 300)




# full growing season
ggplot(emp, aes(x = fgsGDD, y = lengthCM, 
                color = sppfull, 
                fill = sppfull)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  facet_wrap(~sppfull) +
  labs(y = "Ring width (cm)", x = "Growing degree days (GDD)", color = "Tree Species") +
  theme_minimal() +
    theme(strip.text = element_blank(),         
          legend.key.height = unit(1.5, "lines")) +
  guides(fill = "none") 
ggsave("figures/empiricalData/sppLinearRegressions.jpeg", 
       width = 9, height = 6, units = "in", dpi = 300)

# just the number of days, without gdd
emp$pgsNgrowingdays <- emp$budset - emp$leafout
emp$fgsNgrowingdays <- emp$leafcolor - emp$budburst

# number of days  in pgs
ggplot(emp) +
  # geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(aes(x = pgsNgrowingdays, y = lengthCM, 
                  color = sppfull, 
                  fill = sppfull),
              method = "lm", se = TRUE, alpha = 0.2) +

  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  facet_wrap(~sppfull) +
  labs(y = "Ring width (cm)", x = "", color = "Tree Species") +
  theme_minimal() +
  theme(strip.text = element_blank(),         
          legend.key.height = unit(1.5, "lines")) +
  guides(fill = "none") 
ggsave("figures/empiricalData/sppLinearRegressions_pgsNdays.jpeg", width = 9, height = 6, units = "in", dpi = 300)
# number of days  in fgs
ggplot(emp, aes(x = fgsNgrowingdays, y = lengthCM, 
                color = sppfull, 
                fill = sppfull)) +
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  facet_wrap(~sppfull) +
  labs(y = "Ring width (cm)", x = "", color = "Tree Species") +
  theme_minimal() +   theme(strip.text = element_blank(),                    legend.key.height = unit(1.5, "lines")) +
  theme(strip.text = element_blank(),
        legend.key.height = unit(1.5, "lines")) +
  guides(fill = "none") 
# # plot gdd
# gddstats <- aggregate(GDD_10 ~ doy, FUN = mean, data = gdd)
# colnames(gddstats)[2] <-  "mean"
# gddstats2 <- aggregate(GDD_10 ~ doy, FUN = min, data = gdd)
# colnames(gddstats2)[2] <-  "min"
# gddstats3 <- aggregate(GDD_10 ~ doy, FUN = max, data = gdd)
# colnames(gddstats3)[2] <-  "max"
# 
# gddstats <- merge(gddstats, gddstats2, by = "doy")
# gddstats <- merge(gddstats, gddstats3, by = "doy")
# 
# ggplot(gddstats, aes(x = doy, y = mean)) + 
#   geom_line() +
#   geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2, color = NA) +
#   geom_vline(xintercept = 5) +
#   labs(title = "GDD accumulation")


# recover fit parameters for pgsGDD ####
df_fit <- as.data.frame(fit)

###### Recover treeid ######

# grab treeid 
treeid_cols <- colnames(df_fit)[grepl("treeid", colnames(df_fit))]
# remove sigma_asp for now
treeid_cols <- treeid_cols[1:length(treeid_cols)]

treeid_df <- df_fit[, colnames(df_fit) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub(".*treeid:([^]]+)\\].*", "\\1", colnames(treeid_df))

# empty treeid dataframe
treeid_df2 <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_a_treeid = numeric(ncol(treeid_df)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df2$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df2$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df2$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df2$fit_a_treeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df2$fit_a_treeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df2$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

if (FALSE) {
###### Recover a spp  ######
aspp_cols <- colnames(df_fit)[grepl("spp", colnames(df_fit))]
aspp_cols <- aspp_cols[!grepl("Sigma", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("pgsGDDscaled", aspp_cols)]

aspp_df <- df_fit[, colnames(df_fit) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub(".*spp:([^]]+)\\].*", "\\1", colnames(aspp_df))

#empty aspp df
aspp_df2 <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_a_spp = numeric(ncol(aspp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per25 = NA,
  fit_a_spp_per75 = NA,
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df2$spp[i] <- colnames(aspp_df)[i]         
  aspp_df2$fit_a_spp[i] <- round(mean(aspp_df[[i]]),3)  
  aspp_df2$fit_a_spp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df2$fit_a_spp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df2$fit_a_spp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df2$fit_a_spp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df2

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
}
###### Recover b spp ###### 
bspp_cols <- colnames(df_fit)[grepl("pgsGDDscaled", colnames(df_fit))]
# remove sigma_aspp for now
bspp_cols <- bspp_cols[!grepl("Sigma", bspp_cols)]
bspp_cols <- bspp_cols[2:5]


bspp_df <- df_fit[, colnames(df_fit) %in% bspp_cols]
# change their names
colnames(bspp_df) <- sub(".*spp:([^]]+)\\].*", "\\1", colnames(bspp_df))
#empty spp df
bspp_df2 <- data.frame(
  spp = character(ncol(bspp_df)),
  fit_b_spp = numeric(ncol(bspp_df)),
  fit_b_spp_per5 = NA,
  fit_b_spp_per25 = NA,
  fit_b_spp_per75 = NA,
  fit_b_spp_per95 = NA
)
for (i in 1:ncol(bspp_df)) { # i = 1
  bspp_df2$spp[i] <- colnames(bspp_df)[i]
  bspp_df2$fit_b_spp[i] <- round(mean(bspp_df[[i]]),3)
  bspp_df2$fit_b_spp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
  bspp_df2$fit_b_spp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
  bspp_df2$fit_b_spp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
  bspp_df2$fit_b_spp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
}
bspp_df2

###### Recover site ######
site_cols <- colnames(df_fit)[grepl("prov", colnames(df_fit))]
# remove for now
site_cols <- site_cols[1:length(site_cols)]
site_cols <- site_cols[!grepl("Sigma", site_cols)]
site_df <- df_fit[, colnames(df_fit) %in% site_cols]
# change their names
colnames(site_df) <- sub(".*prov:([^]]+)\\].*", "\\1", colnames(site_df))
# empty site df
site_df2 <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per25 = NA,
  fit_a_site_per75 = NA,
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df2$site[i] <- colnames(site_df)[i]         
  site_df2$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df2$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df2$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df2$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df2$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df2


# recover fit parameters for pgsNgrowingdays ####
df_fit_pgsNgrowingday <- as.data.frame(fit_pgsNgrowingdays)

###### Recover treeid ######

# grab treeid 
treeid_cols <- colnames(df_fit_pgsNgrowingday)[grepl("treeid", colnames(df_fit_pgsNgrowingday))]
# remove sigma_asp for now
treeid_cols <- treeid_cols[1:length(treeid_cols)]

treeid_df <- df_fit_pgsNgrowingday[, colnames(df_fit_pgsNgrowingday) %in% treeid_cols]
# change their names
colnames(treeid_df) <- sub(".*treeid:([^]]+)\\].*", "\\1", colnames(treeid_df))

# empty treeid dataframe
treeid_df_pgsNgrowingdays <- data.frame(
  treeid = character(ncol(treeid_df)),
  fit_a_treeid = numeric(ncol(treeid_df)),  
  fit_a_treeid_per5 = NA, 
  fit_a_treeid_per25 = NA,
  fit_a_treeid_per75 = NA,
  fit_a_treeid_per95 = NA
)
for (i in 1:ncol(treeid_df)) { # i = 1
  treeid_df_pgsNgrowingdays$treeid[i] <- colnames(treeid_df)[i]         
  treeid_df_pgsNgrowingdays$fit_a_treeid[i] <- round(mean(treeid_df[[i]]),3)  
  treeid_df_pgsNgrowingdays$fit_a_treeid_per5[i] <- round(quantile(treeid_df[[i]], probs = 0.05), 3)
  treeid_df_pgsNgrowingdays$fit_a_treeid_per25[i] <- round(quantile(treeid_df[[i]], probs = 0.25), 3)
  treeid_df_pgsNgrowingdays$fit_a_treeid_per75[i] <- round(quantile(treeid_df[[i]], probs = 0.75), 3)
  treeid_df_pgsNgrowingdays$fit_a_treeid_per95[i] <- round(quantile(treeid_df[[i]], probs = 0.95), 3)
}
treeid_df_pgsNgrowingdays

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
if (FALSE) {
###### Recover a spp  ######
aspp_cols <- colnames(df_fit_pgsNgrowingday)[grepl("spp", colnames(df_fit_pgsNgrowingday))]
aspp_cols <- aspp_cols[!grepl("Sigma", aspp_cols)]
aspp_cols <- aspp_cols[!grepl("pgsNgrowingdays", aspp_cols)]

aspp_df <- df_fit_pgsNgrowingday[, colnames(df_fit_pgsNgrowingday) %in% aspp_cols]
# change their names
colnames(aspp_df) <- sub(".*spp:([^]]+)\\].*", "\\1", colnames(aspp_df))

#empty aspp df
aspp_df_pgsNgrowingdays <- data.frame(
  spp = character(ncol(aspp_df)),
  fit_a_spp = numeric(ncol(aspp_df)),  
  fit_a_spp_per5 = NA, 
  fit_a_spp_per25 = NA,
  fit_a_spp_per75 = NA,
  fit_a_spp_per95 = NA
)
for (i in 1:ncol(aspp_df)) { # i = 1
  aspp_df_pgsNgrowingdays$spp[i] <- colnames(aspp_df)[i]         
  aspp_df_pgsNgrowingdays$fit_a_spp[i] <- round(mean(aspp_df[[i]]),5)  
  aspp_df_pgsNgrowingdays$fit_a_spp_per5[i] <- round(quantile(aspp_df[[i]], probs = 0.05), 3)
  aspp_df_pgsNgrowingdays$fit_a_spp_per25[i] <- round(quantile(aspp_df[[i]], probs = 0.25), 3)
  aspp_df_pgsNgrowingdays$fit_a_spp_per75[i] <- round(quantile(aspp_df[[i]], probs = 0.75), 3)
  aspp_df_pgsNgrowingdays$fit_a_spp_per95[i] <- round(quantile(aspp_df[[i]], probs = 0.95), 3)
}
aspp_df_pgsNgrowingdays

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
}
###### Recover b spp ###### 
bspp_cols <- colnames(df_fit_pgsNgrowingday)[grepl("pgsNgrowingdays", colnames(df_fit_pgsNgrowingday))]
# remove sigma_aspp for now
bspp_cols <- bspp_cols[!grepl("Sigma", bspp_cols)]
bspp_cols <- bspp_cols[2:5]

bspp_df <- df_fit_pgsNgrowingday[, colnames(df_fit_pgsNgrowingday) %in% bspp_cols]
# change their names
colnames(bspp_df) <- sub(".*spp:([^]]+)\\].*", "\\1", colnames(bspp_df))
#empty spp df
bspp_df_pgsNgrowingdays <- data.frame(
  spp = character(ncol(bspp_df)),
  fit_b_spp = numeric(ncol(bspp_df)),
  fit_b_spp_per5 = NA,
  fit_b_spp_per25 = NA,
  fit_b_spp_per75 = NA,
  fit_b_spp_per95 = NA
)
for (i in 1:ncol(bspp_df)) { # i = 1
  bspp_df_pgsNgrowingdays$spp[i] <- colnames(bspp_df)[i]
  bspp_df_pgsNgrowingdays$fit_b_spp[i] <- round(mean(bspp_df[[i]]),5)
  bspp_df_pgsNgrowingdays$fit_b_spp_per5[i] <- round(quantile(bspp_df[[i]], probs = 0.05), 3)
  bspp_df_pgsNgrowingdays$fit_b_spp_per25[i] <- round(quantile(bspp_df[[i]], probs = 0.25), 3)
  bspp_df_pgsNgrowingdays$fit_b_spp_per75[i] <- round(quantile(bspp_df[[i]], probs = 0.75), 3)
  bspp_df_pgsNgrowingdays$fit_b_spp_per95[i] <- round(quantile(bspp_df[[i]], probs = 0.95), 3)
}
bspp_df_pgsNgrowingdays

###### Recover site ######
site_cols <- colnames(df_fit_pgsNgrowingday)[grepl("prov", colnames(df_fit_pgsNgrowingday))]
# remove for now
site_cols <- site_cols[1:length(site_cols)]
site_cols <- site_cols[!grepl("Sigma", site_cols)]
site_df <- df_fit_pgsNgrowingday[, colnames(df_fit_pgsNgrowingday) %in% site_cols]
# change their names
colnames(site_df) <- sub(".*prov:([^]]+)\\].*", "\\1", colnames(site_df))
# empty site df
site_df_pgsNgrowingdays <- data.frame(
  site = character(ncol(site_df)),
  fit_a_site = numeric(ncol(site_df)),  
  fit_a_site_per5 = NA, 
  fit_a_site_per25 = NA,
  fit_a_site_per75 = NA,
  fit_a_site_per95 = NA
)
for (i in 1:ncol(site_df)) { # i = 1
  site_df_pgsNgrowingdays$site[i] <- colnames(site_df)[i]         
  site_df_pgsNgrowingdays$fit_a_site[i] <- round(mean(site_df[[i]]),3)  
  site_df_pgsNgrowingdays$fit_a_site_per5[i] <- round(quantile(site_df[[i]], probs = 0.05), 3)
  site_df_pgsNgrowingdays$fit_a_site_per25[i] <- round(quantile(site_df[[i]], probs = 0.25), 3)
  site_df_pgsNgrowingdays$fit_a_site_per75[i] <- round(quantile(site_df[[i]], probs = 0.75), 3)
  site_df_pgsNgrowingdays$fit_a_site_per95[i] <- round(quantile(site_df[[i]], probs = 0.95), 3)
}
site_df_pgsNgrowingdays

if (FALSE) {
###### Plot asp ######
aspp_df2$sppfull[which(aspp_df2$spp == "ALNINC")] <- "Alnus incana"
aspp_df2$sppfull[which(aspp_df2$spp == "BETPOP")] <- "Betula populifolia"
aspp_df2$sppfull[which(aspp_df2$spp == "BETPAP")] <- "Betula papyrifera"
aspp_df2$sppfull[which(aspp_df2$spp == "BETALL")] <- "Betula alleghaniensis"

aspp_df_pgsNgrowingdays$sppfull[which(aspp_df_pgsNgrowingdays$spp == "ALNINC")] <- "Gray alder 
(Alnus incana)"
aspp_df_pgsNgrowingdays$sppfull[which(aspp_df_pgsNgrowingdays$spp == "BETPOP")] <- "Gray birch 
(Betula populifolia)"
aspp_df_pgsNgrowingdays$sppfull[which(aspp_df_pgsNgrowingdays$spp == "BETPAP")] <- "Paper birch 
(Betula papyrifera)"
aspp_df_pgsNgrowingdays$sppfull[which(aspp_df_pgsNgrowingdays$spp == "BETALL")] <- "Yellow birch 
(Betula alleghaniensis)"

# order by spp
aspp_df2 <- aspp_df2[order(aspp_df2$spp), ]
aspp_df_pgsNgrowingdays <- aspp_df_pgsNgrowingdays[order(aspp_df_pgsNgrowingdays$spp), ]

# pgsGDD
ggplot(aspp_df2, aes(x = fit_a_spp, y = sppfull, color = sppfull)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_spp_per5, xmax = fit_a_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_spp_per25, xmax = fit_a_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species")+
  theme_minimal() 
ggsave("figures/empiricalData/empiricalData_asp_lmer.jpeg", width = 6, height = 6, units = "in", dpi = 300)

# pgsNgrowingdays
ggplot(aspp_df_pgsNgrowingdays, aes(x = fit_a_spp, y = sppfull, color = sppfull)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_spp_per5, xmax = fit_a_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_spp_per25, xmax = fit_a_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species")+
  theme_minimal() 

}
###### Plot bsp ######
bspp_df2$sppfull[which(bspp_df2$spp == "ALNINC")] <- "Alnus incana"
bspp_df2$sppfull[which(bspp_df2$spp == "BETPOP")] <- "Betula populifolia"
bspp_df2$sppfull[which(bspp_df2$spp == "BETPAP")] <- "Betula papyrifera"
bspp_df2$sppfull[which(bspp_df2$spp == "BETALL")] <- "Betula alleghaniensis"

bspp_df_pgsNgrowingdays$sppfull[which(bspp_df_pgsNgrowingdays$spp == "ALNINC")] <- "Gray alder 
(Alnus incana)"
bspp_df_pgsNgrowingdays$sppfull[which(bspp_df_pgsNgrowingdays$spp == "BETPOP")] <- "Gray birch 
(Betula populifolia)"
bspp_df_pgsNgrowingdays$sppfull[which(bspp_df_pgsNgrowingdays$spp == "BETPAP")] <- "Paper birch 
(Betula papyrifera)"
bspp_df_pgsNgrowingdays$sppfull[which(bspp_df_pgsNgrowingdays$spp == "BETALL")] <- "Yellow birch 
(Betula alleghaniensis)"

# plot bspp for pgsGDD
bspp_df2 <- bspp_df2[order(bspp_df2$spp), ]

ggplot(bspp_df2, aes(x = fit_b_spp, y = sppfull, color = sppfull)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_b_spp_per5, xmax = fit_b_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_b_spp_per25, xmax = fit_b_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "", x = "Species slope values", color = "Tree Species")+
  theme_bw()+  
  scale_y_discrete(limits = rev) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "italic"),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) 
ggsave("figures/empiricalData/empiricalData_bsp_lmer.jpeg", width = 9, height = 6, units = "in", dpi = 300)

# plot bspp for pgsNgrowingdays
bspp_df_pgsNgrowingdays <- bspp_df_pgsNgrowingdays[order(bspp_df_pgsNgrowingdays$spp), ]

ggplot(bspp_df_pgsNgrowingdays, aes(x = fit_b_spp, y = sppfull, color = sppfull)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_b_spp_per5, xmax = fit_b_spp_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_b_spp_per25, xmax = fit_b_spp_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "Species", x = "Species intercept values", color = "Tree Species")+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  ) +
  theme_bw() +
  scale_y_discrete(limits = rev)   

###### Plot asite ######
# add lat to df
site_df_forplot <- site_df2
site_df_forplot$sitefull <- c("Dartmouth College (NH)", "Harvard Forest (MA)", "St-Hyppolyte (Qc)", "White Mountains (NH)")

site_df_forplot$Latitude <- c(44.92,42.55,45.98,44.11)

site_df_forplot$site <- factor(site_df_forplot$site, levels = site_df_forplot$site)
site_df_forplot <- site_df_forplot[order(-site_df_forplot$Latitude), ]

# latitudinal degree difference from arboretum
site_df_forplot$degreedif <- site_df_forplot$Latitude-42.29601035316377

# latitudinal km difference from arboretum
site_df_forplot$kmdif <- site_df_forplot$degreedif*111

site_df_forplot$sitefull <- factor(site_df_forplot$sitefull,
                                   levels = site_df_forplot$sitefull)


ggplot(site_df_forplot, aes(x = fit_a_site, y = kmdif, color = sitefull)) +
  geom_point(size = 6, alpha = 1) + 
  geom_errorbarh(aes(xmin = fit_a_site_per5, xmax = fit_a_site_per95), 
                 width = 0, alpha = 1, linewidth = 0.7) +
  geom_errorbarh(aes(xmin = fit_a_site_per25, xmax = fit_a_site_per75), 
                 width = 0, alpha = 1, linewidth = 2) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Sites intercept values", y = "Latitude distance (km)", color = "Site") +
  scale_y_continuous(breaks = seq(0, 450, by = 50))+
  theme_bw() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10),                  
    legend.key.size = unit(1.5, "lines"),                   
    legend.position = "right"                               
  )
ggsave("figures/empiricalData/empiricalData_asite_lmer.jpeg", width = 8, height = 6, units = "in", dpi = 300)

# Map ####
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# --- Get the map data ---
world <- ne_countries(scale = "medium", returnclass = "sf")

# --- Define bounding box for northeastern North America ---
# Adjust these coordinates as needed
lat_min <- 41
lat_max <- 48
lon_min <- -78
lon_max <- -63

# --- Create example points along a latitudinal gradient ---
# These are arbitrary example locations
locations <- data.frame(
  name = c("Harvard Forest (MA)", "White Mountains (NH)", "Dartmouth College (NH)", "St-Hyppolyte (Qc)"),
  lon = c(-72.2, -71, -70.66, -74.01),
  lat = c(42.55, 44.11, 44.92, 45.98)
)

locations2 <- locations[order(-locations$lat), ]

locations2$col <-  wes_palettes$Darjeeling1[1:4]

special_point <- data.frame(
  name = "Arnold Arboretum of 
  Harvard University (MA)",
  lon = -71.13358611669867,
  lat = 42.29601035316377
)
special_sf <- st_as_sf(special_point, coords = c("lon", "lat"), crs = 4326)

# Convert to sf object
points_sf <- st_as_sf(locations2, coords = c("lon", "lat"), crs = 4326)

# --- Plot the map ---
ggplot(data = world) +
  geom_sf(fill = "white", color = "gray60") +
  geom_sf(data = points_sf, color = locations2$col, size = 4) +
  geom_sf(data = special_sf, color = "#E54E21", shape = 8, size = 6, stroke = 1.2) +
  geom_text(data = locations2, aes(x = lon, y = lat, label = name),
            nudge_y = 0.35, nudge_x = 0, size = 4.5, fontface = "bold") +
  geom_text(data = special_point,
            aes(x = lon, y = lat, label = name),
            nudge_y = 0.2, nudge_x = 2.5, color = "#E54E21", size = 5, fontface = "bold") +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  theme_minimal() +   
  theme(strip.text = element_blank(),                    
        legend.key.height = unit(1.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted")
  )
ggsave("figures/mapSourcePop.jpeg", width = 9, height = 6, units = "in", dpi = 300)


[1] "#FF0000" "#00A08A" "#F2AD00" "#F98400"

 # Fitting empirical data with stan_lmer ####

