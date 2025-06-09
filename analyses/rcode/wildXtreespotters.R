# Wildchrokie model
# CRD 23 April 2025
# Goal is to figure out whether treespotters phenology data matches the one for the common garden in the objective of building a predictive model that woul estimate phenological dates for the missing years in the commmon garden.


# Goal: build a model to understand the relationship between growth and growing degree days using the tree cookies collected in the wildhell common garden in spring of 2023

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 200) 

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)
library(dplyr)

source("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/cleaning/getAll.R")

setwd("/Users/christophe_rouleau-desrochers/github/coringtreespotters/analyses/output")
# read phenology for treespotters
spotphen <- read.csv2("cleanTS.csv", header= TRUE, sep =",")

# need to run getAll.R first
# read phenology for wildchrokie
wildphen <- cgclean 
# wildphen$species<-substr(wildphen$Name, 0,6)

# select main phenological observations for both
str(wildphen)
col <- c("budburst", "leafout", "flowers", "fruits", "coloredLeaves")

# select only species within betulacea because its the only thing I have in the common garden df
# select only columns in treespotters for the observations that match the common garden
spotphen_sub <- subset(spotphen, select = c("plantNickname", "Common_Name", "year", col))

# same but for the common garden
col <- c("budburst", "leafout", "flowers", "fruit", "leafcolor")
wildphen <- subset(wildphen, select = c("Name", "spp", "year", col))


# change colnames for them to match
wildphen <- wildphen %>%
  rename(
    species = spp,
    Year = year,
  ) 
# change colnames for them to match
spotphen_sub_matched <- spotphen_sub %>%
  rename(
    Name = plantNickname,
    species = Common_Name,
    Year = year,
    leafcolor = coloredLeaves,
    fruit = fruits
  ) 

wildphen$Source <- "WildHell"
spotphen_sub_matched$Source <- "Treespotters"

combined <- rbind(wildphen, spotphen_sub_matched)

# start with species that I have in both datasets
### betall
betall <- subset(combined, species %in% c("BETALL", "Yellow birch"))

# keep only years that we have data for wildhell
betall <- subset(betall, Year %in% c(2018:2020))

summary_df <- betall %>%
  group_by(species, Year, Source) %>%
  summarise(
    mean_budburst = mean(budburst, na.rm = TRUE),
    sd_budburst = sd(budburst, na.rm = TRUE),
    n = n(),
    se_budburst = sd_budburst / sqrt(n)
  )

# convert to wide
summary_df$species <- rep("betall", 6)
summary_df$Year <- as.character(summary_df$Year)

ts <- subset(summary_df, Source == "Treespotters")
colnames(ts) <- c("species", "Year", "Source", "mean_budburst_ts", "sd_budburst_ts","n","se_budburst_ts")
# do the same for wh
wh <- subset(summary_df, Source == "WildHell")
colnames(wh) <- c("species", "Year", "Source", "mean_budburst_wh", "sd_budburst_wh","n","se_budburst_wh")

merged <- merge(ts, wh, by =c("species", "Year"))

betall <- ggplot(merged, aes(x = mean_budburst_ts, y = mean_budburst_wh, color = Year)) +
  geom_point(size = 3) +
  
  # Vertical SE bars (WildHell)
  geom_segment(aes(x = mean_budburst_ts,
                   xend = mean_budburst_ts,
                   y = mean_budburst_wh - se_budburst_wh,
                   yend = mean_budburst_wh + se_budburst_wh)) +
  
  # Horizontal SE bars (Treespotters)
  geom_segment(aes(x = mean_budburst_ts - se_budburst_ts,
                   xend = mean_budburst_ts + se_budburst_ts,
                   y = mean_budburst_wh,
                   yend = mean_budburst_wh)) +
  
  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  
  labs(x = "Budburst DOY (Treespotters)",
       y = "Budburst DOY (WildHell)",
       title = "yellow birch",
       color = "Year") +
  theme_minimal()

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")
ggsave("figures/buburst_yellow_birch.jpeg", plot = betall, width = 8, height = 6)


### === === === === === ===
# now with betula
### === === === === === ===
betula<- subset(combined, species %in% c("BETPAP", "BETPOP", "BETALL", "River birch", "Yellow birch"))
betula$genus <- "betula"

betula <- subset(betula, Year %in% c(2018:2020))

sum_betula <- betula %>%
  group_by(genus, Year, Source) %>%
  summarise(
    mean_budburst = mean(budburst, na.rm = TRUE),
    sd_budburst = sd(budburst, na.rm = TRUE),
    n = n(),
    se_budburst = sd_budburst / sqrt(n)
  )
sum_betula

sum_betula$Year <- as.character(sum_betula$Year)

ts_betula <- subset(sum_betula, Source == "Treespotters")
colnames(ts_betula) <- c("species", "Year", "Source", "mean_budburst_ts", "sd_budburst_ts","n","se_budburst_ts")
# do the same for wh
wh_betula <- subset(sum_betula, Source == "WildHell")
colnames(wh_betula) <- c("species", "Year", "Source", "mean_budburst_wh", "sd_budburst_wh","n","se_budburst_wh")

merged <- merge(ts_betula, wh_betula, by =c("species", "Year"))

betula <- ggplot(merged, aes(x = mean_budburst_ts, y = mean_budburst_wh, color = Year)) +
  geom_point(size = 3) +
  
  # Vertical SE bars (WildHell)
  geom_segment(aes(x = mean_budburst_ts,
                   xend = mean_budburst_ts,
                   y = mean_budburst_wh - se_budburst_wh,
                   yend = mean_budburst_wh + se_budburst_wh)) +
  
  # Horizontal SE bars (Treespotters)
  geom_segment(aes(x = mean_budburst_ts - se_budburst_ts,
                   xend = mean_budburst_ts + se_budburst_ts,
                   y = mean_budburst_wh,
                   yend = mean_budburst_wh)) +
  
  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  
  labs(x = "Budburst DOY (Treespotters)",
       y = "Budburst DOY (WildHell)",
       title = "betula",
       color = "Year") +
  theme_minimal()
ggsave("figures/budburst_betula.jpeg", plot = betula, width = 8, height = 6)
