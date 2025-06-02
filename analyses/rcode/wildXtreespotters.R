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


summary_stats <- betall %>%
  group_by(species, Source, Year) %>%
  summarise(
    across(
      c(budburst, leafout, leafcolor),
      list(mean = ~mean(., na.rm = TRUE),
           sd = ~sd(., na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# go with only budburst for now
buburst <- summary_stats[,1:4]
# convert to wide
buburst$species <- rep("betall", 6)
buburst$Year <- as.character(buburst$Year)
buburst_wide <- buburst %>%
  pivot_wider(names_from = Source, values_from = budburst_mean)


buburst <- ggplot(buburst_wide, aes(x = Treespotters, y = WildHell, color = Year)) +
  geom_point(size = 3) +
  labs(x = "Treespotters Budburst Mean", 
       y = "WildHell Budburst Mean",
       title = "budburst dates for yellow birch") +
  theme_minimal()
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/")
ggsave("figures/buburst_yellow_birch.jpeg", plot = buburst, width = 8, height = 6)

#try with uncertainty
buburst_wide <- buburst %>%
  pivot_wider(names_from = Source, 
              values_from = c(budburst_mean, budburst_sd))

ggplot(buburst_wide, aes(x = budburst_mean_Treespotters, y = budburst_mean_WildHell)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = budburst_mean_WildHell - budburst_sd_WildHell, 
                    ymax = budburst_mean_WildHell + budburst_sd_WildHell), width = 0.2) +
  geom_errorbarh(aes(xmin = budburst_mean_Treespotters - budburst_sd_Treespotters,
                     xmax = budburst_mean_Treespotters + budburst_sd_Treespotters), height = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Budburst DOY (Treespotters)",
       y = "Budburst DOY (WildHell)",
       title = "Budburst DOY Comparison with Uncertainty (Yellow birch)") +
  theme_minimal()
