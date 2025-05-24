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

# need to run getAll.R first
# read phenology for wildchrokie
wildphen <- obsdata 
wildphen$species<-substr(wildphen$Name, 0,6)

# read phenology for treespotters
spotphen <- read.csv2("cleanTS.csv", header= TRUE, sep =",")

# select main phenological observations for both
str(wildphen)
col <- c("budburst", "leafout", "flowers", "fruits", "coloredLeaves")

# select only species within betulacea because its the only thing I have in the common garden df
spotphen_bet <- subset(spotphen, genus == "Betula")

coltoselect <- c( phentomatch)
# select only columns in treespotters for the observations that match the common garden
spotphen_bet <- subset(spotphen_bet, select = c("plantNickname", "Common_Name", "year", col))

# same but for the common garden
col <- c("budburst", "leafout", "flowers", "fruit", "leafcolor")
wildphen <- subset(wildphen, select = c("Name", "species", "Year", col))



# change colnames for them to match
spotphen_bet_matched <- spotphen_bet %>%
  rename(
    Name = plantNickname,
    Year = year,
    leafcolor = coloredLeaves,
    fruit = fruits
  ) %>%
  select(Name, species, Year, budburst, leafout, flowers, fruit, leafcolor)

wildphen$Source <- "WildHell"
spotphen_bet_matched$Source <- "Treespotters"

combined <- rbind(wildphen, spotphen_bet_matched)

combined_long <- pivot_longer(
  combined,
  cols = c("budburst", "leafout", "flowers", "fruit", "leafcolor"),
  names_to = "Phenophase",
  values_to = "DOY"
)
combined_long$Year <- as.character(combined_long$Year)
# open device
quartz()

# all species within the betulacea family
betulacea <- ggplot(combined_long, aes(x = Year, y = DOY, color = Source)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.2)) +
  facet_wrap(~ Phenophase, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("WildHell" = "#1b9e77", "Treespotters" = "#d95f02")) +
  labs(
    title = "",
    x = "Year",
    y = "doy"
  )
ggsave("figures/betulacea.png", plot = betulacea, width = 10, height = 8, dpi = 300)

# all species within the betula genus
betula <- subset(combined_long, species != "ALNINC")
betulafig <- ggplot(betula, aes(x = Year, y = DOY, color = Source)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.2)) +
  facet_wrap(~ Phenophase, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("WildHell" = "#1b9e77", "Treespotters" = "#d95f02")) +
  labs(
    title = "",
    x = "Year",
    y = "doy"
  )
ggsave("figures/betula.png", plot = betulafig, width = 10, height = 8, dpi = 300)

# the only overlapping species
all <- c("alleghaniensis","BETALL")
alleng <- subset(combined_long, species %in% all)

allengfig <- ggplot(alleng, aes(x = Year, y = DOY, color = Source)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.2)) +
  facet_wrap(~ Phenophase, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("WildHell" = "#1b9e77", "Treespotters" = "#d95f02")) +
  labs(
    title = "",
    x = "Year",
    y = "DOY"
  )
ggsave("figures/alleng.png", plot = allengfig, width = 10, height = 8, dpi = 300)

dev.off()
