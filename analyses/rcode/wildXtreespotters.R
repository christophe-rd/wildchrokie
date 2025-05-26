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
wildphen <- cgclean 
wildphen$species<-substr(wildphen$Name, 0,6)

# read phenology for treespotters
spotphen <- read.csv2("cleanTS.csv", header= TRUE, sep =",")

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
  ) %>%
  select(Name, species, Year, budburst, leafout, flowers, fruit, leafcolor)

# change colnames for them to match
spotphen_sub_matched <- spotphen_sub %>%
  rename(
    Name = plantNickname,
    species = Common_Name,
    Year = year,
    leafcolor = coloredLeaves,
    fruit = fruits
  ) %>%
  select(Name, species, Year, budburst, leafout, flowers, fruit, leafcolor)

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

wildhellsub <- subset(summary_stats, Source == "WildHell")[, c("Year", "budburst_mean","budburst_sd")]
colnames(wildhellsub) <- c("Year", "wild_budburst_mean", "wild_budburst_sd")

treespotsub <- subset(summary_stats, Source == "Treespotters")[, c("Year", "budburst_mean", "budburst_sd")]
colnames(treespotsub) <- c("Year", "treespot_budburst_mean", "treespot_budburst_sd")

binded <- cbind(wildhellsub, treespotsub)[, c(1:3,5:6)]

ggplot(binded)+
  geom_point(aes(x=wild_budburst_mean, y=treespot_budburst_mean))
