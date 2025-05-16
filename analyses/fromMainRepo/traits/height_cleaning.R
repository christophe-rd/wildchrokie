# code written by Tolu A., summer of 2022
# edited by Deirdre on Sept 1, 2022
rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyr)
library(ggplot2)
library(readr)
library(dplyr)

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/wildhellgarden")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/projects/others/deirdre/synchrony")
} else{
  setwd("/home/deirdre/wildhellgarden") # for midge
}
#Could you clean and combine all the different height data files into one?
# 1. There is height data from 2019 which is clean and in m: analyses/2019_data/2019_heightanddbh.csv
# Could you convert this file to just be ind, plot, and height columns without the blank rows?
# 2. There is also height data in the 2020_data that has not been cleaned and the values will need to be 
# corrected for the measuring rod issues and converted to m
# 3. And of course the 2022 data needs to be cleaned and corrected as well
# 4-7. Same instructions to look for outliers and make plots as Grace above.

#Notes on the measure stick can be found here: misc/Jun2021measuringstick/measuringstick.txt
# first compare the actual conversions with our rod conversions in m
# 1ft = 3.2787m

# 2019 Cleaning
wildhell_2019 <- read_csv("Documents/WL2022/WildHell/wildhellgarden/analyses/2019_data/2019_heightanddbh.csv")

wildhell_2019 <- wildhell_2019 %>%
  select(IND, Plot, Height)

wildhell_2019 <- wildhell_2019[complete.cases(wildhell_2019),]
colnames(wildhell_2019)[which(names(wildhell_2019) == "IND")] <- "ID"
wildhell_2019$Year <- "2019"

wildhell_2019$Height <- as.numeric(as.character(wildhell_2019$Height)) / 100
colnames(wildhell_2019)[which(names(wildhell_2019) == "Height")] <- "Height_m"

# 2021 Cleaning
wildhell_2021 <- read_csv("Documents/WL2022/WildHell/wildhellgarden/analyses/2021_data/heights2021.csv")

wildhell_2021$Height <- paste0(wildhell_2021$Feet, ".", wildhell_2021$Inches2)
wildhell_2021$Height[which(wildhell_2021$Height == "NA.NA")] <- NA

# cleaning the 0.5 inches
temp2021 <- wildhell_2021 %>% dplyr::filter(nchar(Height) >= 5)

wildhell_2021$Height[which(wildhell_2021$Height == "3.0.5")] <- 3.05
wildhell_2021$Height[which(wildhell_2021$Height == "12.0.5")] <- 12.05
wildhell_2021$Height[which(wildhell_2021$Height == "6.0.5")] <- 6.05
wildhell_2021$Height[which(wildhell_2021$Height == "13.0.5")] <- 13.05
wildhell_2021$Height[which(wildhell_2021$Height == "4.0.5")] <- 4.05
wildhell_2021$Height[which(wildhell_2021$Height == "3.0.5")] <- 3.05
wildhell_2021$Height[which(wildhell_2021$Height == "3.NA")] <- 3.0
wildhell_2021$Height[which(wildhell_2021$Height == "4.NA")] <- 4.0

# converting heights

wildhell_2021$Height_m <- as.numeric(as.character(wildhell_2021$Height)) / 3.2787
cleaned_2021 <- wildhell_2021 %>%
  select(ID, Plot, Height_m)
cleaned_2021[,'Height_m']=round(cleaned_2021[,'Height_m'],2)
cleaned_2021 <- cleaned_2021[!is.na(cleaned_2021$ID),]
cleaned_2021$Year <- "2021"

# 2022 Cleaning
wildhell_2022 <- read_csv("Documents/WL2022/WildHell/wildhellgarden/analyses/2022_data/height2022.csv")

wildhell_2022$Height <- paste0(wildhell_2022$Feet, ".", wildhell_2022$Inches2)
wildhell_2022$Height[which(wildhell_2022$Height == "NA.NA")] <- NA

temp2022 <- wildhell_2022 %>% dplyr::filter(nchar(Height) >= 5)
wildhell_2022$Height[which(wildhell_2022$Height == "3.0.5")] <- 3.05
wildhell_2022$Height[which(wildhell_2022$Height == "4.0.5")] <- 4.05
wildhell_2022$Height[which(wildhell_2022$Height == "5.0.5")] <- 5.05

# converting heights

wildhell_2022$Height_m <- as.numeric(as.character(wildhell_2022$Height)) / 3.2787
cleaned_2022 <- wildhell_2022 %>%
  select(id, Plot, Height_m)
cleaned_2022[,'Height_m']=round(cleaned_2022[,'Height_m'],2)
colnames(cleaned_2022)[which(names(cleaned_2022) == "id")] <- "ID"
cleaned_2022$Year <- "2022"

# combining all three data frames
wildhell_height <- rbind(wildhell_2019, cleaned_2021, cleaned_2022)

plot <- wildhell_height[which(wildhell_height$Plot == "3"),] 
ggplot(wildhell_height, aes(x = Plot, y = Height_m)) + geom_point()



