library(tidyr)
library(ggplot2)
library(readr)
library(dplyr)


#Could you clean and combine all the different height data files into one?
# 1. There is height data from 2019 which is clean and in m: analyses/2019_data/2019_heightanddbh.csv
# Could you convert this file to just be ind, plot, and height columns without the blank rows?
# 2. There is also height data in the 2020_data that has not been cleaned and the values will need to be 
# corrected for the measuring rod issues and converted to m
# 3. And of course the 2022 data needs to be cleaned and corrected as well
# 4-7. Same instructions to look for outliers and make plots as Grace above.

#Notes on the measure stick can be found here: misc/Jun2021measuringstick/measuringstick.txt

wildhell_2019 <- read_csv("Documents/WL2022/WildHell/wildhellgarden/analyses/2019_data/2019_heightanddbh.csv")
View(wildhell_2019)

wildhell_2019 <- wildhell_2019 %>%
  select(IND, Plot, Height)

wildhell_2019 <- wildhell_2019[complete.cases(wildhell_2019),]



#converting 2019 height data to just be ind, plot, and height columns without the blank rows?

