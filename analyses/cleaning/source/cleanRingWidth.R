## Started 3 April 2025 ##
# Christophe

# goal of this script is to import the cookie measurements, clean them and change unit from inches to cm

## housekeeping
# rm(list=ls()) 
# options(stringsAsFactors = FALSE)
options(max.print = 150) 

# verification steps
crossverification <- FALSE

## Load Libraries
library(ggplot2)
library(xlsx)

# Set main directory
directory <- "/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses/input"

# List all items and keep only directories
items <- list.files(path = directory, full.names = TRUE)
folders <- items[file.info(items)$isdir]

# Get folder names and exclude undesired ones
folder_names <- basename(folders)
valid_folders <- folder_names[!folder_names %in% c("_notcookies", "_readme.md")]

# Initialize list to store all data
all_data <- list()

# Loop through each folder
for (folder in valid_folders) {
  folder_path <- file.path(directory, folder)
  file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Read all CSVs and tag each row with the folder name
  data_list <- lapply(file_list, function(file) {
    df <- read.csv(file)
    df$sourceFolder <- folder
    return(df)
  })
  
  all_data <- c(all_data, data_list)
}

# Combine all into a single data frame
c <- do.call(rbind, all_data)

### === === === === === ###
##### Read initial table that compiled what tree we had cookies/cores ##### 
# it involves a couple of mistakes when I entered the data that I will fix bellow #
### === === === === === ###
# read table
setwd(directory)
list.files()
d <- read.xlsx("_notcookies/treecookies.xlsx", sheetName = "Sheet1")
# remove _ between provenance and number
d$id <- gsub("(_)([A-Z]+)_([0-9]+[A-Z]?)$", "\\1\\2\\3", d$id)
# paste id and plot
d$idfull <- paste(d$id, paste0("P", d$Plot), sep = "_")
# subset only for species I want
vec <- c("ALNINC", "BETALL", "BETPAP", "BETPOP", "QUERUB")
dsub <- subset(d, species %in% vec )

# Clean cookie and cores when I made mistakes in entering them
### Remove cookie for ALNINC HF9 P6
dsub$cookie.[which(dsub$idfull == "ALNINC_HF_9_6")] <- "0" 
### Remove cookie for ALNINC WM2B P1
dsub$cookie.[which(dsub$idfull == "ALNINC_WM_2B_1")] <- "0"
### Change number of cookies to 2 for BETALL GR13 P15
dsub$cookie.[which(dsub$idfull == "BETALL_GR_13_15")] <- "2"
### Change number of cookies to 2 for BETALL GR13 P1
dsub$cookie.[which(dsub$idfull == "BETALL_GR_13_1")] <- "0"
### Change number of cookie for BETALL_GR9_P3
dsub$cookie.[which(dsub$idfull == "BETALL_GR9_P3")] <- "1"
### Change number of cookie for ALNINC_GR8B_P5
dsub$core[which(dsub$idfull == "ALNINC_GR8B_P5")] <- "1"
### ALNINC_SH5_P2
dsub$cookie.[which(dsub$idfull == "ALNINC_SH5_P2")] <- "1"
### BETPAP_GR5B_P2
dsub$cookie.[which(dsub$idfull == "BETPAP_GR5B_P2")] <- "1"
### BETPOP_GR5_P6
dsub$cookie.[which(dsub$idfull == "BETPOP_GR5_P6")] <- "1"
### ALNINC_HF9_P6
dsub$cookie.[which(dsub$idfull == "ALNINC_HF9_P6")] <- "0"
dsub$core[which(dsub$idfull == "ALNINC_HF9_P6")] <- "1"

cookiesOG <- subset(dsub, cookie. == "1")

coresOG <- subset(dsub, core %in% c("1","2"))

### === === === === === ###
##### Clean labels and years #####
### === === === === === ###
# cleaned name column
c$name <- gsub("[-_]?guides([-_][0-9]+)?\\.tif$", "", sub(":.*", "", c$Label))
# cleaned year column
c$year <- sub(".*:", "", c$Label)
# Create a new table with only year values in the year column i.e. excluding comments
e <- c[grepl("^\\d{2,4}$", c$year), ] 
# Convert 2 digit year into 4 digit years:
e$year <- ifelse(nchar(e$year) == 2, paste0("20", e$year), e$year)
# Convert year column to numeric
e$year <- as.numeric(e$year)
# Create a new column with the length in cm with a conversion factor of 2.54cm because 1 inch is 2.540005 cm
e$lengthCM <- e$Length*2.54
# Create a new df with only necessary columns
d <- e[, c("name", "year", "lengthCM", "sourceFolder")]

#### Extra cleaning steps ####
# clean names that I entered wrong when scanning
# BETALL_WM8B_P95: decision made in CoreXCookies.md
d$name[which(d$name == "BETALL_WM8B_P95")] <- "BETALL_WM8B_P5" 
# BETALL_WM8D_P5: decision made in CoreXCookies.md
d$name[which(d$name == "BETALL_WM8D_P5")] <- "BETALL_WM8B_P5" 
# BETALL_WM_9_P5
d$name[which(d$name == "BETALL_WM_9_P5")] <- "BETALL_WM9_P5" 
# BETALL_WM_8B_P12
d$name[which(d$name == "BETALL_WM_8B_P12")] <- "BETALL_WM8B_P12" 
# BETPOP_WM-7_P2
d$name[which(d$name == "BETPOP_WM-7_P2")] <- "BETPOP_WM7_P2" 
# BETPOP_GRB_P12: decision made in CoreXCookies.md
d$name[which(d$name == "BETPOP_GRB_P12")] <- "BETPOP_GR5B_P12" 
# BETPOP_GR_5B_P12
d$name[which(d$name == "BETPOP_GR_5B_P12")] <- "BETPOP_GR5B_P12" 
# ALNINC_WM_2A_P1
d$name[which(d$name == "ALNINC_WM_2A_P1")] <- "ALNINC_WM2A_P1" 
# ALNINC_WM_5B_P2
d$name[which(d$name == "ALNINC_WM_5B_P2")] <- "ALNINC_WM5B_P2" 
# BETPAP_HF_16A_P6
d$name[which(d$name == "BETPAP_HF_16A_P6")] <- "BETPAP_HF16A_P6" 
# BETPOP_GR_5A_P3
d$name[which(d$name == "BETPOP_GR_5A_P3")] <- "BETPOP_GR5A_P3" 
# BETPOP_GR_5B_P11
d$name[which(d$name == "BETPOP_GR_5B_P11")] <- "BETPOP_GR5B_P11" 
# BETPOP_GR_5C_P12
d$name[which(d$name == "BETPOP_GR_5C_P12")] <- "BETPOP_GR5C_P12" 
# BETPOP_HF__P6
d$name[which(d$name == "BETPOP_HF__P6")] <- "BETPOP_HFNA_P6"
# BETAL_GR9_P3
d$name[which(d$name == "BETAL_GR9_P3")] <- "BETALL_GR9_P3"
# BETALL_SH9_PNA: decision made in CoreXCookies.md
d$name[which(d$name == "BETALL_SH9_PNA")] <- "BETALL_SH9_P9"
# BETALL_WM8_P12: 
d$name[which(d$name == "BETALL_WM8_P12")] <- "BETALL_WM8B_P12" 
# ALNINC_WM8_P1: decision in cookie scraping. Not WM8, but WM5.
d$name[which(d$name == "ALNINC_WM8_P1" & d$sourceFolder == "cookies")] <- "ALNINC_WM5_P1" 
# BETALL_SH9_P6
d$name[which(d$name == "BETALL_SH5_P6" & d$sourceFolder == "coresUnconfident")] <- "BETALL_SH9_P6"
# ALNINC_SH5_P3
d$name[which(d$name == "ALNINC_SH3_P3" & d$sourceFolder == "cookies")] <- "ALNINC_SH5_P3" 
# BETPAP_HF16_P2
d$name[which(d$name == "BETPAP_HF16_B2" & d$sourceFolder == "coresWithoutCookies")] <- "BETPAP_HF16_P2" 
# ALNINC_HF1_P16: ALNINC_HF3_P16 on cookie, but likely a transcription mistake as we have HF1 in wildhell, but not HF3
d$name[which(d$name == "ALNINC_HF3_P16" & d$sourceFolder == "cookies")] <- "ALNINC_HF1_P16" 


### === === === === === ###
##### Verification steps #####
### === === === === === ###
# compare if I have the data from all cookies in the og dataset
if(crossverification) {
  vcook <- c("cookies", "cookiesUnconfident", "coresUnconfident", "coresWithoutCookies")
  coresandcookies <- subset(d, sourceFolder %in% vcook)
  listCookieNames <- unique(coresandcookies$name)
  cookieInOG <- listCookieNames[which(!listCookieNames%in%obsdata$name)]
  # compare if I have scanned the cookies from the data from og dataset 
  OGtoCookies <- unique(obsdata$name[which(!obsdata$name%in%listCookieNames)])
}



### === === === === === ###
# Start playing with the data #
### === === === === === ###
# Take the mean of each year and each tree replicate
agg <- aggregate(d$lengthCM, by = list(d$name, d$year, d$sourceFolder), FUN = mean)

# rename the columns
colnames(agg) <- c("name", "year", "sourceFolder", "lengthCM")
agg$Yearcor <- agg$year

### === === === === === ###
##### Cross date cookie/core #####
### === === === === === ###
# vec of all cores with cookies
corewcookie <- subset(agg, sourceFolder == "coresWithCookies")
corevec <- unique(corewcookie$name)
# Start by selecting cookies for which we have cores
# agg$Yearcor <- agg$year
cookiewcore <- subset(agg, name %in% corevec)

# histograms of every cookie X core 
# start with first 5 names
# onetoten <- subset(agg, name %in% unique(cookiewcore$name)[1:10])
alninc <- subset(cookiewcore, grepl("ALNINC", name))
betall <- subset(cookiewcore, grepl("BETALL", name))
betpap <- subset(cookiewcore, grepl("BETPAP", name))
betpop <- subset(cookiewcore, grepl("BETPOP", name))

# quartz()
ggplot(betpap, aes(x = factor(Yearcor), y = lengthCM, fill = sourceFolder)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ name, scales = "fixed") +
  theme_minimal() +
  labs(
    title = "Cookie Length by year and Core name",
    x = "year",
    y = "Length (cm)"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# dev.off()

### === === === === === === === === === === ###
###### Change 2022 to 2023 when necessary for cookies with cores ######
### === === === === === === === === === === ###
# add yearCor column to fix year patterns after verification 
d$Yearcor <- d$year
# BETALL_GR12_P1
d$Yearcor[d$name == "BETALL_GR12_P1" & d$sourceFolder == "cookies"] <- 
  agg$year[agg$name == "BETALL_GR12_P1" & agg$sourceFolder == "cookies"] + 1
# BETALL_GR12_P1
d$Yearcor[d$name == "BETALL_GR12_P1" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETALL_GR12_P1" & d$sourceFolder == "cookies"] + 1
# BETALL_WM8A_P5
d$Yearcor[d$name == "BETALL_WM8A_P5" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETALL_WM8A_P5" & d$sourceFolder == "cookies"] + 1
# BETPAP_HF16_P12
d$Yearcor[d$name == "BETPAP_HF16_P12" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPAP_HF16_P12" & d$sourceFolder == "cookies"] + 1
# BETPAP_HF16A_P6
d$Yearcor[d$name == "BETPAP_HF16A_P6" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPAP_HF16A_P6" & d$sourceFolder == "cookies"] + 1
# BETPAP_SH1A_P1
d$Yearcor[d$name == "BETPAP_SH1A_P1" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPAP_SH1A_P1" & d$sourceFolder == "cookies"] + 1
# BETPOP_GR5A_P12
d$Yearcor[d$name == "BETPOP_GR5A_P12" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_GR5A_P12" & d$sourceFolder == "cookies"] + 1
# BETPOP_GR5B_P12
d$Yearcor[d$name == "BETPOP_GR5B_P12" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_GR5B_P12" & d$sourceFolder == "cookies"] + 1
# BETPOP_HF1_P3
d$Yearcor[d$name == "BETPOP_HF1_P3" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_HF1_P3" & d$sourceFolder == "cookies"] + 1
# BETPOP_WM7_P6
d$Yearcor[d$name == "BETPOP_WM7_P6" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_WM7_P6" & d$sourceFolder == "cookies"] + 1
# BETPOP_XX_P3
d$Yearcor[d$name == "BETPOP_XX_P3" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_XX_P3" & d$sourceFolder == "cookies"] + 1

### === === === === === === === === === === ###
######Change 2022 to 2023 when necessary for cookies without cores ######
### === === === === === === === === === === ###
# BETALL_GR12_P1
d$Yearcor[d$name == "BETALL_WM8_P1" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETALL_WM8_P1" & d$sourceFolder == "cookies"] + 1
# BETPOP_HF16A_P2
d$Yearcor[d$name == "BETPOP_HF16A_P2" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_HF16A_P2" & d$sourceFolder == "cookies"] + 1
# BETPOP_HF3A_P6
d$Yearcor[d$name == "BETPOP_HF3A_P6" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_HF3A_P6" & d$sourceFolder == "cookies"] + 1
# BETPOP_WM8_PNA
d$Yearcor[d$name == "BETPOP_WM8_PNA" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_WM8_PNA" & d$sourceFolder == "cookies"] + 1
# BETPOP_GR5C_P12
d$Yearcor[d$name == "BETPOP_GR5C_P12" & d$sourceFolder == "cookies"] <- 
  d$year[d$name == "BETPOP_GR5C_P12" & d$sourceFolder == "cookies"] + 1
    
# select cookies for which we have cores
cookienocore <- subset(agg, !(name %in% unique(corewcookie$name)))
alninc <- subset(cookienocore, grepl("ALNINC", name))
betall <- subset(cookienocore, grepl("BETALL", name))
betpap <- subset(cookienocore, grepl("BETPAP", name))
betpop <- subset(cookienocore, grepl("BETPOP", name))

### === === === === === === === === === === ###
###### Change 2022 to 2023 when necessary for cores without cookies ######
### === === === === === === === === === === ###
# Start by selecting cookies for which we have cores
corenocookie <- subset(agg, sourceFolder == "coresWithoutCookies")
alninc <- subset(corenocookie, grepl("ALNINC", name))
betall <- subset(corenocookie, grepl("BETALL", name))
betpap <- subset(corenocookie, grepl("BETPAP", name))
betpop <- subset(corenocookie, grepl("BETPOP", name))



# rename df and remove old year col
rw <- d[, !(names(d) %in% "year")]
# rename year column using base r code
# add species col
names(rw)[names(rw) == "Yearcor"] <- "year"
names(rw)[names(rw) == "sourceFolder"] <- "sampleType"
rw$spp <- sub("_.*", "", rw$name)

library(dplyr)
library(stringr)
# remove duplicates 
rwnodup <- rw[!duplicated(rw$name),]
rwnodup %>%
  count(spp)

rwnodup2 <- rwnodup %>%
  mutate(
    provenance = str_extract(name, "_([A-Z]{2})") %>% str_remove("_")
  )

rwnodup2 %>%
  count(spp, provenance)
