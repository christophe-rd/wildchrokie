# Data Cleaning by Grace
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(ggplot2)
library(gridExtra)

wildhell_2019 <- read.csv("analyses/input/WildHell_2019_traits.csv")
wildhell_2022 <- read.csv("analyses/input/2022_WildHellTraits_Original.csv")

# Checking for typos
wildhell_2022$species <- toupper(wildhell_2022$species)
unique(wildhell_2022$species)

unique(wildhell_2022$site)
unique(wildhell_2022$leaf)
unique(wildhell_2022$Plot)

# Removing ID data with multiple entries (indicated by a "/")
wildhell_2022 <- wildhell_2022[!grepl("/", wildhell_2022$wet.wt),]

wildhell_2022$wet.wt <- as.numeric(wildhell_2022$wet.wt)

# Conversions cm2 to mm2
wildhell_2022$area <- wildhell_2022$area * 100
wildhell_2019$area <- wildhell_2019$area * 100

# dry mass g to mg
wildhell_2019$dr.wt <- wildhell_2019$dr.wt * 1000
wildhell_2022$dr.wt <- wildhell_2022$dr.wt * 1000


# SLA Calculations (area / dry weight)
wildhell_2019$SLA <- (wildhell_2019$area/wildhell_2019$dr.wt)
wildhell_2022$SLA <- (wildhell_2022$area/wildhell_2022$dr.wt)

# LDMC Calculations (dry / fresh)
wildhell_2019$LDMC <- (wildhell_2019$dr.wt / wildhell_2019$wet.wt)
wildhell_2022$LDMC <- (wildhell_2022$dr.wt / wildhell_2022$wet.wt)


# Looking at the data
plot_2022 <- ggplot(data = wildhell_2022, aes(x= LDMC, y= SLA, color = species))+
    geom_point()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SLA vs. LDMC (2022 WildHell Dataset)")
plot_2022

# An SLA of >400 doesn't seem right - it looks like the surface area was incorrectly measured.
# I went to measure this leaf in ImageJ and it looks like the image might have a different
# scale then the other PDFs which might have caused the weird measurement (and it's also the only jpeg file...)
# I will be omitting the individual for simplicity

# Regarding the outliers of LDMC > 1000 ... these all appear to be due to
# missing a decimal place reported in the dry weight column
# So, the following individuals were shifted by one decimal place in a new csv file:
# 3_AMECAN_GR_10B_1, 8_SPITOM_HF_13_2, 11_MYRGAL_WM_6B_1, 5_AMECAN_GR_8_1,
# 4_SPITOM_WM_5D_1, 4_SPIALB_SH_7_1, 1_BETPOP_WM_7_1, 11_MYRGAL_WM_6A_1, 13_DIELON_GR_8_2
# 7_DIELON_HF_1A_1, and 7_DIELON_HF_1A_2

wildhell_2022_updated <- read.csv("analyses/input/wildHellTraits2022.csv")
wildhell_2022_updated$species <- toupper(wildhell_2022_updated$species)
wildhell_2022_updated <- wildhell_2022_updated[!grepl("/", wildhell_2022_updated$wet.wt),]
wildhell_2022_updated$area <- wildhell_2022_updated$area * 100
wildhell_2022_updated$dr.wt <- wildhell_2022_updated$dr.wt * 1000

wildhell_2022_updated$wet.wt <- as.numeric(wildhell_2022_updated$wet.wt)
wildhell_2022_updated$SLA <- (wildhell_2022_updated$area/wildhell_2022_updated$dr.wt)
wildhell_2022_updated$LDMC <- (wildhell_2022_updated$dr.wt / wildhell_2022_updated$wet.wt)
wildhell_2022_updated <- subset(wildhell_2022_updated, SLA < 400)

wildhell_2022_updated <- subset(wildhell_2022_updated, LDMC != "NA")

# Plotting the corrected dataset
updated_plot_2022 <- ggplot(data = wildhell_2022_updated, aes(x= LDMC, y= SLA, color = species))+
  geom_point()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SLA vs. LDMC (2022 WildHell Dataset)")
updated_plot_2022

# Additional outliers:
# 8_MYRGAL_SH_2_2 --> I believe the wet weight was incorrectly entered in the field (it should be 0.0988 instead of 0.988),
# Based on the measurements of other members of the same species
# 4_SPIALB_GR_12_1 --> I think dr.wt was incorrectly entered.

# Plotting the 2019 dataset
wildhell_2019 <- subset(wildhell_2019, LDMC != "NA")
plot_2019 <- ggplot(data = wildhell_2019, aes(x= LDMC, y= SLA, color = species))+
  geom_point()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SLA vs. LDMC (2019 WildHell Dataset)")
plot_2019

# Comparing the two datasets:

p1 <- updated_plot_2022+xlim(0, 600)+ylim(0, 175)
p2 <- plot_2019 +xlim(0, 600)+ylim(0, 175)
comparison_plot <- grid.arrange(p1, p2, nrow = 1)


# 
# pdf(file = "SLA_LDMC_Scatterplot_2019vs2022.pdf", width = 8, height = 6)
# plot(comparison_plot)
# dev.off()

# SLA Distributions

SLA_dist_2022 <- ggplot( data = wildhell_2022_updated, aes(x = SLA, fill = species))+
  geom_histogram()+
  facet_wrap(~species)+theme(axis.line = element_line(colour = "black"),
                                            panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank(),
                                            panel.border = element_blank(),
                                            panel.background = element_blank())+
  ggtitle("SLA Distribution (2022 WildHell Dataset)")


# pdf(file = "WildHell_2022_SLA_Distribution.pdf", width = 8, height = 6)
# plot(SLA_dist_2022)
# dev.off()

SLA_dist_2019 <- ggplot( data = wildhell_2019, aes(x = SLA, fill = species))+
  geom_histogram()+
  facet_wrap(~species)+theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank())+
  ggtitle("SLA Distribution (2019 WildHell Dataset)")
SLA_dist_2019

# pdf(file = "WildHell_2019_SLA_Distribution.pdf", width = 8, height = 6)
# plot(SLA_dist_2019)
# dev.off()

SLA_dist_comparison <- grid.arrange(SLA_dist_2019, SLA_dist_2022, nrow = 1)
# 
# pdf(file = "SLA_Distribution_2019vs2022.pdf", width = 12, height = 6)
# plot(SLA_dist_comparison)
# dev.off()

# Histogram of ALL species:
# I had this colour palette and I never get the chance to use it
colour_palette <-c("#510142", "#9E0142", "#D53E4F",
             "#F46D43",
             "#FDAE61",
             "#FEE08B",
             "#FFFFBF", "#E6F598",
             "#ABDDA4",
             "#66C2A5", "#3288BD",
             "#5E4FA2",
             "#372F5B",
             "#225677")


# LDMC Distribution
LDMC_dist_2022 <- ggplot( data = wildhell_2022_updated, aes(x = LDMC, fill = species))+
  geom_histogram()+
  facet_wrap(~species)+theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank())+
  ggtitle("LDMC Distribution (2022 WildHell Dataset)")

LDMC_dist_2019 <- ggplot( data = wildhell_2019, aes(x = LDMC, fill = species))+
  geom_histogram()+
  facet_wrap(~species)+theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank())+
  ggtitle("LDMC Distribution (2019 WildHell Dataset)")

# pdf(file = "WildHell_2022_LDMC_Distribution.pdf", width = 8, height = 6)
# plot(LDMC_dist_2022)
# dev.off()
# 
# pdf(file = "WildHell_2019_LDMC_Distribution.pdf", width = 8, height = 6)
# plot(LDMC_dist_2019)
# dev.off()

# Stem Density Analysis
stem_density <- read.csv("analyses/2022_data/height_stem_2022.csv")

sort(unique(stem_density$species)) # 
stem_density$species[stem_density$species == "BETPOPX"] <- "BETPOP"
stem_density$species[stem_density$species == "VACMR"] <- "VACMYR"

unique(stem_density$site) # Possible error: there are sites "X" and "XX" - should they both exist?
stem_density$site[stem_density$site == "GB"] <- "GR"

unique(stem_density$Plot) # No errors

# Calculating SSD (mass/vol); converting g to mg
# Removing duplicate weight entries:
stem_density <- stem_density[!grepl("/", stem_density$Weight),]
stem_density <- stem_density[!grepl("/", stem_density$StemVol),]

stem_density$Weight <- as.numeric(stem_density$Weight)
stem_density$StemVol <- as.numeric(stem_density$StemVol)
stem_density$Weight <- stem_density$Weight * 1000
stem_density$ssd <- (stem_density$Weight / stem_density$StemVol)

# Checking for typos
unique(stem_density$species)

# Checking the plot:
ssd_2022_test <- ggplot(data = stem_density, aes(x= Weight, y= ssd, color = species))+
  geom_point()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SSD Scatterplot (2022 WildHell Dataset)")
ssd_2022_test

# Issues with data: there are two entries (in the CSV template) for 11_SPITOM_GR_1A
# i.e. the ID "11_SPITOM_GR_1A" shows up twice
# I will be removing these for graphing purposes
stem_density2 <- subset(stem_density, id != "SPITOM_GR_1A")

# Outlier 6_BETPOP_HF_3A (The highest point) -> volume data entry was correct, weight entered is comparable to other BETPOP entries
betpop_subset <- subset(stem_density2, species == "BETPOP")

# Outlier 8_SPITOM_GR_1 (Second highest point) -> volume data entry was correct, weight measurement seems a little big but still possible
spitom_subset <- subset(stem_density2, species == "SPITOM")

# Outlier 8_SPIALB_HF_2 --> volume data entry was correct, SSD looks high but
# but is actually pretty similar to other SPIALB values:
spialb_subset <- subset(stem_density2, species == "SPIALB")

# There appears to be a high-weight outlier on the plot - I don't think this is a
# mistake, it was just a heavy, high volume stick

# Another note: there are both "BETPOPX_X" and "BETPOP_XX" in plot 3
# Should this only be one entry? I am not sure so I will leave it alone for now.

# Outlier at the bottom left (4_MYRGAL_GR_8D): data entry volume was correct, it's possible that weight was 
# incorrectly entered
myrgal_subset <- subset(stem_density, species == "MYRGAL")

#Plotting stem density once again
stem_density2 <- subset(stem_density, ssd != "NA")
ssd_2022 <- ggplot(data = stem_density2, aes(x= Weight, y= ssd, color = species))+
  geom_point()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SSD Scatterplot (2022 WildHell Dataset)")
ssd_2022

# pdf(file = "StemDensity_Scatterplot_2022.pdf", width = 8, height = 6)
# plot(ssd_2022)
# dev.off()

# SSD distribution

SSD_dist_2022 <- ggplot(data = stem_density2, aes(x = ssd, fill = species))+
  geom_histogram()+
  facet_wrap(~species)+theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank())+
  ggtitle("SSD Distribution (2022 WildHell Dataset)")
SSD_dist_2022

dist_1 <- ggplot( data = wildhell_2022_updated, aes(x = SLA, fill = species))+
  geom_histogram()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SLA Distribution (2022 WildHell Dataset)")+
  scale_fill_manual(values = colour_palette)

dist_2 <- ggplot( data = wildhell_2022_updated, aes(x = LDMC, fill = species))+
  geom_histogram()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("LDMC Distribution (2022 WildHell Dataset)")+
  scale_fill_manual(values = colour_palette)


dist_3 <- ggplot( data = stem_density2, aes(x = ssd, fill = species))+
  geom_histogram()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SSD Distribution (2022 WildHell Dataset)")+
  scale_fill_manual(values = colour_palette)
dist_3

Comparison_Histogram <- grid.arrange(dist_1, dist_2, dist_3, nrow = 1)

# pdf(file = "2022_Trait_Histogram.pdf", width = 12, height = 6)
# plot(Comparison_Histogram)
# dev.off()

SLA_Boxplot <- ggplot(wildhell_2022_updated, aes(x=species, y=SLA, fill=species)) + 
  geom_boxplot()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SLA Distribution (2022 WildHell Dataset)")+
  stat_boxplot(geom ='errorbar')
SLA_Boxplot

LDMC_Boxplot <- ggplot(wildhell_2022_updated, aes(x=species, y=LDMC, fill=species)) + 
  geom_boxplot()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("LDMC Distribution (2022 WildHell Dataset)")+
  stat_boxplot(geom ='errorbar')
LDMC_Boxplot

SSD_Boxplot <- ggplot(stem_density2, aes(x=species, y=ssd, fill=species)) + 
  geom_boxplot()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  ggtitle("SSD Distribution (2022 WildHell Dataset)")+
  stat_boxplot(geom ='errorbar')
SSD_Boxplot

# pdf(file = "2022_SLA_Boxplot.pdf", width = 12, height = 6)
# plot(SLA_Boxplot)
#   dev.off()

# pdf(file = "2022_LDMC_Boxplot.pdf", width = 12, height = 6)
# plot(LDMC_Boxplot)
#   dev.off()
# 
# pdf(file = "2022_SSD_Boxplot.pdf", width = 12, height = 6)
# plot(SSD_Boxplot)
# dev.off()

# Comparing 2019 vs. 2022 traits with boxplots:
wildhell_2022_updated$species_year <- paste(wildhell_2022_updated$species, " (2022)", sep = "")
wildhell_2019$species_year <- paste(wildhell_2019$species, " (2019)", sep =  "")

boxplot_ALNINC <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "ALNINC",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "ALNINC",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_ALNINC

boxplot_BETALL <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "BETALL",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "BETALL",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_BETALL

boxplot_BETPAP <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "BETPAP",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "BETPAP",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_BETPAP

boxplot_BETPOP <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "BETPOP",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "BETPOP",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_BETPOP

boxplot_AMECAN <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "AMECAN",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "AMECAN",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_AMECAN

boxplot_VIBCAS <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "VIBCAS",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "VIBCAS",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_VIBCAS

boxplot_SPIALB <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "SPIALB",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "SPIALB",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_SPIALB

boxplot_MYRGAL <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "MYRGAL",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "MYRGAL",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_MYRGAL

boxplot_AROMEL <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "AROMEL",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "AROMEL",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_AROMEL

boxplot_SPITOM <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "SPITOM",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "SPITOM",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_SPITOM

boxplot_DIELON <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "DIELON",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "DIELON",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_DIELON

boxplot_ACEPEN <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "ACEPEN",],
               mapping = aes(species_year, y = SLA, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "ACEPEN",],
               mapping = aes(x=species_year, y=SLA, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot_ACEPEN

Species_SLA_Boxplot <- grid.arrange(boxplot_ALNINC, boxplot_ACEPEN, boxplot_AMECAN, 
                                    boxplot_AROMEL, boxplot_BETALL, boxplot_BETPAP,
                                    boxplot_BETPOP, boxplot_DIELON, boxplot_MYRGAL,
                                    boxplot_SPIALB, boxplot_SPITOM, boxplot_VIBCAS,
                                    nrow = 3)



boxplot2_ALNINC <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "ALNINC",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "ALNINC",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_ALNINC

boxplot2_BETALL <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "BETALL",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "BETALL",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_BETALL

boxplot2_BETPAP <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "BETPAP",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "BETPAP",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_BETPAP

boxplot2_BETPOP <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "BETPOP",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "BETPOP",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_BETPOP

boxplot2_AMECAN <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "AMECAN",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "AMECAN",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_AMECAN

boxplot2_VIBCAS <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "VIBCAS",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "VIBCAS",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_VIBCAS

boxplot2_SPIALB <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "SPIALB",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "SPIALB",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_SPIALB

boxplot2_MYRGAL <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "MYRGAL",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "MYRGAL",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_MYRGAL

boxplot2_AROMEL <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "AROMEL",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "AROMEL",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_AROMEL

boxplot2_SPITOM <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "SPITOM",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "SPITOM",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_SPITOM

boxplot2_DIELON <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "DIELON",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "DIELON",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_DIELON

boxplot2_ACEPEN <- ggplot() + 
  geom_boxplot(wildhell_2019[wildhell_2019$species %in% "ACEPEN",],
               mapping = aes(species_year, y = LDMC, fill = species_year))+
  geom_boxplot(wildhell_2022_updated[wildhell_2022_updated$species %in% "ACEPEN",],
               mapping = aes(x=species_year, y=LDMC, fill = species_year))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_boxplot(geom ='errorbar')+
  labs(x = "Species + Year", color = "Species")+
  theme(legend.position = "none")
boxplot2_ACEPEN

Species_LDMC_Boxplot <- grid.arrange(boxplot2_ALNINC, boxplot2_ACEPEN, boxplot2_AMECAN, 
                                    boxplot2_AROMEL, boxplot2_BETALL, boxplot2_BETPAP,
                                    boxplot2_BETPOP, boxplot2_DIELON, boxplot2_MYRGAL,
                                    boxplot2_SPIALB, boxplot2_SPITOM, boxplot2_VIBCAS,
                                    nrow = 3)
Species_LDMC_Boxplot

# pdf(file = "LDMC_Boxplot_2019vs2022.pdf", width = 12, height = 6)
# plot(Species_SLA_Boxplot)
# dev.off()
# 
# pdf(file = "SLA_Boxplot_2019vs2022.pdf", width = 12, height = 6)
# plot(Species_LDMC_Boxplot)
# dev.off()

## Save the data output
names(wildhell_2019)
wildhell_2019_sub <- wildhell_2019[,c("id","species","site","ind","leaf", "SLA", "LDMC")]
wildhell_2019_sub$year <- "2019"

names(wildhell_2022)
wildhell_2022_sub <- wildhell_2022[,c("id","species","site","ind","leaf", "SLA", "LDMC")]
colnames(wildhell_2022)[colnames(wildhell_2022) == "Plot"] <- "plot"
wildhell_2022_sub$year <- "2022"

leafT <- rbind(wildhell_2019_sub, wildhell_2022_sub)

write.csv(leafT, "analyses/input/SLA_LDMC_trait.csv", row.names = F)

write.csv(stem_density, "analyses/input/SSD_trait.csv", row.names = F)
