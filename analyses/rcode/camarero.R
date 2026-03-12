setwd("/Users/christophe_rouleau-desrochers/Downloads/chronicle-of-nature-calendar.v1.0.5")

p <- read.csv("phenology.csv")

# which phenoevents are present in this dataset:
message(print(unique(p$eventtype)))

spp <- c("Pinus halepensis", "Pinus sylvestris")




p$taxon[grepl("halepensis", p$taxon)]

p2 <- subset(p, taxon %in% spp)

p3 <- subset(p2, dataset %in% p2$dataset[grepl("Stolby", p2$dataset)])

unique(p3$eventtype)
View(p3)
