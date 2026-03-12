setwd("/Users/christophe_rouleau-desrochers/Downloads/chronicle-of-nature-calendar.v1.0.5")

p <- read.csv("phenology.csv")

# which phenoevents are present in this dataset:


p2 <- subset(p, taxon %in% "Pinus sylvestris" & dataset %in% dataset[grepl("Stolby", p$dataset)])
print(unique(p2$eventtype))
message("Number years with of budburst observations:", 
        length(unique(p2$year[which(p2$eventtype %in% "onset of budburst")])))
message("Number years with of leaf colouring observations:", 
        length(unique(p2$year[which(p2$eventtype %in% "onset of autumn colouring")])))

length(unique (p2$))
p3 <- subset(p2, dataset %in% )

unique(p3$eventtype)
View(p3)
