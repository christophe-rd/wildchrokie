#### Let's make a new datasheet for the Common Garden
## 27 March 2019 - Cat


# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
graphics.off()

# Set Working Directory
setwd("~/Documents/git/wildhellgarden/analyses/2019_data")
cg <-read.csv("input/2019_CG_datasheet.csv", header=TRUE)


cg <- cg[rep(seq_len(nrow(cg)), each=6),]
cg$Phase <- rep(c("Leaves", "Flowers", "Fruits", "Fall Colors", "Senescence", "Budset"), 340)

write.csv(cg, file="2019_CG_datasheet.csv", row.names = FALSE)
