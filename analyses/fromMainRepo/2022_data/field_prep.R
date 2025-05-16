# Started March 25, 2022 by DL

require(dplyr)
# This script is to help me understand how many individuals of each species is alive:
setwd("~/Documents/github/wildhellgarden/analyses/2022_data")

ht <- read.csv("..//2020_data/heights.csv", header=T, na.strings=c("","NA"))

ht$temp <- ht$ID
temp <- str_split_fixed(ht$temp, "_", 3)
ht$speciesName <- temp[,1]
ht$site <- temp[,2]
ht$indiv <- temp[,3]

head(ht)

dead <- c("dead", "Dead","dead sampled?", "dead)","Dead)")
rip <- ht[ht$Notes %in% dead,]

alive <- ht[!ht$Notes %in% dead,] #276

alive$count <- 1
noSp <- aggregate(alive["ID"], alive["speciesName"], FUN = length)

sum(noSp$ID) # 276
# no longer have any VACMR or VACMYR - all dead

# Really don't have enough to sample: AcePen, AceSpi, BetPopx, QueAlb, QueRub, assuming our threshold for measuring leaf
#traits is say 5 individuals - leave 267 individuals to sample

# ACEPEN	2
#	ACESPI	3
# ALNINC	35
# AMECAN	5
# AROMEL	13
# BETALL	23
# BETPAP	11
# BETPOP	25
# BETPOPX	1
# DIELON	29
# MYRGAL	14
# QUEALB	1
# QUERUB	2
# SAMRAC	15
# SORAME	7
# SPIALB	35
# SPITOM	28
# VIBCAS	27


