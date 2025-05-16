# Started April 2023 by DL

# aim of this code is to plot the model output of the wildhell common garden trait models:
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)

library(gridExtra)
library(ggplot2)
library(viridis)
library(rstan)
library(ggdist)
require(cowplot)

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/wildhellgarden/analyses")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/projects/others/deirdre/synchrony")
} else{
  setwd("/home/deirdre/WildHellCommonGarden") # for midge
}

# Height:
pop <- c("GR", "HF", "SH", "WM")
spp <- c("ACEPEN", "ACESPI", "ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "QUEALB", "QUERUB", "SAMRAC", "SORAME", "SPIALB", "SPITOM", "VACMYR", "VIBCAS")

sppHt14 <- c("ACEPEN", "ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "SAMRAC", "SORAME", "SPIALB", "SPITOM", "VIBCAS")

load("output/mdlHtFullTrait.Rda")
sumHtFull <- summary(mdlHtFull)$summary

mu_grand <-sumHtFull[grep("mu_grand", rownames(sumHtFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumHtFull[grep("mu_sp", rownames(sumHtFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- spp
mu_a_pop <-sumHtFull[grep("mu_pop", rownames(sumHtFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumHtFull[grep("sigma_sp", rownames(sumHtFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumHtFull[grep("sigma_pop", rownames(sumHtFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumHtFull[grep("sigma_y", rownames(sumHtFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

HtFull <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
HtFull$rowName <- rownames(HtFull)

# 14 spp only
load("output/mdlHtFull14.Rda")

sumHtFull14 <- summary(mdlHtFull14)$summary

mu_grand <-sumHtFull14[grep("mu_grand", rownames(sumHtFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumHtFull14[grep("mu_sp", rownames(sumHtFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- sppHt14
mu_a_pop <-sumHtFull14[grep("mu_pop", rownames(sumHtFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumHtFull14[grep("sigma_sp", rownames(sumHtFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumHtFull14[grep("sigma_pop", rownames(sumHtFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumHtFull14[grep("sigma_y", rownames(sumHtFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

HtFull14 <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
HtFull14$rowName <- rownames(HtFull14)

# 1. Compare that the two models are not that different
HtDiff <- merge(HtFull, HtFull14, all = T, by = "rowName")
HtDiff$diff <- HtDiff$mean.x - HtDiff$mean.y

# there are some pretty big diff in spp estimates, even some that switch signs, but a lot of the 14 spp model estimates cross zero, but most of the removed spp are trees (Acer and three Quercus, one is a shrub). Of the 7 spp that have esti switch sign, most are shrubs (2 are trees). But the grand mean, pop esti, and sigma are all very similar


# 2. Do we see populations diff?

postHt <- rstan::extract(mdlHtFull)

HtPop <- data.frame(postHt$mu_pop)
names(HtPop) <- pop

HtLong <- melt(HtPop)
popOrder <- c("HF","GR","WM","SH")

HtPopPlot <- ggplot() +
  stat_eye(data = HtLong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species-level height intercept", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 1, label = "a)", cex =5)

# 3. Do we see clear spp diff? Shrub vs tree?
HtSp <- data.frame(postHt$mu_sp)
names(HtSp) <- spp

HtSpLong <- melt(HtSp)

HtSpPlot <- ggplot() +
  stat_eye(data = HtSpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Height effect", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 1, label = "a)", cex =5)


# SLA:
load("output/mdlSLAFull.Rda")
sumSLAFull <- summary(mdlSLAFull)$summary

mu_grand <-sumSLAFull[grep("mu_grand", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumSLAFull[grep("mu_sp", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- spp
mu_a_pop <-sumSLAFull[grep("mu_pop", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumSLAFull[grep("sigma_sp", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumSLAFull[grep("sigma_pop", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumSLAFull[grep("sigma_y", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

SLAFull <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
SLAFull$rowName <- rownames(SLAFull)

# 14 spp only
load("output/mdlSLAFull14.Rda")

sumSLAFull14 <- summary(mdlSLAFull14)$summary

sppSLA13 <- c("ACEPEN", "ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "SAMRAC", "SPIALB", "SPITOM", "VIBCAS")
mu_grand <-sumSLAFull14[grep("mu_grand", rownames(sumSLAFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumSLAFull14[grep("mu_sp", rownames(sumSLAFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- sppSLA13
mu_a_pop <-sumSLAFull14[grep("mu_pop", rownames(sumSLAFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumSLAFull14[grep("sigma_sp", rownames(sumSLAFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumSLAFull14[grep("sigma_pop", rownames(sumSLAFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumSLAFull14[grep("sigma_y", rownames(sumSLAFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

SLAFull14 <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
SLAFull14$rowName <- rownames(SLAFull14)

# 1. Compare that the two models are not that different
SLADiff <- merge(SLAFull, SLAFull14, all = T, by = "rowName")
SLADiff$diff <- SLADiff$mean.x - SLADiff$mean.y

# there are some pretty big diff in spp estimates, even some that switch signs, but most of the removed spp are trees (Acer and three Quercus, one is a shrub). Of the 7 spp that have esti switch sign, most are shrubs (2 are trees). But the grand mean, pop esti, and sigma are all very similar


# 2. Do we see populations diff?

postSLA <- rstan::extract(mdlSLAFull)

SLAPop <- data.frame(postSLA$mu_pop)
names(SLAPop) <- pop

SLALong <- melt(SLAPop)
popOrder <- c("HF","GR","WM","SH")

SLAPopPlot <- ggplot() +
  stat_eye(data = SLALong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill ="#CC6677") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species-level SLA intercept", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "b)", cex =5)

# 3. Do we see clear spp diff? Shrub vs tree?
SLASp <- data.frame(postSLA$mu_sp)
names(SLASp) <- spp

SLASpLong <- melt(SLASp)

SLASpPlot <- ggplot() +
  stat_eye(data = SLASpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#CC6677") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species level variation in height", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 2, label = "b)", cex =5)

# LDMC
load("output/mdlLDMCFull.Rda")
sumLDMCFull <- summary(mdlLDMC)$summary

mu_grand <-sumLDMCFull[grep("mu_grand", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumLDMCFull[grep("mu_sp", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- spp
mu_a_pop <-sumLDMCFull[grep("mu_pop", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumLDMCFull[grep("sigma_sp", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumLDMCFull[grep("sigma_pop", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumLDMCFull[grep("sigma_y", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

LDMCFull <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
LDMCFull$rowName <- rownames(LDMCFull)

# 14 spp only
load("output/mdlLDMCFull14.Rda")

sumLDMCFull14 <- summary(mdlLDMCFull14)$summary

mu_grand <-sumLDMCFull14[grep("mu_grand", rownames(sumLDMCFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumLDMCFull14[grep("mu_sp", rownames(sumLDMCFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- sppSLA13
mu_a_pop <-sumLDMCFull14[grep("mu_pop", rownames(sumLDMCFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumLDMCFull14[grep("sigma_sp", rownames(sumLDMCFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumLDMCFull14[grep("sigma_pop", rownames(sumLDMCFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumLDMCFull14[grep("sigma_y", rownames(sumLDMCFull14)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

LDMCFull14 <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
LDMCFull14$rowName <- rownames(LDMCFull14)

# 1. Compare that the two models are not that different
LDMCDiff <- merge(LDMCFull, LDMCFull14, all = T, by = "rowName")
LDMCDiff$diff <- LDMCDiff$mean.x - LDMCDiff$mean.y

# there are some pretty big diff in spp estimates, even some that switch signs, but a lot of the 14 spp model estimates cross zero, but most of the removed spp are trees (Acer and three Quercus, one is a shrub). Of the 7 spp that have esti switch sign, most are shrubs (2 are trees). But the grand mean, pop esti, and sigma are all very similar


# 2. Do we see populations diff?

postLDMC <- rstan::extract(mdlLDMCFull)

LDMCPop <- data.frame(postLDMC$mu_pop)
names(LDMCPop) <- pop

LDMCLong <- melt(LDMCPop)
popOrder <- c("HF","GR","WM","SH")

LDMCPopPlot <- ggplot() +
  stat_eye(data = LDMCLong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "#44AA99") +
  theme_classic() +
  ylim(200,375) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "LDMC effect", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 375, label = "c)", cex =5)

# 3. Do we see clear spp diff? Shrub vs tree?
LDMCSp <- data.frame(postLDMC$mu_sp)
names(LDMCSp) <- spp

LDMCSpLong <- melt(LDMCSp)

LDMCSpPlot <- ggplot() +
  stat_eye(data = LDMCSpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#44AA99") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species level variation in height", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 2, label = "b)", cex =5)

# SSD
load("output/mdlSSD2022.Rda")

sumSSDFull <- summary(mdlSSD)$summary

mu_grand <-sumSSDFull[grep("mu_grand", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
mu_a_sp <-sumSSDFull[grep("mu_sp", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_sp) <- spp
mu_a_pop <-sumSSDFull[grep("mu_pop", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]; rownames(mu_a_pop) <- pop

sigma_a_sp <- sumSSDFull[grep("sigma_sp", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_a_pop <- sumSSDFull[grep("sigma_pop", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]
sigma_y <- sumSSDFull[grep("sigma_y", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%", "n_eff", "Rhat")]

SSDFull <- data.frame(rbind( mu_grand, mu_a_sp, mu_a_pop, sigma_a_sp, sigma_a_pop, sigma_y))
SSDFull$rowName <- rownames(SSDFull)

# 2. Do we see populations diff?

postSSD <- rstan::extract(mdlSSD)

SSDPop <- data.frame(postSSD$mu_pop)
names(SSDPop) <- pop

SSDLong <- melt(SSDPop)
popOrder <- c("HF","GR","WM","SH")

SSDPopPlot <- ggplot() +
  stat_eye(data = SSDLong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "purple4") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "SSD effect", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 0.2, label = "d)", cex =5)

# 3. Do we see clear spp diff? Shrub vs tree?
SSDSp <- data.frame(postSSD$mu_sp)

sppSSD12 <- c("ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "SAMRAC", "SPIALB", "SPITOM", "VIBCAS")
names(SSDSp) <- sppSSD12

SSDSpLong <- melt(SSDSp)

SSDSpPlot <- ggplot() +
  stat_eye(data = SSDSpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species level variation in height", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 2, label = "b)", cex =5)

# Pop plot:
HtPopPlot <- ggplot() +
  stat_eye(data = HtLong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species-level height intercept", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 2, label = "a)", cex =5)

SLAPopPlot <- ggplot() +
  stat_eye(data = SLALong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill ="#CC6677") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species-level SLA intercept", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 10, label = "b)", cex =5)

LDMCPopPlot <- ggplot() +
  stat_eye(data = LDMCLong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "#44AA99") +
  theme_classic() +
  ylim(200,375) +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species-level LDMC intercept", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 375, label = "c)", cex =5)

SSDPopPlot <- ggplot() +
  stat_eye(data = SSDLong, aes(x = factor(variable, level = popOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "purple4") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Population", y = "Species-level SSD intercept", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 0.1, label = "d)", cex =5)

# Plot the trait by pop figures
pdf("traitPop.pdf", height = 8, width = 8)
plot_grid(HtPopPlot,SLAPopPlot,LDMCPopPlot,SSDPopPlot, nrow = 2, ncol =2)
dev.off()

# Species plot

HtSpPlot <- ggplot() +
  stat_eye(data = HtSpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  theme_classic() +
  theme(axis.title.x= element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x= element_text(size = 12, angle = 78, 
                                  hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in height", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 2.5, label = "a)", cex =5)

SLASpPlot <- ggplot() +
  stat_eye(data = SLASpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#CC6677") +
  theme_classic() +
  theme(axis.title.x= element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x= element_text(size = 12, angle = 78, 
                                  hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in SLA", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 10, label = "b)", cex =5)

LDMCSpPlot <- ggplot() +
  stat_eye(data = LDMCSpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "#44AA99") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12, angle = 78, 
                                 hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in LDMC", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 150, label = "c)", cex =5)

SSDSpPlot <- ggplot() +
  stat_eye(data = SSDSpLong, aes(x = factor(variable), y = value), .width = c(.90, .5), cex = 0.75, fill = "purple4") +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12,angle = 78, 
                                 hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in SSD", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 0.2, label = "d)", cex =5)

pdf("traitSpecies.pdf", height = 15, width = 8)
plot_grid(HtSpPlot,SLASpPlot,LDMCSpPlot,SSDSpPlot, nrow = 4)
dev.off()

# Shrub vs tree:
spInfo <- read.csv("input/species_list.csv")

names(HtSpLong) <- c("species","value"); HtSpLong$species <- tolower(HtSpLong$species)
names(SLASpLong) <- c("species","value"); SLASpLong$species <- tolower(SLASpLong$species)
names(LDMCSpLong) <- c("species","value"); LDMCSpLong$species <- tolower(LDMCSpLong$species)
names(SSDSpLong) <- c("species","value"); SSDSpLong$species <- tolower(SSDSpLong$species)

HtSpLong <- merge(HtSpLong, spInfo, by = "species")
SLASpLong <- merge(SLASpLong, spInfo, by = "species")
LDMCSpLong <- merge(LDMCSpLong, spInfo, by = "species")
SSDSpLong <- merge(SSDSpLong, spInfo, by = "species")

# Species plot coloured by tree/shrub

HtSpPlot <- ggplot() +
  stat_eye(data = HtSpLong, aes(x = factor(species), y = value, fill = type), .width = c(.90, .5), cex = 0.75) +
  theme_classic() +
  theme(axis.title.x= element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x= element_text(size = 12, angle = 78, 
                                  hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in height", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 2.5, label = "a)", cex =5)

SLASpPlot <- ggplot() +
  stat_eye(data = SLASpLong, aes(x = factor(species), y = value, fill = type), .width = c(.90, .5), cex = 0.75) +
  theme_classic() +
  theme(axis.title.x= element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x= element_text(size = 12, angle = 78, 
                                  hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in SLA", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 10, label = "b)", cex =5)

LDMCSpPlot <- ggplot() +
  stat_eye(data = LDMCSpLong, aes(x = factor(species), y = value, fill = type), .width = c(.90, .5), cex = 0.75) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12, angle = 78, 
                                 hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in LDMC", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 150, label = "c)", cex =5)

SSDSpPlot <- ggplot() +
  stat_eye(data = SSDSpLong, aes(x = factor(species), y = value, fill = type), .width = c(.90, .5), cex = 0.75) +
  theme_classic() +
  theme(axis.title.x=element_text(size = 12),
        axis.title.y=element_text(size = 12),
        axis.text.x=element_text(size = 12,angle = 78, 
                                 hjust=1),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  labs( x = "Species", y = "Species level variation in SSD", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 0.2, label = "d)", cex =5)

pdf("traitSpeciesType.pdf", height = 15, width = 8)
plot_grid(HtSpPlot,SLASpPlot,LDMCSpPlot,SSDSpPlot, nrow = 4)
dev.off()