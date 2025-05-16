# started May 3, 2023 by DL

# Aim of this code is to extract the values to reference in the output document. 
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)

library(gridExtra)
library(ggplot2)
library(viridis)
library(rstan)
#library(ggdist)
require(cowplot)

# if(length(grep("deirdreloughnan", getwd())>0)) {
#   setwd("~/Documents/github/wildhellgarden/analyses")
# } else if(length(grep("Lizzie", getwd())>0)) {
#   setwd("~/Documents/git/projects/others/deirdre/synchrony")
# } else{
#   setwd("/home/deirdre/WildHellCommonGarden") # for midge
# }

# Height:
pop <- c("GR", "HF", "SH", "WM")
spp <- c("ACEPEN", "ACESPI", "ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "QUEALB", "QUERUB", "SAMRAC", "SORAME", "SPIALB", "SPITOM", "VACMYR", "VIBCAS")

sppHt14 <- c("ACEPEN", "ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "SAMRAC", "SORAME", "SPIALB", "SPITOM", "VIBCAS")

# HEIGHT
load("..//analyses/output/mdlHtFullTrait.Rda")
sumHtFull <- summary(mdlHtFull)$summary

mu_grandHt <-sumHtFull[grep("mu_grand", rownames(sumHtFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
mu_a_spHt <-sumHtFull[grep("mu_sp", rownames(sumHtFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_spHt) <- spp
mu_a_popHt <-sumHtFull[grep("mu_pop", rownames(sumHtFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_popHt) <- pop

# sigma_a_spHt <- sumHtFull[grep("sigma_sp", rownames(sumHtFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_a_popHt <- sumHtFull[grep("sigma_pop", rownames(sumHtFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_yHt <- sumHtFull[grep("sigma_y", rownames(sumHtFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]

sigSpHt <- round(sumHtFull[grep("sigma_sp", rownames(sumHtFull)), c("mean")],2)
sigSpHt25 <- round(sumHtFull[grep("sigma_sp", rownames(sumHtFull)), c("25%")],2)
sigSpHt75 <- round(sumHtFull[grep("sigma_sp", rownames(sumHtFull)), c("75%")],2)

sigPopHt <- round(sumHtFull[grep("sigma_pop", rownames(sumHtFull)), c("mean")],2)
sigPopHt25 <- round(sumHtFull[grep("sigma_pop", rownames(sumHtFull)), c("25%")],2)
sigPopHt75 <- round(sumHtFull[grep("sigma_pop", rownames(sumHtFull)), c("75%")],2)


grPopMHt <- round(mu_a_popHt["GR","mean"],2)
hfPopMHt <- round(mu_a_popHt["HF","mean"],2)
shPopMHt <- round(mu_a_popHt["SH","mean"],2)
wmPopMHt <- round(mu_a_popHt["WM","mean"],2)

grPop25Ht <- round(mu_a_popHt["GR","25%"],2)
hfPop25Ht <- round(mu_a_popHt["HF","25%"],2)
shPop25Ht <- round(mu_a_popHt["SH","25%"],2)
wmPop25Ht <- round(mu_a_popHt["WM","25%"],2)

grPop75Ht <- round(mu_a_popHt["GR","75%"],2)
hfPop75Ht <- round(mu_a_popHt["HF","75%"],2)
shPop75Ht <- round(mu_a_popHt["SH","75%"],2)
wmPop75Ht <- round(mu_a_popHt["WM","75%"],2)

# SLA
load("..//analyses/output/mdlSLAFull.Rda")
sumSLAFull <- summary(mdlSLAFull)$summary

mu_grandSLA <-sumSLAFull[grep("mu_grand", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
mu_a_spSLA <-sumSLAFull[grep("mu_sp", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_spSLA) <- spp
mu_a_popSLA <-sumSLAFull[grep("mu_pop", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_popSLA) <- pop
# 
# sigma_a_sp <- sumSLAFull[grep("sigma_sp", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_a_pop <- sumSLAFull[grep("sigma_pop", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_y <- sumSLAFull[grep("sigma_y", rownames(sumSLAFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
sigSpSLA <- round(sumSLAFull[grep("sigma_sp", rownames(sumSLAFull)), c("mean")],2)
sigSpSLA25 <- round(sumSLAFull[grep("sigma_sp", rownames(sumSLAFull)), c("25%")],2)
sigSpSLA75 <- round(sumSLAFull[grep("sigma_sp", rownames(sumSLAFull)), c("75%")],2)

sigPopSLA <- round(sumSLAFull[grep("sigma_pop", rownames(sumSLAFull)), c("mean")],2)
sigPopSLA25 <- round(sumSLAFull[grep("sigma_pop", rownames(sumSLAFull)), c("25%")],2)
sigPopSLA75 <- round(sumSLAFull[grep("sigma_pop", rownames(sumSLAFull)), c("75%")],2)

grPopMSLA <- round(mu_a_popSLA["GR","mean"],2)
hfPopMSLA <- round(mu_a_popSLA["HF","mean"],2)
shPopMSLA <- round(mu_a_popSLA["SH","mean"],2)
wmPopMSLA <- round(mu_a_popSLA["WM","mean"],2)

grPop25SLA <- round(mu_a_popSLA["GR","25%"],2)
hfPop25SLA <- round(mu_a_popSLA["HF","25%"],2)
shPop25SLA <- round(mu_a_popSLA["SH","25%"],2)
wmPop25SLA <- round(mu_a_popSLA["WM","25%"],2)

grPop75SLA <- round(mu_a_popSLA["GR","75%"],2)
hfPop75SLA <- round(mu_a_popSLA["HF","75%"],2)
shPop75SLA <- round(mu_a_popSLA["SH","75%"],2)
wmPop75SLA <- round(mu_a_popSLA["WM","75%"],2)
# LDMC

load("..//analyses/output/mdlLDMCFull.Rda")
sumLDMCFull <- summary(mdlLDMC)$summary

mu_grandLDMC <-sumLDMCFull[grep("mu_grand", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
mu_a_spLDMC <-sumLDMCFull[grep("mu_sp", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_spLDMC) <- spp
mu_a_popLDMC <-sumLDMCFull[grep("mu_pop", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_popLDMC) <- pop
# 
# sigma_a_sp <- sumLDMCFull[grep("sigma_sp", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_a_pop <- sumLDMCFull[grep("sigma_pop", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_y <- sumLDMCFull[grep("sigma_y", rownames(sumLDMCFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]

sigSpLDMC <- round(sumLDMCFull[grep("sigma_sp", rownames(sumLDMCFull)), c("mean")],2)
sigSpLDMC25 <- round(sumLDMCFull[grep("sigma_sp", rownames(sumLDMCFull)), c("25%")],2)
sigSpLDMC75 <- round(sumLDMCFull[grep("sigma_sp", rownames(sumLDMCFull)), c("75%")],2)

sigPopLDMC <- round(sumLDMCFull[grep("sigma_pop", rownames(sumLDMCFull)), c("mean")],2)
sigPopLDMC25 <- round(sumLDMCFull[grep("sigma_pop", rownames(sumLDMCFull)), c("25%")],2)
sigPopLDMC75 <- round(sumLDMCFull[grep("sigma_pop", rownames(sumLDMCFull)), c("75%")],2)


grPopMLDMC <- round(mu_a_popLDMC["GR","mean"],2)
hfPopMLDMC <- round(mu_a_popLDMC["HF","mean"],2)
shPopMLDMC <- round(mu_a_popLDMC["SH","mean"],2)
wmPopMLDMC <- round(mu_a_popLDMC["WM","mean"],2)

grPop25LDMC <- round(mu_a_popLDMC["GR","25%"],2)
hfPop25LDMC <- round(mu_a_popLDMC["HF","25%"],2)
shPop25LDMC <- round(mu_a_popLDMC["SH","25%"],2)
wmPop25LDMC <- round(mu_a_popLDMC["WM","25%"],2)

grPop75LDMC <- round(mu_a_popLDMC["GR","75%"],2)
hfPop75LDMC <- round(mu_a_popLDMC["HF","75%"],2)
shPop75LDMC <- round(mu_a_popLDMC["SH","75%"],2)
wmPop75LDMC <- round(mu_a_popLDMC["WM","75%"],2)
# SSD
load("..//analyses/output/mdlSSD2022.Rda")

sumSSDFull <- summary(mdlSSD)$summary
sppSSD12 <- c("ALNINC", "AMECAN", "AROMEL", "BETALL", "BETPAP", "BETPOP", "DIELON", "MYRGAL", "SAMRAC", "SPIALB", "SPITOM", "VIBCAS")

mu_grandSSD <-sumSSDFull[grep("mu_grand", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
mu_a_spSSD <-sumSSDFull[grep("mu_sp", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_spSSD) <- sppSSD12
mu_a_popSSD <-sumSSDFull[grep("mu_pop", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]; rownames(mu_a_popSSD) <- pop

# sigma_a_sp <- sumSSDFull[grep("sigma_sp", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_a_pop <- sumSSDFull[grep("sigma_pop", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]
# sigma_y <- sumSSDFull[grep("sigma_y", rownames(sumSSDFull)), c("mean", "2.5%", "97.5%","25%","75%", "n_eff", "Rhat")]

sigSpSSD <- format(round(sumSSDFull[grep("sigma_sp", rownames(sumSSDFull)), c("mean")],2), nsmall = 2)
sigSpSSD25 <- format(round(sumSSDFull[grep("sigma_sp", rownames(sumSSDFull)), c("25%")],2), nsmall = 2)
sigSpSSD75 <- format(round(sumSSDFull[grep("sigma_sp", rownames(sumSSDFull)), c("75%")],2), nsmall = 2)

sigPopSSD <- format(round(sumSSDFull[grep("sigma_pop", rownames(sumSSDFull)), c("mean")],2), nsmall = 2)
sigPopSSD25 <- format(round(sumSSDFull[grep("sigma_pop", rownames(sumSSDFull)), c("25%")],2), nsmall = 2)
sigPopSSD75 <- format(round(sumSSDFull[grep("sigma_pop", rownames(sumSSDFull)), c("75%")],2), nsmall = 2)


grPopMSSD <- format(round(mu_a_popSSD["GR","mean"],2), nsmall = 2)
hfPopMSSD <- format(round(mu_a_popSSD["HF","mean"],2), nsmall = 2)
shPopMSSD <- format(round(mu_a_popSSD["SH","mean"],2), nsmall = 2)
wmPopMSSD <- format(round(mu_a_popSSD["WM","mean"],2), nsmall = 2)

grPop25SSD <- format(round(mu_a_popSSD["GR","25%"],2), nsmall = 2)
hfPop25SSD <- format(round(mu_a_popSSD["HF","25%"],2), nsmall = 2)
shPop25SSD <- format(round(mu_a_popSSD["SH","25%"],2), nsmall = 2)
wmPop25SSD <- format(round(mu_a_popSSD["WM","25%"],2), nsmall = 2)

grPop75SSD <- format(round(mu_a_popSSD["GR","75%"],2), nsmall = 2)
hfPop75SSD <- format(round(mu_a_popSSD["HF","75%"],2), nsmall = 2)
shPop75SSD <- format(round(mu_a_popSSD["SH","75%"],2), nsmall = 2)
wmPop75SSD <- format(round(mu_a_popSSD["WM","75%"],2), nsmall = 2)
