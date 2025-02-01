# Bayesian Class
# CRD on 22 January 2025

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Load library 
library(rstanarm)
library(ggplot2)
library(arm)
library(readxl)
# wd
setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/data/TransferFromMainRepo")

# read it
d<-read_excel("Seedlings to Raised Beds 2016.edit.xlsx")

# get only data of interest
dcut <- subset(d, ID%in%c("ALNINC","QUERUB", "BETPOP", "BETPAP", "BETALL"))




