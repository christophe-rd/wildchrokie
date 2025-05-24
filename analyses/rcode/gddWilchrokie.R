## 22 May 2025
## CRD

### Goal is to calculate GDD using primary and full growing season and implementing wildhell garden gdd code


# libraries
#library(rstan)
library(dplyr)
library(brms)
library(tidybayes)
library(tidyr)

# Clear workspace
rm(list=ls()) # remove everything currently held in the R memory
options(stringsAsFactors=FALSE)
options(mc.cores = parallel::detectCores())
graphics.off()

setwd("/Users/christophe_rouleau-desrochers/github/wildhellgarden/analyses/")
