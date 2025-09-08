# Posterior recovery tests
# 8 September 2025
# extracting posterior distribution without lme4 functions

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("deirdre", getwd()) > 0)) {
  setwd("~/Documents/github/egret/analyses")
} else if(length(grep("lizzie", getwd()) > 0)) {
  setwd("/Users/lizzie/Documents/git/projects/egret/analyses")
} else if(length(grep("christophe_rouleau-desrochers", getwd())) > 0){
  setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")
} 

fitnested <- readRDS("output/fitnested")

df_fit <- as.data.frame(fitnested)

# recover slope
colnames(df_fit)
# grab ids nested in spp
ids_cols <- colnames(df_fit)[grepl("ids:spp:", colnames(df_fit))]
ids_cols <- ids_cols[1:length(ids_cols)-1]
ids_df <- df_fit[, colnames(df_fit) %in% ids_cols]
# change their names
colnames(ids_df) <- sub(".*ids:spp:(.*)\\]$", "\\1", colnames(ids_df))
# empty ids dataframe
ids_df2 <- data.frame(
  ids_spp = character(ncol(ids_df)),
  fit_a_ids_spp = numeric(ncol(ids_df)),  
  fit_per5 = NA, 
  fit_per95 = NA,
  fit_sd = NA
)
for (i in 1:ncol(ids_df)) { # i = 1
  ids_df2$ids_spp[i] <- colnames(ids_df)[i]         
  ids_df2$fit_a_ids_spp[i] <- round(mean(ids_df[[i]]),3)  
  ids_df2$fit_per5[i] <- round(quantile(ids_df[[i]], probs = 0.055), 3)
  ids_df2$fit_per95[i] <- round(quantile(ids_df[[i]], probs = 0.945), 3)
  ids_df2$fit_sd[i] <- round(sd(ids_df[[i]]), 3)
}
ids_df2

