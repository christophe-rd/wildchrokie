# Checks with new generated quantities 
# Treeid gdd####
df_fitgdd <- as.data.frame(fitgdd)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### extract from stan #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
fullintercept_cols <- grep("^fullintercept", colnames(df_fitgdd), value = TRUE)
fullintercept <- df_fitgdd[, fullintercept_cols]
colnames(fullintercept) <- 1:ncol(fullintercept)

treeid_slope_cols <- grep("^treeid_slope", colnames(df_fitgdd), value = TRUE)
treeid_bspp <- df_fitgdd[, treeid_slope_cols] / gddscale
colnames(treeid_bspp) <- 1:ncol(treeid_bspp)

# sim each tree id for each gddseq n iterations times
y_post_array <- extract(fitgdd, "y_post")$y_post
# dimensions: [n_draws, Ngddseq, Ntreeid]



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### old way #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# full posterior
columns <- colnames(df_fitgdd)[!grepl("prior", colnames(df_fitgdd))]
sigma_df <- df_fitgdd[, columns[grepl("sigma", columns)]]
bspp_df <- df_fitgdd[, columns[grepl("bsp", columns)]]
treeid_df <- df_fitgdd[, grepl("treeid", columns) & !grepl("z|sigma", columns)]
aspp_df <- df_fitgdd[, columns[grepl("aspp", columns)]]
site_df <- df_fitgdd[, columns[grepl("asite", columns)]]

# change colnames
colnames(bspp_df) <- 1:ncol(bspp_df)
colnames(treeid_df) <- 1:ncol(treeid_df)
colnames(aspp_df) <- 1:ncol(aspp_df)
colnames(site_df) <- 1:ncol(site_df)


subyvec <- vector()
for (i in 1:length(unique(emp$treeid_num))) {
  subyvec[i] <- paste("atreeid", "[",i,"]", sep = "")  
}
atreeid_gdd <- subset(df_fitgdd, select = subyvec)
colnames(atreeid_gdd) <- 1:length(subyvec)

# the spp values for each tree id
treeid_aspp_gdd <- data.frame(matrix(ncol = ncol(atreeid_gdd), nrow = nrow(df_fitgdd)))
colnames(treeid_aspp_gdd) <- colnames(atreeid_gdd)

for (i in seq_len(ncol(treeid_aspp_gdd))) { # i = 1
  tree_id <- as.integer(colnames(treeid_aspp_gdd)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_aspp_gdd[, i] <- aspp_df[, spp_id]
}
treeid_aspp_gdd

# the site values for each tree id
treeid_asite_gdd <- data.frame(matrix(ncol = ncol(atreeid_gdd), nrow = nrow(df_fitgdd)))
colnames(treeid_asite_gdd) <- colnames(atreeid_gdd)

for (i in seq_len(ncol(treeid_asite_gdd))) { # i = 1
  tree_id <- as.integer(colnames(treeid_asite_gdd)[i])
  site_id <- treeid_spp_site$site_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_asite_gdd[, i] <- site_df[, site_id]
}
treeid_asite_gdd

# recover a
treeid_a <- data.frame(matrix(ncol = ncol(atreeid_gdd), nrow = nrow(df_fitgdd)))
colnames(treeid_a) <- colnames(atreeid_gdd)

for (i in seq_len(ncol(treeid_a))) { # i = 1
  treeid_a[, i] <- df_fitgdd[, "a"]
}

# sum all 3 dfs together to get the full intercept for each treeid
fullintercept_old <-
  treeid_a + 
  atreeid_gdd +
  treeid_aspp_gdd +
  treeid_asite_gdd
fullintercept_old

# tree id slopes
bspp_df4 <- bspp_df
for (i in 1:ncol(bspp_df4)){
  bspp_df4[[i]] <- bspp_df4[[i]] #/ gddscale
}

# now get the slope for each treeid
treeid_bspp_old <- data.frame(matrix(ncol = ncol(atreeid_gdd), nrow = nrow(df_fitgdd)))
colnames(treeid_bspp_old) <- colnames(atreeid_gdd)

for (i in seq_len(ncol(treeid_bspp_old))) { # i = 30
  tree_id <- as.integer(colnames(treeid_bspp_old)[i])
  spp_id <- treeid_spp_site$spp_num[match(tree_id, treeid_spp_site$treeid_num)]
  treeid_bspp_old[, i] <- bspp_df4[, spp_id]
}
treeid_bspp_old

y_post_list <- list()  # store posterior predictions in a list where each tree id gets matrixad

for (i in seq_along(treeidvecnum)) { # i = 1
  tree_col <- as.character(treeidvecnum[i]) 
  
  y_post <- sapply(1:nrow(df_fitgdd), function(f) {
    rnorm(length(gddseq), 
          fullintercept[f, tree_col] + treeid_bspp_old[f, tree_col] * gddseq,
          sigma_df$sigma_y[f])
  })
  y_post_list[[tree_col]] <- y_post
}


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
##### 1:1 line plot #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# the following should be exactly on the same line
fullintercept_old_means <- colMeans(fullintercept_old)
fullintercept_new_means <- colMeans(fullintercept)

plot(fullintercept_old_means, fullintercept_new_means,
     pch = 1, col = "black",
     xlab = "Old fullintercept", ylab = "New fullintercept",
     main = "1:1 comparison")
abline(0, 1, col = "red", lty = 2)

plot(fullintercept_old$`1`, fullintercept_old$`1`,
     pch = 1, col = "black",
     xlab = "Old fullintercept", ylab = "New fullintercept",
     main = "1:1 comparison")
abline(0, 1, col = "red", lty = 2)

bspp_means <- colMeans(treeid_bspp)
bspp_means_old <- colMeans(treeid_bspp_old)

plot(bspp_means_old, bspp_means,
     pch = 1, col = "black",
     xlab = "Old fullintercept", ylab = "New fullintercept",
     main = "1:1 comparison")
abline(0, 1, col = "red", lty = 2)

treeidvecnum <- 1:ncol(fullintercept)
treeidvecname <- treeid_spp_site$treeid
gddseq <- seq(min(emp$pgsGDD5), max(emp$pgsGDD5), length.out = 100)

# sim at each gddseq and each
# get means across draws for each gddseq point and treeid
y_post_new_means <- apply(y_post_array, c(2, 3), mean)
# dimensions: [Ngddseq, Ntreeid]

# get means from old list
y_post_old_means <- sapply(as.character(treeidvecnum), function(tree_col) {
  rowMeans(y_post_list[[tree_col]])
})
plot(y_post_old_means, y_post_new_means,
     pch = 16, cex = 0.4, col = adjustcolor("steelblue", alpha.f = 0.3),
     xlab = "Old y_post means", ylab = "New y_post means",
     main = "1:1 check: per-tree posterior predictions")
abline(0, 1, col = "red", lty = 2)