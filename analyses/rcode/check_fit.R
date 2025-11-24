rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

file <- 'output/stanOutput/fit'
con <- file(file, 'r')

names <- readLines(con, n = 1)
names
names <- strsplit(names, '\",\"')[[1]]

asp <- grep('^asp\\[', names, value = TRUE)
asite <- grep('asite\\[', names, value = TRUE)
atreeid <- grep('atreeid\\[', names, value = TRUE)
bsp <- grep('bsp\\[', names, value = TRUE)

zasp <- grep('zasp\\[', names, value = TRUE)
ypred <- grep('ypred', names, value = TRUE)

other_var <- setdiff(names[-c(1:8)], c(asp, asite, atreeid, bsp, zasp, ypred)) # lp

a <- rep(NA, 8000)
b <- rep(NA, 8000)
for(i in 1:8000){
  iter <- scan(con, what = 'character', nlines = 1, sep = ',', quiet = TRUE)
  a[i] <- iter[3]
  b[i] <- iter[2]
}
close(con)

a_num <- as.numeric(a)
b_num <- as.numeric(b)

mean(a_num) # 1.42973, ~ 0.7 below simulated a
mean(b_num) # 0.1408094, ~ 0.26 below simulated b

plot(a_num, b_num, col = rep(c('blue', 'red', 'green', 'yellow'), each = 2000),
     xlab = 'a', ylab = 'b')
