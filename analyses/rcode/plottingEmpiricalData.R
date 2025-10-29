# Wildchrokie model
# CRD 28 October April 2025
# Start plotting empirical data for wildchrokie

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(max.print = 150) 
options(digits = 3)
# quartz()

# Load library 
library(ggplot2)

setwd("/Users/christophe_rouleau-desrochers/github/wildchrokie/analyses")

sim <- read.csv("output/simdata.csv")
emp <- read.csv("output/empiricalDataMAIN.csv")

# Plot empirical data
ggplot(emp, aes(x= pgsGDD, y = lengthCM, group = spp, color = spp)) +
  # geom_violin() +
  geom_boxplot() +
  facet_wrap(~spp, ncol = 1, nrow = 4) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) + 
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()
ggsave("figures/empiricalData/sppBoxPlots.jpeg", width = 6, height = 6, units = "in", dpi = 300)

# let's just look at 2018
y18 <- subset(emp, year == "2018")

ggplot(y18, aes(x= pgsGDD, y = lengthCM, group = spp, color = spp)) +
  geom_boxplot(position = "dodge2") +
  facet_wrap(~year) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) + 
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()

ggplot(emp, aes(x = pgsGDD, y = lengthCM, color = spp, fill = spp)) +
  geom_point(size = 2, alpha = 0.7) +  # optional: show raw data points
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  # facet_wrap(~spp)+
  theme_minimal()
ggsave("figures/empiricalData/sppLinearRegressions_nofacet.jpeg", width = 6, height = 6, units = "in", dpi = 300)


test <- data.frame(
  # aspp = rnorm(8000, 0, 1),
  aspp1 = rnorm(8000, 0.2, 1),
  aspp2 = rnorm(8000, -0.5, 1),
  aspp3 = rnorm(8000, -1, 1),
  aspp4 = rnorm(8000, 1, 1)
)
test 
simsub <- sim[sample(sample(1:40000, 1e4)),]
simsub <- subset(sim, spp %in% 1:4)

# get some basic stats
simcoef_summary <- aggregate(ringwidth ~ spp, data = simsub, function(x) {
  mean <- mean(x)
  serror <- sd(x) / sqrt(length(x))
  lower <- mean - 1.96 * serror
  upper <- mean + 1.96 * serror
  c(mean = mean, serror = serror, lower = lower, upper = upper)
})

# split
simcoef_summary <- do.call(data.frame, simcoef_summary)
names(simcoef_summary) <- c("spp", "mean_ringwidth", "serror", "lower", "upper")

# mean gdd and merge!
gdd_mean <- aggregate(gddcons ~ spp, data = simcoef, mean)
simcoef_summary <- merge(simcoef_summary, gdd_mean, by = "spp")

ggplot(simcoef_summary, aes(x = gddcons, y = factor(spp))) +
  geom_vline(xintercept = 0, color = "black") +
  
  # geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, color = "grey50") +
  geom_point(size = 3, color = "black") +
  theme_minimal() +
  labs(
    x = "GDD (centered, gddcons)",
    y = "Species",
    title = "ringwidth by spp "
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


test_long <- data.frame(
  spp = rep(names(test), each = nrow(test)),
  a_spp = as.vector(as.matrix(test))
)
test_long

ggplot(test_long, aes(x = a_spp, y = spp, group = spp, color = spp)) +
  geom_boxplot(position = "dodge2", outlier.shape = NA) +
  geom_vline(xintercept = 0, color = "black") +
  scale_color_manual(values = wes_palette("AsteroidCity1")) +
  scale_fill_manual(values = wes_palette("AsteroidCity1")) +
  theme_minimal()
ggsave("figures/aspp_boxPlot.jpeg", width = 6, height = 6, units = "in", dpi = 300)
  
 


