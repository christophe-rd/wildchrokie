# Started January 2, 2023

# Aim of this code is write test data and test whether Cat's or my new model return values

# We have trait data across different populations

# basing it off code for simulations from traitors

options(mc.cores = 4)

## Set parameters
param <- list(
  Nrep = 20, # rep per trait
  Npop = 10, # number of studies with traits
  Nspp = 15, # number of species with traits
  trait_mu_grand = 20,
  trait_sigma_sp = 4,
  trait_sigma_pop = 2,
  trait_sigma_traity = 3)

## Generate species and study offsets (traits)
mu_sp <- rnorm(n = param[["Nspp"]], mean = 0, sd = param[["trait_sigma_sp"]])
mu_pop <- rnorm(n = param[["Npop"]], mean = 0, sd = param[["trait_sigma_pop"]])


## Make empty table
yTable <- data.frame(Species = c(),
                     Pop = c(),
                     Replicate = c(),
                     mu_grand = c(),
                     mu_sp = c(),
                     mu_pop = c())
## Fill table with parameters

for(i in 1:param[["Nspp"]]){
  for(j in 1:param[["Npop"]]){        
    temp <- data.frame(Species = i,
                       Pop = j,
                       Replicate = 1:param[["Nrep"]],
                       mu_grand = param[["trait_mu_grand"]],
                       mu_sp = mu_sp[i],
                       mu_pop = mu_pop[j])
    yTable <- rbind(yTable, temp)
  }
}

## Calculate latent trait values
yTable$trait_latent <- yTable$mu_grand + yTable$mu_sp + yTable$mu_pop
## Generate trait observation using parameters
yTable$yTraiti <- rnorm(n = nrow(yTable),
                        mean = yTable$trait_latent,
                        sd = param[["trait_sigma_traity"]])

## Calculate response parameters

head(yTable)
hist(yTable$yTraiti)

# Now compare the test data to the different stan models:
all_data <- list(y = yTable$yTraiti,
                 N = nrow(yTable),
                 n_sp = param[["Nspp"]], 
                 sp = yTable$Species,
                 n_pop = param[["Npop"]],
                 pop = yTable$Pop,
                 prior_mu_grand_mu = param[["trait_mu_grand"]],
                 prior_mu_grand_sigma = 5,
                 prior_sigma_sp_mu = param[["trait_sigma_sp"]],
                 prior_sigma_sp_sigma = 5,
                 prior_sigma_study_mu = param[["trait_sigma_pop"]],
                 prior_sigma_study_sigma = 5,
                 prior_sigma_traity_mu = param[["trait_sigma_traity"]],
                 prior_sigma_traity_sigma = 5)

mdl.trait <- stan("stan/trait_3levelwpop.stan",
                      data = all_data,
                      iter = 4000,
                      warmup = 3000,
                      chains = 4,
                      include = FALSE, pars = c("y_hat"))

summary(mdl.trait)$summary
## True value
param


#########################################
# now running it with Cat's model:
mdl.cat <- stan("stan/nointer_3levelwpop_pheno_ncp.stan",
                      data = all_data,
                      iter = 4000,
                      warmup = 3000,
                      chains = 4,
                      include = FALSE, pars = c("y_hat"))

summary(mdl.cat)$summary
