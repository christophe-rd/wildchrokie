//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
// Two-level (1 hierarchical grouping) `random' slope and intercept model
// Partial pooling on intercepts and slopes 

data{
int<lower=0> N;     // number of total observations
int<lower=0> Nspp;  // number of species (grouping factor)
int species[N];     // species identity, coded as int
vector[N] year;     // year of data point (predictor for slope)
real y[N];      // day of year of phenological event (response)
}

parameters{
real mu_a;      // mean intercept across species
real<lower=0> sigma_a;  // variation of intercept across species    
real mu_b;      // mean slope across species
real<lower=0> sigma_b;  // variation of slope across species
real<lower=0> sigma_y;  // measurement error, noise etc.    
real a[Nspp];       //the intercept for each species
real b[Nspp];       //the slope for each species 

}

transformed parameters{
real ypred[N];
for (i in 1:N){
    ypred[i]=a[species[i]]+b[species[i]]*year[i];
}
}

model{  
b ~ normal(mu_b, sigma_b); // this creates the partial pooling on slopes
a ~ normal(mu_a, sigma_a); // this creates the partial pooling on intercepts
y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
// Priors ...
mu_a ~ normal(100,30);
sigma_a ~ normal(0,20);
mu_b ~ normal(0,5);
sigma_b ~ normal(0,15);
sigma_y ~ normal(0,15);
}   
