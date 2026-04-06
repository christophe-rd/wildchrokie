// Two-level (1 hierarchical grouping) `random' intercept model
// Partial pooling on intercepts 
// Updated for new version of Stan (2025!)

data{
int<lower=0> N; 	// number of total observations
int<lower=0> Nspp; 	// number of species (grouping factor)
array[N] int species; 	// species identity, coded as int
int<lower=0> Nyear;  
array[N] int year;  
vector[N] climpredictor; 	// climpredictor (predictor changes according to loop position)
array[N] real y;
}

parameters{
real a;		// mean intercept across everything
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
vector[Nspp] aspp;
vector[Nyear] ayear;
vector[Nspp] bsp;
}

transformed parameters{
array[N] real ypred;
for (i in 1:N){ // don't change this for reparameterization
    ypred[i]=
        a + 
        aspp[species[i]] + 
        ayear[year[i]] + 
        bsp[species[i]]*climpredictor[i];
    }
}

model{	
  a ~ normal(0, 20);
  aspp ~ normal(0, 20);
  ayear ~ normal(0, 15);
  bsp ~ normal(0, 5);
  sigma_y ~ normal(0, 5);
  
  y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
}	

generated quantities {
  // posterior predictive samples
  array[N] real y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(
        a +
        aspp[species[i]] +
        ayear[year[i]] + 
        bsp[species[i]]*climpredictor[i], sigma_y);
  }

  // prior predictive samples
  real a_prior = normal_rng(0, 20);
  real aspp_prior = normal_rng(0, 20);
  real ayear_prior = normal_rng(0, 15);
  real bsp_prior = normal_rng(0, 5);
  real sigma_y_prior = abs(normal_rng(0, 5));
}
