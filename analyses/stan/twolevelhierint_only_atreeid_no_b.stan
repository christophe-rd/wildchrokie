// Two-level (1 hierarchical grouping) `random' intercept model
// Partial pooling on intercepts 
// Updated for new version of Stan (2025!)

data{
int<lower=0> N; 	// number of total observations
int<lower=0> Ntreeid;
array[N] int treeid; 
int<lower=0> Nspp; 	
array[N] int species;
array[N] real y; 		// day of year of phenological event (response)
}

parameters{
real a;		// mean intercept across everything
real<lower=0> sigma_atreeid;
real<lower=0> sigma_asp;
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
vector[Ntreeid] zatreeid;
vector[Nspp] asp;
}

transformed parameters{
vector[Ntreeid] atreeid;
atreeid = 0 + sigma_atreeid*zatreeid; // non-centered a_sp

// vector[Nspp] asp;
// asp = 0 + sigma_asp*zasp;

array[N] real ypred;
for (i in 1:N){ // don't change this for reparameterization
    ypred[i]=
        a + 
        atreeid[treeid[i]] +
        asp[species[i]];
    }
}

model{	
  zatreeid ~ normal(0, 1);
  a ~ normal(5, 1);
  asp ~ normal(0, sigma_asp);
  sigma_atreeid ~ normal(0, 0.3);
  sigma_asp ~ normal(0, 0.3);
  sigma_y ~ normal(0, 1);
  
  y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
}	
