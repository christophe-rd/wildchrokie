// Two-level (1 hierarchical grouping) `random' intercept model
// Partial pooling on intercepts 
// Updated for new version of Stan (2025!)

data{
int<lower=0> N; 	// number of total observations
int<lower=0> Nspp; 	// number of species (grouping factor)
array[N] int species;
int<lower=0> Nsite;
array[N] int site; 
vector[N] gdd; 	// gdd (predictor for slope)
array[N] real y; 		// day of year of phenological event (response)
}

parameters{
real b;        // slope
real a;		// mean intercept across everything
real<lower=0> sigma_asp;
real<lower=0> sigma_asite;// variation of intercept across species	
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
vector[Nspp] zasp;
vector[Nsite] zasite;
}

transformed parameters{
vector[Nspp] asp;
asp = 0 + sigma_asp*zasp; // non-centered a_sp

vector[Nsite] asite;
asite = 0 + sigma_asite*zasite; // non-centered a_sp

array[N] real ypred;
for (i in 1:N){ // don't change this for reparameterization
    ypred[i]=
        a + 
        asp[species[i]]+
        asite[site[i]] +
        b*gdd[i];
    }
}

model{	
  zasp ~ normal(0, 1);
  zasite ~ normal(0, 1);
  a ~ normal(5, 1);
  b ~ normal(0.5, 1);
  sigma_asp ~ normal(0, 0.3);
  sigma_asite ~ normal(0, 0.3);
  sigma_y ~ normal(0, 1);
  
  y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
}	
