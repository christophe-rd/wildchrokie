// Two-level (1 hierarchical grouping) `random' intercept model
// Partial pooling on intercepts 
// Updated for new version of Stan (2025!)

data{
int<lower=0> N; 	// number of total observations
int<lower=0> Nspp; 	// number of species (grouping factor)
array[N] int species; 	// species identity, coded as int
int<lower=0> Nsite;  // number of sites (grouping factor)
array[N] int site;   // site identity, coded as int
int<lower=0> Ntreeid;  // number of tree ids (grouping factor)
array[N] int treeid;   // tree id identity, coded as int
vector[N] gdd; 	// gdd (predictor for slope)
array[N] real y;
}

parameters{
real a;		// mean intercept across everything
real<lower=0> sigma_atreeid;
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
vector[Ntreeid] zatreeid; // variation of intercept across tree ids, no-centered
vector[Nspp] aspp;
vector[Nsite] asite;
vector[Nspp] bsp;
}

transformed parameters{
vector[Ntreeid] atreeid;
atreeid = 0 + sigma_atreeid*zatreeid; // non-centered parameterization on atreeid

array[N] real ypred;
for (i in 1:N){ // don't change this for reparameterization
    ypred[i]=
        a + 
        aspp[species[i]] + 
        asite[site[i]] + 
        atreeid[treeid[i]] + 
        bsp[species[i]]*gdd[i];

    }
}

model{	
  a ~ normal(5, 3);
  zatreeid ~ normal(0, 1); // this creates the partial pooling on intercepts for tree ids, standard sigma for non-centered parameterization
  aspp ~ normal(0, 6);
  asite ~ normal(0, 2);
  bsp ~ normal(0, 0.5);
  sigma_atreeid ~ normal(0, 0.5); 
  sigma_y ~ normal(0, 3);
  
  y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
}	

generated quantities {
  array[N] real y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(
        a + 
        aspp[species[i]] + 
        asite[site[i]] +
        atreeid[treeid[i]] + 
        bsp[species[i]]*gdd[i], sigma_y);
  }
}
