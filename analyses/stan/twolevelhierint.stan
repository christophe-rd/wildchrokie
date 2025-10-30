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
array[N] real y; 		// day of year of phenological event (response)
}

parameters{
real b;        // slope
//real a;		// mean intercept across everything
real<lower=0> sigma_bsp;
real<lower=0> sigma_asp;	// variation of intercept across species	
real<lower=0> sigma_asite;    // variation of intercept across sites
real<lower=0> sigma_atreeid;    // variation of intercept across tree ids
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
vector[Nspp] bsp;
vector[Nspp] zasp; 		// defining transformed asp
vector[Nsite] asite;       //the intercept for each sites
vector[Ntreeid] atreeid;       //the intercept for each tree id
}

transformed parameters{
vector[Nspp] asp;
asp = 0 + sigma_asp*zasp; // non-centered a_sp
array[N] real ypred;
for (i in 1:N){ // don't change this for reparameterization
    ypred[i]=
        //a + 
        asp[species[i]] + 
        asite[site[i]] + 
        atreeid[treeid[i]] + 
        b*gdd[i] +
      bsp[species[i]]*gdd[i];
    }
}

model{	
  bsp ~ normal(0, sigma_bsp); // I guess partial pooling on slopes for species
  asite ~ normal(0, sigma_asite); // this creates the partial pooling on intercepts for sites
  atreeid ~ normal(0, sigma_atreeid); // this creates the partial pooling on intercepts for tree ids
  //a ~ normal(4, 1);
  b ~ normal(0, 0.2);
  sigma_bsp ~ normal(0, 1);
  sigma_asp ~ normal(0, 0.3);
  zasp ~ normal(0, 1); // here i put the standard centered prior on zasp
  sigma_asite ~ normal(0, 1);
  sigma_atreeid ~ normal(0, 0.1);
  sigma_y ~ normal(0, 1);
  
  y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
}	
