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
int<lower=0> Nyear;
array[N] int year; 
vector[N] covariate; 	// covariate (predictor for slope)
array[N] real y;
}

parameters{
real a;		// mean intercept across everything
real<lower=0> sigma_atreeid;
real<lower=0> sigma_asite;
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
vector[Ntreeid] zatreeid; // variation of intercept across tree ids, no-centered
vector[Nspp] aspp;
vector[Nsite] zasite;
vector[Nyear] ayear;
vector[Nspp] bsp;
}

transformed parameters{
vector[Ntreeid] atreeid;
atreeid = 0 + sigma_atreeid*zatreeid; // non-centered parameterization on atreeid

vector[Nsite] asite;
asite = 0 + sigma_asite*zasite; // non-centered parameterization on atreeid


array[N] real ypred;
for (i in 1:N){ // don't change this for reparameterization
    ypred[i]=
        a + 
        aspp[species[i]] + 
        asite[site[i]] + 
        atreeid[treeid[i]] + 
        ayear[year[i]] +
        bsp[species[i]]*covariate[i];

    }
}

model{	
  a ~ normal(2, 8);
  zatreeid ~ normal(0, 1); // this creates the partial pooling on intercepts for tree ids, standard sigma for non-centered parameterization
  zasite ~ normal(0, 1);
  aspp ~ normal(0, 10);
  ayear ~ normal(0, 2);
  bsp ~ normal(0, 1.6);
  sigma_atreeid ~ normal(0, 2); 
  sigma_asite ~ normal(0, 2); 
  sigma_y ~ normal(0, 2);
  
  y ~ normal(ypred, sigma_y); // this creates an error model where error is normally distributed
}	

generated quantities {
  // posterior predictive samples
  array[N] real y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(
        a + 
        aspp[species[i]] + 
        asite[site[i]] +
        atreeid[treeid[i]] + 
        ayear[year[i]] +
        bsp[species[i]]*covariate[i], sigma_y);
  }

  // prior predictive samples
  real a_prior = normal_rng(2, 8);
  real sigma_asite_prior = abs(normal_rng(0, 2));
  real sigma_atreeid_prior = abs(normal_rng(0, 2));  
  real sigma_y_prior = abs(normal_rng(0, 2));    
  real aspp_prior = normal_rng(0, 10);
  real ayear_prior = normal_rng(0, 1.5);
  real bsp_prior = normal_rng(0, 1.6);

  real zatreeid_prior = normal_rng(0, 2);
  real atreeid_prior = abs(normal_rng(0, 1)) * zatreeid_prior;
  real zasite_prior = normal_rng(0, 2);
  real asite_prior = abs(normal_rng(0, 1)) * zasite_prior;
}

