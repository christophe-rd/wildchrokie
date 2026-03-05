// Two-level (1 hierarchical grouping) `random' intercept model
// Partial pooling on intercepts 
// Updated for new version of Stan (2025!)

data{
int<lower=0> N; 	 
int<lower=0> Nspp; 
array[N] int species; 	
int<lower=0> Ntreeid; 
array[N] int treeid;  
array[N] real y;
}

parameters{
real a;		
real<lower=0> sigma_atreeid;
real<lower=0> sigma_y; 	  
vector[Ntreeid] zatreeid; 
vector[Nspp] aspp;
}

transformed parameters{
vector[Ntreeid] atreeid;
atreeid = 0 + sigma_atreeid*zatreeid;
    
array[N] real ypred;
for (i in 1:N){ 
    ypred[i] =
        a + 
        aspp[species[i]] + 
        atreeid[treeid[i]];  
        }
}

model{	
  a ~ normal(5, 3);
  zatreeid ~ normal(0, 1);
  aspp ~ normal(0, 6);
  sigma_atreeid ~ normal(0, 0.5); 
  sigma_y ~ normal(0, 3);
  y ~ normal(ypred, sigma_y); 
}	

generated quantities {
  array[N] real y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(
        a + 
        aspp[species[i]] + 
        atreeid[treeid[i]], sigma_y);  
  }
}
