// Estimate GDD with leafout data for wildchrokie
// CRD on 17 November 2025
// see issue # 13 for more details

data{
int<lower=0> N;
int<lower=0> Nspp; 
array[N] int species;
int<lower=0> Nsite;
array[N] int site;
int<lower=0> Ntreeid;
array[N] int treeid;  
array[N] real y; 		
}

parameters{
real a;
real<lower=0> sigma_atreeid;  
real<lower=0> sigma_y;
vector[Nspp] aspp; 		
vector[Nsite] asite;
vector[Ntreeid] zatreeid;
}

transformed parameters{
vector[Ntreeid] atreeid;
atreeid = 0 + sigma_atreeid*zatreeid;
  
array[N] real ypred;
for (i in 1:N){ 
    ypred[i]=
        a + 
        aspp[species[i]] + 
        asite[site[i]] +
        atreeid[treeid[i]];
    }
}

model{	
  // aspp ~ normal(0, sigma_aspp);
  // asite ~ normal(0, sigma_asite);
  zatreeid ~ normal(0, 1); 
  a ~ normal(8, 6);
  aspp ~ normal(0, 3);
  asite ~ normal(0, 3);
  sigma_atreeid ~ normal(0, 1.5);
  sigma_y ~ normal(0, 3);
  y ~ normal(ypred, sigma_y); 
}	

generated quantities{
  array[N] real y_tilde;
  for (i in 1:N){
    y_tilde[i] = normal_rng(ypred[i], sigma_y);
  }
}
