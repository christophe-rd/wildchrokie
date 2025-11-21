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
real<lower=0> sigma_asp;	 
real<lower=0> sigma_asite;    
real<lower=0> sigma_atreeid;  
real<lower=0> sigma_y;
vector[Nspp] asp; 		
vector[Nsite] asite;       
vector[Ntreeid] atreeid;
}

transformed parameters{
array[N] real ypred;
for (i in 1:N){ 
    ypred[i]=
        a + 
        asp[species[i]] + 
        asite[site[i]] + 
        atreeid[treeid[i]];
    }
}

model{	
  asp ~ normal(0, sigma_asp);
  asite ~ normal(0, sigma_asite);
  atreeid ~ normal(0, sigma_atreeid); 
  a ~ normal(1.5, 1);
  sigma_asp ~ normal(0, 0.3);
  sigma_asite ~ normal(0, 1);
  sigma_atreeid ~ normal(0, 0.1);
  sigma_y ~ normal(0, 1);
  y ~ normal(ypred, sigma_y); 
}	
