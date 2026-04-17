// year model
data{
int<lower=0> N; 	// number of total observations
int<lower=0> Nyear;
array[N] int year19;
array[N] int year20;
// array[N] int year; 
array[N] real y;
}

parameters{
real a;		
real byear19;
real byear20;
real<lower=0> sigma_y; 	// measurement error, noise etc. 	
}

transformed parameters{

array[N] real ypred;
for (i in 1:N){ 
    ypred[i]=
        a + 
        byear19 * year19[i] +
        byear20 * year20[i];
    }
}

model{	
  a ~ normal(2, 4);
  byear19 ~ normal(0, 3);
  byear20 ~ normal(0, 3);
  sigma_y ~ normal(0, 3);
  y ~ normal(ypred, sigma_y); 
}	

generated quantities {
  // posterior predictive samples
  array[N] real y_rep;
  for (i in 1:N) {
    y_rep[i] = normal_rng(
        a + 
        byear19 * year19[i] +
        byear20 * year20[i],
        sigma_y);
  }

  // prior predictive samples
  real a_prior = normal_rng(2, 4);
  real byear19_prior = normal_rng(0, 3);
  real byear20_prior = normal_rng(0, 3);
  real sigma_y_prior = abs(normal_rng(0, 1));

  // // For LOO cross-validation
  // vector[N] log_lik;
  // for (i in 1:N) {
  //   log_lik[i] = normal_lpdf(y[i] | ypred[i], sigma_y);
  // }
}
