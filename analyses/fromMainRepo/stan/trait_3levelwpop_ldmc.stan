

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  
  int<lower=0> n_sp;
   
   // Number of level-3 clusters
   int<lower=0> n_pop;
 
   // Cluster IDs
   int<lower=1, upper=n_sp> sp[N];
   int<lower=1, upper=n_pop> pop[N]; 
 
   // Continuous outcome
   real y[N];
   
  real prior_mu_grand_mu; 
  real prior_mu_grand_sigma;
  real prior_sigma_sp_mu;
  real prior_sigma_sp_sigma;
  real prior_sigma_pop_mu;
  real prior_sigma_pop_sigma;
  real prior_sigma_y_mu;
  real prior_sigma_y_sigma;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu_grand;
  vector[n_sp] mu_sp;
  vector[n_pop] mu_pop;
  
  real<lower=0> sigma_y;
  real<lower=0> sigma_pop;
  real<lower=0> sigma_sp;
}

transformed parameters{
  vector[N] y_hat;
  
  for (i in 1:N){
    y_hat[i] = mu_grand + mu_sp[sp[i]] + mu_pop[pop[i]];
  }
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  mu_sp ~ normal(0, sigma_sp);
  mu_pop ~ normal(0, sigma_pop);
  
  mu_grand ~ normal(prior_mu_grand_mu, prior_mu_grand_sigma);
  sigma_sp ~ normal(prior_sigma_sp_mu, prior_sigma_sp_sigma);
  sigma_pop ~ normal(prior_sigma_pop_mu, prior_sigma_pop_sigma);
  sigma_y ~ normal(prior_sigma_y_mu, prior_sigma_y_sigma);

  y ~ normal(y_hat, sigma_y);
}

