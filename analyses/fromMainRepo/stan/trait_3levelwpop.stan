

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
  
  mu_grand ~ normal(0, 10);
  sigma_sp ~ normal(0, 10);
  sigma_pop ~ normal(0, 10);
  sigma_y ~ normal(0, 10);

  y ~ normal(y_hat, sigma_y);
}

