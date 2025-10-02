// 1 October 2025
// By Christophe
// Goal is to try to attempt to do some kind of prior predictive checks using a new stan model that simulates only from the priors and without data


// https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html#prior-predictive-checks

data {
  int<lower=0> N;
  vector[N] x;
}
generated quantities {
  real b = normal(0, 1);
  real a  = normal(0, 1);
  real sigma_bsp; = normal(0, 1);
  real sigma_asp; = normal(0, 1);
  real sigma_asite;  = normal(0, 1);
  real sigma_atreeid; = normal(0, 1);
  real sigma_y;  = normal(0, 1);
  array[N] real y_sim = normal(alpha + beta * x);
}

