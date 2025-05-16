// OSPREE analysis
// 3-level model for budburst day a function of forcing temperature, chilling units, photoperiod in a meta-analysis of 100+ studies
// Level: Species and population on INTERCEPTS and SLOPES, just forcing and photoperiod for now since chilling covarys with site

data {
   // Define variables in data
   // Number of level-1 observations (an integer)
   int<lower=0> N;
   
   // Number of level-2 clusters
   int<lower=0> n_sp;
   
   // Number of level-3 clusters
   int<lower=0> n_pop;
   
   //Cue parameters
   vector[N] mst;
   //vector[N] mwt;
 
   // Cluster IDs
   int<lower=1, upper=n_sp> sp[N];
   int<lower=1, upper=n_pop> pop[N]; 
 
   // Continuous outcome
   real y[N];
   
 }
 
 parameters {
   // Define parameters to estimate
   // Population intercept (a real number)
   real mu_a_sp;
   // Population slope
   real mu_b_mst_sp;
   //real mu_b_mwt_sp;
 
   // Level-1 errors
   real<lower=0> sigma_y;
   real<lower=0> sigma_a_sp;
 
   // Level-2 random effect
   real<lower=0> sigma_b_mst_sp;
   //real<lower=0> sigma_b_mwt_sp;
 
   // Level-3 random effect
   // Population slope
   real<lower=0> sigma_b_mst_sppop;
   //real<lower=0> sigma_b_mwt_sppop;
   real<lower=0> sigma_a_pop;
   
   // Varying intercepts
   vector[n_pop] a_sppop_raw; // reparameterize here
   vector[n_pop] b_mst_sppop_raw; // reparameterize here
   //vector[n_pop] b_mwt_sppop_raw; // reparameterize here
 
   // Individual mean
   vector[n_sp] a_sp_raw; // reparameterize here
   vector[n_sp] b_mst_raw; // reparameterize here
   //vector[n_sp] b_mwt_raw; // reparameterize here
   
   real alpha; // 'grand mean' ... needed when you have more than one level
   
 }
 
 transformed parameters  {
   vector[n_sp] b_mst = mu_b_mst_sp + sigma_b_mst_sp * b_mst_raw;
   //vector[n_sp] b_mwt = mu_b_mwt_sp + sigma_b_mwt_sp * b_mwt_raw;
   
   vector[n_sp] a_sp = mu_a_sp + sigma_a_sp * a_sp_raw;

   vector[n_pop] a_sppop0 = sigma_a_pop * a_sppop_raw; // You need to seperate out the population effect from the species effect mean to make things easier for the indexing
   vector[n_pop] a_sppop;
   
   vector[n_pop] b_mst_sppop0 = sigma_b_mst_sppop * b_mst_sppop_raw; 
   vector[n_pop] b_mst_sppop;
   
   //vector[n_pop] b_mwt_sppop0 = sigma_b_mwt_sppop * b_mwt_sppop_raw; 
   //vector[n_pop] b_mwt_sppop;
   
   for (j in 1:n_pop){ //
     a_sppop[j] = a_sp[sp[j]] + a_sppop0[j];
     b_mst_sppop[j] = b_mst[sp[j]] + b_mst_sppop0[j];
     ///b_mwt_sppop[j] = b_mwt[sp[j]] + b_mwt_sppop0[j];
   }
   
   
 }
 
 model {
   vector[N] yhat;
   // Individual mean
   for(i in 1:N){
            yhat[i] = alpha + a_sppop[pop[i]] + // indexed with population 
		b_mst_sppop[pop[i]] * mst[i];// +
		//b_mwt_sppop[pop[i]] * mwt[i];;
			     	}
   
   
   // Random effects distribution of raw (ncp) priors, they should always be (0, 1)
   target += normal_lpdf(to_vector(a_sp_raw) | 0, 1);
   target += normal_lpdf(to_vector(a_sppop_raw) | 0, 1);
  
   target += normal_lpdf(to_vector(b_mst_raw) | 0, 1);
   //target += normal_lpdf(to_vector(b_mwt_raw) | 0, 1);
   
   target += normal_lpdf(to_vector(b_mst_sppop_raw) | 0, 1);
   //target += normal_lpdf(to_vector(b_mwt_sppop_raw) | 0, 1);
   
   // Random effects distribution of remaining priors
   target += normal_lpdf(to_vector(a_sp) | 0, 20);
   target += normal_lpdf(to_vector(a_sppop) | 0, 20);
   
	target += normal_lpdf(to_vector(b_mst) | 0, 20);
   target += normal_lpdf(to_vector(b_mst_sppop) | 0, 20);
   //target += normal_lpdf(to_vector(b_mwt) | 0, 20);
   //target += normal_lpdf(to_vector(b_mwt_sppop) | 0, 20);
   
   target += normal_lpdf(sigma_y | 0, 10);
 
   // Likelihood part of Bayesian inference
   for (i in 1:N) {
     target += normal_lpdf(y[i] | yhat[i], sigma_y);
   }
 }
