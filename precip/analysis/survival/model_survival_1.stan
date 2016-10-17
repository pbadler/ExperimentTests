data{
  int<lower=0> N; // observations
  int<lower=0,upper=1> Y[N]; // observation vector
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  vector[N] X; // size vector
  int<lower=0> Wcovs; // number of crowding effects 
  matrix[N,Wcovs] W; // crowding matrix
  int<lower=0> spp; // focal species 
}
parameters{
  real gint_mu; 
  vector[G] gint_raw;
  vector[Yrs] a_raw;
  real b1_mu;
  vector[Yrs] b1_raw;
  real w;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  vector[Yrs] a;
  vector[Yrs] b1;
  vector[G] gint; 
  real mu[N];
  vector[N] crowdEff;
  vector[N] W_intra; 

  W_intra <- W[, spp];
  crowdEff <- W_intra*w;
  
  // reparamaterize the hierarchical parameters  
  a <- 0 + sig_a*a_raw;
  b1 <- b1_mu + sig_b1*b1_raw;
  gint <- gint_mu + sig_G*gint_raw;
  
  for(n in 1:N){
    mu[n] <- inv_logit(gint[gid[n]] + a[yid[n]]  + b1[yid[n]]*X[n] + crowdEff[n]);
  }
  
}
model{
  // Priors
  gint_mu ~ normal(0,10);
  w ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,2);
  sig_G ~ cauchy( 0,5);
  gint_raw ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  b1_raw ~ normal(0, 1);
  
  // Likelihood
  Y ~ binomial(1,mu);

}
generated quantities {

  // Section for calculating log_lik of fitted data 
  
  vector[N] log_lik; 
  
  for(n in 1:N){
    log_lik[n] <- bernoulli_log(Y[n], mu[n]); 
  }
  
}



