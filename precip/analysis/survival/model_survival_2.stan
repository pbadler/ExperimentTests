data{
  int<lower=0> N; // observations
  int<lower=0,upper=1> Y[N]; // observation vector
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  vector[N] X; // size vector
  int<lower=0> gid[N]; // group id
  int<lower=0> Wcovs; // number of crowding effects 
  matrix[N,Wcovs] W; // crowding matrix
  int<lower=0>Covs; // number of climate effects 
  matrix[N,Covs] C; // climate matrix
  real tau_beta;
  int<lower=0> spp; // focal species number 
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  real w;
  real gint[G];
  real<lower=1e-7> sig_a;
  real<lower=1e-7> sig_b1;
  real<lower=1e-7> sig_G;
  vector[Covs] b2;
}
transformed parameters{
  real mu[N];
  vector[N] crowdEff;
  vector[N] climEff;
  vector[N] W_intra; 
  
  W_intra <- W[,spp];
  
  crowdEff <- W_intra*w;
  climEff <- C*b2;
  
  for(n in 1:N){
    mu[n] <- inv_logit( a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n]);
  }
}
model{
  // Priors
  a_mu ~ normal(0,10);
  w ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  for(g in 1:G)
    gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }
  b2 ~ normal(0, tau_beta);
  
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


