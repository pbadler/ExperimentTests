// Intra-specific competition model for survival: includes intraspecific effects only 
data{
  int<lower=0> N; // observations
  int<lower=0> npreds;
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> gid_out[npreds]; // group id holdout
  
  int<lower=0,upper=1> Y[N]; // observation vector
  
  vector[npreds] y_holdout;
  vector[N] X; // size vector
  vector[npreds] Xhold;
  vector[N] W; // intraspecific crowding 
  vector[npreds] Whold; // intraspecific crowding matrix for holdout data 
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  real w;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  real mu[N];

  vector[N] crowdEff;
  crowdEff <- W*w;
  
  for(n in 1:N){
    mu[n] <- inv_logit(a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n]);
  }
  
}
model{
  // Priors
  a_mu ~ normal(0,10);
  w ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,2);
  sig_b1 ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  for(g in 1:G)
    gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

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



