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
  int<lower=0> Covs; // climate covariates 
  matrix[N,Covs] C;  // climate matrix
  real tau_beta;     // climate prior sd 
  
  // for out of sample predictions
  int<lower=0> npreds;
  int<lower=0,upper=1> y_holdout[npreds];
  int<lower=0> nyrs_out; // years out
  int<lower=0> yid_out[npreds]; //year out id
  int<lower=0> gid_out[npreds]; // group id holdout
  vector[npreds] Xhold;
  matrix[npreds,Wcovs] Whold; // crowding matrix for holdout data
  matrix[npreds,Covs] Chold;  // climate matrix

  
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
  vector[Covs] b2;
}
transformed parameters{
  real mu[N];
  vector[N] crowdEff;
  vector[N] W_intra; 
  vector[N] climEff; 
  
  W_intra <- W[, spp];
  crowdEff <- W_intra*w;
  climEff  <- C*b2;
  
  for(n in 1:N){
    mu[n] <- inv_logit(a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n]);
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
  b2 ~ normal(0, tau_beta);

  // Likelihood
  Y ~ binomial(1,mu);

}
generated quantities {
  vector[npreds] crowdhat;
  vector[nyrs_out] a_out;
  vector[nyrs_out] b1_out;
  real muhat[npreds];
  int<lower=0,upper=1> y_hat[npreds]; // pointwise predictions  
  vector[npreds] log_lik; // vector for computing log pointwise predictive density  
  vector[npreds] Whold_intra; //
  vector[npreds] climhat; 
  
  // 1. Holdout data predictions 
  Whold_intra <- Whold[, spp];
  
  crowdhat <- Whold_intra*w;
  climhat  <- Chold*b2;
  
  for( i in 1:nyrs_out){
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1); //draw random year x size effect 
  }
  
  for(n in 1:npreds){
    muhat[n] <- inv_logit(a_out[yid_out[n]] + gint[gid_out[n]] + b1_out[yid_out[n]]*Xhold[n] + crowdhat[n]);
    y_hat[n] <- bernoulli_rng(muhat[n]);
    log_lik[n] <- bernoulli_log(y_holdout[n], muhat[n]);
  }
}
