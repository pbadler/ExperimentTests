// Single species climate model for growth: includes climate + intraspecific effects 
data{
  int<lower=0> N; // observations
  int<lower=0,upper=1> Y[N]; // observation vector
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  vector[N] X; // size vector
  int<lower=0> Wcovs;
  matrix[N, Wcovs] W; // crowding effects
  int<lower=0> Covs; // climate covariates
  matrix[N,Covs] C; // climate matrix
  real<lower=0> tau_beta; // prior sdev for climate effects

  // for out of sample predictions 
  int<lower=0> npreds;
  int<lower=0, upper=1> y_holdout[npreds];  
  int<lower=0> nyrs_out; // years out 
  int<lower=0> yid_out[npreds]; //year out id
  int<lower=0> gid_out[npreds]; // group id holdout
  vector[npreds] Xhold;
  matrix[npreds,Wcovs] Whold; // crowding matrix for holdout data 
  matrix[npreds,Covs] Chold;

  // for year effect estimation using entire dataset 
  int<lower=0> N2; // all observations
  int<lower=0,upper=1> Y2[N2]; // observation vector
  int<lower=0> Yrs2; // all years
  int<lower=0> yid2[N2]; // year id
  int<lower=0> gid2[N2]; // group id 
  vector[N2] X2; // size vector
  matrix[N2,Wcovs] W2; // crowding matrix

}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  vector[Wcovs] w;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
  
  // for year effects model  
  real a_mu2;
  vector[Yrs2] a2;
  real b1_mu2;
  vector[Yrs2] b12;
  real<lower=0> sig_a2;
  real<lower=0> sig_b12;
  real gint2[G];
  vector[Wcovs] w2;
  real<lower=0> sig_G2;
  
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  vector[N] crowdEff;
  real mu2[N2] ; 
  vector[N2] crowdEff2; 

  climEff <- C*b2;
  crowdEff <- W*w;

  for(n in 1:N){
    mu[n] <- inv_logit( a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n]);
  }
  
  // for year effects model 
  
  crowdEff2 <- W2*w2;
  
  for(n in 1:N2){
    mu2[n] <- inv_logit(a2[yid2[n]] + gint2[gid2[n]] + b12[yid2[n]]*X2[n] + crowdEff2[n]);
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
  b2 ~ normal(0, tau_beta);
  for(g in 1:G)
    gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ binomial(1,mu);
  
  //for year effects model 
  a_mu2 ~ normal(0,10);
  w2 ~ normal(0,10);
  b1_mu2 ~ normal(0,10);
  sig_a2 ~ cauchy(0,2);
  sig_b12 ~ cauchy(0,2);
  sig_G2 ~ cauchy(0,2);
  for(g in 1:G)
    gint2[g] ~ normal(0, sig_G2);
  for(y in 1:Yrs){
    a2[y] ~ normal(a_mu2, sig_a2);
    b12[y] ~ normal(b1_mu2, sig_b12);
  }

  // Likelihood
  Y2 ~ binomial(1,mu2);
}
generated quantities {
  vector[npreds] climpred;
  vector[npreds] crowdhat;
  vector[nyrs_out] a_out;
  vector[nyrs_out] b1_out;
  real muhat[npreds];
  int<lower=0,upper=1> y_hat[npreds]; // pointwise predictions  
  vector[npreds] log_lik; // vector for computing log pointwise predictive density  
  
  // for year predictions from year effects model 
  real muhat2[npreds];
  int<lower=0,upper=1> y_hat2[npreds]; // pointwise predictions  
  vector[npreds] log_lik2; // vector for computing log pointwise predictive density  
  int<lower=0> yid_out2[npreds]; //integer for modern year effects  
  
  // 1. Holdout data predictions 

  climpred <- Chold*b2;
  crowdhat <- Whold*w;
  
  for( i in 1:nyrs_out){
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1); //draw random year x size effect 
  }
  
  for(n in 1:npreds){
    muhat[n] <- inv_logit(a_out[yid_out[n]] + gint[gid_out[n]] + b1_out[yid_out[n]]*Xhold[n] + crowdhat[n]);
    y_hat[n] <- bernoulli_rng(muhat[n]);
    log_lik[n] <- bernoulli_log(y_holdout[n], muhat[n]);
  }
  
  // 2. Predictions for holdout data with KNOWN year effects.  
  //    Simulate predictions as if year effects in the out of sample data are known. 
  
  for( n in 1:npreds){ 
    yid_out2[n] <- yid_out[n] + Yrs;  // add number of training years to get correct index for a2 and b12
    muhat2[n] <- inv_logit(a2[yid_out2[n]] + gint[gid_out[n]] + b12[yid_out2[n]]*Xhold[n] + crowdhat[n]);
    y_hat2[n] <- bernoulli_rng(muhat2[n]);
    log_lik2[n] <- bernoulli_log(y_holdout[n], muhat2[n]);
  }
}
