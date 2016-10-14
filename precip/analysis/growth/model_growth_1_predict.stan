// Single species growth model: includes intraspecific effects only 
data{
  int<lower=0> N; // observations
  vector[N] Y; // observation vector
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  vector[N] X; // size vector
  vector[npreds] Xhold;
  int<lower=0> Wcovs; // number of crowding effects 
  matrix[N,Wcovs] W; // crowding matrix
  int<lower=0>Covs; // number of climate effects 
  matrix[N,Covs] C; // climate matrix
  real tau_beta;

  // For out of sample prediction
  int<lower=0> npreds;
  vector[npreds] y_holdout;
  int<lower=0> nyrs_out; // years holdout 
  int<lower=0> yid_out[npreds];// year id holdout
  int<lower=0> gid_out[npreds]; // group id holdout
  matrix[npreds,Wcovs] Whold; // crowding matrix for holdout data 
  matrix[npreds,Covs] Chold;
  
  // For predicting overall cover 
  
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Wcovs] w;
  vector[Covs] b2;
  real gint[G];
  real tau;
  real tauSize;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  vector[N] mu;
  real<lower=0> sigma[N];
  vector[N] crowdEff;
  vector[N] climEff;
  
  
  climEff <- C*b2;
  crowdEff <- W*w;

  for(n in 1:N){
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n];
    sigma[n] <- sqrt((fmax(tau*exp(tauSize*mu[n]), 0.0000001)));  
  }
}
model{
  // Priors
  a_mu ~ normal(0,10);
  b1_mu ~ normal(0,10);
  w ~ normal(0, 10);
  tau ~ normal(0,10);
  tauSize ~ normal(0,10);
  sig_a ~ cauchy(0,2);
  sig_b1 ~ cauchy(0,2);
  b2 ~ normal(0, tau_beta);
  sig_G ~ cauchy(0,2);
  for(g in 1:G)
    gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {
  vector[nyrs_out] a_out;
  vector[nyrs_out] b1_out;
  vector[npreds] crowdhat;
  vector[npreds] sigmahat;
  vector[npreds] muhat;
  vector[npreds] y_hat;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  
  crowdhat <- Whold*w;
  
  for( i in 1:nyrs_out){
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 
    b1_out[i] <- normal_rng(b1_mu, sig_b1); //draw random year x size effect 
  }
  for(n in 1:npreds){
    muhat[n] <- a_out[yid_out[n]] + gint[gid_out[n]] + b1_out[yid_out[n]]*Xhold[n] + crowdhat[n];
    sigmahat[n] <- sqrt((fmax(tau*exp(tauSize*muhat[n]), 0.0000001))); 
    y_hat[n] <- normal_rng(muhat[n], sigmahat[n]);
    log_lik[n] <- normal_log(y_holdout[n], muhat[n], sigmahat[n]);
  }
  
  
  
}

