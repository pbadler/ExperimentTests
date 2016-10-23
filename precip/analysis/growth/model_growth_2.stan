// Intra-specific competition model for survival: includes intraspecific effects only 
data{
  int<lower=0> N; // observations
  vector[N] Y; // observation vector
  int<lower=0> nyrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  vector[N] X; // size vector
  int<lower=0> Wcovs; // number of crowding effects 
  vector[N] W; // crowding vector
  int<lower=0>Covs; // number of climate effects 
  matrix[N,Covs] C; // climate matrix
  real tau_beta;
}
parameters{
  real a_mu;
  vector[nyrs] a;
  real b1_mu;
  vector[nyrs] b1;
  real w;
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
  for(y in 1:nyrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, sigma);
}
generated quantities {

  // Section for calculating log_lik of fitted data 
  
  vector[N] log_lik; 
  
  for(n in 1:N){
    log_lik[n] <- normal_log(Y[n], mu[n] , sigma[n]); 
  }
  
}

