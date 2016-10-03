// Full model for recruitment: includes climate + intraspecific + interspecific competition effects 
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> Nspp;
  int<lower=0> spp;
  
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  
  int<lower=0> Covs; // climate covariates
  matrix[N,Covs] C; // climate matrix
  real tau_beta;
  
}parameters{
  real a_mu;
  vector[Yrs] a;
  vector[Nspp] dd;
  vector[Covs] b2;

  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
  real mu[N];
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] q;
  vector[N] climEff;
  vector[N] parentEff; 
  
  climEff <- C*b2;
  
  for(n in 1:N){
    for( j in 1:Nspp) { 
      trueP1[n, j] <- parents1[n, j]*u + parents2[n, j]*(1-u);
      
      trueP2[n, j] <- sqrt(trueP1[n, j]);
    }
  }
  
  parentEff <- trueP1*dd;
  
  for(n in 1:N){
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + parentEff[n] + climEff[n]);
    
    lambda[n] <- trueP1[n, spp]*mu[n];
      
    q[n] <- lambda[n]*theta;
  }
  
}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,10);
  a_mu ~ normal(0,10);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  dd ~ normal(0, 10);
  
  b2 ~ normal(0, tau_beta);
  
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
