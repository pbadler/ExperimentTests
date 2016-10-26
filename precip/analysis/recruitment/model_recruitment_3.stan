// Full multi-species model with climate: includes climate + competition effects
data{
  int<lower=0> N;             // observations
  int<lower=0> Y[N];          // observation vector
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  int<lower=0> G;             // groups
  int<lower=0> gid[N];        // group id
  int<lower=0> Nspp;          // number of species 
  int<lower=0> spp;           // focal species id
  matrix[N, Nspp] parents1;   // parents in plot
  matrix[N, Nspp] parents2;   // parents in group
  int<lower=0> Covs;          // climate covariates
  matrix[N,Covs] C;           // climate matrix
  real tau_beta;              // prior standard deviation
  matrix[N, G] gm;

}parameters{
  vector[nyrs] a_raw;
  vector[Nspp] w;
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0, upper=1> u;
  vector[Covs] b2;
  vector[G] bg; 
}
transformed parameters{
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] coverEff;
  vector[N] climEff;
  vector[N] gint; 
  vector[nyrs] a; 
  
  climEff <- C*b2;
  trueP1 <- parents1*u + parents2*(1-u);
  
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  gint     <- gm*bg;
  coverEff <- trueP2*w;
  a  <- 0 + a_raw*sig_a; 

  for(n in 1:N){
    mu[n] <- exp(gint[n] + a[yid[n]] + coverEff[n] + climEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  
  }
}

model{
  // Priors
  u ~ uniform(0,1);
  theta ~ cauchy(0,2);
  sig_a ~ cauchy(0,2);
  a_raw ~ normal(0, 1);
  bg ~ normal(0, 10);
  w ~ normal(0, 5);
  b2 ~ normal(0, tau_beta);

  // Likelihood
  Y ~ neg_binomial_2(lambda, theta);
}
generated quantities{
  
  vector[N] log_lik; // vector for computing log pointwise predictive density
  vector[N] y_hat; 
  
  for(n in 1:N){ 
    y_hat[n]   <- neg_binomial_rng( lambda[n], theta);
    log_lik[n] <- neg_binomial_2_log(Y[n], lambda[n], theta); 
  }
}
