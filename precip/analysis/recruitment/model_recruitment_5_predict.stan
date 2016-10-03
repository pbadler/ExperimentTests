// Full model for recruitment: includes climate + intraspecific + interspecific competition effects 
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> Nspp; // number of species 
  int<lower=0> spp; // focal species id
  
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  
  int<lower=0> Covs; // climate covariates
  matrix[N,Covs] C; // climate matrix
  real tau_beta;
  
  int<lower=0> npreds;
  int<lower=0> nyrs_out; // years out 
  int<lower=0> yid_out[npreds]; //year out id
  int<lower=0> gid_out[npreds]; // group id

  int<lower=0> y_holdout[npreds]; // observation vector
  matrix[npreds,Covs] Chold; // climate matrix, holdout

  matrix[npreds, Nspp] parents1_out; // hold out parents in plot 
  matrix[npreds, Nspp] parents2_out; // hold out parents in group
  
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
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] q;
  vector[N] climEff;
  vector[N] coverEff;

  climEff <- C*b2;
    
  trueP1 <- parents1*u + parents2*(1-u);

  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  coverEff <- trueP2*dd;

  for(n in 1:N)
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + coverEff[n] + climEff[n]);
    
  lambda <- trueP1[, spp] .* mu;  // note use of elementwise multiplication operator 
  q <- lambda*theta;
}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,5);
  a_mu ~ normal(0,5);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  dd ~ normal(0, 5);
  
  b2 ~ normal(0, tau_beta);
  
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
generated quantities{
  
  vector[nyrs_out] a_out;
  vector[npreds] climEffpred;
  vector[npreds] coverEffpred;
  matrix[npreds, Nspp] trueP1_pred;
  matrix[npreds, Nspp] trueP2_pred;
  
  vector[npreds] mu_pred;
  vector[npreds] lambda_hat;
  vector[npreds] qpred;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  int<lower=0> y_hat[npreds]; // pointwise predictions  

  climEffpred <- Chold*b2;
    
  trueP1_pred <- parents1_out*u + parents2_out*(1-u);

  for(n in 1:npreds)
    for( j in 1:Nspp)
      trueP2_pred[n, j] <- sqrt(trueP1_pred[n, j]);
  
  coverEffpred <- trueP2_pred*dd;

  for( i in 1:nyrs_out)
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 

  for(n in 1:npreds)
    mu_pred[n] <- exp(a_out[yid_out[n]] + gint[gid_out[n]] + coverEffpred[n] + climEffpred[n]);
    
  lambda_hat <- trueP1_pred[, spp] .* mu_pred;  // note use of elementwise multiplication operator 
  qpred <- lambda_hat*theta;
  
  for(n in 1:npreds){
    y_hat[n] <- neg_binomial_2_rng(qpred[n],  theta);
    log_lik[n] <- neg_binomial_2_log(y_holdout[n], qpred[n], theta);
  }
}
