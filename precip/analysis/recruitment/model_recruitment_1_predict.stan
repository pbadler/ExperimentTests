// Null model for recruitment
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> Nspp; // number of species 
  int<lower=0> spp; // focal species id
  vector[N] parents1; // parents in plot
  vector[N] parents2; // parents in group

  // for out of sample prediction 
  int<lower=0> npreds;
  int<lower=0> nyrs_out; // years out 
  int<lower=0> yid_out[npreds]; // year out id
  int<lower=0> gid_out[npreds]; // group id
  int<lower=0> y_holdout[npreds]; // observation vector
  vector[npreds] parents1_out; // hold out parents in plot 
  vector[npreds] parents2_out; // hold out parents in group
  
}parameters{
  real a_mu;
  vector[Yrs] a;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
  vector[N] mu;
  vector[N] trueP1;
  vector[N] trueP2;
  vector[N] lambda;
  vector[N] q;

  trueP1 <- parents1*u + parents2*(1-u);

  for(n in 1:N){
    mu[n] <- exp(a[yid[n]] + gint[gid[n]]);
    lambda[n] <- trueP1[n]*mu[n];  // elementwise multiplication  
  } 
  
  q <- lambda*theta;

}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,5);
  a_mu ~ normal(0,5);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
generated quantities{
  
  vector[nyrs_out] a_out;
  vector[npreds] trueP1_pred;
  vector[npreds] trueP2_pred;
  
  vector[npreds] mu_pred;
  vector[npreds] lambda_hat;
  vector[npreds] qpred;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  int<lower=0> y_hat[npreds]; // pointwise predictions  

  trueP1_pred <- parents1_out*u + parents2_out*(1-u);

  for( i in 1:nyrs_out)
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 

  for(n in 1:npreds){
    mu_pred[n] <- exp(a_out[yid_out[n]] + gint[gid_out[n]]);
    lambda_hat[n] <- trueP1_pred[n]*mu_pred[n];  // elementwise multiplication 
  }
  
  qpred <- lambda_hat*theta;
  
  for(n in 1:npreds){
    y_hat[n] <- neg_binomial_2_rng(qpred[npreds],  theta);
    log_lik[n] <- neg_binomial_2_log(y_holdout[npreds], qpred[npreds], theta);
  }
}
