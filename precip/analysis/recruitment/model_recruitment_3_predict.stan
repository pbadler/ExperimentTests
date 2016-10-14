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
  
  // for out of sample prediction 
  int<lower=0> npreds;
  int<lower=0> nyrs_out; // years out 
  int<lower=0> yid_out[npreds]; // year out id
  int<lower=0> gid_out[npreds]; // group id
  int<lower=0> y_holdout[npreds]; // observation vector
  matrix[npreds,Covs] Chold; // climate matrix, holdout
  matrix[npreds, Nspp] parents1_out; // hold out parents in plot 
  matrix[npreds, Nspp] parents2_out; // hold out parents in group
  
  // for year effect estimation using entire dataset 
  int<lower=0> N2; // all observations
  int<lower=0> Y2[N2]; // observation vector
  int<lower=0> Yrs2; // all years
  int<lower=0> yid2[N2]; // year id
  int<lower=0> gid2[N2]; // group id 
  matrix[N2, Nspp] parents1_2; // parents in plot
  matrix[N2, Nspp] parents2_2; // parents in group
  
}parameters{
  real a_mu;
  vector[Yrs] a;
  vector[Nspp] w;
  vector[Covs] b2;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
  
  // for year effects model  
  real a_mu2;
  vector[Yrs2] a2;
  vector[Nspp] w2;
  real gint2[G];
  real<lower=0> sig_a2;
  real<lower=0> theta2;
  real<lower=0> sig_G2;
  real<lower=0, upper=1> u2;
}
transformed parameters{
  vector[N] mu;
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] lambda;
  vector[N] q;
  vector[N] climEff;
  vector[N] coverEff;

  // for year effects model 
  vector[N2] mu2; 
  matrix[N2, Nspp] trueP1_2;
  matrix[N2, Nspp] trueP2_2;
  vector[N2] lambda2;
  vector[N2] q2; 
  vector[N2] coverEff2;

  climEff <- C*b2;
  trueP1 <- parents1*u + parents2*(1-u);
  
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j]);
  
  coverEff <- trueP2*w;

  for(n in 1:N){
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + coverEff[n] + climEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  // elementwise multiplication  
  } 
  
  q <- lambda*theta;
  
  // for year effects model 
  
  trueP1_2 <- parents1_2*u2 + parents2_2*(1-u2);

  for(n in 1:N2)
    for(j in 1:Nspp)
      trueP2_2[n, j] <- sqrt(trueP1_2[n, j]);
  
  coverEff2 <- trueP2_2*w2;

  for(n in 1:N2){
    mu2[n] <- exp(a2[yid2[n]] + gint2[gid2[n]] + coverEff2[n]);
    lambda2[n] <- trueP1_2[n, spp]*mu2[n];  // elementwise multiplication  
  } 
  
  q2 <- lambda2*theta2;

}
model{
  // Priors
  u ~ uniform(0,1);
  theta ~ uniform(0,5);
  a_mu ~ normal(0,5);
  sig_a ~ cauchy(0,2);
  sig_G ~ cauchy(0,2);
  w ~ normal(0, 5);
  b2 ~ normal(0, tau_beta);
  gint ~ normal(0, sig_G);
  a ~ normal(a_mu, sig_a);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
  
  // For year effects model 
  u2 ~ uniform(0,1);
  theta2 ~ uniform(0,5);
  a_mu2 ~ normal(0,5);
  sig_a2 ~ cauchy(0,2);
  sig_G2 ~ cauchy(0,2);
  w2 ~ normal(0, 5);
  gint2 ~ normal(0, sig_G2);
  a2 ~ normal(a_mu2, sig_a2);

  // Likelihood
  Y2 ~ neg_binomial_2(q2, theta2);
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

  // for year predictions from year effects model 
  vector[npreds] qpred2;
  vector[npreds] mu_pred2;
  vector[npreds] lambda_hat2;
  int<lower=0> y_hat2[npreds]; // pointwise predictions  
  vector[npreds] log_lik2; // vector for computing log pointwise predictive density  
  int<lower=0> yid_out2[npreds]; //integer for modern year effects  

  // 1. Holdout data predictions 
  climEffpred <- Chold*b2;
    
  trueP1_pred <- parents1_out*u + parents2_out*(1-u);

  for(n in 1:npreds)
    for( j in 1:Nspp)
      trueP2_pred[n, j] <- sqrt(trueP1_pred[n, j]);
  
  coverEffpred <- trueP2_pred*w;

  for( i in 1:nyrs_out)
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 

  for(n in 1:npreds){
    mu_pred[n] <- exp(a_out[yid_out[n]] + gint[gid_out[n]] + coverEffpred[n] + climEffpred[n]);
    lambda_hat[n] <- trueP1_pred[n, spp]*mu_pred[n];  // elementwise multiplication 
  }
  
  qpred <- lambda_hat*theta;
  
  for(n in 1:npreds){
    y_hat[n] <- neg_binomial_2_rng(qpred[npreds],  theta);
    log_lik[n] <- neg_binomial_2_log(y_holdout[npreds], qpred[npreds], theta);
  }
  
  // 2. Predictions for holdout data with KNOWN year effects.  
  //    Simulate predictions as if year effects in the out of sample data are known. 
  
  for(n in 1:npreds){
    yid_out2[n] <- yid_out[n] + Yrs;  // add number of training years to get correct index for a2 and b12
    mu_pred2[n] <- exp(a2[yid_out2[n]] + gint[gid_out[n]] + coverEffpred[n]);
    lambda_hat2[n] <- trueP1_pred[n, spp]*mu_pred2[n];  // elementwise multiplication 
  }
  
  qpred2 <- lambda_hat2*theta;
  
  for(n in 1:npreds){
    y_hat2[n] <- neg_binomial_2_rng(qpred2[npreds],  theta);
    log_lik2[n] <- neg_binomial_2_log(y_holdout[npreds], qpred2[npreds], theta);
  }  
  
}
