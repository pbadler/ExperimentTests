// Full model for recruitment: includes climate + intraspecific + interspecific competition effects 
data{
  int<lower=0> N; // observations
  int<lower=0> Nspp;
  int<lower=0> npreds;
  int<lower=0> spp; // species id  
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  
  int<lower=0> nyrs_out; // years out 
  int<lower=0> yid_out[npreds]; //year out id
  
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> gid_out[npreds]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> y_holdout[npreds]; // observation vector
  matrix[N,Covs] C; // climate matrix
  matrix[npreds,Covs] Chold; // climate matrix, holdout
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  matrix[npreds, Nspp] parents1_out; // hold out parents in plot 
  matrix[npreds, Nspp] parents2_out; // hold out parents in group
  real<lower=0> tau_beta; // prior sd
}parameters{
  real a_mu;
  vector[Yrs] a;
  vector[Covs] b2;
  vector[Nspp] dd;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
real mu[N];
vector[N] climEff;
vector[N] coverEff;
vector[N] parentsEff;
matrix[N, Nspp] trueP1;
matrix[N, Nspp] trueP2;
vector[N] lambda;
vector[N] q;

climEff <- C*b2;    

for(n in 1:N){
  
  trueP1[n, ] <- parents1[n, ]*u + parents2[n, ]*(1-u);
  
  for( j in 1:Nspp){
    trueP2[n, j] <- sqrt(trueP1[n, j]);
  }
  
  coverEff[n] <- dot_product(trueP2[n, ], dd);
  
  mu[n] <- exp(a[yid[n]] + gint[gid[n]] + coverEff[n] + climEff[n]);
  
  lambda[n] <- trueP1[n, spp]*mu[n];
  
  q[n] <- lambda[n]*theta;
  
  }
  
}
model{
// Priors
  u ~ uniform(0,1);
  theta ~ uniform(0.001,100);
  a_mu ~ normal(0,1000);
  sig_a ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  for(g in 1:G){
    gint[g] ~ normal(0, sig_G);
  }
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
  }
  for(j in 1:Covs){
    b2[j] ~ normal(0, tau_beta); 
  }
  for(j in 1:Nspp){
    dd[j] ~ uniform(-100, 100);
  }
  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
generated quantities {
  
  vector[nyrs_out] a_out;
  vector[npreds] climpred;
  vector[npreds] coverEffpred;
  vector[npreds] parentsEffpred;
  matrix[npreds, Nspp] trueP1_pred;
  matrix[npreds, Nspp] trueP2_pred;
  
  real mu_pred[npreds];
  vector[npreds] lambda_hat;
  vector[npreds] qpred;
  vector[npreds] log_lik; // vector for computing log pointwise predictive density
  int<lower=0> y_hat[npreds]; // pointwise predictions  

  climpred <- Chold*b2;    
  
  for( i in 1:nyrs_out){
    a_out[i] <- normal_rng(a_mu, sig_a); // draw random year intercept 
  }
  
  for(n in 1:npreds){
    
    trueP1_pred[n, ] <- parents1_out[n, ]*u + parents2_out[n, ]*(1-u);
    
    for( j in 1:Nspp){  
      trueP2_pred[n, j] <- sqrt(trueP1_pred[n, j]);  // apply sqrt transformation to vector
    }
    
    coverEffpred[n] <- trueP2_pred[n, ]*dd;
    
    mu_pred[n] <- exp(a_out[yid_out[n]] + gint[gid_out[n]] + coverEffpred[n] + climpred[n]);
    
    lambda_hat[n] <- trueP1_pred[n, spp]*mu_pred[n];
    
    qpred[n] <- lambda_hat[n]*theta;
    
    y_hat[n] <- neg_binomial_2_rng(qpred[n],  theta);
    
    log_lik[n] <- neg_binomial_2_log(y_holdout[n], qpred[n], theta);
  }
}
